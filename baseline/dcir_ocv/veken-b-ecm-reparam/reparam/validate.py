"""Whole-file simulation of the RPT (held-out cell) and HPPC records with a parameter set.

One simulator for every parameter set: the 1RC model of ecm.py, row by row over the entire file.
- SOC: coulomb counting. SOC = 1 at the end of each rested full charge (a rest that follows a charge
  ending at 3.45 V) and counted forward from there; rows before the first full charge are counted
  back from it. Never inferred from voltage. The jump at each re-anchoring is reported (it is the
  coulomb-counting error accumulated since the previous full charge).
- OCV: the branch of the most recent current direction (|I| > 1 A), held through rests.
- Parameters at the file's temperature: table rows used directly; 60 C (outside both sets' data)
  from the 25/45 C OCV rows extended linearly in T and per-SOC Arrhenius through the 15/25/45 C R/C rows.
- The first logged row of every step (captured mid-switch) is left out of the error statistics.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

from . import ecm, load

V_MAX = 3.45


def rc_at(tab, T):
    if T in set(tab["T"]):
        return tab[tab["T"] == T].sort_values("SOC")[["SOC", "R0", "R1", "tau"]].reset_index(drop=True)
    rows = []
    for s, g in tab[tab["T"].isin(ecm.MEASURED)].sort_values(["SOC", "T"]).groupby("SOC"):
        a, b = ecm.arrhenius_fit(g["T"].to_numpy(float), g[["R0", "R1", "tau"]].to_numpy())
        v = ecm.arrhenius_eval(a, b, [T])[0]
        rows.append(dict(SOC=s, R0=v[0], R1=v[1], tau=v[2]))
    return pd.DataFrame(rows)


def ocv_at(ocv_tab, T):
    def row(TT):
        g = ocv_tab[ocv_tab.Temperature == TT].sort_values("SOC")
        return g.SOC.to_numpy(), g["OCV-Charge"].to_numpy(), g["OCV-Discharge"].to_numpy()
    if T in set(ocv_tab.Temperature):
        z, c, d = row(T)
    else:
        z, c25, d25 = row(25)
        _, c45, d45 = row(45)
        k = (T - 25.0) / 20.0
        c, d = c25 + k * (c45 - c25), d25 + k * (d45 - d25)
    return PchipInterpolator(z, c, extrapolate=True), PchipInterpolator(z, d, extrapolate=True)


def soc_full_file(df, steps, Q, soc_at_anchor=1.0):
    """(SOC per row, table of anchors with the SOC jump found at each)."""
    R = load.merged_rests(df, steps)
    anchors = R[(R.prev_kind == "chg") & (R.prev_Vend >= V_MAX - 2e-3) & (R.prev_dur > 15)]   # not 10 s pulses
    q = df.q_ah.to_numpy()
    seg = df.seg.to_numpy()
    end_idx = np.array([np.flatnonzero(seg == s)[-1] for s in anchors.seg_last])
    q_anc = q[end_idx]
    k = np.searchsorted(end_idx, np.arange(len(q)), side="left") - 1    # last anchor ending before the row
    k_use = np.where(k < 0, 0, k)                                        # before the first anchor: count back from it
    z = soc_at_anchor - (q - q_anc[k_use]) / Q
    jumps = [dict(anchor_end_h=float(df.t.iloc[e] / 3600),
                  soc_before_reanchoring=float(soc_at_anchor - (q[e] - q_anc[j - 1]) / Q))
             for j, e in enumerate(end_idx) if j > 0]
    return z, pd.DataFrame(jumps)


def simulate(df, z, rc, ocvs):
    """Terminal voltage for every row (row-by-row exact update; current may change within a step)."""
    t, i = df.t.to_numpy(), df.current.to_numpy()
    zc = np.clip(z, 0.0, 1.0)
    R0, R1, tau = (np.interp(zc, rc.SOC, rc[c]) for c in ["R0", "R1", "tau"])
    s = pd.Series(np.where(i > 1.0, 1.0, np.where(i < -1.0, -1.0, np.nan))).ffill().fillna(1.0).to_numpy()
    U = np.where(s < 0, ocvs[0](z), ocvs[1](z))                      # charging -> charge branch
    a = np.exp(-np.diff(t) / tau[:-1])
    x = np.zeros_like(i)
    for k in range(len(a)):
        x[k + 1] = a[k] * x[k] + (1 - a[k]) * i[k]
    return U - R0 * i - R1 * x


def phases(df, steps, test):
    """Phase name per row, for labelling and statistics."""
    seg = df.seg.to_numpy()
    out = np.empty(len(df), dtype=object)
    if test == "RPT":
        p60 = steps[(steps.kind != "rest") & steps.dur_s.between(50, 70)].index
        first_chg_after = steps[(steps.index > p60[-1]) & (steps.kind == "chg") & (steps.dur_s > 3600)].index[0]
        start = steps[(steps.index < p60[0]) & (steps.kind == "rest")].index[-1]
        out[:] = "capacity cycles to 1.5 V"
        out[seg >= start] = "pulses and CP SOC moves"
        out[seg >= first_chg_after] = "cycles to 1.8 V"
    else:
        top = steps[(steps.kind == "rest") & (steps.dur_s >= 7000) & (steps.dur_s < 10000)].index[0]
        out[:] = "conditioning"
        out[seg >= top] = "HPPC pulses and SOC moves"
    return out


def step_type(steps):
    lab = {}
    for sid, s in steps.iterrows():
        if s.kind == "rest":
            lab[sid] = "rest"
        elif s.dur_s <= 70:
            lab[sid] = "pulses (10-60 s)"
        elif s.dur_s <= 200:
            lab[sid] = "150 s pulses"
        else:
            lab[sid] = f"{s['mode']} {'charge' if s.kind == 'chg' else 'discharge'}"
    return lab


def metrics(e):
    a = np.abs(e)
    return dict(n=len(e), rmse_mV=1e3 * np.sqrt(np.mean(e ** 2)), p99_mV=1e3 * np.quantile(a, 0.99),
                max_mV=1e3 * a.max(), within_50mV=float(np.mean(a <= 0.05)))


def run_file(test, T, path, sets, Q, soc_starts):
    """Simulate one file with every parameter set. Returns (metrics rows, trace dict, anchor jumps)."""
    df, steps = load.load_cached(path)
    first = np.r_[False, df.seg.to_numpy()[1:] != df.seg.to_numpy()[:-1]]
    ph = phases(df, steps, test)
    st = pd.Series(df.seg.map(step_type(steps))).to_numpy()
    rows, trace, jumps = [], dict(t_h=df.t.to_numpy() / 3600, V=df.voltage.to_numpy(), I=df.current.to_numpy(),
                                  phase=ph, first=first), None
    for name, (ocv_tab, rc_tab) in sets.items():
        for start_name, z0 in soc_starts.get(name, {"SOC 1 at full charge": 1.0}).items():
            z, jumps_ = soc_full_file(df, steps, Q, z0)
            if jumps is None:
                jumps = jumps_
            v = simulate(df, z, rc_at(rc_tab, T), ocv_at(ocv_tab, T))
            e = v - df.voltage.to_numpy()
            key = dict(test=test, T=T, file=path.name, set=name, soc_start=start_name)
            keep = ~first
            rows.append(dict(key, part="whole file", **metrics(e[keep])))
            for p in pd.unique(ph):
                rows.append(dict(key, part=p, **metrics(e[keep & (ph == p)])))
            for lab, m in [("SOC < 10 %", z < 0.1), ("SOC 10-90 %", (z >= 0.1) & (z <= 0.9)), ("SOC > 90 %", z > 0.9)]:
                if (keep & m).sum() > 5:
                    rows.append(dict(key, part=lab, **metrics(e[keep & m])))
            for lab in pd.unique(st):
                if (keep & (st == lab)).sum() > 5:
                    rows.append(dict(key, part=lab, **metrics(e[keep & (st == lab)])))
            trace[(name, start_name)] = (v, z)
    return rows, trace, jumps
