"""OCV from GITT (plus HPPC 2 h rests to fill gaps), and the candidate dU/dT.

Steps (all choices are explained in README.md):
1. Equilibrium voltage of each rest = average of every sample where |dV/dt| < theta, with dV/dt
   from a trailing straight-line fit whose length makes the slope's standard error <= theta/5.
2. SOC = 1 - dQ/Q(T), dQ counted from the rested state after the slow (C/20) charge to 3.45 V.
3. Discharge curve per cell = GITT discharge points + HPPC 2 h rest voltages of the same cell (15/45 C);
   at 25 C (10 % GITT spacing) the 15/45 C discharge shape bent through the 25 C GITT points.
   Charge curve = discharge curve + charge/discharge gap (measured by GITT at 15/45 C; at 25 C taken
   from 15/45 C and corrected to the measured 25 C gap at 10-90 % SOC).
4. Shape-preserving (PCHIP) interpolation onto 0.01 SOC; branches meet at their mean at SOC 1; at SOC 0
   both equal the rested voltage after the normal-rate discharge to 1.5 V; mean of cells.
   dU/dT = least-squares slope of the per-cell curves vs (T - 25 C).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize_scalar

from . import load

MV_PER_H = 1e-3 / 3600.0            # 1 mV/h in V/s
GRID = np.round(np.arange(0, 1.0 + 1e-9, 0.01), 2)


# ---------------------------------------------------------------- relaxation criterion
def trailing_slope(t, v, W, min_pts=3):
    """Least-squares slope of v(t) over the samples in [t_k - W, t_k]; NaN until the window is full."""
    t = np.asarray(t, float) - t[0]
    v = np.asarray(v, float) - v[0]
    c = lambda x: np.r_[0.0, np.cumsum(x)]
    St, Sv, Stt, Stv = c(t), c(v), c(t * t), c(t * v)
    hi = np.arange(1, len(t) + 1)
    lo = np.searchsorted(t, t - W - 1e-6, side="left")
    n = hi - lo
    st, sv, stt, stv = St[hi] - St[lo], Sv[hi] - Sv[lo], Stt[hi] - Stt[lo], Stv[hi] - Stv[lo]
    den = n * stt - st * st
    with np.errstate(invalid="ignore", divide="ignore"):
        slope = (n * stv - st * sv) / den
    span = t - t[np.minimum(lo, len(t) - 1)]
    slope[(n < min_pts) | (span < 0.9 * W) | (den <= 0)] = np.nan
    return slope


def window_for(theta, sigma_v, dt, w_max=4e5):
    """Shortest window (multiple of dt) whose slope standard error is <= theta/5 at noise sigma_v."""
    for n in range(3, int(w_max / dt) + 2):
        if sigma_v / (dt * np.sqrt(n * (n * n - 1) / 12.0)) <= theta / 5.0:
            return (n - 1) * dt
    return np.inf


def relaxed_point(t, v, theta, W):
    """Equilibrium voltage = average of every sample with |slope| < theta.
    Returns (index of the first such sample, average voltage, number of samples averaged);
    (None, nan, 0) if no sample qualifies (unrelaxed)."""
    s = trailing_slope(t, v, W)
    ok = np.flatnonzero(np.abs(s) < theta)
    if len(ok) == 0:
        return None, np.nan, 0
    return ok[0], float(np.mean(np.asarray(v)[ok])), len(ok)


def logging_class(t):
    d = np.diff(t)
    dt = np.median(d[len(d) // 2:]) if len(d) > 2 else np.nan
    return "300s" if dt > 100 else ("1s" if dt > 0.5 else "0.1s")


def window_s(logging, cfg):
    rc = cfg["relaxation"]
    dt = {"0.1s": 0.1, "1s": 1.0, "300s": 300.0}[logging]
    W = window_for(rc["theta_mV_per_h"] * MV_PER_H, rc["noise_uV"][logging] * 1e-6, dt)
    return max(W, (rc["min_samples_300s"] - 1) * 300.0) if logging == "300s" else W


def temp_rise(t, Tc, block_s):
    """Range of block means (the probes have ~0.3-0.4 K sample-to-sample noise)."""
    m = pd.Series(Tc).groupby((t // block_s).astype(int)).mean()
    return float(m.max() - m.min())


def equilibrium(df, r, cfg, T_setpoint):
    t, v, d = load.rest_trace(df, r)
    lg = logging_class(t)
    W = window_s(lg, cfg)
    k, V, n_avg = relaxed_point(t, v, cfg["relaxation"]["theta_mV_per_h"] * MV_PER_H, W)
    s = trailing_slope(t, v, W)[-1]
    out = dict(seg=int(r.seg_first), logging=lg, W_s=W, rest_h=t[-1] / 3600, paused_h=d.pause_s.sum() / 3600,
               relaxed=k is not None, V_eq=V if k is not None else float(v[-1]),
               t_eq_h=(t[k] if k is not None else t[-1]) / 3600,
               samples_averaged=n_avg, slope_at_end_mV_h=s / MV_PER_H if np.isfinite(s) else np.nan, V_end=float(v[-1]),
               q_ah=float(r.q_ah), prev_I=float(r.prev_I))
    if "aux_T1a" in d:
        Tc = d[["aux_T1a", "aux_T2a"]].mean(axis=1).to_numpy()
        rise = temp_rise(t, Tc, cfg["temperature"]["rise_block_s"])
        out.update(T_label=float(Tc.mean()), T_source="logged", T_rise_K=rise,
                   rejected_T=bool(rise > cfg["temperature"]["rise_cap_K"]))
    else:
        out.update(T_label=float(T_setpoint), T_source="setpoint", T_rise_K=np.nan, rejected_T=False)
    return out


# ---------------------------------------------------------------- GITT branches and HPPC rests
def gitt_branches(df, steps, top="c20"):
    """Charge and discharge branch rests of a GITT file and the 0 % / 100 % anchors.

    The conditioning leg is the last CP step before the first CC step. A conditioning discharge to
    1.5 V means the charge branch comes first (15/45 C files); a conditioning charge means the
    discharge branch comes first (25 C files). Each branch = its start rest + every rest after a CC
    step in that direction. 100 % ('c20') = rest after the last CC charge that reached 3.45 V;
    'cp' = rest after the last CP (conditioning-rate) charge to 3.45 V (used only to measure where a
    normal-rate charge ends)."""
    R = load.merged_rests(df, steps)
    cc = steps[(steps["mode"] == "CC") & (steps.kind != "rest")]
    first_cc = cc.index[0]
    cond = steps[(steps["mode"] == "CP") & (steps.index < first_cc)].iloc[-1]
    cond_rest = R[R.seg_first == cond.name + 1].iloc[0]
    first_dir = "chg" if cond.kind == "dchg" else "dchg"
    second_dir = "dchg" if first_dir == "chg" else "chg"
    gR = R[R.index >= cond_rest.name]
    b1 = [cond_rest.name] + gR[(gR.prev_kind == first_dir) & (gR.prev_seg > cond.name)].index.tolist()
    after = gR[gR.index > b1[-1]]
    b2 = [b1[-1]] + after[after.prev_kind == second_dir].index.tolist()
    branch = {first_dir: R.loc[b1].copy(), second_dir: R.loc[b2].copy()}
    if cond.kind == "dchg":
        anchor0 = cond_rest
    else:
        pre = steps[(steps["mode"] == "CP") & (steps.kind == "dchg") & (steps.index < cond.name)].iloc[-1]
        anchor0 = R[R.seg_first == pre.name + 1].iloc[0]
    v_max = steps[steps.kind == "chg"].V_end.max()
    c20 = steps[(steps["mode"] == "CC") & (steps.kind == "chg") & (steps.V_end >= v_max - 2e-3)].index[-1]
    if top == "c20":
        top_seg = c20
    else:
        top_seg = steps[(steps["mode"] == "CP") & (steps.kind == "chg") & (steps.V_end >= v_max - 2e-3)
                        & (steps.index < first_cc)].index[-1]
    anchor100 = R[R.seg_first == top_seg + 1].iloc[0]
    q_c20 = float(R[R.seg_first == c20 + 1].iloc[0].q_ah)
    return dict(branch=branch, anchor0=anchor0, anchor100=anchor100, Q_test=float(anchor0.q_ah - q_c20))


def hppc_rests(df, steps, cfg, T):
    """Equilibrium points of the HPPC 2 h rests. 'post-move' = after a CC discharge SOC move.
    q100 = charge count at the first 2 h rest (after the top-up charge to 3.45 V)."""
    R = load.merged_rests(df, steps)
    long = R[(R.dur_s >= 7000) & (R.dur_s < 10000)]
    rows = []
    for _, r in long.iterrows():
        e = equilibrium(df, r, cfg, T)
        e["kind"] = "top" if r.name == long.index[0] else ("post-move" if r.prev_dur > 15 else "post-pulse")
        rows.append(e)
    return pd.DataFrame(rows), float(long.iloc[0].q_ah)


# ---------------------------------------------------------------- curves
def grid_branch(z, v, grid=GRID):
    """Shape-preserving interpolation through (z, v), evaluated on grid (extrapolates at the ends)."""
    z = np.asarray(z, float); v = np.asarray(v, float)
    o = np.argsort(z)
    z, v = z[o], v[o]
    keep = np.r_[True, np.diff(z) > 1e-9]
    return PchipInterpolator(z[keep], v[keep], extrapolate=True)(grid)


def monotonic_violations(u):
    return [float(GRID[i + 1]) for i in np.flatnonzero(np.diff(u) <= 0)]


def slope_fit(T, U, T_ref=25.0):
    """Per-column least squares U = a + b (T - T_ref). Returns (a, b, standard error of b)."""
    X = np.c_[np.ones_like(T), T - T_ref]
    coef, *_ = np.linalg.lstsq(X, U, rcond=None)
    s2 = ((U - X @ coef) ** 2).sum(axis=0) / max(len(T) - 2, 1)
    return coef[0], coef[1], np.sqrt(s2 * np.linalg.inv(X.T @ X)[1, 1])


def align_offset(z_raw, v, z_ref, v_ref):
    """SOC offset that puts (z_raw, v) onto the curve through (z_ref, v_ref), fitted at 15-85 % SOC."""
    u = grid_branch(z_ref, v_ref)
    m = (z_raw > 0.15) & (z_raw < 0.85)
    f = lambda dz: float(np.sum((v[m] - np.interp(z_raw[m] + dz, GRID, u)) ** 2))
    r = minimize_scalar(f, bounds=(-0.03, 0.03), method="bounded")
    return r.x, 1e3 * np.sqrt(f(0) / m.sum()), 1e3 * np.sqrt(r.fun / m.sum())


# ---------------------------------------------------------------- build
def build(cfg):
    """Return (ocv table, entropy table, info) in the datapack format."""
    Q = {int(k): v for k, v in cfg["soc"]["capacity_ah"].items()}
    pts, anchor0 = [], []
    for T, names in load.GITT.items():
        for name in names:
            df, steps = load.load_cached(load.gitt_path(name))
            b = gitt_branches(df, steps, "c20")
            a0 = equilibrium(df, b["anchor0"], cfg, T)                  # rested after the normal-rate discharge to 1.5 V
            anchor0.append(dict(T=T, cell=load.cell_id(name), V=a0["V_eq"],
                                z=1.0 - (a0["q_ah"] - float(b["anchor100"].q_ah)) / Q[T]))
            for br, R in b["branch"].items():
                for _, r in R.iterrows():
                    e = equilibrium(df, r, cfg, T)
                    e.update(T=T, cell=load.cell_id(name), file=name, branch=br, Q_test_ah=b["Q_test"],
                             z=1.0 - (e["q_ah"] - float(b["anchor100"].q_ah)) / Q[T])
                    pts.append(e)
    pts = pd.DataFrame(pts)
    anchor0 = pd.DataFrame(anchor0)
    use = pts[~pts.rejected_T]
    cells = {T: sorted(pts[pts["T"] == T].cell.unique()) for T in load.GITT}
    gitt_only = {(T, c, br): grid_branch(g.z, g.V_eq) for (T, c, br), g in use.groupby(["T", "cell", "branch"])}

    # HPPC 2 h rests of the GITT cells, SOC axis shifted so they sit on the cell's own GITT discharge curve
    hp, align = [], []
    for T, names in load.HPPC.items():
        for name in names:
            cid = load.cell_id(name)
            if cid not in cells[T]:
                continue
            df, steps = load.load_cached(load.hppc_path(name))
            h, q100 = hppc_rests(df, steps, cfg, T)
            h["z_raw"] = 1.0 - (h.q_ah - q100) / Q[T]
            g = use[(use["T"] == T) & (use.cell == cid) & (use.branch == "dchg")]
            pm = h[h.kind == "post-move"]
            dz, before, after = align_offset(pm.z_raw.to_numpy(), pm.V_eq.to_numpy(), g.z, g.V_eq)
            h["z"] = h.z_raw + dz
            h["T"], h["cell"], h["file"] = T, cid, name
            align.append(dict(T=T, cell=cid, soc_offset_pct=100 * dz, rms_before_mV=before, rms_after_mV=after))
            hp.append(h)
    hp = pd.concat(hp, ignore_index=True)

    # 25 C charge/discharge gap: from 15/45 C (interpolated in T), corrected to the measured 25 C gap at 10-90 %
    mean_c = lambda T, br: np.mean([gitt_only[(T, c, br)] for c in cells[T]], axis=0)
    T15 = float(use[use["T"] == 15].T_label.mean())
    w = (25.0 - T15) / (45.0 - T15)
    gap_T = (mean_c(15, "chg") - mean_c(15, "dchg")) * (1 - w) + (mean_c(45, "chg") - mean_c(45, "dchg")) * w
    lo, hi = cfg["ocv"]["gap_correction_soc"]
    zc = np.round(np.arange(lo, hi + 1e-9, 0.1), 2)

    curves, used_hp, ends = {}, [], []
    for T in [15, 25, 45]:
        for cid in cells[T]:
            gd = use[(use["T"] == T) & (use.cell == cid) & (use.branch == "dchg")]
            hh = hp[(hp["T"] == T) & (hp.cell == cid) & (hp.kind == "post-move")]
            hh = hh[[np.min(np.abs(gd.z.to_numpy() - z)) > cfg["ocv"]["hppc_min_gap_soc"] for z in hh.z]]
            used_hp.append(hh)
            zd, vd = np.r_[gd.z, hh.z], np.r_[gd.V_eq, hh.V_eq]
            if T == 25 and cfg["ocv"]["method_25C"] == "shape_15_45":
                zd, vd = gd.z.to_numpy(), gd.V_eq.to_numpy()          # HPPC rests not used at 25 C
            if cfg["ocv"]["soc0_rule"] == "normal_discharge_rest":
                # SOC 0 = rested voltage after the normal-rate discharge to 1.5 V; the slow discharge's points
                # past 0 % are dropped so the discharge curve starts from that state
                v0 = float(anchor0[(anchor0["T"] == T) & (anchor0.cell == cid)].V.iloc[0])
                zd, vd = np.r_[0.0, zd[zd > 0]], np.r_[v0, vd[zd > 0]]
            if T == 25 and cfg["ocv"]["method_25C"] == "shape_15_45":
                # 15/45 C discharge shape, bent through the 25 C points
                shape = mean_c(15, "dchg") * (1 - w) + mean_c(45, "dchg") * w
                ud = shape + grid_branch(zd, vd - np.interp(zd, GRID, shape))
            else:
                ud = grid_branch(zd, vd)
            measured_gap = gitt_only[(T, cid, "chg")] - gitt_only[(T, cid, "dchg")]
            if T == 25:
                resid = np.interp(zc, GRID, measured_gap - gap_T)
                gap = gap_T + PchipInterpolator(zc, resid)(np.clip(GRID, lo, hi))
            else:
                gap = measured_gap
            uc = ud + gap
            ends.append(dict(T=T, cell=cid, chg_soc0=uc[0], dchg_soc0=ud[0], chg_soc1=uc[-1], dchg_soc1=ud[-1]))
            for i in (0, -1):                       # branches meet at their mean at SOC 0 and 1
                uc[i] = ud[i] = 0.5 * (uc[i] + ud[i])
            if cfg["ocv"]["soc0_rule"] == "normal_discharge_rest":
                uc[0] = ud[0] = v0
            curves[(T, cid, "chg")], curves[(T, cid, "dchg")] = uc, ud

    rows, mono = [], []
    for T in [15, 25, 45]:
        uc = np.mean([curves[(T, c, "chg")] for c in cells[T]], axis=0)
        ud = np.mean([curves[(T, c, "dchg")] for c in cells[T]], axis=0)
        for nm, u in [("OCV-Charge", uc), ("OCV-Discharge", ud), ("Average", 0.5 * (uc + ud))]:
            if monotonic_violations(u):
                mono.append(dict(T=T, column=nm, soc=monotonic_violations(u)))
        rows += [[T, z, uc[i], ud[i], 0.5 * (uc[i] + ud[i])] for i, z in enumerate(GRID)]
    ocv_tab = pd.DataFrame(rows, columns=load.OCV_COLS)

    ent_soc = load.read_table(load.DATAPACK / "entropy.csv").SOC.to_numpy()
    ent = {"SOC": ent_soc}
    slopes = {}
    for br, col in [("chg", "Charge"), ("dchg", "Discharge")]:
        ks = [(T, c) for T in [15, 25, 45] for c in cells[T]]
        Tl = np.array([use[(use["T"] == T) & (use.cell == c) & (use.branch == br)].T_label.mean() for T, c in ks])
        a, b, se = slope_fit(Tl, np.array([curves[(T, c, br)] for T, c in ks]))
        ent[col] = np.interp(ent_soc, GRID, b)
        slopes[br] = dict(T_labels={f"{T}C {c}": round(x, 2) for (T, c), x in zip(ks, Tl)},
                          dUdT_V_per_K=b.tolist(), std_error_V_per_K=se.tolist())
    ent_tab = pd.DataFrame(ent)[load.ENTROPY_COLS]

    cols = ["T", "cell", "file", "branch", "seg", "z", "V_eq", "relaxed", "t_eq_h", "samples_averaged", "slope_at_end_mV_h", "V_end",
            "rest_h", "paused_h", "logging", "W_s", "T_label", "T_source", "T_rise_K", "rejected_T"]
    info = dict(
        gitt_points=pts[cols].round(6), hppc_points_used=pd.concat(used_hp)[["T", "cell", "file", "seg", "z_raw", "z", "V_eq", "relaxed"]].round(6),
        hppc_alignment=pd.DataFrame(align), branch_ends_before_meeting=pd.DataFrame(ends), monotonicity_violations=mono,
        q_test_ah={f"{T}C {c}": round(q, 2) for (T, c), q in pts.groupby(["T", "cell"]).Q_test_ah.first().items()},
        entropy_slopes=slopes, T15_label=T15)
    return ocv_tab, ent_tab, info
