"""1RC model from HPPC: R0 from the pulses, one tau1 for every temperature and R1 per SOC from the
SOC moves, gridded onto the team SOC/temperature grid with an Arrhenius fill for the cold rows.

Model (current discharge-positive, exact for piecewise-constant current):
    V_k = OCV(SOC_k, T) - R0 i_k - R1 x_k,     x_{k+1} = a_k x_k + (1 - a_k) i_k,   a_k = exp(-dt_k / tau1)
The current on row k is held until row k+1. The last row of one step and the first row of the next
are the same instant (dt = 0).

HPPC schedule at each SOC point: 2 h rest after the SOC move, then for each amplitude a discharge
pulse, 60 s rest, charge pulse, 900 s rest (2 h after the last charge pulse); the last point runs
charge-first. Fit windows: a pulse window is the rest sample before the pulse + the 10 s pulse + its
rest (capped at 900 s); a move window is the rest sample before the SOC move + the move + its 2 h rest.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.optimize import least_squares, minimize_scalar, nnls

from . import load, ocv

KB_EV = 8.617333262e-5
MEASURED = [15, 25, 45]
SOC_GRID = np.array(load.R_SOC_GRID, float)


# ---------------------------------------------------------------- model
def rc_current(t, i, tau, x0=0.0, step=None):
    """RC-branch current x at every sample. With `step` ids the current is held at each step's mean
    and the recursion is evaluated in closed form (fast, exact for constant-current steps); without,
    row by row (needed for constant-power steps)."""
    t = np.asarray(t, float); i = np.asarray(i, float)
    x = np.empty_like(i)
    if step is None:
        a = np.exp(-np.diff(t) / tau)
        x[0] = x0
        for k in range(len(a)):
            x[k + 1] = a[k] * x[k] + (1.0 - a[k]) * i[k]
        return x
    step = np.asarray(step)
    starts = np.flatnonzero(np.r_[True, step[1:] != step[:-1]])
    xs = x0
    for s0, s1 in zip(starts, np.r_[starts[1:], len(t)]):
        I = i[s0:s1].mean()
        x[s0:s1] = I + (xs - I) * np.exp(-(t[s0:s1] - t[s0]) / tau)
        xs = x[s1 - 1]
    x[0] = x0
    return x


def fit_window(t, i, eta, mask, bounds, weights, n_hist=0, step=None):
    """Fit (R0, R1, tau, c) to one window's overpotential eta = V - OCV.

    Rows [0, n_hist) are history since the last long rest: they drive the RC state into the window
    with the same parameters but carry no residual. c is the window's zero level. Residuals are
    weighted by sqrt(local dt) and taken on `mask` rows. The start is global: for fixed tau the model
    is linear in (R0, R1, c), so tau is scanned on a log grid with an exact solve at each value."""
    t = np.asarray(t, float); i = np.asarray(i, float); eta = np.asarray(eta, float)
    m = np.asarray(mask, bool).copy(); m[:n_hist] = False
    w = np.asarray(weights, float)[m]

    def model(p):
        R0, R1, tau, c = p
        e = -R0 * i - R1 * rc_current(t, i, tau, 0.0, step)
        return e - e[n_hist] + c

    def vp(tau):
        x = rc_current(t, i, tau, 0.0, step)
        A = np.c_[-(i - i[n_hist]), -(x - x[n_hist]), np.ones_like(t), -np.ones_like(t)][m] * w[:, None] * 1e3
        coef, rn = nnls(A, eta[m] * w * 1e3)
        return rn, np.r_[coef[:2], coef[2] - coef[3]]

    lo = np.r_[np.asarray(bounds[0], float), -np.inf]
    hi = np.r_[np.asarray(bounds[1], float), np.inf]
    taus = np.geomspace(lo[2], hi[2], 60)
    k = int(np.argmin([vp(tt)[0] for tt in taus]))
    lt = np.log(taus)
    r = minimize_scalar(lambda q: vp(np.exp(q))[0], bounds=(lt[max(k - 1, 0)], lt[min(k + 1, len(lt) - 1)]),
                        method="bounded", options=dict(xatol=1e-4))
    R0, R1, c = vp(np.exp(r.x))[1]
    best = least_squares(lambda p: (model(p) - eta)[m] * w * 1e3, np.clip([R0, R1, np.exp(r.x), c], lo, hi),
                         bounds=(lo, hi), x_scale="jac", max_nfev=4000)
    p = best.x
    try:
        se = np.sqrt(np.diag(np.linalg.inv(best.jac.T @ best.jac) * (2 * best.cost / max(m.sum() - 4, 1))))
    except np.linalg.LinAlgError:
        se = np.full(4, np.nan)
    resid = model(p) - eta
    return dict(R0=p[0], R1=p[1], tau=p[2], c=p[3], se=se, resid=resid,
                rmse_mV=1e3 * np.sqrt(np.mean(resid[m] ** 2)), max_mV=1e3 * np.max(np.abs(resid[m])))


def solve_given_tau(W, tau):
    """For fixed tau and R0 a move window is linear in (R1, c): returns (weighted SSE, R1, c)."""
    m = W["mask"].copy(); m[:W["nh"]] = False
    w = W["w"][m]
    x = rc_current(W["t"], W["i"], tau, 0.0, W["step"])
    y = W["eta"] + W["R0"] * (W["i"] - W["i"][W["nh"]])
    A = np.c_[-(x - x[W["nh"]]), np.ones_like(x), -np.ones_like(x)][m] * w[:, None] * 1e3
    coef, rn = nnls(A, y[m] * w * 1e3)
    return rn ** 2, coef[0], coef[1] - coef[2]


# ---------------------------------------------------------------- HPPC windows
def structure(steps):
    """One row per pulse / SOC move: SOC-point index, direction, amplitude class."""
    top = steps[(steps.kind == "rest") & (steps.dur_s >= 7000) & (steps.dur_s < 10000)].index[0]
    rows, point = [], 0
    for sid, s in steps[(steps.kind != "rest") & (steps.index > top)].iterrows():
        if s.dur_s > 15:
            point += 1
        rows.append(dict(seg=sid, kind="move" if s.dur_s > 15 else "pulse", point=point, dir=s.kind,
                         I=abs(s.I_mean)))
    st = pd.DataFrame(rows)
    lv = np.sort(st[st.kind == "pulse"].I.to_numpy())
    cuts = lv[:-1][np.diff(np.log(lv)) > np.log(1.2)]          # amplitude classes: gaps wider than 20 %
    st["amp_A"] = st.assign(c=np.where(st.kind == "pulse", np.searchsorted(cuts, st.I, side="left"), -1)) \
        .groupby(["kind", "c"]).I.transform("median").round(1)
    return st


def soc_axis(df, steps, Q, ocv_dchg):
    """SOC per row, counted from the HPPC's first 2 h rest (after its top-up charge), then shifted by
    one offset so the post-move 2 h rest voltages at 15-85 % SOC sit on the OCV discharge curve."""
    R = load.merged_rests(df, steps)
    long = R[(R.dur_s >= 7000) & (R.dur_s < 10000)]
    q100 = float(long.iloc[0].q_ah)
    pm = long[long.prev_dur > 15]
    zr = 1.0 - (pm.q_ah.to_numpy() - q100) / Q
    m = (zr > 0.15) & (zr < 0.85)
    f = lambda dz: float(np.sum((pm.V_end.to_numpy()[m] - ocv_dchg(zr[m] + dz)) ** 2))
    r = minimize_scalar(f, bounds=(-0.05, 0.05), method="bounded")
    return 1.0 - (df.q_ah.to_numpy() - q100) / Q + r.x, dict(soc_offset_pct=100 * r.x,
                                                             rms_before_mV=1e3 * np.sqrt(f(0) / m.sum()),
                                                             rms_after_mV=1e3 * np.sqrt(r.fun / m.sum()))


def window(df, steps, seg_id, kind, rest_cap_s=900.0, long_rest_s=3000.0):
    """(rows, n_hist): history rows since the last long rest, then the window rows."""
    seg = df.seg.to_numpy()
    i0 = int(np.flatnonzero(seg == seg_id - 1)[-1])
    nxt = last = seg_id + 1
    while last + 1 in steps.index and steps.loc[last + 1, "kind"] == "rest":
        last += 1
    win = np.r_[i0, np.flatnonzero(seg == seg_id), np.flatnonzero((seg >= nxt) & (seg <= last))]
    t = df.t.to_numpy()
    if kind == "pulse" and steps.loc[nxt:last, "dur_s"].sum() < long_rest_s:
        win = win[t[win] <= t[np.flatnonzero(seg == nxt)[0]] + rest_cap_s + 1e-9]
    before = steps[(steps.kind == "rest") & (steps.dur_s >= long_rest_s) & (steps.index < seg_id)]
    h0 = int(np.flatnonzero(seg == before.index[-1])[-1]) if len(before) else i0
    return df.iloc[np.r_[np.arange(h0, i0), win]], i0 - h0


def window_arrays(d, nh, z, ocv_fn):
    """Time from window start, current, overpotential, mask (first row of every step is logged
    mid-switch and excluded), sqrt(dt) weights."""
    t = d.t.to_numpy() - d.t.iloc[nh]
    V, ov = d.voltage.to_numpy(), ocv_fn(z)
    eta = (V - ov) - (V[nh] - ov[nh])
    seg = d.seg.to_numpy()
    mask = ~np.r_[False, seg[1:] != seg[:-1]]
    w = np.sqrt(np.gradient(t))
    w[~np.isfinite(w) | (w == 0)] = np.sqrt(0.1)
    return t, d.current.to_numpy(), eta, mask, w, seg


def window_temperature(d, nh, T, t):
    w = d.iloc[nh:]
    if "aux_T1a" not in w:
        return float(T), "setpoint", np.nan
    Tc = w[["aux_T1a", "aux_T2a"]].mean(axis=1).to_numpy()
    return float(Tc.mean()), "logged", ocv.temp_rise(t[nh:], Tc, 60.0 if t[-1] <= 900 else 600.0)


def excluded(fits):
    """Windows that start in a state no other SOC point has: the 100 % point's pulses and the first SOC
    move (fresh off the charge to 3.45 V), and the last point's pulses (charge-first order); plus moves
    left out by name (see fit_hppc)."""
    last = fits.groupby("file").point.transform("max")
    named = fits["left_out"].fillna(False).astype(bool) if "left_out" in fits else False
    return (((fits.kind == "pulse") & ((fits.point == 0) | (fits.point == last)))
            | ((fits.kind == "move") & (fits.point == 1)) | named)


def fit_hppc(cfg, ocv_tab, leave_out=()):
    """Per-window results for every HPPC file: discharge pulses (R0) and SOC moves (R1 at one tau1).
    leave_out: [{T, soc}] moves ending within 1 % of that SOC are reported but kept out of tau1 and the tables."""
    Q = {int(k): v for k, v in cfg["soc"]["capacity_ah"].items()}
    bounds = (cfg["hppc"]["bounds"]["R0_ohm"][0], cfg["hppc"]["bounds"]["R1_ohm"][0], cfg["hppc"]["bounds"]["tau_s"][0]), \
             (cfg["hppc"]["bounds"]["R0_ohm"][1], cfg["hppc"]["bounds"]["R1_ohm"][1], cfg["hppc"]["bounds"]["tau_s"][1])
    heat = 1.0 / (cfg["temperature"]["mass_kg"] * cfg["temperature"]["cp_J_kgK"])
    rows, moves, axes = [], [], []
    for T, names in load.HPPC.items():
        o = ocv_tab[ocv_tab.Temperature == T]
        f_ocv = PchipInterpolator(o.SOC.to_numpy(), o["OCV-Discharge"].to_numpy(), extrapolate=True)
        for name in names:
            df, steps = load.load_cached(load.hppc_path(name))
            z, ax = soc_axis(df, steps, Q[T], f_ocv)
            axes.append(dict(T=T, cell=load.cell_id(name), **ax))
            st = structure(steps)
            r0_at = {}
            for _, r in st[(st.kind == "pulse") & (st.dir == "dchg")].iterrows():
                d, nh = window(df, steps, r.seg, "pulse")
                t, i, eta, mask, w, seg = window_arrays(d, nh, z[d.index.to_numpy()], f_ocv)
                f = fit_window(t, i, eta, mask, bounds, w, nh, seg)
                Tl, Ts, rise = window_temperature(d, nh, T, t)
                rows.append(dict(T=T, cell=load.cell_id(name), file=name, seg=int(r.seg), kind="pulse",
                                 point=int(r.point), amp_A=r.amp_A, soc=float(z[d.index[nh]]), R0=f["R0"],
                                 R1=f["R1"], tau=f["tau"], rmse_mV=f["rmse_mV"], T_label=Tl, T_source=Ts,
                                 T_rise_K=rise, heat_bound_K=heat * (f["R0"] + f["R1"]) * float(np.sum(i[nh:-1] ** 2 * np.diff(t[nh:])))))
                r0_at.setdefault(int(r.point), []).append(f["R0"])
            for _, r in st[(st.kind == "move") & (st.point > 1)].iterrows():
                d, nh = window(df, steps, r.seg, "move")
                t, i, eta, mask, w, seg = window_arrays(d, nh, z[d.index.to_numpy()], f_ocv)
                Tl, Ts, rise = window_temperature(d, nh, T, t)
                moves.append(dict(t=t, i=i, eta=eta, mask=mask, w=w, nh=nh, step=seg,
                                  R0=float(np.median(r0_at[int(r.point)])),
                                  info=dict(T=T, cell=load.cell_id(name), file=name, seg=int(r.seg), kind="move",
                                            point=int(r.point), amp_A=float(r.I), soc=float(z[d.index[-1]]),
                                            T_label=Tl, T_source=Ts, T_rise_K=rise,
                                            left_out=any(T == x["T"] and abs(float(z[d.index[-1]]) - x["soc"]) < 0.01
                                                         for x in leave_out),
                                            i2dt=float(np.sum(i[nh:-1] ** 2 * np.diff(t[nh:]))))))
    # one tau1 for all temperatures: minimise the summed weighted error of every move window
    total = lambda lt: sum(solve_given_tau(W, np.exp(lt))[0] for W in moves if not W["info"]["left_out"])
    grid = np.log(np.geomspace(5, 5000, 40))
    k = int(np.argmin([total(q) for q in grid]))
    tau = float(np.exp(minimize_scalar(total, bounds=(grid[max(k - 1, 0)], grid[min(k + 1, len(grid) - 1)]),
                                       method="bounded").x))
    for W in moves:
        _, R1, c = solve_given_tau(W, tau)
        m = W["mask"].copy(); m[:W["nh"]] = False
        e = -W["R0"] * W["i"] - R1 * rc_current(W["t"], W["i"], tau, 0.0, W["step"])
        res = ((e - e[W["nh"]] + c) - W["eta"])[m]
        on = (np.abs(W["i"]) > 0.5)[m]
        tt = W["t"][m]
        t_off = W["t"][np.flatnonzero(np.abs(W["i"]) > 0.5)[-1]]
        info = dict(W["info"])
        info.update(R0=W["R0"], R1=R1, tau=tau, rmse_mV=1e3 * np.sqrt(np.mean(res ** 2)),
                    resid_during_move_mV=1e3 * res[on].mean(),
                    resid_first_10s_of_rest_mV=1e3 * res[(~on) & (tt > t_off) & (tt <= t_off + 10)].mean(),
                    heat_bound_K=heat * (W["R0"] + R1) * info.pop("i2dt"))
        rows.append(info)
    return pd.DataFrame(rows), tau, pd.DataFrame(axes)


# ---------------------------------------------------------------- tables
def _interp_extrap(x, y, xg):
    """Linear inside the data; linear extrapolation from the two nearest points outside it
    (straight line in ln(y) if a linear one would go <= 0; returns whether that happened)."""
    o = np.argsort(x)
    x, y = np.asarray(x)[o], np.asarray(y)[o]
    out, fell_back = np.interp(xg, x, y), False
    for sel, (i0, i1) in [(xg < x[0], (0, 1)), (xg > x[-1], (-1, -2))]:
        if sel.any():
            lin = y[i0] + (xg[sel] - x[i0]) * (y[i1] - y[i0]) / (x[i1] - x[i0])
            if (lin <= 0).any():
                ly = np.log(y)
                lin, fell_back = np.exp(ly[i0] + (xg[sel] - x[i0]) * (ly[i1] - ly[i0]) / (x[i1] - x[i0])), True
            out[sel] = lin
    return out, fell_back


def arrhenius_fit(T_c, y):
    """ln y = a + b / T_K per column; activation energy = b * kB [eV]."""
    X = np.c_[np.ones(len(T_c)), 1.0 / (np.asarray(T_c, float) + 273.15)]
    coef, *_ = np.linalg.lstsq(X, np.log(np.asarray(y, float)), rcond=None)
    return coef[0], coef[1]


def arrhenius_eval(a, b, T_c):
    return np.exp(a[None, :] + b[None, :] / (np.asarray(T_c, float)[:, None] + 273.15))


def tables(fits):
    """(full table with every T row, per-cell table, Arrhenius activation energies)."""
    f = fits[~excluded(fits)]
    p = f[f.kind == "pulse"].groupby(["T", "cell", "point"]).agg(R0=("R0", "median"), soc=("soc", "median")).reset_index()
    mv = f[f.kind == "move"]
    rows, fallback = [], []
    for (T, cell), g in mv.groupby(["T", "cell"]):
        gp = p[(p["T"] == T) & (p.cell == cell)]
        R0, fb0 = _interp_extrap(gp.soc, gp.R0, SOC_GRID)
        R1, fb1 = _interp_extrap(g.soc, g.R1, SOC_GRID)
        tau, _ = _interp_extrap(g.soc, g.tau, SOC_GRID)
        fallback += [f"{T}C {cell} R0"] * fb0 + [f"{T}C {cell} R1"] * fb1
        rows += [dict(T=T, cell=cell, SOC=s, R0=R0[k], R1=R1[k], tau=tau[k], T_label=g.T_label.mean())
                 for k, s in enumerate(SOC_GRID)]
    cg = pd.DataFrame(rows)
    cg["C1"] = cg.tau / cg.R1
    g = cg.groupby(["T", "SOC"])
    cons = g[["R0", "R1", "tau", "T_label"]].median().join(g[["R0", "R1", "C1"]].std(ddof=1).add_prefix("err_")).reset_index()
    cons["C1"] = cons.tau / cons.R1
    cons["source"] = "measured"
    filled, ea = [], []
    for s, gg in cons.groupby("SOC"):
        a, b = arrhenius_fit(gg.T_label.to_numpy(), gg[["R0", "R1", "tau"]].to_numpy())
        ea.append(dict(SOC=s, Ea_R0_eV=b[0] * KB_EV, Ea_R1_eV=b[1] * KB_EV))
        rel = np.nanmedian(gg[["err_R0", "err_R1", "err_C1"]].to_numpy() / gg[["R0", "R1", "C1"]].to_numpy(), axis=0)
        fill_T = [T for T in load.R_T_GRID if T not in MEASURED]
        for T, (R0, R1, tau) in zip(fill_T, arrhenius_eval(a, b, fill_T)):
            filled.append(dict(T=T, SOC=s, R0=R0, R1=R1, tau=tau, C1=tau / R1, err_R0=R0 * rel[0],
                               err_R1=R1 * rel[1], err_C1=tau / R1 * rel[2], T_label=float(T), source="Arrhenius"))
    full = pd.concat([cons, pd.DataFrame(filled)], ignore_index=True).sort_values(["T", "SOC"]).reset_index(drop=True)
    return full, cg, pd.DataFrame(ea), fallback


def to_datapack(full):
    r0 = pd.DataFrame({"Temperature": full["T"], "SOC": full.SOC, "Mean (mOhms)": 1e3 * full.R0, "Error (mOhms)": 1e3 * full.err_R0})
    r1 = pd.DataFrame({"Temperature": full["T"], "SOC": full.SOC, "Mean (mOhms)": 1e3 * full.R1, "Error (mOhms)": 1e3 * full.err_R1})
    c1 = pd.DataFrame({"Temperature": full["T"], "SOC": full.SOC, "Mean (Farad)": full.C1, "Error (Farad)": full.err_C1})
    return r0, r1, c1


def from_datapack(r0, r1, c1):
    out = pd.DataFrame({"T": r0.Temperature, "SOC": r0.SOC, "R0": r0.iloc[:, 2] * 1e-3,
                        "R1": r1.iloc[:, 2] * 1e-3, "C1": c1.iloc[:, 2]})
    out["tau"] = out.R1 * out.C1
    return out
