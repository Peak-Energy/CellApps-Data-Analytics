"""RC models from HPPC, gridded onto the team SOC/temperature grid with an Arrhenius fill for the cold rows.

1RC: R0 from the pulses, one tau1 for every SOC and temperature and R1 per SOC from the SOC moves.
2RC: the same R0; one fast and one slow time constant for every SOC and temperature, found together on
     all windows (discharge pulses, charge pulses and SOC moves); the two resistances fitted jointly per
     SOC point, the fast one kept where the point has pulses and the slow one where it has a move.

Model (current discharge-positive, exact for piecewise-constant current), one term per RC branch:
    V_k = OCV(SOC_k, T) - R0 i_k - sum_b R_b x_{b,k},   x_{b,k+1} = a_k x_{b,k} + (1 - a_k) i_k,   a_k = exp(-dt_k / tau_b)
The current on row k is held until row k+1. The last row of one step and the first row of the next
are the same instant (dt = 0).

HPPC schedule at each SOC point: 2 h rest after the SOC move, then for each amplitude a discharge
pulse, 60 s rest, charge pulse, 900 s rest (2 h after the last charge pulse); the last point runs
charge-first. Fit windows: a pulse window is the rest sample before the pulse + the 10 s pulse + its
rest (capped at 900 s); a move window is the rest sample before the SOC move + the move + its 2 h rest.
The first row of every step (mid-switch) and a second row that repeats its voltage are never fitted.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.optimize import least_squares, minimize_scalar, nnls
from pathlib import Path

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


def rc_currents(t, i, taus, step):
    """rc_current from x0 = 0 with step-mean current for several time constants at once: n x K."""
    t = np.asarray(t, float); i = np.asarray(i, float); taus = np.asarray(taus, float)
    step = np.asarray(step)
    starts = np.flatnonzero(np.r_[True, step[1:] != step[:-1]])
    x = np.empty((len(t), len(taus)))
    xs = np.zeros(len(taus))
    for s0, s1 in zip(starts, np.r_[starts[1:], len(t)]):
        I = i[s0:s1].mean()
        x[s0:s1] = I + (xs - I) * np.exp(-(t[s0:s1, None] - t[s0]) / taus)
        xs = x[s1 - 1]
    x[0] = 0.0
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


def unusable_rows(seg, V):
    """Rows the logger did not capture properly: the first row of every step (logged mid-switch) and a
    second row that repeats the first row's voltage exactly (the logger had not refreshed the voltage; it
    happens whenever the second row falls within ~0.5 s on a step logged coarser than 0.1 s)."""
    seg, V = np.asarray(seg), np.asarray(V)
    first = np.r_[False, seg[1:] != seg[:-1]]
    stale = np.r_[False, first[:-1]] & np.r_[False, V[1:] == V[:-1]]
    return first | stale


def window_arrays(d, nh, z, ocv_fn):
    """Time from window start, current, overpotential, mask (see unusable_rows), sqrt(dt) weights."""
    t = d.t.to_numpy() - d.t.iloc[nh]
    V, ov = d.voltage.to_numpy(), ocv_fn(z)
    eta = (V - ov) - (V[nh] - ov[nh])
    seg = d.seg.to_numpy()
    mask = ~unusable_rows(seg, V)
    w = np.sqrt(np.gradient(t))
    w[~np.isfinite(w) | (w == 0)] = np.sqrt(0.1)
    return t, d.current.to_numpy(), eta, mask, w, seg


def window_temperature(d, nh, T, t):
    w = d.iloc[nh:]
    if "aux_T1a" not in w:
        return float(T), "setpoint", np.nan
    Tc = w[["aux_T1a", "aux_T2a"]].mean(axis=1).to_numpy()
    return float(Tc.mean()), "logged", ocv.temp_rise(t[nh:], Tc, 60.0 if t[-1] <= 900 else 600.0)


def is_excluded(kind, point, last_point, left_out=False):
    """Windows that start in a state no other SOC point has: the 100 % point's pulses and the first SOC
    move (fresh off the charge to 3.45 V), and the last point's pulses (charge-first order); plus moves
    left out by name (see prepare)."""
    return (kind == "pulse" and point in (0, last_point)) or (kind == "move" and point == 1) or bool(left_out)


def excluded(fits):
    last = fits.groupby("file").point.transform("max")
    named = fits["left_out"].fillna(False).astype(bool) if "left_out" in fits else False
    return (((fits.kind == "pulse") & ((fits.point == 0) | (fits.point == last)))
            | ((fits.kind == "move") & (fits.point == 1)) | named)


KEEP_COLS = ["t", "current", "voltage", "seg", "aux_T1a", "aux_T2a"]


def prepare(cfg, ocv_tab, leave_out=(), files=None):
    """Load every HPPC file, build its SOC axis and fit the discharge pulses one by one (R0; the free 1RC of
    each pulse is kept for the record). Returns (file records for `windows`, per-pulse rows, SOC-axis rows).
    files: {T: [path, ...]}, default the -B HPPC files.
    leave_out: [{T, soc}] moves ending within 1 % of that SOC are reported but kept out of every fit and table."""
    files = files or {T: [load.hppc_path(n) for n in names] for T, names in load.HPPC.items()}
    Q = {int(k): v for k, v in cfg["soc"]["capacity_ah"].items()}
    bounds = (cfg["hppc"]["bounds"]["R0_ohm"][0], cfg["hppc"]["bounds"]["R1_ohm"][0], cfg["hppc"]["bounds"]["tau_s"][0]), \
             (cfg["hppc"]["bounds"]["R0_ohm"][1], cfg["hppc"]["bounds"]["R1_ohm"][1], cfg["hppc"]["bounds"]["tau_s"][1])
    heat = 1.0 / (cfg["temperature"]["mass_kg"] * cfg["temperature"]["cp_J_kgK"])
    recs, rows, axes = [], [], []
    for T, paths in files.items():
        o = ocv_tab[ocv_tab.Temperature == T]
        f_ocv = PchipInterpolator(o.SOC.to_numpy(), o["OCV-Discharge"].to_numpy(), extrapolate=True)
        for path in paths:
            path = Path(path)
            df, steps = load.load_cached(path)
            name, cell = path.name, load.cell_id(path.name)
            z, ax = soc_axis(df, steps, Q[T], f_ocv)
            axes.append(dict(T=T, cell=cell, **ax))
            st = structure(steps)
            last = int(st.point.max())
            r0_at = {}
            for _, r in st[(st.kind == "pulse") & (st.dir == "dchg")].iterrows():
                d, nh = window(df, steps, r.seg, "pulse")
                t, i, eta, mask, w, seg = window_arrays(d, nh, z[d.index.to_numpy()], f_ocv)
                f = fit_window(t, i, eta, mask, bounds, w, nh, seg)
                Tl, Ts, rise = window_temperature(d, nh, T, t)
                rows.append(dict(T=T, cell=cell, file=name, seg=int(r.seg), kind="pulse", dir="dchg",
                                 point=int(r.point), amp_A=r.amp_A, soc=float(z[d.index[nh]]), R0=f["R0"],
                                 R1=f["R1"], tau1=f["tau"], rmse_mV=f["rmse_mV"], T_label=Tl, T_source=Ts,
                                 T_rise_K=rise, left_out=False, excluded=is_excluded("pulse", int(r.point), last),
                                 heat_bound_K=heat * (f["R0"] + f["R1"]) * float(np.sum(i[nh:-1] ** 2 * np.diff(t[nh:])))))
                r0_at.setdefault(int(r.point), []).append(f["R0"])
            recs.append(dict(T=T, cell=cell, name=name, df=df[[c for c in KEEP_COLS if c in df]], steps=steps, z=z,
                             f_ocv=f_ocv, st=st, last=last, heat=heat, leave_out=list(leave_out),
                             r0_point={p: float(np.median(v)) for p, v in r0_at.items()}))
    return recs, pd.DataFrame(rows), pd.DataFrame(axes)


def windows(rec, kinds=("pulse", "move"), dirs=("dchg", "chg")):
    """The fit windows of one prepared file, with the SOC point's R0 attached (median of its discharge
    pulses). The first SOC move is never a window (it starts fresh off the charge)."""
    df, steps, z, f_ocv = rec["df"], rec["steps"], rec["z"], rec["f_ocv"]
    for _, r in rec["st"].iterrows():
        if r.kind not in kinds or (r.kind == "pulse" and r.dir not in dirs) or (r.kind == "move" and r.point <= 1):
            continue
        d, nh = window(df, steps, r.seg, r.kind)
        t, i, eta, mask, w, seg = window_arrays(d, nh, z[d.index.to_numpy()], f_ocv)
        Tl, Ts, rise = window_temperature(d, nh, rec["T"], t)
        soc = float(z[d.index[nh]]) if r.kind == "pulse" else float(z[d.index[-1]])
        left = r.kind == "move" and any(rec["T"] == x["T"] and abs(soc - x["soc"]) < 0.01 for x in rec["leave_out"])
        info = dict(T=rec["T"], cell=rec["cell"], file=rec["name"], seg=int(r.seg), kind=r.kind, dir=r.dir,
                    point=int(r.point), amp_A=float(r.amp_A if r.kind == "pulse" else r.I), soc=soc, T_label=Tl,
                    T_source=Ts, T_rise_K=rise, left_out=left,
                    excluded=is_excluded(r.kind, int(r.point), rec["last"], left),
                    i2dt=float(np.sum(i[nh:-1] ** 2 * np.diff(t[nh:]))))
        yield dict(t=t, i=i, eta=eta, mask=mask, w=w, nh=nh, step=seg, R0=rec["r0_point"][int(r.point)], info=info)


# ---------------------------------------------------------------- 1RC fit
def fit_1rc(prep):
    """Rows for every discharge pulse (R0) and SOC move (R1 at one tau1 for all temperatures), and tau1.
    tau1 minimises the summed weighted error of every move window that is not left out."""
    recs, pulse_rows, _ = prep
    moves = [W for rec in recs for W in windows(rec, kinds=("move",))]
    total = lambda lt: sum(solve_given_tau(W, np.exp(lt))[0] for W in moves if not W["info"]["left_out"])
    grid = np.log(np.geomspace(5, 5000, 40))
    k = int(np.argmin([total(q) for q in grid]))
    tau = float(np.exp(minimize_scalar(total, bounds=(grid[max(k - 1, 0)], grid[min(k + 1, len(grid) - 1)]),
                                       method="bounded").x))
    rows = []
    for W in moves:
        _, R1, c = solve_given_tau(W, tau)
        m = W["mask"].copy(); m[:W["nh"]] = False
        e = -W["R0"] * W["i"] - R1 * rc_current(W["t"], W["i"], tau, 0.0, W["step"])
        res = ((e - e[W["nh"]] + c) - W["eta"])[m]
        on = (np.abs(W["i"]) > 0.5)[m]
        tt = W["t"][m]
        t_off = W["t"][np.flatnonzero(np.abs(W["i"]) > 0.5)[-1]]
        info = dict(W["info"])
        info.update(R0=W["R0"], R1=R1, tau1=tau, rmse_mV=1e3 * np.sqrt(np.mean(res ** 2)),
                    resid_during_move_mV=1e3 * res[on].mean(),
                    resid_first_10s_of_rest_mV=1e3 * res[(~on) & (tt > t_off) & (tt <= t_off + 10)].mean(),
                    heat_bound_K=recs[0]["heat"] * (W["R0"] + R1) * info.pop("i2dt"))
        rows.append(info)
    return pd.concat([pulse_rows, pd.DataFrame(rows)], ignore_index=True), tau


# ---------------------------------------------------------------- 2RC fit
def summarize(W, taus):
    """What the joint fit needs from one window at the given time constants, with the window's zero level
    projected out: the weighted normal equations (G, b, yy) and, for the residual statistics, the unweighted
    ones (Gr, br, yyr) plus sums over the rows during the current step and the first 10 s of the rest."""
    m = W["mask"].copy(); m[:W["nh"]] = False
    w = W["w"][m]
    X = rc_currents(W["t"], W["i"], taus, W["step"])[m]
    y = -(W["eta"] + W["R0"] * (W["i"] - W["i"][W["nh"]]))[m]          # eta = -R0 i - sum R x + c
    w2 = w ** 2
    Xc, yc = X - (w2 @ X) / w2.sum(), y - (w2 @ y) / w2.sum()
    Xw, yw = Xc * w[:, None], yc * w
    on = (np.abs(W["i"]) > 0.5)[m]
    tm, t_off = W["t"][m], W["t"][np.flatnonzero(np.abs(W["i"]) > 0.5)[-1]]
    first10 = (~on) & (tm > t_off) & (tm <= t_off + 10)
    return dict(G=Xw.T @ Xw, b=Xw.T @ yw, yy=float(yw @ yw), Gr=Xc.T @ Xc, br=Xc.T @ yc, yyr=float(yc @ yc),
                n=int(m.sum()), sel={k: (Xc[s].sum(0), float(yc[s].sum()), int(s.sum()))
                                     for k, s in [("during", on), ("first10", first10)]})


def nnls2(G11, G12, G22, b1, b2, yy):
    """Least squares in (R1, R2) >= 0 from 2x2 normal equations, elementwise over arrays: (R1, R2, SSE)."""
    G11, G12, G22, b1, b2, yy = np.broadcast_arrays(*(np.asarray(a, float) for a in (G11, G12, G22, b1, b2, yy)))
    sse = lambda a, c: yy - 2 * (a * b1 + c * b2) + a * a * G11 + 2 * a * c * G12 + c * c * G22
    with np.errstate(divide="ignore", invalid="ignore"):
        det = G11 * G22 - G12 ** 2
        r1, r2 = (G22 * b1 - G12 * b2) / det, (G11 * b2 - G12 * b1) / det
        c1 = np.where(G11 > 0, np.clip(b1 / np.where(G11 > 0, G11, 1), 0, None), 0.0)
        c2 = np.where(G22 > 0, np.clip(b2 / np.where(G22 > 0, G22, 1), 0, None), 0.0)
    ok = (det > 0) & (r1 >= 0) & (r2 >= 0)
    one = sse(c1, 0.0) <= sse(0.0, c2)
    R1 = np.where(ok, r1, np.where(one, c1, 0.0))
    R2 = np.where(ok, r2, np.where(one, 0.0, c2))
    return R1, R2, sse(R1, R2)


def _pair_sse(P):
    """SSE of the best non-negative (R_a, R_b) for every pair of grid time constants: K x K."""
    d = np.diag(P["G"])
    return nnls2(d[:, None], P["G"], d[None, :], P["b"][:, None], P["b"][None, :], P["yy"])[2]


def _accumulate(pts, key, S, kind):
    P = pts.setdefault(key, dict(G=0.0, b=0.0, yy=0.0, n=0, kinds=set()))
    P["G"] = P["G"] + S["G"]; P["b"] = P["b"] + S["b"]; P["yy"] += S["yy"]; P["n"] += S["n"]; P["kinds"].add(kind)


def _refine(lt, sse, kf, ks):
    """Quadratic refinement of the grid minimum (kf, ks) in log tau from its 3 x 3 neighbourhood; falls back
    to the grid point when the fitted minimum leaves the cell."""
    u, v = np.meshgrid([-1, 0, 1], [-1, 0, 1], indexing="ij")
    u, v = u.ravel(), v.ravel()
    f = sse[kf - 1:kf + 2, ks - 1:ks + 2].ravel()
    A = np.c_[np.ones(9), u, v, u * u, u * v, v * v]
    c = np.linalg.lstsq(A, f, rcond=None)[0]
    H = np.array([[2 * c[3], c[4]], [c[4], 2 * c[5]]])
    if np.all(np.linalg.eigvalsh(H) > 0):
        du, dv = np.linalg.solve(H, -c[1:3])
        if abs(du) <= 1 and abs(dv) <= 1:
            h = lt[1] - lt[0]
            return float(np.exp(lt[kf] + du * h)), float(np.exp(lt[ks] + dv * h))
    return float(np.exp(lt[kf])), float(np.exp(lt[ks]))


def _best_pair(sse, taus):
    upper = np.triu(np.ones_like(sse, bool), 1)
    kf, ks = np.unravel_index(np.argmin(np.where(upper, sse, np.inf)), sse.shape)
    return int(kf), int(ks)


def fit_2rc(prep, two_rc, fixed=None):
    """Joint fast + slow RC fit. two_rc: {tau_grid_s: [lo, hi, n]}.

    Pass 1: for every window not excluded, the normal equations for every grid tau, summed per SOC point;
    the summed error of the best (R_fast, R_slow) >= 0 per point for every tau pair picks the pair, refined
    by a quadratic through its neighbours. Only points with both pulses and a move count.
    Pass 2: at the chosen pair, (R_fast, R_slow) per point. A point that has no usable move (point 1) or no
    usable pulses (the last point) holds the branch it cannot see at the nearest complete point's value.
    fixed = (tau_fast, tau_slow) skips pass 1 and fits at that pair.
    Returns (rows per window, {fast, slow} time constants, diagnostics)."""
    recs, pulse_rows, _ = prep
    lo, hi, n = two_rc["tau_grid_s"]
    taus = np.geomspace(lo, hi, int(n)); lt = np.log(taus)
    diag = dict(tau_grid_s=[lo, hi, int(n)])
    if fixed is None:
        pts, by_kind = {}, {"pulse": {}, "move": {}}
        for rec in recs:
            for W in windows(rec):
                if not W["info"]["excluded"]:
                    S, key, kind = summarize(W, taus), (rec["name"], W["info"]["point"]), W["info"]["kind"]
                    _accumulate(pts, key, S, kind)
                    _accumulate(by_kind[kind], key, S, kind)
        complete = {k: P for k, P in pts.items() if P["kinds"] == {"pulse", "move"}}
        sse = sum(_pair_sse(P) for P in complete.values())
        n_rows = sum(P["n"] for P in complete.values())
        kf, ks = _best_pair(sse, taus)
        at_edge = kf == 0 or ks == len(taus) - 1 or ks == kf + 1
        tf, ts = (float(taus[kf]), float(taus[ks])) if at_edge else _refine(lt, sse, kf, ks)
        parts = {}
        for kind, P_k in by_kind.items():
            s_k = sum(_pair_sse(P_k[k]) for k in complete)
            a, b = _best_pair(s_k, taus)
            parts[kind] = dict(best_pair_s=dict(fast=float(taus[a]), slow=float(taus[b])),
                               best_fast_s_at_chosen_slow=float(taus[int(np.argmin(s_k[:ks, ks]))]),
                               rms_mV_at_chosen_pair=1e3 * np.sqrt(s_k[kf, ks] / n_rows))
        d = np.diag(sse)
        upper = np.triu(np.ones_like(sse, bool), 1)
        diag.update(grid_best=dict(fast=float(taus[kf]), slow=float(taus[ks])), refined=dict(fast=tf, slow=ts),
                    at_grid_edge=bool(at_edge), by_window_kind=parts,
                    rms_mV_best_single_tau=1e3 * np.sqrt(d.min() / n_rows), best_single_tau_s=float(taus[int(np.argmin(d))]),
                    sse_surface_mV=[[None if not np.isfinite(v) else round(float(v), 4) for v in row]
                                    for row in 1e3 * np.sqrt(np.where(upper, sse, np.nan) / n_rows)])
        t4 = [tf, ts, float(taus[kf]), float(taus[ks])]
    else:
        tf, ts = map(float, fixed)
        t4 = [tf, ts, tf, ts]
        diag.update(fixed=dict(fast=tf, slow=ts))
    # pass 2: exact solve at the refined and the grid pair, keep the better
    pts, wins = {}, []
    for rec in recs:
        for W in windows(rec):
            S = summarize(W, t4); S["info"] = W["info"]; S["R0"] = W["R0"]
            wins.append(S)
            if not W["info"]["excluded"]:
                _accumulate(pts, (rec["name"], W["info"]["point"]), S, W["info"]["kind"])
    complete = {k: P for k, P in pts.items() if P["kinds"] == {"pulse", "move"}}
    tot = {}
    for a, b in [(0, 1), (2, 3)]:
        tot[(a, b)] = sum(float(nnls2(P["G"][a, a], P["G"][a, b], P["G"][b, b], P["b"][a], P["b"][b], P["yy"])[2])
                          for P in complete.values())
    a, b = min(tot, key=tot.get)
    tf, ts = t4[a], t4[b]
    sol = {}
    for k, P in complete.items():
        Rf, Rs, _ = nnls2(P["G"][a, a], P["G"][a, b], P["G"][b, b], P["b"][a], P["b"][b], P["yy"])
        sol[k] = dict(R_fast=float(Rf), R_slow=float(Rs), held=None)
    for k, P in pts.items():
        if k in sol:
            continue
        near = min((kk for kk in complete if kk[0] == k[0]), key=lambda kk: abs(kk[1] - k[1]), default=None)
        if near is None:
            continue
        if "move" not in P["kinds"]:                       # hold the slow branch, fit the fast one
            Rs = sol[near]["R_slow"]
            Rf = max((P["b"][a] - Rs * P["G"][a, b]) / P["G"][a, a], 0.0)
            sol[k] = dict(R_fast=float(Rf), R_slow=Rs, held="slow")
        else:                                              # hold the fast branch, fit the slow one
            Rf = sol[near]["R_fast"]
            Rs = max((P["b"][b] - Rf * P["G"][a, b]) / P["G"][b, b], 0.0)
            sol[k] = dict(R_fast=Rf, R_slow=float(Rs), held="fast")
    rows = []
    for S in wins:
        info = dict(S["info"])
        key = (info["file"], info["point"])
        s = sol.get(key, dict(R_fast=np.nan, R_slow=np.nan, held="none"))
        R = np.array([s["R_fast"], s["R_slow"]])
        idx = [a, b]
        Gr, br = S["Gr"][np.ix_(idx, idx)], S["br"][idx]
        rss = S["yyr"] - 2 * R @ br + R @ Gr @ R if np.all(np.isfinite(R)) else np.nan
        mean = lambda sel: 1e3 * (sel[1] - R @ sel[0][idx]) / sel[2] if sel[2] and np.all(np.isfinite(R)) else np.nan
        own = max((S["b"][b] - s["R_fast"] * S["G"][a, b]) / S["G"][b, b], 0.0) if np.isfinite(s["R_fast"]) else np.nan
        info.update(R0=S["R0"], R_fast=s["R_fast"], R_slow=s["R_slow"], tau_fast=tf, tau_slow=ts, held=s["held"],
                    rmse_mV=1e3 * np.sqrt(rss / S["n"]), resid_during_mV=mean(S["sel"]["during"]),
                    resid_first_10s_of_rest_mV=mean(S["sel"]["first10"]),
                    R_slow_this_window_only=own if info["kind"] == "move" else np.nan,
                    heat_bound_K=recs[0]["heat"] * (S["R0"] + np.nansum(R)) * info.pop("i2dt"))
        rows.append(info)
    n_rows = sum(P["n"] for P in complete.values())
    diag.update(used=("refined" if (a, b) == (0, 1) else "grid"), n_points_complete=len(complete),
                n_points_held=sum(v["held"] is not None for v in sol.values()),
                rms_mV_best_pair=1e3 * np.sqrt(tot[(a, b)] / n_rows))
    return pd.DataFrame(rows), dict(fast=tf, slow=ts), diag


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


BRANCHES_1RC = [("1", "R1", "tau1", "move")]
BRANCHES_2RC = {"slow": [("1", "R_slow", "tau_slow", "move"), ("2", "R_fast", "tau_fast", "pulse")],
                "fast": [("1", "R_fast", "tau_fast", "pulse"), ("2", "R_slow", "tau_slow", "move")]}


def tables(fits, branches=BRANCHES_1RC):
    """(full table with every T row, per-cell table, Arrhenius activation energies, notes on points that were
    interpolated or log-extrapolated). A point whose resistance fitted to exactly 0 (the non-negative solve's
    boundary) is interpolated from its neighbours.
    branches: [(number, R column, tau column, window kind that defines it)] -> output columns R<n>, tau<n>, C<n>."""
    f = fits[~excluded(fits)]
    src = {"R0": ("R0", "pulse", "dchg")} | {f"R{n}": (rc, kind, None) for n, rc, _, kind in branches}
    rows, fallback = [], []
    for (T, cell), g in f.groupby(["T", "cell"]):
        vals = {}
        for q, (col, kind, d) in src.items():
            gk = g[(g.kind == kind) & ((g.dir == d) if d else True)]
            gq = gk.groupby("point").agg(soc=("soc", "median"), v=(col, "median")).dropna()
            zero = gq[gq.v <= 0]
            fallback += [f"{T}C {cell} {q} = 0 at SOC {z:.3f}: interpolated from the neighbours" for z in zero.soc]
            gq = gq[gq.v > 0]
            vals[q], fb = _interp_extrap(gq.soc, gq.v, SOC_GRID)
            fallback += [f"{T}C {cell} {q}: log extrapolation at a SOC end"] * fb
        for n, _, tc, kind in branches:
            vals[f"tau{n}"] = np.full(len(SOC_GRID), float(g[g.kind == kind][tc].median()))
        T_label = g[g.kind == "move"].T_label.mean()                 # logged cell temperature over the move windows
        rows += [dict(T=T, cell=cell, SOC=s, T_label=T_label, **{q: v[k] for q, v in vals.items()})
                 for k, s in enumerate(SOC_GRID)]
    cg = pd.DataFrame(rows)
    nums = [n for n, *_ in branches]
    for n in nums:
        cg[f"C{n}"] = cg[f"tau{n}"] / cg[f"R{n}"]
    R_cols = ["R0"] + [f"R{n}" for n in nums]
    tau_cols = [f"tau{n}" for n in nums]
    C_cols = [f"C{n}" for n in nums]
    g = cg.groupby(["T", "SOC"])
    cons = g[R_cols + tau_cols + ["T_label"]].median().join(g[R_cols + C_cols].std(ddof=1).add_prefix("err_")).reset_index()
    for n in nums:
        cons[f"C{n}"] = cons[f"tau{n}"] / cons[f"R{n}"]
    cons["source"] = "measured"
    filled, ea = [], []
    fill_T = [T for T in load.R_T_GRID if T not in MEASURED]
    for s, gg in cons.groupby("SOC"):
        a, b = arrhenius_fit(gg.T_label.to_numpy(), gg[R_cols + tau_cols].to_numpy())
        ea.append(dict(SOC=s, **{f"Ea_{c}_eV": b[k] * KB_EV for k, c in enumerate(R_cols)}))
        rel = np.nanmedian(gg[[f"err_{c}" for c in R_cols + C_cols]].to_numpy() / gg[R_cols + C_cols].to_numpy(), axis=0)
        for T, v in zip(fill_T, arrhenius_eval(a, b, fill_T)):
            row = dict(T=T, SOC=s, T_label=float(T), source="Arrhenius", **dict(zip(R_cols + tau_cols, v)))
            for n in nums:
                row[f"C{n}"] = row[f"tau{n}"] / row[f"R{n}"]
            for k, c in enumerate(R_cols + C_cols):
                row[f"err_{c}"] = row[c] * rel[k]
            filled.append(row)
    full = pd.concat([cons, pd.DataFrame(filled)], ignore_index=True).sort_values(["T", "SOC"]).reset_index(drop=True)
    return full, cg, pd.DataFrame(ea), fallback


def to_datapack(full):
    """{'r0': table, 'r1': table, 'c1': table[, 'r2', 'c2']} in the team file format."""
    out = {"r0": pd.DataFrame({"Temperature": full["T"], "SOC": full.SOC, "Mean (mOhms)": 1e3 * full.R0,
                               "Error (mOhms)": 1e3 * full.err_R0})}
    for n in [c[1:] for c in full.columns if c.startswith("R") and c != "R0"]:
        out[f"r{n}"] = pd.DataFrame({"Temperature": full["T"], "SOC": full.SOC, "Mean (mOhms)": 1e3 * full[f"R{n}"],
                                     "Error (mOhms)": 1e3 * full[f"err_R{n}"]})
        out[f"c{n}"] = pd.DataFrame({"Temperature": full["T"], "SOC": full.SOC, "Mean (Farad)": full[f"C{n}"],
                                     "Error (Farad)": full[f"err_C{n}"]})
    return out


def from_datapack(r0, *rc):
    """Table of R0 and (R<n>, tau<n>) per T and SOC from datapack tables r0, r1, c1[, r2, c2, ...]."""
    out = pd.DataFrame({"T": r0.Temperature, "SOC": r0.SOC, "R0": r0.iloc[:, 2] * 1e-3})
    for n, (r, c) in enumerate(zip(rc[0::2], rc[1::2]), start=1):
        out[f"R{n}"] = r.iloc[:, 2].to_numpy() * 1e-3
        out[f"tau{n}"] = out[f"R{n}"] * c.iloc[:, 2].to_numpy()
    return out


def branches_of(tab):
    """[(R column, tau column)] of a parameter table."""
    return [(f"R{n}", f"tau{n}") for n in range(1, 10) if f"R{n}" in tab]
