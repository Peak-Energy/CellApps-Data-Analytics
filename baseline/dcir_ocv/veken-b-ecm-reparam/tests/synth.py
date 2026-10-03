"""Synthetic Neware-format HPPC record with known OCV and 1RC or 2RC parameters (verification only)."""
from __future__ import annotations

import datetime as dt
from pathlib import Path

import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from reparam import ecm, load  # noqa: E402

HEADER = load.CORE_HEADER


def _hms(s: float) -> str:
    s = int(np.floor(s + 1e-9))
    return f"{s // 3600:02d}:{(s % 3600) // 60:02d}:{s % 60:02d}"


def schedule(amps=(56.0, 85.0, 170.0), n_points=4, move_I=56.0, move_s=557.0):
    """List of (step_type, current discharge-positive, duration_s, logging_s)."""
    steps = [("Rest", 0.0, 7200.0, 300.0)]
    for p in range(n_points):
        for a in amps:
            steps += [("CC DChg", a, 10.0, 0.1), ("Rest", 0.0, 60.0, 0.1), ("CC Chg", -a, 10.0, 0.1)]
            steps += [("Rest", 0.0, 900.0, 0.1)] if a != amps[-1] else [("Rest", 0.0, 7200.0, 300.0)]
        if p < n_points - 1:
            steps += [("CC DChg", move_I, move_s, 10.0), ("Rest", 0.0, 7200.0, 300.0)]
    return steps


def write(path, ocv_fn, R0, R1, tau, Q_ah=172.4, z0=1.0, noise_uV=37.0, transition_frac=0.3,
          seed=0, steps=None, R2=0.0, tau2=1.0):
    """Simulate the schedule exactly at the logged instants and write a Neware CSV.

    A second RC branch (R2, tau2) is added when R2 > 0. Returns a dict of the true values. The first
    row of every current step shows full current but only `transition_frac` of the step's instantaneous
    voltage change (as in the real logs).
    """
    rng = np.random.default_rng(seed)
    steps = steps or schedule()
    rows_t, rows_i, rows_step, rows_trel = [], [], [], []
    t_abs = 0.0
    for k, (_, I, dur, lg) in enumerate(steps):
        n = int(round(dur / lg))
        tr = np.arange(n + 1) * lg
        rows_t.append(t_abs + tr); rows_i.append(np.full(n + 1, I)); rows_step.append(np.full(n + 1, k))
        rows_trel.append(tr)
        t_abs += dur
    t = np.concatenate(rows_t); i = np.concatenate(rows_i)
    step = np.concatenate(rows_step); trel = np.concatenate(rows_trel)
    # coulomb counting (current held over (t_k, t_k+1])
    dq = np.r_[0.0, i[:-1] * np.diff(t) / 3600.0]
    q = np.cumsum(dq)
    z = z0 - q / Q_ah
    v = ocv_fn(z) - R0 * i - R1 * ecm.rc_current(t, i, tau) - R2 * ecm.rc_current(t, i, tau2)
    first = np.r_[True, step[1:] != step[:-1]]
    cur_first = first & (i != 0)
    v[cur_first] = v[np.flatnonzero(cur_first) - 1] + transition_frac * (v[cur_first] - v[np.flatnonzero(cur_first) - 1])
    off_first = first & (i == 0) & (np.arange(len(i)) > 0)
    v[off_first] = v[np.flatnonzero(off_first) - 1] + transition_frac * (v[off_first] - v[np.flatnonzero(off_first) - 1])
    v_noisy = v + rng.normal(0, noise_uV * 1e-6, len(v))
    # per-step capacity counter
    cap = np.zeros_like(t)
    for k in np.unique(step):
        m = step == k
        cap[m] = np.abs(np.r_[0.0, np.cumsum(i[m][:-1] * np.diff(t[m]))]) / 3600.0
    t0 = dt.datetime(2026, 1, 1, 0, 0, 0)
    phase = 0.37                                     # wall clock ticks at a different sub-second phase
    lines = [",".join(HEADER)]
    for r in range(len(t)):
        st = steps[step[r]][0]
        chg = cap[r] if i[r] < 0 else 0.0
        dch = cap[r] if i[r] > 0 else 0.0
        date = (t0 + dt.timedelta(seconds=int(np.floor(t[r] + phase)))).strftime("%Y-%m-%d %H:%M:%S")
        lines.append(",".join([str(r + 1), "1", str(step[r] + 1), st, _hms(trel[r]), _hms(t[r]),
                               f"{-i[r]:.4f}", f"{v_noisy[r]:.6f}", f"{cap[r]:.4f}", f"{chg:.4f}", f"{dch:.4f}",
                               "0.0000", "0.0000", "0.0000", date, f"{-i[r] * v_noisy[r]:.4f}"]))
    Path(path).write_text("﻿" + "\n".join(lines) + "\n", encoding="utf-8")
    return dict(t=t, i=i, v=v, z=z, R0=R0, R1=R1, tau=tau, R2=R2, tau2=tau2, Q_ah=Q_ah, z0=z0)
