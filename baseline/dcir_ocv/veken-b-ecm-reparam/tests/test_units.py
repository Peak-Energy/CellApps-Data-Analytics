"""Unit tests for the loader, OCV and 1RC building blocks, and the datapack file format."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from reparam import ecm, load, ocv  # noqa: E402

AUX = ",V1(V),V2(V),T1(℃),T2(℃),T1(℃),T2(℃),H2(%),H22(%),Aux. ΔV(V),Aux. ΔT(℃)"


def _csv(tmp_path, rows, aux=False):
    p = tmp_path / "CV_0010_CS2D0001_test.csv"
    p.write_text("﻿" + ",".join(load.CORE_HEADER) + (AUX if aux else "") + "\n" + "\n".join(rows) + "\n",
                 encoding="utf-8")
    return p


def _row(n, step, typ, tl, I, V, cap, aux=""):
    return f"{n},1,{step},{typ},{tl},{tl},{I:.4f},{V:.6f},{cap:.4f},0,0,0,0,0,2026-01-01 00:00:00,0" + aux


def test_discharge_current_becomes_positive(tmp_path):
    rows = [_row(1, 1, "Rest", "00:00:00", 0, 3.0, 0), _row(2, 2, "CC DChg", "00:00:00", -50, 2.9, 0),
            _row(3, 2, "CC DChg", "00:00:01", -50, 2.9, 0.0139), _row(4, 3, "CC Chg", "00:00:00", 50, 3.1, 0),
            _row(5, 3, "CC Chg", "00:00:01", 50, 3.1, 0.0139)]
    df, _ = load.load(_csv(tmp_path, rows))
    assert (df[df.kind == "dchg"].current > 0).all() and (df[df.kind == "chg"].current < 0).all()
    assert df.q_ah.iloc[2] == pytest.approx(0.0139) and df.q_ah.iloc[-1] == pytest.approx(0.0)


def test_rows_sharing_a_printed_second_are_spaced_evenly(tmp_path):
    rows = [_row(1, 1, "Rest", "00:00:00", 0, 3.0, 0)]
    for k in range(101):                                  # 10 s pulse logged every 0.1 s, time floored
        rows.append(_row(k + 2, 2, "CC DChg", f"00:00:{int(k * 0.1 + 1e-9):02d}", -56.1, 2.9, 56.1 * k * 0.1 / 3600))
    df, _ = load.load(_csv(tmp_path, rows))
    p = df[df.seg == 1]
    assert np.allclose(p.t_rel.to_numpy(), np.arange(101) * 0.1, atol=1e-9)
    assert np.nanmax(np.abs(p.t_cap_dev)) < 0.004         # capacity counter agrees within its 6 ms resolution


def test_duplicated_temperature_headers_read_by_position(tmp_path):
    aux = ",0.0,0.0,16.1,16.3,16.5,16.5,0,0,0.0,0.2"
    rows = [_row(1, 1, "Rest", "00:00:00", 0, 3.0, 0, aux), _row(2, 1, "Rest", "00:00:01", 0, 3.0, 0, aux)]
    df, _ = load.load(_csv(tmp_path, rows, aux=True))
    assert (df.aux_T1a.iloc[0], df.aux_T2a.iloc[0], df.aux_T1b.iloc[0]) == (16.1, 16.3, 16.5)


def test_rc_update_matches_the_analytic_step_response():
    R0, R1, tau, I = 3e-4, 2e-4, 50.0, 100.0
    t = np.r_[0.0, np.sort(np.random.default_rng(1).uniform(0, 400, 300))]
    i = np.full_like(t, I)
    v = 3.0 - R0 * i - R1 * ecm.rc_current(t, i, tau)
    assert np.max(np.abs(v - (3.0 - R0 * I - R1 * I * (1 - np.exp(-t / tau))))) < 1e-12


def test_step_mean_shortcut_matches_row_by_row_update():
    t = np.r_[np.arange(0, 10.0001, 0.1), 10 + np.arange(0, 60.0001, 0.1)]       # step boundary repeats t = 10
    step = np.r_[np.zeros(101, int), np.ones(601, int)]
    i = np.where(step == 0, 56.0, 0.0)
    assert np.max(np.abs(ecm.rc_current(t, i, 30.0, 2.0) - ecm.rc_current(t, i, 30.0, 2.0, step))) < 1e-9


def test_relaxation_criterion_on_an_exponential_with_known_tau():
    tau, A, theta, W = 2000.0, 0.01, ocv.MV_PER_H, 173.0
    t = np.arange(0, 25200.0, 1.0)
    v = 3.0 + A * np.exp(-t / tau)
    k, V, n = ocv.relaxed_point(t, v, theta, W)
    assert t[k] == pytest.approx(tau * np.log(A / (tau * theta)) + W / 2, abs=5.0)   # trailing fit lags ~W/2
    assert n == len(t) - k                                                         # every later sample qualifies
    assert V == pytest.approx(v[k:].mean(), abs=1e-12)                             # equilibrium = their average
    assert ocv.relaxed_point(t[:1000], v[:1000], theta, W)[0] is None              # cut short -> unrelaxed


def test_ocv_interpolation_is_monotone_without_overshoot():
    z = np.array([0, 0.01, 0.02, 0.04, 0.1, 0.3, 0.6, 0.9, 0.98, 1.0])
    v = 2.0 + 1.4 * z ** 0.25 + 4 * np.clip(z - 0.95, 0, None)
    u = ocv.grid_branch(z, v)
    assert ocv.monotonic_violations(u) == [] and u.min() >= v.min() - 1e-12 and u.max() <= v.max() + 1e-12


def test_ocv_is_subtracted_per_sample_along_soc():
    Q, I, R0 = 172.4, 56.0, 3e-4
    t = np.arange(0, 600.0, 1.0)
    z = 0.5 - I * t / 3600 / Q
    f = lambda zz: 3.0 + 0.4 * zz
    V = f(z) - R0 * I
    assert np.allclose(V - f(z), -R0 * I)                       # per sample: flat overpotential
    assert abs((V - f(z[0]))[-1] - (V - f(z[0]))[0]) > 1e-3     # single start value: 22 mV of false drift


def test_arrhenius_fit_recovers_known_activation_energy():
    Ea, R25 = np.array([0.30, 0.45]), np.array([3e-4, 1e-3])
    T = np.array([15.0, 25.0, 45.0])
    y = R25 * np.exp(Ea / ecm.KB_EV * (1 / (T[:, None] + 273.15) - 1 / 298.15))
    a, b = ecm.arrhenius_fit(T, y)
    assert np.allclose(b * ecm.KB_EV, Ea)
    assert np.allclose(ecm.arrhenius_eval(a, b, [-30.0])[0], R25 * np.exp(Ea / ecm.KB_EV * (1 / 243.15 - 1 / 298.15)))


@pytest.mark.parametrize("name,kind", [("ocv.csv", "ocv"), ("r0.csv", "r"), ("r1.csv", "r"), ("c1.csv", "c"),
                                       ("entropy.csv", "entropy")])
def test_team_files_round_trip_through_reader_and_writer(tmp_path, name, kind):
    a = load.read_table(load.DATAPACK / name)
    load.write_table(a, tmp_path / name, kind)
    b = load.read_table(tmp_path / name)
    keys = ["SOC"] if kind == "entropy" else ["Temperature", "SOC"]
    assert list(a.columns) == list(b.columns) and a[keys].equals(b[keys])
    assert np.allclose(a.drop(columns=keys).to_numpy(), b.drop(columns=keys).to_numpy(), rtol=1e-6, atol=1e-9)
    assert (tmp_path / name).read_bytes()[:3] != b"\xef\xbb\xbf"


def test_simulation_with_two_branches_is_the_sum_of_two_step_responses():
    import pandas as pd
    from reparam import validate
    R0, R1, tau1, R2, tau2, I = 3e-4, 2e-4, 300.0, 1e-4, 8.0, 100.0
    t = np.arange(0, 600.0, 0.5)
    df = pd.DataFrame(dict(t=t, current=np.full_like(t, I)))
    rc = pd.DataFrame(dict(SOC=[0.0, 1.0], R0=R0, R1=R1, tau1=tau1, R2=R2, tau2=tau2))
    flat = lambda z: np.full_like(np.asarray(z, float), 3.0)
    v = validate.simulate(df, np.full_like(t, 0.5), rc, (flat, flat))
    want = 3.0 - R0 * I - R1 * I * (1 - np.exp(-t / tau1)) - R2 * I * (1 - np.exp(-t / tau2))
    assert np.max(np.abs(v - want)) < 1e-12
    one = validate.simulate(df, np.full_like(t, 0.5), rc.drop(columns=["R2", "tau2"]), (flat, flat))
    assert np.max(np.abs(one - (want + R2 * I * (1 - np.exp(-t / tau2))))) < 1e-12


def test_stale_second_row_of_a_step_is_masked_with_the_mid_switch_row():
    seg = np.array([0, 0, 1, 1, 1, 1, 2, 2, 2])
    V = np.array([3.0, 3.0, 3.0, 3.0, 2.9, 2.9, 2.9, 2.95, 2.95])
    bad = ecm.unusable_rows(seg, V)
    assert bad.tolist() == [False, False, True, True, False, False, True, False, False]
