"""Main correctness test: a synthetic Neware HPPC file with known answers goes through the same loader,
time reconstruction, SOC axis, windows and fitter as the real data.

Truth: R0 = 0.35 mOhm, R1 = 0.25 mOhm, tau1 = 360 s, Q = 172.4 Ah, OCV = team 25 C discharge curve.
Realistic logging: 37 uV voltage noise, 0.1 s pulse logging with the time printed to the whole second,
300 s logging on long rests, 10 s on the SOC moves, the first row of each step captured mid-switch,
and the top rest at 99.5 % SOC so the SOC-axis shift has a known answer.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.interpolate import PchipInterpolator

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from reparam import ecm, load  # noqa: E402
import synth  # noqa: E402

TRUE = dict(R0=0.35e-3, R1=0.25e-3, tau=360.0)
Q, Z0 = 172.4, 0.995
BOUNDS = ([0.0, 0.0, 0.05], [5e-3, 5e-2, 2e4])


def _make(tmp_path_factory, noise_uV):
    t = load.read_table(load.DATAPACK / "ocv.csv")
    t = t[t.Temperature == 25]
    ocv_fn = PchipInterpolator(t.SOC.to_numpy(), t["OCV-Discharge"].to_numpy(), extrapolate=True)
    p = tmp_path_factory.mktemp("synth") / "CV_0010_CS2D0000_Pr_HPPC_25C_.csv"
    truth = synth.write(p, ocv_fn, Q_ah=Q, z0=Z0, seed=3, steps=synth.schedule(n_points=4, move_s=2230.0),
                        noise_uV=noise_uV, **TRUE)
    df, steps = load.load(p)
    return df, steps, ocv_fn, truth


@pytest.fixture(scope="module")
def noisy(tmp_path_factory):
    return _make(tmp_path_factory, 37.0)


@pytest.fixture(scope="module")
def clean(tmp_path_factory):
    return _make(tmp_path_factory, 0.0)


def _fits(rec, kind, rest_after=None):
    df, steps, ocv_fn, _ = rec
    z, _ = ecm.soc_axis(df, steps, Q, ocv_fn)
    out = []
    for _, r in ecm.structure(steps).query("kind == @kind").iterrows():
        if rest_after and not rest_after[0] <= steps.loc[r.seg + 1, "dur_s"] <= rest_after[1]:
            continue
        d, nh = ecm.window(df, steps, r.seg, kind)
        t, i, eta, mask, w, seg = ecm.window_arrays(d, nh, z[d.index.to_numpy()], ocv_fn)
        f = ecm.fit_window(t, i, eta, mask, BOUNDS, w, nh, seg)
        out.append(dict(R0=f["R0"], R1=f["R1"], tau=f["tau"], se_R1=f["se"][1], se_tau=f["se"][2]))
    return pd.DataFrame(out)


def test_loader_reproduces_the_true_time_and_current(noisy):
    df, _, _, truth = noisy
    assert len(df) == len(truth["t"]) and np.max(np.abs(df.t.to_numpy() - truth["t"])) < 1e-6
    assert np.allclose(df.current.to_numpy(), truth["i"], atol=1e-3)


def test_soc_axis_shift_is_recovered(noisy):
    df, steps, ocv_fn, truth = noisy
    z, info = ecm.soc_axis(df, steps, Q, ocv_fn)
    assert info["soc_offset_pct"] == pytest.approx(100 * (Z0 - 1.0), abs=0.05)
    assert np.max(np.abs(z - truth["z"])) < 5e-4


@pytest.mark.parametrize("rest", [(50, 70), (800, 8000)])
def test_without_noise_every_pulse_window_returns_the_truth(clean, rest):
    """Left-over error is only the 1 uV print resolution of the voltage (worth 1-2 % in R1 and tau)."""
    f = _fits(clean, "pulse", rest)
    assert len(f) == 12 and np.allclose(f.R0, TRUE["R0"], rtol=1e-3)
    assert np.allclose(f.R1, TRUE["R1"], rtol=0.03) and np.allclose(f.tau, TRUE["tau"], rtol=0.03)


def test_with_noise_pulse_errors_stay_within_three_standard_errors(noisy):
    """A 10 s pulse charges the 360 s branch only 2.7 %, so R1 and tau are noise-limited; the fit must
    report that honestly. R0 stays within 1 %."""
    f = _fits(noisy, "pulse", (800, 8000))
    assert np.allclose(f.R0, TRUE["R0"], rtol=0.01)
    for k in ["R1", "tau"]:
        assert ((f[k] - TRUE[k]).abs() / f[f"se_{k}"] < 3).all(), f.to_string()


def test_soc_move_windows_recover_r1_and_tau(noisy):
    f = _fits(noisy, "move")
    assert len(f) == 3
    assert np.allclose(f.R1, TRUE["R1"], rtol=0.03) and np.allclose(f.tau, TRUE["tau"], rtol=0.05)


def test_shared_tau_solve_recovers_r1_at_the_true_tau(noisy):
    df, steps, ocv_fn, _ = noisy
    z, _ = ecm.soc_axis(df, steps, Q, ocv_fn)
    for _, r in ecm.structure(steps).query("kind == 'move'").iterrows():
        d, nh = ecm.window(df, steps, r.seg, "move")
        t, i, eta, mask, w, seg = ecm.window_arrays(d, nh, z[d.index.to_numpy()], ocv_fn)
        _, R1, _ = ecm.solve_given_tau(dict(t=t, i=i, eta=eta, mask=mask, w=w, nh=nh, step=seg, R0=TRUE["R0"]),
                                       TRUE["tau"])
        assert R1 == pytest.approx(TRUE["R1"], rel=0.02)


def test_mid_switch_row_gives_only_a_third_of_r0(noisy):
    """The team's R0 convention (voltage change on the first logged row) recovers ~30 % of R0 here."""
    df, steps, _, _ = noisy
    r = ecm.structure(steps).query("kind == 'pulse'").iloc[0]
    d, nh = ecm.window(df, steps, r.seg, "pulse")
    v = d.voltage.to_numpy()[nh:]
    assert (v[0] - v[1]) / r.amp_A / TRUE["R0"] == pytest.approx(0.3, abs=0.05)
