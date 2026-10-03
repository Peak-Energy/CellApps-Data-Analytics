"""2RC correctness test: a synthetic HPPC file with a fast and a slow branch of known size goes through the
same loader, SOC axis, windows and joint fit as the real data.

Truth: R0 = 0.35 mOhm; slow branch 0.25 mOhm at 360 s; fast branch 0.12 mOhm at 8 s. Logging as in the
real files: 0.1 s on the pulses and their rests (60 s after a discharge pulse, 900 s after a charge pulse),
10 s on the SOC moves, 300 s on the 2 h rests, the first row of each step captured mid-switch.
"""
import sys
from pathlib import Path

import numpy as np
import pytest
import yaml
from scipy.interpolate import PchipInterpolator

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from reparam import ecm, load  # noqa: E402
import synth  # noqa: E402

CFG = yaml.safe_load(open(Path(__file__).resolve().parents[1] / "config.yaml"))
TRUE = dict(R0=0.35e-3, R1=0.25e-3, tau=360.0, R2=0.12e-3, tau2=8.0)
Q, Z0 = 172.4, 0.995


def _fit(tmp_path_factory, noise_uV):
    t = load.read_table(load.DATAPACK / "ocv.csv")
    t = t[t.Temperature == 25]
    ocv_fn = PchipInterpolator(t.SOC.to_numpy(), t["OCV-Discharge"].to_numpy(), extrapolate=True)
    p = tmp_path_factory.mktemp("synth2") / "CV_0010_CS2D0000_Pr_HPPC_25C_.csv"
    synth.write(p, ocv_fn, Q_ah=Q, z0=Z0, seed=5, steps=synth.schedule(n_points=5, move_s=2230.0),
                noise_uV=noise_uV, **TRUE)
    cfg = dict(CFG, soc={"capacity_ah": {25: Q}})
    prep = ecm.prepare(cfg, t, files={25: [p]})
    fits, taus, diag = ecm.fit_2rc(prep, cfg["two_rc"])
    return fits, taus, diag


@pytest.fixture(scope="module")
def clean(tmp_path_factory):
    return _fit(tmp_path_factory, 0.0)


@pytest.fixture(scope="module")
def noisy(tmp_path_factory):
    return _fit(tmp_path_factory, 37.0)


def _branches(fits):
    """(R_fast, R_slow) per SOC point from the rows the tables use: fast from pulses, slow from moves."""
    f = fits[~ecm.excluded(fits)]
    fast = f[f.kind == "pulse"].groupby("point").R_fast.median()
    slow = f[f.kind == "move"].groupby("point").R_slow.median()
    return fast, slow


def test_without_noise_both_time_constants_and_resistances_come_back(clean):
    fits, taus, diag = clean
    assert taus["fast"] == pytest.approx(TRUE["tau2"], rel=0.03)
    assert taus["slow"] == pytest.approx(TRUE["tau"], rel=0.03)
    fast, slow = _branches(fits)
    assert len(fast) >= 3 and len(slow) >= 3
    assert np.allclose(fast, TRUE["R2"], rtol=0.03) and np.allclose(slow, TRUE["R1"], rtol=0.03)
    assert np.allclose(fits.R0, TRUE["R0"], rtol=0.02)      # R0 is the 1RC pulse fit, kept unchanged on purpose
    assert not diag["at_grid_edge"]


def test_with_noise_the_answers_stay_within_five_percent(noisy):
    fits, taus, _ = noisy
    assert taus["fast"] == pytest.approx(TRUE["tau2"], rel=0.05)
    assert taus["slow"] == pytest.approx(TRUE["tau"], rel=0.05)
    fast, slow = _branches(fits)
    assert np.allclose(fast, TRUE["R2"], rtol=0.05) and np.allclose(slow, TRUE["R1"], rtol=0.05)


def test_two_branches_beat_one_on_the_same_windows(noisy):
    """The best single time constant on the same windows must leave more error than the pair."""
    _, _, diag = noisy
    assert diag["rms_mV_best_pair"] < 0.5 * diag["rms_mV_best_single_tau"]


def test_the_two_branches_are_fitted_jointly_per_point(noisy):
    """Every window of a SOC point carries the same pair of resistances."""
    fits, _, _ = noisy
    for _, g in fits[~ecm.excluded(fits)].groupby("point"):
        assert g.R_fast.nunique() == 1 and g.R_slow.nunique() == 1
