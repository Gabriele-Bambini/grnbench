"""Calibration gates use finite Monte Carlo bounds and expose model violations."""
import json

import pytest

from grnbench.calibration import binomial_bounds, calibration_gate, run_calibration


def config(reps=8):
    return {"seed": 9281, "replicates": reps, "alpha": .05, "margin": 0.,
            "tolerance": .025, "mc_confidence": .95,
            "scenarios": [{"name": "iid", "n_units": 8, "means": [0., 0., .8],
                           "shared_noise": .5, "independent_noise": 1., "repeat_units": 1}]}


def test_exact_binomial_bound_zero_errors_is_not_zero_risk():
    bounds = binomial_bounds(0, 20, .95)
    assert bounds["upper"] == pytest.approx(1 - .05 ** (1 / 20))
    assert bounds["lower"] == 0.
    assert calibration_gate(0, 20, 20, alpha=.05, tolerance=.025)["status"] == "inconclusive"
    assert calibration_gate(0, 1000, 1000, alpha=.05, tolerance=.025)["status"] == "pass"
    assert calibration_gate(500, 500, 1000, alpha=.05, tolerance=.025)["status"] == "fail"


def test_pilot_reproducible_and_known_truth_recorded(tmp_path):
    a = run_calibration(config(), tmp_path / "a")
    b = run_calibration(config(), tmp_path / "b")
    assert a["scenarios"] == b["scenarios"]
    assert a["scope"] == "score_level_gaussian_pilot"
    assert a["biological_validation"] is False
    assert (tmp_path / "a" / "experiments.csv").exists()
    assert json.loads((tmp_path / "a" / "manifest.json").read_text())["config_sha256"]
    with pytest.raises(FileExistsError):
        run_calibration(config(), tmp_path / "a")


def test_duplicated_units_are_declared_as_assumption_violation(tmp_path):
    c = config()
    c["scenarios"][0]["repeat_units"] = 4
    result = run_calibration(c, tmp_path / "out")
    row = result["scenarios"][0]
    assert row["independence_assumption_holds"] is False
    assert row["n_independent_units"] == 2


@pytest.mark.parametrize("field,value", [("replicates", 0), ("alpha", 1), ("tolerance", -.01),
                                         ("mc_confidence", .1), ("margin", -1)])
def test_invalid_config_rejected_before_creating_output(tmp_path, field, value):
    c = config()
    c[field] = value
    out = tmp_path / "bad"
    with pytest.raises(ValueError):
        run_calibration(c, out)
    assert not out.exists()


def test_failed_execution_records_failure_manifest(tmp_path, monkeypatch):
    import grnbench.calibration as module
    def broken(*args, **kwargs):
        raise RuntimeError("injected simulator failure")
    monkeypatch.setattr(module, "compare_units", broken)
    with pytest.raises(RuntimeError):
        run_calibration(config(), tmp_path / "failed")
    manifest = json.loads((tmp_path / "failed" / "manifest.json").read_text())
    assert manifest["status"] == "failed"
    assert "injected simulator failure" in manifest["error"]


def test_truth_is_invariant_to_score_units(tmp_path):
    """Changing units must not turn a true small effect into a false call."""
    c = config(reps=20)
    c["scenarios"] = [{"name": "small_effect", "n_units": 10,
                       "means": [0., 5e-13], "shared_noise": 0.,
                       "independent_noise": 1e-14}]
    small = run_calibration(c, tmp_path / "small")["scenarios"][0]
    c["scenarios"][0].update(means=[0., 50.], independent_noise=1.)
    large = run_calibration(c, tmp_path / "large")["scenarios"][0]
    assert small["mean_power"] == large["mean_power"] == 1.
    assert small["gate"]["familywise_false_superiority"]["events"] == 0
