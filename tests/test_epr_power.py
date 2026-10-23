"""benchmarks/epr_power.py: C4's EPR comparison re-read from the committed rows (EXP-20261020-GB-1).

The expected values are the preliminary estimate of the Phase 0 plan (point 15), computed
before this record from the same committed CSV; the gate lines are copied from C4's
committed REPORT_REAL.md.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

from grnbench.eval.power import paired_t_power

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "benchmarks"))
import epr_power  # noqa: E402

C4 = ROOT / "benchmarks" / "results" / "EXP-20260925-LM-2"
AUTO_LINE = "| pearson_motif | pearson_signed | 20 | +0.120 | [-0.046, +0.286] | 0.219 | 1 | unresolved |"
BASE_LINE = "| pearson_motif | pearson_signed | 20 | +0.152 | [-0.014, +0.318] | 0.105 | 1 | unresolved |"
# the registered report code on the 19 other datasets (found by the independent review of PR #129)
BASE_LOO_CHANGES = {
    "10x-3-v3__homo-sapiens__prostate-gland__20de8570":
        "| pearson_motif | pearson_signed | 19 | +0.190 | [+0.035, +0.345] | 0.0288 | 1 | pearson_motif |",
    "10x-3-v3__homo-sapiens__respiratory-airway__eec804b9":
        "| pearson_motif | pearson_signed | 19 | +0.188 | [+0.030, +0.345] | 0.0331 | 1 | pearson_motif |"}


@pytest.fixture(scope="module")
def rows():
    return epr_power.load_rows(C4 / "all_rows_real.csv")


@pytest.fixture(scope="module")
def report_text():
    return (C4 / "REPORT_REAL.md").read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def result(rows, report_text):
    return epr_power.analyse(rows, report_text)


def test_gate_reprints_the_registered_epr_comparisons_of_c4(rows, report_text):
    lines = epr_power.registered_epr_lines(rows)
    epr_power.check_against_report(lines, report_text)
    assert AUTO_LINE in lines and BASE_LINE in lines
    assert sum(line.startswith("#### Declared comparisons, EPR") for line in lines) == 2


def test_gate_refuses_rows_that_no_longer_give_the_registered_report(rows, report_text):
    changed = rows.copy()
    changed.loc[changed.index[0], "m_pearson_motif_auto_epr"] += 0.5
    with pytest.raises(epr_power.GateError, match="pearson_motif"):
        epr_power.check_against_report(epr_power.registered_epr_lines(changed), report_text)


def test_gate_refuses_a_report_with_an_epr_comparison_the_code_does_not_print(rows, report_text):
    extra = AUTO_LINE.replace("pearson_signed", "spearman_signed")
    with pytest.raises(epr_power.GateError):
        epr_power.check_against_report(epr_power.registered_epr_lines(rows),
                                       report_text.replace(AUTO_LINE, AUTO_LINE + "\n" + extra))


def test_inputs_are_hashed_as_committed_and_match_c4_sha256sums():
    sums = epr_power.verify_inputs(C4 / "all_rows_real.csv", C4 / "REPORT_REAL.md")
    assert sums == {"all_rows_real.csv": "08a79ccae6ca18e47fa14c8ef4b66f66f27c17c531cd71beb14c009d4a3eed14",
                    "REPORT_REAL.md": "b7f4963f79614fa1b4c91d6a6c9ab9983805d3fbbd547b0f1d4a0a154a34979c"}


def test_an_input_that_differs_from_c4_sha256sums_is_refused(tmp_path):
    for name in ("all_rows_real.csv", "REPORT_REAL.md", "SHA256SUMS"):
        (tmp_path / name).write_bytes((C4 / name).read_bytes())
    with open(tmp_path / "all_rows_real.csv", "ab") as fh:
        fh.write(b"\n")
    with pytest.raises(epr_power.GateError, match="all_rows_real.csv"):
        epr_power.verify_inputs(tmp_path / "all_rows_real.csv", tmp_path / "REPORT_REAL.md")


def test_a_zero_epr_drops_only_the_log_ratio_diagnostic(rows):
    zero = rows.copy()
    zero.loc[zero.index[0], "m_pearson_signed_auto_epr"] = 0.0
    diag = epr_power.diagnose(zero, "auto")
    assert diag["n"] == 20 and diag["log_ratio"] is None


def test_an_unreachable_sample_size_is_reported_not_raised():
    assert "n_max" in epr_power._n_or_reason(1e-4, 0.0125, 0.8)
    assert epr_power._n_or_reason(-0.1, 0.0125, 0.8).startswith("none")


def test_other_power_errors_are_not_turned_into_a_reason():
    with pytest.raises(ValueError):
        epr_power._n_or_reason(float("nan"), 0.0125, 0.8)


def test_a_missing_auprc_score_names_the_dataset(rows):
    gap = rows.copy()
    gap.loc[gap.index[3], "m_pearson_signed_auto_auprc_ratio"] = float("nan")
    with pytest.raises(ValueError, match=gap.loc[gap.index[3], "task"].replace("cxg__", "")):
        epr_power.diagnose(gap, "auto")


def test_a_duplicate_dataset_id_is_refused(rows):
    twice = rows.copy()
    twice.loc[twice.index[1], "task"] = twice.loc[twice.index[0], "task"]
    with pytest.raises(ValueError, match="unique"):
        epr_power.diagnose(twice, "auto")


def test_the_diagnosed_pair_matches_its_gate_line(result):
    for arm, d in result["arms"].items():
        prefix = f"| pearson_motif | pearson_signed | {d['n']} | {d['mean_diff']:+.3f} |"
        assert any(line.startswith(prefix) for line in result["gate"]["lines"])
        assert d["dz_bootstrap_zero_sd_resamples"] == 0


def test_a_missing_score_column_is_an_error_not_a_smaller_table(rows, tmp_path):
    path = tmp_path / "rows.csv"
    rows.drop(columns=["m_pearson_signed_baseline_epr"]).to_csv(path, index=False)
    with pytest.raises(ValueError, match="m_pearson_signed_baseline_epr"):
        epr_power.load_rows(path)


def test_diagnosis_reproduces_the_preliminary_estimate(result):
    auto, base = result["arms"]["auto"], result["arms"]["baseline"]
    assert auto["n"] == base["n"] == 20
    assert auto["mean_diff"] == pytest.approx(0.12025588513894907, abs=1e-12)
    assert auto["sd_diff"] == pytest.approx(0.35484266253376456, abs=1e-12)
    assert auto["dz"] == pytest.approx(0.33889917373592665, abs=1e-12)
    assert auto["p_one_sided"] == pytest.approx(0.07304013339009496, abs=1e-12)
    assert auto["loo_t_range"] == pytest.approx([1.135225444991594, 2.0472590538313535], abs=1e-9)
    assert auto["n_positive"] == 13 and base["n_positive"] == 14
    assert auto["dz_bootstrap_10_90"] == pytest.approx([0.057031221235723274, 0.6843057632785443], abs=1e-12)
    assert base["dz"] == pytest.approx(0.4293240240230651, abs=1e-12)
    assert base["sd_diff"] == pytest.approx(0.3546578579360856, abs=1e-12)
    assert auto["auprc_ratio"]["dz"] == pytest.approx(1.899681823524661, abs=1e-9)
    assert auto["log_ratio"]["dz"] == pytest.approx(0.4635062129003433, abs=1e-9)


def test_planning_tables_reproduce_the_preliminary_numbers(result):
    observed = {arm: result["planning"][arm]["observed"] for arm in ("auto", "baseline")}
    assert observed["auto"]["alpha/3"]["n80"] == 80 and observed["auto"]["alpha/4"]["n80"] == 86
    assert observed["baseline"]["alpha/3"]["n80"] == 51 and observed["baseline"]["alpha/4"]["n80"] == 55
    smallest = result["planning"]["auto"]["effects"]["0.10"]["alpha/2"]
    assert (smallest["n80"], smallest["n90"]) == (101, 135)
    grid = {row["n"]: row for row in result["grid"]["auto"]["alpha/2"]}
    assert grid[60]["power_at_0.10"] == pytest.approx(0.57, abs=0.01)
    assert grid[60]["mde80_epr"] == pytest.approx(0.13, abs=0.005)
    assert grid[100]["mde80_epr"] == pytest.approx(0.10, abs=0.005)
    for arm in ("auto", "baseline"):
        for level in ("alpha/4", "alpha/2"):
            at = {row["n"]: row for row in result["grid"][arm][level]}
            assert 0.085 <= at[120]["mde80_epr"] <= 0.105
            assert 0.065 <= at[200]["mde80_epr"] <= 0.085


def test_c4_as_run_is_described_by_its_smallest_detectable_effect(result):
    # sensitivity, not observed power: what n = 20 could detect at 80 %
    at20 = result["sensitivity_n20"]["auto"]["alpha/2"]
    assert at20["mde80_dz"] == pytest.approx(0.66, abs=0.01)
    assert at20["mde80_epr"] == pytest.approx(0.23, abs=0.01)


def _floats(obj):
    if isinstance(obj, dict):
        for v in obj.values():
            yield from _floats(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from _floats(v)
    elif isinstance(obj, float):
        yield obj


def test_no_output_reports_power_computed_from_c4s_own_estimate(tmp_path):
    # record section 1: power is never computed from C4's estimate; the numbers at C4's estimate
    # are only the datasets a new study would need (with the bootstrap band of dz)
    epr_power.main(["--out", str(tmp_path)])
    data = json.loads((tmp_path / "power.json").read_text(encoding="utf-8"))
    report = (tmp_path / "POWER.md").read_text(encoding="utf-8")
    numbers = np.array(list(_floats(data)))
    for arm, d in data["arms"].items():
        for alpha in epr_power.LEVELS.values():
            for sides in (1, 2):
                at_estimate = np.array([paired_t_power(d["dz"], n, alpha, sides=sides) for n in range(3, 251)])
                assert not np.isclose(numbers[:, None], at_estimate[None, :], rtol=0, atol=1e-9).any()
            assert f"{paired_t_power(d['dz'], d['n'], alpha):.3f}" not in report
    assert "power at c4's estimate" not in report.lower()
    curves = epr_power.power_curves(data)
    assert curves and all(c["dz"] not in {d["dz"] for d in data["arms"].values()} for c in curves)


def test_leave_one_out_reruns_the_registered_report_code_and_names_each_verdict_change(result):
    # C4's own procedure (models_section: paired t, Holm over 4 directed hypotheses) on the 19 others
    loo = result["leave_one_out"]
    for arm in ("auto", "baseline"):
        assert loo[arm]["verdict_all"] == "unresolved" and loo[arm]["drops"] == 20
    assert loo["auto"]["changed"] == []
    assert {c["dropped"]: c["line"] for c in loo["baseline"]["changed"]} == BASE_LOO_CHANGES


def test_main_writes_the_report_json_and_figure(tmp_path):
    epr_power.main(["--out", str(tmp_path)])
    report = (tmp_path / "POWER.md").read_text(encoding="utf-8")
    data = json.loads((tmp_path / "power.json").read_text(encoding="utf-8"))
    assert AUTO_LINE in report and BASE_LINE in report
    assert "not resolved is not zero" in report.lower()
    assert "auto: verdict with all 20 datasets `unresolved`; it changes in 0 of 20 single-dataset drops" in report
    assert "baseline: verdict with all 20 datasets `unresolved`; it changes in 2 of 20 single-dataset drops" in report
    assert all(line in report for line in BASE_LOO_CHANGES.values())
    assert data["arms"]["auto"]["n"] == 20 and data["gate"]["matches_report"] is True
    assert data["inputs"]["benchmarks/results/EXP-20260925-LM-2/all_rows_real.csv"].startswith("08a79cca")
    assert (tmp_path / "power_epr.png").stat().st_size > 0
