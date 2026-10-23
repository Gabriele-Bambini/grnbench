"""benchmarks/compare_runs.py: judging a rerun of a registered experiment against its committed
results from files (rows, report, the record's primary comparison), not by eye."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
from compare_runs import MIN_UNITS, compare_report, main  # noqa: E402

import pytest  # noqa: E402
from grnbench.eval.comparison import compare_units  # noqa: E402

PRIMARY = "pm > ps | auprc_ratio | auto"
SCORES = {"d1": (1.50, 1.30), "d2": (1.42, 1.31), "d3": (1.61, 1.40), "d4": (1.33, 1.29)}


def _report(rows: dict, seconds: float, eol: str = "\n") -> str:
    lines = ["# Engine mass test on real datasets - 4 datasets", "",
             "#### Declared comparisons, AUPRC ratio, auto pipeline (paired t, Holm over 2 directed hypotheses)", "",
             "| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |",
             "|---|---|---|---|---|---|---|---|",
             "| pm | ps | 4 | +0.140 | [+0.010, +0.270] | 0.0221 | 1 | pm |", "",
             "## Per dataset", "",
             "| dataset | outcome | s |", "|---|---|---|"]
    lines += [f"| {t} | {r['outcome']} | {seconds:.1f} |" for t, r in rows.items()]
    return eol.join(lines) + eol


def _run(root: Path, *, seconds: float = 10.0, eol: str = "\n", change: dict | None = None) -> Path:
    rows = {}
    for t, (pm, ps) in SCORES.items():
        rows[t] = {"task": t, "outcome": "ok, scored", "scored": True, "reference_edges": 900,
                   "m_pm_auto_auprc_ratio": pm, "m_ps_auto_auprc_ratio": ps,
                   "m_pm_auto_epr": 2.0, "m_pm_auto_seconds": seconds, "seconds": seconds}
    for (t, key), value in (change or {}).items():
        if value is None:
            rows.pop(t)
        else:
            rows[t][key] = value
    (root / "rows").mkdir(parents=True)
    for t, r in rows.items():
        (root / "rows" / f"{t}.json").write_text(json.dumps(r), encoding="utf-8")
    (root / "REPORT_REAL.md").write_bytes(_report(rows, seconds, eol).encode())
    return root


def _record(tmp_path: Path) -> Path:
    rec = tmp_path / "EXP-x.md"
    rec.write_text("## 1. Registration\n\n- primary comparison: " + PRIMARY + "\n  (restated before the run)\n",
                   encoding="utf-8")
    return rec


def test_identical_rerun_matches_across_line_endings_and_timings(tmp_path, capsys):
    committed = _run(tmp_path / "committed", seconds=10.0, eol="\r\n")      # written on Windows
    rerun = _run(tmp_path / "rerun", seconds=37.5)                             # Linux, slower machine
    code = main([str(committed), str(rerun), "--record", str(_record(tmp_path)),
                 "--json", str(tmp_path / "diff.json")])
    out = capsys.readouterr().out
    assert code == 0, out
    assert "tasks: same (4 committed, 4 rerun)" in out
    # 4 rows x (reference_edges + 3 scores) numeric; 4 bool; task + outcome; 2 timing keys per row
    assert "compared 16 numeric, 4 bool and 8 other fields; 8 timing fields not compared; exact" in out
    assert "report: identical (line endings and timing columns ignored)" in out
    assert f"primary: MATCH - {PRIMARY}" in out
    assert "n=4 mean +0.140" in out
    assert out.rstrip().endswith("verdict: MATCH")
    assert json.loads((tmp_path / "diff.json").read_text())["verdict"] == "MATCH"


def test_a_changed_secondary_value_fails_the_verdict_but_not_the_primary(tmp_path, capsys):
    committed = _run(tmp_path / "committed")
    rerun = _run(tmp_path / "rerun", change={("d2", "m_pm_auto_epr"): 2.5})
    code = main([str(committed), str(rerun), "--primary", PRIMARY])
    out = capsys.readouterr().out
    assert code == 1
    assert "rows: 1 DIFFERENT" in out and "d2\tm_pm_auto_epr\t2.0 -> 2.5" in out
    assert f"primary: MATCH - {PRIMARY}" in out
    assert "verdict: DIFFERENT" in out


def test_a_changed_primary_input_is_named(tmp_path, capsys):
    committed = _run(tmp_path / "committed")
    rerun = _run(tmp_path / "rerun", change={("d3", "m_ps_auto_auprc_ratio"): 1.41})
    code = main([str(committed), str(rerun), "--primary", PRIMARY])
    out = capsys.readouterr().out
    assert code == 1
    assert f"primary: DIFFERENT - {PRIMARY}" in out
    assert "d3\tps\t1.4 -> 1.41" in out


def test_tolerance_is_opt_in(tmp_path, capsys):
    committed = _run(tmp_path / "committed")
    rerun = _run(tmp_path / "rerun", change={("d1", "m_pm_auto_auprc_ratio"): 1.50 * (1 + 1e-12)})
    assert main([str(committed), str(rerun), "--primary", PRIMARY]) == 1
    assert main([str(committed), str(rerun), "--primary", PRIMARY, "--rtol", "1e-9"]) == 0
    assert "tolerance rtol=1e-09" in capsys.readouterr().out


def test_missing_task_and_type_change_are_differences(tmp_path, capsys):
    committed = _run(tmp_path / "committed")
    rerun = _run(tmp_path / "rerun", change={("d4", "task"): None, ("d1", "scored"): 1})
    code = main([str(committed), str(rerun), "--primary", PRIMARY])
    out = capsys.readouterr().out
    assert code == 1
    assert "tasks: DIFFERENT (4 committed, 3 rerun; only committed: d4)" in out
    assert "d1\tscored\tTrue -> 1" in out                       # a bool is not the number 1
    assert "unit d4 only in committed" in out


def test_a_unit_without_a_numeric_score_leaves_the_pair_on_both_sides(tmp_path, capsys):
    committed = _run(tmp_path / "committed", change={("d2", "m_ps_auto_auprc_ratio"): "error"})
    rerun = _run(tmp_path / "rerun", change={("d2", "m_ps_auto_auprc_ratio"): "error"})
    assert main([str(committed), str(rerun), "--primary", PRIMARY]) == 0
    assert "committed: n=3 mean" in capsys.readouterr().out    # as real_bench drops it from the pair


def test_report_differences_outside_timing_columns_are_reported():
    old = "| dataset | outcome | s |\n|---|---|---|\n| d1 | ok | 3.0 |\n"
    assert compare_report(old, old.replace("3.0", "9.9")) == []
    assert compare_report(old, old.replace("ok", "error")) == ["- line 3: '| d1 | ok | 3.0 |'",
                                                               "+ line 3: '| d1 | error | 3.0 |'"]
    assert compare_report(old, old + "extra\n") == ["+ line 4: 'extra'"]


def test_a_newer_report_format_is_information_not_a_difference(tmp_path, capsys):
    committed = _run(tmp_path / "committed")
    rerun = _run(tmp_path / "rerun")
    report = rerun / "REPORT_REAL.md"                 # same rows, report written by newer code
    lines = report.read_text(encoding="utf-8").splitlines()
    report.write_text("\n".join(lines[:1] + ["", "- 4 planned datasets", "- permanent errors: 0"] + lines[1:]) + "\n",
                      encoding="utf-8")
    assert main([str(committed), str(rerun), "--primary", PRIMARY]) == 0
    out = capsys.readouterr().out
    assert "report: 3 lines differ (information only" in out     # only the inserted lines, aligned
    assert f"primary: MATCH - {PRIMARY}" in out and out.rstrip().endswith("verdict: MATCH")


def test_the_primary_report_row_is_shown_but_the_inputs_decide(tmp_path, capsys):
    committed = _run(tmp_path / "committed")
    rerun = _run(tmp_path / "rerun")
    report = rerun / "REPORT_REAL.md"                 # e.g. a different Holm family in newer code
    report.write_text(report.read_text(encoding="utf-8").replace("| 0.0221 |", "| 0.0442 |"), encoding="utf-8")
    assert main([str(committed), str(rerun), "--primary", PRIMARY]) == 0
    out = capsys.readouterr().out
    assert f"primary: MATCH - {PRIMARY}" in out
    assert "report row differs (information)" in out and "0.0442" in out


def test_no_primary_given_is_an_error_not_a_match(tmp_path, capsys):
    committed, rerun = _run(tmp_path / "committed"), _run(tmp_path / "rerun")
    assert main([str(committed), str(rerun)]) == 2
    captured = capsys.readouterr()
    assert "no primary comparison given" in captured.err and "--rows-only" in captured.err
    assert "MATCH" not in captured.out                      # nothing certified


def test_rows_only_is_explicit_and_says_so(tmp_path, capsys):
    committed, rerun = _run(tmp_path / "committed"), _run(tmp_path / "rerun")
    assert main([str(committed), str(rerun), "--rows-only", "--json", str(tmp_path / "d.json")]) == 0
    out = capsys.readouterr().out
    assert "primary: not checked (--rows-only)" in out
    assert out.rstrip().endswith("verdict: MATCH (rows only: no primary comparison checked)")
    assert json.loads((tmp_path / "d.json").read_text())["mode"] == "rows-only"
    changed = _run(tmp_path / "changed", change={("d2", "m_pm_auto_epr"): 2.5})   # rows still decide
    assert main([str(committed), str(changed), "--rows-only"]) == 1
    with pytest.raises(SystemExit) as exc:                 # a primary and --rows-only contradict
        main([str(committed), str(rerun), "--rows-only", "--primary", PRIMARY])
    assert exc.value.code == 2


@pytest.mark.parametrize("given", ["not a comparison", "pm > ps | auprc_ratio", "",
                                   PRIMARY + "; broken second comparison"])
def test_an_unparseable_primary_fails_loudly(tmp_path, capsys, given):
    committed, rerun = _run(tmp_path / "committed"), _run(tmp_path / "rerun")
    assert main([str(committed), str(rerun), "--primary", given]) == 2
    captured = capsys.readouterr()
    assert "no primary comparison could be parsed from --primary" in captured.err
    assert "MATCH" not in captured.out


def test_a_record_without_a_parsable_primary_fails_loudly(tmp_path, capsys):
    committed, rerun = _run(tmp_path / "committed"), _run(tmp_path / "rerun")
    for text in ("## 1. Registration\n\n- type: exploratory\n",                        # field absent
                 "## 1. Registration\n\n- primary comparison: pearson beats clr on AUPRC\n"):   # not A > B | metric | arm
        record = tmp_path / "EXP-y.md"
        record.write_text(text, encoding="utf-8")
        assert main([str(committed), str(rerun), "--record", str(record)]) == 2
        captured = capsys.readouterr()
        assert "could be parsed from the record" in captured.err and "MATCH" not in captured.out

    assert main([str(committed), str(rerun), "--record", str(tmp_path / "missing.md")]) == 2
    captured = capsys.readouterr()
    assert "experiment record does not exist" in captured.err and "MATCH" not in captured.out


def test_record_primary_allows_explanatory_prose_and_multiple_valid_comparisons(tmp_path, capsys):
    committed, rerun = _run(tmp_path / "committed"), _run(tmp_path / "rerun")
    record = tmp_path / "EXP-z.md"
    record.write_text("- primary comparison: " + PRIMARY + "; ps > pm | auprc_ratio | auto\n"
                      "  (machine-readable restatement; no content change)\n", encoding="utf-8")
    assert main([str(committed), str(rerun), "--record", str(record)]) == 0
    out = capsys.readouterr().out
    assert f"primary: MATCH - {PRIMARY}" in out
    assert "primary: MATCH - ps > pm | auprc_ratio | auto" in out


def test_a_primary_whose_columns_are_missing_on_both_sides_is_not_evaluable(tmp_path, capsys):
    committed, rerun = _run(tmp_path / "committed"), _run(tmp_path / "rerun")
    code = main([str(committed), str(rerun), "--primary", "pm > ps | auprc_ratio | no_such_arm",
                 "--json", str(tmp_path / "d.json")])
    out = capsys.readouterr().out
    assert code == 1
    assert "primary: NOT EVALUABLE - pm > ps | auprc_ratio | no_such_arm" in out
    assert "committed: n=0 (too few units for the test)" in out
    assert out.rstrip().endswith("verdict: NOT EVALUABLE") and "primary: MATCH" not in out
    assert json.loads((tmp_path / "d.json").read_text())["verdict"] == "NOT EVALUABLE"


def test_a_primary_with_fewer_than_the_minimum_units_is_not_evaluable(tmp_path, capsys):
    gone = {(t, "m_ps_auto_auprc_ratio"): "error" for t in ("d3", "d4")}     # two usable units left
    committed, rerun = _run(tmp_path / "committed", change=gone), _run(tmp_path / "rerun", change=gone)
    assert main([str(committed), str(rerun), "--primary", PRIMARY]) == 1
    out = capsys.readouterr().out
    assert "primary: NOT EVALUABLE" in out and "committed n=2, rerun n=2" in out
    assert out.rstrip().endswith("verdict: NOT EVALUABLE")


def test_the_minimum_matches_compare_units():
    scores = pd.DataFrame({"a": [1.0, 2.0, 4.0], "b": [0.5, 1.0, 1.0]})
    compare_units(scores.iloc[:MIN_UNITS])                  # exactly the minimum runs
    with pytest.raises(ValueError):
        compare_units(scores.iloc[:MIN_UNITS - 1])          # one fewer is refused: keep MIN_UNITS in step


def test_differences_outrank_an_unevaluable_primary(tmp_path, capsys):
    gone = {(t, "m_ps_auto_auprc_ratio"): "error" for t in ("d3", "d4")}
    committed = _run(tmp_path / "committed", change=gone)
    rerun = _run(tmp_path / "rerun", change={**gone, ("d1", "m_pm_auto_epr"): 2.5})
    assert main([str(committed), str(rerun), "--primary", PRIMARY]) == 1
    out = capsys.readouterr().out
    assert "primary: NOT EVALUABLE" in out and out.rstrip().endswith("verdict: DIFFERENT")


@pytest.mark.parametrize("field", ["v2." + PRIMARY, "not " + PRIMARY, PRIMARY + " | two-sided",
                                   PRIMARY + "-v2", PRIMARY + ".v2"])
def test_a_malformed_record_primary_is_not_read_as_the_comparison_inside_it(tmp_path, capsys, field):
    """The record path must not accept a valid comparison found inside a malformed one."""
    committed, rerun = _run(tmp_path / "committed"), _run(tmp_path / "rerun")
    record = tmp_path / "EXP-m.md"
    record.write_text("- primary comparison: " + field + "\n", encoding="utf-8")
    assert main([str(committed), str(rerun), "--record", str(record)]) == 2
    captured = capsys.readouterr()
    assert "could be parsed from the record" in captured.err and "MATCH" not in captured.out


@pytest.mark.parametrize("field", [PRIMARY + ".", "`" + PRIMARY + "`", PRIMARY + " - restated before the run",
                                   PRIMARY + " (restated before the run)"])
def test_prose_after_a_record_primary_still_parses(tmp_path, capsys, field):
    committed, rerun = _run(tmp_path / "committed"), _run(tmp_path / "rerun")
    record = tmp_path / "EXP-p.md"
    record.write_text("- primary comparison: " + field + "\n", encoding="utf-8")
    assert main([str(committed), str(rerun), "--record", str(record)]) == 0
    assert f"primary: MATCH - {PRIMARY}" in capsys.readouterr().out


@pytest.mark.parametrize("flag,value", [("--rtol", "-1"), ("--atol", "-0.1"), ("--atol", "inf"), ("--rtol", "nan")])
def test_a_negative_or_non_finite_tolerance_is_refused(tmp_path, flag, value):
    committed, rerun = _run(tmp_path / "committed"), _run(tmp_path / "rerun")
    with pytest.raises(SystemExit) as exc:
        main([str(committed), str(rerun), "--primary", PRIMARY, flag, value])
    assert exc.value.code == 2


def test_a_primary_naming_the_same_model_twice_exits_2(tmp_path, capsys):
    committed, rerun = _run(tmp_path / "committed"), _run(tmp_path / "rerun")
    assert main([str(committed), str(rerun), "--primary", "pm > pm | auprc_ratio | auto"]) == 2
    captured = capsys.readouterr()
    assert "same model twice" in captured.err and "MATCH" not in captured.out


@pytest.mark.parametrize("content", ["null", "{not json", "[1, 2]"])
def test_an_unreadable_row_file_exits_2_and_names_it(tmp_path, capsys, content):
    committed, rerun = _run(tmp_path / "committed"), _run(tmp_path / "rerun")
    (rerun / "rows" / "d2.json").write_text(content, encoding="utf-8")
    assert main([str(committed), str(rerun), "--primary", PRIMARY]) == 2
    captured = capsys.readouterr()
    assert "unreadable row file" in captured.err and "d2.json" in captured.err and "MATCH" not in captured.out
