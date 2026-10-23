"""benchmarks/compare_confirm.py: judge a candidate-study rerun's primary lines from CONFIRM.json."""
import copy
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("compare_confirm", ROOT / "benchmarks" / "compare_confirm.py")
cc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cc)

PRIMARY = ("- primary comparison: baseline > noref_shipped | eval_mean | clean; "
           "ref_geneset_preserving > baseline | eval_mean | clean\n")


def _pair(a, b, diff, p, supported):
    return {"a": a, "b": b, "n": 20, "excluded": 0, "diff": diff, "ci": [diff - 0.07, diff + 0.07],
            "p_a_gt_b": p, "p_b_gt_a": 1.0, "supported": supported}


CONFIRM = {
    "clean": [_pair("baseline", "noref_shipped", 0.08898038017836948, 0.04347084550116562, "baseline"),
              _pair("ref_geneset_preserving", "baseline", 0.0226, 0.201, "unresolved")],
    "doublets_20": [_pair("baseline", "noref_shipped", 0.1119, 0.0054, "baseline"),
                    _pair("ref_geneset_preserving", "baseline", 0.024, 0.0061, "ref_geneset_preserving")],
}


def _run(tmp_path, rerun, *extra, record_text=PRIMARY):
    a, b, rec = tmp_path / "a.json", tmp_path / "b.json", tmp_path / "rec.md"
    a.write_text(json.dumps(CONFIRM), encoding="utf-8")
    b.write_text(json.dumps(rerun), encoding="utf-8")
    rec.write_text("# EXP\n\n" + record_text, encoding="utf-8")
    return cc.main([str(a), str(b), "--record", str(rec), *extra])


def test_identical_files_match(tmp_path, capsys):
    assert _run(tmp_path, CONFIRM) == 0
    out = capsys.readouterr().out
    assert out.count("primary: MATCH") == 2
    assert out.strip().endswith("verdict: MATCH")


def test_last_digit_drift_is_different_exactly_but_matches_within_rtol(tmp_path, capsys):
    rerun = copy.deepcopy(CONFIRM)
    rerun["clean"][0]["diff"] *= 1 + 1e-12
    assert _run(tmp_path, rerun) == 1
    assert "primary: DIFFERENT - baseline > noref_shipped | eval_mean | clean" in capsys.readouterr().out
    assert _run(tmp_path, rerun, "--rtol", "1e-9") == 0


def test_a_changed_decision_is_different_whatever_the_tolerance(tmp_path, capsys):
    rerun = copy.deepcopy(CONFIRM)
    rerun["clean"][1]["supported"] = "ref_geneset_preserving"
    assert _run(tmp_path, rerun, "--rtol", "1") == 1
    out = capsys.readouterr().out
    assert "supported: 'unresolved' -> 'ref_geneset_preserving'" in out


def test_interval_bound_is_compared(tmp_path):
    rerun = copy.deepcopy(CONFIRM)
    rerun["clean"][0]["ci"][1] += 1e-6
    assert _run(tmp_path, rerun, "--rtol", "1e-9") == 1


def test_secondary_conditions_are_information_only(tmp_path, capsys):
    rerun = copy.deepcopy(CONFIRM)
    rerun["doublets_20"][0]["diff"] += 0.01
    assert _run(tmp_path, rerun) == 0
    out = capsys.readouterr().out
    assert "secondary (information only): DIFFERENT - baseline > noref_shipped | doublets_20" in out


def test_missing_primary_pair_is_an_error_not_a_pass(tmp_path):
    rerun = copy.deepcopy(CONFIRM)
    rerun["clean"] = rerun["clean"][:1]
    assert _run(tmp_path, rerun) == 2


def test_record_without_primary_line_is_an_error(tmp_path):
    assert _run(tmp_path, CONFIRM, record_text="- type: confirmatory\n") == 2


def test_parses_the_registered_record():
    rec = ROOT / "experiments" / "EXP-20261003-LM-1.md"
    assert cc.read_primaries(rec) == [("baseline", "noref_shipped", "eval_mean", "clean"),
                                      ("ref_geneset_preserving", "baseline", "eval_mean", "clean")]


def test_committed_file_against_itself_matches(capsys):
    committed = ROOT / "benchmarks" / "results" / "EXP-20261003-LM-1" / "CONFIRM.json"
    if not committed.exists():
        pytest.skip("committed CONFIRM.json not present")
    rec = ROOT / "experiments" / "EXP-20261003-LM-1.md"
    assert cc.main([str(committed), str(committed), "--record", str(rec)]) == 0
    assert capsys.readouterr().out.count("primary: MATCH") == 2
