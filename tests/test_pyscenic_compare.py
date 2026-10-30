from __future__ import annotations

import json
import sys
from argparse import Namespace
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "benchmarks"))
import pyscenic_compare as pc  # noqa: E402


def write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload), encoding="utf-8")


def test_regulon_parser_strips_direction_and_self_pairs(tmp_path):
    path = tmp_path / "regulons.json"
    write_json(path, {"TF1(+ )": ["G1"], "TF1(+ ) ": ["G2"],
                      "TF2(-)": ["TF2", "G3"]})
    # The official names are '(+)'/'(-)'; surrounding whitespace is accepted.
    write_json(path, {"TF1 (+)": ["G1", "TF1"], "TF2(-)": ["G2"]})
    assert pc.load_regulon_pairs(path) == {("TF1", "G1"), ("TF2", "G2")}


def test_set_metrics_and_module_union_have_defined_empty_case():
    modules = [{"tf": "TF", "genes": ["TF", "A", "B"]},
               {"tf": "TF", "genes": ["B", "C"]}]
    assert pc.pair_union(modules) == {("TF", "A"), ("TF", "B"), ("TF", "C")}
    assert pc.set_metrics(set(), set()) == {
        "candidate": 0, "reference": 0, "intersection": 0, "union": 0,
        "exact": True, "jaccard": 1.0, "precision": 1.0, "recall": 1.0,
    }
    metric = pc.set_metrics({1, 2}, {2, 3})
    assert metric["jaccard"] == pytest.approx(1 / 3)
    assert metric["precision"] == metric["recall"] == 0.5


def test_rerank_preserves_full_support_and_base_order():
    edges = pd.DataFrame({"source": ["A", "A", "B", "B"],
                          "target": ["x", "y", "x", "y"],
                          "weight": [4.0, 3.0, 2.0, 1.0]})
    got = pc.rerank_by_support(edges, {("A", "y"), ("B", "y")})
    assert list(zip(got.source, got.target)) == [
        ("A", "y"), ("B", "y"), ("A", "x"), ("B", "x")]
    assert set(zip(got.source, got.target)) == set(zip(edges.source, edges.target))


def _portable_fixture(tmp_path: Path):
    units = tmp_path / "units.tsv"
    units.write_text("task\tarm\tcontext\tspecies\n"
                     "human_task\tbaseline\tctx0\thuman\n"
                     "mouse_task\tbaseline\tctx0\tmouse\n", encoding="utf-8")
    prepared = tmp_path / "prepared"
    records = []
    for task, species in (("human_task", "human"), ("mouse_task", "mouse")):
        uid = f"{task}__baseline__ctx0"
        records.append({"id": uid, "task": task, "arm": "baseline", "context": "ctx0",
                        "species": species})
    write_json(prepared / "PREPARED.json", {"status": "prepared", "units": records})
    complete = prepared / "units" / records[0]["id"]
    write_json(complete / "run.json", {"status": "completed", "seconds": 2.0})
    write_json(complete / "grnbench_modules.json",
               [{"tf": "TF1", "genes": ["G1", "G2"]}])
    write_json(complete / "pyscenic_modules.json",
               [{"tf": "TF1", "genes": ["TF1", "G1", "G3"]}])
    pd.DataFrame([("TF1", "G1")], columns=["source", "target"]).to_csv(
        complete / "grnbench_support.tsv", sep="\t", index=False)
    write_json(complete / "pyscenic_regulons.json", {"TF1 (+)": ["G1", "G2"]})
    pd.DataFrame([("TF1", "G1", 3.0), ("TF1", "G2", 2.0), ("TF1", "G3", 1.0)],
                 columns=["TF", "target", "importance"]).to_csv(
        complete / "adjacencies.tsv", sep="\t", index=False)
    pd.DataFrame([("TF1", "G2")], columns=["source", "target"]).to_csv(
        complete / "reference_evaluation.tsv", sep="\t", index=False)
    failed = prepared / "units" / records[1]["id"]
    write_json(failed / "run.json", {"status": "permanent_error", "seconds": 1.0,
                                      "error": "ImportError: reference package failed"})
    return units, prepared


def test_report_accounts_for_completed_and_permanent_error(tmp_path):
    units, prepared = _portable_fixture(tmp_path)
    output = tmp_path / "portable"
    assert pc.report(Namespace(prepared=str(prepared), units=str(units), out=str(output))) == 0
    result = json.loads((output / "COMPARISON.json").read_text(encoding="utf-8"))
    assert (result["expected"], result["terminal"], result["completed"],
            result["permanent_errors"]) == (2, 2, 1, 1)
    assert result["units"][0]["modules"]["jaccard"] == pytest.approx(1 / 3)
    assert result["units"][0]["support"]["jaccard"] == pytest.approx(0.5)
    assert (output / "COMPLETE.json").is_file()
    assert (output / "viz" / "implementation_overlap.png").is_file()
    assert (output / "viz" / "reference_recovery.png").is_file()
    assert (output / "viz" / "index.md").is_file()
    assert (output / "latex" / "report.tex").is_file()
    assert "NaN" not in (output / "COMPARISON.json").read_text(encoding="utf-8")


def test_report_refuses_silent_missing_unit(tmp_path):
    units, prepared = _portable_fixture(tmp_path)
    missing = prepared / "units" / "mouse_task__baseline__ctx0" / "run.json"
    missing.unlink()
    with pytest.raises(RuntimeError, match="not terminal"):
        pc.report(Namespace(prepared=str(prepared), units=str(units), out=str(tmp_path / "out")))


def test_recreated_top_edges_must_match_kept_artifact(tmp_path):
    full = pd.DataFrame({"source": ["A", "A"], "target": ["x", "y"],
                         "weight": [2.0, 1.0]})
    kept = tmp_path / "grnboost2.tsv.gz"
    full.iloc[:1].to_csv(kept, sep="\t", index=False, compression="gzip")
    pc._same_kept_edges(full, kept)
    changed = full.iloc[:1].assign(target="z")
    changed.to_csv(kept, sep="\t", index=False, compression="gzip")
    with pytest.raises(RuntimeError, match="pairs"):
        pc._same_kept_edges(full, kept)
