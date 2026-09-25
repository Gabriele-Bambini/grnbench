"""Reproduce audit findings; these are audit fixtures, not repair/regression tests.

Run: python audit/2026-09-13/reproduce_eval_findings.py
Loads the real evaluation modules without executing optional scanpy-dependent
package initializers. Requires only numpy, pandas, and scipy. No source changes.
"""
from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path
import sys
import types

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def load_eval_modules():
    for name, path in (("grnbench", ROOT / "src/grnbench"),
                       ("grnbench.eval", ROOT / "src/grnbench/eval")):
        package = types.ModuleType(name)
        package.__path__ = [str(path)]
        sys.modules[name] = package
    loaded = {}
    for short in ("metrics", "nulls", "significance"):
        name = f"grnbench.eval.{short}"
        spec = importlib.util.spec_from_file_location(
            name, ROOT / "src/grnbench/eval" / f"{short}.py")
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
        loaded[short] = module
    return loaded


def edges(rows):
    return pd.DataFrame(rows, columns=["source", "target", "weight"])


def clean_json(value):
    if isinstance(value, dict):
        return {k: clean_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean_json(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def main():
    modules = load_eval_modules()
    metrics, nulls, significance = (modules[n] for n in
                                    ("metrics", "nulls", "significance"))
    findings = []

    # 1. Same candidate universe, but the default resampling changes EPR's k.
    gold = [("T", f"P{i}") for i in range(100)]
    negative = [("T", f"N{i}") for i in range(900)]

    def ranking(positive_positions):
        positions = set(positive_positions)
        positives, negatives = iter(gold), iter(negative)
        return edges([(next(positives) if i in positions else next(negatives))
                      + (1000 - i,) for i in range(1000)])

    a = ranking(list(range(50, 90)) + list(range(200, 260)))
    b = ranking(list(range(25)) + list(range(300, 375)))
    default = significance.paired_bootstrap(
        {"A": a, "B": b}, set(gold), metric="epr", n_boot=200,
        frac=0.5, seed=42)
    fixed = significance.paired_bootstrap(
        {"A": a, "B": b}, set(gold), metric="epr", n_boot=200,
        frac=0.5, seed=42, top_k=100)
    findings.append({
        "case": "default_cutoff_inversion",
        "input": {"n_candidates": 1000, "n_gold": 100,
                  "A_positive_rank_ranges_inclusive": [[51, 90], [201, 260]],
                  "B_positive_rank_ranges_inclusive": [[1, 25], [301, 375]]},
        "default": default, "fixed_k_100": fixed,
        "finding_reproduced": bool(
            default["observed"]["A"] > default["observed"]["B"]
            and default["pairs"]["A|B"]["diff_ci"][1] < 0
            and fixed["pairs"]["A|B"]["diff_ci"][0] > 0),
        "interpretation": "Default resampling changes k from 100 to 50, reversing the estimand's ordering.",
    })

    # 2. Marginal row counts survive, but unique directed degrees do not.
    original = edges([("A", "x", 4), ("A", "y", 3),
                      ("B", "x", 2), ("B", "y", 1)])
    shuffled = nulls.degree_preserving_nulls(original, 1, seed=0)[0]
    findings.append({
        "case": "degree_null_duplicates",
        "original": original.to_dict("records"),
        "null": shuffled.to_dict("records"),
        "original_unique_outdegree": original.groupby("source").target.nunique().to_dict(),
        "null_unique_outdegree": shuffled.groupby("source").target.nunique().to_dict(),
        "finding_reproduced": bool(shuffled.duplicated(["source", "target"]).any()),
        "interpretation": "Target permutation preserves multigraph row margins, not simple-graph degrees.",
    })

    # 3. Self-loops are introduced despite being absent from the candidate set.
    original = edges([("A", "B", 6), ("A", "C", 5),
                      ("B", "A", 4), ("B", "C", 3),
                      ("C", "A", 2), ("C", "B", 1)])
    shuffled = nulls.degree_preserving_nulls(original, 1, seed=0)[0]
    random = nulls.random_null(original, 1, seed=0)[0]
    loops = int((shuffled.source == shuffled.target).sum())
    findings.append({
        "case": "null_support_and_random_margins",
        "original": original.to_dict("records"),
        "degree_null": shuffled.to_dict("records"),
        "n_null_self_loops": loops,
        "random_null_source_row_counts_preserved": (
            original.source.value_counts().to_dict() == random.source.value_counts().to_dict()),
        "random_null_target_row_counts_preserved": (
            original.target.value_counts().to_dict() == random.target.value_counts().to_dict()),
        "finding_reproduced": loops > 0,
        "interpretation": "Degree null leaves the original loop-free support; random mode also preserves both marginal row counts.",
    })

    # 4. A leader-only tier comparison misses significant differences inside it.
    pairs = {"A|B": {"significant": False, "diff_mean": 1.0},
             "A|C": {"significant": False, "diff_mean": 2.0},
             "B|C": {"significant": True, "diff_mean": 1.0}}
    tiers = significance.significant_tiers(
        ["A", "B", "C"], {"A": 3.0, "B": 2.0, "C": 1.0}, pairs)
    findings.append({
        "case": "tier_nontransitivity", "input_pairs": pairs, "tiers": tiers,
        "finding_reproduced": any("B" in tier and "C" in tier for tier in tiers),
        "interpretation": "The same tier can contain a pair explicitly declared significantly different.",
    })

    # 5. The denominator is only the positives surviving in each method's list.
    original = edges([("T", f"G{i}", 100 - i) for i in range(100)])
    gold = {("T", "G0"), ("T", "G99")}
    full = metrics.evaluate_ranking(original, gold, with_curve=False)
    truncated = metrics.evaluate_ranking(original.head(1), gold, with_curve=False)
    findings.append({
        "case": "candidate_list_truncation", "full": full,
        "truncated_to_first_edge": truncated,
        "finding_reproduced": full["auprc"] < truncated["auprc"] == 1.0,
        "interpretation": "Truncation silently removes a missed positive from recall's denominator; cross-list AP is not a common-universe comparison.",
    })

    # 6. Exact finite-sample AP expectation differs from prevalence.
    scores = []
    for positive_index in range(4):
        labels = np.zeros(4, dtype=int)
        labels[positive_index] = 1
        scores.append(metrics.average_precision(labels))
    expectation = float(np.mean(scores))
    findings.append({
        "case": "finite_random_ap_baseline", "n_edges": 4, "n_positives": 1,
        "all_equally_likely_random_APs": scores,
        "exact_expected_random_AP": expectation, "prevalence": 0.25,
        "finding_reproduced": not np.isclose(expectation, 0.25),
        "interpretation": "Prevalence is not the exact finite-sample expectation of random-ranking AP.",
    })

    result = clean_json({
        "artifact_type": "read-only source audit fixtures, not repair tests",
        "source_loading": "importlib; package initializers bypassed, real evaluation modules executed",
        "seed": 42,
        "findings": findings,
        "all_six_findings_reproduced": all(x["finding_reproduced"] for x in findings),
    })
    output = Path(__file__).with_name("eval_findings.json")
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(output),
                      "all_six_findings_reproduced": result["all_six_findings_reproduced"]}))
    return 0 if result["all_six_findings_reproduced"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
