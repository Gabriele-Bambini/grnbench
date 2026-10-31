"""random_baseline: the negative control of new benchmark runs (D-037).

It scores the candidate universe the other models score with a seeded uniform weight, so its AUPRC
ratio and EPR are an empirical null (about 1 with many evaluation positives); it is recorded in the task lists that new plans write (old lists keep
their exact behaviour) and reported apart, never in a comparison family.
"""
import inspect
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))

import candidate_study  # noqa: E402
import real_bench  # noqa: E402
from grnbench import seeds  # noqa: E402
from grnbench.infer import base as infer_base  # noqa: E402
from grnbench.infer import run_method  # noqa: E402


def _adata(n_cells=200, n_genes=300, seed=3):
    rng = np.random.default_rng(seed)
    genes = [f"TF{i}" for i in range(20)] + [f"G{i}" for i in range(n_genes - 20)]
    X = rng.poisson(2.0, size=(n_cells, n_genes)).astype(np.float32)
    return ad.AnnData(X, var=pd.DataFrame(index=genes)), [f"TF{i}" for i in range(20)]


def test_registered_as_a_control_outside_the_default_set():
    spec = infer_base.REGISTRY["random_baseline"]
    assert spec.origin == "builtin" and spec.default is False
    assert real_bench.NEGATIVE_CONTROLS == ("random_baseline",)


def test_seed_is_fixed_in_seeds_and_is_the_runner_default():
    from grnbench.infer.random_baseline import random_baseline
    assert isinstance(seeds.RANDOM_BASELINE_SEED, int)
    assert inspect.signature(random_baseline).parameters["seed"].default == seeds.RANDOM_BASELINE_SEED
    assert seeds.LEGACY["library"]["model.random_baseline"] == seeds.RANDOM_BASELINE_SEED


def test_same_universe_as_the_correlation_baseline_and_deterministic():
    adata, tfs = _adata()
    rand, _ = run_method("random_baseline", adata, tf_list=tfs)
    pear, _ = run_method("pearson_signed", adata, tf_list=tfs)
    pairs = lambda e: set(zip(e["source"], e["target"]))  # noqa: E731
    assert pairs(rand) == pairs(pear)
    again, _ = run_method("random_baseline", adata, tf_list=tfs)
    pd.testing.assert_frame_equal(rand, again)
    other, _ = run_method("random_baseline", adata, tf_list=tfs, seed=seeds.RANDOM_BASELINE_SEED + 1)
    assert not np.array_equal(rand["weight"].to_numpy(), other["weight"].to_numpy())
    w = rand["weight"].to_numpy()
    assert ((w >= 0) & (w < 1)).all()


def test_auprc_ratio_and_epr_are_about_one():
    from grnbench.eval.metrics import GoldIndex, evaluate_ranking
    adata, tfs = _adata(n_genes=520)
    edges, _ = run_method("random_baseline", adata, tf_list=tfs)
    universe = edges[["source", "target"]]
    gold = universe.sample(n=600, random_state=11)            # 600 of ~10,000 candidate pairs
    r = evaluate_ranking(edges, GoldIndex.coerce(set(zip(gold["source"], gold["target"]))),
                         with_curve=False, universe=universe)
    assert 0.8 < r["auprc_ratio"] < 1.25, r["auprc_ratio"]
    assert 0.7 < r["epr"] < 1.35, r["epr"]


def test_never_in_a_comparison_family():
    members = {m for pair in real_bench.DECLARED_PAIRS for m in pair}
    assert not members & set(real_bench.NEGATIVE_CONTROLS)
    arms = {m for pair in candidate_study.CONFIRM_FAMILY for m in pair}
    assert not arms & set(real_bench.NEGATIVE_CONTROLS)
    assert "random_baseline" not in real_bench.EVAL_METHODS          # the screening proxies stay as they were


def test_new_plans_record_the_control_old_task_lists_do_not(tmp_path):
    root = tmp_path / "data"
    d = root / "cellxgene" / "study" / "ds"
    d.mkdir(parents=True)
    (d / "truth.json").write_text(json.dumps({"organism": "human"}), encoding="utf-8")
    tasks = real_bench.make_tasks(root)
    assert tasks and all(t["negative_controls"] == ["random_baseline"] for t in tasks)
    study = candidate_study.make_tasks(root, conditions=("clean",))
    assert study and all(t["negative_controls"] == ["random_baseline"] for t in study)
    assert real_bench.task_controls(tasks[0]) == ["random_baseline"]
    legacy = {k: v for k, v in tasks[0].items() if k != "negative_controls"}
    assert real_bench.task_controls(legacy) == []


def _rows(with_control: bool) -> pd.DataFrame:
    rows = []
    for i in range(4):
        r = {"task": f"cxg__d{i}", "m_pearson_signed_auto_auprc_ratio": 1.5 + i / 10,
             "m_pearson_signed_baseline_auprc_ratio": 1.4, "m_pearson_signed_auto_epr": 1.3,
             "m_pearson_signed_baseline_epr": 1.2, "m_pearson_signed_auto_seconds": 1.0,
             "m_clr_auto_auprc_ratio": 1.2, "m_clr_baseline_auprc_ratio": 1.1, "m_clr_auto_epr": 1.1,
             "m_clr_baseline_epr": 1.0, "m_clr_auto_seconds": 1.0}
        if with_control:
            r.update({"nc_random_baseline_auto_auprc_ratio": 0.98 + i / 100,
                      "nc_random_baseline_baseline_auprc_ratio": 1.01,
                      "nc_random_baseline_auto_epr": 0.97, "nc_random_baseline_baseline_epr": 1.02,
                      "nc_random_baseline_auto_seconds": 0.1,
                      "nc_random_baseline_auto_n_pos_min": 25 + i, "nc_random_baseline_auto_n_pos_max": 300 + i,
                      "nc_random_baseline_baseline_n_pos_min": 31, "nc_random_baseline_baseline_n_pos_max": 412})
        rows.append(r)
    return pd.DataFrame(rows)


def test_report_shows_the_control_apart_and_only_when_present():
    plain = real_bench.models_section(_rows(False))
    assert not any("egative control" in line for line in plain)
    lines = real_bench.models_section(_rows(True))
    text = "\n".join(lines)
    assert "Negative control (random)" in text
    head = text.split("Negative control (random)")[0]
    assert "random_baseline" not in head                         # not in the model tables or comparisons
    assert lines[:len(plain)] == plain                           # everything before it is unchanged
    assert "| random_baseline | 4 |" in text


def test_report_prints_the_reading_rule_and_the_evaluation_positives():
    text = "\n".join(real_bench.models_section(_rows(True)))
    control = text.split("Negative control (random)")[1]
    assert "not a test against 1" in control and "number of evaluation positives" in control
    assert "one fixed draw per universe" in control
    assert "positives auto (min-max) | positives baseline (min-max) |" in control
    # auto: min over datasets of the per-dataset min (25), max of the per-dataset max (303); baseline 31-412
    assert "| 0.995 | 0.995 | 1.010 | 0.970 | 1.020 | 25-303 | 31-412 |" in control


def test_score_models_writes_the_control_under_its_own_prefix(monkeypatch):
    seen = {}

    def fake(records, plans, arms, methods, tfs, pairs, artifacts_dir=None, support=False):
        seen["support"] = support
        return {f"m_{m}_{a}_auprc_ratio": 1.0 for m in methods for a in arms}
    monkeypatch.setattr(real_bench, "score_models", fake)
    out = real_bench.score_controls({}, {}, {"baseline": "baseline", "auto": "x"}, ["random_baseline"], None, set())
    assert set(out) == {"nc_random_baseline_baseline_auprc_ratio", "nc_random_baseline_auto_auprc_ratio"}
    assert seen["support"] is True                                 # the control carries its positives


class _Plan:
    def value(self, key):
        return "raw_counts"


def _processed(tmp_path, n_genes=120):
    adata, tfs = _adata(n_genes=n_genes)
    path = tmp_path / "processed.h5ad"
    adata.write_h5ad(path)
    genes = list(adata.var_names)
    gold = {(tfs[i % len(tfs)], genes[20 + 3 * i]) for i in range(30)}
    return {"x": [{"processed": str(path), "context": "all"}]}, {"x": _Plan()}, tfs, gold


def test_scorers_count_the_evaluation_positives_only_when_asked(tmp_path):
    import engine_bench
    records, plans, tfs, gold = _processed(tmp_path)
    plain = engine_bench.score_candidates(records, plans, tfs, gold, methods=["random_baseline"])["x"]
    assert not any(k.endswith(("_n_pos_min", "_n_pos_max")) for k in plain)
    sc = engine_bench.score_candidates(records, plans, tfs, gold, methods=["random_baseline"], support=True)["x"]
    assert sc["eval_random_baseline_n_pos_min"] == sc["eval_random_baseline_n_pos_max"] == 30
    assert sc["eval_random_baseline"] == plain["eval_random_baseline"]
    arms = {"baseline": "x", "auto": "x"}
    models = real_bench.score_models(records, plans, arms, ["pearson_signed"], tfs, gold)
    assert not any("n_pos" in k for k in models)                  # model rows gain no key
    nc = real_bench.score_controls(records, plans, arms, ["random_baseline"], tfs, gold)
    assert nc["nc_random_baseline_auto_n_pos_min"] == nc["nc_random_baseline_baseline_n_pos_max"] == 30
    assert nc["nc_random_baseline_auto_auprc_ratio"] == nc["nc_random_baseline_baseline_auprc_ratio"]


def test_draw_rule_is_what_the_code_does():
    from grnbench.infer.random_baseline import DRAW_RULE
    adata, tfs = _adata(n_genes=60)
    edges, _ = run_method("random_baseline", adata, tf_list=tfs[::-1])      # TF-list order does not matter
    genes = list(adata.var_names)
    regs = [g for g in genes if g in set(tfs)]                              # expression-column order
    w = np.random.default_rng(seeds.RANDOM_BASELINE_SEED).random(len(regs) * len(genes))
    want = pd.DataFrame({"source": np.repeat(regs, len(genes)), "target": genes * len(regs), "weight": w})
    want = want[want["source"] != want["target"]]
    key = ["source", "target"]
    got = edges[key + ["weight"]].sort_values(key).reset_index(drop=True)
    want = want.sort_values(key).reset_index(drop=True)
    pd.testing.assert_frame_equal(got, want, check_dtype=False)
    assert "default_rng(seed)" in DRAW_RULE and "regulator-major" in DRAW_RULE


def test_run_config_records_the_control_seed_and_draw_only_for_new_task_lists():
    from grnbench.infer.random_baseline import DRAW_RULE
    old = real_bench.scientific_config(["pearson_signed"], [], None)
    assert set(old) == {"models", "screen_methods", "max_tfs", "max_score_cells", "cell_sampling",
                        "min_reference_edges", "reference_split", "reference_sha256"}
    new = real_bench.scientific_config(["pearson_signed"], ["random_baseline"], None)
    assert {k: v for k, v in new.items() if not k.startswith("negative_control")} == old
    assert new["negative_controls"] == ["random_baseline"]
    assert new["negative_control_draws"] == {"random_baseline": {"seed": 20261027, "draw": DRAW_RULE}}
    assert seeds.RANDOM_BASELINE_SEED == 20261027
    json.dumps(new)                                                          # RUN_CONFIG.json must serialise


def test_candidate_study_summary_shows_the_control_only_when_present():
    base = {"task": "d__clean", "dataset": "d", "condition": "clean", "seed": 0, "error": None,
            "candidates": [{"name": "baseline", "source": "baseline", "eval_mean": 1.2}]}
    plain = candidate_study.summary([dict(base)])
    assert "egative control" not in plain
    row = dict(base, negative_control={"candidate": "baseline", "eval_random_baseline": 1.01,
                                       "eval_random_baseline_epr": 0.99, "eval_random_baseline_n_pos_min": 12,
                                       "eval_random_baseline_n_pos_max": 640,
                                       "draws": dict(real_bench.NEGATIVE_CONTROL_DRAWS)})
    text = candidate_study.summary([row])
    assert "Negative control (random)" in text and text.startswith(plain.rstrip("\n"))
    control = text.split("Negative control (random)")[1]
    assert "not a test against 1" in control and "one fixed draw per universe" in control
    assert "| random_baseline | clean | 1 | 1.010 | 0.990 | 12-640 |" in control
    assert "n_pos" not in control                                   # the counts are columns, not phantom controls
