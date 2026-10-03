"""Task manifests, completeness checks and the controlled stress design."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest


def _engine_bench():
    path = Path(__file__).parents[1] / "benchmarks" / "engine_bench.py"
    spec = importlib.util.spec_from_file_location("engine_bench_runtime_tests", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_controlled_stress_grid_has_every_declared_curve_and_740_tasks():
    bench = _engine_bench()
    tasks = bench.make_stress_tasks()
    assert len(tasks) == 740
    assert len({t["id"] for t in tasks}) == 740
    assert sum(bench.is_sentinel_task(t) for t in tasks) == 12
    assert {t["shared_edge_fraction"] for t in tasks} == {1.0, 0.75, 0.5, 0.25}
    for shared in (1.0, 0.75, 0.5, 0.25):
        regime = [t for t in tasks if t["shared_edge_fraction"] == shared]
        assert len(regime) == 37 * 5
        assert {t["stressor"] for t in regime} == {"clean", "mixed", *bench.STRESS_FACTORS}
        assert sum(t["stressor"] == "clean" for t in regime) == 5


def test_task_manifest_is_hashed_and_rejects_mutation(tmp_path):
    bench = _engine_bench()
    path = tmp_path / "tasks.json"
    tasks = [{"id": "a", "seed": 0}, {"id": "b", "seed": 1}]
    written = bench.write_task_manifest(path, tasks, task_type="test")
    loaded, meta = bench.load_task_manifest(path)
    assert loaded == tasks and meta["tasks_sha256"] == written["tasks_sha256"]
    payload = json.loads(path.read_text())
    payload["tasks"][0]["seed"] = 99
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="hash mismatch"):
        bench.load_task_manifest(path)


def test_context_truth_is_macro_scored_and_shared_edges_cannot_leak_between_halves():
    bench = _engine_bench()
    from grnbench.screen import split_reference

    by_context = {"type1": {("TF", "shared"), ("TF", "only1")},
                  "type2": {("TF", "shared"), ("TF", "only2")}}
    union = set().union(*by_context.values())
    selection, evaluation = split_reference(union, seed=4)
    select_by = {k: v & selection for k, v in by_context.items()}
    eval_by = {k: v & evaluation for k, v in by_context.items()}
    assert selection.isdisjoint(evaluation)
    assert all(select_by[k].isdisjoint(eval_by[k]) for k in by_context)
    pooled = bench._truth_groups({"subset_value": None}, evaluation, eval_by)
    one_type = bench._truth_groups({"subset_value": "type1"}, evaluation, eval_by)
    assert set(pooled) == {"type1", "type2"}
    assert one_type == {"type1": eval_by["type1"]}


def test_run_configuration_is_immutable(tmp_path):
    bench = _engine_bench()
    meta = {"tasks_sha256": "abc"}
    first = bench.ensure_run_config(tmp_path, meta, runner="engine_bench",
                                    scientific_config={"models": ["pearson_signed"]})
    assert bench.ensure_run_config(tmp_path, meta, runner="engine_bench",
                                   scientific_config={"models": ["pearson_signed"]}) == first
    with pytest.raises(ValueError, match="configuration differs"):
        bench.ensure_run_config(tmp_path, meta, runner="engine_bench",
                                scientific_config={"models": ["clr"]})


def test_completeness_refuses_missing_or_retryable_rows_and_allows_watermark(tmp_path):
    bench = _engine_bench()
    tasks_path = tmp_path / "tasks.json"
    tasks = [{"id": "done"}, {"id": "retry"}, {"id": "missing"}]
    meta = bench.write_task_manifest(tasks_path, tasks, task_type="test")
    bench.ensure_run_config(tmp_path, meta, runner="test", scientific_config={})
    rows = tmp_path / "rows"
    rows.mkdir()
    (rows / "done.json").write_text(json.dumps({"task": "done", "status": "completed"}))
    (rows / "retry.json").write_text(json.dumps(
        {"task": "retry", "status": "transient_error", "error": "network"}))
    with pytest.raises(RuntimeError, match="missing=1, retryable=1"):
        bench.enforce_completeness(tmp_path, tasks_path)
    audit = bench.enforce_completeness(tmp_path, tasks_path, partial=True)
    assert not audit["complete"]
    assert (tmp_path / "PARTIAL.json").exists()
    assert not (tmp_path / "COMPLETE.json").exists()


def test_permanent_error_is_terminal_but_remains_visible(tmp_path):
    bench = _engine_bench()
    tasks_path = tmp_path / "tasks.json"
    meta = bench.write_task_manifest(tasks_path, [{"id": "bad"}], task_type="test")
    bench.ensure_run_config(tmp_path, meta, runner="test", scientific_config={})
    rows = tmp_path / "rows"
    rows.mkdir()
    (rows / "bad.json").write_text(json.dumps(
        {"task": "bad", "status": "permanent_error", "error": "deterministic"}))
    audit = bench.enforce_completeness(tmp_path, tasks_path)
    assert audit["complete"] and audit["unexpected_errors"] == ["bad"]
    assert json.loads((tmp_path / "COMPLETE.json").read_text())["unexpected_errors"] == ["bad"]


def test_stress_breaking_points_write_generated_figures(tmp_path):
    import numpy as np
    import pandas as pd

    bench = _engine_bench()
    rows = []
    for shared in (1.0, 0.5):
        for seed in range(5):
            rows.append({"shared_edge_fraction": shared, "stressor": "clean", "severity": 0,
                         "severity_label": "clean", "severity_value": 0,
                         "m_pearson_signed_auto_auprc_ratio": 2.0 + seed / 100})
            for severity, label in ((1, "mild"), (2, "extreme")):
                rows.append({"shared_edge_fraction": shared, "stressor": "depth",
                             "severity": severity, "severity_label": label,
                             "severity_value": severity,
                             "m_pearson_signed_auto_auprc_ratio":
                                 np.nan if severity == 2 and seed == 0 else 2.0 - severity / 2})

    lines = bench.stress_breaking_points(pd.DataFrame(rows), tmp_path)

    assert any("Controlled stress-test breaking points" in line for line in lines)
    assert (tmp_path / "stress_levels.csv").exists()
    assert (tmp_path / "breaking_points.csv").exists()
    assert (tmp_path / "viz" / "stress_depth.png").exists()
    assert "stress_depth.png" in (tmp_path / "viz" / "index.md").read_text()


def test_report_rebuilds_committed_lm1_run_without_tasks_flag(tmp_path):
    """Runs planned before D-019 saved tasks.json, not real_tasks.json: `real_bench.py report
    --out <RUN>` must still regenerate them, and values.tex must not move (AGENTS.md section 7)."""
    import shutil
    import sys

    from grnbench import texreport

    root = Path(__file__).parents[1]
    committed = root / "benchmarks" / "results" / "EXP-20260926-LM-1"
    run = tmp_path / "EXP-20260926-LM-1"
    run.mkdir()
    shutil.copytree(committed / "rows", run / "rows")
    shutil.copy(committed / "tasks.json", run / "tasks.json")
    sys.path.insert(0, str(root / "benchmarks"))
    spec = importlib.util.spec_from_file_location("real_bench_runtime_tests",
                                                  root / "benchmarks" / "real_bench.py")
    real_bench = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(real_bench)

    real_bench.report(run)
    audit = json.loads((run / "COMPLETE.json").read_text(encoding="utf-8"))
    assert audit["complete"] and audit["task_count"] == 20
    texreport.build_report(str(run), exp="EXP-20260926-LM-1", out=run / "latex", compile_pdf=False)
    assert ((run / "latex" / "values.tex").read_text(encoding="utf-8")
            == (committed / "latex" / "values.tex").read_text(encoding="utf-8"))


def test_real_report_outcome_does_not_call_scoreless_rows_scored():
    import sys
    import numpy as np
    import pandas as pd

    root = Path(__file__).parents[1]
    sys.path.insert(0, str(root / "benchmarks"))
    spec = importlib.util.spec_from_file_location("real_bench_outcome_tests",
                                                  root / "benchmarks" / "real_bench.py")
    real_bench = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(real_bench)

    assert real_bench._outcome_label(pd.Series({
        "m_pearson_signed_auto_auprc_ratio": 1.2,
    })) == "ok, scored"
    assert real_bench._outcome_label(pd.Series({
        "m_pearson_signed_auto_auprc_ratio": np.nan,
        "contexts_ok": 3,
        "reference_edges": 0,
    })) == "completed, unscored (3 contexts; no reference edges)"
    assert real_bench._outcome_label(pd.Series({
        "m_pearson_signed_auto_auprc_ratio": np.nan,
        "contexts_ok": np.nan,
        "reference_edges": 3778,
    })) == "completed, no model score"
    assert real_bench._outcome_label(pd.Series({
        "m_pearson_signed_auto_auprc_ratio": np.nan,
        "n_failed_contexts": 4,
    })) == "completed, screening failures (4 candidate-contexts)"


def test_real_rows_record_screening_failures_without_personal_paths(tmp_path, monkeypatch):
    import sys

    root = Path(__file__).parents[1]
    sys.path.insert(0, str(root / "benchmarks"))
    spec = importlib.util.spec_from_file_location("real_bench_failure_tests",
                                                  root / "benchmarks" / "real_bench.py")
    real_bench = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(real_bench)
    home = tmp_path / "home"
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("USERPROFILE", str(home))   # Path.home() on Windows
    row = {}
    records = {
        "recommended": [
            {"context": "type_a", "error": f"PrepError: failed at {home}/run/a"},
            {"context": "type_b", "error": None},
        ],
        "baseline": [{"context": "all_cells", "error": "PrepError: QC gate failure(s)"}],
    }

    real_bench.record_screening_failures(row, records)

    assert row["n_failed_contexts"] == 2
    assert row["n_gate_refusals"] == 1
    assert str(home) not in row["first_context_error"]
    failures = json.loads(row["context_failures"])
    assert failures[0] == {
        "candidate": "recommended", "context": "type_a",
        "error": "PrepError: failed at ~/run/a",
    }


def test_real_task_manifest_portabilizes_home_paths(tmp_path, monkeypatch):
    import sys

    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("USERPROFILE", str(home))   # Path.home() on Windows
    root = Path(__file__).parents[1]
    sys.path.insert(0, str(root / "benchmarks"))
    spec = importlib.util.spec_from_file_location("real_bench_portable_tests",
                                                  root / "benchmarks" / "real_bench.py")
    real_bench = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(real_bench)
    engine_bench = _engine_bench()
    source = tmp_path / "source.json"
    original = [{"id": "a", "path": str(home / "grnbench_data" / "a.h5ad")}]
    written = engine_bench.write_task_manifest(source, original, task_type="real")
    portable = tmp_path / "portable.json"

    payload = real_bench.portable_task_manifest(source, portable)
    loaded, meta = engine_bench.load_task_manifest(portable)

    assert payload["source_tasks_sha256"] == written["tasks_sha256"]
    assert meta["source_tasks_sha256"] == written["tasks_sha256"]
    assert loaded[0]["path"] == "~/grnbench_data/a.h5ad"
    assert real_bench._expand_task_paths(loaded)[0]["path"] == str(home / "grnbench_data" / "a.h5ad")
