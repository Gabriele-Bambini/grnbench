"""grnbench.viz: figures + index.md built from a real_bench.py-shaped run dir."""
from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from grnbench import cli_viz


def _write_edges(path: Path, rng: np.random.Generator, tfs, targets) -> None:
    rows = [(tf, g, float(rng.normal())) for tf in tfs for g in targets if g != tf]
    df = pd.DataFrame(rows, columns=["source", "target", "weight"])
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, compression="gzip")


def _small_h5ad(rng: np.random.Generator, n_obs=200, n_vars=50) -> ad.AnnData:
    genes = [f"G{i}" for i in range(n_vars - 4)] + ["MT-CO1", "MT-ND1", "TF1", "TF2"]
    X = rng.poisson(2.0, size=(n_obs, n_vars)).astype(np.float32)
    obs = pd.DataFrame({"cell_type": rng.choice(["Tcell", "Bcell", "NK"], size=n_obs)},
                       index=[f"cell{i}" for i in range(n_obs)])
    return ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))


def _build_run_dir(tmp_path: Path) -> Path:
    run_dir = tmp_path / "run"
    rows_dir = run_dir / "rows"
    rows_dir.mkdir(parents=True)
    rng = np.random.default_rng(0)

    tasks = ["task_a", "task_b", "task_c"]
    for i, task in enumerate(tasks):
        row = {
            "task": task,
            "m_pearson_signed_auto_auprc_ratio": 1.2 + 0.1 * i,
            "m_pearson_signed_motif_auto_auprc_ratio": 1.5 + 0.1 * i,
            "m_pearson_signed_baseline_auprc_ratio": 1.0 + 0.05 * i,
            "m_pearson_signed_motif_baseline_auprc_ratio": 1.1 + 0.05 * i,
            "m_pearson_signed_auto_epr": 2.0 + i,
            "m_pearson_signed_motif_auto_epr": 2.3 + i,
            "rules_choices": json.dumps({"normalization": "cpm_log1p", "batch": "none"}),
        }
        (rows_dir / f"{task}.json").write_text(json.dumps(row), encoding="utf-8")

    # artifacts for one task only, exercising edges + reference + processed h5ad
    art = run_dir / "artifacts" / "task_a"
    art.mkdir(parents=True, exist_ok=True)
    tfs, targets = ["TF1", "TF2"], [f"G{i}" for i in range(46)] + ["TF1", "TF2"]
    ref_pairs = [("TF1", "G0"), ("TF1", "G1"), ("TF2", "G2"), ("TF2", "G3")]
    pd.DataFrame(ref_pairs, columns=["source", "target"]).to_csv(
        art / "reference_evaluation.tsv", sep="\t", index=False)
    pd.DataFrame(ref_pairs, columns=["source", "target"]).to_csv(
        art / "reference_selection.tsv", sep="\t", index=False)
    for ctx in ("ctx0", "ctx1"):
        ctx_dir = art / "auto" / ctx
        _write_edges(ctx_dir / "pearson_signed.tsv.gz", rng, tfs, targets)
        _write_edges(ctx_dir / "pearson_signed_motif.tsv.gz", rng, tfs, targets)
        h5 = _small_h5ad(rng)
        h5.write_h5ad(ctx_dir / "processed.h5ad", compression="gzip")
        (ctx_dir / "context.json").write_text(json.dumps(
            {"context": "all_cells", "n_obs": int(h5.n_obs), "n_vars": int(h5.n_vars)}),
            encoding="utf-8")
    decisions = {
        "rules_choices": {"normalization": "cpm_log1p", "batch": "none"},
        "candidates": {"baseline": "baseline", "auto": "candidate_1"},
        "chosen_plans": {
            "auto": {"candidate": "candidate_1", "selection": {}, "warnings": ["mito gate relaxed"],
                     "choices": {"normalization": "cpm_log1p"}},
        },
    }
    (art / "decisions.json").write_text(json.dumps(decisions), encoding="utf-8")
    return run_dir


def test_build_run_writes_index_and_expected_figures(tmp_path):
    run_dir = _build_run_dir(tmp_path)
    out_dir = tmp_path / "viz_out"
    rc = cli_viz.main([str(run_dir), "--out", str(out_dir)])
    assert rc == 0

    index_path = out_dir / "index.md"
    assert index_path.exists()
    index_text = index_path.read_text(encoding="utf-8")

    # scores across datasets (rows/*.json, no artifacts needed)
    assert (out_dir / "scores_auprc_ratio_auto.png").exists()
    assert (out_dir / "lift_pearson_signed_motif_vs_pearson_signed.png").exists()

    # per-dataset figures that need artifacts/task_a
    assert (out_dir / "pr_task_a.png").exists()
    assert (out_dir / "tf_recovery_task_a.png").exists()
    assert (out_dir / "tf_recovery_task_a.tsv").exists()
    assert (out_dir / "calibration_task_a.png").exists()
    assert (out_dir / "network_task_a.png").exists()
    assert (out_dir / "motif_lift_task_a_pearson_signed.png").exists()
    assert (out_dir / "decisions_task_a.md").exists()

    # datasets without artifacts are explicitly listed as not produced, with a reason
    assert "task_b" in index_text
    assert "no edge artifacts" in index_text
    assert "--keep-artifacts" in index_text

    # every figure noted in the index carries an "n = ..." count
    for line in index_text.splitlines():
        if line.startswith("- [") and ".png]" in line:
            assert "n = " in line, line


def test_build_run_accepts_versioned_task_manifest(tmp_path):
    run_dir = _build_run_dir(tmp_path)
    raw_path = tmp_path / "raw.h5ad"
    _small_h5ad(np.random.default_rng(1)).write_h5ad(raw_path, compression="gzip")
    tasks_path = tmp_path / "real_tasks.json"
    tasks_path.write_text(json.dumps({
        "schema_version": 1,
        "task_type": "real",
        "task_count": 1,
        "tasks_sha256": "test-hash",
        "tasks": [{"id": "task_a", "path": str(raw_path)}],
    }), encoding="utf-8")

    out_dir = tmp_path / "viz_versioned"
    from grnbench.viz import build_run
    build_run(run_dir, out_dir=out_dir, tasks_json=tasks_path, max_datasets=1)

    assert (out_dir / "qc_task_a.png").exists()
    assert (out_dir / "embed_task_a.png").exists()


def test_build_run_can_write_portable_global_summary_only(tmp_path):
    run_dir = _build_run_dir(tmp_path)
    out_dir = tmp_path / "viz_summary"
    from grnbench.viz import build_run
    build_run(run_dir, out_dir=out_dir, max_datasets=0)

    assert (out_dir / "scores_auprc_ratio_auto.png").exists()
    assert not list(out_dir.glob("decisions_*"))
    assert "omitted from this summary" in (out_dir / "index.md").read_text()


def test_decisions_markdown_handles_missing_choices_from_dataframe():
    from grnbench.viz import decisions_markdown

    text = decisions_markdown({"rules_choices": float("nan")})
    assert "(none found)" in text


def test_score_distributions_and_paired_lift_title_has_n(tmp_path):
    """The figure functions embed a literal "n = ..." in the suptitle; capture it
    before the figure is closed (build_run's index.md is checked separately)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from grnbench import viz

    run_dir = _build_run_dir(tmp_path)
    rows_df = pd.DataFrame([json.loads(p.read_text(encoding="utf-8"))
                            for p in sorted((run_dir / "rows").glob("*.json"))])

    titles = []
    orig_close = plt.close
    orig_suptitle = matplotlib.figure.Figure.suptitle

    def capture_suptitle(self, t, *a, **kw):
        titles.append(t)
        return orig_suptitle(self, t, *a, **kw)

    matplotlib.figure.Figure.suptitle = capture_suptitle
    try:
        path, reason = viz.score_distributions(rows_df, "auprc_ratio", "auto", tmp_path / "s.png")
        assert reason is None and path is not None
        path, reason = viz.paired_lift(rows_df, "pearson_signed_motif", "pearson_signed",
                                       "auprc_ratio", "auto", tmp_path / "p.png")
        assert reason is None and path is not None
    finally:
        matplotlib.figure.Figure.suptitle = orig_suptitle

    assert titles, "no figure suptitle was set"
    assert all("n = " in t for t in titles)


def test_score_distributions_missing_columns_returns_reason(tmp_path):
    from grnbench import viz

    empty = pd.DataFrame({"task": ["x"]})
    path, reason = viz.score_distributions(empty, "auprc_ratio", "auto", tmp_path / "none.png")
    assert path is None
    assert reason


def test_viz_on_a_path_that_is_not_a_run_creates_nothing_and_says_so(tmp_path, capsys):
    """`grnbench viz /wrong/path` used to create /wrong/path/viz/ and an empty index.md before saying nothing."""
    empty = tmp_path / "not_a_run"
    empty.mkdir()
    assert cli_viz.main([str(empty)]) == 2
    assert not (empty / "viz").exists()
    assert "rows" in capsys.readouterr().err
    missing = tmp_path / "gone" / "deeper"
    assert cli_viz.main([str(missing)]) == 2
    assert not missing.exists() and not missing.parent.exists()
    assert cli_viz.main([str(empty), "--out", str(tmp_path / "chosen_out")]) == 2
    assert not (tmp_path / "chosen_out").exists(), "an explicit --out is not created for a wrong run either"


def test_decisions_markdown_reads_a_candidate_study_row():
    """A candidate_study.py row keeps `candidates` as a list of dicts, not the arm -> candidate mapping."""
    from grnbench import viz

    row = {"candidates": [{"name": "recommended", "source": "engine", "n_contexts": 3},
                          {"name": "baseline", "source": "baseline", "n_contexts": 1}]}
    text = viz.decisions_markdown(row)
    assert "| candidate | recommended | screened (engine, 3 context(s)) |" in text
    assert "| candidate | baseline | screened (baseline, 1 context(s)) |" in text
    assert "(none found)" not in text
    assert "(none found)" in viz.decisions_markdown({"candidates": []})
    assert "candidate (auto) | c1 |" in viz.decisions_markdown({"candidates": {"auto": "c1"}})
