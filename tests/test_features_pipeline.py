"""Feature selection (TF-union rule) and end-to-end pipeline tests."""
from __future__ import annotations

import collections
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import anndata as ad
from scipy import sparse

from grnbench.features import select_features
from grnbench.pipeline import PrepError, run_prep
from grnbench.qc import Thresholds


def _weaken_gene(adata: ad.AnnData, gene: str) -> ad.AnnData:
    """Make ``gene`` barely expressed so HVG selection will ignore it."""
    out = adata.copy()
    X = out.X.toarray() if sparse.issparse(out.X) else np.asarray(out.X)
    j = list(out.var_names).index(gene)
    X[:, j] = (X[:, j] > 0) * 1.0
    out.X = sparse.csr_matrix(X.astype(np.float32))
    return out


def test_tf_union_rule_rescues_weak_tf(synthetic_adata):
    weak = _weaken_gene(synthetic_adata, "GENE100")
    out, rep = select_features(weak, n_hvg=50, tf_list=["GENE100"])
    assert rep["tfs_matched"] == 1
    assert rep["tfs_added_by_union"] >= 1
    assert "GENE100" in out.var_names
    assert rep["n_features_kept"] <= 51


def test_select_features_without_tf_list(synthetic_adata):
    out, rep = select_features(synthetic_adata, n_hvg=100)
    assert out.n_vars <= 101
    assert rep["tfs_requested"] == 0
    assert rep["n_hvg_selected"] >= 1


def test_end_to_end_run(tmp_path, synthetic_adata):
    src = tmp_path / "in.h5ad"
    synthetic_adata.write_h5ad(src)
    res = run_prep(src, tmp_path / "run", n_hvg=150, tf_list=["GENE100"], seed=1)

    assert res.processed_path.exists()
    assert res.manifest_path.exists()
    assert res.report_path.exists()
    proc = ad.read_h5ad(res.processed_path)
    assert proc.n_vars <= 151  # 150 HVGs + rescued TF at most
    assert res.manifest["strategy_chosen"]
    assert len(res.manifest["steps"]) >= 5
    text = res.report_path.read_text(encoding="utf-8")
    assert "## QC gates" in text and "## Steps" in text


def test_auto_strategy_default_is_cpm_log1p(tmp_path, synthetic_adata):
    src = tmp_path / "in.h5ad"
    synthetic_adata.write_h5ad(src)
    res = run_prep(src, tmp_path / "run_auto")  # auto, no group key
    assert res.manifest["strategy_chosen"] == "cpm_log1p"


def test_auto_recommends_pseudobulk_when_groups_tiny(tmp_path, synthetic_adata):
    a = synthetic_adata.copy()
    a.obs["group"] = [f"g{i % 40}" for i in range(a.n_obs)]  # 40 groups x 10 cells
    src = tmp_path / "in.h5ad"
    a.write_h5ad(src)
    res = run_prep(src, tmp_path / "run_pb", group_key="group")
    top = res.manifest["recommendations"][0]["strategy"]
    assert top == "pseudobulk"
    assert res.strategy == "pseudobulk"
    assert res.processed_path.exists()


def test_gate_failure_aborts_then_force_proceeds(tmp_path, extreme_mito_adata):
    src = tmp_path / "bad.h5ad"
    extreme_mito_adata.write_h5ad(src)
    with pytest.raises(PrepError, match="QC gate failure"):
        run_prep(src, tmp_path / "r1")
    # extreme_mito_adata has only 10 genes; the default min_genes_per_cell=200
    # would filter out every cell regardless of --force (a separate, later
    # failure the pipeline now also refuses loudly - see the dedicated
    # empty-matrix test below), so relax it here to isolate what --force
    # actually does: proceed past the mito gate to a real, non-empty result.
    res = run_prep(src, tmp_path / "r2", force=True,
                   thresholds=Thresholds(min_genes_per_cell=1))
    assert res.processed_path.exists()
    assert res.manifest.get("forced_past_gates") is True
    assert ad.read_h5ad(res.processed_path).n_obs > 0


def test_force_cannot_rescue_a_result_every_filter_empties(tmp_path, extreme_mito_adata):
    """--force bypasses the QC *gate* abort, but a downstream filter
    (min_genes_per_cell) can still remove every cell on its own; that must
    fail loudly rather than write a 0-cell processed.h5ad as success."""
    src = tmp_path / "bad.h5ad"
    extreme_mito_adata.write_h5ad(src)
    with pytest.raises(PrepError, match="empty matrix"):
        run_prep(src, tmp_path / "r3", force=True)  # default min_genes_per_cell=200 > 10 genes total


def test_manifest_has_environment_and_fingerprint(tmp_path, synthetic_adata):
    src = tmp_path / "in.h5ad"
    synthetic_adata.write_h5ad(src)
    res = run_prep(src, tmp_path / "run_env")
    env = res.manifest["environment"]
    assert "scanpy" in env and "python" in env
    fp = res.manifest["input"]
    assert fp["size_bytes"] > 0


# --------------------------------------------------------------------------- #
# stratified subsampling and filter arithmetic
# --------------------------------------------------------------------------- #

def _labelled_adata(labels, n_genes=400, seed=0):
    # >= Thresholds.min_genes_per_cell (200), else filter_cells empties
    # the matrix and the subsample branch never runs
    import anndata as ad
    rng = np.random.default_rng(seed)
    n = len(labels)
    return ad.AnnData(
        X=rng.poisson(3.0, size=(n, n_genes)).astype("float32"),
        obs=pd.DataFrame({"ct": labels}, index=[f"c{i}" for i in range(n)]),
        var=pd.DataFrame(index=[f"G{i}" for i in range(n_genes)]))


def test_stratified_subsample_keeps_the_majority_class():
    """A class holding half the cells must not be truncated out of existence.

    Bumping every zero-quota class to one cell overshot n_keep, and the
    `picked[:n_keep]` slice that repaired it cut in alphabetical class order -
    so 200 singleton classes crowded out a class holding 50% of the data, which
    ended up with ZERO cells while the manifest still said "stratified by".
    """
    from grnbench.pipeline import _stratified_indices

    labels = np.array([f"rare{i:03d}" for i in range(200)] + ["zzz_major"] * 200)
    idx, info = _stratified_indices(labels, 20, np.random.default_rng(42))
    kept = collections.Counter(labels[idx])
    assert len(idx) == 20, "the allocation must spend exactly its budget"
    assert kept["zzz_major"] == 10, "50% of the cells earns 50% of the subsample"
    assert info["classes_dropped"] == 190
    assert info["min_per_class_honoured"] is False


def test_stratified_subsample_is_proportional_when_it_fits():
    from grnbench.pipeline import _stratified_indices

    labels = np.array(["a"] * 60 + ["b"] * 30 + ["c"] * 10)
    idx, info = _stratified_indices(labels, 20, np.random.default_rng(0))
    kept = collections.Counter(labels[idx])
    assert (kept["a"], kept["b"], kept["c"]) == (12, 6, 2)
    assert info["classes_dropped"] == 0 and info["min_per_class_honoured"]


def test_dropped_strata_are_surfaced_as_a_warning(tmp_path):
    from grnbench.pipeline import run_prep

    labels = [f"rare{i:03d}" for i in range(60)] + ["major"] * 140
    p = tmp_path / "in.h5ad"
    _labelled_adata(labels).write_h5ad(p)
    res = run_prep(p, tmp_path / "out", strategy="cpm_log1p", force=True,
                   subsample_cells=10, stratify_key="ct", seed=1)
    assert any("class(es) are absent from the subsample" in w
               for w in res.manifest["warnings"])
    step = [s for s in res.manifest["steps"] if s["step"] == "subsample_cells"][0]
    assert step["mode"] == "stratified"
    assert step["stratification"]["classes_dropped"] > 0


def test_all_missing_stratification_labels_fail_one_context_not_the_dataset(tmp_path):
    """A bad context is recorded as a PrepError; later contexts still run (#169)."""
    from grnbench.autoprep import PrepPlan, execute_plan

    src = tmp_path / "in.h5ad"
    _labelled_adata(["unknown"] * 120 + ["beta"] * 120).write_h5ad(src)
    plan = PrepPlan(contexts=[
        {"name": "ct=unknown", "subset_key": "ct", "subset_value": "unknown"},
        {"name": "ct=beta", "subset_key": "ct", "subset_value": "beta"},
    ])
    plan.add("strategy", "Normalization", "cpm_log1p", "test")
    plan.add("stratify_key", "Stratification column", "ct", "test")
    plan.add("subsample_cells", "Subsample cells", 50, "test")
    plan.add("subsample_mode", "Subsample mode", "stratified", "test")

    records = execute_plan(plan, src, tmp_path / "run")

    assert records[0]["error"] == (
        "PrepError: cannot stratify a subsample by 'ct': all 120 cells in this context "
        "have missing labels; use random subsampling or a different stratification column"
    )
    assert records[0]["processed"] is None
    assert records[1]["error"] is None and Path(records[1]["processed"]).exists()
    saved = json.loads((tmp_path / "run" / "contexts.json").read_text(encoding="utf-8"))
    assert saved == records


def test_explicit_random_subsample_does_not_silently_stratify(tmp_path):
    p = tmp_path / "in.h5ad"
    _labelled_adata(["major"] * 300 + ["minor"] * 100).write_h5ad(p)
    res = run_prep(
        p, tmp_path / "out", strategy="cpm_log1p", force=True,
        subsample_cells=80, subsample_mode="random", stratify_key="ct", seed=1)
    step = [s for s in res.manifest["steps"] if s["step"] == "subsample_cells"][0]
    assert step["mode"] == "random"
    assert step["stratification"] is None
    assert step["detail"].startswith("random (requested)")
    assert res.manifest["args"]["subsample_mode"] == "random"


def test_legacy_missing_stratify_key_records_the_random_fallback(tmp_path):
    p = tmp_path / "in.h5ad"
    _labelled_adata(["major"] * 300 + ["minor"] * 100).write_h5ad(p)
    res = run_prep(
        p, tmp_path / "out", strategy="cpm_log1p", force=True,
        subsample_cells=80, stratify_key="not_a_column", seed=1)
    step = [s for s in res.manifest["steps"] if s["step"] == "subsample_cells"][0]
    assert step["mode"] == "random"
    assert step["stratification"] is None
    assert step["detail"].startswith("random (no stratify key)")


def test_explicit_stratified_subsample_requires_a_key_before_writing_output(tmp_path, synthetic_adata):
    p = tmp_path / "in.h5ad"
    synthetic_adata.write_h5ad(p)
    out = tmp_path / "out"
    with pytest.raises(PrepError, match="needs --stratify-key"):
        run_prep(p, out, subsample_cells=50, subsample_mode="stratified")
    assert not out.exists()


def test_explicit_stratified_subsample_refuses_a_missing_obs_key(tmp_path, synthetic_adata):
    p = tmp_path / "in.h5ad"
    synthetic_adata.write_h5ad(p)
    out = tmp_path / "out"
    with pytest.raises(PrepError, match="is not in obs"):
        run_prep(
            p, out, subsample_cells=50, subsample_mode="stratified",
            stratify_key="not_a_column")
    assert not (out / "manifest.json").exists()


def test_filter_cells_counts_add_up(tmp_path):
    """Zero-count cells are a SUBSET of low-gene cells, never a separate group."""
    from grnbench.pipeline import run_prep

    rng = np.random.default_rng(1)
    n, g = 300, 400
    X = rng.poisson(1.0, size=(n, g)).astype("float32")
    X[:10, :] = 0                       # zero-count
    X[10:25, :] = 0; X[10:25, :5] = 3   # low-gene but non-zero
    p = tmp_path / "f.h5ad"
    _labelled_adata(["a"] * n, n_genes=g).copy()
    import anndata as ad
    ad.AnnData(X=X, obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]),
               var=pd.DataFrame(index=[f"G{i}" for i in range(g)])).write_h5ad(p)

    steps = run_prep(p, tmp_path / "o", strategy="cpm_log1p",
                     force=True).manifest["steps"]
    fc = [s for s in steps if s["step"] == "filter_cells"][0]
    assert fc["n_cells_removed"] == 25
    assert fc["n_cells_removed"] + fc["n_cells_after"] == n
    assert "of which 10 had zero counts" in fc["detail"]
