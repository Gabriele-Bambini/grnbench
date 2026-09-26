"""Issue #6: per-batch scrublet on tiny batches.

In EXP-20260925-LM-2, 4/20 real datasets had donor/batch levels of 1-25 cells.
Scrublet's 30-component PCA refused them and the whole dataset fell back to the
joint-extreme heuristic. The rule under test: batches below
``SCRUBLET_MIN_CELLS`` are pooled and scored together; a pool still below the
minimum is left unscored (cells kept, score NaN) and recorded. Batch levels that
are missing count as their own level.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import anndata as ad
import scipy.sparse as sp

from grnbench import hygiene as hyg
from grnbench.pipeline import run_prep

pytest.importorskip("skimage")


def _data(tiny_sizes, n_per=250, n_genes=400, n_doublets=30, seed=3):
    """Two programs + synthetic doublets, split into two large batches plus
    tiny batches of the given sizes (cells drawn at random)."""
    rng = np.random.default_rng(seed)
    base = rng.gamma(2.0, 1.0, size=n_genes)
    prog_a = base.copy(); prog_a[: n_genes // 2] *= 6
    prog_b = base.copy(); prog_b[n_genes // 2:] *= 6
    A = rng.poisson(prog_a, size=(n_per, n_genes))
    B = rng.poisson(prog_b, size=(n_per, n_genes))
    D = A[:n_doublets] + B[:n_doublets]
    X = np.vstack([A, B, D]).astype(np.float32)
    n = len(X)
    labels = np.empty(n, dtype=object)
    order = rng.permutation(n)
    start = 0
    for i, size in enumerate(tiny_sizes):
        labels[order[start:start + size]] = f"tiny{i}"
        start += size
    rest = order[start:]
    labels[rest[: len(rest) // 2]] = "big0"
    labels[rest[len(rest) // 2:]] = "big1"
    obs = pd.DataFrame({"batch": labels}, index=[f"c{i}" for i in range(n)])
    return ad.AnnData(X=sp.csr_matrix(X), obs=obs,
                      var=pd.DataFrame(index=[f"G{i}" for i in range(n_genes)]))


def test_tiny_batches_are_pooled_and_scored_not_a_fallback():
    adata = _data([13, 20, 25, 25, 25])  # 108 cells in tiny batches
    method, is_dbl, warning, groups = hyg.detect_doublets(adata, batch_key="batch", seed=0)
    assert method == "scrublet"
    assert np.isfinite(adata.obs["doublet_score"].to_numpy(float)).all()
    assert groups["min_cells"] == hyg.SCRUBLET_MIN_CELLS
    assert groups["per_batch"] == 2
    assert groups["pooled_batches"] == 5 and groups["pooled_cells"] == 108
    assert groups["unscored_cells"] == 0
    assert "pooled" in warning


def test_pool_below_minimum_is_left_unscored_and_kept():
    adata = _data([13, 20, 27])  # 60 cells: even pooled, below the minimum
    tiny = adata.obs["batch"].str.startswith("tiny").to_numpy()
    method, is_dbl, warning, groups = hyg.detect_doublets(adata, batch_key="batch", seed=0)
    assert method == "scrublet"
    assert groups["unscored_cells"] == 60 and groups["pooled_cells"] == 60
    assert np.isnan(adata.obs["doublet_score"].to_numpy(float)[tiny]).all()
    assert np.isfinite(adata.obs["doublet_score"].to_numpy(float)[~tiny]).all()
    assert not is_dbl[tiny].any()
    assert adata.obs["predicted_doublet"].dtype == bool
    assert "unscored" in warning and "60" in warning


def test_missing_batch_labels_are_their_own_level():
    adata = _data([])
    adata.obs.loc[adata.obs.index[:40], "batch"] = np.nan
    method, _is_dbl, _warning, groups = hyg.detect_doublets(adata, batch_key="batch", seed=0)
    assert method == "scrublet"
    assert groups["pooled_batches"] == 1 and groups["pooled_cells"] == 40


def test_all_batches_tiny_pool_everything():
    adata = _data([60] * 8 + [50])  # 530 cells, every batch below the minimum
    method, _is_dbl, _warning, groups = hyg.detect_doublets(adata, batch_key="batch", seed=0)
    assert method == "scrublet"
    assert groups["per_batch"] == 0 and groups["pooled_cells"] == adata.n_obs
    assert np.isfinite(adata.obs["doublet_score"].to_numpy(float)).all()


def test_dataset_below_minimum_falls_back_loudly():
    adata = _data([])[:80].copy()
    del adata.obs["batch"]
    with pytest.warns(RuntimeWarning, match="DOUBLET FALLBACK"):
        method, _is_dbl, warning, groups = hyg.detect_doublets(adata, seed=0)
    assert method == "heuristic"
    assert str(hyg.SCRUBLET_MIN_CELLS) in warning
    assert groups["unscored_cells"] == 80
    assert groups["batch_key"] is None and groups["pooled_batches"] == 0


def test_large_batches_only_match_plain_per_batch_scrublet():
    """No small batch: the rule changes nothing, scores equal scanpy's own per-batch run."""
    import scanpy as sc
    adata = _data([])
    ref = adata.copy()
    sc.pp.scrublet(ref, batch_key="batch", random_state=0, verbose=False)
    method, is_dbl, warning, groups = hyg.detect_doublets(adata, batch_key="batch", seed=0)
    assert method == "scrublet" and warning is None and groups["pooled_batches"] == 0
    np.testing.assert_allclose(adata.obs["doublet_score"].to_numpy(float),
                               ref.obs["doublet_score"].to_numpy(float))
    assert (is_dbl == ref.obs["predicted_doublet"].to_numpy(bool)).all()
    assert "_scrublet_group" not in adata.obs and "batches" in adata.uns["scrublet"]


def test_manifest_records_the_grouping(tmp_path):
    # 1230 cells: on this two-program toy scrublet over-calls small groups (the
    # 108-cell pool calls ~20%), so the large batches must carry the dataset
    # under the 20% ceiling. On real data pooled calls matched per-batch ones
    # (0-3% in the four LM-2 datasets that used to fall back; PR for issue #6).
    src = tmp_path / "in.h5ad"
    _data([13, 20, 25, 25, 25], n_per=600).write_h5ad(src)
    res = run_prep(src, tmp_path / "run", drop_doublets=True, batch_key="batch", seed=0)
    step = next(s for s in res.manifest["steps"] if s["step"] == "drop_doublets")
    assert step["method"] == "scrublet" and step["fallback"] is False
    assert step["scrublet_groups"]["pooled_batches"] == 5
    assert "pooled" in step["detail"]
