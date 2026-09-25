"""Strategy registry and normalization behavior tests."""
from __future__ import annotations

import numpy as np
import pytest

from grnbench.ingest import looks_like_counts
from grnbench.strategies import (
    StrategyError,
    apply_strategy,
    get_strategy,
    list_strategies,
)


def test_cpm_log1p_normalizes_rows_to_target(synthetic_adata):
    out, prov = apply_strategy("cpm_log1p", synthetic_adata)
    recovered = np.asarray(np.expm1(out.X).sum(axis=1)).ravel()
    assert np.allclose(recovered, 1e4, rtol=1e-3)
    assert out.X.min() >= 0
    assert "counts" in out.layers
    assert prov["strategy"] == "cpm_log1p"
    # input never mutated
    ok, _ = looks_like_counts(synthetic_adata.X)
    assert ok


def test_pseudobulk_sums_are_exact_and_ordered(synthetic_adata):
    out, prov = apply_strategy("pseudobulk", synthetic_adata, group_key="group")
    assert out.n_obs == 3
    counts_layer = out.layers["counts"]
    orig = synthetic_adata.layers["counts"]
    groups = list(out.obs["group"])
    for gi, grp in enumerate(groups):
        mask = (synthetic_adata.obs["group"] == grp).to_numpy()
        expected = np.asarray(orig[mask].sum(axis=0)).ravel()
        row = counts_layer[gi]
        got = np.asarray(row.todense()).ravel() if hasattr(row, "todense") \
            else np.asarray(row).ravel()
        assert np.allclose(expected, got), f"group {grp}"
    assert list(out.obs["n_cells"]) == [200, 195, 5]


def test_pseudobulk_missing_key_raises(synthetic_adata):
    with pytest.raises(StrategyError, match="group"):
        apply_strategy("pseudobulk", synthetic_adata, group_key="nope")


def test_pseudobulk_warns_on_tiny_groups(synthetic_adata):
    with pytest.warns(UserWarning, match="< 10 cells"):
        apply_strategy("pseudobulk", synthetic_adata, group_key="group")


def test_unavailable_strategy_raises_with_reason():
    with pytest.raises(StrategyError, match="scran"):
        get_strategy("scran_pooled")


def test_registry_contains_expected_entries():
    entries = {s["name"] for s in list_strategies()}
    assert {"cpm_log1p", "pearson_residuals", "pseudobulk",
            "raw_counts", "scran_pooled", "scvi_latent"} <= entries
    avail = {s["name"] for s in list_strategies() if s["available"]}
    assert {"cpm_log1p", "pseudobulk", "raw_counts"} <= avail


def test_no_imputation_strategy_exists():
    """Hard rule: the registry must never grow an imputation entry."""
    names = {s["name"] for s in list_strategies()}
    forbidden = {"magic", "impute", "knn_smooth", "alchemy"}
    assert not (names & forbidden)


def test_pearson_residuals_smoke(synthetic_adata):
    try:
        out, prov = apply_strategy("pearson_residuals", synthetic_adata)
    except StrategyError:
        pytest.skip("scanpy experimental API unavailable")
    assert out.shape == synthetic_adata.shape
    dense = out.X.toarray() if hasattr(out.X, "toarray") else np.asarray(out.X)
    assert np.isfinite(dense).all()
