"""Regressions for audited inference and TF-identifier failures."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from sklearn.linear_model import Ridge

from examples.plugins.ridge_model import RidgeGRN
from grnbench.infer.base import run_method
from grnbench.infer.classical import clr, mutual_info
from grnbench import tfs


def _edge_weight(edges: pd.DataFrame, source: str, target: str) -> float:
    row = edges.loc[(edges["source"] == source) & (edges["target"] == target)]
    assert len(row) == 1
    return float(row.iloc[0]["weight"])


def test_mutual_information_of_constant_variables_is_zero():
    expr = pd.DataFrame(np.zeros((16, 3)), columns=["A", "B", "C"])

    edges = mutual_info(expr, tf_names=["A"], n_bins=8)

    assert _edge_weight(edges, "A", "B") == pytest.approx(0.0, abs=1e-12)


def test_mutual_information_is_invariant_to_common_cell_permutation():
    rng = np.random.default_rng(3)
    expr = pd.DataFrame(
        rng.binomial(1, 0.1, (1000, 3)), columns=["A", "B", "C"]
    )

    original = mutual_info(expr, tf_names=["A"], n_bins=8)
    permuted = mutual_info(
        expr.sample(frac=1.0, random_state=5), tf_names=["A"], n_bins=8
    )

    np.testing.assert_allclose(
        original.sort_values(["source", "target"])["weight"],
        permuted.sort_values(["source", "target"])["weight"],
        atol=1e-12,
        rtol=0.0,
    )


def test_clr_zero_variance_background_is_finite_zero():
    expr = pd.DataFrame(np.ones((20, 4)), columns=["A", "B", "C", "D"])

    edges = clr(expr, tf_names=["A"], n_bins=8)

    assert np.isfinite(edges["weight"]).all()
    assert np.count_nonzero(edges["weight"].to_numpy()) == 0


def test_ridge_excludes_each_target_and_matches_independent_fits():
    rng = np.random.default_rng(19)
    a = rng.normal(size=80)
    b = 0.7 * a + rng.normal(scale=0.4, size=80)
    c = 1.8 * a - 0.9 * b + rng.normal(scale=0.1, size=80)
    values = np.column_stack([a, b, c])
    adata = AnnData(values)
    adata.var_names = ["A", "B", "C"]
    alpha = 0.6

    edges = RidgeGRN(alpha=alpha).infer(adata, tf_list=["A", "B"])
    observed = {
        (str(row.source), str(row.target)): float(row.weight) * int(row.direction)
        for row in edges.itertuples(index=False)
    }

    genes = ["A", "B", "C"]
    regulators = ["A", "B"]
    for target_index, target in enumerate(genes):
        predictors = [r for r in regulators if r != target]
        expected = {r: 0.0 for r in regulators}
        if predictors:
            predictor_indices = [genes.index(r) for r in predictors]
            design = values[:, predictor_indices]
            design = (design - design.mean(axis=0)) / np.where(
                design.std(axis=0) > 0, design.std(axis=0), 1.0
            )
            fitted = Ridge(alpha=alpha).fit(design, values[:, target_index])
            expected.update(dict(zip(predictors, fitted.coef_, strict=True)))
        for regulator in regulators:
            assert observed[(regulator, target)] == pytest.approx(
                expected[regulator], rel=1e-10, abs=1e-12
            )


def test_tf_resolution_preserves_mouse_expression_names_and_order():
    genes = ["Nanog", "Pou5f1", "Sox2", "Klf4"]

    resolved = tfs.resolve_tf_names(genes, ["SOX2", "POU5F1"])

    assert resolved == ["Pou5f1", "Sox2"]
    assert tfs.match_tfs_in_data(["sox2", "POU5F1"], genes) == (
        ["SOX2", "POU5F1"],
        2,
    )


def test_tf_resolution_rejects_ambiguous_case_collisions():
    with pytest.raises(tfs.TFListError, match="ambiguous.*case"):
        tfs.resolve_tf_names(["Sox2", "SOX2", "Nanog"], ["sox2"])
    with pytest.raises(tfs.TFListError, match="ambiguous.*case"):
        tfs.match_tfs_in_data(["sox2"], ["Sox2", "SOX2", "Nanog"])


@pytest.mark.parametrize("requested", [[], ["NOT_A_GENE"]])
def test_tf_resolution_rejects_empty_or_zero_match_lists(requested):
    with pytest.raises(tfs.TFListError, match="no.*TF|none.*TF|empty"):
        tfs.resolve_tf_names(["Pou5f1", "Sox2"], requested)
    with pytest.raises(tfs.TFListError, match="no.*TF|none.*TF|empty"):
        tfs.match_tfs_in_data(requested, ["Pou5f1", "Sox2"])


def test_run_method_uses_resolved_names_and_never_expands_zero_matches():
    rng = np.random.default_rng(7)
    adata = AnnData(rng.normal(size=(24, 3)))
    adata.var_names = ["Pou5f1", "Sox2", "Nanog"]

    edges, _ = run_method(
        "spearman_signed", adata, tf_list=["POU5F1"], assume_normalized=True
    )
    assert set(edges["source"]) == {"Pou5f1"}

    with pytest.raises(tfs.TFListError, match="no.*TF|none.*TF"):
        run_method(
            "spearman_signed", adata, tf_list=["TP53"], assume_normalized=True
        )
