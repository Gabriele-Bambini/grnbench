"""Regressions for audited inference and TF-identifier failures."""
from __future__ import annotations

import sys
import types

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from sklearn.linear_model import Ridge

from examples.plugins.ridge_model import RidgeGRN
from grnbench.infer import tree_methods
from grnbench.infer.base import run_method
from grnbench.infer.classical import _fit_lasso, clr, mutual_info
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


def test_lasso_default_path_scores_an_orthogonal_regulator_zero():
    """An exactly uncorrelated pair is no evidence of regulation, at any offset.

    ``lasso_path`` fits no intercept, so uncentered columns let the means buy a
    coefficient: this orthogonal pair scored 0.487, and 0.971 once 10 was added
    to both genes.
    """
    x = np.array([[0.0], [0.0], [1.0], [1.0]])
    y = np.array([0.0, 1.0, 0.0, 1.0])

    assert _fit_lasso(x, y) == pytest.approx([0.0], abs=1e-12)
    assert _fit_lasso(x + 10, y + 10) == pytest.approx([0.0], abs=1e-12)


def test_lasso_provenance_marks_the_centered_default_path():
    """Old (uncentered) and new default-path lasso runs must be told apart (D-005)."""
    expr = pd.DataFrame(np.random.default_rng(0).normal(size=(30, 3)), columns=["A", "B", "C"])
    adata = AnnData(expr.to_numpy(np.float32), var=pd.DataFrame(index=expr.columns))

    _, default = run_method("lasso", adata, tf_list=["A", "B"], assume_normalized=True)
    _, explicit = run_method("lasso", adata, tf_list=["A", "B"], assume_normalized=True, alpha=0.1)

    assert default["params"]["path"] == "centered"
    assert "path" not in explicit["params"]  # the explicit-alpha branch did not change


def test_random_forest_uses_the_requested_thread_budget(monkeypatch):
    """n_jobs=-1 takes every core on the node, whatever --threads said."""
    import sklearn.ensemble

    seen = []

    class _RecordingForest:
        def __init__(self, **kwargs):
            seen.append(kwargs)

        def fit(self, X, y):
            self.feature_importances_ = np.zeros(X.shape[1])
            return self

    monkeypatch.setattr(sklearn.ensemble, "RandomForestRegressor", _RecordingForest)
    expr = pd.DataFrame(np.random.default_rng(0).normal(size=(20, 3)), columns=["A", "B", "C"])

    run_method("random_forest", AnnData(expr.to_numpy(np.float32), var=pd.DataFrame(index=expr.columns)),
               tf_list=["A"], assume_normalized=True, threads=3)

    assert seen and all(kwargs["n_jobs"] == 3 for kwargs in seen)


def test_arboreto_client_gets_the_requested_thread_budget(monkeypatch):
    """A dask client with no budget sizes itself from the machine, not from --threads."""
    seen = {}

    class _RecordingClient:
        def __init__(self, **kwargs):
            seen.update(kwargs)

    fake = types.ModuleType("dask.distributed")
    fake.Client = _RecordingClient
    monkeypatch.setitem(sys.modules, "dask", types.ModuleType("dask"))
    monkeypatch.setitem(sys.modules, "dask.distributed", fake)

    assert tree_methods._quiet_client(threads=3) is not None
    assert seen["n_workers"] == 1 and seen["threads_per_worker"] == 3


@pytest.mark.parametrize("runner", ["genie3_runner", "grnboost2_runner"])
def test_arboreto_runners_forward_the_thread_budget(monkeypatch, runner):
    captured = {}

    def _fake(expr, tf_names, seed, **kwargs):
        captured.update(kwargs)
        return pd.DataFrame(columns=["TF", "target", "importance"])

    monkeypatch.setattr(tree_methods, "_run_arboreto", _fake)
    expr = pd.DataFrame(np.zeros((4, 2)), columns=["A", "B"])

    getattr(tree_methods, runner)(expr, tf_names=["A"], seed=1, threads=3)

    assert captured["threads"] == 3


def _fake_arboreto(monkeypatch):
    """Stand-in for arboreto's core: records create_graph's kwargs, returns no edges."""
    seen = {}

    def create_graph(**kwargs):
        seen.update(kwargs)
        return "graph", None

    core = types.SimpleNamespace(create_graph=create_graph, DEMON_SEED=666)
    monkeypatch.setattr(tree_methods, "_load", lambda: (core, {}, {}))
    return seen


def test_run_arboreto_passes_the_thread_budget_to_its_client(monkeypatch):
    seen = _fake_arboreto(monkeypatch)
    budgets = []

    class _Client:
        def compute(self, graph, sync):
            return pd.DataFrame({"TF": [], "target": [], "importance": []})

        def close(self):
            pass

    monkeypatch.setattr(tree_methods, "_quiet_client", lambda threads=None: budgets.append(threads) or _Client())
    expr = pd.DataFrame(np.zeros((4, 2)), columns=["A", "B"])

    tree_methods._run_arboreto(expr, ["A"], 1, regressor_type="RF",
                               early_stop_window_length=0, threads=3)

    assert budgets == [3] and isinstance(seen["client"], _Client)


@pytest.mark.parametrize("threads", [3, None])
def test_arboreto_client_construction_failure_is_raised_with_its_cause(monkeypatch, threads):
    """A Client that cannot start used to become None, then `None.compute` (AttributeError)."""
    _fake_arboreto(monkeypatch)

    def _refuse(**kwargs):
        raise OSError("no free port")

    fake = types.ModuleType("dask.distributed")
    fake.Client = _refuse
    monkeypatch.setitem(sys.modules, "dask", types.ModuleType("dask"))
    monkeypatch.setitem(sys.modules, "dask.distributed", fake)
    expr = pd.DataFrame(np.zeros((4, 2)), columns=["A", "B"])

    with pytest.raises(RuntimeError, match="dask Client") as err:
        tree_methods._run_arboreto(expr, ["A"], 1, regressor_type="RF",
                                   early_stop_window_length=0, threads=threads)
    assert isinstance(err.value.__cause__, OSError)


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


def test_default_lasso_path_raises_no_future_warning():
    """scikit-learn 1.9 deprecates lasso_path(n_alphas=) and 1.11 removes it."""
    import warnings

    rng = np.random.default_rng(0)
    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        _fit_lasso(rng.normal(size=(30, 4)), rng.normal(size=30))


def test_extra_trees_uses_the_requested_thread_budget(monkeypatch):
    import sklearn.ensemble

    seen = []

    class _RecordingForest:
        def __init__(self, **kwargs):
            seen.append(kwargs)

        def fit(self, X, y):
            self.feature_importances_ = np.zeros(X.shape[1])
            return self

    monkeypatch.setattr(sklearn.ensemble, "ExtraTreesRegressor", _RecordingForest)
    expr = pd.DataFrame(np.random.default_rng(0).normal(size=(20, 3)), columns=["A", "B", "C"])

    run_method("extra_trees", AnnData(expr.to_numpy(np.float32), var=pd.DataFrame(index=expr.columns)),
               tf_list=["A"], assume_normalized=True, threads=3)

    assert seen and all(kwargs["n_jobs"] == 3 for kwargs in seen)
