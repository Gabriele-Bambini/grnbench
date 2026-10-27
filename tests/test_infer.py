"""Inference adapter tests: contract, baseline recovery, tree methods (ARACNe3: test_aracne3.py)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from grnbench.infer import (
    EDGE_COLUMNS,
    MethodError,
    list_methods,
    load_edges,
    run_method,
    save_edges,
)

CASCADE_PAIRS = {("TF_A", "G_B"), ("G_B", "G_C"), ("TF_A", "G_C")}


def test_list_methods_reports_unavailability():
    table = list_methods()
    names = {m["name"] for m in table}
    assert {"spearman_signed", "genie3", "grnboost2", "aracne3"} <= names
    arboreto_entry = next(m for m in table if m["name"] == "grnboost2")
    assert arboreto_entry["available"] or "arboreto" in arboreto_entry["unavailable_reason"]
    aracne3_entry = next(m for m in table if m["name"] == "aracne3")
    # one prerequisite: the native ARACNe3 program (GRNBENCH_ARACNE3_BIN, or ARACNe3_app_release on PATH)
    assert aracne3_entry["available"] or "GRNBENCH_ARACNE3_BIN" in aracne3_entry["unavailable_reason"]


def test_spearman_recovers_cascade_and_sign(cascade_adata):
    edges, prov = run_method("spearman_signed", cascade_adata)
    assert list(edges.columns[:3]) == EDGE_COLUMNS
    assert "direction" in edges.columns
    assert edges["weight"].is_monotonic_decreasing
    top_pairs = set(map(tuple, edges.head(6)[["source", "target"]].to_numpy()))
    assert top_pairs & CASCADE_PAIRS, f"cascade not recovered: {top_pairs}"
    ab = edges[(edges.source == "TF_A") & (edges.target == "G_B")].iloc[0]
    assert ab["weight"] > 0.5 and ab["direction"] == 1
    assert all(e["source"] != e["target"] for _, e in edges.iterrows())
    assert prov["normalization"].startswith("raw counts detected")


def test_tf_list_restricts_regulators(cascade_adata):
    edges, prov = run_method("spearman_signed", cascade_adata, tf_list=["TF_A"])
    assert (edges["source"] == "TF_A").all()
    assert prov["regulators"] == "tf_list"


def test_unknown_and_unavailable_methods_raise(cascade_adata):
    with pytest.raises(MethodError, match="Unknown"):
        run_method("nope", cascade_adata)
    table = {m["name"]: m for m in list_methods()}
    unavailable = [n for n, m in table.items() if not m["available"]]
    if unavailable:
        with pytest.raises(MethodError, match="unavailable"):
            run_method(unavailable[0], cascade_adata)


def test_edges_roundtrip(tmp_path, cascade_adata):
    edges, _ = run_method("spearman_signed", cascade_adata)
    path = save_edges(edges, tmp_path, "demo")
    back = load_edges(path)
    assert set(EDGE_COLUMNS) <= set(back.columns)  # direction may travel along
    assert len(back) == len(edges)
    assert back["weight"].is_monotonic_decreasing


@pytest.mark.skipif(
    __import__("importlib.util", fromlist=["find_spec"]).find_spec("arboreto") is None,
    reason="arboreto not installed",
)
class TestTreeMethods:
    def test_genie3_schema_and_signal(self, cascade_adata):
        edges, prov = run_method("genie3", cascade_adata, n_estimators=150)
        assert list(edges.columns) == EDGE_COLUMNS
        assert len(edges) > 0
        top = set(map(tuple, edges.head(15)[["source", "target"]].to_numpy()))
        assert top & CASCADE_PAIRS, f"weak signal? top15={top}"
        assert prov["method"] == "genie3"

    def test_grnboost2_schema(self, cascade_adata):
        edges, prov = run_method("grnboost2", cascade_adata, seed=7)
        assert list(edges.columns) == EDGE_COLUMNS
        assert len(edges) > 0
        assert edges["weight"].is_monotonic_decreasing


@pytest.mark.parametrize("sparse", [False, True])
@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
def test_a_matrix_with_nan_or_inf_is_refused_not_turned_into_zeros(cascade_adata, sparse, bad):
    """One NaN used to give pearson_signed/spearman_signed an all-zero network and exit code 0."""
    import scipy.sparse as sp
    adata = cascade_adata.copy()
    X = np.asarray(adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X, dtype=float)
    X[3, 1] = bad
    X[7, 1] = bad
    X[5, 2] = bad
    adata.X = sp.csr_matrix(X) if sparse else X
    for method in ("spearman_signed", "pearson_signed"):
        with pytest.raises(MethodError, match=r"3 non-finite value\(s\).*in 2 gene\(s\) and 3 cell\(s\)"):
            run_method(method, adata)
    with pytest.raises(MethodError, match="non-finite"):       # also when the caller vouches for the scale
        run_method("pearson_signed", adata, assume_normalized=True)
    clean = cascade_adata.copy()
    edges, _ = run_method("pearson_signed", clean)               # the guard is silent on good data
    assert len(edges) and np.isfinite(edges["weight"]).all()


def test_a_plugin_that_edits_its_input_cannot_corrupt_the_next_method(cascade_adata):
    """With --assume-normalized the plugin received the caller's own AnnData, so an in-place edit leaked on."""
    from grnbench.infer.base import REGISTRY, MethodSpec, register_method

    def vandal(adata, tf_list, **kw):
        (adata.X.data if hasattr(adata.X, "data") and hasattr(adata.X, "tocsr") else adata.X)[:] = 0
        return pd.DataFrame({"source": ["TF_A"], "target": ["G_B"], "weight": [1.0]})

    register_method(MethodSpec("vandal_test", "edits X in place", vandal, takes="adata", origin="python"),
                    replace=True)
    try:
        before = np.asarray(cascade_adata.X.toarray() if hasattr(cascade_adata.X, "toarray") else cascade_adata.X).copy()
        run_method("vandal_test", cascade_adata, assume_normalized=True)
        after = np.asarray(cascade_adata.X.toarray() if hasattr(cascade_adata.X, "toarray") else cascade_adata.X)
        assert before.any() and np.array_equal(before, after), "the caller's matrix was edited by the plugin"
    finally:
        REGISTRY.pop("vandal_test", None)
