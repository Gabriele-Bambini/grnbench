"""Inference adapter tests: contract, baseline recovery, tree methods, aracne3 parser."""
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
from grnbench.infer.base import get_method

CASCADE_PAIRS = {("TF_A", "G_B"), ("G_B", "G_C"), ("TF_A", "G_C")}


def test_list_methods_reports_unavailability():
    table = list_methods()
    names = {m["name"] for m in table}
    assert {"spearman_signed", "genie3", "grnboost2", "aracne3"} <= names
    arboreto_entry = next(m for m in table if m["name"] == "grnboost2")
    assert arboreto_entry["available"] or "arboreto" in arboreto_entry["unavailable_reason"]
    aracne3_entry = next(m for m in table if m["name"] == "aracne3")
    # Two prerequisites, probed in order: java on PATH, then ARACNE3_JAR.
    assert aracne3_entry["available"] or any(
        piece in aracne3_entry["unavailable_reason"] for piece in ("java", "ARACNE3_JAR"))


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


def test_aracne3_parser_handles_format_drift():
    from grnbench.infer.aracne3 import parse_aracne3_text

    text = (
        "# comment line that must be ignored\n"
        "Regulator\tTarget\tScore\tSign\n"
        "TP53\tCDKN1A\t12.4\t+\n"
        "TP53\tBAX\t8.1\tactivation\n"
        "MYC\tMDM2\t-0.5\t-\n"   # negative weights allowed; sign maps to -1
        "broken line here\n"
    )
    net = parse_aracne3_text(text)
    assert list(net.columns) == ["source", "target", "weight", "direction"]
    assert len(net) == 3
    cdk = net[net.target == "CDKN1A"].iloc[0]
    assert cdk["source"] == "TP53" and cdk["weight"] == 12.4 and cdk["direction"] == 1
    mdm = net[net.target == "MDM2"].iloc[0]
    assert mdm["direction"] == -1


def test_aracne3_runner_raises_without_jar(monkeypatch, cascade_adata):
    from grnbench.infer import base as infer_base
    from grnbench.infer.aracne3 import aracne3_availability

    ok, why = aracne3_availability()
    if ok:
        pytest.skip("java + ARACNE3_JAR available on this host")
    monkeypatch.delenv("ARACNE3_JAR", raising=False)
    with pytest.raises(MethodError, match="unavailable"):
        expr_df, _ = infer_base.prepare_expression(cascade_adata)
        spec = get_method("aracne3")
        spec.runner(expr_df, tf_names=None)


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
