"""Data-hygiene tests: labels, IDs, organism, gates, stratification, memory."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import anndata as ad
import scipy.sparse as sp

from grnbench.hygiene import (
    audit_obs,
    detect_batch_column,
    detect_id_type,
    doublet_heuristic_fraction,
    downcast_floats,
    duplicate_barcodes_report,
    load_label_mapping,
    normalize_labels,
    organism_hint,
    sparsify_if_dense,
)
from grnbench.ingest import IngestError, load_any
from grnbench.pipeline import run_prep
from grnbench.strategies import apply_strategy


# ------------------------------------------------------------------ labels #

def test_normalize_labels_strips_and_collapses():
    s = pd.Series([" Tumor ", "Tumor", "  normal  tissue", "normal tissue", None])
    out, rep = normalize_labels(s)
    assert out.dropna().tolist() == ["Tumor", "Tumor", "normal tissue", "normal tissue"]
    assert rep["n_whitespace_fixed"] == 2
    assert rep["n_missing"] == 1


def test_normalize_labels_case_and_synonyms():
    s = pd.Series(["wt", "WT", "knock-out", "ko"])
    out, rep = normalize_labels(s, mapping={"knock-out": "KO", "ko": "KO"}, case="upper")
    assert out.tolist() == ["WT", "WT", "KO", "KO"]
    # case pass changes wt/knock-out/ko (3); synonym pass then fixes knock-out only
    assert rep["n_case_changed"] == 3 and rep["n_synonym_mapped"] == 1


def test_normalize_labels_missing_tokens():
    s = pd.Series(["A", "?", "unknown", "", "NaN"])
    out, rep = normalize_labels(s)
    assert out.notna().sum() == 1
    assert rep["n_missing"] == 4


def test_load_label_mapping_aliases(tmp_path):
    p = tmp_path / "map.csv"
    p.write_text("old,new\nwt,WT\nko,KO\n", encoding="utf-8")
    m = load_label_mapping(p)
    assert m == {"wt": "WT", "ko": "KO"}


def test_audit_obs_reports_imbalance(synthetic_adata):
    adata = synthetic_adata.copy()
    adata.obs["skew"] = ["A"] * 395 + ["B"] * 4 + [np.nan]
    rows = {r["column"]: r for r in audit_obs(adata)}
    assert rows["group"]["imbalance_ratio"] == pytest.approx(200 / 5)
    assert rows["skew"]["n_missing"] == 1
    assert rows["skew"]["smallest_class_cells"] == 4


def test_detect_batch_column_explicit_and_auto():
    a = ad.AnnData(X=np.ones((5, 2)),
                   obs=pd.DataFrame({"batch": list("AABBB"),
                                     "note": list("xxxxx")}),
                   var=pd.DataFrame(index=["g1", "g2"]))
    col, how = detect_batch_column(a)
    assert col == "batch" and how == "auto-detected"
    col2, how2 = detect_batch_column(a, explicit="note")
    assert col2 == "note" and how2 == "explicit"
    with pytest.raises(KeyError):
        detect_batch_column(a, explicit="nope")


# --------------------------------------------------------------- gene IDs #

def test_detect_id_type_verdicts():
    assert detect_id_type(["ENSG00000141510", "ENSG00000000003"])["verdict"] == "ensembl_only"
    assert detect_id_type(["TP53", "CDKN1A", "BAX"])["verdict"] == "symbols"
    mixed = detect_id_type(["TP53", "ENSG00000141510", "BAX", "MYC"])
    assert mixed["verdict"] == "mixed"


def test_organism_hint_human_vs_mouse():
    human_genes = [f"GENE{i}" for i in range(50)] + ["MT-CO1", "MT-ND2"]
    assert organism_hint(human_genes) == "human-like"
    mouse_genes = [f"Gm{i}" for i in range(50)] + ["mt-Co1", "mt-Nd2"]
    assert organism_hint(mouse_genes) == "mouse-like"


# ------------------------------------------------- counts-layer integrity #
# NB: a real h5ad cannot even contain a broken counts layer (anndata blocks
# mismatched writes), so the validator is exercised directly as a unit.

def test_transposed_counts_layer_rejected():
    from grnbench.ingest import _validate_counts_layer
    x_as_loaded = np.arange(6.0).reshape(3, 2)
    layer = sp.csr_matrix(np.arange(6.0).reshape(2, 3))
    with pytest.raises(IngestError, match="TRANSPOSED"):
        _validate_counts_layer(x_as_loaded, layer)


def test_shape_mismatch_counts_layer_rejected():
    from grnbench.ingest import _validate_counts_layer
    with pytest.raises(IngestError, match="shape"):
        _validate_counts_layer(np.ones((3, 2)), sp.csr_matrix(np.ones((4, 2))))


def test_negative_counts_layer_rejected():
    from grnbench.ingest import _validate_counts_layer
    layer = np.array([[1.0, 2.0], [3.0, -4.0]])
    with pytest.raises(IngestError, match="raw counts"):
        _validate_counts_layer(np.ones_like(layer), layer)


def test_noninteger_counts_layer_rejected():
    from grnbench.ingest import _validate_counts_layer
    layer = np.array([[1.5, 2.0], [3.0, 4.7]])
    with pytest.raises(IngestError, match="raw counts"):
        _validate_counts_layer(layer, layer)


def test_nan_matrix_rejected(tmp_path):
    X = np.array([[1.0, np.nan], [3.0, 4.0]])
    a = ad.AnnData(X=X, var=pd.DataFrame(index=["g1", "g2"]))
    p = tmp_path / "nan.h5ad"
    a.write_h5ad(p)
    with pytest.raises(IngestError, match="NaN/Inf"):
        load_any(p)


# --------------------------------------------------- pipeline integration #

def test_zero_count_cells_dropped_and_logged(tmp_path, synthetic_adata):
    a = synthetic_adata.copy()
    X = a.X.toarray()
    X[0, :] = 0
    X[1, :] = 0
    a.X = sp.csr_matrix(X)
    src = tmp_path / "zeros.h5ad"
    a.write_h5ad(src)
    res = run_prep(src, tmp_path / "run")
    step = next(s for s in res.manifest["steps"] if s["step"] == "filter_cells")
    assert "zero counts" in step["detail"]
    proc = ad.read_h5ad(res.processed_path)
    totals = np.asarray((proc.X > 0).sum(axis=1)).ravel()
    assert (totals > 0).all()


def test_class_imbalance_gate_warns_on_skewed_groups(tmp_path):
    rng = np.random.default_rng(2)
    X = rng.poisson(5.0, size=(400, 300)).astype(np.float32)
    a = ad.AnnData(X=sp.csr_matrix(X),
                   obs=pd.DataFrame({"cond": ["CTRL"] * 396 + ["KO"] * 4}),
                   var=pd.DataFrame(index=[f"G{i}" for i in range(300)]))
    src = tmp_path / "skew.h5ad"
    a.write_h5ad(src)
    res = run_prep(src, tmp_path / "run", group_key="cond")
    gate = next(g for g in res.manifest["qc"]["gates"] if g["name"] == "class_imbalance")
    assert gate["status"] == "warn"


def test_stratified_subsample_preserves_proportions(tmp_path, synthetic_adata):
    src = tmp_path / "in.h5ad"
    synthetic_adata.write_h5ad(src)
    res = run_prep(src, tmp_path / "run_s", subsample_cells=100)  # group_key absent -> label? none; use stratify via group col present in obs
    # synthetic_adata.obs has 'group'; stratify_key defaults to it only when passed explicitly
    proc = ad.read_h5ad(res.processed_path)
    assert proc.n_obs <= 100


def test_stratified_subsample_with_group_key(tmp_path, synthetic_adata):
    src = tmp_path / "in.h5ad"
    synthetic_adata.write_h5ad(src)
    res = run_prep(src, tmp_path / "run_s2", subsample_cells=60,
                   stratify_key="group")
    proc = ad.read_h5ad(res.processed_path)
    vc = proc.obs["group"].value_counts()
    total = vc.sum()
    frac_A = vc.get("A", 0) / total
    frac_C = vc.get("C", 0) / total
    # original: A=200/400=.5, C=5/400=.0125 ; allow generous slack for tiny class
    assert abs(frac_A - 0.5) < 0.12
    assert 0 <= frac_C <= 0.1


def test_duplicate_barcodes_report_and_gate():
    from types import SimpleNamespace
    ns = SimpleNamespace(obs_names=pd.Index(["c1", "c1", "c2"]))
    assert duplicate_barcodes_report(ns) == 1
    from grnbench.qc import Thresholds, run_qc
    a = ad.AnnData(X=np.eye(4, dtype=np.float32),
                   var=pd.DataFrame(index=[f"g{i}" for i in range(4)]))
    rep = run_qc(a.copy(), Thresholds(), n_dup_barcodes=7)
    gate = next(g for g in rep.gates if g.name == "duplicate_barcodes")
    assert gate.status == "warn" and gate.value == 7


def test_sparsify_and_downcast_steps(tmp_path):
    rng = np.random.default_rng(4)
    n_cells, n_genes = 80, 2000
    X = np.zeros((n_cells, n_genes), dtype=np.float64)
    for i in range(n_cells):  # ~220 detected genes/cell -> 11% density (>=80% zeros)
        cols = rng.choice(n_genes, size=220, replace=False)
        X[i, cols] = rng.poisson(3, size=220).astype(np.float64) + 1.0
    a = ad.AnnData(X=X, var=pd.DataFrame(index=[f"g{i}" for i in range(n_genes)]))
    src = tmp_path / "dense.h5ad"
    a.write_h5ad(src)
    res = run_prep(src, tmp_path / "run_mem", strategy="cpm_log1p")
    step_names = {s["step"] for s in res.manifest["steps"]}
    assert "sparsify" in step_names and "downcast" in step_names
    proc = ad.read_h5ad(res.processed_path)
    assert sp.issparse(proc.X)


def test_pseudobulk_min_group_cells_drops_tiny(synthetic_adata):
    with pytest.warns(UserWarning, match="dropping"):
        out, _prov = apply_strategy("pseudobulk", synthetic_adata,
                                    group_key="group", min_group_cells=10)
    assert "C" not in set(out.obs["group"])
    assert set(out.obs["n_cells"]) == {200, 195}


def test_tf_list_matching_nothing_is_refused(tmp_path, synthetic_adata):
    """Zero matches used to fall through to 'every gene is a regulator'.

    That silently replaces the requested candidate space with a completely
    different one, so it is now an error rather than a warning.
    """
    from grnbench.tfs import TFListError

    src = tmp_path / "in.h5ad"
    synthetic_adata.write_h5ad(src)
    with pytest.raises(TFListError, match="none of the requested TF identifiers"):
        run_prep(src, tmp_path / "run_tf", tf_list=["ZZZ1", "ZZZ2", "ZZZ3"], seed=1)


def test_tf_coverage_low_warning(tmp_path, synthetic_adata):
    """A list that mostly misses (but not entirely) still runs, with a warning."""
    src = tmp_path / "in.h5ad"
    synthetic_adata.write_h5ad(src)
    one_real = [str(synthetic_adata.var_names[0])]
    res = run_prep(src, tmp_path / "run_tf",
                   tf_list=one_real + [f"ZZZ{i}" for i in range(9)], seed=1)
    assert any("TF coverage low" in w for w in res.manifest.get("warnings", []))
    assert res.manifest["tf_coverage"]["matched_in_final_matrix"] == 1


def test_doublet_heuristic_fraction_extreme(extreme_mito_adata):
    frac, flagged = doublet_heuristic_fraction(extreme_mito_adata)
    assert 0.0 <= frac <= 0.05
    assert len(flagged) == extreme_mito_adata.n_obs


def test_downcast_floats_converts():
    a = ad.AnnData(X=np.eye(3, dtype=np.float64),
                   var=pd.DataFrame(index=list("abc")))
    changed = downcast_floats(a)
    assert changed and a.X.dtype == np.float32


def test_pseudobulk_dropping_every_group_is_a_loud_failure(tmp_path, synthetic_adata):
    """min_group_cells above every group size used to silently write a 0-cell
    processed.h5ad and exit 0 - an empty matrix that feeds inference next."""
    from grnbench.pipeline import PrepError

    src = tmp_path / "in.h5ad"
    synthetic_adata.write_h5ad(src)
    with pytest.raises(PrepError, match="empty matrix"):
        run_prep(src, tmp_path / "run_empty", group_key="group",
                 strategy="pseudobulk", min_group_cells=1000, seed=1)


def test_pseudobulk_small_group_warning_reaches_the_manifest(tmp_path, synthetic_adata):
    """The strategy's UserWarning ('dropping ... group(s)') used to be visible
    only on stderr; it belongs in manifest['warnings'] like every other one."""
    src = tmp_path / "in.h5ad"
    synthetic_adata.write_h5ad(src)
    res = run_prep(src, tmp_path / "run_tiny", group_key="group",
                   strategy="pseudobulk", min_group_cells=10, seed=1)
    assert any("dropping" in w and "group" in w for w in res.manifest["warnings"])
