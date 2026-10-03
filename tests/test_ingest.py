"""Ingest tests: format loading, orientation detection, counts integrity."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import anndata as ad

from grnbench.ingest import (
    IngestError,
    IngestReport,
    _clean_gene_names,
    load_any,
    looks_like_counts,
)


def test_load_h5ad_roundtrip(tmp_path, synthetic_adata):
    p = tmp_path / "syn.h5ad"
    synthetic_adata.write_h5ad(p)
    adata, rep = load_any(p)
    assert adata.shape == synthetic_adata.shape
    assert "counts" in adata.layers
    ok, why = looks_like_counts(adata.X)
    assert ok, why
    assert rep.fmt == "h5ad"


def test_csv_orientation_detection_both_ways(tmp_path, synthetic_adata):
    df = pd.DataFrame(
        synthetic_adata.X.toarray(),
        index=synthetic_adata.obs_names.astype(str),
        columns=synthetic_adata.var_names.astype(str),
    )
    p_cells_x_genes = tmp_path / "cxg.csv"
    df.to_csv(p_cells_x_genes)
    a1, r1 = load_any(p_cells_x_genes)
    assert a1.shape == (400, 600)
    assert not r1.transposed

    p_genes_x_cells = tmp_path / "gxc.csv"
    df.T.to_csv(p_genes_x_cells)
    a2, r2 = load_any(p_genes_x_cells)
    assert a2.shape == (400, 600)
    assert r2.transposed


def test_rejects_normalized_input_without_counts_layer(tmp_path, synthetic_adata):
    from grnbench.strategies import cpm_log1p
    norm = cpm_log1p(synthetic_adata)
    del norm.layers["counts"]  # simulate a dataset saved post-normalization
    p = tmp_path / "norm.h5ad"
    norm.write_h5ad(p)
    with pytest.raises(IngestError, match="counts"):
        load_any(p)


def test_rejects_negative_scaled_input():
    X = np.array([[1.0, -0.5], [0.3, 2.0]])
    ok, why = looks_like_counts(X)
    assert not ok and "negative" in why


def test_duplicate_genes_aggregated_by_sum():
    a = ad.AnnData(
        X=np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),
        var=pd.DataFrame(index=["v0", "v1", "v2"]),  # unique at construction
    )
    a.var.index = pd.Index(["A", "A", "B"])  # introduce duplicates post-hoc
    rep = IngestReport(path="test", fmt="synthetic")
    a = _clean_gene_names(a, rep)
    assert list(a.var_names) == ["A", "B"]
    dense = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
    assert np.allclose(dense, [[3.0, 3.0], [9.0, 6.0]])
    assert rep.n_duplicates_merged == 1


def test_case_insensitive_duplicate_genes_are_aggregated_and_reported():
    """Regression for #30: case-fold collisions must be resolved before TF matching."""
    a = ad.AnnData(
        X=np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),
        var=pd.DataFrame(index=["LPCAT4", "Lpcat4", "Gata1"]),
    )
    rep = IngestReport(path="test", fmt="synthetic")

    a = _clean_gene_names(a, rep)

    assert list(a.var_names) == ["LPCAT4", "Gata1"]
    dense = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
    assert np.allclose(dense, [[3.0, 3.0], [9.0, 6.0]])
    assert rep.n_case_duplicates_merged == 1
    assert any("case-insensitive" in note for note in rep.notes)


def _write_10x_h5(path, symbols, ids, dense_genes_x_cells):
    """Minimal CellRanger v3 .h5 (CSC, features x barcodes) for loader tests."""
    h5py = pytest.importorskip("h5py")
    from scipy import sparse

    m = sparse.csc_matrix(np.asarray(dense_genes_x_cells, dtype=np.int32))
    n_genes, n_cells = m.shape
    with h5py.File(path, "w") as f:
        g = f.create_group("matrix")
        g.create_dataset("data", data=m.data.astype(np.int32))
        g.create_dataset("indices", data=m.indices.astype(np.int64))
        g.create_dataset("indptr", data=m.indptr.astype(np.int64))
        g.create_dataset("shape", data=np.array([n_genes, n_cells], dtype=np.int32))
        g.create_dataset("barcodes",
                         data=np.array([f"CELL{i}".encode() for i in range(n_cells)]))
        ft = g.create_group("features")
        ft.create_dataset("id", data=np.array([i.encode() for i in ids]))
        ft.create_dataset("name", data=np.array([s.encode() for s in symbols]))
        ft.create_dataset("feature_type",
                          data=np.array([b"Gene Expression"] * n_genes))
        ft.create_dataset("genome", data=np.array([b"GRCh38"] * n_genes))
    return path


def test_10x_h5_duplicate_symbols_are_aggregated_like_the_mtx_path(tmp_path):
    """Regression: the .h5 loader used to call var_names_make_unique() first,
    which renamed TBCE -> TBCE-1 and left _clean_gene_names nothing to merge -
    the same biological data then gave a different matrix than the MTX path,
    and the manifest recorded n_duplicates_merged=0."""
    p = _write_10x_h5(
        tmp_path / "filtered_feature_bc_matrix.h5",
        symbols=["TBCE", "TBCE", "ACTB"],
        ids=["ENSG00000284770", "ENSG00000285053", "ENSG00000075624"],
        dense_genes_x_cells=[[1, 4], [2, 5], [3, 6]],
    )
    adata, rep = load_any(p)

    assert list(adata.var_names) == ["TBCE", "ACTB"], "duplicates must be merged"
    dense = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    assert np.allclose(dense, [[3.0, 3.0], [9.0, 6.0]]), "counts must be summed"
    assert rep.n_duplicates_merged == 1
    assert rep.fmt == "10x_h5"
    assert adata.layers["counts"].shape == adata.X.shape


def test_10x_h5_without_duplicates_is_untouched(tmp_path):
    p = _write_10x_h5(
        tmp_path / "clean.h5",
        symbols=["TP53", "ACTB"], ids=["ENSG1", "ENSG2"],
        dense_genes_x_cells=[[1, 2], [3, 4]],
    )
    adata, rep = load_any(p)
    assert list(adata.var_names) == ["TP53", "ACTB"]
    assert rep.n_duplicates_merged == 0


def test_ensembl_version_suffixes_stripped():
    a = ad.AnnData(
        X=np.array([[1.0, 2.0]]),
        var=pd.DataFrame(index=["ENSG00000141510.12", "ENSG00000000003.19"]),
    )
    rep = IngestReport(path="test", fmt="synthetic")
    a = _clean_gene_names(a, rep)
    assert list(a.var_names) == ["ENSG00000141510", "ENSG00000000003"]
    assert rep.version_suffixes_stripped == 2


# --------------------------------------------------------------------------- #
# orientation: 10x barcodes match the gene-symbol regex, so gene-likeness alone
# could not tell a GEO-layout matrix (genes as rows) from a cells-as-rows one.
# The matrix went in transposed and every QC metric was computed on the wrong
# axis, including n_mito_genes = 0, which silently disables the mito gate.
# --------------------------------------------------------------------------- #

def _geo_layout_csv(tmp_path, n_genes=300, n_cells=80, seed=0):
    rng = np.random.default_rng(seed)
    genes = ([f"GENE{i}" for i in range(n_genes - 4)]
             + ["MT-CO1", "MT-ND1", "RPS6", "ACTB"])
    barcodes = ["".join(rng.choice(list("ACGT"), 16)) + "-1" for _ in range(n_cells)]
    X = rng.poisson(1.5, size=(n_genes, n_cells))
    p = tmp_path / "geo.csv"
    pd.DataFrame(X, index=genes, columns=barcodes).to_csv(p)
    return p, n_genes, n_cells


def test_barcode_columns_identify_the_cell_axis(tmp_path):
    p, n_genes, n_cells = _geo_layout_csv(tmp_path)
    adata, rep = load_any(p)
    assert rep.transposed is True
    assert adata.shape == (n_cells, n_genes)
    assert "MT-CO1" in set(adata.var_names)      # genes really are on the var axis


def test_orientation_can_be_declared_and_is_obeyed(tmp_path):
    p, n_genes, n_cells = _geo_layout_csv(tmp_path)
    a, rep = load_any(p, orientation="genes_x_cells")
    assert a.shape == (n_cells, n_genes) and rep.transposed
    # the override must be able to *disagree* with the detector, not just echo it
    b, rep_b = load_any(p, orientation="cells_x_genes")
    assert b.shape == (n_genes, n_cells) and not rep_b.transposed
    with pytest.raises(IngestError, match="orientation must be"):
        load_any(p, orientation="sideways")


def test_ambiguous_orientation_is_refused_not_guessed(tmp_path):
    """Neither axis looks like barcodes or carries markers: refuse to guess."""
    names_a = [f"AAA{i}" for i in range(40)]
    names_b = [f"BBB{i}" for i in range(40)]
    p = tmp_path / "ambiguous.csv"
    pd.DataFrame(np.ones((40, 40), dtype=int), index=names_a, columns=names_b).to_csv(p)
    with pytest.raises(IngestError, match="cannot tell whether rows are cells"):
        load_any(p)


def test_counts_are_taken_from_raw_when_x_is_normalized(tmp_path):
    counts = np.array([[1.0, 0.0, 3.0], [2.0, 5.0, 0.0]])
    a = ad.AnnData(X=np.log1p(counts), var=pd.DataFrame(index=["A", "B", "C"]))
    a.raw = ad.AnnData(X=counts, var=pd.DataFrame(index=["A", "B", "C"]))
    p = tmp_path / "cxg.h5ad"
    a.write_h5ad(p)
    adata, rep = load_any(p)
    assert rep.counts_source == "raw.X"
    assert np.allclose(adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X, counts)


def test_ensembl_ids_are_mapped_to_symbols_from_a_var_column():
    a = ad.AnnData(
        X=np.array([[1.0, 2.0, 3.0]]),
        var=pd.DataFrame({"feature_name": ["TP53", "MYC", "MYC"]},
                         index=["ENSG00000141510", "ENSG00000136997", "ENSG00000999999"]),
    )
    rep = IngestReport(path="test", fmt="synthetic")
    a = _clean_gene_names(a, rep)
    assert list(a.var_names) == ["TP53", "MYC"]
    assert np.allclose(np.asarray(a.X.todense()), [[1.0, 5.0]])
    assert any("mapped to symbols" in n for n in rep.notes)
