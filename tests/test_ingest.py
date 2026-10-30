"""Ingest tests: format loading, orientation detection, counts integrity."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import anndata as ad
import scipy.sparse as sp

from grnbench.ingest import (
    IngestError,
    IngestReport,
    _clean_gene_names,
    load_gene_mapping,
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
    first = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    a = ad.AnnData(
        X=np.hstack([first, np.zeros((2, 999))]),
        var=pd.DataFrame(index=["LPCAT4", "Lpcat4", "Gata1",
                                *[f"Gene{i}" for i in range(999)]]),
    )
    rep = IngestReport(path="test", fmt="synthetic")

    a = _clean_gene_names(a, rep)

    assert list(a.var_names[:2]) == ["LPCAT4", "Gata1"]
    dense = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
    assert np.allclose(dense[:, :2], [[3.0, 3.0], [9.0, 6.0]])
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


def test_human_and_mouse_ensembl_version_suffixes_are_stripped_consistently():
    a = ad.AnnData(
        X=np.ones((1, 12)),
        var=pd.DataFrame(index=[
            "ENSG00000141510.12", "ENSMUSG00000059552.7",
            *[f"GENE{i}" for i in range(10)],
        ]),
    )
    rep = IngestReport(path="test", fmt="synthetic")
    a = _clean_gene_names(a, rep)
    assert list(a.var_names[:2]) == ["ENSG00000141510", "ENSMUSG00000059552"]
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


def test_ensembl_mapping_chooses_the_input_column_with_most_usable_symbols():
    ids = [f"ENSG{i:011d}" for i in range(20)]
    sparse = [f"OLD{i}" if i < 3 else "" for i in range(20)]
    complete = [f"GENE{i}" for i in range(20)]
    a = ad.AnnData(
        X=np.ones((1, len(ids))),
        var=pd.DataFrame(
            {"feature_name": sparse, "gene_symbols": complete}, index=ids),
    )
    rep = IngestReport(path="test", fmt="synthetic")

    out = _clean_gene_names(a, rep)

    assert list(out.var_names) == complete
    assert rep.identifier_provenance["mapping_source"] == {
        "kind": "input_column", "column": "gene_symbols"}
    assert rep.identifier_provenance["mapping_outcome"]["mapped_features"] == 20


def test_ensembl_dominant_input_without_symbol_source_is_refused():
    names = [f"ENSG{i:011d}" for i in range(19)] + ["TP53"]
    a = ad.AnnData(X=np.ones((1, len(names))), var=pd.DataFrame(index=names))
    rep = IngestReport(path="test", fmt="synthetic")

    with pytest.raises(IngestError, match=r"more than 90%.*--gene-mapping"):
        _clean_gene_names(a, rep)


def test_user_gene_mapping_is_hashed_applied_and_fully_recorded(tmp_path):
    mapping_path = tmp_path / "ids.tsv"
    mapping_path.write_text(
        "id\tsymbol\nENSG00000141510\tTP53\nENSMUSG00000059552\tTrp53\n",
        encoding="utf-8",
    )
    mapping, metadata = load_gene_mapping(mapping_path)
    a = ad.AnnData(
        X=np.array([[2.0, 3.0]]),
        var=pd.DataFrame(index=["ENSG00000141510.18", "ENSMUSG00000059552.7"]),
    )
    rep = IngestReport(path="test", fmt="synthetic")

    a = _clean_gene_names(a, rep, gene_mapping=mapping, mapping_source=metadata)

    assert list(a.var_names) == ["TP53", "Trp53"]
    provenance = rep.identifier_provenance
    assert provenance["mapping_source"] == {
        "kind": "user_table",
        "sha256": metadata["sha256"],
        "rows": 2,
        "unique_ids": 2,
    }
    assert provenance["mapping_outcome"]["mapped_features"] == 2
    assert provenance["source_to_representative"] == [
        {
            "as_loaded": "ENSG00000141510.18", "representative": "TP53",
            "occurrences": 1, "transformations": ["user_mapping"],
        },
        {
            "as_loaded": "ENSMUSG00000059552.7", "representative": "Trp53",
            "occurrences": 1, "transformations": ["user_mapping"],
        },
    ]


def test_conflicting_user_gene_mapping_is_refused(tmp_path):
    mapping_path = tmp_path / "bad.tsv"
    mapping_path.write_text("id\tsymbol\nENSG1\tTP53\nENSG1\tTrp53\n", encoding="utf-8")
    with pytest.raises(IngestError, match="maps one identifier to multiple symbols"):
        load_gene_mapping(mapping_path)


def test_paul15_sized_single_case_group_is_accepted_and_recorded():
    names = ["LPCAT4", "Lpcat4", *[f"Gene{i}" for i in range(3449)]]
    a = ad.AnnData(X=np.ones((1, len(names))), var=pd.DataFrame(index=names))
    rep = IngestReport(path="test", fmt="synthetic")

    out = _clean_gene_names(a, rep)

    assert out.n_vars == 3450
    assert rep.identifier_provenance["case_fold_groups"] == [{
        "representative": "LPCAT4",
        "spellings": ["LPCAT4", "Lpcat4"],
        "method": "sum",
        "cross_convention": True,
    }]


def test_many_exact_duplicates_are_accepted_summed_and_recorded():
    names = ["Y_RNA"] * 733 + [f"GENE{i}" for i in range(51_976 - 733)]
    a = ad.AnnData(X=np.ones((1, len(names))), var=pd.DataFrame(index=names))
    rep = IngestReport(path="test", fmt="synthetic")

    out = _clean_gene_names(a, rep)

    assert out.n_vars == 51_244
    assert float(out[:, "Y_RNA"].X[0, 0]) == 733.0
    assert rep.identifier_provenance["exact_duplicate_groups"] == [{
        "representative": "Y_RNA",
        "occurrences": 733,
        "method": "sum",
    }]
    assert rep.identifier_provenance["collision_summary"]["case_fold_groups"] == 0


def test_tp53_and_trp53_are_not_a_case_fold_collision():
    a = ad.AnnData(X=np.ones((1, 2)), var=pd.DataFrame(index=["TP53", "Trp53"]))
    rep = IngestReport(path="test", fmt="synthetic")
    out = _clean_gene_names(a, rep)
    assert list(out.var_names) == ["TP53", "Trp53"]
    assert rep.identifier_provenance["case_fold_groups"] == []


def test_human_mouse_mixed_symbol_fixture_is_refused():
    human = ["ACTB", "GAPDH", "MYC", "SOX2", "JUN", "FOS"]
    mouse = ["Actb", "Gapdh", "Myc", "Sox2", "Jun", "Fos"]
    a = ad.AnnData(
        X=np.ones((1, 12)),
        obs=pd.DataFrame({"organism": ["Homo sapiens"]}, index=["cell"]),
        var=pd.DataFrame(index=human + mouse),
    )
    with pytest.raises(IngestError, match="mixed human/mouse annotation"):
        _clean_gene_names(a, IngestReport(path="test", fmt="synthetic"))


def test_measured_orf_token_does_not_create_a_false_cross_convention_group():
    names = ["C1orf112", "C1ORF112", *[f"GENE{i}" for i in range(1_000)]]
    a = ad.AnnData(X=np.ones((1, len(names))), var=pd.DataFrame(index=names))
    rep = IngestReport(path="test", fmt="synthetic")
    _clean_gene_names(a, rep)
    assert rep.identifier_provenance["case_fold_groups"][0]["cross_convention"] is False


def _cross_convention_fixture(n_genes, n_groups):
    names = []
    for i in range(n_groups):
        names.extend([f"GENE{i}", f"Gene{i}"])
    names.extend(f"FILLER{i}" for i in range(n_genes - len(names)))
    return ad.AnnData(
        X=np.ones((1, n_genes)),
        obs=pd.DataFrame({"organism": ["Homo sapiens"]}, index=["cell"]),
        var=pd.DataFrame(index=names),
    )


def test_cross_convention_boundary_accepts_two_of_2000_and_refuses_three():
    accepted = _cross_convention_fixture(2_000, 2)
    rep = IngestReport(path="test", fmt="synthetic")
    _clean_gene_names(accepted, rep)
    assert rep.identifier_provenance["collision_summary"]["cross_convention_groups"] == 2

    refused = _cross_convention_fixture(2_000, 3)
    with pytest.raises(IngestError, match=r"3 cross-convention.*0.1%"):
        _clean_gene_names(refused, IngestReport(path="test", fmt="synthetic"))


def test_cross_convention_boundary_has_a_two_group_floor_for_small_panels():
    accepted = _cross_convention_fixture(200, 2)
    rep = IngestReport(path="test", fmt="synthetic")
    _clean_gene_names(accepted, rep)
    summary = rep.identifier_provenance["collision_summary"]
    assert summary["cross_convention_groups"] == 2
    assert summary["cross_convention_min_groups"] == 2

    refused = _cross_convention_fixture(200, 3)
    with pytest.raises(IngestError, match=r"3 cross-convention.*max\(2 groups, 0.1%\)"):
        _clean_gene_names(refused, IngestReport(path="test", fmt="synthetic"))


def test_distinct_spelling_case_groups_over_one_percent_are_refused():
    names = ["Abc", "aBc", "Def", "dEf", *[f"GENE{i}" for i in range(96)]]
    a = ad.AnnData(
        X=np.ones((1, 100)),
        obs=pd.DataFrame({"organism": ["Homo sapiens"]}, index=["cell"]),
        var=pd.DataFrame(index=names),
    )
    with pytest.raises(IngestError, match=r"2 distinct-spelling.*1%"):
        _clean_gene_names(a, IngestReport(path="test", fmt="synthetic"))


def test_third_organism_preserves_case_variants_and_sums_only_exact_duplicates():
    names = ["NTM", "ntm", "zgc:1", "zgc:1", *[f"gene{i}" for i in range(100)]]
    a = ad.AnnData(
        X=np.ones((1, len(names))),
        obs=pd.DataFrame({"organism": ["Danio rerio"]}, index=["cell"]),
        var=pd.DataFrame(index=names),
    )
    rep = IngestReport(path="test", fmt="synthetic")

    out = _clean_gene_names(a, rep)

    assert list(out.var_names[:3]) == ["NTM", "ntm", "zgc:1"]
    assert float(out[:, "zgc:1"].X[0, 0]) == 2.0
    assert rep.n_duplicates_merged == 1
    assert rep.n_case_duplicates_merged == 0
    assert rep.identifier_provenance["policy"] == "D-028"
    assert rep.identifier_provenance["organism_policy"] == {
        "scope": "other",
        "label": "danio rerio",
        "source": "metadata:obs['organism']",
        "case_fold_action": "preserve_and_warn",
    }
    assert rep.identifier_provenance["case_fold_groups"] == [{
        "representative": None,
        "spellings": ["NTM", "ntm"],
        "method": "preserve",
        "cross_convention": True,
    }]
    assert any("preserved" in warning for warning in rep.warnings)


def test_zebrafish_scale_case_groups_are_preserved_instead_of_refused():
    symbols = []
    for i in range(91):
        symbols.extend([f"ZGC{i}", f"zgc{i}"])
    symbols.extend(f"gene{i}" for i in range(32_031 - len(symbols)))
    identifiers = [f"ENSDARG{i:011d}" for i in range(32_031)]
    a = ad.AnnData(
        X=sp.csr_matrix(np.ones((1, len(symbols)))),
        var=pd.DataFrame({"feature_name": symbols}, index=identifiers),
    )
    rep = IngestReport(path="test", fmt="synthetic")

    out = _clean_gene_names(a, rep)

    assert out.n_vars == 32_031
    assert rep.identifier_provenance["mapping_source"] == {
        "kind": "input_column", "column": "feature_name"}
    assert rep.identifier_provenance["organism_policy"] == {
        "scope": "other",
        "label": "non-human/mouse Ensembl",
        "source": "gene_identifiers",
        "case_fold_action": "preserve_and_warn",
    }
    assert rep.identifier_provenance["collision_summary"]["cross_convention_groups"] == 91
    assert rep.identifier_provenance["collision_summary"]["case_fold_action"] == \
        "preserve_and_warn"
    assert rep.n_case_duplicates_merged == 0


def test_unknown_organism_preserves_case_variants_conservatively():
    a = ad.AnnData(X=np.ones((1, 2)), var=pd.DataFrame(index=["ABC", "abc"]))
    rep = IngestReport(path="test", fmt="synthetic")

    out = _clean_gene_names(a, rep)

    assert list(out.var_names) == ["ABC", "abc"]
    assert rep.identifier_provenance["organism_policy"]["scope"] == "unknown"
    assert rep.identifier_provenance["collision_summary"]["case_fold_action"] == \
        "preserve_and_warn"


def test_human_metadata_keeps_casefold_merge_and_guard_policy():
    names = ["ABC", "abc", *[f"FILLER{i}" for i in range(998)]]
    a = ad.AnnData(
        X=np.array([[2.0, 3.0, *([1.0] * 998)]]),
        obs=pd.DataFrame({"organism": ["NCBITaxon:9606"]}, index=["cell"]),
        var=pd.DataFrame(index=names),
    )
    rep = IngestReport(path="test", fmt="synthetic")

    out = _clean_gene_names(a, rep)

    assert out.var_names[0] == "ABC"
    assert out.n_vars == 999
    assert float(out.X[0, 0]) == 5.0
    assert rep.identifier_provenance["organism_policy"]["scope"] == "human"
    assert rep.identifier_provenance["collision_summary"]["case_fold_action"] == \
        "merge_and_guard"


# --------------------------------------------------------------------------- #
# defects found by the end-to-end test of 2026-10-09 (issue #74)
# --------------------------------------------------------------------------- #

def test_a_repeated_gene_in_a_cells_by_genes_header_is_aggregated_like_every_other_format(tmp_path):
    """pandas renames the second 'A' to 'A.1'; the two copies must be summed, as for h5ad, MTX and genes-x-cells."""
    p = tmp_path / "dup.csv"
    p.write_text("cell,A,B,A,C\nc1,1,2,3,4\nc2,5,6,7,8\n")
    adata, rep = load_any(p, orientation="cells_x_genes")
    assert list(adata.var_names) == ["A", "B", "C"]
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    assert X.tolist() == [[4.0, 2.0, 4.0], [12.0, 6.0, 8.0]]
    assert rep.n_duplicates_merged == 1 and any("duplicate gene names aggregated" in n for n in rep.notes)


def test_loom_without_loompy_says_what_to_install(tmp_path, monkeypatch):
    import sys
    p = tmp_path / "x.loom"
    p.write_bytes(b"not a real loom, the import fails first")
    monkeypatch.setitem(sys.modules, "loompy", None)
    with pytest.raises(IngestError, match="loompy"):
        load_any(p)


@pytest.mark.parametrize("name,payload,fragment", [
    ("empty.h5ad", b"", "empty (0 bytes)"),
    ("broken.h5ad", b"this is not an hdf5 file at all" * 20, "not a valid .h5ad file"),
    ("binary.csv", b"\xc3\x28\xff\xfe\x00\x81" * 50, "binary file"),
])
def test_an_unreadable_file_gets_a_message_naming_the_file_not_library_jargon(tmp_path, name, payload, fragment):
    p = tmp_path / name
    p.write_bytes(payload)
    import re
    with pytest.raises(IngestError, match=re.escape(fragment)) as exc:
        load_any(p)
    assert name in str(exc.value)


def test_a_folder_that_is_not_a_10x_matrix_says_which_files_it_needs(tmp_path):
    (tmp_path / "somefolder").mkdir()
    with pytest.raises(IngestError, match="matrix.mtx"):
        load_any(tmp_path / "somefolder")


def test_gene_names_with_stray_whitespace_are_stripped_and_reported():
    """' TP53' matches no TF list: the TF coverage dropped without a word (HI-12)."""
    from grnbench.ingest import IngestReport, _clean_gene_names
    a = ad.AnnData(X=np.array([[1.0, 2.0, 3.0]], dtype=np.float32),
                   var=pd.DataFrame(index=[" TP53", "MYC\xa0", "STAT1"]))
    rep = IngestReport(path="test", fmt="synthetic")
    a = _clean_gene_names(a, rep)
    assert list(a.var_names) == ["TP53", "MYC", "STAT1"]
    assert any("2 gene name(s) had leading or trailing whitespace" in n for n in rep.notes), rep.notes
    clean = ad.AnnData(X=np.array([[1.0, 2.0]], dtype=np.float32), var=pd.DataFrame(index=["A", "B"]))
    rep2 = IngestReport(path="test", fmt="synthetic")
    _clean_gene_names(clean, rep2)
    assert not any("whitespace" in n for n in rep2.notes)
