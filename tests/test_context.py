"""Dataset context: facts from matrix, metadata, README files and URLs."""
import anndata as ad
import numpy as np
import pandas as pd
import pytest

from grnbench.context import (build_context, discover_text_files, extract_from_text,
                              geo_text_url, html_to_text)

GEO_TEXT = """!Series_title = Single-cell atlas of human bone marrow
!Series_summary = We profiled 5,000 cells with 10x Genomics Chromium v3.
!Series_overall_design = Raw UMI counts are provided; reads aligned to GRCh38.
!Sample_characteristics_ch1 = tissue: bone marrow
!Series_geo_accession = GSE123456
"""


def test_geo_like_text_gives_every_fact_with_its_snippet():
    ctx = extract_from_text(GEO_TEXT, "url:geo")
    s = ctx.summary()
    assert s["organism"] == "human"
    assert s["technology"] == "umi" and s["platform"] == "10x Chromium"
    assert s["value_type"] == "raw_counts"
    assert s["tissue"] == "bone marrow"
    assert s["n_cells_declared"] == "5000"
    assert s["accessions"] == ["GSE123456"]
    assert all(f.source == "url:geo" and f.snippet for f in ctx.facts)


def test_full_length_mouse_protocol():
    ctx = extract_from_text("Mus musculus cortex, Smart-seq2 full-length, TPM values.", "readme")
    s = ctx.summary()
    assert s["organism"] == "mouse"
    assert s["technology"] == "full_length" and s["platform"] == "Smart-seq2"
    assert s["value_type"] == "normalized"


def test_passing_mention_does_not_outvote_the_species_name():
    ctx = extract_from_text("Mus musculus samples; relevant to human disease.", "readme")
    assert ctx.get("organism") == "mouse"
    assert "organism" in ctx.conflicts()


def _adata(var_names, obs=None, uns=None):
    n = 4 if obs is None else len(obs)
    return ad.AnnData(X=np.ones((n, len(var_names)), dtype=np.float32),
                      obs=obs if obs is not None else pd.DataFrame(index=[f"c{i}" for i in range(n)]),
                      var=pd.DataFrame(index=var_names), uns=uns or {})


def test_matrix_gene_ids_beat_a_readme_that_disagrees():
    a = _adata([f"ENSG{i:011d}" for i in range(20)])
    ctx = build_context(adata=a, notes="Mus musculus, Mus musculus", discover=False)
    assert ctx.get("organism") == "human"
    assert ctx.conflicts()["organism"] == {"human": 10, "mouse": 6}


def test_cellxgene_style_metadata():
    obs = pd.DataFrame({
        "assay": ["10x 3' v3"] * 6, "organism": ["Homo sapiens"] * 6,
        "tissue": ["lung"] * 6, "donor_id": ["d1", "d1", "d2", "d2", "d3", "d3"],
        "cell_type": ["T", "B", "T", "B", "NK", "T"],
    }, index=[f"c{i}" for i in range(6)])
    a = _adata(["CD3E", "MS4A1", "MT-CO1"], obs=obs, uns={"title": "Lung atlas"})
    s = build_context(adata=a, discover=False).summary()
    assert s["organism"] == "human" and s["technology"] == "umi"
    assert s["tissue"] == "lung"
    assert s["donor_key"] == "donor_id" and s["cell_type_key"] == "cell_type"
    assert s["batch_key"] == "donor_id"


def test_readme_discovery_skips_tables_binaries_and_desktop_ini(tmp_path):
    data = tmp_path / "counts.h5ad"
    data.write_bytes(b"\x89HDF")
    (tmp_path / "README.md").write_text("Drop-seq of mouse retina", encoding="utf-8")
    (tmp_path / "desktop.ini").write_text("[x]", encoding="utf-8")
    header = "\t".join(f"cell{i}" for i in range(30))
    (tmp_path / "matrix.txt").write_text((header + "\n") * 5, encoding="utf-8")
    files = discover_text_files(data)
    assert data not in files and tmp_path / "desktop.ini" not in files
    ctx = build_context(data)
    assert ctx.sources_read == ["readme:README.md"]
    assert ctx.get("platform") == "Drop-seq" and ctx.get("organism") == "mouse"


def test_urls_are_fetched_and_failures_are_recorded_not_raised():
    def fetch(url):
        if "bad" in url:
            raise OSError("offline")
        return GEO_TEXT
    ctx = build_context(urls=["https://good.example/GSE1", "https://bad.example"],
                        fetch=fetch, discover=False)
    assert ctx.get("platform") == "10x Chromium"
    assert ctx.sources_failed == [{"source": "https://bad.example", "reason": "OSError: offline"}]


@pytest.mark.parametrize("given", ["GSE123456", "gse123456",
                                   "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123456"])
def test_geo_links_and_accessions_map_to_the_text_record(given):
    assert geo_text_url(given) == ("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"
                                   "acc=GSE123456&targ=self&form=text&view=brief")
    assert geo_text_url("https://zenodo.org/records/1") is None


def test_html_is_reduced_to_text():
    page = "<html><script>var human=1</script><p>Mus musculus</p><br>10x&nbsp;Genomics</html>"
    text = html_to_text(page)
    assert "var human" not in text and "Mus musculus" in text and "10x" in text


def test_curated_cell_type_column_beats_cluster_assignments():
    obs = pd.DataFrame({"initialClusterAssignments": ["c1", "c2"] * 5,
                        "cell_type": ["T", "B"] * 5}, index=[f"c{i}" for i in range(10)])
    a = _adata(["CD3E", "MS4A1", "MT-CO1"], obs=obs)
    assert build_context(adata=a, discover=False).get("cell_type_key") == "cell_type"


@pytest.mark.parametrize("label,expect", [("Macaca mulatta", "macaca mulatta"),
                                          ("Homo sapiens", "human"), ("Mus musculus", "mouse")])
def test_curated_species_beats_gene_casing(label, expect):
    obs = pd.DataFrame({"organism": [label] * 4}, index=[f"c{i}" for i in range(4)])
    a = _adata(["CD3E", "MS4A1", "MT-CO1"], obs=obs)      # human-looking symbols
    assert build_context(adata=a, discover=False).get("organism") == expect
