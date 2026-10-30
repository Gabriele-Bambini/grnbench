"""Silent behaviours made loud (issues #149 and #116): each one now says what happened, nothing else changes.

A planning fact is a note in the reason of the decision it concerns (no value and no ``plan.warnings`` entry
changes); a run-time fact is a ``run_prep`` warning (manifest, report). The one new plan warning is #116's:
the knee found no cliff on an unfiltered matrix, which never happens on the registered datasets.
"""
from __future__ import annotations

import copy
import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from grnbench.autoprep import diagnose, plan_preprocessing
from grnbench.cli import main as prep_main
from grnbench.demo import write_demo
from grnbench.ingest import IdentifierPolicyError
from grnbench.pipeline import run_prep
from grnbench.qc import Thresholds
from grnbench.simulate import Scenario, write_scenario

MITO_NOTE = "mitochondrial genes are found only by symbol (MT- prefix): none was found"


@pytest.fixture(scope="module")
def demo(tmp_path_factory):
    paths = write_demo(tmp_path_factory.mktemp("demo"), n_cells=1200, n_genes=900, n_tfs=15, targets_per_tf=10)
    return paths, paths["tfs"].read_text().split(), diagnose(paths["counts"])


# -- #149 item 1: Ensembl identifiers switch mito QC off --------------------------------------------------- #

def test_unmapped_ensembl_ids_are_refused_before_mito_qc(tmp_path):
    """D-028 supersedes the old pure-Ensembl warning: symbol references cannot score this input."""
    paths = write_scenario(tmp_path, Scenario("ens", gene_ids="ensembl", n_cells=600))
    with pytest.raises(IdentifierPolicyError, match="more than 90%.*Ensembl IDs"):
        diagnose(paths["counts"])


def test_unmapped_ensembl_ids_are_refused_before_per_batch_bounds(tmp_path):
    paths = write_scenario(tmp_path, Scenario("ens_b", gene_ids="ensembl", n_cells=600, n_batches=2))
    with pytest.raises(IdentifierPolicyError, match="more than 90%.*Ensembl IDs"):
        diagnose(paths["counts"])


def test_mixed_identifiers_without_mito_genes_give_the_real_ensembl_share(demo, tmp_path):
    """'mixed' (5-90% Ensembl) fires too: the mitochondrial genes can be among the Ensembl IDs."""
    paths, tfs, _ = demo
    diag = _with_uns(paths, tmp_path, {}, rename=lambda i, g: f"ENSG{i:011d}" if i >= 450 else g)
    assert diag.id_verdict["verdict"] == "mixed" and diag.qc.metrics["n_mito_genes"] == 0
    reason = plan_preprocessing(diag, tf_list=tfs, methods=[]).decisions["max_pct_mito_cell"].reason
    assert "50% of the gene identifiers are Ensembl IDs" in reason and MITO_NOTE in reason, reason


def test_symbols_with_mito_genes_get_no_mito_note(demo):
    _, tfs, diag = demo
    plan = plan_preprocessing(diag, tf_list=tfs, methods=[])
    assert all(MITO_NOTE not in d.reason for d in plan.decisions.values())


def test_mixed_identifiers_without_mito_genes_warn_at_run_time(demo, tmp_path):
    paths, _, _ = demo
    adata = ad.read_h5ad(paths["counts"])
    adata.var_names = [f"ENSG{i:011d}" if i >= 450 else g for i, g in enumerate(adata.var_names)]
    mixed = tmp_path / "mixed.h5ad"
    adata.write_h5ad(mixed)
    res = run_prep(mixed, tmp_path / "run", force=True)
    notes = [w for w in res.manifest["warnings"] if MITO_NOTE in w]
    assert len(notes) == 1 and "per-cell mito filter removes no cell" in notes[0], res.manifest["warnings"]
    assert notes[0] in (tmp_path / "run" / "report.md").read_text(encoding="utf-8")


def test_symbols_with_mito_genes_get_no_run_time_mito_note(demo, tmp_path):
    paths, _, _ = demo
    res = run_prep(paths["counts"], tmp_path / "run", force=True)
    assert not any(MITO_NOTE in w for w in res.manifest["warnings"])


# -- #149 item 3: a missing-label group that becomes a context is named -------------------------------- #

def test_a_missing_label_group_that_becomes_a_context_is_named(demo):
    _, tfs, diag = demo
    unk = copy.copy(diag)
    unk.adata = diag.adata.copy()
    unk.adata.obs["cell_type"] = unk.adata.obs["cell_type"].astype(str).replace({"beta": "unknown"})
    unk.context = copy.deepcopy(diag.context)
    n = int((unk.adata.obs["cell_type"] == "unknown").sum())
    assert n >= 200
    plan = plan_preprocessing(unk, tf_list=tfs, methods=[], choices={"context_scope": "per_cell_type"})
    assert "cell_type=unknown" in [c["name"] for c in plan.contexts], "still a context: excluding it is GB's call"
    reason = plan.decisions["context_scope"].reason
    assert f"the group 'unknown' ({n:,} cells) is a missing-label value, but it is inferred like a cell type" \
        in reason, reason
    assert not any("unknown" in w for w in plan.warnings)
    plain = plan_preprocessing(diag, tf_list=tfs, methods=[], choices={"context_scope": "per_cell_type"})
    assert "missing-label value" not in plain.decisions["context_scope"].reason


# -- #149 item 4: a context that ends below the minimum after QC ------------------------------------------- #

def _two_types(path, n_small=210, n_low=60):
    rng = np.random.default_rng(0)
    n_big, n_genes = 400, 300
    X = rng.poisson(3.0, size=(n_big + n_small, n_genes)).astype(np.float32)
    low = np.arange(n_big, n_big + n_low)                      # damaged cells of the small type: 10 genes each
    X[low] = 0
    X[np.ix_(low, np.arange(10))] = 5
    obs = pd.DataFrame({"ct": ["big"] * n_big + ["small"] * n_small},
                       index=[f"c{i}" for i in range(n_big + n_small)])
    ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=pd.DataFrame(index=[f"G{i}" for i in range(n_genes)])).write_h5ad(path)


SMALL_CONTEXT_NOTE = ("context ct == 'small' has 150 cells after cell QC (210 before), fewer than 200 cells, the "
                      "minimum the planner uses for a per-cell-type context")


def test_a_user_subset_that_ends_below_the_minimum_after_qc_gets_the_same_warning(tmp_path):
    """--subset-key/--subset-value is the user's choice, not the planner's: the warning must be true for it too
    (stress test: "fewer than the 200 cells the planner asks of a per-cell-type context")."""
    import json

    from grnbench.cli import main as cli_main
    src = tmp_path / "in.h5ad"
    _two_types(src)
    out = tmp_path / "cli"
    assert cli_main([str(src), "--out", str(out), "--subset-key", "ct", "--subset-value", "small",
                     "--min-genes-per-cell", "50", "--force"]) == 0
    warnings = json.loads((out / "manifest.json").read_text(encoding="utf-8"))["warnings"]
    assert [w for w in warnings if "after cell QC" in w] == [SMALL_CONTEXT_NOTE], warnings
    assert not any("the planner asks" in w for w in warnings)


def test_a_context_that_ends_below_the_minimum_after_qc_is_warned(tmp_path):
    src = tmp_path / "in.h5ad"
    _two_types(src)
    t = Thresholds(min_genes_per_cell=50)
    res = run_prep(src, tmp_path / "small", subset_key="ct", subset_value="small", thresholds=t, force=True)
    notes = [w for w in res.manifest["warnings"] if "after cell QC" in w]
    assert notes == [SMALL_CONTEXT_NOTE], res.manifest["warnings"]
    big = run_prep(src, tmp_path / "big", subset_key="ct", subset_value="big", thresholds=t, force=True)
    assert not any("after cell QC" in w for w in big.manifest["warnings"])


# -- #149 item 5: no --n-hvg on the command line ------------------------------------------------------------ #

def test_prep_without_n_hvg_says_all_genes_are_used(demo, tmp_path, capsys):
    from grnbench import autoprep, kb
    assert autoprep.DEFAULT_N_HVG == autoprep._N_HVG[kb.DECISION_INDEX["hvg"].default] == 1000
    paths, _, _ = demo
    out = tmp_path / "run"
    assert prep_main([str(paths["counts"]), "--out", str(out), "--force"]) == 0
    n = ad.read_h5ad(out / "processed.h5ad").n_vars
    note = (f"no --n-hvg: all {n:,} filtered genes are used; the planner and the interface default to the "
            "TFs plus 1,000 highly variable genes")
    assert note in json.loads((out / "manifest.json").read_text())["warnings"]
    assert f"[NO --n-hvg] {note}" in capsys.readouterr().err
    assert prep_main([str(paths["counts"]), "--out", str(tmp_path / "hvg"), "--force", "--n-hvg", "100"]) == 0
    assert not any("no --n-hvg" in w for w in json.loads((tmp_path / "hvg" / "manifest.json").read_text())["warnings"])
    assert "no --n-hvg" not in capsys.readouterr().err
    direct = run_prep(paths["counts"], tmp_path / "direct", force=True)      # the planner's path: no note
    assert not any("no --n-hvg" in w for w in direct.manifest["warnings"])


def test_n_hvg_all_is_an_explicit_choice_the_args_tell_from_no_n_hvg(demo, tmp_path, capsys):
    """LD's review of #162: the note must follow from manifest['args'], so an explicit every-gene choice (a plan's
    n_hvg = None, `--n-hvg all`) and an omitted --n-hvg are told apart there, not only by their warnings."""
    paths, _, _ = demo
    runs = {}
    for name, extra in (("omitted", []), ("all", ["--n-hvg", "all"]), ("upper", ["--n-hvg", " ALL "])):
        assert prep_main([str(paths["counts"]), "--out", str(tmp_path / name), "--force", *extra]) == 0
        runs[name] = json.loads((tmp_path / name / "manifest.json").read_text(encoding="utf-8"))
    printed = capsys.readouterr().err
    assert printed.count("[NO --n-hvg]") == 1
    direct = run_prep(paths["counts"], tmp_path / "direct", force=True).manifest
    api_all = run_prep(paths["counts"], tmp_path / "api_all", force=True, n_hvg="all").manifest
    omitted = runs.pop("omitted")
    assert omitted["args"]["note_all_genes"] is True and omitted["args"]["n_hvg"] is None
    assert any(w.startswith("no --n-hvg:") for w in omitted["warnings"])
    for man in (*runs.values(), direct, api_all):
        assert man["args"]["note_all_genes"] is False and man["args"]["n_hvg"] is None
        assert man["warnings"] == [w for w in omitted["warnings"] if not w.startswith("no --n-hvg:")]
        assert man["steps"] == omitted["steps"] and man["effective_genes"] == omitted["effective_genes"]


# -- #149 comment item 6: organism from gene-name casing only ------------------------------------------------ #

def _with_uns(paths, tmp_path, uns: dict, rename=None):
    """The demo with ``uns`` entries added, next to the same README and TF list (the same text sources)."""
    import shutil
    adata = ad.read_h5ad(paths["counts"])
    if rename is not None:
        adata.var_names = [rename(i, g) for i, g in enumerate(adata.var_names)]
    adata.uns.update(uns)
    folder = tmp_path / paths["counts"].parent.name
    folder.mkdir()
    for f in paths["counts"].parent.iterdir():
        if f.name != "counts.h5ad":
            shutil.copy(f, folder / f.name)
    adata.write_h5ad(folder / "counts.h5ad")
    return diagnose(folder / "counts.h5ad")


def test_organism_from_casing_only_names_a_disagreeing_uns_organism(demo, tmp_path):
    paths, tfs, diag = demo
    other = _with_uns(paths, tmp_path, {"organism": "Danio rerio", "organism_ontology_term_id": "NCBITaxon:7955"})
    for tf_list in (tfs, None):
        base = plan_preprocessing(diag, tf_list=tf_list, methods=[])
        plan = plan_preprocessing(other, tf_list=tf_list, methods=[])
        assert plan.profile["organism"] == base.profile["organism"] == "human"     # the organism is not changed
        assert plan.choices == base.choices and plan.value("tf_source") == base.value("tf_source")
        assert plan.warnings == base.warnings
        reason = plan.decisions["tf_source"].reason
        assert reason.startswith(base.decisions["tf_source"].reason)
        assert ("the organism 'human' comes only from gene-name casing, but uns['organism'] says 'Danio rerio': "
                "the plan does not follow it") in reason, reason
        assert "gene-name casing" not in base.decisions["tf_source"].reason


def test_no_detected_organism_names_the_uns_organism(demo, tmp_path):
    paths, tfs, _ = demo
    diag = _with_uns(paths, tmp_path, {"organism": "Macaca mulatta"}, rename=lambda i, g: f"g{i}")
    assert diag.context.get("organism") is None
    plan = plan_preprocessing(diag, tf_list=None, methods=[])
    assert "no organism was detected, but uns['organism'] says 'Macaca mulatta': the plan does not follow it" \
        in plan.decisions["tf_source"].reason


def test_an_agreeing_uns_organism_adds_no_note(demo, tmp_path):
    paths, tfs, _ = demo
    plan = plan_preprocessing(_with_uns(paths, tmp_path, {"organism": "Homo sapiens"}), tf_list=tfs, methods=[])
    assert "uns['organism']" not in plan.decisions["tf_source"].reason


# -- #116: the knee finds no cliff on an unfiltered matrix --------------------------------------------------- #

def _no_cliff(path):
    """30,000 barcodes, 120 with 150 counts, the rest with 1-40: unfiltered, but the knee finds no cliff."""
    rng = np.random.default_rng(0)
    n, n_genes = 30_000, 40
    totals = np.concatenate([np.full(120, 150), rng.integers(1, 41, n - 120)])
    X = np.stack([rng.multinomial(t, np.full(n_genes, 1 / n_genes)) for t in totals]).astype(np.float32)
    ad.AnnData(X=sp.csr_matrix(X), obs=pd.DataFrame(index=[f"bc{i}" for i in range(n)]),
               var=pd.DataFrame(index=[f"G{i}" for i in range(n_genes)])).write_h5ad(path)


def test_a_knee_without_a_cliff_on_an_unfiltered_matrix_is_recorded_and_warned(tmp_path):
    src = tmp_path / "raw.h5ad"
    _no_cliff(src)
    diag = diagnose(src)
    assert diag.raw_droplets and diag.cell_calling is None and diag.adata.n_obs == 30_000
    assert diag.knee_no_cliff is True
    plan = plan_preprocessing(diag, methods=[])
    assert plan.choices["cell_calling"] == "knee"
    reason = plan.decisions["cell_calling"].reason
    assert "?" not in reason and "nan" not in reason, reason
    assert ("the barcode-rank curve shows no clear cell/empty-droplet cliff, so the knee called no cells; "
            "running this plan stops at cell calling") in reason
    warned = [w for w in plan.warnings if "no clear cell/empty-droplet cliff" in w]
    assert warned == ["The matrix looks unfiltered (30,000 barcodes, most with < 100 counts) but the barcode-rank "
                      "knee found no clear cell/empty-droplet cliff: every barcode was kept, so the profile, the QC "
                      "and the cell bounds describe droplets, not cells."], plan.warnings


def test_knee_chosen_on_a_filtered_matrix_says_why_diagnose_called_no_cells(demo):
    _, tfs, diag = demo
    assert diag.knee_no_cliff is False
    plan = plan_preprocessing(diag, tf_list=tfs, methods=[], choices={"cell_calling": "knee"})
    reason = plan.decisions["cell_calling"].reason
    assert "?" not in reason and "nan" not in reason, reason
    assert "the matrix does not look unfiltered, so no cells were called at diagnosis" in reason
    assert not any("cliff" in w for w in plan.warnings)
