"""Cleaning engine: knowledge base, profile, selector, plan, execution."""
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from grnbench import kb
from grnbench.autoprep import (PrepPlan, diagnose, execute_plan, mad_bounds,
                               plan_preprocessing, resolve_strategy)
from grnbench.demo import make_demo, write_demo
from grnbench.profile import Profile, cramers_v, technology_from_counts
from grnbench.selector import RULES, select
from grnbench.strategies import cpm_log1p_center_batch, metacells


# -- knowledge base ------------------------------------------------------------------ #

def test_kb_is_consistent_and_every_decision_has_a_rule():
    assert kb.validate_kb() == []
    assert set(RULES) == {d.id for d in kb.DECISIONS}


def test_kb_document_is_generated_from_the_code():
    doc = Path(__file__).resolve().parents[1] / "docs" / "knowledge" / "preprocessing_kb.md"
    assert doc.read_text(encoding="utf-8") == kb.render_markdown(), \
        "regenerate docs/knowledge/preprocessing_kb.md from grnbench.kb.render_markdown()"


def test_embedding_only_integration_is_never_available():
    ok, why = kb.DECISION_INDEX["batch"].option("harmony").availability()
    assert not ok and "not applicable" in why


# -- profile ----------------------------------------------------------------------------- #

def test_count_values_reveal_the_protocol():
    rng = np.random.default_rng(0)
    umi = sparse.csr_matrix(rng.poisson(0.5, (200, 300)).astype(np.float32))
    reads = sparse.csr_matrix((rng.poisson(0.3, (200, 300)) * rng.integers(50, 5000, (200, 300)))
                              .astype(np.float32))
    assert technology_from_counts(umi)[0] == "umi"
    assert technology_from_counts(reads)[0] == "full_length"


def test_cramers_v_separates_confounded_from_balanced_designs():
    batch = ["a"] * 50 + ["b"] * 50
    assert cramers_v(batch, ["x"] * 50 + ["y"] * 50) == pytest.approx(1.0)
    assert cramers_v(batch, ["x", "y"] * 50) == pytest.approx(0.0)


def test_mad_bounds_are_nan_without_spread():
    lo, hi, med = mad_bounds(np.full(50, 300), 3)
    assert np.isnan(lo) and np.isnan(hi) and med == pytest.approx(300)


# -- selector ------------------------------------------------------------------------------ #

def _profile(**kw) -> Profile:
    base = dict(n_cells=3000, n_genes=15000, sparsity=0.9, median_genes_per_cell=1500,
                median_counts_per_cell=4000, library_cv=0.5, median_pct_mito=3.0, n_mito_genes=13,
                frac_ones_among_nonzero=0.5, max_count=500, data_technology_hint="umi",
                raw_droplets_suspected=False, organism="human", technology="umi",
                platform="10x Chromium", time_course=False, id_verdict="symbols")
    base.update(kw)
    return Profile(**base)


def test_full_length_data_get_no_doublet_removal_and_no_residuals():
    sel = select(_profile(technology="full_length", platform="Smart-seq2"))
    assert sel.results["doublets"].chosen == "none"
    assert sel.results["normalization"].chosen == "cpm_log1p"
    assert sel.results["ambient_rna"].chosen == "none"


def test_imputation_is_never_chosen():
    for sparsity in (0.5, 0.99):
        sel = select(_profile(sparsity=sparsity, n_cells=50_000))
        assert sel.results["sparsity"].chosen != "imputation"


def test_confounded_batches_are_not_corrected():
    p = _profile(batch_key="donor", n_batches=3, min_cells_per_batch=800,
                 batch_celltype_cramers_v=0.9, batch_silhouette=0.4)
    assert select(p).results["batch"].chosen == "none"


def test_separated_but_balanced_batches_are_handled():
    p = _profile(batch_key="run", n_batches=2, min_cells_per_batch=1500,
                 batch_celltype_cramers_v=0.05, batch_silhouette=0.35)
    assert select(p).results["batch"].chosen in {"center_per_batch", "combat", "per_batch_consensus"}


def test_several_large_cell_types_mean_one_network_each():
    p = _profile(cell_type_key="cell_type", n_cell_types=4, n_cell_types_with_enough_cells=3,
                 min_cells_per_type=150)
    assert select(p).results["context_scope"].chosen == "per_cell_type"


def test_ensembl_ids_are_flagged_even_while_mapping_is_planned():
    r = select(_profile(id_verdict="ensembl_only")).results["identifier_mapping"]
    assert r.best_overall == "ensembl_to_symbol"


def test_candidates_change_one_uncertain_decision_at_a_time():
    sel = select(_profile(batch_key="run", n_batches=2, min_cells_per_batch=1500,
                          batch_celltype_cramers_v=0.05, batch_silhouette=0.35, library_cv=1.5))
    cands = sel.candidates(max_candidates=4)
    assert cands[0]["name"] == "recommended" and cands[0]["changes"] == {}
    for c in cands[1:]:
        diff = {k for k in c["choices"] if c["choices"][k] != cands[0]["choices"][k]}
        assert diff == set(c["changes"]) and len(diff) == 1


def test_strategy_resolution_priorities():
    assert resolve_strategy({"sparsity": "metacells", "batch": "combat"})[0] == "metacells"
    assert resolve_strategy({"batch": "center_per_batch", "normalization": "pearson_residuals"})[0] \
        == "cpm_log1p_center_batch"
    assert resolve_strategy({"normalization": "pearson_residuals"})[0] == "pearson_residuals"


# -- strategies ------------------------------------------------------------------------------ #

def test_center_per_batch_removes_between_batch_shift_only():
    rng = np.random.default_rng(1)
    X = rng.poisson(5, (200, 30)).astype(np.float32)
    X[100:] *= 3
    a = ad.AnnData(X=sparse.csr_matrix(X), obs=pd.DataFrame({"b": ["x"] * 100 + ["y"] * 100}))
    out = cpm_log1p_center_batch(a, batch_key="b")
    means = pd.DataFrame(out.X).groupby(a.obs["b"].to_numpy()).mean()
    np.testing.assert_allclose(means.loc["x"], means.loc["y"], atol=1e-4)


def test_metacells_sum_integer_counts_within_batch_and_report_purity():
    adata, _, _ = make_demo(n_cells=600, n_genes=400, n_tfs=10, targets_per_tf=8, seed=3)
    out = metacells(adata, cells_per_metacell=20, batch_key="batch", label_key="cell_type", seed=0)
    assert out.obs["n_cells"].sum() == adata.n_obs
    counts = out.layers["counts"].toarray()
    assert np.allclose(counts, np.round(counts)) and counts.sum() == adata.X.sum()
    assert set(out.obs["batch"]) == set(adata.obs["batch"])
    assert out.obs["purity"].between(0, 1).all()


# -- plan and execution --------------------------------------------------------------------- #

@pytest.fixture(scope="module")
def demo(tmp_path_factory):
    d = tmp_path_factory.mktemp("demo")
    paths = write_demo(d, n_cells=1200, n_genes=900, n_tfs=15, targets_per_tf=10)
    tfs = paths["tfs"].read_text().split()
    return paths, tfs, diagnose(paths["counts"])


def test_plan_reads_the_readme_and_explains_every_decision(demo):
    paths, tfs, diag = demo
    plan = plan_preprocessing(diag, tf_list=tfs)
    assert "readme:README.md" in diag.context.sources_read
    assert plan.choices["doublets"] == "scrublet"            # 10x from the README
    assert all(d.reason for d in plan.decisions.values())
    again = PrepPlan.from_dict(json.loads(json.dumps(plan.to_dict(), default=str)))
    assert again.choices == plan.choices and again.contexts == plan.contexts


def test_per_context_execution_calls_doublets_before_subsetting(demo, tmp_path):
    paths, tfs, diag = demo
    plan = plan_preprocessing(diag, tf_list=tfs, choices={"context_scope": "per_cell_type"})
    recs = execute_plan(plan, paths["counts"], tmp_path / "run", tf_list=tfs)
    assert len(recs) >= 2 and all(r["error"] is None for r in recs)
    steps = [s["step"] for s in json.loads(Path(recs[0]["manifest"]).read_text())["steps"]]
    assert steps.index("drop_doublets") < steps.index("subset_context")
    assert "--subset-key" in recs[0]["cli"] and (tmp_path / "run" / "plan.json").exists()


def test_small_contexts_with_metacells_are_warned_before_running(demo):
    paths, tfs, diag = demo
    plan = plan_preprocessing(diag, tf_list=tfs,
                              choices={"context_scope": "per_cell_type", "sparsity": "metacells"})
    assert any("metacells: the smallest context" in w for w in plan.warnings)


def test_per_batch_bounds_and_qc_gene_exclusion_run(demo, tmp_path):
    paths, tfs, diag = demo
    plan = plan_preprocessing(diag, tf_list=tfs, choices={
        "context_scope": "pooled", "cell_qc_thresholds": "mad_per_batch",
        "gene_filter": "min_cells_drop_qc_genes"})
    assert set(plan.value("per_batch_bounds")) == set(diag.adata.obs["batch"].astype(str))
    rec = execute_plan(plan, paths["counts"], tmp_path / "run", tf_list=tfs)[0]
    assert rec["error"] is None
    out = ad.read_h5ad(rec["processed"])
    assert not out.var_names.str.upper().str.startswith("MT-").any()
    assert "--per-batch-bounds" in rec["cli"] and "--exclude-qc-genes" in rec["cli"]


# -- screening -------------------------------------------------------------------------------- #

def test_reference_split_is_disjoint_complete_and_per_tf_balanced():
    from grnbench.screen import split_reference
    pairs = {(f"TF{i}", f"G{j}") for i in range(5) for j in range(9)}
    a, b = split_reference(pairs, seed=1)
    assert a.isdisjoint(b) and a | b == pairs
    for i in range(5):
        na = sum(s == f"TF{i}" for s, _ in a)
        assert abs(na - sum(s == f"TF{i}" for s, _ in b)) <= 1


def test_screening_ranks_candidates_on_the_selection_half_only(demo, tmp_path):
    from grnbench.screen import screen_candidates
    paths, tfs, diag = demo
    gold = pd.read_csv(paths["gold"], sep="\t")
    res = screen_candidates(diag, paths["counts"], tmp_path / "scr", tf_list=tfs,
                            reference_pairs=set(zip(gold.source, gold.target)), max_candidates=1,
                            extra=[{"name": "pooled", "choices": {"context_scope": "pooled"}}])
    assert set(res["rows"]["candidate"]) == {"recommended", "pooled"}
    assert res["best"] in {"recommended", "pooled"} and "AUPRC ratio" in res["rule"]
    sel, ev = res["reference"]["selection"], res["reference"]["evaluation"]
    assert sel.isdisjoint(ev)
    assert (tmp_path / "scr" / "reference_split.json").exists()
    assert (tmp_path / "scr" / "screening_summary.json").exists()


def test_context_folder_names_are_valid_on_every_os():
    from grnbench.autoprep import safe_folder_name
    for label in ["cell_type=EC:c1", "cell_type=CD4+ T|naive", r'a/b\c?*<>"']:
        name = safe_folder_name(label)
        assert name and not set(name) & set(r':/\|?*<>"')


def test_knee_calls_cells_and_rejects_curves_without_a_cliff():
    import numpy as np
    from grnbench.hygiene import knee_cell_threshold
    rng = np.random.default_rng(0)
    totals = np.concatenate([rng.lognormal(np.log(3000), 0.4, 1000), rng.lognormal(np.log(15), 0.8, 30000)])
    call = knee_cell_threshold(totals)
    assert call is not None and 900 <= call["n_called"] <= 1100
    assert knee_cell_threshold(rng.lognormal(np.log(3000), 0.4, 2000)) is None


def test_high_mito_tissue_relaxes_the_dataset_gate_only(demo):
    from copy import deepcopy
    from grnbench.qc import Gate
    paths, tfs, diag = demo
    d = deepcopy(diag)
    d.context.add("tissue", "liver", "test", "tissue: liver", 99)
    d.qc.metrics["median_pct_mito"] = 30.0
    d.qc.gates.append(Gate("mito_fraction", "fail", "median % mitochondrial counts per cell: 30.0%", 30.0))
    plan = plan_preprocessing(d, tf_list=tfs)
    assert plan.value("mito_fail") == 60.0 and plan.value("force") is False
    assert not any("mito_fraction" in w for w in plan.warnings)
    assert plan.value("max_pct_mito_cell") is not None       # per-cell bound still applies
    plain = plan_preprocessing(diag, tf_list=tfs)
    assert "mito_fail" not in plain.decisions


def test_real_bench_records_which_doublet_caller_ran(tmp_path):
    import json
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
    from real_bench import doublet_methods

    recs = []
    for i, steps in enumerate([[{"step": "drop_doublets", "method": "scrublet"}],
                               [{"step": "drop_doublets", "method": "heuristic"}],
                               [{"step": "filter_cells"}]]):
        m = tmp_path / f"m{i}.json"
        m.write_text(json.dumps({"steps": steps}))
        recs.append({"manifest": str(m)})
    assert doublet_methods(recs) == "heuristic,not run,scrublet"
    assert doublet_methods([]) == "unknown"


def test_real_bench_model_matrix_keeps_datasets_with_a_common_prefix_apart():
    import sys
    from pathlib import Path
    import pandas as pd
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
    from real_bench import _model_matrix

    long = "cxg__patch-seq__mus-musculus__primary-motor-cortex__"
    df = pd.DataFrame({"task": [long + "8f98c236", long + "341a0702"],
                       "m_pearson_signed_auto_auprc_ratio": [1.1, 1.2]})
    assert _model_matrix(df, "auto", "auprc_ratio").index.is_unique
