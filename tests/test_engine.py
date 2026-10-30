"""Cleaning engine: knowledge base, profile, selector, plan, execution."""
import copy
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


def test_heuristic_evidence_level_does_not_claim_validation():
    """#112: the level said "validated on the grnbench bench"; no confirmatory experiment has tested a heuristic."""
    text = kb.EVIDENCE_LEVELS["heuristic"]
    assert "not yet confirmed by a confirmatory experiment" in text, text
    if kb.__doc__ is not None:                       # python -OO strips docstrings
        assert text in " ".join(kb.__doc__.split()), "module docstring and EVIDENCE_LEVELS say the same thing"


def test_knee_option_describes_the_threshold_the_code_uses():
    """#112: hygiene.knee_cell_threshold keeps barcodes at or above the inflection after the knee."""
    import inspect
    from grnbench.hygiene import knee_cell_threshold
    floor = inspect.signature(knee_cell_threshold).parameters["min_counts"].default
    summary = kb.DECISION_INDEX["cell_calling"].option("knee").summary
    assert "inflection" in summary and "above the knee" not in summary, summary
    assert str(floor) in summary, summary


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


def test_pseudobulk_reason_names_ignored_compatibility_choices_without_new_warning():
    strategy, reason, warnings = resolve_strategy({
        "sparsity": "pseudobulk", "batch": "combat",
        "normalization": "pearson_residuals",
    })
    assert strategy == "pseudobulk"
    assert "batch option 'combat' is not applied" in reason
    assert "normalized with cpm_log1p, not 'pearson_residuals'" in reason
    assert warnings == []


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


def _without_context_fields(diag, *, obs_columns=(), fact_keys=()):
    """A diagnosis variant for planner fallbacks, without re-running ingestion."""
    out = copy.copy(diag)
    out.adata = diag.adata.copy()
    out.adata.obs = out.adata.obs.drop(columns=list(obs_columns))
    out.context = copy.deepcopy(diag.context)
    out.context.facts = [f for f in out.context.facts if f.key not in set(fact_keys)]
    return out


def test_plan_reads_the_readme_and_explains_every_decision(demo):
    paths, tfs, diag = demo
    plan = plan_preprocessing(diag, tf_list=tfs)
    assert "readme:README.md" in diag.context.sources_read
    assert plan.choices["doublets"] == "scrublet"            # 10x from the README
    assert all(d.reason for d in plan.decisions.values())
    again = PrepPlan.from_dict(json.loads(json.dumps(plan.to_dict(), default=str)))
    assert again.choices == plan.choices and again.contexts == plan.contexts


def test_unknown_decision_key_is_refused_instead_of_ignored(demo):
    _, tfs, diag = demo
    with pytest.raises(KeyError, match="unknown preprocessing decision"):
        plan_preprocessing(diag, tf_list=tfs, choices={"not_a_decision": "none"})


def test_per_batch_qc_without_batch_explains_global_mad_fallback(demo):
    _, tfs, diag = demo
    no_batch = _without_context_fields(
        diag, obs_columns=("batch",), fact_keys=("batch_key",))
    plan = plan_preprocessing(
        no_batch, tf_list=tfs, choices={"cell_qc_thresholds": "mad_per_batch"})
    assert plan.value("per_batch_bounds") is None
    reasons = [plan.decisions[k].reason for k in (
        "min_genes_per_cell", "max_genes_per_cell", "max_counts_per_cell",
        "max_pct_mito_cell")]
    assert all("no batch column was detected" in reason for reason in reasons)
    assert not any("no batch column was detected" in warning for warning in plan.warnings)


def test_hvg_reason_explains_when_small_matrix_keeps_all_genes(demo):
    _, tfs, diag = demo
    plan = plan_preprocessing(diag, tf_list=tfs, choices={"hvg": "tf_plus_1000"})
    assert plan.value("n_hvg") is None
    reason = plan.decisions["n_hvg"].reason
    assert "900 genes" in reason and "all genes are kept instead" in reason


def test_batch_correction_is_not_planned_inside_its_own_cell_type_context(demo, tmp_path, ran):
    """#164: subsetting on the batch column leaves one batch, so centring used to fail every context."""
    paths, tfs, diag = demo
    same_key = copy.copy(diag)
    same_key.adata = diag.adata.copy()
    same_key.context = copy.deepcopy(diag.context)
    same_key.context.facts = [f for f in same_key.context.facts if f.key != "batch_key"]
    same_key.context.add("batch_key", "cell_type", "test", "cell type is the reported batch", 99)

    plan = plan_preprocessing(
        same_key, tf_list=tfs, choices={"context_scope": "per_cell_type"})

    assert len(plan.contexts) >= 2
    assert plan.choices["batch"] == "center_per_batch"  # requested choice remains provenance
    assert plan.value("strategy") == "cpm_log1p"
    assert "also the cell-type context key" in plan.decisions["strategy"].reason
    assert any("every context contains only one batch" in warning for warning in plan.warnings)

    records = execute_plan(plan, paths["counts"], tmp_path / "run", tf_list=tfs)
    assert all(call["strategy"] == "cpm_log1p" for call in ran)
    assert all(record["error"] is None for record in records)


def test_context_fallback_without_cell_type_is_explicit(demo):
    _, tfs, diag = demo
    no_cell_type = _without_context_fields(
        diag, obs_columns=("cell_type",), fact_keys=("cell_type_key",))
    plan = plan_preprocessing(
        no_cell_type, tf_list=tfs, choices={"context_scope": "per_cell_type"})
    assert plan.contexts == [{"name": "all_cells", "subset_key": None, "subset_value": None}]
    assert "no cell-type column was detected" in plan.decisions["context_scope"].reason
    assert "pooled all-cells context" in plan.decisions["context_scope"].reason


def test_context_fallback_without_an_eligible_cell_type_is_explicit(demo):
    _, tfs, diag = demo
    too_small = copy.copy(diag)
    too_small.adata = diag.adata.copy()
    too_small.adata.obs["cell_type"] = [f"type_{i}" for i in range(too_small.adata.n_obs)]
    too_small.context = copy.deepcopy(diag.context)
    plan = plan_preprocessing(
        too_small, tf_list=tfs, choices={"context_scope": "per_cell_type"})
    assert plan.contexts == [{"name": "all_cells", "subset_key": None, "subset_value": None}]
    assert "no cell type has at least 200 cells" in plan.decisions["context_scope"].reason


def test_selected_type_names_the_automatic_largest_type(demo):
    _, tfs, diag = demo
    largest = str(diag.adata.obs["cell_type"].dropna().astype(str).value_counts().index[0])
    plan = plan_preprocessing(
        diag, tf_list=tfs, choices={"context_scope": "selected_type"})
    assert plan.contexts == [{
        "name": f"cell_type={largest}", "subset_key": "cell_type", "subset_value": largest}]
    assert f"selected the largest eligible type '{largest}'" in plan.decisions["context_scope"].reason


def test_requested_subsample_mode_is_recorded_as_the_actual_runtime_mode(demo):
    paths, tfs, diag = demo
    random_plan = plan_preprocessing(
        diag, tf_list=tfs, choices={"subsample": "random"})
    assert random_plan.value("subsample_mode") == "random"
    argv = random_plan.cli_args(paths["counts"], "out", paths["tfs"], seed=42)
    assert argv[argv.index("--subsample-mode") + 1] == "random"

    no_cell_type = _without_context_fields(
        diag, obs_columns=("cell_type",), fact_keys=("cell_type_key",))
    fallback = plan_preprocessing(
        no_cell_type, tf_list=tfs, choices={"subsample": "stratified"})
    assert fallback.choices["subsample"] == "stratified"
    assert fallback.value("subsample_mode") == "random"
    assert "uniform random subsampling runs instead" in fallback.decisions["subsample_mode"].reason


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


# -- user overrides keep what the engine chose (#113) ------------------------------------------ #

_TRIO = ("auto_value", "auto_reason", "auto_source")


def _trio(d):
    return d.auto_value, d.auto_reason, d.auto_source


def _written(plan) -> str:
    return json.dumps(plan.to_dict(), indent=2, default=str)


def test_a_plan_nobody_edited_is_written_as_before(demo):
    paths, tfs, diag = demo
    plan = plan_preprocessing(diag, tf_list=tfs)
    data = plan.to_dict()
    assert all(list(d) == ["key", "label", "value", "reason", "source", "options"]
               for d in data["decisions"].values())
    assert data["decisions"]["min_cells_per_gene"] == {
        "key": "min_cells_per_gene", "label": "Min cells expressing a gene", "value": 3,
        "reason": "genes seen in fewer than max(3, 0.1% of 1200 cells) = 3 cells", "source": "data",
        "options": None}
    for same in ({"doublets": plan.choices["doublets"]}, dict(plan.choices)):   # the selector's own choices
        assert _written(plan_preprocessing(diag, tf_list=tfs, choices=same)) == _written(plan)


def test_a_user_value_keeps_the_data_derived_value_and_its_reason(demo):
    paths, tfs, diag = demo
    plan = plan_preprocessing(diag, tf_list=tfs)
    d = plan.decisions["max_pct_mito_cell"]
    auto = (d.value, d.reason, d.source)
    assert "MADs" in d.reason and d.source == "data"
    plan.set("max_pct_mito_cell", 12.5, "set by the user in the interface")
    assert (d.value, d.reason, d.source) == (12.5, "set by the user in the interface", "user")
    assert _trio(d) == auto
    plan.set("max_pct_mito_cell", 12.5, "another reason")                    # the same value again: no change
    assert d.reason == "set by the user in the interface" and _trio(d) == auto
    plan.set("max_pct_mito_cell", 30.0, "set by the user in the interface")
    assert d.value == 30.0 and _trio(d) == auto                               # the first engine trio stays
    out = plan.to_dict()["decisions"]
    assert list(out["max_pct_mito_cell"]) == ["key", "label", "value", "reason", "source", "options", *_TRIO]
    assert [out["max_pct_mito_cell"][k] for k in _TRIO] == list(auto)
    assert not any(k in v for k in _TRIO for key, v in out.items() if key != "max_pct_mito_cell")


def test_the_engine_value_set_again_is_no_override(demo):
    paths, tfs, diag = demo
    plan = plan_preprocessing(diag, tf_list=tfs)
    before = _written(plan)
    d = plan.decisions["min_genes_per_cell"]
    engine = (d.value, d.reason, d.source)
    plan.set("min_genes_per_cell", engine[0], "set by the user in the interface")
    assert _written(plan) == before
    plan.set("min_genes_per_cell", engine[0] + 100, "set by the user in the interface")
    plan.set("min_genes_per_cell", engine[0], "set by the user in the interface")     # back to the engine's
    assert _written(plan) == before                                                    # no trace, by design
    plan.set("min_genes_per_cell", engine[0] + 50, "set by the user in the interface")
    assert d.value == engine[0] + 50 and _trio(d) == engine                            # captured again
    plan.set("group_key", "batch", "set by the user in the interface")                 # the engine had it off
    g = plan.decisions["group_key"]
    assert _trio(g) == (None, "only used by 'pseudobulk'", "default")
    assert plan.to_dict()["decisions"]["group_key"]["auto_value"] is None
    plan.set("group_key", None, "set by the user in the interface")
    assert g.value is None and _trio(g) == (None, None, None) and g.source == "default"


def test_plan_files_load_with_and_without_the_engine_trio(demo, tmp_path):
    paths, tfs, diag = demo
    old = {"decisions": {"n_hvg": {"key": "n_hvg", "label": "Highly variable genes", "value": 1000,
                                   "reason": "[hvg = tf_plus_1000] default", "source": "knowledge base",
                                   "options": [500, 1000]}}}
    d = PrepPlan.from_dict(old).decisions["n_hvg"]
    assert d.value == 1000 and _trio(d) == (None, None, None)
    assert PrepPlan.from_dict(old).to_dict()["decisions"] == old["decisions"]
    plan = plan_preprocessing(diag, tf_list=tfs)
    plan.set("max_pct_mito_cell", 12.5, "set by the user in the interface")
    plan.set("group_key", "batch", "set by the user in the interface")
    path = plan.save(tmp_path / "plan.json")
    again = PrepPlan.load(path)
    for key in ("max_pct_mito_cell", "group_key"):
        a, b = again.decisions[key], plan.decisions[key]
        assert (a.value, a.source, *_trio(a)) == (b.value, b.source, *_trio(b))
    assert _written(again) == path.read_text(encoding="utf-8")


def test_a_broken_engine_trio_or_an_unknown_key_is_refused_by_name():
    base = {"key": "n_hvg", "label": "HVG", "value": 500, "reason": "set by the user", "source": "user",
            "options": None}
    bad = [({**base, "auto_value": 1000}, "n_hvg"),                                   # no reason
           ({**base, "auto_source": "knowledge base"}, "n_hvg"),
           ({**base, "auto_value": 1000, "auto_reason": "engine"}, "n_hvg"),          # no source
           ({**base, "source": "data", "auto_value": 1000, "auto_reason": "engine",
             "auto_source": "data"}, "n_hvg"),                                         # trio, source not user
           ({**base, "auto_vlaue": 1000}, "auto_vlaue")]
    for entry, name in bad:
        with pytest.raises(ValueError, match=name):
            PrepPlan.from_dict({"decisions": {"n_hvg": entry}})
    plan = PrepPlan.from_dict({"decisions": {"n_hvg": base}})
    plan.decisions["n_hvg"].auto_value = 1000                     # half a trio, made after loading
    with pytest.raises(ValueError, match="n_hvg"):
        plan.to_dict()


def test_a_knowledge_base_override_keeps_the_engine_value_and_reason(demo):
    paths, tfs, diag = demo
    auto = plan_preprocessing(diag, tf_list=tfs)
    assert auto.choices["doublets"] == "scrublet" and auto.choices["cell_qc_thresholds"] == "mad"
    plan = plan_preprocessing(diag, tf_list=tfs, choices={"doublets": "none", "cell_qc_thresholds": "fixed"})
    m, e = plan.decisions["min_genes_per_cell"], auto.decisions["min_genes_per_cell"]
    assert m.value == 200 and m.reason.startswith("[cell_qc_thresholds = fixed]") and m.source == "user"
    assert _trio(m) == (e.value, e.reason, "data") and "MAD" in e.reason      # the MAD number, not "mad"
    d = plan.decisions["drop_doublets"]
    assert d.value is False and _trio(d) == (True, auto.decisions["drop_doublets"].reason, "knowledge base")
    for k, x in plan.decisions.items():         # the others are exactly as the engine wrote them
        if x.auto_reason is None:
            b = auto.decisions[k]
            assert (x.value, x.reason, x.source) == (b.value, b.reason, b.source), k
        else:
            assert x.source == "user", k
    plan.set("min_genes_per_cell", 300, "set by the user in the interface")
    assert m.value == 300 and _trio(m) == (e.value, e.reason, "data")         # a later edit keeps the trio
    plan.set("min_genes_per_cell", e.value, "set by the user in the interface")
    assert (m.value, m.reason, m.source, *_trio(m)) == (e.value, e.reason, "data", None, None, None)
    plan.set("drop_doublets", True)                                             # revert on the KB path too
    assert (d.value, d.source, d.auto_reason) == (True, "knowledge base", None)


def test_strategy_cell_calling_and_global_mad_bounds_keep_the_engine_values(demo):
    import copy
    paths, tfs, diag = demo
    sel = copy.deepcopy(select(Profile(**plan_preprocessing(diag, tf_list=tfs).profile)))
    sel.results["cell_qc_thresholds"].chosen = "mad_per_batch"       # an engine that bounds per batch
    engine = plan_preprocessing(diag, tf_list=tfs, selection=sel)
    plan = plan_preprocessing(diag, tf_list=tfs, selection=sel, choices={
        "cell_qc_thresholds": "mad", "normalization": "pearson_residuals", "cell_calling": "knee"})
    g, e = plan.decisions["min_genes_per_cell"], engine.decisions["min_genes_per_cell"]
    assert "MAD" in g.reason and g.source == "user" and _trio(g) == (e.value, e.reason, "knowledge base")
    s, es = plan.decisions["strategy"], engine.decisions["strategy"]
    assert s.value == "pearson_residuals" and _trio(s) == (es.value, es.reason, "knowledge base")
    for key in ("cell_calling", "max_genes_per_cell"):                # not in the engine's own plan
        assert key not in engine.decisions and plan.decisions[key].source == "user"
        assert _trio(plan.decisions[key]) == (None, "not in the engine's own plan", "knowledge base")
    assert "per_batch_bounds" in engine.decisions and "per_batch_bounds" not in plan.decisions


def test_every_decision_a_choice_changed_names_the_option_behind_it(demo):
    from grnbench.autoprep import KB_DECISIONS_OF
    paths, tfs, diag = demo
    many = {"doublets": "none", "cell_qc_thresholds": "fixed", "gene_filter": "min_cells_drop_qc_genes",
            "normalization": "pearson_residuals", "hvg": "tf_plus_500", "subsample": "stratified",
            "regulators": "all_genes", "context_scope": "pooled", "cell_calling": "knee"}
    other = {"sparsity": "metacells", "cell_qc_thresholds": "mad_per_batch"}
    assert {k for ids in KB_DECISIONS_OF.values() for k in ids} <= set(kb.DECISION_INDEX)
    for choices in (many, other):
        plan = plan_preprocessing(diag, tf_list=tfs, choices=choices)
        changed = [k for k, d in plan.decisions.items() if d.auto_reason is not None]
        assert len(changed) >= 3
        for k in changed:
            assert plan.chosen_options(k) and set(plan.chosen_options(k)) <= set(choices.values()), k
        assert all(not plan.chosen_options(k) for k, d in plan.decisions.items() if d.auto_reason is None)
def test_reproduction_command_always_names_the_seed():
    plan = PrepPlan()
    assert plan.cli_args("in.h5ad", "out", seed=0) == ["in.h5ad", "--out", "out", "--n-hvg", "all", "--seed", "0"]
    with pytest.raises(TypeError):
        plan.cli_args("in.h5ad", "out")           # no command without the seed it replays at


def test_gene_mapping_reaches_plan_kwargs_and_reproduction_command(tmp_path):
    plan = PrepPlan()
    mapping = tmp_path / "gene map.tsv"
    assert plan.run_prep_kwargs(gene_mapping_path=mapping)["gene_mapping_path"] == mapping
    argv = plan.cli_args("in.h5ad", "out", gene_mapping_path=mapping, seed=0)
    assert argv[argv.index("--gene-mapping") + 1] == str(mapping)


def test_one_default_seed_for_the_plan_the_prep_command_and_run_prep():
    import inspect

    from grnbench.autoprep import DEFAULT_SEED
    from grnbench.cli import _build_parser
    from grnbench.pipeline import run_prep
    assert _build_parser().parse_args([]).seed == DEFAULT_SEED == 42
    assert inspect.signature(run_prep).parameters["seed"].default == DEFAULT_SEED
    assert inspect.signature(execute_plan).parameters["seed"].default == DEFAULT_SEED


def test_reproduction_command_names_force_when_the_run_used_it():
    plan = PrepPlan()
    assert "--force" in plan.cli_args("in.h5ad", "out", seed=0, force=True)
    assert "--force" not in plan.cli_args("in.h5ad", "out", seed=0)


@pytest.fixture
def ran(monkeypatch):
    """``run_prep`` replaced by a stub that keeps the keyword arguments of every call."""
    from types import SimpleNamespace

    import grnbench.pipeline
    calls = []

    def fake_run_prep(input_path, out_dir, **kwargs):
        calls.append(kwargs)
        return SimpleNamespace(processed_path=Path(out_dir) / "processed.h5ad",
                               manifest_path=Path(out_dir) / "manifest.json")

    monkeypatch.setattr(grnbench.pipeline, "run_prep", fake_run_prep)
    return calls


@pytest.mark.parametrize("kw, seed", [({"seed": 0}, 0), ({"seed": 42}, 42), ({}, 42)])
def test_contexts_json_records_the_seed_each_context_ran_with(demo, tmp_path, ran, kw, seed):
    """Screening preps at seed 0; the recorded command replayed it at the CLI default 42 (audit 2026-09-27)."""
    paths, tfs, diag = demo
    plan = plan_preprocessing(diag, tf_list=tfs, choices={"context_scope": "per_cell_type"})
    execute_plan(plan, paths["counts"], tmp_path / "run", tf_list=tfs, **kw)
    recorded = json.loads((tmp_path / "run" / "contexts.json").read_text(encoding="utf-8"))
    assert len(recorded) == len(plan.contexts) >= 2 and [c["seed"] for c in ran] == [seed] * len(recorded)
    for rec in recorded:
        assert rec["seed"] == seed and rec["cli"][rec["cli"].index("--seed") + 1] == str(seed), rec


@pytest.mark.parametrize("plan_force, run_force, forced", [(False, True, True), (False, False, False),
                                                           (True, False, True)])
def test_contexts_json_records_force_when_the_run_used_it(demo, tmp_path, ran, plan_force, run_force, forced):
    """execute_plan(force=True) ran past the QC gates but recorded a command without --force (issue #117)."""
    paths, tfs, diag = demo
    plan = plan_preprocessing(diag, tf_list=tfs, choices={"context_scope": "pooled"})
    plan.set("force", plan_force)
    rec = execute_plan(plan, paths["counts"], tmp_path / "run", tf_list=tfs, force=run_force)[0]
    assert [c["force"] for c in ran] == [forced] and ("--force" in rec["cli"]) is forced, rec["cli"]


def _replay_gives_the_same_matrix(rec, out) -> None:
    """Run the recorded command into ``out``: it must succeed and write the matrix the run wrote."""
    from grnbench.cli import main
    assert rec["error"] is None and rec["cli"][:2] == ["grnbench", "prep"]
    assert main([str(out) if a == rec["out_dir"] else a for a in rec["cli"][2:]]) == 0
    a, b = ad.read_h5ad(rec["processed"]), ad.read_h5ad(Path(out) / Path(rec["processed"]).name)
    assert list(a.obs_names) == list(b.obs_names) and list(a.var_names) == list(b.var_names)
    dense = [m.toarray() if sparse.issparse(m) else np.asarray(m) for m in (a.X, b.X)]
    np.testing.assert_array_equal(*dense)


def test_replaying_the_recorded_command_rebuilds_the_same_matrix(demo, tmp_path):
    paths, tfs, diag = demo
    plan = plan_preprocessing(diag, tf_list=tfs, choices={"context_scope": "pooled"})
    plan.set("subsample_cells", 800)          # a random step, so a replay at another seed keeps other cells
    rec = execute_plan(plan, paths["counts"], tmp_path / "run", tf_list=tfs, seed=0)[0]
    _replay_gives_the_same_matrix(rec, tmp_path / "replay")


def test_replaying_context_label_that_starts_with_dash(demo, tmp_path):
    """A recorded value such as ``-alpha`` must not be parsed as another CLI option (#172)."""
    paths, tfs, _ = demo
    adata = ad.read_h5ad(paths["counts"])
    adata.obs["cell_type"] = adata.obs["cell_type"].astype(str).replace({"alpha": "-alpha"})
    counts = tmp_path / "leading-dash.h5ad"
    adata.write_h5ad(counts)

    plan = plan_preprocessing(diagnose(counts), tf_list=tfs,
                              choices={"context_scope": "per_cell_type"})
    records = execute_plan(plan, counts, tmp_path / "run", tf_list=tfs)
    rec = next(r for r in records if r["context"] == "cell_type=-alpha")
    _replay_gives_the_same_matrix(rec, tmp_path / "replay")


def test_replaying_a_forced_run_passes_the_gate_it_passed(demo, tmp_path):
    paths, tfs, diag = demo
    plan = plan_preprocessing(diag, tf_list=tfs, choices={"context_scope": "pooled"})
    plan.add("mito_warn", "Dataset mito gate: warn above median %", 5.0, "test: below the demo's median")
    plan.add("mito_fail", "Dataset mito gate: fail above median %", 10.0, "test: below the demo's median")
    gated = execute_plan(plan, paths["counts"], tmp_path / "gated", tf_list=tfs)[0]
    assert gated["error"].startswith("PrepError: QC gate failure"), gated["error"]
    rec = execute_plan(plan, paths["counts"], tmp_path / "run", tf_list=tfs, force=True)[0]
    _replay_gives_the_same_matrix(rec, tmp_path / "replay")


def _report_warnings(run_dir) -> list[str]:
    text = (Path(run_dir) / "report.md").read_text(encoding="utf-8")
    return text.split("## Warnings\n\n", 1)[1].split("\n\n", 1)[0].splitlines() if "## Warnings" in text else []


def test_replaying_a_plan_that_keeps_every_gene_records_the_same_steps_and_warnings(demo, tmp_path, capsys):
    """LD's review of #162: execute_plan ran a plan's own n_hvg = None without the 'no --n-hvg' note, but the
    recorded command had no --n-hvg, so its replay added the note: the same matrix and manifest['args'], other
    warnings in manifest.json, report.md and the terminal."""
    paths, tfs, diag = demo
    plan = plan_preprocessing(diag, tf_list=tfs, choices={"context_scope": "pooled"})
    assert plan.value("n_hvg") is None                     # 900 genes: every gene is kept
    rec = execute_plan(plan, paths["counts"], tmp_path / "run", tf_list=tfs)[0]
    _replay_gives_the_same_matrix(rec, tmp_path / "replay")
    printed = capsys.readouterr().err
    ran = json.loads(Path(rec["manifest"]).read_text(encoding="utf-8"))
    replayed = json.loads((tmp_path / "replay" / "manifest.json").read_text(encoding="utf-8"))
    assert replayed["warnings"] == ran["warnings"]
    assert _report_warnings(tmp_path / "replay") == _report_warnings(rec["out_dir"])
    assert "--n-hvg" not in printed
    assert replayed["steps"] == ran["steps"]
    assert replayed["args"] == ran["args"] and replayed["feature_selection"] == ran["feature_selection"] is None
    assert rec["cli"][rec["cli"].index("--n-hvg") + 1] == "all"


# -- screening -------------------------------------------------------------------------------- #

def test_reference_split_is_disjoint_complete_and_per_tf_balanced():
    from grnbench.screen import split_reference
    pairs = {(f"TF{i}", f"G{j}") for i in range(5) for j in range(9)}
    a, b = split_reference(pairs, seed=1)
    assert a.isdisjoint(b) and a | b == pairs
    for i in range(5):
        na = sum(s == f"TF{i}" for s, _ in a)
        assert abs(na - sum(s == f"TF{i}" for s, _ in b)) <= 1


@pytest.mark.parametrize("kw, seed", [({}, 0), ({"seed": 3}, 3)])
def test_screening_cleans_every_candidate_at_its_own_seed(demo, tmp_path, monkeypatch, kw, seed):
    import grnbench.autoprep
    from grnbench.screen import screen_candidates
    seen = []

    def fake_execute_plan(plan, input_path, out_dir, *, tf_list=None, seed=None, force=False):
        seen.append(seed)
        return [{"context": c["name"], "out_dir": str(out_dir), "processed": None, "error": "stub", "seed": seed,
                 "cli": []} for c in plan.contexts]

    monkeypatch.setattr(grnbench.autoprep, "execute_plan", fake_execute_plan)
    paths, tfs, diag = demo
    screen_candidates(diag, paths["counts"], tmp_path / "scr", tf_list=tfs, max_candidates=1, **kw)
    assert len(seen) >= 2 and set(seen) == {seed}          # recommended and baseline at least


def test_screening_ranks_candidates_on_the_selection_half_only(demo, tmp_path):
    from grnbench.screen import screen_candidates
    paths, tfs, diag = demo
    gold = pd.read_csv(paths["gold"], sep="\t")
    res = screen_candidates(diag, paths["counts"], tmp_path / "scr", tf_list=tfs,
                            reference_pairs=set(zip(gold.source, gold.target)), max_candidates=1,
                            extra=[{"name": "pooled", "choices": {"context_scope": "pooled"}}])
    assert set(res["rows"]["candidate"]) == {"recommended", "pooled", "baseline"}
    assert res["best"] in {"recommended", "pooled", "baseline"} and "AUPRC ratio" in res["rule"]
    sel, ev = res["reference"]["selection"], res["reference"]["evaluation"]
    assert sel.isdisjoint(ev)
    assert (tmp_path / "scr" / "reference_split.json").exists()
    assert (tmp_path / "scr" / "screening_summary.json").exists()
    diagnostics = json.loads((tmp_path / "scr" / "diagnostics.json").read_text())
    assert diagnostics["selected_candidate"] == res["best"]
    assert set(diagnostics["screen"]) == {"recommended", "pooled", "baseline"}
    for candidate, contexts in diagnostics["prep"].items():
        for entry in contexts.values():
            assert entry["choices"]["strategy"] == res["plans"][candidate].value("strategy")
    assert not [c for c in res["rows"].columns if c.startswith("auto_")]      # screening rows unchanged (#113)
    recommended = json.loads((tmp_path / "scr" / "recommended" / "plan.json").read_text())
    assert not any(k in d for k in _TRIO for d in recommended["decisions"].values())
    pooled = json.loads((tmp_path / "scr" / "pooled" / "plan.json").read_text())["decisions"]["context_scope"]
    assert (pooled["value"], pooled["source"], pooled["auto_value"]) == ("pooled", "user", "per_cell_type")


def _rank_rows(**cands):
    """cands: name -> (edge_stability, error or None)."""
    return pd.DataFrame([{"candidate": n, "context": "all", "error": e, "edge_stability": st,
                          "depth_leakage": 0.0, "batch_leakage": 0.0, "reference_auprc_ratio": r}
                         for n, (st, e, r) in cands.items()])


def test_rank_without_reference_puts_baseline_first_even_if_less_stable():
    from grnbench.screen import rank_candidates
    rows = _rank_rows(baseline=(0.2, None, np.nan), recommended=(0.9, None, np.nan))
    ranking, rule = rank_candidates(rows)
    assert list(ranking.index) == ["baseline", "recommended"] and "D-021" in rule


def test_rank_without_reference_falls_back_when_baseline_failed():
    from grnbench.screen import rank_candidates
    rows = _rank_rows(baseline=(0.9, "refused", np.nan), a=(0.3, None, np.nan), b=(0.5, None, np.nan))
    ranking, _ = rank_candidates(rows)
    assert list(ranking.index) == ["b", "a"]


def test_rank_with_reference_ignores_the_baseline_rule():
    from grnbench.screen import rank_candidates
    rows = _rank_rows(baseline=(0.9, None, 1.0), a=(0.1, None, 2.0))
    ranking, rule = rank_candidates(rows)
    assert ranking.index[0] == "a" and "AUPRC" in rule


def _ctx_rows(**cands):
    """cands: name -> {context: auprc_ratio}; one row per scored context."""
    return pd.DataFrame([{"candidate": n, "context": c, "error": None, "edge_stability": 0.5,
                          "depth_leakage": 0.0, "batch_leakage": 0.0, "reference_auprc_ratio": r}
                         for n, ctxs in cands.items() for c, r in ctxs.items()])


def test_rank_with_reference_never_picks_a_candidate_scored_on_fewer_contexts():
    from grnbench.screen import rank_candidates
    rows = _ctx_rows(recommended={"A": 2.0, "B": 1.0, "C": 1.0}, selected_type={"A": 3.0},
                     other={"A": 1.0, "B": 1.0, "C": 1.0})
    ranking, rule = rank_candidates(rows)
    assert list(ranking.index) == ["recommended", "other", "selected_type"]
    assert ranking["fewer_contexts"].to_dict() == {"recommended": False, "other": False,
                                                   "selected_type": True} and "D-023" in rule


def test_rank_fewer_contexts_flag_ignores_a_failed_superset_so_the_subset_is_not_flagged():
    from grnbench.screen import rank_candidates
    rows = pd.concat([_ctx_rows(subset={"A": 5.0}),
                      _ctx_rows(failed_superset={"A": 1.0, "B": 1.0, "C": 1.0}).assign(
                          error=[None, None, "boom"])], ignore_index=True)
    ranking, _ = rank_candidates(rows)
    assert list(ranking.index) == ["subset", "failed_superset"]
    assert not ranking.loc["subset", "fewer_contexts"] and ranking.loc["subset", "score"] > 0
    assert ranking.loc["failed_superset", "score"] == -np.inf


def test_rank_non_finite_score_is_unrankable_before_the_context_tiers():
    from grnbench.screen import rank_candidates
    ranking, _ = rank_candidates(_ctx_rows(nan_superset={"A": np.nan, "B": np.nan}, finite_subset={"A": 2.0}))
    assert list(ranking.index) == ["finite_subset", "nan_superset"]
    assert not ranking.loc["finite_subset", "fewer_contexts"] and np.isfinite(ranking["score"].iloc[0])
    partial, _ = rank_candidates(_ctx_rows(half_nan={"A": 9.0, "B": np.nan}, finite_subset={"A": 2.0}))
    assert list(partial.index) == ["finite_subset", "half_nan"]


def test_rank_with_reference_keeps_failed_last_and_orders_by_score_on_equal_or_disjoint_contexts():
    from grnbench.screen import rank_candidates
    same = rank_candidates(_ctx_rows(a={"A": 1.0, "B": 1.0}, b={"A": 2.0, "B": 2.0}))[0]
    assert list(same.index) == ["b", "a"]
    overlap = rank_candidates(_ctx_rows(a={"A": 3.0, "B": 3.0}, b={"B": 1.0, "C": 1.0}))[0]
    assert list(overlap.index) == ["a", "b"]
    rows = pd.concat([_ctx_rows(a={"A": 1.0, "B": 1.0}, s={"A": 5.0}),
                      _ctx_rows(f={"A": 9.0, "B": 9.0}).assign(error=[None, "boom"])], ignore_index=True)
    assert list(rank_candidates(rows)[0].index) == ["a", "s", "f"]


def test_screen_candidates_adds_one_baseline_unless_caller_passes_it(demo, tmp_path, monkeypatch):
    import grnbench.autoprep as autoprep
    from grnbench.screen import BASELINE_CHOICES, screen_candidates
    paths, tfs, diag = demo

    def stub(plan, input_path, out_dir, **kw):
        return [{"context": plan.contexts[0]["name"], "processed": None, "error": "stub",
                 "out_dir": str(out_dir), "cli": []}]

    monkeypatch.setattr(autoprep, "execute_plan", stub)
    res = screen_candidates(diag, paths["counts"], tmp_path / "a", tf_list=tfs, max_candidates=1)
    assert [n for n in res["plans"] if n == "baseline"] == ["baseline"]
    assert all(res["plans"]["baseline"].choices[k] == v for k, v in BASELINE_CHOICES.items())
    res = screen_candidates(diag, paths["counts"], tmp_path / "b", tf_list=tfs, max_candidates=1,
                            extra=[{"name": "baseline", "choices": {"hvg": "tf_plus_500"}}])
    assert list(res["plans"]).count("baseline") == 1
    assert res["plans"]["baseline"].choices["hvg"] == "tf_plus_500"  # the caller's own baseline wins


def test_screening_preserves_raw_droplet_detection(demo, tmp_path, monkeypatch):
    """Regression for #30: screening must not re-profile called cells as filtered input."""
    from copy import deepcopy
    from grnbench import autoprep
    from grnbench.screen import screen_candidates

    paths, tfs, diag = demo
    raw = deepcopy(diag)
    raw.raw_droplets = True
    raw.cell_calling = {"threshold": 100.0, "n_called": 1200, "n_barcodes": 30_000}

    def stop_after_planning(plan, input_path, out_dir, **kwargs):
        return [{"context": plan.contexts[0]["name"], "processed": None,
                 "error": "test stop", "out_dir": str(out_dir), "cli": []}]

    monkeypatch.setattr(autoprep, "execute_plan", stop_after_planning)
    result = screen_candidates(raw, paths["counts"], tmp_path / "scr", tf_list=tfs,
                               max_candidates=1)

    assert result["plans"]["recommended"].value("cell_calling") == "knee"


def test_screening_fails_one_candidate_not_all_when_a_context_has_nan(demo, tmp_path, monkeypatch):
    import grnbench.autoprep as autoprep
    from grnbench.screen import screen_candidates
    paths, tfs, diag = demo
    real = autoprep.execute_plan

    def poison_pooled(plan, *a, **kw):
        recs = real(plan, *a, **kw)
        if plan.value("context_scope") == "pooled":
            for rec in recs:
                if rec["processed"]:
                    adata = ad.read_h5ad(rec["processed"])
                    X = adata.X.toarray() if sparse.issparse(adata.X) else np.asarray(adata.X)
                    X[:, 0] = np.nan
                    adata.X = X
                    adata.write_h5ad(rec["processed"])
        return recs

    monkeypatch.setattr(autoprep, "execute_plan", poison_pooled)
    res = screen_candidates(diag, paths["counts"], tmp_path / "scr", tf_list=tfs, max_candidates=1,
                            extra=[{"name": "pooled", "choices": {"context_scope": "pooled"}}])
    rows = res["rows"]
    bad = rows[rows["candidate"] == "pooled"]
    assert bad["error"].notna().all() and bad["error"].str.contains("non-finite").all()
    assert rows[rows["candidate"] == "recommended"]["error"].isna().all()
    assert res["best"] == "recommended"


def test_context_folder_names_are_valid_on_every_os():
    from grnbench.autoprep import safe_folder_name
    for label in ["cell_type=EC:c1", "cell_type=CD4+ T|naive", r'a/b\c?*<>"', "cell_type=" + "A" * 69 + ".1"]:
        name = safe_folder_name(label)
        assert name and not set(name) & set(r':/\|?*<>"') and not name.endswith("."), name
    # the 80-character cut falls just after a '.': Windows would drop it and make 'X.' the folder 'X'
    assert safe_folder_name("cell_type=" + "A" * 69 + ".1") == "cell_type-" + "A" * 69


def test_context_folders_stay_distinct_when_a_label_already_spells_another_labels_fallback(tmp_path):
    """LD's review of #87: the hashed fallback was never checked, so a real label spelling it shared a folder.
    A verifier then found 'X.' and 'X' (one folder on Windows) after the 80-character cut."""
    import hashlib
    import random
    import unicodedata
    from grnbench.autoprep import unique_folder_names

    def same_folder(n):            # how Windows and macOS compare names: case, trailing dots, Unicode normalization
        return unicodedata.normalize("NFD", n).casefold().rstrip(". ")
    sha = hashlib.sha1("cell_type=CD4+T".encode("utf-8")).hexdigest()
    ld = ["cell_type-CD4_T-462ea0", "cell_type=CD4 T", "cell_type=CD4+T"]
    assert sha[:6] == "462ea0" and unique_folder_names(ld)[:2] == ["cell_type-CD4_T-462ea0", "cell_type-CD4_T"]
    cases = [ld,
             [f"cell_type-CD4_T-{sha[:12]}", f"cell_type-CD4_T-{sha}"] + ld,  # the longer fallbacks are labels too
             ["x", "x", "x"],                                                  # one label twice: one hash
             ["T", "t"],                                                       # one folder on Windows and macOS
             ["a" * 100, "a" * 100 + "b", "a" * 100 + "c"],                    # truncated to the same 80 characters
             ["cell_type=" + "A" * 69 + ".1", "cell_type=" + "A" * 69],        # 'X.' after the cut is 'X' on Windows
             ["cell_type=\ud55c", unicodedata.normalize("NFD", "cell_type=\ud55c")]]  # one folder on macOS
    rng = random.Random(0)
    for _ in range(30):            # labels that are other labels' folder names, in any order
        labels = rng.sample(ld, 3)
        labels += rng.sample(unique_folder_names(labels), 2)
        labels += unique_folder_names(labels)
        rng.shuffle(labels)
        cases.append(labels)
    for i, labels in enumerate(cases):
        names = unique_folder_names(labels)
        assert len({same_folder(n) for n in names}) == len(labels), (labels, names)
        assert names == unique_folder_names(labels), "deterministic"
        assert all(n and len(n) <= 80 and not set(n) & set(r':/\|?*<>"') and not n.endswith(".") for n in names)
        for n in names:            # and on this machine's file system: mkdir refuses a folder that already exists
            (tmp_path / str(i) / n).mkdir(parents=True)
        assert len(list((tmp_path / str(i)).iterdir())) == len(labels), (labels, names)
    # folders of past runs keep their names: no collision -> plain name, first collision -> the same 6-character hash
    assert unique_folder_names(["cell_type=CD4+T", "cell_type=CD4 T", "cell_type=CD4:T"]) == [
        "cell_type-CD4_T", "cell_type-CD4_T-be7bf9", "cell_type-CD4_T-10a5ca"]


def test_cell_types_whose_labels_differ_only_in_punctuation_get_their_own_folder(demo, tmp_path):
    """'CD4+T', 'CD4 T' and 'CD4:T' were one folder name: the later contexts failed with 'choose a new path'."""
    from grnbench.autoprep import unique_folder_names
    labels = ["cell_type=CD4+T", "cell_type=CD4 T", "cell_type=CD4:T"]
    names = unique_folder_names(labels)
    assert len(set(names)) == 3 and names[0] == "cell_type-CD4_T"
    assert names == unique_folder_names(labels), "deterministic"
    assert all(not set(n) & set(r':/\|?*<>"') and len(n) <= 80 for n in names)
    assert unique_folder_names(["a", "b"]) == ["a", "b"], "labels that do not collide keep their plain name"

    paths, tfs, diag = demo
    adata = ad.read_h5ad(paths["counts"])
    rename = dict(zip(sorted(adata.obs["cell_type"].astype(str).unique()), ["CD4+T", "CD4 T", "CD4:T"]))
    adata.obs["cell_type"] = adata.obs["cell_type"].astype(str).map(rename)
    src = tmp_path / "collide.h5ad"
    adata.write_h5ad(src)
    plan = plan_preprocessing(diagnose(src), tf_list=tfs, choices={"context_scope": "per_cell_type"})
    recs = execute_plan(plan, src, tmp_path / "run", tf_list=tfs)
    assert len(recs) >= 2 and all(r["error"] is None for r in recs), [r["error"] for r in recs]
    assert len({r["out_dir"] for r in recs}) == len(recs) and all(Path(r["processed"]).exists() for r in recs)


def test_cells_without_a_cell_type_label_are_not_a_cell_type_called_nan(demo, tmp_path):
    paths, tfs, _ = demo
    adata = ad.read_h5ad(paths["counts"])
    labels = adata.obs["cell_type"].astype(object)
    labels.iloc[:300] = np.nan
    adata.obs["cell_type"] = labels
    src = tmp_path / "nolabel.h5ad"
    adata.write_h5ad(src)
    diag = diagnose(src)
    plan = plan_preprocessing(diag, tf_list=tfs, choices={"context_scope": "per_cell_type"})
    names = [c["name"] for c in plan.contexts]
    assert names and not any("nan" in n.lower() for n in names), names
    assert any("300 cells have no label in 'cell_type'" in w for w in plan.warnings), plan.warnings
    assert plan.profile["n_cell_types"] == adata.obs["cell_type"].nunique(), "NaN is not a cell type"


def test_knee_calls_cells_and_rejects_curves_without_a_cliff():
    import numpy as np
    from grnbench.hygiene import knee_cell_threshold
    rng = np.random.default_rng(0)
    totals = np.concatenate([rng.lognormal(np.log(3000), 0.4, 1000), rng.lognormal(np.log(15), 0.8, 30000)])
    call = knee_cell_threshold(totals)
    assert call is not None and 900 <= call["n_called"] <= 1100
    assert knee_cell_threshold(rng.lognormal(np.log(3000), 0.4, 2000)) is None


def test_knee_threshold_respects_floor_knee_and_count_band():
    """#112 review: floor <= threshold <= knee, n_called counts totals >= threshold, None under 200 at the floor."""
    import inspect
    import numpy as np
    from grnbench.hygiene import knee_cell_threshold
    floor = inspect.signature(knee_cell_threshold).parameters["min_counts"].default
    rng = np.random.default_rng(0)
    totals = np.concatenate([rng.lognormal(np.log(3000), 0.4, 1000), rng.lognormal(np.log(15), 0.8, 30000)])
    call = knee_cell_threshold(totals)
    assert call is not None and floor <= call["threshold"] <= call["knee_counts"]
    assert call["n_called"] == int((totals >= call["threshold"]).sum())
    cells, cliff = np.geomspace(5000, 1000, 150), np.geomspace(900, 120, 50)
    at_200 = np.concatenate([cells, cliff, np.full(5000, floor - 1.0)])   # exactly 200 barcodes reach the floor
    assert knee_cell_threshold(at_200) is not None
    assert knee_cell_threshold(at_200[1:]) is None                       # 199


# Barcode-rank curve of 10x's public heart_10k_v3 raw matrix (rank, total UMIs), sampled at
# log-spaced ranks; Cell Ranger called 7,713 cells. EXP-20261004-LM-2 / #51.
HEART_10K_CURVE = [(1, 139767), (11, 79956), (55, 54654), (228, 40284), (939, 25184), (1907, 15958),
                   (3057, 9462), (4901, 3607), (6206, 1989), (7858, 610), (9950, 257), (12599, 181),
                   (20199, 124), (32386, 95), (51924, 76), (65747, 67), (83250, 46), (105413, 8),
                   (133476, 4), (214002, 2), (343110, 1), (1116796, 1)]


def test_knee_ignores_the_ambient_to_empty_cliff_on_heart_10k():
    """With a 10-UMI floor the steepest drop on this curve is ambient -> empty droplets (~75 -> 8
    UMIs) and 92k barcodes were called; DropletUtils' default floor of 100 UMIs finds the cells."""
    import numpy as np
    from grnbench.hygiene import knee_cell_threshold
    ranks, totals = np.array(HEART_10K_CURVE, dtype=float).T
    t = 10 ** np.interp(np.log10(np.arange(1, ranks[-1] + 1)), np.log10(ranks), np.log10(totals))
    assert knee_cell_threshold(t, min_counts=10)["n_called"] > 50_000       # the old floor
    call = knee_cell_threshold(t)
    assert call is not None and 6_000 <= call["n_called"] <= 9_500


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


def _with_tissue_and_high_mito(diag, tissue):
    from copy import deepcopy
    from grnbench.qc import Gate
    d = deepcopy(diag)
    d.context.add("tissue", tissue, "test", f"tissue: {tissue}", 99)
    d.qc.metrics["median_pct_mito"] = 30.0
    d.qc.gates.append(Gate("mito_fraction", "fail", "median % mitochondrial counts per cell: 30.0%", 30.0))
    return d


@pytest.mark.parametrize("tissue", ["liver", "caudate lobe of liver", "hepatic lobule", "heart",
                                    "heart left ventricle", "left ventricle", "interventricular septum",
                                    "left cardiac atrium", "right atrium auricular region",
                                    "atrioventricular node", "sinoatrial node", "myocardium", "kidney",
                                    "renal medulla", "cortex of kidney", "nephron tubule", "mesonephros",
                                    "endocardium", "intrahepatic bile duct",
                                    "skeletal muscle tissue", "brown adipose tissue"])
def test_high_mito_tissues_are_matched_by_whole_words(demo, tissue):
    paths, tfs, diag = demo
    plan = plan_preprocessing(_with_tissue_and_high_mito(diag, tissue), tf_list=tfs)
    assert plan.value("mito_fail") == 60.0, tissue


@pytest.mark.parametrize("tissue", ["adrenal gland", "suprarenal gland", "pediatric bone marrow",
                                    "psychiatric cohort brain", "lateral ventricle", "subventricular zone",
                                    "choroid plexus of fourth ventricle", "lateral septum", "cardia of stomach",
                                    "extracellular matrix", "hepatic flexure of colon", "interrenal primordium"])
def test_a_tissue_that_only_contains_a_high_mito_fragment_keeps_the_mito_gate(demo, tissue):
    """The fragments 'renal', 'atri', 'ventric', 'septum', 'cardi' used to match inside other words, which
    relaxed the gate to 60% and wrote a false reason ("naturally rich in mitochondrial RNA")."""
    paths, tfs, diag = demo
    plan = plan_preprocessing(_with_tissue_and_high_mito(diag, tissue), tf_list=tfs)
    assert "mito_fail" not in plan.decisions, tissue
    assert any("mito_fraction" in w for w in plan.warnings), tissue


# Spellings of high-mito tissues that the substring rule of main relaxed and the first whole-word rule did not
# (stress test of #161): every snake_case, hyphenated or plural spelling of the CELLxGENE listing's tissue labels.
HIGH_MITO_SPELLINGS = [
    "left-ventricle", "anterior-wall-of-left-ventricle", "apical-region-of-left-ventricle",
    "apex of hearts", "basal zone of hearts", "capsule of livers", "caudate lobe of livers", "endocardiums", "hearts",
    "left lobe of livers", "livers", "outflow tract myocardiums", "proepicardiums", "right lobe of livers",
    "anterior_wall_of_left_ventricle", "apex_of_heart", "apical_region_of_left_ventricle", "atrioventricular_node",
    "basal_zone_of_heart", "capsule_of_liver", "cardiac_atrium", "cardiac_muscle_tissue", "cardiac_ventricle",
    "caudate_lobe_of_liver", "cortex_of_kidney", "craniocervical_muscle", "fast_muscle_tissue", "heart_left_ventricle",
    "heart_primordium", "heart_right_ventricle", "hindlimb_muscle", "inner_medulla_of_kidney", "intercostal_muscle",
    "interventricular_septum", "kidney_blood_vessel", "left_atrium_auricular_region", "left_cardiac_atrium",
    "left_lobe_of_liver", "limb_muscle", "liver_primordium", "muscle_of_abdomen", "muscle_of_pelvic_diaphragm",
    "muscle_of_posterior_compartment_of_hindlimb_stylopod", "muscle_organ", "muscle_tissue", "outer_medulla_of_kidney",
    "rectus_abdominis_muscle", "renal_glomerulus", "renal_medulla", "renal_papilla", "renal_pelvis",
    "right_atrium_auricular_region", "right_cardiac_atrium", "right_lobe_of_liver", "sinoatrial_node",
    "skeletal_muscle_organ,_vertebrate", "skeletal_muscle_tissue", "skeletal_muscle_tissue_of_tibialis_anterior",
    "slow_muscle_tissue", "stomach_smooth_muscle", "tongue_muscle"]


@pytest.mark.parametrize("tissue", HIGH_MITO_SPELLINGS)
def test_snake_case_hyphenated_and_plural_high_mito_tissues_are_matched(tissue):
    from grnbench.autoprep import is_high_mito_tissue
    assert is_high_mito_tissue(tissue), tissue
    assert is_high_mito_tissue(tissue.upper() + " "), tissue


@pytest.mark.parametrize("tissue", ["skeletal_muscle_tissue", "cortex_of_kidney", "heart_left_ventricle", "hearts",
                                    "left-ventricle", "endocardiums"])
def test_the_planner_relaxes_the_gate_for_those_spellings(demo, tissue):
    paths, tfs, diag = demo
    plan = plan_preprocessing(_with_tissue_and_high_mito(diag, tissue), tf_list=tfs)
    assert plan.value("mito_fail") == 60.0, tissue


@pytest.mark.parametrize("tissue", ["hepatic  flexure of colon", "hepatic\tflexure of colon", "hepatic-flexure",
                                    "Hepatic_Flexure_Of_Colon"])
def test_the_hepatic_flexure_exclusion_holds_for_any_spacing(demo, tissue):
    """The flexure of the colon is not liver: two spaces, a tab, a hyphen or an underscore used to skip the
    exclusion (the look-ahead expected one space)."""
    paths, tfs, diag = demo
    plan = plan_preprocessing(_with_tissue_and_high_mito(diag, tissue), tf_list=tfs)
    assert "mito_fail" not in plan.decisions, tissue
    from grnbench.autoprep import is_high_mito_tissue
    assert not is_high_mito_tissue(tissue), tissue


@pytest.mark.parametrize("tissue", ["SkeletalMuscle", "HeartLeftVentricle", "LiverTissue", "KidneyCortex", "heart1",
                                    "liver2", "Kidney2Cortex", "LIVER3"])
def test_camel_case_and_digit_suffixed_high_mito_tissues_are_matched(demo, tissue):
    """Stress test round 2: main's substring rule relaxed these; the words are split at lower/upper and
    letter/digit boundaries before lower-casing."""
    paths, tfs, diag = demo
    plan = plan_preprocessing(_with_tissue_and_high_mito(diag, tissue), tf_list=tfs)
    assert plan.value("mito_fail") == 60.0, tissue
    from grnbench.autoprep import is_high_mito_tissue
    assert is_high_mito_tissue(tissue), tissue


@pytest.mark.parametrize("tissue", ["skeletalmuscle", "heartleftventricle", "kidneycortex", "liverwort", "deliver",
                                    "heartwood", "CA1 field of hippocampus"])
def test_run_together_lower_case_words_are_not_matched(tissue):
    """On purpose: finding a word inside a lower-case run of letters is what relaxed liverwort, deliver and heartwood."""
    from grnbench.autoprep import is_high_mito_tissue
    assert not is_high_mito_tissue(tissue), tissue


def _obs_tissue_diag(paths, tmp_path, label, uns=None):
    from copy import deepcopy

    from grnbench.qc import Gate
    a = ad.read_h5ad(paths["counts"])
    a.obs["tissue"] = pd.Categorical([label] * a.n_obs)
    a.uns.update(uns or {})
    a.write_h5ad(tmp_path / "tissue.h5ad")
    diag = deepcopy(diagnose(tmp_path / "tissue.h5ad"))
    diag.qc.metrics["median_pct_mito"] = 30.0
    diag.qc.gates.append(Gate("mito_fraction", "fail", "median % mitochondrial counts per cell: 30.0%", 30.0))
    return diag


def test_a_camel_case_obs_tissue_reaches_the_matcher_through_the_context(demo, tmp_path):
    """Stress test round 2: obs['tissue'] = 'SkeletalMuscle' reached the planner as 'skeletalmuscle'. Round 3: the
    context still stores and shows the label in lower case, as main does; the planner gives the matcher the label
    as written (``DatasetContext.spelling``), and only the matcher splits its words."""
    from grnbench.context import extract_from_text
    paths, tfs, _ = demo
    diag = _obs_tissue_diag(paths, tmp_path, "SkeletalMuscle")
    assert diag.context.get("tissue") == "skeletalmuscle"
    assert [f.value for f in diag.context.evidence("tissue")] == ["skeletalmuscle", "skeletalmuscle"]
    plan = plan_preprocessing(diag, tf_list=tfs)
    assert plan.value("mito_fail") == 60.0
    assert "tissue 'skeletalmuscle' is naturally rich in mitochondrial RNA" in plan.decisions["mito_warn"].reason
    readme = extract_from_text("tissue: KidneyCortex\n", "readme:x")
    assert readme.get("tissue") == "kidneycortex" and readme.spelling("tissue", "kidneycortex") == "KidneyCortex"
    assert extract_from_text("tissue: heart left ventricle\n", "readme:x").get("tissue") == "heart left ventricle"


def test_an_earlier_uns_spelling_does_not_hide_the_curated_obs_spelling(demo, tmp_path):
    """Review of #161 (LD): uns is read before obs, and the first spelling used to be kept whatever its weight, so
    'tissue: skeletalmuscle' in uns (weight 2) hid obs['tissue'] = 'SkeletalMuscle' (weight 3) and kept the 25%
    gate. The spelling kept is now the one with the strongest evidence."""
    paths, tfs, _ = demo
    diag = _obs_tissue_diag(paths, tmp_path, "SkeletalMuscle", uns={"description": "tissue: skeletalmuscle"})
    ctx = diag.context
    assert [(f.source, f.weight) for f in ctx.evidence("tissue")] == [
        ("metadata:uns", 2), ("metadata:obs['tissue']", 2), ("metadata:obs['tissue']", 3)]
    assert ctx.get("tissue") == "skeletalmuscle"
    assert ctx.spelling("tissue", "skeletalmuscle") == "SkeletalMuscle"
    plan = plan_preprocessing(diag, tf_list=tfs)
    assert plan.value("mito_fail") == 60.0
    assert "tissue 'skeletalmuscle' is naturally rich in mitochondrial RNA" in plan.decisions["mito_warn"].reason


def test_the_context_shows_a_tissue_label_as_main_does(demo, tmp_path):
    """Round 3: 'iPSC-derived cardiomyocytes' was stored as 'i psc-derived cardiomyocytes'."""
    paths, _, _ = demo
    diag = _obs_tissue_diag(paths, tmp_path, "iPSC-derived cardiomyocytes")
    assert diag.context.get("tissue") == "ipsc-derived cardiomyocytes"
    assert diag.context.to_dict()["summary"]["tissue"] == "ipsc-derived cardiomyocytes"


@pytest.mark.parametrize("text, words", [("SkeletalMuscle", "skeletal muscle"), ("Liver2Kidney", "liver 2 kidney"),
                                         ("iPSC-derived cardiomyocytes", "ipsc derived cardiomyocytes"),
                                         ("HEK293", "hek 293"), ("heart_Left-Ventricle", "heart left ventricle"),
                                         ("MuscleInvasiveBladderCancer", "muscle invasive bladder cancer")])
def test_the_matcher_splits_camel_case_only_before_a_lower_case_letter(text, words):
    """A lower->Upper boundary followed by a lower-case letter is a word break; 'iPSC' stays one word."""
    from grnbench.autoprep import tissue_words
    assert tissue_words(text) == words


@pytest.mark.parametrize("tissue", ["hepaticflexure", "hepaticflexureofcolon", "HepaticFlexure",
                                    "HepaticFlexureOfColon", "hepatic-flexure", "hepatic_flexure_2"])
def test_the_hepatic_flexure_is_not_liver_in_any_spelling(demo, tissue):
    paths, tfs, diag = demo
    plan = plan_preprocessing(_with_tissue_and_high_mito(diag, tissue), tf_list=tfs)
    assert "mito_fail" not in plan.decisions, tissue
    from grnbench.autoprep import is_high_mito_tissue
    assert not is_high_mito_tissue(tissue)
    assert is_high_mito_tissue("intrahepaticbileduct") and is_high_mito_tissue("hepatocyte")


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
