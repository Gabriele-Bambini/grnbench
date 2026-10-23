"""An option is available only with its dependencies AND a grnbench executor (issue #128).

An option grnbench cannot run is never chosen by the selector, and forcing one fails loudly. Closes
the audit finding "external options silently dropped" (audit/2026-09-27 findings.json:513-525, AUDIT
section 7.5). The A/B tests compare with golden values written from a read-only copy of main 6917c88
(tests/golden/plan_values.py). Known gaps, not fixed here: #131 (regulators recorded only), #132.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse

from grnbench import executors, kb
from grnbench.autoprep import PrepPlan, plan_preprocessing, resolve_strategy
from grnbench.profile import Profile, build_profile
from grnbench.selector import select

sys.path.insert(0, str(Path(__file__).resolve().parent / "golden"))
import plan_values as golden  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
GOLDEN = json.loads((ROOT / "tests" / "golden" / "option_availability_main.json").read_text(encoding="utf-8"))
NO_EXECUTOR = "no grnbench executor for {} yet"


@pytest.fixture()
def registered_env(monkeypatch):
    """The environment of the registered runs: no R, no scvi, no magic, no CellBender."""
    monkeypatch.setattr(kb, "_have", lambda req: req in golden.REGISTERED_ENV)


@pytest.fixture()
def everything_installed(monkeypatch):
    """Every dependency reported as installed: R packages, scvi, magic, CellBender, ..."""
    monkeypatch.setattr(kb, "_have", lambda req: True)


@pytest.fixture()
def fresh_registry(monkeypatch):
    """Registrations made by a test stay in that test."""
    monkeypatch.setattr(executors, "_REGISTRY", {})
    monkeypatch.setattr(executors, "PLAN_STEPS", dict(executors.PLAN_STEPS))


@pytest.fixture(scope="module")
def datasets(tmp_path_factory):
    return golden.datasets(tmp_path_factory.mktemp("presets"))


def _availability(decision, option):
    return kb.DECISION_INDEX[decision].availability(option)


def _profile(**kw) -> Profile:
    base = dict(n_cells=3000, n_genes=15000, sparsity=0.9, median_genes_per_cell=1500,
                median_counts_per_cell=4000, library_cv=0.5, median_pct_mito=3.0, n_mito_genes=13,
                frac_ones_among_nonzero=0.5, max_count=500, data_technology_hint="umi",
                raw_droplets_suspected=False, organism="human", technology="umi",
                platform="10x Chromium", time_course=False, id_verdict="symbols")
    base.update(kw)
    return Profile(**base)


def _mismatch(a, b, path="$"):
    """First place where two JSON values differ (floats to 1e-9 relative: BLAS differs by OS), or None."""
    if isinstance(a, float) or isinstance(b, float):
        ok = isinstance(a, (int, float)) and isinstance(b, (int, float)) and (
            a == b or math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12))
        return None if ok else f"{path}: {a!r} != {b!r}"
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b):
            return f"{path}: keys {sorted(set(a) ^ set(b))}"
        return next((m for k in a if (m := _mismatch(a[k], b[k], f"{path}.{k}"))), None)
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return f"{path}: length {len(a)} != {len(b)}"
        return next((m for i, (x, y) in enumerate(zip(a, b)) if (m := _mismatch(x, y, f"{path}[{i}]"))), None)
    return None if a == b else f"{path}: {a!r} != {b!r}"


# -- availability: reasons, order, repeated ids ------------------------------------------- #

def test_reasons_come_in_order_and_keep_the_historical_text(registered_env, monkeypatch):
    assert _availability("batch", "harmony") == (False, "not applicable to GRN inference (see note)")
    assert _availability("regulators", "collectri") == (False, "planned in grnbench, not implemented yet")
    assert _availability("normalization", "scran") == (False, "missing: r:scran")
    assert _availability("batch", "scvi_denoised") == (False, "missing: py:scvi")
    assert _availability("doublets", "scrublet") == (True, "")
    monkeypatch.setattr(kb, "_have", lambda req: True)
    assert _availability("normalization", "scran") == (False, NO_EXECUTOR.format("normalization=scran"))
    assert _availability("regulators", "collectri") == (False, "planned in grnbench, not implemented yet")
    assert _availability("batch", "harmony") == (False, "not applicable to GRN inference (see note)")
    monkeypatch.setattr(kb, "_have", lambda req: False)
    full = kb.DECISION_INDEX["doublets"].availability_full("scrublet")
    assert (full.ok, full.reason, full.executor) == (False, "missing: py:scanpy", "python")
    assert full.install_hint == "pip install scanpy"
    scran = kb.DECISION_INDEX["normalization"].availability_full("scran")
    assert scran.install_hint == "" and "would not make it available" in scran.detail


def test_forbidden_options_with_their_dependency_say_d002(everything_installed):      # D-002
    for decision, option in executors.FORBIDDEN:
        full = kb.DECISION_INDEX[decision].availability_full(option)
        assert (full.ok, full.reason) == (False, "not offered: D-002, no imputation or smoothing")
        assert "D-002" in full.detail


def test_repeated_option_ids_are_resolved_by_their_decision(everything_installed, fresh_registry):
    nones = [d for d in kb.DECISIONS if any(o.id == "none" for o in d.options)]
    assert len(nones) == 6
    for d in nones:
        full = d.availability_full("none")
        assert full.ok and full.executor == "python", d.id
    for decision in ("hvg", "regulators"):
        assert _availability(decision, "all_genes") == (True, "")
    with pytest.raises(KeyError, match="doublets"):
        _availability("doublets", "all_genes")
    # One requirement (py:scvi), two pairs: an executor for one does not make the other available.
    executors.register_executor("doublets", "solo", lambda *a: None, kind="python")
    assert _availability("doublets", "solo") == (True, "")
    assert "solo" in executors.PLAN_STEPS["doublets"]
    assert _availability("batch", "scvi_denoised") == (False, "not offered: D-002, no imputation or smoothing")


def test_option_availability_still_reports_dependencies_only(everything_installed):
    assert kb.DECISION_INDEX["doublets"].option("scdblfinder").availability() == (True, "")
    assert _availability("doublets", "scdblfinder") == (False, NO_EXECUTOR.format("doublets=scdblfinder"))


def test_option_table_shows_the_decision_availability(everything_installed):
    row = next(r for r in kb.option_table() if (r["decision"], r["option"]) == ("doublets", "scdblfinder"))
    assert not row["available"] and row["unavailable_reason"] == NO_EXECUTOR.format("doublets=scdblfinder")


# -- executors: PLAN_STEPS, registry, D-002, D-001 -------------------------------------------- #

def test_knowledge_base_labels_agree_with_plan_steps():
    for d in kb.DECISIONS:
        for o in d.options:
            assert executors.has_executor(d.id, o.id) == (o.implementation == "grnbench"), (d.id, o.id)
    assert executors.RECORDED_ONLY <= {(d, o) for d, opts in executors.PLAN_STEPS.items() for o in opts}


def test_every_plan_step_that_is_not_recorded_only_changes_what_runs(datasets, registered_env):
    """Against the decision's default: run_prep's arguments, the contexts or the strategy."""
    diag, tfs = datasets["demo"]

    def runs(decision, option):
        plan = plan_preprocessing(diag, tf_list=tfs, choices={decision: option}, methods=[])
        return plan.run_prep_kwargs(tfs), plan.contexts, plan.value("strategy")

    for decision, options in executors.PLAN_STEPS.items():
        default = runs(decision, kb.DECISION_INDEX[decision].default)
        for option in options:
            if option == kb.DECISION_INDEX[decision].default:
                continue
            changed = runs(decision, option) != default
            assert changed != ((decision, option) in executors.RECORDED_ONLY), (decision, option)


def test_register_executor_checks_kind_fn_and_the_knowledge_base(fresh_registry):
    fn = lambda *a: None  # noqa: E731
    with pytest.raises(ValueError, match="kind"):
        executors.register_executor("doublets", "scdblfinder", fn, kind="java")
    with pytest.raises(TypeError):
        executors.register_executor("doublets", "scdblfinder", None, kind="r")
    with pytest.raises(KeyError):
        executors.register_executor("doublets", "no_such_option", fn, kind="r")
    ex = executors.register_executor("doublets", "scdblfinder", fn, kind="r")
    assert executors.executor_for("doublets", "scdblfinder") is ex and ex.kind == "r"
    with pytest.raises(ValueError, match="already"):
        executors.register_executor("doublets", "scdblfinder", fn, kind="r")
    with pytest.raises(ValueError, match="already"):
        executors.register_executor("doublets", "scrublet", fn, kind="python")
    with pytest.raises(ValueError, match="kind"):
        executors.Executor("doublets", "scrublet", "java", fn)
    assert not hasattr(executors, "EXECUTORS")      # the registry is private


def test_executor_for_refuses_an_unknown_pair():
    with pytest.raises(KeyError):
        executors.executor_for("doublets", "no_such_option")
    with pytest.raises(KeyError):
        executors.has_executor("no_such_decision", "none")


def test_smoothing_and_imputation_never_get_an_executor(everything_installed, fresh_registry):  # D-002
    assert executors.FORBIDDEN == {("sparsity", "imputation"), ("batch", "scvi_denoised")}
    for decision, option in executors.FORBIDDEN:
        with pytest.raises(ValueError, match="D-002"):
            executors.register_executor(decision, option, lambda *a: None, kind="python")
        assert not executors.has_executor(decision, option)


@pytest.mark.parametrize("as_sparse", [False, True])
def test_executor_output_must_be_integer_counts(as_sparse):                             # D-001
    wrap = sparse.csr_matrix if as_sparse else np.asarray
    good = np.array([[0, 3, 0], [1, 7, 2]], dtype=np.float32)
    executors.check_integer_counts(wrap(good))
    for bad in (good + 0.5, -good, np.where(good > 2, np.nan, good)):
        with pytest.raises(ValueError, match="D-001"):
            executors.check_integer_counts(wrap(bad))


# -- the audit regression -------------------------------------------------------------------- #

def test_audit_installed_dependencies_no_longer_select_options_grnbench_cannot_run(everything_installed):
    r = select(_profile(raw_droplets_suspected=True)).results
    assert (r["doublets"].chosen, r["doublets"].best_overall) == ("scrublet", "scdblfinder")
    assert (r["cell_calling"].chosen, r["cell_calling"].best_overall) == ("knee", "emptydrops")
    assert r["ambient_rna"].chosen == "none"
    best = next(o for o in r["doublets"].options if o.option == "scdblfinder")
    assert not best.available and best.unavailable_reason == NO_EXECUTOR.format("doublets=scdblfinder")
    profiles = [_profile(raw_droplets_suspected=True),
                _profile(batch_key="run", n_batches=2, min_cells_per_batch=1500,
                         batch_celltype_cramers_v=0.05, batch_silhouette=0.35),
                _profile(batch_key="donor", n_batches=3, min_cells_per_batch=800,
                         batch_celltype_cramers_v=0.9, batch_silhouette=0.4),
                _profile(sparsity=0.99, n_cells=50_000, time_course=True)]
    for p in profiles:
        for res in select(p).results.values():
            for option in (res.chosen, res.runner_up):    # runner-ups become screened candidates
                if option:
                    assert _availability(res.decision, option)[0], (res.decision, option)
    scvi = next(o for o in select(profiles[1]).results["batch"].options if o.option == "scvi_denoised")
    assert not scvi.available


def test_audit_plan_keeps_removing_doublets_and_calling_cells(datasets, everything_installed):
    diag, tfs = datasets["raw_droplets"]
    plan = plan_preprocessing(diag, tf_list=tfs, methods=[])
    assert plan.value("drop_doublets") is True and plan.value("cell_calling") == "knee"
    assert plan.choices["doublets"] == "scrublet" and plan.choices["cell_calling"] == "knee"
    blocked = plan.option_notes["blocked"]
    assert (blocked["doublets"]["option"], blocked["doublets"]["reason"]) == \
        ("scdblfinder", NO_EXECUTOR.format("doublets=scdblfinder"))
    assert blocked["cell_calling"]["option"] == "emptydrops"
    assert any("scDblFinder" in w and NO_EXECUTOR.format("doublets=scdblfinder") in w for w in plan.warnings)


# -- A/B against main 6917c88 ------------------------------------------------------------------ #

def test_option_table_is_that_of_main(registered_env):
    assert _mismatch(golden.option_table(), GOLDEN["option_table"]) is None


@pytest.mark.parametrize("env", ["registered_env", "everything_installed"])
def test_plans_and_candidates_are_those_of_main(datasets, env, request):
    """Without R: the whole plan (warnings, reasons, sources, selection, profile) and the screening
    candidates as on main. With every dependency installed: the same choices, decision values and
    warning count (main picked scDblFinder, EmptyDrops and DecontX there and ran none of them); only the
    unavailability reasons in the selection and the warnings name the missing executor instead."""
    request.getfixturevalue(env)
    got = golden.plan_values(datasets)
    assert set(got) == set(GOLDEN["plans"])
    for case, want in GOLDEN["plans"].items():
        have = got[case]
        if env == "registered_env":
            assert _mismatch(have, want) is None, (case, _mismatch(have, want))
        else:
            for key in ("choices", "contexts"):
                assert have["plan"][key] == want["plan"][key], (case, key)
            values = {k: d["value"] for k, d in have["plan"]["decisions"].items()}
            assert values == {k: d["value"] for k, d in want["plan"]["decisions"].items()}, case
            assert len(have["plan"]["warnings"]) == len(want["plan"]["warnings"]), case
            assert have["candidates"] == want["candidates"], case


def test_candidate_study_variants_are_those_of_main_and_of_the_c7_rows(registered_env):
    variants = golden.variant_candidates()
    assert json.dumps(variants, sort_keys=True) == json.dumps(GOLDEN["variant_candidates"], sort_keys=True)
    rows = sorted((ROOT / "benchmarks" / "results" / "EXP-20261003-LM-1" / "rows").glob("*.json"))
    assert len(rows) == 60
    for f in rows:   # byte-identical, key order included
        assert json.dumps(json.loads(f.read_text(encoding="utf-8"))["skipped_variants"]) == json.dumps(variants[1])


def test_with_every_dependency_installed_scran_is_skipped_for_want_of_an_executor(everything_installed):
    _, skipped = golden.variant_candidates()
    assert skipped == [{"variant": "normalization=scran", "reason": NO_EXECUTOR.format("normalization=scran")},
                       {"variant": "batch=harmony", "reason": "not applicable to GRN inference (see note)"}]


def test_every_option_in_the_committed_rows_can_run_in_the_registered_environment(registered_env):
    seen = set()
    for f in ROOT.glob("benchmarks/results/**/rows/*.json"):
        row = json.loads(f.read_text(encoding="utf-8"))
        rules = row.get("rules_choices")
        sets = [json.loads(rules) if isinstance(rules, str) else rules] if rules else []
        sets += [c.get("choices") or {} for c in row.get("candidates") or [] if isinstance(c, dict)]
        seen |= {(d, o) for choices in sets for d, o in choices.items()}
    assert len(seen) > 20
    from grnbench.screen import BASELINE_CHOICES
    seen |= set(BASELINE_CHOICES.items())                  # the fixed baseline every benchmark screens
    bad = sorted((d, o) for d, o in seen
                 if d not in kb.DECISION_INDEX or o not in executors.PLAN_STEPS[d] or not _availability(d, o)[0])
    assert bad == []


# -- plan side: an option grnbench cannot run fails loudly ---------------------------------- #

@pytest.mark.parametrize("decision,option,reason", [
    ("doublets", "scdblfinder", "missing: r:scDblFinder"),
    ("cell_calling", "emptydrops", "missing: r:DropletUtils"),
    ("normalization", "scran", "missing: r:scran"),
    ("regulators", "collectri", "planned in grnbench, not implemented yet"),
    ("batch", "harmony", "not applicable to GRN inference (see note)"),
    ("sparsity", "imputation", "missing: py:magic")])
def test_a_forced_option_grnbench_cannot_run_fails_loudly(datasets, registered_env, decision, option, reason):
    diag, tfs = datasets["demo"]
    with pytest.raises(ValueError) as err:
        plan_preprocessing(diag, tf_list=tfs, choices={decision: option}, methods=[])
    assert f"{decision}={option} was set by the user" in str(err.value) and reason in str(err.value)


def test_a_forced_option_with_a_missing_python_dependency_names_the_install(datasets, monkeypatch):
    monkeypatch.setattr(kb, "_have", lambda req: req == "py:scanpy")       # no decoupler
    diag, tfs = datasets["demo"]
    with pytest.raises(ValueError, match="missing: py:decoupler. Install: pip install decoupler"):
        plan_preprocessing(diag, tf_list=tfs, choices={"regulators": "dorothea"}, methods=[])


def test_with_every_dependency_installed_forcing_scdblfinder_still_fails(datasets, everything_installed):
    diag, tfs = datasets["demo"]
    with pytest.raises(ValueError, match=NO_EXECUTOR.format("doublets=scdblfinder")):
        plan_preprocessing(diag, tf_list=tfs, choices={"doublets": "scdblfinder"}, methods=[])


def test_a_selected_option_grnbench_cannot_run_is_an_error(datasets, registered_env):
    diag, tfs = datasets["demo"]
    sel = select(build_profile(diag.adata, diag.qc.metrics, diag.context, id_verdict=diag.id_verdict,
                               tf_list=tfs, raw_droplets=diag.raw_droplets))
    sel.results["doublets"].chosen = "scdblfinder"
    with pytest.raises(ValueError, match="doublets=scdblfinder was chosen by the selector"):
        plan_preprocessing(diag, tf_list=tfs, selection=sel, methods=[])


def test_screening_an_extra_candidate_grnbench_cannot_run_fails_before_anything_runs(
        datasets, registered_env, tmp_path):
    from grnbench.demo import write_demo
    from grnbench.screen import screen_candidates

    diag, tfs = datasets["demo"]
    paths = write_demo(tmp_path / "data")
    with pytest.raises(ValueError, match="doublets=scdblfinder was set by the user"):
        screen_candidates(diag, paths["counts"], tmp_path / "screen", tf_list=tfs, max_candidates=1,
                          extra=[{"name": "dbl", "choices": {"doublets": "scdblfinder"}}])
    assert not [p for p in (tmp_path / "screen").glob("*") if p.is_dir()]


@pytest.mark.parametrize("choices", [{"normalization": "scran"}, {"normalization": "sctransform"},
                                     {"batch": "harmony"}, {"batch": "scvi_denoised"},
                                     {"batch": "per_batch_consensus"}, {"sparsity": "imputation"},
                                     {"sparsity": "metacells", "normalization": "scran"}])
def test_resolve_strategy_has_no_silent_fallback(choices):
    with pytest.raises(ValueError, match="no plan step"):
        resolve_strategy(choices)


def test_the_unavailable_best_option_warning_names_what_runs(datasets, registered_env):
    diag, tfs = datasets["demo"]
    default = plan_preprocessing(diag, tf_list=tfs, methods=[])
    assert any(w.endswith("using 'scrublet'.") for w in default.warnings if w.startswith("doublets:"))
    override = plan_preprocessing(diag, tf_list=tfs, choices={"doublets": "none"}, methods=[])
    assert [w for w in override.warnings if w.startswith("doublets:")] == [
        "doublets: best option 'scDblFinder' is unavailable (missing: r:scDblFinder); using 'none'."]
    assert len(override.warnings) == len(default.warnings)


def test_option_notes_carry_detail_and_hint_and_old_plans_load(datasets, registered_env):
    from grnbench.ui.logic import option_notes

    diag, tfs = datasets["demo"]
    plan = plan_preprocessing(diag, tf_list=tfs, methods=[])
    note = plan.option_notes["blocked"]["doublets"]
    assert note["option"] == "scdblfinder" and "would not make it available" in note["detail"]
    lines = option_notes(plan)
    assert any(line.startswith("scDblFinder: grnbench cannot run it yet") for line in lines)
    old = plan.to_dict()
    old.pop("option_notes")
    assert PrepPlan.from_dict(json.loads(json.dumps(old, default=str))).option_notes == {"blocked": {}}
    assert option_notes(PrepPlan()) == []
