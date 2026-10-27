"""grnbench.seeds: derived child seeds and the three SeedPlan policies (issue #126, step 1a)."""
import inspect
import os
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from grnbench import seeds
from grnbench.seeds import DEFAULT_SEED, LEGACY, MASTER, SEED_MAX, TASK, UNSEEDED, SeedPlan, derive

ROOT = Path(__file__).resolve().parents[1]

COMPONENTS = [
    "prep.sampling", "prep.doublets", "prep.metacells", "profile.nonzero_sample",
    "profile.batch_silhouette", "screen.split_reference", "screen.cell_cap", "screen.batch_leakage",
    "screen.edge_stability", "score.cell_cap", "model.grnboost2", "model.genie3", "model.regdiffusion",
    "model.random_forest", "model.extra_trees", "model.gradient_boosting", "model.tigress",
    "model.aracne3", "eval.null", "eval.stability", "eval.significance", "degrade", "viz", "simulate",
] + [f"extra.component_{i}" for i in range(26)]
PREP = ("prep.sampling", "prep.doublets", "prep.metacells")
SCREEN = ("screen.split_reference", "screen.cell_cap", "screen.batch_leakage", "screen.edge_stability")
SCREENING = ("interface_screening", "real_bench", "engine_bench", "candidate_study")

# models whose cell depends on the runner: a seeded one, one defaulting to 0, a wrapper, one with no seed
MODEL_SAMPLE = ("model.grnboost2", "model.regdiffusion", "model.lasso", "model.grnboost2_motif",
                "model.community", "model.aracne3", "model.pearson_signed")


def _components(entry):
    table = LEGACY[entry]
    names = [c for c in table if c != "model.*"]
    return names + (list(MODEL_SAMPLE) if "model.*" in table else [])


def _task(entry, component):
    table = LEGACY[entry]
    return 3 if table.get(component, table.get("model.*")) is TASK else None


# -- derive ---------------------------------------------------------------------------- #

def test_default_seed_is_42():
    assert DEFAULT_SEED == 42


def test_derive_golden_values():
    # frozen: changing the key construction would move every derived seed of every run
    assert derive(42, "prep.sampling") == 1448362217
    assert derive(7, "model.grnboost2", 3) == 1213207632
    assert derive(0, "viz") == 2045888699


def test_derive_is_deterministic_and_a_python_int():
    a, b = derive(42, "prep.sampling"), derive(42, "prep.sampling")
    assert a == b and type(a) is int


def test_derive_differs_per_component_replicate_and_master():
    assert len(COMPONENTS) == 50 and len(set(COMPONENTS)) == 50
    assert len({derive(42, c) for c in COMPONENTS}) == 50
    reps = [derive(42, "screen.edge_stability", r) for r in range(10)]
    assert len(set(reps)) == 10
    assert derive(42, "screen.edge_stability") not in reps      # no replicate is its own stream
    assert derive(0, "viz") != derive(42, "viz")


def test_component_words_are_distinct():
    names = set(COMPONENTS) | {c for table in LEGACY.values() for c in table}
    names |= {k for table in seeds.KNOBS.values() for k in table.values()}
    assert len({seeds._component_word(n) for n in names}) == len(names)


def test_derive_stays_in_31_bits():
    values = [derive(m, c, r) for m in (0, 1, 42, 2**32 + 5) for c in COMPONENTS[:10]
              for r in (None, 0, 3)]
    assert all(0 <= v < 2**31 for v in values)


def test_derive_does_not_depend_on_python_hash_randomization():
    code = ("from grnbench.seeds import derive; "
            "print(derive(42, 'prep.sampling'), derive(7, 'model.grnboost2', 3))")
    path = os.pathsep.join([str(ROOT / "src"), os.environ.get("PYTHONPATH", "")])
    out = {subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True,
                          env={**os.environ, "PYTHONHASHSEED": h, "PYTHONPATH": path}).stdout.strip()
           for h in ("0", "123", "987654")}
    assert out == {f"{derive(42, 'prep.sampling')} {derive(7, 'model.grnboost2', 3)}"}


@pytest.mark.parametrize("bad", [-1, 1.5, "42", None, True, np.bool_(True)])
def test_derive_refuses_a_bad_master(bad):
    with pytest.raises(ValueError, match="master"):
        derive(bad, "viz")


@pytest.mark.parametrize("bad", ["", None, 3])
def test_derive_refuses_a_bad_component(bad):
    with pytest.raises(ValueError, match="component"):
        derive(42, bad)


@pytest.mark.parametrize("bad", [-1, 0.5, True])
def test_derive_refuses_a_bad_replicate(bad):
    with pytest.raises(ValueError, match="replicate"):
        derive(42, "viz", bad)


# -- the tables ------------------------------------------------------------------------ #

def test_legacy_is_read_only():
    with pytest.raises(TypeError):
        LEGACY["cli_prep"] = {}
    with pytest.raises(TypeError):
        LEGACY["cli_prep"]["prep.sampling"] = 7
    with pytest.raises(TypeError):
        seeds.KNOBS["cli_prep"]["prep.sampling"] = "x"


def test_legacy_cells_are_ints_markers_or_a_none_default():
    for entry, table in LEGACY.items():
        assert table, entry
        for component, cell in table.items():
            if entry == "library":
                assert cell is None or (type(cell) is int and cell >= 0), (component, cell)
            else:
                assert cell in (TASK, MASTER, UNSEEDED) or (type(cell) is int and cell >= 0), \
                    (entry, component, cell)


def test_knobs_follow_the_call_site():
    for entry in ("cli_prep", "cli_run", "interface", "real_bench_noref"):        # run_prep(seed=)
        assert {seeds.knob(c, entry) for c in PREP} == {"prep"}, entry
    for entry in SCREENING:                                       # screen_candidates(seed=) runs the prep too
        assert {seeds.knob(c, entry) for c in PREP + SCREEN} == {"screen"}, entry
    assert len({seeds.knob(c, "library") for c in PREP + SCREEN}) == len(PREP + SCREEN)
    assert seeds.knob("model.community", "cli_infer") == "model.community"    # it seeds its members
    assert seeds.knob("model.lasso", "cli_infer") == "model.lasso"
    assert set(seeds.KNOBS) == set(LEGACY)


@pytest.mark.parametrize("policy", ["master", "derived"])
def test_components_fed_by_one_int_get_one_seed(policy):
    for entry in LEGACY:
        plan = SeedPlan(master=7, policy=policy, entry_point=entry)
        groups = {}
        for component in _components(entry):
            groups.setdefault(seeds.knob(component, entry), []).append(component)
        for name, members in groups.items():
            values = {plan.seed_for(c, task=_task(entry, c)) for c in members
                      if seeds._has_seed(c)}
            assert len(values) <= 1, (entry, name, members, values)


def test_not_covered_is_listed_in_the_docstring():
    assert seeds.NOT_COVERED
    for path, why in seeds.NOT_COVERED:
        assert Path(path).name in seeds.__doc__ and why, path


# -- SeedPlan construction ------------------------------------------------------------- #

def test_plan_fields():
    plan = SeedPlan(entry_point="cli_prep")
    assert (plan.master, plan.policy, plan.entry_point) == (DEFAULT_SEED, "legacy", "cli_prep")


def test_as_dict_records_the_effective_seeds():
    d = SeedPlan(master=7, entry_point="cli_prep").as_dict()
    assert d == {"master": 7, "policy": "legacy", "entry_point": "cli_prep", "task": None,
                 "seeds": {"prep.sampling": 7, "prep.doublets": 7, "prep.metacells": 7}}
    d = SeedPlan(master=7, entry_point="real_bench").as_dict()
    assert d["master"] is None                     # no cell takes the user's seed there
    assert d["seeds"]["screen.split_reference"] == 0 and d["seeds"]["model.grnboost2"] == 42
    assert d["seeds"]["model.regdiffusion"] == 0 and d["seeds"]["model.community"] is None
    d = SeedPlan(entry_point="engine_bench").as_dict(task=3)
    assert d["task"] == 3 and d["seeds"]["simulate"] == 3 and d["seeds"]["score.cell_cap"] == 0
    assert SeedPlan(entry_point="engine_bench").as_dict()["seeds"]["simulate"] is None
    assert SeedPlan(master=np.int64(7), policy="master", entry_point="cli_prep").as_dict()["master"] == 7


def test_master_is_coerced_to_int():
    assert type(SeedPlan(master=np.int64(7), policy="master", entry_point="cli_prep").master) is int
    assert SeedPlan(master=7.0, policy="master", entry_point="cli_prep").master == 7


@pytest.mark.parametrize("bad", [-1, 1.5, True, np.bool_(False), "7", None, SEED_MAX + 1])
def test_bad_master_is_refused(bad):
    with pytest.raises(ValueError, match="master"):
        SeedPlan(master=bad, policy="master", entry_point="cli_prep")


def test_seeds_are_bounded_to_32_bits():
    assert SeedPlan(master=SEED_MAX, policy="master", entry_point="cli_prep").seed_for("prep.sampling") == SEED_MAX
    plan = SeedPlan(entry_point="engine_bench")
    with pytest.raises(ValueError, match="task"):
        plan.seed_for("simulate", task=SEED_MAX + 1)
    with pytest.raises(ValueError, match="master \\+ task"):
        SeedPlan(master=SEED_MAX, policy="master", entry_point="engine_bench").seed_for("simulate", task=1)


@pytest.mark.parametrize("policy", ["legacy", "master", "derived"])
def test_every_policy_needs_an_entry_point(policy):
    with pytest.raises(ValueError, match="needs an entry point"):
        SeedPlan(policy=policy)


def test_construction_validates():
    with pytest.raises(ValueError, match="entry point"):
        SeedPlan(entry_point="no_such_runner")
    with pytest.raises(ValueError, match="policy"):
        SeedPlan(policy="random", entry_point="cli_prep")


def test_plan_is_frozen():
    with pytest.raises(AttributeError):
        SeedPlan(entry_point="cli_prep").master = 1


# -- legacy ------------------------------------------------------------------------------ #

def test_legacy_returns_todays_effective_seeds():
    assert SeedPlan(entry_point="real_bench").seed_for("screen.split_reference") == 0
    assert SeedPlan(entry_point="real_bench_noref").seed_for("prep.doublets") == 42
    assert SeedPlan(entry_point="interface").seed_for("prep.sampling") == 42
    # no seed is passed: the adapter default applies, never arboreto's DEMON_SEED
    assert SeedPlan(entry_point="real_bench").seed_for("model.grnboost2") == 42
    assert SeedPlan(entry_point="real_bench").seed_for("model.regdiffusion") == 0
    assert SeedPlan(entry_point="interface").seed_for("model.grnboost2_motif") == 42
    assert SeedPlan(entry_point="library").seed_for("model.regdiffusion") == 0
    assert SeedPlan(entry_point="real_bench").seed_for("score.cell_cap") == 0
    assert SeedPlan(entry_point="real_bench").seed_for("viz") == 0


def test_cli_cells_follow_the_users_seed():
    assert SeedPlan(master=7, entry_point="cli_prep").seed_for("prep.sampling") == 7
    assert SeedPlan(master=7, entry_point="cli_infer").seed_for("model.grnboost2") == 7
    assert SeedPlan(master=7, entry_point="cli_infer").seed_for("model.lasso") == 7
    assert SeedPlan(master=7, entry_point="cli_eval").seed_for("eval.null") == 7
    assert SeedPlan(master=7, entry_point="cli_eval").seed_for("stability.model") == 42
    plan = SeedPlan(master=7, entry_point="cli_run")
    assert {plan.seed_for(c) for c in ("prep.sampling", "model.genie3", "eval.significance")} == {7}
    assert plan.seed_for("stability.model") == 42
    # the master does not reach cells that do not take the user's seed
    assert SeedPlan(master=7, entry_point="interface").seed_for("prep.sampling") == 42


def test_kwargs_mirror_what_is_passed_today_under_the_real_parameter_name():
    assert SeedPlan(entry_point="real_bench").kwargs("model.grnboost2") == {}
    assert SeedPlan(entry_point="interface").kwargs("model.regdiffusion") == {}
    assert SeedPlan(entry_point="interface").kwargs("prep.sampling") == {"seed": 42}
    assert SeedPlan(master=7, entry_point="cli_infer").kwargs("model.regdiffusion") == {"seed": 7}
    assert SeedPlan(entry_point="engine_bench").kwargs("screen.cell_cap", task=3) == {"seed": 3}
    assert SeedPlan(entry_point="library").kwargs("prep.doublets") == {}
    assert SeedPlan(master=7, entry_point="cli_eval").kwargs("stability.model") == {"algorithm_seed": 42}
    assert SeedPlan(master=7, policy="master", entry_point="cli_eval").kwargs("stability.model") == \
        {"algorithm_seed": 7}
    assert SeedPlan(entry_point="interface").kwargs("viz") == {"seed": 0}      # the ui/logic figure functions
    # aracne3 takes a seed (the program's --seed): the interface passes none, so its own 42 applies
    assert SeedPlan(entry_point="interface").kwargs("model.aracne3") == {}
    assert SeedPlan(entry_point="interface").seed_for("model.aracne3") == 42
    assert SeedPlan(master=7, entry_point="cli_infer").kwargs("model.aracne3") == {"seed": 7}
    with pytest.raises(ValueError, match="pearson_signed"):
        SeedPlan(entry_point="cli_infer").seed_for("model.pearson_signed")    # no seed parameter
    with pytest.raises(ValueError, match="community"):
        SeedPlan(entry_point="interface").seed_for("model.community")     # default None: no seed


@pytest.mark.parametrize("component, entry", [("score.cell_cap", "real_bench"), ("viz", "real_bench"),
                                              ("score.cell_cap", "engine_bench"), ("viz", "library")])
def test_kwargs_refuse_a_component_without_a_parameter(component, entry):
    for policy in ("legacy", "master"):
        with pytest.raises(ValueError, match="no seed parameter"):
            SeedPlan(policy=policy, entry_point=entry).kwargs(component)


def test_kwargs_never_pass_none():
    for entry in LEGACY:
        plan = SeedPlan(entry_point=entry)
        for component in _components(entry):
            if seeds.param(component, entry) is not None:
                assert None not in plan.kwargs(component, task=_task(entry, component)).values(), \
                    (entry, component)


def test_legacy_task_marker_needs_the_task():
    plan = SeedPlan(entry_point="engine_bench")
    assert LEGACY["engine_bench"]["screen.split_reference"] is TASK
    assert plan.seed_for("screen.split_reference", task=3) == 3
    assert plan.seed_for("score.cell_cap") == 0                      # fixed, whatever the task
    with pytest.raises(ValueError, match="task"):
        plan.seed_for("screen.split_reference")
    with pytest.raises(ValueError, match="engine_bench"):
        plan.seed_for("score.cell_cap", task=3)                      # not a task cell
    for bad in (-1, 2.5, True, "3"):
        with pytest.raises(ValueError, match="task"):
            plan.seed_for("simulate", task=bad)
    assert type(plan.seed_for("simulate", task=np.int64(3))) is int


def test_legacy_refuses_replicates():
    with pytest.raises(ValueError, match="fan-out"):
        SeedPlan(entry_point="cli_eval").seed_for("stability.model", replicate=2)


def test_legacy_sweep_reads_the_table_and_never_derives(monkeypatch):
    def boom(*a, **k):
        raise AssertionError("legacy must not derive")
    monkeypatch.setattr(seeds, "derive", boom)
    for entry, table in LEGACY.items():
        plan = SeedPlan(master=7, entry_point=entry)
        for component in _components(entry):
            cell = table.get(component, table.get("model.*"))
            effective = cell
            if entry == "library" or cell is UNSEEDED:
                effective = LEGACY["library"].get(component)
            elif cell is MASTER:
                effective = 7
            elif cell is TASK:
                effective = 3
            if component.startswith("model.") and component not in LEGACY["library"]:
                effective = None                                  # the runner has no seed parameter
            if effective is None:
                with pytest.raises(ValueError):
                    plan.seed_for(component, task=_task(entry, component))
            else:
                assert plan.seed_for(component, task=_task(entry, component)) == effective, (entry, component)


# -- master and derived ---------------------------------------------------------------- #

def test_master_gives_the_master_everywhere():
    plan = SeedPlan(master=7, policy="master", entry_point="real_bench")
    assert {plan.seed_for(c) for c in ("prep.sampling", "model.grnboost2", "viz")} == {7}
    assert plan.kwargs("model.grnboost2") == {"seed": 7}


def test_task_cells_under_master_and_derived():
    master = SeedPlan(master=7, policy="master", entry_point="engine_bench")
    assert [master.seed_for("simulate", task=t) for t in (0, 1, 2)] == [7, 8, 9]
    derived = SeedPlan(master=7, policy="derived", entry_point="engine_bench")
    assert derived.seed_for("simulate", task=2) == derive(7, "simulate", 2)
    for plan in (master, derived):
        with pytest.raises(ValueError, match="task"):
            plan.seed_for("simulate")


def test_derived_seeds_each_knob():
    plan = SeedPlan(master=7, policy="derived", entry_point="cli_run")
    assert plan.seed_for("prep.sampling") == plan.seed_for("prep.doublets") == derive(7, "prep")
    assert plan.seed_for("eval.null") == derive(7, "eval.null")
    assert plan.seed_for("model.grnboost2") == derive(7, "model.grnboost2")
    assert plan.seed_for("eval.null") != plan.seed_for("prep.sampling")
    screening = SeedPlan(master=7, policy="derived", entry_point="real_bench")
    assert screening.seed_for("prep.sampling") == screening.seed_for("screen.cell_cap") == derive(7, "screen")


def test_replicates_are_derived_under_master_and_derived_only():
    for policy in ("master", "derived"):
        plan = SeedPlan(master=7, policy=policy, entry_point="cli_prep")
        assert plan.seed_for("prep.sampling", replicate=2) == derive(7, "prep", 2)
    with pytest.raises(ValueError, match="fan-out"):
        SeedPlan(master=7, entry_point="cli_prep").seed_for("prep.sampling", replicate=2)
    assert "replicate" not in inspect.signature(SeedPlan.kwargs).parameters


# -- names ------------------------------------------------------------------------------- #

@pytest.mark.parametrize("policy", ["legacy", "master", "derived"])
def test_unknown_components_are_refused_in_every_policy(policy):
    plan = SeedPlan(master=7, policy=policy, entry_point="cli_prep")
    for bad in ("no_such.component", "prep.sampleing", "prep.sampling ", "model.no_such_model", ""):
        with pytest.raises(ValueError):
            plan.seed_for(bad)
        with pytest.raises(ValueError):
            plan.kwargs(bad)


@pytest.mark.parametrize("policy", ["legacy", "master", "derived"])
def test_a_component_the_entry_point_does_not_run_is_refused(policy):
    plan = SeedPlan(master=7, policy=policy, entry_point="cli_prep")
    with pytest.raises(ValueError, match="cli_prep"):
        plan.seed_for("model.grnboost2")                          # prep runs no model
    with pytest.raises(ValueError, match="cli_prep"):
        plan.kwargs("screen.cell_cap")


def test_has_seed_reads_the_runner_signature(monkeypatch):
    from grnbench.infer import base
    from grnbench.infer.base import MethodSpec

    def seeded(expr, tf_names=None, seed: int = 5, **_):
        return pd.DataFrame()

    def unseeded(expr, tf_names=None, **_):
        return pd.DataFrame()
    monkeypatch.setitem(base.REGISTRY, "my_seeded", MethodSpec("my_seeded", "x", seeded, origin="python"))
    monkeypatch.setitem(base.REGISTRY, "my_unseeded", MethodSpec("my_unseeded", "x", unseeded, origin="python"))
    plan = SeedPlan(entry_point="interface")
    assert plan.seed_for("model.my_seeded") == 5                  # its own default, from its signature
    assert plan.kwargs("model.my_seeded") == {}
    with pytest.raises(ValueError, match="my_unseeded"):
        plan.seed_for("model.my_unseeded")
    assert SeedPlan(master=7, entry_point="cli_infer").seed_for("model.my_seeded") == 7


# -- kwargs fit the real signatures ------------------------------------------------------ #

def _functions(entry, component):
    """The real callables that take this component's seed at this entry point."""
    from grnbench import autoprep, hygiene, pipeline, plugins, profile, screen, simulate, strategies
    from grnbench.eval import nulls, significance, stability
    from grnbench.infer.base import REGISTRY
    from grnbench.ui import logic
    if component.startswith("model."):
        runner = REGISTRY[component[len("model."):]].runner
        base = inspect.getclosurevars(runner).nonlocals.get("base")     # the motif and cisTarget wrappers
        return [base if base is not None and "seed" not in inspect.signature(runner).parameters else runner]
    if component in PREP + SCREEN and entry in SCREENING:
        return [screen.screen_candidates]
    if component in PREP and entry in ("cli_prep", "cli_run"):
        return [pipeline.run_prep]
    if component in PREP:
        return [autoprep.execute_plan] if entry != "library" else [
            {"prep.sampling": pipeline.run_prep, "prep.doublets": hygiene.detect_doublets,
             "prep.metacells": strategies.metacells}[component]]
    if component == "degrade":
        if str(ROOT / "benchmarks") not in sys.path:
            sys.path.insert(0, str(ROOT / "benchmarks"))
        import candidate_study
        return [candidate_study.degrade]
    return {
        "screen.split_reference": [screen.split_reference], "screen.cell_cap": [screen.screen_processed],
        "screen.batch_leakage": [screen.batch_leakage], "screen.edge_stability": [screen.edge_stability],
        "profile.nonzero_sample": [profile._nonzero_values], "profile.batch_silhouette": [profile.batch_silhouette],
        "eval.null": [nulls.null_metrics], "eval.stability": [stability.subsample_stability],
        "eval.significance": [significance.gold_thinning_sensitivity],
        "stability.model": [stability.subsample_stability],
        "viz": [logic.subsample_cells, logic.embedding_frame, logic.expression_frame, logic.network_figure],
        "simulate": [simulate.Scenario, simulate.random_scenario, simulate.PRESETS[next(iter(simulate.PRESETS))].with_seed],
        "plugin.command": [plugins.command_method_spec({"name": "x", "command": "true"}).runner],
    }[component]


def test_kwargs_bind_to_the_real_signatures():
    checked = 0
    for entry in LEGACY:
        for policy in ("legacy", "master"):
            plan = SeedPlan(master=7, policy=policy, entry_point=entry)
            for component in _components(entry):
                if seeds.param(component, entry) is None:
                    continue
                kwargs = plan.kwargs(component, task=_task(entry, component))
                for fn in _functions(entry, component):
                    params = inspect.signature(fn).parameters
                    inspect.signature(fn).bind_partial(**kwargs)
                    if seeds._has_seed(component):
                        assert set(kwargs) <= set(params), (entry, component, fn, kwargs)
                    checked += 1
    assert checked > 100


# -- through run_method's real call shape ----------------------------------------------- #

def _spy(runner, seen):
    sig = inspect.signature(runner)

    def spy(*args, **kwargs):
        bound = sig.bind(*args, **kwargs)
        bound.apply_defaults()
        seen.append(bound.arguments["seed"])
        return pd.DataFrame({"source": ["g0"], "target": ["g1"], "weight": [1.0]})
    spy.__signature__ = sig
    return spy


SEEDED = ("grnboost2", "genie3", "regdiffusion", "random_forest", "lasso", "tigress", "extra_trees",
          "gradient_boosting")


@pytest.mark.parametrize("entry", ["interface", "real_bench", "engine_bench", "cli_infer", "cli_run"])
def test_kwargs_through_run_method_give_the_registered_seed(entry, monkeypatch):
    """The seed a runner sees, with its real signature and defaults, is ``seed_for``; prov is as today."""
    from grnbench.infer import base
    from grnbench.infer.base import MethodSpec, run_method
    rng = np.random.default_rng(1)
    adata = ad.AnnData(rng.random((12, 4)).astype(np.float32),
                       var=pd.DataFrame(index=[f"g{i}" for i in range(4)]))
    plan = SeedPlan(master=7, entry_point=entry)
    for name in SEEDED:
        seen = []
        monkeypatch.setitem(base.REGISTRY, name, MethodSpec(name, "spy", _spy(base.REGISTRY[name].runner, seen)))
        kwargs = plan.kwargs(f"model.{name}")
        _, prov = run_method(name, adata, assume_normalized=True, **kwargs)
        assert seen == [plan.seed_for(f"model.{name}")], (entry, name)
        assert ("seed" in prov["params"]) == bool(kwargs), (entry, name)
