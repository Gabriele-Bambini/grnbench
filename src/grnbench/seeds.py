"""Seeds: one master seed, named child seeds, and the seeds each entry point uses today.

Issue #126, step 1a. This module only records (:data:`LEGACY`) and derives (:func:`derive`,
:class:`SeedPlan`); nothing calls it yet, so no seed and no output changes.

``tests/test_seeds_legacy.py`` reads every LEGACY value back from the code: the default, constant, parser
default or call-site literal it comes from; the path that hands it on; and, in every function that holds
it, that each RNG constructor and each known RNG consumer (PCA, KMeans, silhouette_score, scanpy's pca and
scrublet, UMAP, the scikit-learn estimators, spring_layout, RegDiffusionTrainer) is seeded with it. The
model internals, the simulator and the task enumerations are included. It does not see a random call
through an object outside that list, nor arithmetic on a seed inside a callee it does not trace.

A LEGACY cell is an int (that seed is passed), ``MASTER`` (the user's ``--seed``: ``SeedPlan.master``),
``TASK`` (the benchmark task's seed) or ``UNSEEDED`` (no seed is passed, so the callee's own default
applies). The ``library`` row holds those defaults; ``None`` there means the default is None (nothing is
passed on). A model without a seed parameter (aracne3, the correlation methods) has no library row.
``model.*`` is the row for every registered model at that entry point.

:data:`KNOBS` follow the call site: components fed by one ``seed=`` int share a knob, and so one seed under
every policy. :data:`PARAMS` give the parameter name each component's seed goes in.

One run can span two entry points: ``quick_run(screen=True)`` in the interface is
``interface_screening`` (the screening) and then ``interface`` (the cleaning, the models, the figures).

Not covered (their seeds are not in LEGACY, see :data:`NOT_COVERED`): beeline_backtest.py,
motif_sensitivity.py, calibration.py, demo.py, plugins.py (Python and entry-point plugin runners),
explore.py, rbridge.py, epr_power.py.
"""
from __future__ import annotations

import enum
import hashlib
import inspect
import numbers
from dataclasses import dataclass
from types import MappingProxyType
from typing import Final, Literal, Mapping, Union

import numpy as np

DEFAULT_SEED: Final = 42
SEED_MAX: Final = 2**32 - 1          # numpy, scikit-learn and torch all take a seed in [0, 2**32)
POLICIES = ("legacy", "master", "derived")
Policy = Literal["legacy", "master", "derived"]


class Marker(enum.Enum):
    """A LEGACY cell that is not a fixed number."""
    TASK = "task"            # the benchmark task's own seed, task["seed"]
    MASTER = "master"        # the user's --seed, i.e. SeedPlan.master
    UNSEEDED = "unseeded"    # no seed is passed: the callee's own default (the library row) applies

    def __repr__(self) -> str:
        return f"seeds.{self.name}"


TASK: Final = Marker.TASK
MASTER: Final = Marker.MASTER
UNSEEDED: Final = Marker.UNSEEDED
Cell = Union[int, None, Literal[Marker.TASK, Marker.MASTER, Marker.UNSEEDED]]

# Seeds this module does not record, with where they live (paths from the repository root).
NOT_COVERED: tuple[tuple[str, str], ...] = (
    ("benchmarks/beeline_backtest.py", "BEELINE back-test: --seed 42 to the runners and the nulls"),
    ("benchmarks/motif_sensitivity.py", "EXP-20261003-GB-1: its permutation seeds range(n_perm) and "
                                        "bootstrap seed 0 stay in the script"),
    ("src/grnbench/calibration.py", "grnbench calibrate: config seed 42, split by SeedSequence.spawn"),
    ("src/grnbench/demo.py", "grnbench demo: make_demo seed 0"),
    ("src/grnbench/rbridge.py", "R bridge: seed and RNG kind fixed by r/_wrapper.R per call"),
    ("src/grnbench/explore.py", "grnbench.explore (exploratory): SEED = 0 for its subsamples"),
    ("benchmarks/epr_power.py", "EPR power analysis: BOOTSTRAP_SEED = 42 for its bootstrap"),
    ("src/grnbench/plugins.py", "Python and entry-point plugin runners keep their own defaults "
                                "(command plugins are the library row plugin.command)"),
)


def _is_count(x) -> bool:
    return isinstance(x, (int, np.integer)) and not isinstance(x, bool) and x >= 0


def _as_seed(x, what: str) -> int:
    """``x`` as a Python int in ``[0, SEED_MAX]``: integers and integral floats only, never a bool."""
    if isinstance(x, (bool, np.bool_)):
        raise ValueError(f"{what} must be an integer in [0, 2**32 - 1], got the bool {x!r}")
    if isinstance(x, numbers.Integral):
        value = int(x)
    elif isinstance(x, (float, np.floating)) and float(x).is_integer():
        value = int(x)
    else:
        raise ValueError(f"{what} must be an integer in [0, 2**32 - 1], got {x!r}")
    if not 0 <= value <= SEED_MAX:
        raise ValueError(f"{what} must be an integer in [0, 2**32 - 1], got {x!r}")
    return value


def _component_word(component: str) -> int:
    """A 32-bit word from the name, stable across processes (no ``hash()``, so no PYTHONHASHSEED)."""
    return int.from_bytes(hashlib.blake2b(component.encode("utf-8"), digest_size=4).digest(), "big")


def derive(master: int, component: str, replicate: int | None = None) -> int:
    """Child seed of ``master`` for ``component`` (and ``replicate``), in ``[0, 2**31)``.

    Keyed by name, so a new component never moves the others; 31 bits suit scikit-learn,
    torch and R's ``set.seed``.
    """
    if not _is_count(master):
        raise ValueError(f"master seed must be a non-negative integer, got {master!r}")
    if not isinstance(component, str) or not component:
        raise ValueError(f"component must be a non-empty name, got {component!r}")
    if replicate is not None and not _is_count(replicate):
        raise ValueError(f"replicate must be a non-negative integer or None, got {replicate!r}")
    key = (_component_word(component),) + ((int(replicate),) if replicate is not None else ())
    state = np.random.SeedSequence(int(master), spawn_key=key).generate_state(1, np.uint32)[0]
    return int(state) & 0x7FFFFFFF


# The seed each component gets today, per entry point (origin/main 6917c88). Paths are under
# src/grnbench/, except real_bench.py, engine_bench.py and candidate_study.py, which are in benchmarks/.
_LEGACY: dict[str, dict[str, Cell]] = {
    # Each function's own default when called with no seed: src/grnbench, plus viz.SEED,
    # simulate.Scenario.seed and benchmarks/candidate_study.degrade. model.aracne3 has no row:
    # run_aracne3 (infer/aracne3.py:33) has no seed parameter and the jar gets none.
    "library": {
        "prep.sampling": 42,                  # pipeline.py:128, its one rng at :153
        "prep.doublets": 0,                   # hygiene.py:259, scrublet random_state at :288
        "prep.metacells": 42,                 # strategies.py:245, PCA :268, KMeans :281
        "profile.nonzero_sample": 0,          # profile.py:69, rng :73
        "profile.batch_silhouette": 0,        # profile.py:110, rng :120, PCA :128, silhouette :132
        "screen.split_reference": 0,          # screen.py:36, rng :38
        "screen.cell_cap": 0,                 # screen.py:120 screen_processed, cap :135
        "screen.batch_leakage": 0,            # screen.py:80, PCA and silhouette :91-92
        "screen.edge_stability": 0,           # screen.py:101, permutation :109
        "model.genie3": 42,                   # infer/tree_methods.py:82; None would give arboreto's
        "model.grnboost2": 42,                #   DEMON_SEED (:72), not 42; :93
        "model.random_forest": 42,            # infer/classical.py:245
        "model.lasso": 42,                    # infer/classical.py:251 (only with alpha; the default path
                                              #   is deterministic)
        "model.regdiffusion": 0,              # infer/regdiffusion_adapter.py:23
        "model.tigress": 42,                  # infer/regression.py:43
        "model.extra_trees": 42,              # infer/regression.py:294
        "model.gradient_boosting": 42,        # infer/regression.py:303
        "model.community": None,              # infer/community.py:40; None, so members get none (:60)
        "model.grnboost2_motif": 42,          # infer/cistarget.py:330-332 hands **params to grnboost2
        "model.regdiffusion_motif": 0,        # infer/cistarget.py:330-332 -> regdiffusion
        "model.grnboost2_cistarget": 42,      # infer/cistarget.py:255-257 -> grnboost2
        "plugin.command": 42,                 # plugins.py:374, the runner of every command plugin
        "eval.null": 42,                      # eval/nulls.py:82, SeedSequence(seed).spawn at :119
        "eval.stability": 42,                 # eval/stability.py:27, rng :37
        "eval.significance": 42,              # eval/significance.py:41, rng :68
        "stability.model": 42,                # eval/stability.py:28 algorithm_seed, to run_method at :50
        "viz": 0,                             # viz.py:34 SEED
        "simulate": 0,                        # simulate.py:59 Scenario.seed; random_scenario draws from
                                              #   10_000 + seed (simulate.py:361)
        "degrade": 0,                         # benchmarks/candidate_study.py:117, rng :123
    },
    # grnbench prep -> run_prep
    "cli_prep": {
        "prep.sampling": MASTER,              # cli.py:79 --seed (default 42), run_prep(seed=args.seed) :200
        "prep.doublets": MASTER,              # run_prep passes its seed on: pipeline.py:346
        "prep.metacells": MASTER,             # pipeline.py:447
    },
    # grnbench infer: run_method(seed=args.seed) for every runner, so regdiffusion gets 42, not its 0,
    # and aracne3 gets a seed that it drops (**_)
    "cli_infer": {
        "model.*": MASTER,                    # cli_infer.py:53 --seed (default 42), passed at :179
    },
    # grnbench eval
    "cli_eval": {
        "eval.null": MASTER,                  # cli_eval.py:44 --seed (default 42), passed at :374
        "eval.stability": MASTER,             # cli_eval.py:44, passed at :249
        "eval.significance": MASTER,          # cli_eval.py:44, passed at :450
        "stability.model": 42,                # cli_eval.py:65 --algorithm-seed, :250; eval/stability.py:50
                                              #   adds i only with --vary-algorithm-seed
    },
    # grnbench run: --seed goes to prep, infer and eval; --algorithm-seed is never passed
    "cli_run": {
        "prep.sampling": MASTER,              # cli_run.py:79 --seed (default 42), to prep at :112
        "prep.doublets": MASTER,              # cli_run.py:112 -> pipeline.py:346
        "prep.metacells": MASTER,             # cli_run.py:112 -> pipeline.py:447
        "model.*": MASTER,                    # cli_run.py:145 to infer
        "eval.null": MASTER,                  # cli_run.py:184 to eval
        "eval.stability": MASTER,             # cli_run.py:184
        "eval.significance": MASTER,          # cli_run.py:184
        "stability.model": 42,                # not passed, so cli_eval.py:65 applies
    },
    # Streamlit interface: profile, cleaning, models, figures
    "interface": {
        "profile.nonzero_sample": 0,          # profile.py:69, called without a seed at :79
        "profile.batch_silhouette": 0,        # profile.py:110, called without a seed at :180
        "prep.sampling": 42,                  # ui/app.py:699, ui/logic.py:1125 execute_plan without seed
        "prep.doublets": 42,                  #   -> autoprep.py:555 seed=42, passed at :583
        "prep.metacells": 42,                 # autoprep.py:555, :583
        "model.*": UNSEEDED,                  # ui/logic.py:311 run_method without seed
        "viz": 0,                             # ui/logic.py:568, :576, :696, :775 seed=0; ui/app.py:1004,
                                              #   :1093, :1169 pass none
    },
    # Streamlit interface, "Run and compare" and the quick run's screening
    "interface_screening": {
        "screen.split_reference": 0,          # ui/app.py:643, ui/logic.py:1119 screen_candidates without
        "screen.cell_cap": 0,                 #   seed -> screen.py:249 seed=0, passed at :280, :331
        "screen.batch_leakage": 0,            # screen.py:249, :331 -> :141
        "screen.edge_stability": 0,           # screen.py:249, :331 -> :142
        "prep.sampling": 0,                   # screen.py:311 execute_plan(seed=seed), so not 42
        "prep.doublets": 0,                   # screen.py:311
        "prep.metacells": 0,                  # screen.py:311
        "model.*": UNSEEDED,                  # screen.py:113, :152 run_method without seed
    },
    # benchmarks/real_bench.py, datasets with a reference (all registered rows of C4-C6)
    "real_bench": {
        "profile.nonzero_sample": 0,          # profile.py:69, :79
        "profile.batch_silhouette": 0,        # profile.py:110, :180
        "screen.split_reference": 0,          # real_bench.py:347 screen_candidates without seed
                                              #   -> screen.py:249 seed=0, passed at :280
        "screen.cell_cap": 0,                 # screen.py:249, :331
        "screen.batch_leakage": 0,            # screen.py:249, :331
        "screen.edge_stability": 0,           # screen.py:249, :331
        "prep.sampling": 0,                   # screen.py:311
        "prep.doublets": 0,                   # screen.py:311
        "prep.metacells": 0,                  # screen.py:311
        "score.cell_cap": 0,                  # real_bench.py:253 default_rng(0); engine_bench.py:339 via :352
        "model.*": UNSEEDED,                  # real_bench.py:274 run_method without seed; screening
                                              #   screen.py:113, :152; proxies engine_bench.py:348
        "viz": 0,                             # viz.py:34, drawn by real_bench.py:618 build_run
    },
    # benchmarks/real_bench.py, datasets without a reference (the row's "seed": 0 is only a label)
    "real_bench_noref": {
        "profile.nonzero_sample": 0,          # profile.py:69, :79
        "profile.batch_silhouette": 0,        # profile.py:110, :180
        "prep.sampling": 42,                  # real_bench.py:397 execute_plan without seed -> autoprep.py:555
        "prep.doublets": 42,                  # autoprep.py:555, :583
        "prep.metacells": 42,                 # autoprep.py:555, :583
    },
    # benchmarks/engine_bench.py, simulated tasks
    "engine_bench": {
        "simulate": TASK,                     # engine_bench.py:233-236 Scenario from task["seed"];
                                              #   random_scenario adds 10_000 (simulate.py:361)
        "profile.nonzero_sample": 0,          # profile.py:69, :79
        "profile.batch_silhouette": 0,        # profile.py:110, :180
        "screen.split_reference": TASK,       # engine_bench.py:455 screen_candidates(seed=task["seed"])
        "screen.cell_cap": TASK,              # engine_bench.py:455
        "screen.batch_leakage": TASK,         # engine_bench.py:455
        "screen.edge_stability": TASK,        # engine_bench.py:455
        "prep.sampling": TASK,                # engine_bench.py:455 -> screen.py:311
        "prep.doublets": TASK,                # screen.py:311
        "prep.metacells": TASK,               # screen.py:311
        "score.cell_cap": 0,                  # engine_bench.py:339 default_rng(0), whatever the task seed
        "model.*": UNSEEDED,                  # engine_bench.py:348 run_method; BENCH_METHOD_PARAMS (:43)
                                              #   has no seed, passed at :408
    },
    # benchmarks/candidate_study.py (C7); make_tasks gives every task seed 0 (candidate_study.py:180)
    "candidate_study": {
        "degrade": TASK,                      # candidate_study.py:253 -> :151 degrade(..., seed)
        "profile.nonzero_sample": 0,          # profile.py:69, :79
        "profile.batch_silhouette": 0,        # profile.py:110, :180
        "screen.split_reference": TASK,       # candidate_study.py:291-292
                                              #   screen_candidates(seed=task["seed"])
        "screen.cell_cap": TASK,              # candidate_study.py:292
        "screen.batch_leakage": TASK,         # candidate_study.py:292
        "screen.edge_stability": TASK,        # candidate_study.py:292
        "prep.sampling": TASK,                # candidate_study.py:292 -> screen.py:311
        "prep.doublets": TASK,                # screen.py:311
        "prep.metacells": TASK,               # screen.py:311
        "score.cell_cap": 0,                  # engine_bench.py:339 via candidate_study.py:302
        "model.*": UNSEEDED,                  # engine_bench.py:348 via candidate_study.py:302
                                              #   (EVAL_METHODS only, all deterministic)
    },
}
LEGACY: Mapping[str, Mapping[str, Cell]] = MappingProxyType(
    {entry: MappingProxyType(table) for entry, table in _LEGACY.items()})
del _LEGACY

_PREP = ("prep.sampling", "prep.doublets", "prep.metacells")
_SCREEN = ("screen.split_reference", "screen.cell_cap", "screen.batch_leakage", "screen.edge_stability")
# One call takes one int and uses it for several components; a caller cannot split it without changing
# the callee, so those components are one knob. run_prep(seed=) (directly or through execute_plan) feeds
# prep.*; screen_candidates(seed=) feeds screen.* and, through execute_plan(seed=seed) at screen.py:311,
# the prep it runs. community(seed=) hands its seed to every member (community.py:60), so the members it
# runs share model.community's knob by construction. Any other component is its own knob.
_FAN_OUT = {
    "cli_prep": {"prep": _PREP}, "cli_run": {"prep": _PREP}, "interface": {"prep": _PREP},
    "real_bench_noref": {"prep": _PREP},
    "interface_screening": {"screen": _SCREEN + _PREP}, "real_bench": {"screen": _SCREEN + _PREP},
    "engine_bench": {"screen": _SCREEN + _PREP}, "candidate_study": {"screen": _SCREEN + _PREP},
}
KNOBS: Mapping[str, Mapping[str, str]] = MappingProxyType({
    entry: MappingProxyType({c: next((k for k, members in _FAN_OUT.get(entry, {}).items() if c in members), c)
                             for c in table if c != "model.*"})
    for entry, table in LEGACY.items()})


def knob(component: str, entry_point: str) -> str:
    """The knob of ``component`` at ``entry_point``: a model not listed is its own knob."""
    return KNOBS[entry_point].get(component, component)


# The parameter a component's seed goes in; None: no parameter, so kwargs() has nothing to pass.
_PARAM = {
    "stability.model": "algorithm_seed",          # eval/stability.py:28
    "score.cell_cap": None,                       # a literal at the call site: default_rng(0)
    "viz": None,                                  # the module constant viz.SEED
}
_PARAM_AT = {("interface", "viz"): "seed"}        # the ui/logic figure functions take seed= (logic.py:568)
PARAMS: Mapping[str, Mapping[str, str | None]] = MappingProxyType({
    entry: MappingProxyType({c: _PARAM_AT.get((entry, c), _PARAM.get(c, "seed")) for c in table})
    for entry, table in LEGACY.items()})


def param(component: str, entry_point: str) -> str | None:
    """The parameter name of ``component``'s seed at ``entry_point`` (``seed`` for a model)."""
    table = PARAMS[entry_point]
    return table[component] if component in table else table.get("model.*", "seed")


_KNOWN = frozenset(c for table in LEGACY.values() for c in table if c != "model.*")
_REPLICATE_RULE = (
    "the legacy policy does not derive replicates: the parent seed comes from the plan and the fan-out "
    "stays in its owner (stability algorithm_seed + i, nulls and calibration SeedSequence.spawn, "
    "motif_sensitivity range(n_perm), engine_bench stress range(n_seeds))")


def _models() -> dict:
    from grnbench.infer.base import REGISTRY
    from grnbench.plugins import ensure_plugins_loaded
    ensure_plugins_loaded()
    return REGISTRY


def _check_component(component) -> None:
    if not isinstance(component, str) or not component:
        raise ValueError(f"component must be a non-empty name, got {component!r}")
    if component in _KNOWN:
        return
    if component.startswith("model.") and component[len("model."):] in _models():
        return
    raise ValueError(f"unknown seed component {component!r}; known: {', '.join(sorted(_KNOWN))} "
                     "and model.<a registered method>")


def _runner_seed(name: str) -> tuple[bool, object]:
    """(takes a seed, its default) for a registered method, read from the runner's signature, or from
    the base a wrapper hands ``**params`` to (the motif and cisTarget wrappers)."""
    runner = _models()[name].runner
    candidates = [runner]
    try:
        base = inspect.getclosurevars(runner).nonlocals.get("base")
    except TypeError:
        base = None
    if callable(base):
        candidates.append(base)
    for fn in candidates:
        try:
            p = inspect.signature(fn).parameters.get("seed")
        except (TypeError, ValueError):
            continue
        if p is not None:
            return True, p.default
    return False, None


def _has_seed(component: str) -> bool:
    """False for a model whose runner takes no seed (aracne3, the correlation methods)."""
    return not component.startswith("model.") or _runner_seed(component[len("model."):])[0]


def _library_default(component: str) -> int:
    library = LEGACY["library"]
    if component in library:
        value = library[component]
    elif component.startswith("model."):
        value = _runner_seed(component[len("model."):])[1]
    else:
        value = None
    if not _is_count(value):
        raise ValueError(f"{component}: no seed is passed here and its default is {value!r}, so it has "
                         "no effective seed (kwargs() is {})")
    return int(value)


@dataclass(frozen=True)
class SeedPlan:
    """Which seed each component gets, at one entry point (required under every policy).

    ``legacy``: today's seeds, from :data:`LEGACY`; ``master``: the master seed everywhere, ``master + task``
    for a task cell; ``derived``: :func:`derive` of the component's knob (and of the task for a task cell).
    A ``replicate`` is derived under ``master`` and ``derived``; under ``legacy`` it is refused, because
    today the loop that fans out owns its seeds.
    """
    master: int = DEFAULT_SEED
    policy: Policy = "legacy"
    entry_point: str | None = None

    def __post_init__(self):
        if self.policy not in POLICIES:
            raise ValueError(f"unknown seed policy {self.policy!r}; choose from {', '.join(POLICIES)}")
        object.__setattr__(self, "master", _as_seed(self.master, "master seed"))
        if self.entry_point is None:
            raise ValueError(f"a seed plan needs an entry point (a task cell needs it to know it is one); "
                             f"known: {', '.join(LEGACY)}")
        if self.entry_point not in LEGACY:
            raise ValueError(f"unknown entry point {self.entry_point!r}; known: {', '.join(LEGACY)}")

    def as_dict(self, *, task: int | None = None) -> dict:
        """For manifests: the plan and the effective seed of every component at the entry point.

        ``master`` is null under ``legacy`` when no cell takes the user's seed; a component without an
        effective seed (or a task cell without ``task``) records null. ``model.*`` lists the seeded models.
        """
        table = LEGACY[self.entry_point]
        names = [c for c in table if c != "model.*"]
        if "model.*" in table:
            names += [c for c in LEGACY["library"] if c.startswith("model.")]
        seeds = {}
        for component in names:
            is_task = self._cell(component) is TASK
            try:
                seeds[component] = None if is_task and task is None else self.seed_for(
                    component, task=task if is_task else None)
            except ValueError:
                seeds[component] = None
        uses_master = self.policy != "legacy" or MASTER in table.values()
        return {"master": self.master if uses_master else None, "policy": self.policy,
                "entry_point": self.entry_point, "task": None if task is None else _as_seed(task, "task"),
                "seeds": seeds}

    def seed_for(self, component: str, *, task: int | None = None, replicate: int | None = None) -> int:
        """The effective seed of ``component``: what its random calls really use.

        A cell that passes no seed gives the callee's default (grnboost2 42, regdiffusion 0). A model
        whose runner takes no seed (aracne3) raises: it has no effective seed. A loop that fans out
        (replicates) asks here with ``replicate=``; :meth:`kwargs` takes none.
        """
        cell = self._cell(component)
        task = self._task(component, cell, task)
        if replicate is not None:
            replicate = _as_seed(replicate, "replicate")
        if not _has_seed(component):
            raise ValueError(f"{component} takes no seed (its runner has no seed parameter), so it has no "
                             "effective seed; kwargs() passes what the entry point passes today")
        if self.policy != "legacy":
            return self._fresh(component, cell, task, replicate)
        if replicate is not None:
            raise ValueError(_REPLICATE_RULE)
        if self.entry_point == "library" or cell is UNSEEDED:
            return _library_default(component)
        if cell is MASTER:
            return self.master
        if cell is TASK:
            return task
        return cell

    def kwargs(self, component: str, *, task: int | None = None) -> dict:
        """What to pass to the call, under the component's real parameter name (:data:`PARAMS`):
        ``{}`` where no seed is passed today, else ``{name: value}``. No replicate: a loop that fans
        out uses :meth:`seed_for`. A component without a parameter (a literal at the call site, a module
        constant) raises."""
        cell = self._cell(component)
        task = self._task(component, cell, task)
        name = param(component, self.entry_point)
        if name is None:
            raise ValueError(f"{component} has no seed parameter at {self.entry_point} (a literal at the call "
                             "site or a module constant), so there is nothing to pass; seed_for() gives its value")
        if self.policy != "legacy":
            return {name: self._fresh(component, cell, task, None)}
        if self.entry_point == "library" or cell is UNSEEDED:
            return {}
        return {name: self.master if cell is MASTER else task if cell is TASK else cell}

    def _cell(self, component: str) -> Cell:
        """The LEGACY cell of ``component`` at the entry point; refuses unknown names."""
        _check_component(component)
        table = LEGACY[self.entry_point]
        if component in table:
            return table[component]
        if component.startswith("model."):
            if self.entry_point == "library":
                return None                    # a runner without a library row: its signature decides
            if "model.*" in table:
                return table["model.*"]
        raise ValueError(f"entry point {self.entry_point!r} does not run {component!r}; "
                         f"it runs: {', '.join(table)}")

    def _task(self, component: str, cell: Cell, task) -> int | None:
        if cell is TASK:
            if task is None:
                raise ValueError(f"{component} uses the task's seed at {self.entry_point}: pass task=")
            return _as_seed(task, "task")
        if task is not None:
            raise ValueError(f"{component} has a fixed seed at {self.entry_point}, not the task's: "
                             "do not pass task=")
        return None

    def _fresh(self, component: str, cell: Cell, task: int | None, replicate: int | None) -> int:
        """The seed under ``master`` or ``derived``."""
        name = knob(component, self.entry_point)
        if cell is TASK:
            if replicate is not None:
                raise ValueError(f"{component} fans out by task here: pass task or replicate, not both")
            if self.policy == "derived":
                return derive(self.master, name, task)
            if self.master + task > SEED_MAX:
                raise ValueError(f"master + task = {self.master + task} is above 2**32 - 1")
            return self.master + task
        if replicate is not None:
            return derive(self.master, name, replicate)
        return self.master if self.policy == "master" else derive(self.master, name)
