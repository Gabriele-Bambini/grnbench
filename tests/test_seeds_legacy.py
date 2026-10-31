"""Every value of ``seeds.LEGACY`` is read back from the code that sets it today (issue #126).

What is traced, per entry point:

1. the value: a default (``inspect.signature``), a module constant, a parser default, or a literal at
   the call site (AST; ``ui/app.py`` and ``benchmarks/*.py`` are only parsed, never imported);
2. the path: each call passes no seed, or exactly ``args.seed`` / ``int(args.seed)`` / ``task["seed"]``,
   and each function that takes a seed hands exactly that name on (run_prep, execute_plan,
   screen_candidates, screen_processed, the model runners, community, the stability reruns);
3. the leaves: in every function holding a traced seed, every RNG constructor (``default_rng``,
   ``SeedSequence``, ``RandomState``) has an argument and every known RNG consumer (``CONSUMERS``) has its
   seed keyword, and all of them use that function's seed; the model internals (arboreto, the
   scikit-learn fits, TIGRESS, RegDiffusion, community), the simulator (``default_rng(scn.seed)``,
   ``10_000 + seed``) and the task enumerations (seeds from 0) included;
4. coverage: every seeded file under ``src/grnbench`` and ``benchmarks`` is read by a tracer or listed in
   ``seeds.NOT_COVERED``, whose registered values are pinned too, and nothing seeds a global RNG.

Known misses: a random call through an object these checks do not name (a new RNG consumer outside
``CONSUMERS``), and arithmetic on the seed inside a callee that is not traced. Line numbers live in the
comments of ``seeds.py``; these checks match structure, and evaluate the stability fan-out instead of
matching its shape, so a harmless rewrite passes.
"""
import ast
import dataclasses
import inspect
from pathlib import Path

import pytest

from grnbench import seeds
from grnbench.seeds import MASTER, TASK, UNSEEDED

ROOT = Path(__file__).resolve().parents[1]
APP = "src/grnbench/ui/app.py"
LOGIC = "src/grnbench/ui/logic.py"
REAL = "benchmarks/real_bench.py"
ENGINE = "benchmarks/engine_bench.py"
STUDY = "benchmarks/candidate_study.py"
MOTIF = "benchmarks/motif_sensitivity.py"
BEELINE = "benchmarks/beeline_backtest.py"

PREP = ("prep.sampling", "prep.doublets", "prep.metacells")
SCREEN = ("screen.split_reference", "screen.cell_cap", "screen.batch_leakage", "screen.edge_stability")
EVAL = ("eval.null", "eval.stability", "eval.significance")
RNG_CTORS = ("default_rng", "SeedSequence", "RandomState")
# callables that draw random numbers, and the keyword that seeds them
CONSUMERS = {
    "PCA": "random_state", "KMeans": "random_state", "MiniBatchKMeans": "random_state",
    "silhouette_score": "random_state", "pca": "random_state", "scrublet": "random_state",
    "Scrublet": "random_state", "UMAP": "random_state", "umap": "random_state", "leiden": "random_state",
    "louvain": "random_state", "train_test_split": "random_state",
    "RandomForestRegressor": "random_state", "RandomForestClassifier": "random_state",
    "ExtraTreesRegressor": "random_state", "GradientBoostingRegressor": "random_state",
    "Lasso": "random_state", "spring_layout": "seed", "RegDiffusionTrainer": "seed",
    "stable_subsample": "seed", "stable_cell_order": "seed",   # screen.py: cells by a hash of seed and barcode (D-024)
}
TRACED: set[str] = set()          # every source file a tracer read


# -- helpers --------------------------------------------------------------------------- #

def _tree(rel: str) -> ast.Module:
    TRACED.add(rel)
    return ast.parse((ROOT / rel).read_text(encoding="utf-8"))


def _src(module) -> str:
    path = Path(inspect.getsourcefile(module)).resolve()
    try:
        rel = path.relative_to(ROOT).as_posix()
    except ValueError:
        pytest.fail(f"{module.__name__} is imported from {path}, not from this checkout ({ROOT}): the trace "
                    "would compare this checkout's sources with another one's values. Install this checkout "
                    "editable (pip install -e .) or set PYTHONPATH=src.")
    TRACED.add(rel)
    return rel


def _func(rel: str, name: str) -> ast.FunctionDef:
    found = [n for n in ast.walk(_tree(rel)) if isinstance(n, ast.FunctionDef) and n.name == name]
    assert len(found) == 1, f"{rel}: expected one function '{name}', found {len(found)}"
    return found[0]


def _callee(call: ast.Call) -> str | None:
    f = call.func
    return f.id if isinstance(f, ast.Name) else f.attr if isinstance(f, ast.Attribute) else None


def _calls(node, name: str, *, required: bool = True) -> list[ast.Call]:
    found = [n for n in ast.walk(node) if isinstance(n, ast.Call) and _callee(n) == name]
    assert found or not required, f"no call to {name}() where one is expected"
    return found


def _kw(call: ast.Call, name: str):
    return next((k.value for k in call.keywords if k.arg == name), None)


def _is_name(node, name: str) -> bool:
    return isinstance(node, ast.Name) and node.id == name


def _names(node) -> set[str]:
    return {n.id for n in ast.walk(node) if isinstance(n, ast.Name)} if node is not None else set()


def _strip_int(node):
    if isinstance(node, ast.Call) and _callee(node) == "int" and len(node.args) == 1 and not node.keywords:
        return node.args[0]
    return node


def _is_param(node, name: str = "seed") -> bool:
    """Exactly ``name`` or ``int(name)``."""
    return _is_name(_strip_int(node), name)


def _reads_args(node, attr: str) -> bool:
    """Exactly ``args.<attr>`` or ``int(args.<attr>)``: no arithmetic on the user's value."""
    node = _strip_int(node)
    return isinstance(node, ast.Attribute) and node.attr == attr and _is_name(node.value, "args")


def _positional(fn) -> list[str]:
    """Names that a positional argument can fill: from a function, or from its AST definition."""
    if isinstance(fn, ast.FunctionDef):
        return [a.arg for a in fn.args.posonlyargs + fn.args.args]
    return [p.name for p in inspect.signature(fn).parameters.values()
            if p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD)]


def _seed_arg(call: ast.Call, fn, param: str = "seed"):
    """What the call passes as ``param`` of ``fn`` (keyword or positional), else None."""
    node = _kw(call, param)
    pos = _positional(fn)
    if node is None and param in pos and len(call.args) > pos.index(param):
        node = call.args[pos.index(param)]
    return node


def _omits_seed(call: ast.Call, fn=None) -> bool:
    """No seed keyword, no ``**`` splat that could carry one, and no positional seed."""
    if any(k.arg in (None, "seed") for k in call.keywords):
        return False
    return fn is None or _seed_arg(call, fn) is None


def _passes(call: ast.Call, fn, name: str = "seed") -> bool:
    return _is_name(_seed_arg(call, fn), name)


def _default(fn, param: str = "seed"):
    value = inspect.signature(fn).parameters[param].default
    assert value is not inspect.Parameter.empty, f"{fn.__qualname__}: '{param}' has no default"
    return value


def _ast_default(fdef: ast.FunctionDef, param: str = "seed"):
    args = fdef.args.posonlyargs + fdef.args.args
    defaults = dict(zip([a.arg for a in args][len(args) - len(fdef.args.defaults):], fdef.args.defaults))
    defaults.update({a.arg: d for a, d in zip(fdef.args.kwonlyargs, fdef.args.kw_defaults) if d is not None})
    return ast.literal_eval(defaults[param])


def _literal_cap_seed(fdef: ast.FunctionDef):
    """The one literal seed of every ``stable_subsample(..., seed=<literal>)`` in ``fdef`` (D-024 cell caps)."""
    values = {ast.literal_eval(_kw(c, "seed")) for c in _calls(fdef, "stable_subsample")}
    assert len(values) == 1, f"{fdef.name}: stable_subsample seeds {values}"
    return values.pop()


def _literal_rng_seed(fdef: ast.FunctionDef):
    """The one literal seed of every ``default_rng(<literal>)`` in ``fdef``."""
    values = {ast.literal_eval(c.args[0]) for c in _calls(fdef, "default_rng")}
    assert len(values) == 1, f"{fdef.name}: default_rng seeds {values}"
    return values.pop()


def _leaf(node, *, param=("seed",), keywords=(), allow=()) -> None:
    """Every RNG constructor in ``node`` has an argument, every known consumer its seed keyword, and
    each of them uses exactly one of ``param`` (or ``int(...)`` of it).

    ``keywords`` adds other keywords to check (``seed=`` handed on); ``allow`` names other variables
    that may seed an RNG (``child_seed`` from a spawn)."""
    label = getattr(node, "name", "module")
    uses = []
    for call in (n for n in ast.walk(node) if isinstance(n, ast.Call)):
        name = _callee(call)
        if name in RNG_CTORS:
            assert call.args or call.keywords, f"{label}: {name}() without a seed draws from the OS"
            uses.append(call.args[0] if call.args else call.keywords[0].value)
        elif name in CONSUMERS:
            value = _kw(call, CONSUMERS[name])
            assert value is not None, f"{label}: {ast.unparse(call.func)}() without {CONSUMERS[name]}="
            uses.append(value)
        uses += [k.value for k in call.keywords if k.arg in keywords]
    uses = [u for u in uses if not any(_is_name(u, a) for a in allow)]
    assert uses, f"{label}: no RNG uses {param}"
    bad = [ast.unparse(u) for u in uses if not any(_is_param(u, p) for p in param)]
    assert not bad, f"{label}: RNG seeded with {bad}, not {param}"


def _is_task_seed(node) -> bool:
    """``task["seed"]``, also inside ``int(...)``."""
    node = _strip_int(node)
    return (isinstance(node, ast.Subscript) and _is_name(node.value, "task")
            and isinstance(node.slice, ast.Constant) and node.slice.value == "seed")


def _argv_seed(node: ast.List, attr: str = "seed") -> bool:
    """An argv list holds ``"--seed", str(args.seed)`` exactly."""
    items = node.elts
    for a, b in zip(items, items[1:]):
        if isinstance(a, ast.Constant) and a.value == f"--{attr}":
            return (isinstance(b, ast.Call) and _callee(b) == "str" and len(b.args) == 1
                    and _reads_args(b.args[0], attr))
    return False


def _module_assign(rel: str, name: str):
    return next(n.value for n in _tree(rel).body if isinstance(n, ast.Assign)
                and any(_is_name(t, name) for t in n.targets))


def _hands_on(fdef: ast.FunctionDef, callee: str, name: str = "seed") -> None:
    """Every call to ``callee`` in ``fdef`` gets ``seed=<name>``."""
    assert all(_is_name(_kw(c, "seed"), name) for c in _calls(fdef, callee)), f"{fdef.name} -> {callee}"


def _seeds_from_zero(fdef: ast.FunctionDef) -> None:
    """Each task's "seed" comes from a ``range(n)``: seeds 0..n-1, as registered."""
    seed_vars = {v.id for d in ast.walk(fdef) if isinstance(d, ast.Dict)
                 for k, v in zip(d.keys, d.values)
                 if isinstance(k, ast.Constant) and k.value == "seed" and isinstance(v, ast.Name)}
    loops = [n.iter for n in ast.walk(fdef) if isinstance(n, (ast.comprehension, ast.For))
             and isinstance(n.target, ast.Name) and n.target.id in seed_vars]
    assert loops, f"{fdef.name}: no seed loop"
    for it in loops:
        assert (isinstance(it, ast.Call) and _callee(it) == "range" and len(it.args) == 1
                and not it.keywords), f"{fdef.name}: seeds from {ast.unparse(it)}, not range(n)"


# -- how a seed travels through the shared code ------------------------------------------ #

def _prep_from(seed) -> dict:
    """run_prep hands its seed to the sampling rng, scrublet and metacells."""
    from grnbench import hygiene, pipeline, strategies
    rel = _src(pipeline)
    fdef = _func(rel, "run_prep")
    rngs = _calls(fdef, "default_rng")
    assert len(rngs) == 1 and _is_name(rngs[0].args[0], "seed"), "run_prep: one rng = default_rng(seed)"
    _leaf(fdef, keywords=("random_state", "seed"))
    assert all(_passes(c, hygiene.detect_doublets) for c in _calls(fdef, "detect_doublets")), \
        "run_prep: detect_doublets(..., seed=seed)"
    metacells_branch = [n for n in ast.walk(fdef) if isinstance(n, ast.If)
                        and isinstance(n.test, ast.Compare) and any(
                            isinstance(c, ast.Constant) and c.value == "metacells" for c in n.test.comparators)]
    assert len(metacells_branch) == 1, "run_prep: one `if chosen == \"metacells\"` branch"
    assert any(_is_name(_kw(c, "seed"), "seed") for c in _calls(metacells_branch[0], "update")), \
        "run_prep: params.update(..., seed=seed) in the metacells branch"
    _leaf(_func(_src(hygiene), "detect_doublets"))
    _leaf(_func(_src(strategies), "metacells"))
    return {c: seed for c in PREP}


def _execute_plan_from(seed) -> dict:
    """execute_plan hands its seed to every run_prep."""
    from grnbench import autoprep
    fdef = _func(_src(autoprep), "execute_plan")
    assert any(_is_name(_kw(c, "seed"), "seed") for c in _calls(fdef, "update")), \
        "execute_plan: kwargs.update(..., seed=seed)"
    assert all(any(k.arg is None for k in c.keywords) for c in _calls(fdef, "run_prep")), \
        "execute_plan: run_prep(..., **kwargs)"
    return _prep_from(seed)


def _screening_from(seed) -> dict:
    """screen_candidates hands its seed to the split, to execute_plan and to every measurement;
    the models it runs (edge stability, proxies) get no seed."""
    from grnbench import autoprep, screen
    rel = _src(screen)
    fdef = _func(rel, "screen_candidates")
    for fn in (screen.split_reference, autoprep.execute_plan, screen.screen_processed):
        assert all(_passes(c, fn) for c in _calls(fdef, fn.__name__)), f"screen_candidates -> {fn.__name__}"
    proc = _func(rel, "screen_processed")
    caps = _calls(proc, "stable_subsample")
    assert len(caps) == 1 and _is_name(_kw(caps[0], "seed"), "seed"), "screen_processed: cap stable_subsample(seed=seed)"
    for fn in (screen.batch_leakage, screen.edge_stability):
        assert all(_passes(c, fn) for c in _calls(proc, fn.__name__)), f"screen_processed -> {fn.__name__}"
    _leaf(proc, keywords=("seed",))
    for name in ("split_reference", "batch_leakage", "edge_stability"):
        _leaf(_func(rel, name))
    for name in ("edge_stability", "screen_processed"):
        assert all(_omits_seed(c) for c in _calls(_func(rel, name), "run_method")), f"{name}: model seeded"
    return {**{c: seed for c in SCREEN}, **_execute_plan_from(seed), "model.*": UNSEEDED}


def _profile() -> dict:
    """The profile is called without a seed, so its own defaults apply, and its leaves use them."""
    from grnbench import profile
    rel = _src(profile)
    assert all(_omits_seed(c, profile._nonzero_values)
               for c in _calls(_func(rel, "technology_from_counts"), "_nonzero_values"))
    assert all(_omits_seed(c, profile.batch_silhouette)
               for c in _calls(_func(rel, "build_profile"), "batch_silhouette"))
    _leaf(_func(rel, "_nonzero_values"))
    _leaf(_func(rel, "batch_silhouette"))
    return {"profile.nonzero_sample": _default(profile._nonzero_values),
            "profile.batch_silhouette": _default(profile.batch_silhouette)}


def _unseeded_models(fdef: ast.FunctionDef) -> dict:
    assert all(_omits_seed(c) for c in _calls(fdef, "run_method")), f"{fdef.name}: run_method gets a seed"
    return {"model.*": UNSEEDED}


def _engine_scoring() -> dict:
    """engine_bench.score_candidates: literal cell cap; models get only method_params, never a seed."""
    score = _func(ENGINE, "score_candidates")
    for call in _calls(score, "run_method"):
        splats = [k.value for k in call.keywords if k.arg is None]
        assert _kw(call, "seed") is None and len(splats) == 1 and _names(splats[0]) == {"method_params", "m"}
    params = ast.literal_eval(_module_assign(ENGINE, "BENCH_METHOD_PARAMS"))
    assert not any("seed" in inner for inner in params.values())
    return {"score.cell_cap": _literal_cap_seed(score), "model.*": UNSEEDED}


def _viz_seed() -> int:
    """viz.SEED; every RNG and consumer in viz.py uses SEED or a ``seed`` parameter defaulting to SEED."""
    from grnbench import viz
    tree = _tree(_src(viz))
    for f in (n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)):
        if "seed" in _positional(f):
            assert _is_name(_default_node(f, "seed"), "SEED"), f.name
    _leaf(tree, param=("SEED", "seed"), keywords=("random_state", "seed"))
    return viz.SEED


def _default_node(fdef: ast.FunctionDef, param: str):
    args = fdef.args.posonlyargs + fdef.args.args
    return dict(zip([a.arg for a in args][len(args) - len(fdef.args.defaults):], fdef.args.defaults))[param]


def _simulator() -> int:
    """simulate() draws from default_rng(scn.seed); random_scenario from default_rng(10_000 + seed)."""
    from grnbench import simulate
    rel = _src(simulate)
    rngs = sorted(ast.unparse(c.args[0]) for c in _calls(_tree(rel), "default_rng"))
    assert rngs == ["10000 + seed", "scn.seed"], rngs
    assert [ast.unparse(c.args[0]) for c in _calls(_func(rel, "simulate"), "default_rng")] == ["scn.seed"]
    assert [ast.unparse(c.args[0]) for c in _calls(_func(rel, "random_scenario"), "default_rng")] == \
        ["10000 + seed"]
    return next(f.default for f in dataclasses.fields(simulate.Scenario) if f.name == "seed")


def _model_internals() -> None:
    """The seed reaches the random calls inside the adapters."""
    from grnbench.infer import aracne3, classical, community, random_baseline, regdiffusion_adapter, regression, tree_methods
    rel = _src(tree_methods)
    arboreto = _func(rel, "_run_arboreto")
    for call in _calls(arboreto, "create_graph"):
        value = _kw(call, "seed")      # int(seed); arboreto's DEMON_SEED only for an explicit None
        assert (isinstance(value, ast.IfExp) and _is_param(value.body)
                and ast.unparse(value.test) == "seed is not None") or _is_param(value), ast.unparse(value)
    for name in ("genie3_runner", "grnboost2_runner"):
        assert all(_passes(c, arboreto) for c in _calls(_func(rel, name), "_run_arboreto")), name
    rel = _src(classical)
    for name in ("_fit_rf", "_fit_lasso"):
        _leaf(_func(rel, name))
    loop = _func(rel, "_per_target_importance")
    assert all(any(k.arg is None and _is_name(k.value, "kw") for k in c.keywords) for c in _calls(loop, "fit_one"))
    for name in ("random_forest", "lasso"):
        _hands_on(_func(rel, name), "_per_target_importance")
    rel = _src(regression)
    for name in ("tigress", "_fit_extra_trees", "_fit_gradient_boosting"):
        _leaf(_func(rel, name))
    for name in ("extra_trees", "gradient_boosting"):
        _hands_on(_func(rel, name), "_per_target_importance")
    _leaf(_func(_src(regdiffusion_adapter), "regdiffusion_runner"))
    _leaf(_func(_src(random_baseline), "random_baseline"))          # the negative control (D-037)
    argv = [n for n in ast.walk(_func(_src(aracne3), "run_aracne3")) if isinstance(n, ast.List)
            and any(isinstance(e, ast.Constant) and e.value == "--seed" for e in n.elts)]
    assert len(argv) == 1, "run_aracne3: one command line with --seed"
    flag = next(i for i, e in enumerate(argv[0].elts) if isinstance(e, ast.Constant) and e.value == "--seed")
    assert ast.unparse(argv[0].elts[flag + 1]) == "str(int(seed))", "run_aracne3 hands its seed to --seed"
    comm = _func(_src(community), "community")
    extra = [n.value for n in ast.walk(comm) if isinstance(n, ast.Assign)
             and any(_is_name(t, "extra") for t in n.targets)]
    assert len(extra) == 1 and any(isinstance(d, ast.Dict) and [getattr(k, "value", None) for k in d.keys] == ["seed"]
                                   and _is_name(d.values[0], "seed") for d in ast.walk(extra[0])), \
        "community hands its seed to its members"
    assert all(any(k.arg is None and _is_name(k.value, "extra") for k in c.keywords)
               for c in _calls(comm, "run_method"))


# -- one tracer per entry point ---------------------------------------------------------- #

def _library_models() -> dict:
    """Every registry runner with a seed parameter, and the wrappers that hand ``**params`` to one."""
    from grnbench.infer import cistarget
    from grnbench.infer.base import REGISTRY
    _model_internals()
    for wrapper in ("with_motif_prior", "with_cistarget"):
        inner = _calls(_func(_src(cistarget), wrapper), "base")
        assert all(any(k.arg is None for k in c.keywords) for c in inner), f"{wrapper}: base(..., **params)"
    out = {}
    for name, spec in REGISTRY.items():
        if spec.origin != "builtin":
            continue
        params = inspect.signature(spec.runner).parameters
        if "seed" in params:
            out[f"model.{name}"] = _default(spec.runner)
            continue
        base = inspect.getclosurevars(spec.runner).nonlocals.get("base")
        if base is not None and "seed" in inspect.signature(base).parameters:
            out[f"model.{name}"] = _default(base)
    return out


def _library() -> dict:
    from grnbench import hygiene, pipeline, plugins, profile, screen, strategies
    from grnbench.eval import nulls, significance, stability
    command = [f for f in ast.walk(_func(_src(plugins), "command_method_spec"))
               if isinstance(f, ast.FunctionDef) and f.name == "runner"]
    assert len(command) == 1
    nm = _func(_src(nulls), "null_metrics")
    _leaf(nm, allow=("child_seed",))
    assert all(_is_name(c.args[0], "seed") for c in _calls(nm, "SeedSequence"))
    _leaf(_func(_src(significance), "gold_thinning_sensitivity"))
    _leaf(_func(_src(stability), "subsample_stability"))
    _leaf(_func(STUDY, "degrade"))
    return {
        "prep.sampling": _default(pipeline.run_prep),
        "prep.doublets": _default(hygiene.detect_doublets),
        "prep.metacells": _default(strategies.metacells),
        **_profile(),
        "screen.split_reference": _default(screen.split_reference),
        "screen.cell_cap": _default(screen.screen_processed),
        "screen.batch_leakage": _default(screen.batch_leakage),
        "screen.edge_stability": _default(screen.edge_stability),
        **_library_models(),
        "plugin.command": _ast_default(command[0]),
        "eval.null": _default(nulls.null_metrics),
        "eval.stability": _default(stability.subsample_stability),
        "eval.significance": _default(significance.gold_thinning_sensitivity),
        "stability.model": _default(stability.subsample_stability, "algorithm_seed"),
        "viz": _viz_seed(),
        "simulate": _simulator(),
        "degrade": _ast_default(_func(STUDY, "degrade")),
    }


def _cli_master(module) -> None:
    assert module._build_parser().get_default("seed") == seeds.DEFAULT_SEED, f"{module.__name__} --seed"


def _cli_prep() -> dict:
    from grnbench import cli
    _cli_master(cli)
    calls = _calls(_func(_src(cli), "main"), "run_prep")
    assert all(_reads_args(_kw(c, "seed"), "seed") for c in calls), "cli main: run_prep(seed=args.seed)"
    return _prep_from(MASTER)


def _cli_infer() -> dict:
    from grnbench import cli_infer
    _cli_master(cli_infer)
    calls = _calls(_func(_src(cli_infer), "main"), "run_method")
    assert all(_reads_args(_kw(c, "seed"), "seed") for c in calls), "infer: run_method(seed=args.seed)"
    _library_models()                      # the runners hand the seed on
    return {"model.*": MASTER}


def _stability_fan_out(node) -> bool:
    """The rerun seed is algorithm_seed + i with vary_algorithm_seed, algorithm_seed without."""
    if node is None or not _names(node) <= {"algorithm_seed", "i", "vary_algorithm_seed"}:
        return False
    code = compile(ast.Expression(node), "<stability seed>", "eval")
    return all(eval(code, {"__builtins__": {}}, {"algorithm_seed": 100, "i": i, "vary_algorithm_seed": vary})
               == (100 + i if vary else 100) for i in range(4) for vary in (True, False))


def _cli_eval() -> dict:
    from grnbench import cli_eval
    from grnbench.eval import nulls, significance, stability
    _cli_master(cli_eval)
    parser = cli_eval._build_parser()
    assert parser.get_default("vary_algorithm_seed") is False
    tree = _tree(_src(cli_eval))
    for fn, name in ((nulls.null_metrics, "null_metrics"),
                     (significance.gold_thinning_sensitivity, "gold_thinning_sensitivity"),
                     (stability.subsample_stability, "subsample_stability")):
        assert all(_reads_args(_seed_arg(c, fn), "seed") for c in _calls(tree, name)), name
    assert all(_reads_args(_kw(c, "algorithm_seed"), "algorithm_seed")
               for c in _calls(tree, "subsample_stability"))
    rerun = _calls(_func(_src(stability), "subsample_stability"), "run_method")
    assert all(_stability_fan_out(_kw(c, "seed")) for c in rerun), "stability: algorithm_seed + i"
    return {**{c: MASTER for c in EVAL}, "stability.model": parser.get_default("algorithm_seed")}


def _cli_run() -> dict:
    from grnbench import cli_run
    rel = _src(cli_run)
    _cli_master(cli_run)
    main = _func(rel, "main")
    argvs = {t.id: n.value for n in ast.walk(main) if isinstance(n, ast.Assign)
             for t in n.targets if isinstance(t, ast.Name) and isinstance(n.value, ast.List)}
    for name in ("prep_argv", "infer_argv", "eval_argv"):
        assert _argv_seed(argvs[name]), f"cli_run: {name} holds '--seed', str(args.seed)"
    strings = {n.value for n in ast.walk(_tree(rel)) if isinstance(n, ast.Constant) and isinstance(n.value, str)}
    assert "--algorithm-seed" not in strings, "cli_run now forwards --algorithm-seed"
    return {**_cli_prep(), **_cli_infer(), **_cli_eval()}


def _interface() -> dict:
    from grnbench import autoprep
    from grnbench.ui import logic
    app = _tree(APP)
    for call in _calls(app, "execute_plan") + _calls(_func(LOGIC, "quick_run"), "execute_plan"):
        assert _omits_seed(call, autoprep.execute_plan) or _passes(call, autoprep.execute_plan, "DEFAULT_SEED"), \
            "the interface cleans with execute_plan's default (omitted, or DEFAULT_SEED since #121)"
    figures = (logic.subsample_cells, logic.embedding_frame, logic.expression_frame, logic.network_figure)
    viz_seed = {_default(f) for f in figures}
    assert len(viz_seed) == 1, viz_seed
    for name in ("embedding_frame", "expression_frame"):
        assert all(_passes(c, logic.subsample_cells) for c in _calls(_func(LOGIC, name), "subsample_cells"))
    for name in ("subsample_cells", "embedding_frame", "network_figure"):
        _leaf(_func(LOGIC, name), keywords=("random_state", "seed"))
    for fn in figures[1:]:
        assert all(_omits_seed(c, fn) for c in _calls(app, fn.__name__)), fn.__name__
    return {**_profile(), **_execute_plan_from(_default(autoprep.execute_plan)),
            **_unseeded_models(_func(LOGIC, "_run_models_here")), "viz": viz_seed.pop()}


def _interface_screening() -> dict:
    from grnbench import screen
    calls = _calls(_tree(APP), "screen_candidates") + _calls(_func(LOGIC, "quick_run"), "screen_candidates")
    for call in calls:
        assert _omits_seed(call, screen.screen_candidates), "the interface screens with the default seed"
    return _screening_from(_default(screen.screen_candidates))


def _real_bench() -> dict:
    from grnbench import screen
    task = _func(REAL, "run_real_task")
    assert all(_omits_seed(c, screen.screen_candidates) for c in _calls(task, "screen_candidates"))
    assert all(_kw(c, "method_params") is None for c in _calls(task, "score_candidates"))
    models = _func(REAL, "score_models")
    engine = _engine_scoring()
    cap = {_literal_cap_seed(models), engine["score.cell_cap"]}   # its models, and the proxies via engine_bench
    assert len(cap) == 1, cap
    assert all(_omits_seed(c) for c in _calls(_func(REAL, "final_report"), "build_run")), \
        "the report draws the figures with viz.SEED"
    _library_models()
    return {**_profile(), **_screening_from(_default(screen.screen_candidates)),
            "score.cell_cap": cap.pop(), **_unseeded_models(models), "viz": _viz_seed()}


def _real_bench_noref() -> dict:
    from grnbench import autoprep
    calls = _calls(_func(REAL, "run_real_task"), "execute_plan")
    assert all(_omits_seed(c, autoprep.execute_plan) for c in calls)
    return {**_profile(), **_execute_plan_from(_default(autoprep.execute_plan))}


def _engine_bench() -> dict:
    scenario = _func(ENGINE, "_scenario")
    for call in [c for name in ("Scenario", "random_scenario", "with_seed") for c in _calls(scenario, name)]:
        node = _kw(call, "seed") if _callee(call) == "Scenario" else call.args[0]
        assert _is_task_seed(node), ast.unparse(call)
    _simulator()
    for name in ("make_tasks", "make_stress_tasks"):
        _seeds_from_zero(_func(ENGINE, name))
    calls = _calls(_func(ENGINE, "run_task"), "screen_candidates")
    assert all(_is_task_seed(_kw(c, "seed")) for c in calls)
    assert all(ast.unparse(_kw(c, "method_params")) == "BENCH_METHOD_PARAMS"
               for c in _calls(_func(ENGINE, "score_selected_models"), "score_candidates"))
    _library_models()
    return {"simulate": TASK, **_profile(), **_screening_from(TASK), **_engine_scoring()}


def _candidate_study() -> dict:
    task = _func(STUDY, "run_study_task")
    assert all(_is_task_seed(_kw(c, "seed")) for c in _calls(task, "screen_candidates"))
    copy = _func(STUDY, "write_degraded_copy")
    assert all(_is_task_seed(_seed_arg(c, copy)) for c in _calls(task, "write_degraded_copy"))
    assert all(_passes(c, _func(STUDY, "degrade")) for c in _calls(copy, "degrade"))
    _leaf(_func(STUDY, "degrade"))
    # every task make_tasks writes has seed 0, so C7 ran all of the above at 0
    dicts = [d for d in ast.walk(_func(STUDY, "make_tasks")) if isinstance(d, ast.Dict)]
    seeds_written = [ast.literal_eval(v) for d in dicts for k, v in zip(d.keys, d.values)
                     if isinstance(k, ast.Constant) and k.value == "seed"]
    assert seeds_written == [0], seeds_written
    assert all(_kw(c, "method_params") is None for c in _calls(task, "score_candidates"))
    return {"degrade": TASK, **_profile(), **_screening_from(TASK), **_engine_scoring()}


TRACERS = {
    "library": _library,
    "cli_prep": _cli_prep,
    "cli_infer": _cli_infer,
    "cli_eval": _cli_eval,
    "cli_run": _cli_run,
    "interface": _interface,
    "interface_screening": _interface_screening,
    "real_bench": _real_bench,
    "real_bench_noref": _real_bench_noref,
    "engine_bench": _engine_bench,
    "candidate_study": _candidate_study,
}


def test_every_entry_point_has_a_tracer():
    assert set(TRACERS) == set(seeds.LEGACY)


@pytest.mark.parametrize("entry", sorted(TRACERS))
def test_legacy_matches_the_code(entry):
    assert dict(seeds.LEGACY[entry]) == TRACERS[entry]()


# -- coverage -------------------------------------------------------------------------- #

def _seeded(tree: ast.Module) -> bool:
    """The file takes, passes or hard-codes a seed somewhere, or draws from an RNG."""
    for n in ast.walk(tree):
        if isinstance(n, ast.arg) and "seed" in n.arg:
            return True
        if isinstance(n, ast.keyword) and n.arg and ("seed" in n.arg or n.arg == "random_state"):
            return True
        if isinstance(n, ast.Constant) and isinstance(n.value, str) and n.value.startswith("--") and "seed" in n.value:
            return True
        if isinstance(n, ast.Call) and _callee(n) in RNG_CTORS:
            return True
    return False


def test_seeded_files_are_traced_or_not_covered():
    for tracer in TRACERS.values():
        tracer()
    not_covered = {path for path, _ in seeds.NOT_COVERED}
    assert all((ROOT / p).is_file() for p in not_covered), "stale NOT_COVERED path"
    files = sorted(ROOT.glob("benchmarks/*.py")) + sorted(ROOT.glob("src/grnbench/**/*.py"))
    seeded = {f.relative_to(ROOT).as_posix() for f in files if _seeded(ast.parse(f.read_text(encoding="utf-8")))}
    seeded.discard("src/grnbench/seeds.py")
    assert len(seeded) > 20, "the scan found too few seeded files"
    missing = sorted(seeded - TRACED - not_covered)
    assert not missing, f"seeded files that no tracer reads and NOT_COVERED does not list: {missing}"


def test_not_covered_values_that_registered_runs_used():
    """EXP-20261003-GB-1 and the BEELINE back-test keep their seeds; these pins say when they move."""
    assert ast.literal_eval(_module_assign(MOTIF, "N_PERM")) == 100
    perms = [c for c in _calls(_tree(MOTIF), "range") if "n_perm" in _names(c)]
    assert perms and all(len(c.args) == 1 and _is_name(c.args[0], "n_perm") for c in perms), \
        "permutation seeds 0..n_perm-1"
    interval = _func(MOTIF, "ratio_interval")
    assert _ast_default(interval) == 0
    _leaf(interval)
    assert _literal_rng_seed(_func(MOTIF, "descriptive_tables")) == 0          # the bootstrap seed
    seeds_args = [c for c in _calls(_tree(BEELINE), "add_argument")
                  if c.args and isinstance(c.args[0], ast.Constant) and c.args[0].value == "--seed"]
    assert len(seeds_args) == 1 and ast.literal_eval(_kw(seeds_args[0], "default")) == 42


def test_library_covers_every_seeded_runner():
    from grnbench.infer.base import REGISTRY
    library = seeds.LEGACY["library"]
    for name, spec in REGISTRY.items():
        if spec.origin != "builtin":          # command plugins: plugin.command; others: NOT_COVERED
            continue
        seeded = "seed" in inspect.signature(spec.runner).parameters
        base = inspect.getclosurevars(spec.runner).nonlocals.get("base") if not seeded else None
        seeded = seeded or (base is not None and "seed" in inspect.signature(base).parameters)
        assert (f"model.{name}" in library) == seeded, name
    assert library["model.aracne3"] == 42        # the native program's --seed (#105)


def test_no_global_seeding():
    """No ``np.random.seed``, ``random.seed`` or torch seeding: every stream is explicit."""
    files = sorted(ROOT.glob("src/grnbench/**/*.py")) + sorted(ROOT.glob("benchmarks/*.py"))
    bad = []
    for f in files:
        for n in ast.walk(ast.parse(f.read_text(encoding="utf-8"))):
            if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute) and n.func.attr in (
                    "seed", "manual_seed", "manual_seed_all"):
                bad.append(f"{f.relative_to(ROOT).as_posix()}:{n.lineno} {ast.unparse(n.func)}")
    assert not bad, bad
