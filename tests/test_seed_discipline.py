"""Every random call in the engine and the benchmarks names its seed (issue #111).

An AST scan of ``src/grnbench`` and ``benchmarks`` finds the calls that draw random numbers without
an explicit seed: ``np.random.default_rng()`` with no argument, numpy's global random state
(``np.random.seed``, ``np.random.choice``, ...), the stdlib ``random`` module, scikit-learn
estimators, splitters and utilities without ``random_state``, scipy ``svds`` without ``v0`` and
``bootstrap`` / ``permutation_test`` without ``random_state``, scanpy and umap calls without
``random_state``, networkx layouts without ``seed``, torch draws without a ``generator`` and pandas
``.sample`` without ``random_state``. Names are resolved through the file's imports and plain
assignments (``mk = np.random.default_rng``), so ``rng.choice`` on a seeded generator is not
numpy's global ``choice``.

``ALLOWED`` freezes the offenders found when this test was written. A new offender fails, and so
does an allowed entry that no longer exists: the list can only shrink. A seed variable or a literal
counts as seeded (making literal defaults explicit is a later step).

Known misses (the scan reads code, it does not run it):

- a seed parameter that is None at run time (``def f(seed=None): default_rng(seed)``) counts as seeded;
- time- or entropy-based seeds (``default_rng(time.time_ns())``) count as seeded;
- methods of objects the scan cannot type: ``tensor.uniform_()``, ``stats.norm(0, 1).rvs()``, a torch
  draw that relies on ``torch.manual_seed``;
- literal seeds that ignore the run's seed, e.g. the scored-cell subsample ``default_rng(0)`` in
  ``benchmarks/engine_bench.py`` (``score_candidates``, line 339 today) and ``benchmarks/real_bench.py``
  (``score_models``, line 253 today);
- a seed passed only through ``**kwargs`` is not seen (the call is reported, strictly);
- an alias assigned inside one function is applied to the whole file.
"""
from __future__ import annotations

import ast
import sys
from collections import Counter
from functools import lru_cache
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCANNED = ("src/grnbench", "benchmarks")

# (file, enclosing function, call) of each unseeded random call on 2026-10-18. Only ever remove entries.
ALLOWED: list[tuple[str, str, str]] = []

# numpy generator and bit-generator constructors -> the parameter that seeds them (also the first positional)
NUMPY_CONSTRUCTORS = {"default_rng": "seed", "RandomState": "seed", "SeedSequence": "entropy",
                      "Generator": "bit_generator", "PCG64": "seed", "PCG64DXSM": "seed",
                      "MT19937": "seed", "Philox": "seed", "SFC64": "seed"}
SKLEARN = {"PCA", "KMeans", "MiniBatchKMeans", "TruncatedSVD", "TSNE", "NMF", "FastICA", "IsolationForest",
           "train_test_split", "ShuffleSplit", "StratifiedShuffleSplit", "GroupShuffleSplit", "RepeatedKFold",
           "RepeatedStratifiedKFold", "shuffle", "resample"}
SKLEARN_PREFIXES = ("RandomForest", "ExtraTrees", "GradientBoosting", "HistGradientBoosting")
SKLEARN_IF_SHUFFLED = {"KFold", "StratifiedKFold"}                     # random only with shuffle=True
SCIPY_RESAMPLING = {"bootstrap", "permutation_test"}
SCANPY_MODULES = {"pp", "tl", "preprocessing", "tools", "external"}     # not pl: plots draw nothing random
SCANPY = {"pca", "scrublet", "scrublet_simulate_doublets", "neighbors", "umap", "leiden", "louvain", "tsne",
          "draw_graph", "subsample", "sample", "downsample_counts", "score_genes"}
NETWORKX = {"spring_layout", "fruchterman_reingold_layout", "louvain_communities", "louvain_partitions",
            "asyn_lpa_communities", "asyn_fluidc"}
TORCH = {"rand", "randn", "randint", "randperm", "rand_like", "randn_like", "randint_like", "normal",
         "bernoulli", "multinomial", "poisson"}


def _given(call: ast.Call, key: str, positional: bool = False) -> bool:
    """``key`` is passed (as a keyword, or first positional when ``positional``) and is not a literal None."""
    kw = {k.arg: k.value for k in call.keywords if k.arg}
    v = kw.get(key, call.args[0] if positional and call.args else None)
    return v is not None and not (isinstance(v, ast.Constant) and v.value is None)


def _on(call: ast.Call, key: str) -> bool:
    """Keyword ``key`` is passed and is not a literal False, 0 or None."""
    v = next((k.value for k in call.keywords if k.arg == key), None)
    return v is not None and not (isinstance(v, ast.Constant) and not v.value)


def _check(name: str, call: ast.Call) -> tuple[str, str | None] | None:
    """``(family, why unseeded or None)`` for a random call resolved to ``name``; None for any other call."""
    parts = name.split(".")
    last = parts[-1]
    if parts[:2] == ["numpy", "random"] and len(parts) >= 3:
        if last not in NUMPY_CONSTRUCTORS:
            return "numpy", "numpy's global random state"
        return "numpy", None if _given(call, NUMPY_CONSTRUCTORS[last], positional=True) else "no seed"
    if parts[0] == "random":
        if last == "Random" and _given(call, "x", positional=True):
            return "random", None
        return "random", "the stdlib random module"
    if parts[0] == "sklearn" and (last in SKLEARN or last.startswith(SKLEARN_PREFIXES)):
        return "sklearn", None if _given(call, "random_state") else "no random_state"
    if parts[0] == "sklearn" and last in SKLEARN_IF_SHUFFLED and _on(call, "shuffle"):
        return "sklearn", None if _given(call, "random_state") else "shuffle without random_state"
    if parts[0] == "sklearn" and last == "silhouette_score" and _given(call, "sample_size"):
        return "sklearn", None if _given(call, "random_state") else "sample_size without random_state"
    if parts[0] == "scipy" and last == "svds":
        seeded = any(_given(call, k) for k in ("v0", "random_state", "rng"))
        return "scipy", None if seeded else "no v0 or random_state"
    if parts[0] == "scipy" and last in SCIPY_RESAMPLING:
        seeded = any(_given(call, k) for k in ("random_state", "rng"))
        return "scipy", None if seeded else "no random_state"
    if parts[0] == "scanpy" and len(parts) >= 3 and parts[1] in SCANPY_MODULES and last in SCANPY:
        return "scanpy", None if _given(call, "random_state") else "no random_state"
    if parts[0] == "umap" and last == "UMAP":
        return "umap", None if _given(call, "random_state") else "no random_state"
    if parts[0] == "networkx" and (last in NETWORKX or "random" in last):
        return "networkx", None if _given(call, "seed") else "no seed"
    if parts[0] == "torch" and last == "seed":
        return "torch", "torch.seed() sets a non-deterministic seed"
    if parts[0] == "torch" and last in TORCH:
        return "torch", None if _given(call, "generator") else "no generator"
    return None


def _aliases(tree: ast.AST) -> dict[str, str]:
    """Local name -> module path, from every import of the file (imports inside functions included) and
    from plain assignments of an imported name (``mk = np.random.default_rng``)."""
    names = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                names[a.asname or a.name.split(".")[0]] = a.name if a.asname else a.name.split(".")[0]
        elif isinstance(node, ast.ImportFrom) and node.module and not node.level:
            for a in node.names:
                names[a.asname or a.name] = f"{node.module}.{a.name}"
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            parts = _dotted(node.value)
            if parts and parts[0] in names:
                names[node.targets[0].id] = ".".join([names[parts[0]], *parts[1:]])
    return names


def _dotted(node: ast.AST) -> list[str] | None:
    if isinstance(node, ast.Name):
        return [node.id]
    if isinstance(node, ast.Attribute):
        head = _dotted(node.value)
        return None if head is None else head + [node.attr]
    return None


class _Scan(ast.NodeVisitor):
    def __init__(self, aliases: dict[str, str]):
        self.aliases, self.scope, self.found = aliases, [], []

    def _enter(self, node) -> None:
        self.scope.append(node.name)
        self.generic_visit(node)
        self.scope.pop()

    visit_FunctionDef = visit_AsyncFunctionDef = visit_ClassDef = _enter

    def visit_Call(self, node: ast.Call) -> None:
        parts = _dotted(node.func)
        hit = None
        if parts and parts[0] in self.aliases:
            hit = _check(".".join([self.aliases[parts[0]], *parts[1:]]), node)
        elif isinstance(node.func, ast.Attribute) and node.func.attr == "sample":     # DataFrame / Series
            hit = "pandas", None if _given(node, "random_state") else "no random_state"
        if hit:
            self.found.append((".".join(self.scope) or "<module>", ast.unparse(node.func), *hit))
        self.generic_visit(node)


def scan_source(source: str) -> list[tuple[str, str, str, str | None]]:
    """``(enclosing function, call, family, why unseeded or None)`` of every random call in ``source``."""
    tree = ast.parse(source)
    scan = _Scan(_aliases(tree))
    scan.visit(tree)
    return scan.found


@lru_cache(maxsize=1)
def scan_repository() -> tuple[tuple[str, ...], tuple[tuple[str, str, str, str, str | None], ...]]:
    files, found = [], []
    for root in SCANNED:
        for path in sorted((ROOT / root).rglob("*.py")):
            rel = path.relative_to(ROOT).as_posix()
            files.append(rel)
            found += [(rel, *hit) for hit in scan_source(path.read_text(encoding="utf-8"))]
    return tuple(files), tuple(found)


CONTROL = '''
import random
import numpy as np
import scanpy as sc
import networkx as nx
import torch
import umap
from numpy.random import default_rng
from scipy import stats
from scipy.sparse.linalg import svds
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingRegressor, IsolationForest, RandomForestRegressor
from sklearn.metrics import silhouette_score
from sklearn.model_selection import KFold, ShuffleSplit
from sklearn.utils import resample, shuffle

mk = np.random.default_rng


def unseeded(X, G, df, f):
    np.random.default_rng()
    default_rng(None)
    np.random.seed(0)
    np.random.choice(X, 3)
    random.shuffle(X)
    PCA(n_components=2)
    RandomForestRegressor(random_state=None)
    HistGradientBoostingRegressor()
    IsolationForest()
    KFold(5, shuffle=True)
    ShuffleSplit(3)
    shuffle(X)
    resample(X)
    silhouette_score(X, X, sample_size=100)
    svds(X, k=2)
    stats.bootstrap((X,), np.mean)
    stats.permutation_test((X, X), f)
    sc.pp.pca(X)
    sc.tl.umap(X)
    sc.pp.downsample_counts(X, counts_per_cell=10)
    sc.tl.score_genes(X, ["A"])
    umap.UMAP()
    nx.spring_layout(G)
    torch.rand(3)
    df.sample(5)
    mk()


def seeded(X, G, df, f, seed, stream, v0):
    rng = np.random.default_rng(seed)
    default_rng(stream)
    np.random.default_rng(0).choice(X, 3)
    rng.permutation(X)
    random.Random(seed)
    PCA(n_components=2, random_state=seed)
    RandomForestRegressor(random_state=0)
    HistGradientBoostingRegressor(random_state=seed)
    IsolationForest(random_state=seed)
    KFold(5)
    KFold(5, shuffle=False)
    KFold(5, shuffle=True, random_state=seed)
    ShuffleSplit(3, random_state=seed)
    shuffle(X, random_state=seed)
    resample(X, random_state=seed)
    silhouette_score(X, X)
    silhouette_score(X, X, sample_size=100, random_state=seed)
    svds(X, k=2, v0=v0)
    svds(X, k=2, random_state=seed)
    stats.bootstrap((X,), np.mean, random_state=seed)
    stats.permutation_test((X, X), f, rng=seed)
    sc.pp.pca(X, random_state=seed)
    sc.pp.normalize_total(X)
    sc.pl.umap(X)
    sc.pp.downsample_counts(X, counts_per_cell=10, random_state=seed)
    sc.tl.score_genes(X, ["A"], random_state=seed)
    umap.UMAP(random_state=seed)
    nx.spring_layout(G, seed=seed)
    torch.rand(3, generator=torch.Generator().manual_seed(seed))
    df.sample(5, random_state=seed)
    mk(seed)
'''


def test_the_scan_flags_unseeded_calls_and_passes_seeded_ones():
    found = scan_source(CONTROL)
    assert [(f, c) for f, c, _, why in found if why] == [("unseeded", c) for c in (
        "np.random.default_rng", "default_rng", "np.random.seed", "np.random.choice", "random.shuffle", "PCA",
        "RandomForestRegressor", "HistGradientBoostingRegressor", "IsolationForest", "KFold", "ShuffleSplit",
        "shuffle", "resample", "silhouette_score", "svds", "stats.bootstrap", "stats.permutation_test", "sc.pp.pca",
        "sc.tl.umap", "sc.pp.downsample_counts", "sc.tl.score_genes", "umap.UMAP", "nx.spring_layout", "torch.rand",
        "df.sample", "mk")]
    assert [c for f, c, _, why in found if f == "seeded"] == [
        "np.random.default_rng", "default_rng", "np.random.default_rng", "random.Random", "PCA",
        "RandomForestRegressor", "HistGradientBoostingRegressor", "IsolationForest", "KFold", "ShuffleSplit",
        "shuffle", "resample", "silhouette_score", "svds", "svds", "stats.bootstrap", "stats.permutation_test",
        "sc.pp.pca", "sc.pp.downsample_counts", "sc.tl.score_genes", "umap.UMAP", "nx.spring_layout", "torch.rand",
        "df.sample", "mk"]


def test_the_scan_reads_this_repository():
    """The scan sees both folders and resolves this code's import styles (lazy imports inside functions)."""
    files, found = scan_repository()
    assert all(any(f.startswith(root + "/") for f in files) for root in SCANNED)
    seen = {family for *_, family, why in found if not why}
    assert {"numpy", "sklearn", "scanpy", "umap", "networkx"} <= seen, (
        f"the scan recognises seeded calls of {sorted(seen)} only: a broken resolver, or the code stopped using "
        "a library (then drop it from this set)")


def test_no_new_random_call_without_a_seed():
    _, found = scan_repository()
    why = {(f, func, call): reason for f, func, call, _, reason in found if reason}
    now = Counter((f, func, call) for f, func, call, _, reason in found if reason)
    allowed = Counter(ALLOWED)
    new, gone = now - allowed, allowed - now
    assert not new, ("New random call(s) without an explicit seed. Pass the run's seed (seed=, random_state=, "
                     "generator=, v0=); do not add them to ALLOWED:\n"
                     + "\n".join(f"    {e!r},  # {why[e]}" for e in sorted(new.elements())))
    assert not gone, ("ALLOWED lists call(s) that are no longer in the code (seeded, moved or renamed). "
                      "Delete these entries from ALLOWED:\n"
                      + "\n".join(f"    {e!r}," for e in sorted(gone.elements())))


def test_the_ratchet_fails_on_a_new_offender(monkeypatch):
    files, found = scan_repository()
    offender = ("src/grnbench/new.py", "draw", "np.random.rand", "numpy", "numpy's global random state")
    monkeypatch.setattr(sys.modules[__name__], "scan_repository", lambda: (files, (*found, offender)))
    with pytest.raises(AssertionError, match=r"New random call\(s\) without an explicit seed") as err:
        test_no_new_random_call_without_a_seed()
    assert "('src/grnbench/new.py', 'draw', 'np.random.rand'),  # numpy's global random state" in str(err.value)


def test_the_ratchet_fails_on_an_allowed_entry_that_is_gone(monkeypatch):
    stale = ("src/grnbench/gone.py", "draw", "np.random.rand")
    monkeypatch.setattr(sys.modules[__name__], "ALLOWED", [*ALLOWED, stale])
    with pytest.raises(AssertionError, match="no longer in the code") as err:
        test_no_new_random_call_without_a_seed()
    assert repr(stale) in str(err.value)
