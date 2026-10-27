"""Every registered model runs once through the registry on tiny data (issue #114).

The other method tests check what each method is for; this module checks that no registered
method is left without a run. Each name goes through ``run_method`` (normalization, regulator
resolution, the adapter, ``standardize_edges``) and must return an edge table that

- gives edges to every regulator of the list it was given, and to no other gene;
- carries a ranking: the weights are not all equal (an all-zero grid fails), sorted high to low;
- for the cisTarget methods, reports motif-supported edges (``cistarget_supported`` > 0);
- for every method that reads the motif rankings, ranks TFA's planted motif targets well above
  the genes TFA drives just as strongly but without a motif (a disabled motif step fails).

Only the methods whose dependencies are outside ``requirements/lock-py311.txt`` may skip:
aracne3 (the native ARACNe3 program, named by ``GRNBENCH_ARACNE3_BIN`` or on PATH), regdiffusion and regdiffusion_motif
(the regdiffusion package and torch). ``SKIP_REASONS`` names those dependencies, and a method
skips only when the registry says it cannot run here and one of them is really absent; the skip
message is "not tested here: <the registry's reason> (waits for: <dependency>)". Everything else
fails: a method that cannot run while its declared dependencies are present or that declares
none (arboreto, pyarrow and scikit-learn are in the lock file), a probe that crashed
("availability probe failed: ..."), a registered method this module does not smoke-run, and a
synthetic cisTarget folder the registry cannot read (an error in the fixture).
"""
from __future__ import annotations

import importlib.util
import os
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from grnbench.infer import EDGE_COLUMNS, run_method
from grnbench.infer import base as infer_base
from grnbench.infer import cistarget as ct

LOCK_FILE = Path(__file__).resolve().parents[1] / "requirements" / "lock-py311.txt"

# every method smoke-run here, with settings that keep the run small ({} = its defaults)
SMOKE = {
    "spearman_signed": {}, "pearson_signed": {}, "mutual_info": {}, "clr": {}, "ppcor": {},
    "lasso": {}, "random_forest": {"n_estimators": 10},
    "genie3": {"n_estimators": 10}, "grnboost2": {"n_estimators": 20},
    "aracne3": {"timeout_s": 120},
    "regdiffusion": {"n_steps": 20},
    "motif_prior": {}, "pearson_motif": {}, "pearson_cistarget": {},
    "grnboost2_motif": {"n_estimators": 20}, "grnboost2_cistarget": {"n_estimators": 20},
    "regdiffusion_motif": {"n_steps": 20},
    "bicor": {}, "wgcna_tom": {}, "propr_rho": {}, "aracne_dpi": {}, "mrnet": {},
    "tigress": {"n_resamples": 10}, "bbsr": {}, "pcnet": {},
    "extra_trees": {"n_estimators": 10}, "gradient_boosting": {"n_estimators": 10},
    "community": {},
}
# the only methods allowed to skip, and the dependencies outside the lock file each one waits
# for: ("module", m) = an importable Python module, ("executable", x) = a program on PATH,
# ("env_path", v) = an environment variable naming an existing file. arboreto, pyarrow and
# scikit-learn are in the lock file, so genie3, grnboost2 and the motif methods are not listed:
# they fail when their dependency is missing.
SKIP_REASONS = {
    "aracne3": [("env_path", "GRNBENCH_ARACNE3_BIN")],         # or on PATH: then the registry runs it
    "regdiffusion": [("module", "regdiffusion"), ("module", "torch")],
    "regdiffusion_motif": [("module", "regdiffusion"), ("module", "torch")],
}
# these read the motif rankings, so they run on the genes of the synthetic cisTarget database
MOTIF = {"motif_prior", "pearson_motif", "grnboost2_motif", "regdiffusion_motif",
         "pearson_cistarget", "grnboost2_cistarget"}
# these prune with cisTarget and report how many edges the motifs support
CISTARGET = {"pearson_cistarget", "grnboost2_cistarget"}

REGISTERED = sorted(n for n, s in infer_base.REGISTRY.items() if s.origin == "builtin")


def _absent(dep: tuple[str, str]) -> bool:
    kind, name = dep
    if kind == "module":
        return importlib.util.find_spec(name) is None
    if kind == "executable":
        return shutil.which(name) is None
    if kind == "env_path":
        value = os.environ.get(name, "")
        return not value or not Path(value).exists()
    raise ValueError(f"unknown dependency kind {kind!r}")


def _describe(dep: tuple[str, str]) -> str:
    kind, name = dep
    return {"module": f"Python module {name}", "executable": f"{name} on PATH",
            "env_path": f"{name} naming an existing file"}[kind]


# --------------------------------------------------------------------------- #
# the synthetic cisTarget folder of tests/test_motif_methods.py (copied, not imported: that
# module skips as a whole without pyarrow, and this one must still report every method)
# --------------------------------------------------------------------------- #

N_BACKGROUND = 400        # genome-wide genes in the fake database
TARGETS = [f"T{i}" for i in range(25)]      # TFA's true targets: motif near them
OTHERS = [f"G{i}" for i in range(25)]       # co-expressed with TFA, no motif


def _write_resources(root, rng):
    genes = ["TFA", "TFB"] + TARGETS + OTHERS + [f"BG{i}" for i in range(N_BACKGROUND)]
    motifs = ["mA"] + [f"m{i}" for i in range(60)]
    ann = pd.DataFrame({"#motif_id": ["mA", "m1"], "gene_name": ["TFA", "TFB"]})
    ann.to_csv(root / "motifs-test.hgnc-m0.001-o0.0.tbl", sep="\t", index=False)
    try:
        import pyarrow.feather as pf
    except ImportError:          # the registry's probe then names the missing piece
        return
    ranks = np.array([rng.permutation(len(genes)) for _ in motifs], dtype=np.int16)
    # motif mA ranks TFA's targets at the very top of the genome, OTHERS at the bottom
    order = TARGETS + [g for g in genes if g not in TARGETS and g not in OTHERS] + OTHERS
    ranks[0] = np.array([order.index(g) for g in genes], dtype=np.int16)
    table = pd.DataFrame(ranks, columns=genes)
    table["motifs"] = motifs
    pf.write_feather(table, str(root / "hg38_test.genes_vs_motifs.rankings.feather"))


@pytest.fixture()
def cistarget_resources(tmp_path, monkeypatch):
    _write_resources(tmp_path, np.random.default_rng(0))
    monkeypatch.setenv(ct.ENV_DIR, str(tmp_path))
    ct.load_annotation.cache_clear()
    if importlib.util.find_spec("pyarrow") is not None:     # the module the registry's probe imports
        ok, why = ct.cistarget_availability()
        assert ok, f"the synthetic cisTarget folder is broken, not a missing dependency: {why}"
    return tmp_path


def _coexpressed_adata() -> ad.AnnData:
    """TFA drives TARGETS and OTHERS equally; only TARGETS carry TFA's motif (normalized scale)."""
    rng = np.random.default_rng(1)
    n = 300
    tfa = rng.normal(size=n)
    cols = {"TFA": tfa, "TFB": rng.normal(size=n)}
    for g in TARGETS + OTHERS:
        cols[g] = tfa + rng.normal(scale=0.7, size=n)
    for i in range(30):
        cols[f"BG{i}"] = rng.normal(size=n)
    frame = pd.DataFrame(cols)
    return ad.AnnData(frame.to_numpy(dtype=np.float32) + 5.0, var=pd.DataFrame(index=frame.columns))


# --------------------------------------------------------------------------- #
# the smoke run
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("name", REGISTERED)
def test_every_registered_method_runs_through_the_registry(name, cascade_adata, request):
    if name not in SMOKE:
        pytest.fail(f"'{name}' is registered but this module does not smoke-run it: add it to SMOKE")
    spec = infer_base.REGISTRY[name]
    assert spec.origin == "builtin", f"'{name}' was replaced by a {spec.origin} method"
    if name in MOTIF:            # before the probe: it reads the resource folder from the environment
        request.getfixturevalue("cistarget_resources")
        adata, tfs, normalized = _coexpressed_adata(), ["TFA", "TFB"], True
    else:
        adata, tfs, normalized = cascade_adata, ["TF_A", "G_B"], False
    reason = infer_base.unavailability_reason(spec)
    if reason is not None:
        if reason.startswith("availability probe failed"):
            pytest.fail(f"'{name}': the registry's availability probe crashed ({reason})")
        missing = [d for d in SKIP_REASONS.get(name, []) if _absent(d)]
        if not missing:
            why = ("its declared dependencies are all present" if name in SKIP_REASONS else
                   "it may not skip: its dependencies are in the lock file")
            pytest.fail(f"'{name}' cannot run here ({reason}) and {why}")
        pytest.skip(f"not tested here: {reason} (waits for: {', '.join(map(_describe, missing))})")

    edges, prov = run_method(name, adata, tf_list=tfs, assume_normalized=normalized, **SMOKE[name])

    assert list(edges.columns[:3]) == EDGE_COLUMNS
    assert set(edges.columns) <= set(EDGE_COLUMNS) | {"direction"}
    assert set(edges["source"]) == set(tfs), (
        f"'{name}': every regulator of the list gets edges, and no other gene does "
        f"(sources: {sorted(set(edges['source']))})")
    assert set(edges["target"]) <= set(adata.var_names), "targets must be genes of the data"
    assert (edges["source"] != edges["target"]).all(), "self-loops must be dropped"
    assert not edges.duplicated(["source", "target"]).any(), "each pair is reported once"
    w = edges["weight"].to_numpy()
    assert w.dtype.kind == "f" and np.isfinite(w).all(), "weights are finite floats"
    assert len(np.unique(w)) > 1, f"'{name}' gave every edge the same weight: no ranking"
    assert np.all(np.diff(w) <= 0), "edges must arrive sorted by weight, high to low"
    if "direction" in edges.columns:
        assert set(edges["direction"].unique()) <= {-1, 0, 1}
    assert prov["method"] == name and prov["origin"] == "builtin"
    assert prov["n_edges"] == len(edges)
    assert prov["regulators"] == "tf_list"

    if name in CISTARGET:
        supported = prov.get("method_info", {}).get("cistarget_supported")
        assert supported is not None and supported > 0, (
            f"'{name}': cisTarget pruning supported no edge (cistarget_supported={supported})")
    if name in MOTIF:
        tfa = edges.loc[edges["source"] == "TFA", "target"].tolist()     # in weight order
        position = {g: i for i, g in enumerate(tfa)}
        lost = [g for g in TARGETS + OTHERS if g not in position]
        assert not lost, f"'{name}': TFA has no edge to {lost}"
        at_targets = np.mean([position[g] for g in TARGETS])
        at_others = np.mean([position[g] for g in OTHERS])
        assert at_targets + 10 < at_others, (
            f"'{name}': TFA's motif targets sit at mean position {at_targets:.1f} and the genes it "
            f"drives without a motif at {at_others:.1f}: the motif step did not rank them first")


def test_every_registered_method_is_smoke_run_and_only_off_lock_methods_may_skip():
    """A newly registered method fails here until it gets a smoke run; a skip names dependencies
    that are outside the lock file."""
    registered = {n for n, s in infer_base.REGISTRY.items() if s.origin == "builtin"}
    assert not registered - set(SMOKE), f"registered without a smoke run: {sorted(registered - set(SMOKE))}"
    assert not set(SMOKE) - registered, f"smoke-run but no longer registered: {sorted(set(SMOKE) - registered)}"
    assert set(REGISTERED) == registered, "the parametrized smoke run must see the whole registry"
    assert set(SKIP_REASONS) <= set(SMOKE)
    assert MOTIF <= set(SMOKE) and CISTARGET <= MOTIF
    pinned = {line.split("==")[0].strip().lower().replace("_", "-")
              for line in LOCK_FILE.read_text(encoding="utf-8").splitlines()
              if line.strip() and not line.lstrip().startswith("#")}
    for name, deps in SKIP_REASONS.items():
        assert deps, f"'{name}': a skip needs the dependencies it waits for"
        for kind, dep in deps:
            assert kind in {"module", "executable", "env_path"}, f"'{name}': unknown kind {kind!r}"
            assert not (kind == "module" and dep.lower().replace("_", "-") in pinned), (
                f"'{name}' waits for {dep}, which the lock file pins: it must run, not skip")
