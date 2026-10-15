"""Unified GRN inference adapter contract.

Every method consumes a gene-expression frame (cells x genes, log-normalized
by default) plus an optional regulator list, and returns the canonical edge
list: a pandas DataFrame with columns

    source   regulator (TF) gene symbol
    target   regulated gene symbol
    weight   float, higher = stronger evidence of regulation

plus optional ``direction`` (-1/0/1) when the method infers sign. Adapters are
registered in :data:`REGISTRY` exactly like preprocessing strategies; methods
with missing external dependencies stay listed with an explicit unavailability
reason instead of disappearing.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import anndata as ad
import numpy as np
import pandas as pd

from grnbench.ingest import looks_like_counts
from grnbench.strategies import cpm_log1p
from grnbench.tfs import resolve_tf_names

EDGE_COLUMNS = ["source", "target", "weight"]


class MethodError(ValueError):
    pass


@dataclass(frozen=True)
class MethodSpec:
    """One registered inference method.

    ``takes`` selects the calling convention:

    - ``"expr"`` (default) — ``runner(expr_df, tf_names=..., **params)`` where
      ``expr_df`` is cells x genes, normalized per :func:`prepare_adata`;
    - ``"adata"`` — ``runner(adata, tf_list, **params)`` with the *normalized*
      AnnData, so plugins can reach ``obs``/``var``/``layers["counts"]``.

    ``probe`` is an optional callable returning ``(available, reason)``,
    re-evaluated on every listing so environment changes (a jar appearing on
    disk, a docker daemon starting) are picked up without a restart.

    ``default`` says whether the method belongs to the set that ``--methods available``
    runs. Methods added after the first registered experiments are ``default=False`` so
    that command keeps meaning what it meant; ``--methods all`` runs every available one.
    """
    name: str
    description: str
    runner: Callable[..., pd.DataFrame] | None = None
    unavailable: str | None = None  # static reason when not runnable in this env
    takes: str = "expr"             # "expr" | "adata"
    origin: str = "builtin"         # builtin | entry_point | python | command
    probe: Callable[[], tuple[bool, str]] | None = None
    default: bool = True


def require_finite(X) -> None:
    """Refuse a matrix with NaN or infinite entries.

    Correlation methods turn one NaN into an all-zero network and tree methods into an obscure error;
    neither says the data were the problem. ``MethodError`` names how many entries, genes and cells.
    """
    import scipy.sparse as sp
    if X is None:
        return
    if sp.issparse(X):
        coo = X.tocoo()
        data, cells, genes = coo.data, coo.row, coo.col
    else:
        data = np.asarray(X)
        cells = genes = None
    if data.dtype.kind not in "fc":      # integer and boolean matrices cannot hold NaN or inf
        return
    bad = ~np.isfinite(data)
    if not bad.any():
        return
    if cells is None:
        n_cells, n_genes = int(bad.any(axis=1).sum()), int(bad.any(axis=0).sum())
    else:
        n_cells, n_genes = len(set(cells[bad].tolist())), len(set(genes[bad].tolist()))
    raise MethodError(
        f"the expression matrix has {int(bad.sum()):,} non-finite value(s) (NaN or infinity) in {n_genes:,} gene(s) "
        f"and {n_cells:,} cell(s); a correlation method would return an all-zero network without saying why. "
        "Remove or impute them first (for example drop those genes or cells) and run again.")


def prepare_adata(adata: ad.AnnData,
                  assume_normalized: bool = False) -> tuple[ad.AnnData, str]:
    """Return (normalized AnnData, normalization provenance note).

    Counts-backed inputs are CPM/1e4 + log1p normalized **on a copy** (tree and
    correlation methods expect comparable scales); already-normalized matrices
    pass through untouched. A matrix with NaN or infinite entries is refused
    (:func:`require_finite`).
    """
    require_finite(adata.X)
    if assume_normalized:
        return adata, "assume_normalized=True; used X as-is"
    ok, why = looks_like_counts(adata.X)
    if ok:
        return cpm_log1p(adata), f"raw counts detected ({why}) -> cpm_log1p applied on copy"
    return adata, f"X already normalized ({why}); used as-is"


def prepare_expression(adata: ad.AnnData,
                       assume_normalized: bool = False) -> tuple[pd.DataFrame, str]:
    """Return (cells x genes DataFrame, normalization provenance note)."""
    prepped, note = prepare_adata(adata, assume_normalized=assume_normalized)
    return _to_frame(prepped), note


def _to_frame(adata: ad.AnnData) -> pd.DataFrame:
    dense = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    return pd.DataFrame(dense, index=adata.obs_names.astype(str),
                        columns=adata.var_names.astype(str))


def standardize_edges(df: pd.DataFrame,
                      source_col: str,
                      target_col: str,
                      score_col: str,
                      direction_col: str | None = None) -> pd.DataFrame:
    """Project any method output onto the canonical edge schema.

    Self-loops are dropped (count recorded in ``attrs['self_loops_dropped']``);
    output sorted by descending weight.
    """
    out = pd.DataFrame({
        "source": df[source_col].astype(str),
        "target": df[target_col].astype(str),
        "weight": df[score_col].astype(float),
    })
    if direction_col is not None and direction_col in df.columns:
        out["direction"] = df[direction_col].astype(int)
    n_self = int((out["source"] == out["target"]).sum())
    if n_self:
        out = out[out["source"] != out["target"]]
        out.attrs["self_loops_dropped"] = n_self
    return out.sort_values("weight", ascending=False, kind="stable").reset_index(drop=True)


def save_edges(edges: pd.DataFrame, out_dir: str | Path, name: str) -> Path:
    """Write edges as tab-separated gzip (no parquet dependency needed)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{name}.edges.tsv.gz"
    edges.to_csv(path, sep="\t", index=False, compression="gzip")
    return path


def load_edges(path: str | Path) -> pd.DataFrame:
    """Read back an edge list saved by :func:`save_edges` (or plain .tsv/.csv)."""
    p = Path(path)
    sep = "\t" if p.suffix in {".gz", ".tsv"} or "edges.tsv" in p.name else ","
    df = pd.read_csv(p, sep=sep)
    missing = set(EDGE_COLUMNS) - set(df.columns)
    if missing:
        raise MethodError(f"{p.name}: missing edge columns {sorted(missing)}")
    return df.sort_values("weight", ascending=False, kind="stable").reset_index(drop=True)


# --------------------------------------------------------------------------- #
# registry
# --------------------------------------------------------------------------- #

def register_method(spec: MethodSpec, *, replace: bool = False) -> MethodSpec:
    """Add a method to the registry (used by the plugin loader).

    Refuses to shadow an existing name unless ``replace=True`` — a silently
    overridden built-in would make two runs with the same method name mean
    different things.
    """
    existing = REGISTRY.get(spec.name)
    if existing is not None and not replace:
        raise MethodError(
            f"method '{spec.name}' is already registered (origin: {existing.origin}); "
            "rename the plugin or pass replace=True"
        )
    REGISTRY[spec.name] = spec
    return spec


def unavailability_reason(spec: MethodSpec) -> str | None:
    """None when the method can run here, else the human-readable reason."""
    if spec.probe is not None:
        try:
            ok, why = spec.probe()
        except Exception as e:  # noqa: BLE001 - a broken probe is unavailability, not a crash
            return f"availability probe failed: {type(e).__name__}: {e}"
        if not ok:
            return why or "unavailable (probe returned False without a reason)"
        return None
    return spec.unavailable


def list_methods() -> list[dict]:
    """Every registered method with its current availability (plugins included)."""
    from grnbench.plugins import ensure_plugins_loaded
    ensure_plugins_loaded()
    entries = []
    for spec in REGISTRY.values():
        reason = unavailability_reason(spec)
        entries.append({
            "name": spec.name,
            "description": spec.description,
            "available": reason is None,
            "unavailable_reason": reason,
            "origin": spec.origin,
            "default": spec.default,
        })
    return entries


def get_method(name: str) -> MethodSpec:
    from grnbench.plugins import ensure_plugins_loaded
    ensure_plugins_loaded()
    if name not in REGISTRY:
        raise MethodError(f"Unknown method '{name}'. Registered: {', '.join(REGISTRY)}")
    spec = REGISTRY[name]
    reason = unavailability_reason(spec)
    if reason is not None:
        raise MethodError(f"Method '{name}' unavailable: {reason}")
    return spec


def run_method(name: str,
               adata: ad.AnnData,
               tf_list: list[str] | None = None,
               *,
               assume_normalized: bool = False,
               **params) -> tuple[pd.DataFrame, dict]:
    """Run one registered inference method; returns (canonical edges, provenance).

    ``tf_list`` restricts regulators when provided; otherwise every gene is a
    candidate regulator.
    """
    spec = get_method(name)
    prepped, norm_note = prepare_adata(adata, assume_normalized=assume_normalized)
    if prepped is adata and spec.origin != "builtin":
        prepped = adata.copy()      # a plugin may edit X in place; the next method must not see that
    effective_tfs = resolve_tf_names(prepped.var_names, tf_list)
    if spec.takes == "adata":
        raw = spec.runner(prepped, effective_tfs, **params)
    else:
        raw = spec.runner(_to_frame(prepped), tf_names=effective_tfs, **params)
    if not isinstance(raw, pd.DataFrame) or raw.shape[1] < 3:
        raise MethodError(
            f"method '{name}' must return a DataFrame with at least "
            f"[source, target, weight]; got {type(raw).__name__} "
            f"with shape {getattr(raw, 'shape', None)}"
        )
    edges = standardize_edges(
        raw,
        source_col=raw.columns[0],
        target_col=raw.columns[1],
        score_col=raw.columns[2],
        direction_col="direction" if "direction" in raw.columns else None,
    )
    prov = {
        "method": name,
        "origin": spec.origin,
        "n_edges": int(len(edges)),
        "self_loops_dropped": int(edges.attrs.get("self_loops_dropped", 0)),
        "normalization": norm_note,
        "regulators": "tf_list" if tf_list else "all_genes",
        "params": {k: str(v) for k, v in params.items()},
    }
    info = dict(raw.attrs)                      # what a method says about its own run (WGCNA's power ...)
    if info:
        prov["method_info"] = info
    return edges, prov


# adapters are imported lazily-safe: they must not raise at import time when
# their optional dependencies are absent (unavailability is data, not an error)
from grnbench.infer.correlation import spearman_signed_abs  # noqa: E402
from grnbench.infer.tree_methods import ARBORETO_UNAVAILABLE, genie3_runner, grnboost2_runner  # noqa: E402
from grnbench.infer.aracne3 import aracne3_availability  # noqa: E402
from grnbench.infer.aracne3 import run_aracne3 as aracne3_runner  # noqa: E402
from grnbench.infer.cistarget import (  # noqa: E402
    cistarget_availability,
    motif_prior,
    with_cistarget,
    with_motif_prior,
)
from grnbench.infer.regdiffusion_adapter import REGDIFFUSION_UNAVAILABLE, regdiffusion_runner  # noqa: E402
from grnbench.infer.association import (  # noqa: E402
    aracne_dpi,
    bicor,
    mrnet,
    propr_rho,
    wgcna_tom,
)
from grnbench.infer.community import community, community_availability  # noqa: E402
from grnbench.infer.regression import (  # noqa: E402
    bbsr,
    extra_trees,
    gradient_boosting,
    pcnet,
    tigress,
)
from grnbench.infer.classical import (  # noqa: E402
    SKLEARN_UNAVAILABLE,
    clr,
    lasso,
    mutual_info,
    pearson_signed,
    ppcor,
    random_forest,
)



def _all_available(static_reason: str | None):
    """Probe for a motif-combined method: its base dependency and the cisTarget resources."""
    def probe() -> tuple[bool, str]:
        if static_reason is not None:
            return False, static_reason
        return cistarget_availability()
    return probe


REGISTRY: dict[str, MethodSpec] = {
    spec.name: spec
    for spec in [
        MethodSpec("spearman_signed",
                   "Fast deterministic baseline: |Spearman| TF-target correlations with sign",
                   spearman_signed_abs),
        MethodSpec("genie3",
                   "GENIE3: random-forest feature importance (arboreto)",
                   genie3_runner,
                   unavailable=ARBORETO_UNAVAILABLE),
        MethodSpec("grnboost2",
                   "GRNBoost2: gradient-boosting with early stopping (arboreto)",
                   grnboost2_runner,
                   unavailable=ARBORETO_UNAVAILABLE),
        MethodSpec("aracne3",
                   "ARACNe3 MI-based network (Java jar; java on PATH + ARACNE3_JAR env)",
                   aracne3_runner,
                   probe=aracne3_availability),
        # numpy/scipy only - always available, so a benchmark can always run
        MethodSpec("pearson_signed",
                   "Linear baseline: |Pearson| TF-target correlations with sign",
                   pearson_signed),
        MethodSpec("mutual_info",
                   "Plug-in mutual information on equal-frequency bins",
                   mutual_info),
        MethodSpec("clr",
                   "CLR: MI z-scored against each gene's own background (Faith et al. 2007)",
                   clr),
        # scikit-learn only - no arboreto, no dask cluster
        MethodSpec("ppcor",
                   "Partial correlation from a Ledoit-Wolf shrunk precision matrix",
                   ppcor,
                   unavailable=SKLEARN_UNAVAILABLE),
        MethodSpec("random_forest",
                   "Per-target random-forest importance (GENIE3's estimator, sklearn only)",
                   random_forest,
                   unavailable=SKLEARN_UNAVAILABLE),
        MethodSpec("lasso",
                   "Per-target L1 regression over the regulator set",
                   lasso,
                   unavailable=SKLEARN_UNAVAILABLE),
        MethodSpec("regdiffusion",
                   "RegDiffusion: diffusion-model adjacency (Zhu & Slonim 2024; torch, CPU ok)",
                   regdiffusion_runner,
                   unavailable=REGDIFFUSION_UNAVAILABLE),
        # sequence prior from the public cisTarget motif rankings (GRNBENCH_CISTARGET_DIR);
        # no curated TF-target collection is used, so CollecTRI scoring stays non-circular
        MethodSpec("motif_prior",
                   "Sequence-only control: TF-target pairs ranked by the TF's best motif rank near the gene",
                   motif_prior,
                   probe=cistarget_availability),
        MethodSpec("pearson_motif",
                   "Pearson x motif prior: rank product of |Pearson| and motif rank",
                   with_motif_prior(pearson_signed),
                   probe=cistarget_availability),
        MethodSpec("grnboost2_motif",
                   "GRNBoost2 x motif prior: rank product of importance and motif rank",
                   with_motif_prior(grnboost2_runner),
                   probe=_all_available(ARBORETO_UNAVAILABLE)),
        MethodSpec("regdiffusion_motif",
                   "RegDiffusion x motif prior: rank product of |adjacency| and motif rank",
                   with_motif_prior(regdiffusion_runner),
                   probe=_all_available(REGDIFFUSION_UNAVAILABLE)),
        MethodSpec("pearson_cistarget",
                   "Pearson + cisTarget pruning (pySCENIC step 2): motif-supported edges ranked first",
                   with_cistarget(pearson_signed),
                   probe=cistarget_availability),
        MethodSpec("grnboost2_cistarget",
                   "GRNBoost2 + cisTarget pruning, i.e. pySCENIC without AUCell",
                   with_cistarget(grnboost2_runner),
                   probe=_all_available(ARBORETO_UNAVAILABLE)),

        # --- widely used methods added later (docs/models_catalog.md). default=False: the set
        # --methods available runs is unchanged; ``--methods all`` or an explicit name runs them.
        MethodSpec("bicor",
                   "Biweight midcorrelation (WGCNA, Langfelder & Horvath 2012); Pearson for genes with MAD 0",
                   bicor, default=False),
        MethodSpec("wgcna_tom",
                   "WGCNA topological overlap: soft power from the scale-free fit (Zhang & Horvath 2005)",
                   wgcna_tom, default=False),
        MethodSpec("propr_rho",
                   "Proportionality rho_p (Quinn et al. 2017) on the normalised log expression",
                   propr_rho, default=False),
        MethodSpec("aracne_dpi",
                   "ARACNe-style: binned MI, G-test threshold, DPI pruning over TF triangles (Margolin 2006); not the ARACNe-AP code",
                   aracne_dpi, default=False),
        MethodSpec("mrnet",
                   "MRNET (Meyer et al. 2007): maximum relevance, minimum redundancy per target on binned MI",
                   mrnet, default=False),
        MethodSpec("tigress",
                   "TIGRESS-style (Haury et al. 2012): LARS stability selection with random regulator weights",
                   tigress, unavailable=SKLEARN_UNAVAILABLE, default=False),
        MethodSpec("bbsr",
                   "BBSR-style (Inferelator): CLR preselection and BIC best-subset regression per target",
                   bbsr, default=False),
        MethodSpec("pcnet",
                   "pcNet (scTenifoldNet, Osorio et al. 2020): principal-component regression per gene",
                   pcnet, default=False),
        MethodSpec("extra_trees",
                   "Extremely randomised trees importance per target (GENIE3's ET option; scikit-learn only)",
                   extra_trees, unavailable=SKLEARN_UNAVAILABLE, default=False),
        MethodSpec("gradient_boosting",
                   "Gradient-boosted trees importance per target (GRNBoost2's idea; scikit-learn only, no dask)",
                   gradient_boosting, unavailable=SKLEARN_UNAVAILABLE, default=False),
        MethodSpec("community",
                   "Mean-rank aggregate of spearman_signed, clr, ppcor, lasso, mrnet, bicor (Marbach et al. 2012)",
                   community, takes="adata", probe=community_availability, default=False),
    ]
}
