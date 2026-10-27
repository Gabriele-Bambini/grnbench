"""Preprocessing strategy registry.

Every strategy consumes a counts-backed AnnData and returns a new AnnData with
normalized ``X`` while preserving raw counts in ``layers["counts"]``. Strategies
are explicit and logged — the platform never silently mixes recipes. There is
deliberately NO imputation strategy: smoothing (MAGIC, knn-graph imputation)
manufactures spurious correlations and corrupts GRN inference.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Callable

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from grnbench.ingest import looks_like_counts


class StrategyError(ValueError):
    pass


@dataclass(frozen=True)
class StrategySpec:
    name: str
    description: str
    apply: Callable[..., ad.AnnData] | None = None
    requires_counts: bool = True
    unavailable: str | None = None  # reason when not implemented in this env


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #

def _ensure_counts_layer(adata: ad.AnnData) -> None:
    if "counts" in adata.layers:
        return
    ok, why = looks_like_counts(adata.X)
    if not ok:
        raise StrategyError(
            f"No counts layer and X does not look like raw counts ({why}); "
            "cannot normalize safely."
        )
    adata.layers["counts"] = adata.X.copy()


def _group_sum(X, codes: np.ndarray, n_groups: int):
    """Sum rows of X per group code (dense/sparse safe) -> (n_groups, n_vars)."""
    ind = sp.csr_matrix(
        (np.ones(len(codes), dtype=np.float32),
         (codes, np.arange(len(codes)))),
        shape=(n_groups, len(codes)),
    )
    return ind @ X


# --------------------------------------------------------------------------- #
# strategies
# --------------------------------------------------------------------------- #

def raw_counts(adata: ad.AnnData, **_) -> ad.AnnData:
    """Pass-through: keep raw counts as X (for methods normalizing internally)."""
    out = adata.copy()
    _ensure_counts_layer(out)
    out.X = out.layers["counts"].copy()
    return out


def cpm_log1p(adata: ad.AnnData, target_sum: float = 1e4, **_) -> ad.AnnData:
    """CPM-style depth normalization to ``target_sum`` followed by log1p."""
    out = adata.copy()
    _ensure_counts_layer(out)
    sc.pp.normalize_total(out, target_sum=target_sum)
    sc.pp.log1p(out)
    return out


def pearson_residuals(adata: ad.AnnData, theta: float = 100.0, **_) -> ad.AnnData:
    """Analytic Pearson residuals (depth-stable variance normalization).

    Output values can be negative by design; downstream edge scorers must not
    assume non-negativity. Operates on raw counts.
    """
    try:
        from scanpy.experimental import pp as experimental_pp
    except ImportError as e:  # pragma: no cover
        raise StrategyError("pearson_residuals requires scanpy>=1.9") from e
    out = adata.copy()
    _ensure_counts_layer(out)
    out.X = out.layers["counts"].copy()
    experimental_pp.normalize_pearson_residuals(out, theta=theta)
    # A gene with no counts has expected value 0, so scanpy returns 0/0 = NaN for the whole
    # column. Observed = expected = 0: its residual is 0.
    empty = np.asarray(out.layers["counts"].sum(axis=0)).ravel() == 0
    if empty.any():
        X = out.X.toarray() if sp.issparse(out.X) else np.asarray(out.X)
        X[:, empty] = 0
        out.X = X
    return out


def pseudobulk(adata: ad.AnnData, group_key: str, normalize: bool = True,
               target_sum: float = 1e4, min_group_cells: int | None = None,
               **_) -> ad.AnnData:
    """Aggregate counts per group (cell type / sample / replicate).

    Recommended when per-cell data is too sparse for stable correlations.
    Default policy keeps every group and records ``obs['n_cells']``; with
    ``min_group_cells=N``, smaller groups are DROPPED with a loud warning
    listing exactly what was removed.
    """
    if group_key not in adata.obs.columns:
        raise StrategyError(f"group key '{group_key}' not found in obs")
    out_counts = adata.layers.get("counts", adata.X)
    labels = adata.obs[group_key].astype(str).to_numpy()
    uniq = pd.unique(labels)
    codes = pd.Categorical(labels, categories=uniq).codes
    summed = _group_sum(sp.csr_matrix(out_counts), codes, len(uniq))

    order = np.argsort([str(u) for u in uniq])  # deterministic ordering
    summed = summed[order]
    uniq = [uniq[i] for i in order]

    counts_per_group = np.bincount(codes)[order]
    obs = pd.DataFrame({
        group_key: np.asarray(uniq, dtype=object),
        "n_cells": counts_per_group,
    }, index=pd.Index(np.asarray(uniq, dtype=str), name=group_key))
    bulk = ad.AnnData(
        X=summed.tocsr(),
        obs=obs,
        var=adata.var.copy(),
    )
    bulk.layers["counts"] = bulk.X.copy()
    if min_group_cells is not None:
        keep = bulk.obs["n_cells"].to_numpy() >= int(min_group_cells)
        dropped = bulk.obs_names[~keep].tolist()
        if dropped:
            warnings.warn(
                f"dropping {len(dropped)} pseudobulk group(s) below "
                f"{min_group_cells} cells: {dropped[:8]}...",
                stacklevel=2,
            )
            bulk = bulk[keep].copy()
    else:
        tiny = bulk.obs["n_cells"] < 10
        if tiny.any():
            warnings.warn(
                f"{int(tiny.sum())} pseudobulk group(s) have < 10 cells "
                f"(kept, see obs['n_cells']): {bulk.obs_names[tiny].tolist()[:5]}...",
                stacklevel=2,
            )
    if normalize and bulk.n_obs > 0:
        sc.pp.normalize_total(bulk, target_sum=target_sum)
        sc.pp.log1p(bulk)
    return bulk


def cpm_log1p_combat(adata: ad.AnnData, batch_key: str, target_sum: float = 1e4,
                     **_) -> ad.AnnData:
    """``cpm_log1p`` followed by ComBat location/scale batch correction.

    Never auto-recommended: corrected values are model residuals on the log
    scale, can be negative, and remove *all* batch-aligned covariance -
    including any biology that happens to be confounded with batch. Use it
    when batches are technical replicates of the same biology; do not use it
    to erase a condition that was run one batch per condition. Raw counts stay
    in ``layers["counts"]``.
    """
    if batch_key not in adata.obs.columns:
        raise StrategyError(f"batch key '{batch_key}' not found in obs")
    n_missing = int(adata.obs[batch_key].isna().sum())
    if n_missing:
        raise StrategyError(f"ComBat needs a batch label for every cell: {n_missing:,} cell(s) have none in "
                            f"'{batch_key}' (fill them, drop them, or use cpm_log1p_center_batch)")
    batches = adata.obs[batch_key].astype(str)
    sizes = batches.value_counts()
    tiny = sizes[sizes < 2]
    if len(tiny):   # QC can leave a batch with one cell; ComBat cannot fit it
        if len(sizes) - len(tiny) < 2:
            raise StrategyError(
                f"ComBat needs >= 2 batches with >= 2 cells; '{batch_key}' has "
                f"{len(sizes) - len(tiny)} (batches with 1 cell: {sorted(tiny.index)[:5]})")
        warnings.warn(f"ComBat: dropped {int(tiny.sum())} cell(s) from {len(tiny)} batch(es) "
                      f"of '{batch_key}' with < 2 cells: {sorted(tiny.index)[:5]}", stacklevel=2)
        adata = adata[~batches.isin(tiny.index).to_numpy()].copy()
        batches = adata.obs[batch_key].astype(str)
    elif len(sizes) < 2:
        raise StrategyError(f"ComBat needs >= 2 batches; '{batch_key}' has {len(sizes)}")
    out = cpm_log1p(adata, target_sum=target_sum)
    out.obs[batch_key] = pd.Categorical(batches)
    if sp.issparse(out.X):
        out.X = out.X.toarray()
    out.X = np.asarray(out.X, dtype=np.float32)
    sc.pp.combat(out, key=batch_key)
    out.X = np.asarray(out.X, dtype=np.float32)
    return out


def cpm_log1p_center_batch(adata: ad.AnnData, batch_key: str, target_sum: float = 1e4,
                           **_) -> ad.AnnData:
    """``cpm_log1p``, then each gene's per-batch mean is replaced by its global mean.

    Removes between-batch location shifts - the part of gene-gene correlation
    that comes from batches sitting at different levels - and nothing else:
    within-batch variation and scale are untouched. Transparent alternative to
    ComBat (COBRA 2024 shows batch effects also persist in covariance; this does
    not claim to remove those). Only sensible when batches share cell-type
    composition. Raw counts stay in ``layers["counts"]``.
    """
    if batch_key not in adata.obs.columns:
        raise StrategyError(f"batch key '{batch_key}' not found in obs")
    n_missing = int(adata.obs[batch_key].isna().sum())
    batches = adata.obs[batch_key].astype(str).to_numpy()
    if len(set(batches)) < 2:
        raise StrategyError(f"centering needs >= 2 batches; '{batch_key}' has 1")
    sizes = pd.Series(batches).value_counts()
    if (sizes < 2).all():
        raise StrategyError(f"centering needs batches with >= 2 cells; every batch of '{batch_key}' has 1")
    if (sizes < 2).any():
        warnings.warn(f"{int((sizes < 2).sum())} batch(es) of '{batch_key}' have a single cell: centering replaces "
                      "each such cell by the grand mean of every gene, which erases it; drop those cells or merge "
                      "their batch labels", stacklevel=2)
    if n_missing:
        warnings.warn(f"{n_missing:,} cell(s) have no label in '{batch_key}': they are centred together as one "
                      "batch called 'nan'", stacklevel=2)
    out = cpm_log1p(adata, target_sum=target_sum)
    X = out.X.toarray() if sp.issparse(out.X) else np.array(out.X)
    X = X.astype(np.float32)
    grand = X.mean(axis=0)
    for b in pd.unique(batches):
        rows = batches == b
        X[rows] += grand - X[rows].mean(axis=0)
    out.X = X
    return out


def metacells(adata: ad.AnnData, cells_per_metacell: int = 20, batch_key: str | None = None,
              label_key: str | None = None, n_pcs: int = 30, seed: int = 42,
              target_sum: float = 1e4, **_) -> ad.AnnData:
    """Aggregate raw counts of transcriptionally similar cells into metacells.

    Partition: k-means on the PCA of log-normalized HVGs (a scratch copy),
    ``n_cells / cells_per_metacell`` groups, built inside each batch so a
    metacell never mixes batches. Counts are summed (integer), then
    depth-normalized and log1p like ``cpm_log1p``. ``obs`` records
    ``n_cells`` and, with ``label_key``, ``purity`` (majority-label fraction):
    impure metacells mix states and create fake co-expression (Bilous et al.
    2024; mcRigor 2025).
    """
    from sklearn.cluster import KMeans

    _ensure_counts_layer(adata)
    counts = sp.csr_matrix(adata.layers["counts"])
    scratch = ad.AnnData(X=counts.astype(np.float32), var=adata.var[[]].copy())
    sc.pp.normalize_total(scratch, target_sum=target_sum)
    sc.pp.log1p(scratch)
    n_top = min(2000, scratch.n_vars)
    if n_top < scratch.n_vars:
        sc.pp.highly_variable_genes(scratch, n_top_genes=n_top, subset=True)
    n_comps = max(2, min(n_pcs, scratch.n_obs - 1, scratch.n_vars - 1))
    sc.pp.pca(scratch, n_comps=n_comps, random_state=seed)
    pcs = scratch.obsm["X_pca"]

    groups = (adata.obs[batch_key].astype(str).to_numpy() if batch_key
              else np.full(adata.n_obs, "all", dtype=object))
    codes = np.empty(adata.n_obs, dtype=np.int64)
    batch_of_code: list[str] = []
    for g in pd.unique(groups):
        rows = np.flatnonzero(groups == g)
        k = max(1, int(round(len(rows) / cells_per_metacell)))
        if k == 1:
            labels = np.zeros(len(rows), dtype=np.int64)
        else:
            labels = KMeans(n_clusters=k, n_init=3, random_state=seed).fit_predict(pcs[rows])
        codes[rows] = labels + len(batch_of_code)
        batch_of_code += [str(g)] * k
    n_meta = len(batch_of_code)
    summed = _group_sum(counts, codes, n_meta).tocsr()
    sizes = np.bincount(codes, minlength=n_meta)
    keep = sizes > 0
    obs = pd.DataFrame({"n_cells": sizes, "batch": batch_of_code},
                       index=[f"metacell{i:05d}" for i in range(n_meta)])
    if label_key is not None and label_key in adata.obs.columns:
        lab = pd.Series(adata.obs[label_key].astype(str).to_numpy())
        tab = pd.crosstab(codes, lab.to_numpy())
        obs["majority_label"] = tab.idxmax(axis=1).reindex(range(n_meta)).to_numpy()
        obs["purity"] = (tab.max(axis=1) / tab.sum(axis=1)).reindex(range(n_meta)).to_numpy()
    meta = ad.AnnData(X=summed[keep], obs=obs[keep], var=adata.var.copy())
    meta.layers["counts"] = meta.X.copy()
    if meta.n_obs < 100:
        warnings.warn(f"only {meta.n_obs} metacells: few samples for gene-gene statistics; "
                      "lower cells_per_metacell or use single cells", stacklevel=2)
    if "purity" in meta.obs and float(meta.obs["purity"].median()) < 0.9:
        warnings.warn(f"median metacell purity {meta.obs['purity'].median():.2f} < 0.90: "
                      "metacells mix cell types", stacklevel=2)
    sc.pp.normalize_total(meta, target_sum=target_sum)
    sc.pp.log1p(meta)
    return meta


REGISTRY: dict[str, StrategySpec] = {
    spec.name: spec
    for spec in [
        StrategySpec("cpm_log1p",
                     "Depth-normalize to 1e4 counts/cell, then log1p (safe default)",
                     cpm_log1p),
        StrategySpec("pearson_residuals",
                     "Analytic Pearson residuals; variance-stable across depths",
                     pearson_residuals),
        StrategySpec("pseudobulk",
                     "Sum counts per group (needs group_key); stabilizes sparse data",
                     pseudobulk),
        StrategySpec("raw_counts",
                     "No transformation; for methods that expect raw counts",
                     raw_counts),
        StrategySpec("cpm_log1p_combat",
                     "cpm_log1p then ComBat batch correction (needs batch_key); "
                     "never auto-recommended, removes batch-aligned covariance",
                     cpm_log1p_combat),
        StrategySpec("cpm_log1p_center_batch",
                     "cpm_log1p then per-batch gene centering (needs batch_key); removes "
                     "between-batch shifts only",
                     cpm_log1p_center_batch),
        StrategySpec("metacells",
                     "Sum counts of kNN-similar cells (k-means on PCA, per batch), then cpm_log1p; "
                     "reduces dropout",
                     metacells),
        StrategySpec("scran_pooled",
                     "scran pooled-size-factor normalization (R interop)",
                     unavailable="requires R scran through the optional R bridge (D-032) - planned"),
        StrategySpec("scvi_latent",
                     "scVI latent space as expression surrogate",
                     unavailable="requires scvi-tools - planned"),
    ]
}


def list_strategies() -> list[dict]:
    return [
        {"name": s.name, "description": s.description,
         "available": s.unavailable is None,
         "unavailable_reason": s.unavailable}
        for s in REGISTRY.values()
    ]


def get_strategy(name: str) -> StrategySpec:
    if name not in REGISTRY:
        raise StrategyError(
            f"Unknown strategy '{name}'. Available: {', '.join(REGISTRY)}"
        )
    spec = REGISTRY[name]
    if spec.unavailable is not None:
        raise StrategyError(f"Strategy '{name}' unavailable: {spec.unavailable}")
    return spec


def apply_strategy(name: str, adata: ad.AnnData, **params) -> tuple[ad.AnnData, dict]:
    """Apply a registered strategy; return (processed AnnData, provenance)."""
    spec = get_strategy(name)
    out = spec.apply(adata, **params)
    prov = {"strategy": name, "params": {k: str(v) for k, v in params.items()}}
    out.uns["grnbench_strategy"] = prov
    return out, prov
