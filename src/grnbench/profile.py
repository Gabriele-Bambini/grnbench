"""Dataset profile: every property the selector's rules and the interface use.

Computed from the counts matrix after ``run_qc`` annotated per-cell metrics
(``obs['_n_genes']``, ``obs['_total_counts']``, ``obs['_pct_mito']``) and
from the dataset context. Everything here is descriptive and deterministic;
the selector turns it into decisions.

Data-driven protocol hints complement what README files say:

* UMI counts are dominated by 1s among non-zero entries and rarely exceed a
  few thousand; full-length read counts have few 1s and large values;
* an unfiltered droplet matrix has tens of thousands of barcodes, most of
  them with almost no counts;
* batch/cell-type confounding is measured with Cramer's V, batch separation
  with the silhouette of batch labels on the first principal components.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field

import numpy as np
import pandas as pd
import scipy.sparse as sp

from grnbench.context import DatasetContext

PROFILE_PCA_CELLS = 3000   # silhouette/PCA on a fixed-seed subsample
MIN_CELLS_PER_CONTEXT = 200


@dataclass
class Profile:
    n_cells: int
    n_genes: int
    sparsity: float
    median_genes_per_cell: float
    median_counts_per_cell: float
    library_cv: float
    median_pct_mito: float
    n_mito_genes: int
    frac_ones_among_nonzero: float
    max_count: float
    data_technology_hint: str          # umi | full_length | unknown
    raw_droplets_suspected: bool
    organism: str | None
    technology: str | None             # context (README/metadata) first, data hint second
    platform: str | None
    time_course: bool
    id_verdict: str
    batch_key: str | None = None
    n_batches: int = 0
    min_cells_per_batch: int = 0
    batch_depth_ratio: float = 1.0      # max/min of per-batch median counts
    batch_celltype_cramers_v: float | None = None
    batch_silhouette: float | None = None
    cell_type_key: str | None = None
    n_cell_types: int = 0
    n_cell_types_with_enough_cells: int = 0
    min_cells_per_type: int = 0
    tf_requested: int = 0
    tf_detected: int = 0
    dependencies: dict = field(default_factory=dict)
    notes: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


def _nonzero_values(X, max_values: int = 2_000_000, seed: int = 0) -> np.ndarray:
    vals = X.data if sp.issparse(X) else np.asarray(X)[np.asarray(X) != 0]
    vals = np.asarray(vals).ravel()
    if vals.size > max_values:
        vals = np.random.default_rng(seed).choice(vals, max_values, replace=False)
    return vals


def technology_from_counts(X) -> tuple[str, float, float]:
    """``(hint, fraction of 1s among non-zeros, max value)`` from the matrix alone."""
    vals = _nonzero_values(X)
    if vals.size == 0:
        return "unknown", 0.0, 0.0
    ones = float(np.mean(vals == 1))
    vmax = float(vals.max())
    med = float(np.median(vals))
    # UMI data: small integers dominate (many 1s, median non-zero value of a few);
    # full-length read counts: amplified values, few 1s, large median and maximum.
    if (ones >= 0.3 or med <= 4) and vmax < 1e5:
        return "umi", ones, vmax
    if ones < 0.15 and med >= 10 and vmax >= 1e3:
        return "full_length", ones, vmax
    return "unknown", ones, vmax


def looks_unfiltered(totals) -> bool:
    """Tens of thousands of barcodes, most with almost no counts: an unfiltered droplet matrix."""
    totals = np.asarray(totals)
    return bool(totals.size > 20_000 and np.mean(totals < 100) > 0.5)


def cramers_v(a, b) -> float:
    table = pd.crosstab(pd.Series(a, dtype=str).to_numpy(), pd.Series(b, dtype=str).to_numpy()).to_numpy()
    n = table.sum()
    if n == 0 or min(table.shape) < 2:
        return 0.0
    expected = table.sum(1, keepdims=True) * table.sum(0, keepdims=True) / n
    chi2 = float(((table - expected) ** 2 / np.where(expected > 0, expected, 1)).sum())
    return float(np.sqrt(chi2 / (n * (min(table.shape) - 1))))


def batch_silhouette(adata, batch_key: str, n_cells: int = PROFILE_PCA_CELLS, seed: int = 0) -> float | None:
    """Silhouette of batch labels on 20 PCs of log-normalized HVGs (-1..1; ~0 = mixed)."""
    import scanpy as sc
    from sklearn.metrics import silhouette_score

    labels = adata.obs[batch_key].astype(str)
    if labels.nunique() < 2:
        return None
    idx = np.arange(adata.n_obs)
    if adata.n_obs > n_cells:
        idx = np.sort(np.random.default_rng(seed).choice(adata.n_obs, n_cells, replace=False))
    sub = adata[idx].copy()
    sub.X = sub.X.astype(np.float32)
    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)
    n_top = min(2000, sub.n_vars)
    sc.pp.highly_variable_genes(sub, n_top_genes=n_top, subset=True)
    n_comps = max(2, min(20, sub.n_obs - 1, sub.n_vars - 1))
    sc.pp.pca(sub, n_comps=n_comps, random_state=seed)
    lab = sub.obs[batch_key].astype(str).to_numpy()
    if len(set(lab)) < 2:
        return None
    return float(silhouette_score(sub.obsm["X_pca"], lab, random_state=seed))


def build_profile(adata, qc_metrics: dict, ctx: DatasetContext, *, id_verdict: dict | None = None,
                  tf_list: list[str] | None = None, compute_silhouette: bool = True,
                  raw_droplets: bool | None = None) -> Profile:
    from grnbench.kb import DECISIONS

    obs = adata.obs
    hint, ones, vmax = technology_from_counts(adata.X)
    totals = obs["_total_counts"].to_numpy()
    if raw_droplets is None:
        raw_droplets = looks_unfiltered(totals)
    tech = ctx.get("technology") or (hint if hint != "unknown" else None)
    notes = []
    platform = ctx.get("platform")
    full_length_platforms = {"Smart-seq2", "Smart-seq3", "Fluidigm C1"}
    if platform and tech and (platform in full_length_platforms) != (tech == "full_length"):
        notes.append(f"platform '{platform}' from the text contradicts the {tech} counts; ignored")
        platform = None
    if ctx.get("technology") and hint != "unknown" and ctx.get("technology") != hint:
        notes.append(f"text sources say {ctx.get('technology')} but count values look {hint}")

    p = Profile(
        n_cells=int(adata.n_obs), n_genes=int(adata.n_vars),
        sparsity=float(qc_metrics.get("sparsity", 0.0)),
        median_genes_per_cell=float(qc_metrics.get("median_genes_per_cell", 0.0)),
        median_counts_per_cell=float(qc_metrics.get("median_counts_per_cell", 0.0)),
        library_cv=float(qc_metrics.get("library_cv", 0.0)),
        median_pct_mito=float(qc_metrics.get("median_pct_mito", 0.0)),
        n_mito_genes=int(qc_metrics.get("n_mito_genes", 0)),
        frac_ones_among_nonzero=ones, max_count=vmax, data_technology_hint=hint,
        raw_droplets_suspected=raw_droplets, organism=ctx.get("organism"),
        technology=tech, platform=platform,
        time_course=ctx.get("time_course") == "yes",
        id_verdict=(id_verdict or {}).get("verdict", "symbols"), notes=notes,
    )

    batch = ctx.get("batch_key")
    if batch is not None and batch in obs.columns:
        levels = obs[batch].astype(str)
        counts = levels.value_counts()
        p.batch_key, p.n_batches = batch, int(len(counts))
        p.min_cells_per_batch = int(counts.min())
        med = pd.Series(totals, index=obs.index).groupby(levels.to_numpy()).median()
        p.batch_depth_ratio = float(med.max() / max(med.min(), 1e-9))
        if compute_silhouette and p.n_batches >= 2:
            try:
                p.batch_silhouette = batch_silhouette(adata, batch)
            except Exception as exc:  # noqa: BLE001 - descriptive only
                p.notes.append(f"batch silhouette not computed: {exc}")

    ctype = ctx.get("cell_type_key")
    if ctype is not None and ctype in obs.columns:
        counts = obs[ctype].astype(str).value_counts()
        p.cell_type_key, p.n_cell_types = ctype, int(len(counts))
        p.min_cells_per_type = int(counts.min())
        p.n_cell_types_with_enough_cells = int((counts >= MIN_CELLS_PER_CONTEXT).sum())
        if p.batch_key and p.batch_key != ctype:
            p.batch_celltype_cramers_v = cramers_v(obs[p.batch_key], obs[ctype])

    if tf_list:
        names = set(pd.Index(adata.var_names).astype(str).str.upper())
        p.tf_requested = len(tf_list)
        p.tf_detected = sum(str(t).upper() in names for t in tf_list)

    p.dependencies = {f"{d.id}.{o.id}": o.availability()[0] for d in DECISIONS for o in d.options}
    return p
