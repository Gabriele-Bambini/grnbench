"""Loud-failure QC gates for GRN-bound scRNA-seq datasets.

Philosophy: QC *reports* first (on the unfiltered matrix), then gates decide
whether preprocessing may proceed. A gate ``fail`` aborts the pipeline unless
the caller forces it; every gate carries its numbers so reports are auditable.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
from anndata import AnnData

from grnbench.ingest import MITO_PREFIXES, RIBO_PREFIXES

STATUS_PASS, STATUS_WARN, STATUS_FAIL = "pass", "warn", "fail"


@dataclass
class Thresholds:
    min_genes_per_cell: int = 200
    min_cells_per_gene: int = 3
    mito_warn: float = 10.0
    mito_fail: float = 25.0
    low_genes_cell_frac_warn: float = 0.10
    low_genes_cell_frac_fail: float = 0.50
    library_cv_warn: float = 2.0
    sparsity_warn: float = 0.95
    sparsity_fail: float = 0.995
    min_cells_warn: int = 100
    min_genes_warn: int = 3000
    zero_cells_warn: float = 0.005
    zero_cells_fail: float = 0.05
    imbalance_ratio_warn: float = 50.0
    min_class_cells_warn: int = 30
    doublet_pct_warn: float = 5.0
    doublet_rate_max: float = 0.20   # --drop-doublets refuses above this call rate
    # Per-cell upper filters; None = off. Set by the user or by grnbench.autoprep
    # (MAD-based outlier bounds) - never silently.
    max_genes_per_cell: int | None = None
    max_counts_per_cell: float | None = None
    max_pct_mito_cell: float | None = None


@dataclass
class Gate:
    name: str
    status: str
    detail: str
    value: float | int | str | None = None

    @property
    def failed(self) -> bool:
        return self.status == STATUS_FAIL

    def to_dict(self) -> dict:
        return {"name": self.name, "status": self.status,
                "detail": self.detail, "value": self.value}


@dataclass
class QCReport:
    metrics: dict = field(default_factory=dict)
    gates: list[Gate] = field(default_factory=list)

    @property
    def has_failures(self) -> bool:
        return any(g.failed for g in self.gates)

    def to_dict(self) -> dict:
        return {"metrics": self.metrics,
                "gates": [g.to_dict() for g in self.gates]}


def _colsums(X) -> np.ndarray:
    return np.asarray(X.sum(axis=0)).ravel()


def run_qc(adata: AnnData, thresholds: Thresholds | None = None, *,
           group_key: str | None = None,
           label_key: str | None = None,
           batch_key: str | None = None,
           id_verdict: dict | None = None,
           n_dup_barcodes: int = 0) -> QCReport:
    """Compute per-cell QC metrics (annotated into ``adata.obs``) and gates.

    ``label_key`` (class labels) and ``group_key`` (pseudobulk / recommendation
    grouping) are audited for imbalance separately: ``class_imbalance`` covers
    the labels - or the group key when that is the only column given - and
    ``group_imbalance`` covers the group key when both are present. Passing
    both used to audit only the group key, so a 500/90/6 label split came back
    "balanced" because the donors were.
    """
    t = thresholds or Thresholds()
    X = adata.X
    sparse = sp.issparse(X)

    total_counts = np.asarray(X.sum(axis=1)).ravel()
    n_genes_cell = X.getnnz(axis=1) if sparse else (X > 0).sum(axis=1)
    n_cells_gene = X.getnnz(axis=0) if sparse else (X > 0).sum(axis=0)

    names_upper = pd_index_upper(adata.var_names)
    mito_mask = np.zeros(len(names_upper), dtype=bool)
    ribo_mask = np.zeros(len(names_upper), dtype=bool)
    for pref in MITO_PREFIXES:
        mito_mask |= names_upper.str.startswith(pref)
    for pref in RIBO_PREFIXES:
        ribo_mask |= names_upper.str.startswith(pref)
    mito_counts = np.asarray(X[:, mito_mask].sum(axis=1)).ravel()
    ribo_counts = np.asarray(X[:, ribo_mask].sum(axis=1)).ravel()
    safe = np.where(total_counts > 0, total_counts, 1)

    pct_mito = 100.0 * mito_counts / safe
    pct_ribo = 100.0 * ribo_counts / safe
    sparsity = 1.0 - (X.getnnz() if sparse else np.count_nonzero(X)) / (X.shape[0] * X.shape[1])
    lib_cv = float(np.std(total_counts) / np.mean(total_counts)) if np.mean(total_counts) else np.inf

    adata.obs["_total_counts"] = total_counts
    adata.obs["_n_genes"] = n_genes_cell
    adata.obs["_pct_mito"] = pct_mito
    adata.obs["_pct_ribo"] = pct_ribo

    metrics = {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "median_genes_per_cell": float(np.median(n_genes_cell)),
        "median_counts_per_cell": float(np.median(total_counts)),
        "library_cv": lib_cv,
        "sparsity": float(sparsity),
        "median_pct_mito": float(np.median(pct_mito)),
        "median_pct_ribo": float(np.median(pct_ribo)),
        "n_mito_genes": int(mito_mask.sum()),
        "n_ribo_genes": int(ribo_mask.sum()),
        "frac_cells_below_min_genes": float(np.mean(n_genes_cell < t.min_genes_per_cell)),
        "frac_genes_below_min_cells": float(np.mean(n_cells_gene < t.min_cells_per_gene)),
    }

    gates: list[Gate] = [
        _threshold_gate(
            "cells_below_min_genes", metrics["frac_cells_below_min_genes"],
            t.low_genes_cell_frac_warn, t.low_genes_cell_frac_fail,
            f"fraction of cells with < {t.min_genes_per_cell} detected genes",
            higher_is_worse=True, fmt="{:.1%}",
        ),
        _threshold_gate(
            "mito_fraction", metrics["median_pct_mito"],
            t.mito_warn, t.mito_fail, "median % mitochondrial counts per cell",
            higher_is_worse=True, fmt="{:.1f}%",
        ),
        _threshold_gate(
            "library_size_cv", lib_cv, t.library_cv_warn, None,
            "CV of total counts per cell (uneven depth: consider "
            "pearson_residuals / pseudobulk)",
            higher_is_worse=True, fmt="{:.2f}",
        ),
        _threshold_gate(
            "sparsity", sparsity, t.sparsity_warn, t.sparsity_fail,
            "matrix zero fraction (high dropout hurts correlation/MI methods)",
            higher_is_worse=True, fmt="{:.1%}",
        ),
        Gate(
            "small_dataset", STATUS_WARN if adata.n_obs < t.min_cells_warn else STATUS_PASS,
            f"n_cells={adata.n_obs}" + ("" if adata.n_obs >= t.min_cells_warn
                                        else f" (< {t.min_cells_warn}: GRN inference will be unstable)"),
            adata.n_obs,
        ),
        Gate(
            "small_gene_space", STATUS_WARN if adata.n_vars < t.min_genes_warn else STATUS_PASS,
            f"n_genes={adata.n_vars}" + ("" if adata.n_vars >= t.min_genes_warn
                                         else f" (< {t.min_genes_warn})"),
            adata.n_vars,
        ),
    ]

    # ---- hygiene gates (appended; original six keep their order) ---- #
    n_zero = int((total_counts == 0).sum())
    zero_frac = n_zero / max(adata.n_obs, 1)
    gates.append(Gate(
        "zero_count_cells",
        STATUS_FAIL if zero_frac > t.zero_cells_fail else (
            STATUS_WARN if zero_frac > t.zero_cells_warn else STATUS_PASS),
        f"{n_zero} cells with zero total counts ({zero_frac:.2%}) - dropped "
        "before normalization to avoid divide-by-zero",
        n_zero,
    ))

    audits: list[tuple[str, str]] = []
    if label_key is not None and label_key in adata.obs.columns:
        audits.append(("class_imbalance", label_key))
    if group_key is not None and group_key in adata.obs.columns and group_key != label_key:
        audits.append(("group_imbalance" if audits else "class_imbalance", group_key))
    for gate_name, col in audits:
        gates.append(_imbalance_gate(gate_name, col, adata.obs[col], t))

    if batch_key is None:
        gates.append(Gate(
            "batch_annotation", STATUS_WARN,
            "no batch/donor annotation found: cross-batch GRN edges may be "
            "confounded; provide one if it exists", "",
        ))
    else:
        gates.append(Gate("batch_annotation", STATUS_PASS,
                          f"batch column '{batch_key}' present", batch_key))

    if id_verdict is not None:
        v = id_verdict.get("verdict")
        gates.append(Gate(
            "gene_id_type",
            STATUS_WARN if v in {"ensembl_only", "mixed"} else STATUS_PASS,
            ("gene IDs are Ensembl-only or mixed: TF-list matching will likely "
             "fail - map symbols first (mygene/gget)") if v == "ensembl_only"
            else ("gene IDs are MIXED symbols/Ensembl: harmonize before "
                  "TF matching") if v == "mixed"
            else f"gene IDs look like symbols (ensembl frac "
                 f"{id_verdict.get('ensembl_frac', 0):.1%})",
            v,
        ))

    gates.append(Gate(
        "duplicate_barcodes", STATUS_WARN if n_dup_barcodes > 0 else STATUS_PASS,
        (f"{n_dup_barcodes} duplicated obs_names, made unique with numeric suffixes "
         "(originals in obs['barcode_as_loaded']); check for pooled libraries")
        if n_dup_barcodes else "obs_names unique", n_dup_barcodes,
    ))

    from grnbench.hygiene import doublet_heuristic_fraction
    dbl_frac, _flagged = doublet_heuristic_fraction(adata)
    dbl_pct = dbl_frac * 100.0
    gates.append(Gate(
        "doublet_heuristic",
        STATUS_WARN if dbl_pct > t.doublet_pct_warn else STATUS_PASS,
        (f"{dbl_pct:.1f}% cells jointly top-1% in genes+counts "
         f"(> {t.doublet_pct_warn:.0f}% suspicious - consider scrublet)")
        if dbl_pct > t.doublet_pct_warn
        else f"{dbl_pct:.1f}% joint-extreme cells (within tolerance)",
        round(dbl_pct, 2),
    ))
    return QCReport(metrics=metrics, gates=gates)


def _imbalance_gate(name: str, col: str, series, t: Thresholds) -> Gate:
    from grnbench.hygiene import missing_label_mask
    labelled = series[~missing_label_mask(series)].astype(str)
    vc = labelled.value_counts()
    if len(vc) < 2:
        return Gate(name, STATUS_WARN, f"'{col}': only one non-missing class", len(vc))
    ratio = float(vc.max() / vc.min()) if vc.min() > 0 else float("inf")
    smallest = int(vc.min())
    issues = []
    if ratio > t.imbalance_ratio_warn:
        issues.append(f"max/min ratio {ratio:.0f} > {t.imbalance_ratio_warn:.0f}")
    if smallest < t.min_class_cells_warn:
        issues.append(f"smallest class has {smallest} cells < {t.min_class_cells_warn}")
    return Gate(
        name, STATUS_WARN if issues else STATUS_PASS,
        (f"'{col}': {len(vc)} classes; " + "; ".join(issues)) if issues
        else f"'{col}': {len(vc)} classes, balanced enough (max/min {ratio:.1f})",
        round(ratio, 2),
    )


def _threshold_gate(name, value, warn_at, fail_at, detail, *, higher_is_worse, fmt):
    if higher_is_worse:
        status = STATUS_FAIL if (fail_at is not None and value > fail_at) else (
            STATUS_WARN if value > warn_at else STATUS_PASS)
    else:
        status = STATUS_FAIL if (fail_at is not None and value < fail_at) else (
            STATUS_WARN if value < warn_at else STATUS_PASS)
    return Gate(name, status, f"{detail}: {fmt.format(value)}", value=round(float(value), 6))


def pd_index_upper(idx) -> "pd.Index":
    """Uppercase a var-name index as a pandas Index (keeps .str accessors)."""
    import pandas as pd
    return pd.Index(idx).astype(str).str.upper()
