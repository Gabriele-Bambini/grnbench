"""Feature selection with the TF-union guarantee.

Rule: HVG selection may never silently drop transcription factors. The kept
set is ``HVGs UNION (TFs detected in data)``, and the report states exactly how
many TFs were rescued by the union.
"""
from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from grnbench.ingest import looks_like_counts
from grnbench.tfs import match_tfs_in_data


def select_features(adata: ad.AnnData, n_hvg: int = 3000,
                    tf_list: list[str] | None = None) -> tuple[ad.AnnData, dict]:
    """Select HVGs (union TFs) on a log-normalized scratch copy.

    Works whether ``adata.X`` is raw counts or already normalized: if counts
    are detected, normalization is computed on a scratch copy so the input is
    never mutated. Returns the subset AnnData and a report dict.
    """
    matched, requested = (match_tfs_in_data(tf_list, adata.var_names)
                          if tf_list else ([], 0))

    if n_hvg >= adata.n_vars:
        report = {
            "n_hvg_requested": n_hvg,
            "n_hvg_selected": int(adata.n_vars),
            "tfs_requested": requested,
            "tfs_matched": len(matched),
            "tfs_added_by_union": 0,
            "n_features_kept": int(adata.n_vars),
            "note": (f"requested {n_hvg:,} HVGs but the matrix has {adata.n_vars:,} genes; "
                     f"kept all {adata.n_vars:,} genes"),
        }
        return adata.copy(), report

    scratch = adata.copy()
    ok, why = looks_like_counts(scratch.X)
    # X is rebuilt from counts below, so a log1p marker left by an earlier step is stale: drop it on the
    # copy, or scanpy warns "already log-transformed" about a matrix that is not
    if ok:
        scratch.uns.pop("log1p", None)
        sc.pp.normalize_total(scratch, target_sum=1e4)
        sc.pp.log1p(scratch)
    elif "counts" in scratch.layers:  # normalized X, but counts available
        scratch.uns.pop("log1p", None)
        scratch.X = scratch.layers["counts"].copy()
        sc.pp.normalize_total(scratch, target_sum=1e4)
        sc.pp.log1p(scratch)
    else:  # normalized X without counts: HVG on X as-is (seurat flavor expects log)
        pass

    n_top = min(int(n_hvg), scratch.n_vars)
    sc.pp.highly_variable_genes(scratch, n_top_genes=n_top, flavor="seurat")
    hvg_mask = scratch.var["highly_variable"].to_numpy(dtype=bool)

    var_upper = pd.Index(adata.var_names).astype(str).str.upper()
    # Index.isin returns a numpy bool array directly
    tf_mask = (var_upper.isin(set(matched)) if matched
               else np.zeros(adata.n_vars, dtype=bool))

    keep = hvg_mask | tf_mask
    report = {
        "n_hvg_requested": int(n_top),
        "n_hvg_selected": int(hvg_mask.sum()),
        "tfs_requested": requested,
        "tfs_matched": len(matched),
        "tfs_added_by_union": int(np.sum(tf_mask & ~hvg_mask)),
        "n_features_kept": int(keep.sum()),
    }
    return adata[:, keep].copy(), report
