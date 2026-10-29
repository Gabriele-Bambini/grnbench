"""Deterministic correlation baseline: signed absolute Spearman correlations.

Rank-transform genes once, then a single Pearson correlation on the ranks
(equivalent to Spearman, vectorized, seedless-reproducible). Emits the full
regulator x target block with sign as ``direction``; |rho| as weight.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from grnbench.tfs import resolve_tf_names


def spearman_signed_abs(expr: pd.DataFrame,
                        tf_names: list[str] | None = None,
                        **_) -> pd.DataFrame:
    regulators = resolve_tf_names(expr.columns, tf_names)
    ranks = expr.rank(axis=0, method="average")
    corr = np.corrcoef(ranks.to_numpy(dtype=np.float64), rowvar=False)
    corr = np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)
    # corrcoef divides by the two std vectors one after the other, so r[i,j] and r[j,i] can differ in
    # the last bit (np.cov itself is symmetric); averaging makes A->B == B->A
    corr = (corr + corr.T) / 2

    reg_idx = [expr.columns.get_loc(g) for g in regulators]
    block = corr[np.ix_(reg_idx, list(range(expr.shape[1])))]
    src, tgt = np.meshgrid(regulators, expr.columns, indexing="ij")

    out = pd.DataFrame({
        "source": src.ravel(),
        "target": tgt.ravel(),
        "weight": np.abs(block).ravel(),
        "direction": np.sign(block).ravel().astype(int),
    })
    return out
