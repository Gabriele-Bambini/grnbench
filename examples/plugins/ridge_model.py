"""Example Python plugin: multivariate ridge regression of each gene on the TFs.

A real, usable method (not a stub) that shows the whole contract in ~40 lines:
one class, one ``infer``, canonical edge columns out. Unlike GENIE3/GRNBoost2 it
is linear and deterministic, so it makes a decent cheap third opinion.

Run it without installing anything::

    grnbench-infer processed.h5ad --out runs/infer --tf-list tfs.txt \\
        --plugin examples/plugins/ridge_model.py:RidgeGRN

or via the spec file::

    grnbench-infer processed.h5ad --out runs/infer \\
        --plugin examples/plugins/ridge.json

To ship it as a proper package, declare it in your own ``pyproject.toml``::

    [project.entry-points."grnbench.methods"]
    ridge = "my_pkg.ridge_model:RidgeGRN"

and grnbench discovers it automatically — no flags at all.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from grnbench.tfs import resolve_tf_names


class RidgeGRN:
    """Per-target ridge regression; |coefficient| is the edge weight."""

    name = "ridge"
    description = "Ridge regression of each gene on the TF panel (scikit-learn)"

    def __init__(self, alpha: float = 1.0):
        self.alpha = alpha

    def infer(self, adata, tf_list=None, *, alpha: float | None = None,
              **params) -> pd.DataFrame:
        # `adata` arrives already normalized by grnbench (see prepare_adata).
        from sklearn.linear_model import Ridge

        genes = [str(g) for g in adata.var_names]
        X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
        X = np.asarray(X, dtype=np.float64)

        regulators = resolve_tf_names(genes, tf_list)
        reg_idx = [genes.index(g) for g in regulators]
        R = X[:, reg_idx]
        # standardize regulators so coefficients are comparable across genes
        R = (R - R.mean(0)) / np.where(R.std(0) > 0, R.std(0), 1.0)

        ridge_alpha = float(alpha if alpha is not None else self.alpha)
        coef = np.zeros((len(genes), len(regulators)), dtype=np.float64)
        for target_index in range(len(genes)):
            keep = np.asarray(reg_idx) != target_index
            if not keep.any():
                continue
            fitted = Ridge(alpha=ridge_alpha).fit(R[:, keep], X[:, target_index])
            coef[target_index, keep] = fitted.coef_

        src, tgt = np.meshgrid(regulators, genes, indexing="ij")
        return pd.DataFrame({
            "source": src.ravel(),
            "target": tgt.ravel(),
            "weight": np.abs(coef.T).ravel(),
            "direction": np.sign(coef.T).ravel().astype(int),
        })
