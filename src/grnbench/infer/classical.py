"""Classical GRN inference methods, in numpy/scipy/scikit-learn only.

Every adapter here is dependency-light on purpose. ``genie3`` and ``grnboost2``
need arboreto **and** a live dask cluster, which is exactly the dependency that
tends to be missing on a login node or in a reviewer's fresh environment — so a
benchmark whose only tree-based method is arboreto-backed cannot be run at all
there. ``random_forest`` computes the same quantity GENIE3 does (per-target
feature importance over the regulator set) with scikit-learn alone.

Methods and their sources:

``pearson_signed``  linear correlation baseline; the honest floor that any
                    method claiming to find regulation must clear.
``mutual_info``     plug-in MI on equal-frequency bins; catches monotone and
                    non-monotone dependence that correlation misses.
``clr``             Context Likelihood of Relatedness — Faith et al., *PLoS
                    Biology* 5(1):e8 (2007). MI z-scored against each gene's own
                    MI background, which is what suppresses the hub artifacts
                    raw MI is full of.
``ppcor``           partial correlation via the precision matrix, the quantity
                    BEELINE's PPCOR reports. scRNA-seq has far more genes than
                    cells, so the covariance is singular and must be shrunk;
                    this uses Ledoit-Wolf and says so in the output.
``lasso``           per-target L1 regression over the regulators (the sparse
                    linear-model family TIGRESS belongs to, without its
                    stability selection).
``random_forest``   per-target random-forest importance — GENIE3's estimator,
                    without arboreto or dask.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from grnbench.tfs import resolve_tf_names

try:  # scikit-learn is an optional extra ("infer")
    import sklearn  # noqa: F401
    SKLEARN_UNAVAILABLE: str | None = None
except ImportError:  # pragma: no cover - exercised by the availability probe
    SKLEARN_UNAVAILABLE = ("scikit-learn not installed - "
                           "pip install 'grnbench[infer]'")


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #

def _regulators(expr: pd.DataFrame, tf_names: list[str] | None) -> list[str]:
    return resolve_tf_names(expr.columns, tf_names)


def _long(regulators, targets, block: np.ndarray,
          direction: np.ndarray | None = None) -> pd.DataFrame:
    """Melt a regulators x targets score block into the canonical edge frame."""
    src, tgt = np.meshgrid(np.asarray(regulators, dtype=object),
                           np.asarray(targets, dtype=object), indexing="ij")
    out = pd.DataFrame({"source": src.ravel(), "target": tgt.ravel(),
                        "weight": np.asarray(block, dtype=float).ravel()})
    if direction is not None:
        out["direction"] = np.asarray(direction).ravel().astype(int)
    return out


def _equal_frequency_bins(X: np.ndarray, n_bins: int) -> np.ndarray:
    """Discretize each column into ``n_bins`` roughly equal-count bins.

    Quantile-style bins rather than equal-width: scRNA-seq expression is heavily
    zero-inflated and right-skewed. Equal values always receive the same code;
    this may leave bins empty, but avoids creating information from cell order
    when a large tied block (usually zeros) crosses a quantile boundary.
    """
    n_cells, n_genes = X.shape
    codes = np.empty((n_cells, n_genes), dtype=np.int64)
    for j in range(n_genes):
        col = X[:, j]
        _values, inverse, counts = np.unique(
            col, return_inverse=True, return_counts=True
        )
        starts = np.concatenate(([0], np.cumsum(counts[:-1], dtype=np.int64)))
        # Bin by the first rank occupied by each distinct value. This preserves
        # equal-frequency behavior for continuous data without splitting ties.
        value_codes = np.minimum((starts * n_bins) // n_cells, n_bins - 1)
        codes[:, j] = value_codes[inverse]
    return codes


def _mi_matrix(X: np.ndarray, rows: np.ndarray, n_bins: int) -> np.ndarray:
    """Plug-in mutual information between ``rows`` columns and every column.

    Returns ``(len(rows), n_genes)`` in nats. Computed by one-hot encoding each
    gene's bin assignment and reading joint counts off a matrix product, which
    keeps the whole thing in BLAS instead of a Python double loop.
    """
    n_cells, n_genes = X.shape
    codes = _equal_frequency_bins(X, n_bins)

    # one-hot: (n_cells, n_genes * n_bins), column (j * n_bins + b)
    onehot = np.zeros((n_cells, n_genes * n_bins), dtype=np.float32)
    flat = codes + np.arange(n_genes)[None, :] * n_bins
    np.put_along_axis(onehot, flat, 1.0, axis=1)

    marg = onehot.sum(axis=0).reshape(n_genes, n_bins) / n_cells       # p(b)
    with np.errstate(divide="ignore", invalid="ignore"):
        log_marg = np.where(marg > 0, np.log(marg), 0.0)
    entropy = -(marg * log_marg).sum(axis=1)                           # H(gene)

    out = np.empty((len(rows), n_genes), dtype=np.float64)
    for r, i in enumerate(rows):
        a = onehot[:, i * n_bins:(i + 1) * n_bins]                     # (cells, B)
        joint = (a.T @ onehot).reshape(n_bins, n_genes, n_bins) / n_cells
        joint = np.moveaxis(joint, 1, 0)                               # (genes,B,B)
        with np.errstate(divide="ignore", invalid="ignore"):
            log_j = np.where(joint > 0, np.log(joint), 0.0)
        h_joint = -(joint * log_j).sum(axis=(1, 2))                    # H(i, gene)
        out[r] = entropy[i] + entropy - h_joint                        # I = Hi+Hj-Hij
    return np.maximum(out, 0.0)   # plug-in MI is >=0; clip estimator noise


# --------------------------------------------------------------------------- #
# methods
# --------------------------------------------------------------------------- #

def pearson_signed(expr: pd.DataFrame, tf_names: list[str] | None = None,
                   **_) -> pd.DataFrame:
    regulators = _regulators(expr, tf_names)
    corr = np.corrcoef(expr.to_numpy(dtype=np.float64), rowvar=False)
    corr = np.nan_to_num(corr, nan=0.0)   # genes with zero variance
    idx = [expr.columns.get_loc(g) for g in regulators]
    block = corr[np.ix_(idx, range(expr.shape[1]))]
    return _long(regulators, expr.columns, np.abs(block), np.sign(block))


def mutual_info(expr: pd.DataFrame, tf_names: list[str] | None = None,
                n_bins: int = 8, **_) -> pd.DataFrame:
    regulators = _regulators(expr, tf_names)
    X = expr.to_numpy(dtype=np.float64)
    rows = np.array([expr.columns.get_loc(g) for g in regulators])
    return _long(regulators, expr.columns, _mi_matrix(X, rows, int(n_bins)))


def _clr_scores(X: np.ndarray, n_bins: int) -> np.ndarray:
    """Full gene x gene CLR score matrix (Faith et al. 2007) from plug-in MI.

    Needs the *full* gene x gene MI matrix, not just the regulator rows: the
    score of edge (i, j) is defined against the MI background of both i and j.
    """
    n_genes = X.shape[1]
    mi = _mi_matrix(X, np.arange(n_genes), int(n_bins))
    np.fill_diagonal(mi, 0.0)

    mean = mi.mean(axis=1, keepdims=True)
    sd = mi.std(axis=1, keepdims=True)
    sd = np.where(sd > 0, sd, 1.0)
    z = np.maximum((mi - mean) / sd, 0.0)                # one-sided z per gene
    score = np.sqrt(z ** 2 + (z.T) ** 2)                 # symmetric combination
    return np.nan_to_num(score, nan=0.0, posinf=0.0, neginf=0.0)


def clr(expr: pd.DataFrame, tf_names: list[str] | None = None,
        n_bins: int = 8, **_) -> pd.DataFrame:
    """Context Likelihood of Relatedness (Faith et al. 2007)."""
    regulators = _regulators(expr, tf_names)
    score = _clr_scores(expr.to_numpy(dtype=np.float64), n_bins)
    idx = [expr.columns.get_loc(g) for g in regulators]
    return _long(regulators, expr.columns, score[np.ix_(idx, range(expr.shape[1]))])


def ppcor(expr: pd.DataFrame, tf_names: list[str] | None = None,
          **_) -> pd.DataFrame:
    """Partial correlation from a shrunk precision matrix.

    With more genes than cells the sample covariance is singular and a plain
    inverse is meaningless, so this uses Ledoit-Wolf shrinkage. That is a real
    departure from R's ``ppcor``, which inverts directly and simply fails (or
    returns garbage) in the p > n regime every scRNA-seq dataset is in.
    """
    from sklearn.covariance import LedoitWolf

    regulators = _regulators(expr, tf_names)
    X = expr.to_numpy(dtype=np.float64)
    keep = X.std(axis=0) > 0
    prec = np.zeros((X.shape[1], X.shape[1]), dtype=np.float64)
    sub = LedoitWolf(assume_centered=False).fit(X[:, keep]).precision_
    prec[np.ix_(keep, keep)] = sub

    d = np.sqrt(np.abs(np.diag(prec)))
    d = np.where(d > 0, d, 1.0)
    partial = -prec / np.outer(d, d)
    np.fill_diagonal(partial, 0.0)

    idx = [expr.columns.get_loc(g) for g in regulators]
    block = partial[np.ix_(idx, range(X.shape[1]))]
    return _long(regulators, expr.columns, np.abs(block), np.sign(block))


def _per_target_importance(expr, tf_names, fit_one, **kw) -> pd.DataFrame:
    """Shared loop: regress every target on the regulators, minus itself.

    The self-column is dropped before fitting rather than after: leaving a gene
    in its own predictor set lets it explain itself perfectly, and the resulting
    importances for the real regulators are then all near zero.
    """
    regulators = _regulators(expr, tf_names)
    X_all = expr.to_numpy(dtype=np.float64)
    reg_idx = np.array([expr.columns.get_loc(g) for g in regulators])
    block = np.zeros((len(regulators), expr.shape[1]), dtype=np.float64)

    for t in range(expr.shape[1]):
        mask = reg_idx != t
        if not mask.any():
            continue
        scores = fit_one(X_all[:, reg_idx[mask]], X_all[:, t], **kw)
        block[mask, t] = scores
    return _long(regulators, expr.columns, block)


def _fit_rf(Xr, y, n_estimators=100, max_features="sqrt", seed=42):
    from sklearn.ensemble import RandomForestRegressor

    if y.std() == 0:
        return np.zeros(Xr.shape[1])
    m = RandomForestRegressor(n_estimators=int(n_estimators),
                              max_features=max_features,
                              random_state=int(seed), n_jobs=-1)
    m.fit(Xr, y)
    return m.feature_importances_


def _fit_lasso(Xr, y, alpha=None, seed=42):
    from sklearn.linear_model import LassoCV, lasso_path

    if y.std() == 0:
        return np.zeros(Xr.shape[1])
    if alpha is None:
        # a single shared path is far cheaper than LassoCV per target; take the
        # coefficient magnitude at the knee of the path
        _, coefs, _ = lasso_path(Xr, y, n_alphas=20)
        return np.abs(coefs[:, len(coefs[0]) // 2])
    from sklearn.linear_model import Lasso
    return np.abs(Lasso(alpha=float(alpha), random_state=int(seed)).fit(Xr, y).coef_)


def random_forest(expr: pd.DataFrame, tf_names: list[str] | None = None,
                  n_estimators: int = 100, seed: int = 42, **_) -> pd.DataFrame:
    return _per_target_importance(expr, tf_names, _fit_rf,
                                  n_estimators=n_estimators, seed=seed)


def lasso(expr: pd.DataFrame, tf_names: list[str] | None = None,
          alpha: float | None = None, seed: int = 42, **_) -> pd.DataFrame:
    return _per_target_importance(expr, tf_names, _fit_lasso,
                                  alpha=alpha, seed=seed)
