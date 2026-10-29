"""Regression-based methods from the GRN literature.

Same rules as :mod:`grnbench.infer.association`: native reimplementations from the
published descriptions, deterministic given ``seed``, scikit-learn and scipy only,
and not the original packages (each docstring says where they differ).

``tigress``            LARS stability selection with random regulator weights.
``bbsr``               CLR preselection + BIC best-subset regression (Inferelator style).
``pcnet``              principal-component regression per gene (scTenifoldNet).
``extra_trees``        extremely randomised trees importance (GENIE3's ET variant).
``gradient_boosting``  boosted-tree importance (GRNBoost2's idea, without dask).
"""
from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd

from grnbench.infer.association import _error, _tf_index
from grnbench.infer.classical import (
    SKLEARN_UNAVAILABLE,
    _clr_scores,
    _long,
    _per_target_importance,
)


def _standardize(X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Columns centred and scaled to unit variance; zero-variance columns become 0."""
    centred = X - X.mean(axis=0)
    sd = centred.std(axis=0)
    z = np.divide(centred, sd, out=np.zeros_like(centred), where=sd > 0)
    return z, sd > 0


# --------------------------------------------------------------------------- #
# TIGRESS
# --------------------------------------------------------------------------- #

def tigress(expr: pd.DataFrame, tf_names: list[str] | None = None,
            n_resamples: int = 50, alpha: float = 0.2, n_steps: int = 5,
            subsample: float = 0.5, seed: int = 42, **_) -> pd.DataFrame:
    """TIGRESS-style stability selection (Haury, Mordelet, Vera-Licona & Vert 2012,
    *BMC Syst Biol* 6:145).

    For every target, ``n_resamples`` times: draw ``subsample`` of the cells, multiply
    each regulator by a random weight in ``[alpha, 1]`` (the randomised lasso), run
    ``n_steps`` LARS steps, and record which regulators have a non-zero coefficient after
    each step ``L``. The score is the "area" of the selection-frequency curve: the mean
    over ``L = 1..n_steps`` of the fraction of resamples in which the regulator is selected
    at step ``L``.

    Differences from the R package: ``n_resamples`` defaults to 50 for CPU time; the
    cell subsets are shared by all targets of a resample (the random regulator weights
    are drawn per target); columns are standardised once over all cells. Regulators
    never selected get 0 and tie at the bottom.
    """
    if SKLEARN_UNAVAILABLE is not None:
        raise _error(f"unavailable: {SKLEARN_UNAVAILABLE}")
    from sklearn.linear_model import lars_path_gram

    regulators, rows = _tf_index(expr, tf_names)
    X = expr.to_numpy(dtype=np.float64)
    n_cells, n_genes = X.shape
    z, varying = _standardize(X)
    n_tf = len(rows)
    tf_pos = {int(g): i for i, g in enumerate(rows)}
    n_sub = min(max(int(round(float(subsample) * n_cells)), 3), n_cells)
    n_steps = int(n_steps)
    rng = np.random.default_rng(int(seed))
    score = np.zeros((n_tf, n_genes))

    for _ in range(int(n_resamples)):
        part = z[np.sort(rng.choice(n_cells, n_sub, replace=False))]
        part = part - part.mean(axis=0)
        reg = part[:, rows]
        gram = reg.T @ reg
        xy_all = reg.T @ part
        weights = rng.uniform(float(alpha), 1.0, size=(n_genes, n_tf))
        for g in range(n_genes):
            if not varying[g]:
                continue
            w = weights[g]
            usable = varying[rows].copy()
            if g in tf_pos:
                usable[tf_pos[g]] = False
            idx = np.flatnonzero(usable)
            if idx.size == 0:
                continue
            wu = w[idx]
            _, _, path = lars_path_gram(xy_all[idx, g] * wu, gram[np.ix_(idx, idx)] * np.outer(wu, wu),
                                        n_samples=n_sub, max_iter=n_steps, method="lasso")
            last = path.shape[1] - 1                             # column 0 is the empty model
            for step in range(1, n_steps + 1):
                score[idx[path[:, min(step, last)] != 0.0], g] += 1.0 / n_steps
    return _long(regulators, expr.columns, score / max(int(n_resamples), 1))


# --------------------------------------------------------------------------- #
# BBSR
# --------------------------------------------------------------------------- #

MAX_PRESELECT = 16


def _subset_table(k: int) -> list[np.ndarray]:
    return [np.array(list(combinations(range(k), s)), dtype=int) for s in range(1, k + 1)]


def _best_subset(corr_pp: np.ndarray, corr_py: np.ndarray, n: int,
                 table: list[np.ndarray]) -> tuple[np.ndarray, float]:
    """BIC-best subset of standardised predictors, from correlations only.

    With unit-variance columns the residual variance of subset ``S`` is
    ``1 - r_S' R_SS^-1 r_S`` (``R`` predictor correlations, ``r`` predictor-target
    correlations). Returns (indices of the best subset, its residual variance); the
    empty model (residual variance 1, BIC 0) wins when nothing beats it.
    """
    best_bic, best, best_rss = 0.0, np.array([], dtype=int), 1.0
    for subsets in table:
        s = subsets.shape[1]
        sub = corr_pp[subsets[:, :, None], subsets[:, None, :]] + 1e-9 * np.eye(s)
        r = corr_py[subsets]
        beta = np.linalg.solve(sub, r[..., None])[..., 0]
        rss = np.maximum(1.0 - (beta * r).sum(axis=1), 1e-12)
        bic = n * np.log(rss) + s * np.log(n)
        i = int(np.argmin(bic))
        if bic[i] < best_bic:
            best_bic, best, best_rss = float(bic[i]), subsets[i], float(rss[i])
    return best, best_rss


def bbsr(expr: pd.DataFrame, tf_names: list[str] | None = None,
         n_preselect: int = 10, n_bins: int = 8, **_) -> pd.DataFrame:
    """Bayesian-best-subset-style regression, after the Inferelator (Greenfield et al.
    2013, *Bioinformatics* 29:1060; Skok Gibbs et al. 2022, *Bioinformatics* 38:2519).

    For each target the ``n_preselect`` regulators with the highest CLR score are kept,
    every subset of them is fitted by least squares, and the subset with the lowest BIC
    is chosen (the empty model if nothing beats it). The weight of a chosen regulator
    is the fraction of the residual variance it removes given the others,
    ``1 - rss_full / rss_without_it``; unchosen regulators get 0 and tie at the bottom.

    Reimplemented from that description, not the Inferelator code: no prior network,
    no multitask (AMuSR) mode, and the BIC is the plain Gaussian one. Choosing the best
    ``n_preselect`` of many regulators and then testing them with a BIC that ignores that
    choice overstates inclusion on noise: with 100 regulators and 1000 cells about one pure
    noise target in five still gets a non-empty model, more with fewer cells. It refuses
    ``n_preselect`` above 16 (the subsets are enumerated) and fewer than 2 x ``n_preselect``
    cells.
    """
    regulators, rows = _tf_index(expr, tf_names)
    X = expr.to_numpy(dtype=np.float64)
    n_cells, n_genes = X.shape
    if int(n_preselect) > MAX_PRESELECT:
        raise _error(f"bbsr enumerates every subset of the preselected regulators: "
                     f"n_preselect={int(n_preselect)} exceeds {MAX_PRESELECT}")
    if n_cells < 2 * int(n_preselect):
        raise _error(f"bbsr needs at least 2 x n_preselect cells ({2 * int(n_preselect)}) for the BIC "
                     f"to mean anything; got {n_cells}")
    z, varying = _standardize(X)
    n_tf = len(rows)
    k = min(int(n_preselect), n_tf)
    clr = _clr_scores(X, int(n_bins))[rows]                      # TF x gene
    corr = (z[:, rows].T @ z) / n_cells                          # TF x gene
    corr_tf = corr[:, rows]
    table = _subset_table(k)
    tf_pos = {int(g): i for i, g in enumerate(rows)}
    block = np.zeros((n_tf, n_genes))
    for g in range(n_genes):
        if not varying[g] or k == 0:
            continue
        score = np.where(varying[rows], clr[:, g], -np.inf)
        if g in tf_pos:
            score[tf_pos[g]] = -np.inf
        top = np.argsort(-score, kind="stable")[:k]
        top = top[np.isfinite(score[top])]
        if top.size == 0:
            continue
        sub_table = table if top.size == k else _subset_table(top.size)
        chosen, rss_full = _best_subset(corr_tf[np.ix_(top, top)], corr[top, g], n_cells, sub_table)
        for pos in chosen:
            rest = np.array([c for c in chosen if c != pos], dtype=int)
            if rest.size:
                r = corr[top[rest], g]
                rss_rest = 1.0 - float(r @ np.linalg.solve(
                    corr_tf[np.ix_(top[rest], top[rest])] + 1e-9 * np.eye(rest.size), r))
            else:
                rss_rest = 1.0
            block[top[pos], g] = max(1.0 - rss_full / max(rss_rest, 1e-12), 0.0)
    return _long(regulators, expr.columns, block)


# --------------------------------------------------------------------------- #
# pcNet
# --------------------------------------------------------------------------- #

def pcnet(expr: pd.DataFrame, tf_names: list[str] | None = None,
          n_comp: int = 3, max_genes: int = 4000, **_) -> pd.DataFrame:
    """pcNet: principal-component regression of every gene on all the others.

    The network-construction step of scTenifoldNet (Osorio et al. 2020, *Patterns* 1(9)).
    For each gene ``j`` it is regressed on the top ``n_comp`` principal components of the
    other genes, and the coefficients are mapped back to genes,
    ``beta_j = sum_c v_c (v_c' r_j) / lambda_c`` in terms of the correlation matrix
    ``R`` (``v_c, lambda_c`` its leave-``j``-out eigenpairs, ``r_j`` the correlations
    of ``j`` with the rest). The weight of ``TF i -> j`` is ``|beta_j[i]|``.

    Only that step is implemented: not the tensor decomposition or the comparison of
    two networks that scTenifoldNet adds. One partial eigendecomposition of a genes x genes
    matrix per gene: the cost grows with the cube of the number of genes (about a minute at
    1000 genes, extrapolated ten at 2000), so it refuses more than ``max_genes`` genes; raise it
    knowingly.
    """
    from scipy.sparse.linalg import ArpackNoConvergence, LinearOperator, eigsh

    if expr.shape[1] > int(max_genes):
        raise _error(f"pcnet builds a genes x genes matrix: {expr.shape[1]} genes exceed "
                     f"max_genes={int(max_genes)}; keep fewer genes or raise max_genes")
    regulators, rows = _tf_index(expr, tf_names)
    X = expr.to_numpy(dtype=np.float64)
    n_cells, n_genes = X.shape
    z, varying = _standardize(X)
    corr = (z.T @ z) / n_cells
    n_comp = max(1, min(int(n_comp), n_genes - 2))
    tf_pos = {int(g): i for i, g in enumerate(rows)}
    block = np.zeros((len(rows), n_genes))
    sign = np.zeros_like(block)
    for j in range(n_genes):
        if not varying[j]:
            continue
        others = np.flatnonzero(varying & (np.arange(n_genes) != j))
        if others.size <= n_comp + 1:
            continue
        sub_n = others.size
        if sub_n <= 400:
            values, vectors = np.linalg.eigh(corr[np.ix_(others, others)])
            values, vectors = values[-n_comp:], vectors[:, -n_comp:]
        else:
            op = LinearOperator((sub_n, sub_n), dtype=np.float64,
                                matvec=lambda v, others=others: _masked_matvec(corr, others, v))
            try:
                values, vectors = eigsh(op, k=n_comp, which="LA", v0=np.ones(sub_n) / np.sqrt(sub_n))
            except ArpackNoConvergence:                  # near-degenerate top eigenvalues: go exact
                values, vectors = np.linalg.eigh(corr[np.ix_(others, others)])
                values, vectors = values[-n_comp:], vectors[:, -n_comp:]
        keep = values > 1e-12
        if not keep.any():
            continue
        r = corr[others, j]
        beta = vectors[:, keep] @ ((vectors[:, keep].T @ r) / values[keep])
        for i, g in enumerate(rows):
            pos = np.searchsorted(others, g)
            if pos < sub_n and others[pos] == g:
                block[i, j] = abs(beta[pos])
                sign[i, j] = np.sign(beta[pos])
    return _long(regulators, expr.columns, block, sign)


def _masked_matvec(corr: np.ndarray, others: np.ndarray, v: np.ndarray) -> np.ndarray:
    full = np.zeros(corr.shape[0])
    full[others] = v
    return (corr @ full)[others]


# --------------------------------------------------------------------------- #
# tree ensembles, scikit-learn only
# --------------------------------------------------------------------------- #

def _fit_extra_trees(Xr, y, n_estimators=100, max_features="sqrt", seed=42, threads=None):
    from sklearn.ensemble import ExtraTreesRegressor

    if y.std() == 0:
        return np.zeros(Xr.shape[1])
    m = ExtraTreesRegressor(n_estimators=int(n_estimators), max_features=max_features,
                            random_state=int(seed),
                            n_jobs=-1 if threads is None else int(threads))
    m.fit(Xr, y)
    return m.feature_importances_


def _fit_gradient_boosting(Xr, y, n_estimators=100, learning_rate=0.1, seed=42):
    from sklearn.ensemble import GradientBoostingRegressor

    if y.std() == 0:
        return np.zeros(Xr.shape[1])
    m = GradientBoostingRegressor(n_estimators=int(n_estimators), learning_rate=float(learning_rate),
                                  subsample=0.9, max_features=0.3, random_state=int(seed))
    m.fit(Xr, y)
    return m.feature_importances_


def extra_trees(expr: pd.DataFrame, tf_names: list[str] | None = None,
                n_estimators: int = 100, seed: int = 42, threads: int | None = None, **_) -> pd.DataFrame:
    """Per-target extremely randomised trees importance, the Extra-Trees option of
    GENIE3 (Huynh-Thu et al. 2010, *PLoS ONE* 5:e12776); scikit-learn only."""
    return _per_target_importance(expr, tf_names, _fit_extra_trees,
                                  n_estimators=n_estimators, seed=seed, threads=threads)


def gradient_boosting(expr: pd.DataFrame, tf_names: list[str] | None = None,
                      n_estimators: int = 100, learning_rate: float = 0.1,
                      seed: int = 42, **_) -> pd.DataFrame:
    """Per-target gradient-boosted trees importance: GRNBoost2's idea (Moerman et al.
    2019, *Bioinformatics* 35:2159) with scikit-learn's boosting and no dask cluster.

    Not GRNBoost2: arboreto uses a slower learning rate with early stopping on a
    sliding window; here the number of trees is fixed. Use ``grnboost2`` to reproduce
    arboreto's scores."""
    return _per_target_importance(expr, tf_names, _fit_gradient_boosting,
                                  n_estimators=n_estimators, learning_rate=learning_rate, seed=seed)
