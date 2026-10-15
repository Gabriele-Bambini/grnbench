"""Co-expression and information-theoretic methods from the GRN literature.

Native numpy/scipy reimplementations: no extra dependency, deterministic, CPU
minutes at the sizes the engine works at (a few thousand cells, ~2,000 genes,
~150 regulators). They follow the published description of each method; they are
**not** the original packages, and each docstring says where they differ. None of
them is registered in a pre-declared comparison family: until an experiment record
declares them, scores they produce are exploratory.

``bicor``       biweight midcorrelation (WGCNA), robust to outlier cells.
``wgcna_tom``   WGCNA: soft-thresholded |correlation| and topological overlap.
``propr_rho``   proportionality rho_p (Quinn et al. 2017), a log-ratio association.
``aracne_dpi``  ARACNe-style: binned MI, significance threshold, DPI pruning.
``mrnet``       MRNET (Meyer et al. 2007): maximum relevance, minimum redundancy.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from grnbench.infer.classical import (
    _equal_frequency_bins,
    _long,
    _mi_matrix,
    _regulators,
)


def _error(message: str):
    from grnbench.infer.base import MethodError
    return MethodError(message)


def _tf_index(expr: pd.DataFrame, tf_names) -> tuple[list[str], np.ndarray]:
    regulators = _regulators(expr, tf_names)
    return regulators, np.array([expr.columns.get_loc(g) for g in regulators], dtype=int)


# --------------------------------------------------------------------------- #
# bicor
# --------------------------------------------------------------------------- #

def _bicor_block(X: np.ndarray, rows: np.ndarray, c: float = 9.0) -> np.ndarray:
    """Biweight midcorrelation of ``rows`` columns against every column.

    Langfelder & Horvath, *J Stat Softw* 46(11) 2012: each gene is centred on its
    median, cells are down-weighted by ``(1 - u^2)^2`` with ``u = dev / (c MAD)``
    (zero beyond ``|u| = 1``), and the weighted deviations are correlated. A gene
    whose MAD is zero (at least half the cells at one value, the rule for sparse
    scRNA-seq genes) has no robust scale; every pair involving such a gene uses
    Pearson, WGCNA's ``pearsonFallback = "individual"``.
    """
    med = np.median(X, axis=0)
    dev = X - med
    mad = np.median(np.abs(dev), axis=0)
    robust = mad > 0
    u = np.zeros_like(X)
    u[:, robust] = dev[:, robust] / (c * mad[robust])
    weighted = dev * np.where(np.abs(u) < 1.0, (1.0 - u ** 2) ** 2, 0.0)
    weighted[:, ~robust] = 0.0
    norm = np.linalg.norm(weighted, axis=0)
    zb = np.divide(weighted, norm, out=np.zeros_like(weighted), where=norm > 0)

    centred = X - X.mean(axis=0)
    pnorm = np.linalg.norm(centred, axis=0)
    zp = np.divide(centred, pnorm, out=np.zeros_like(centred), where=pnorm > 0)

    biweight = zb[:, rows].T @ zb
    pearson = zp[:, rows].T @ zp
    both_robust = robust[rows][:, None] & robust[None, :]
    return np.nan_to_num(np.where(both_robust, biweight, pearson), nan=0.0)


def bicor(expr: pd.DataFrame, tf_names: list[str] | None = None, **_) -> pd.DataFrame:
    regulators, rows = _tf_index(expr, tf_names)
    block = _bicor_block(expr.to_numpy(dtype=np.float64), rows)
    return _long(regulators, expr.columns, np.abs(block), np.sign(block))


# --------------------------------------------------------------------------- #
# wgcna
# --------------------------------------------------------------------------- #

WGCNA_POWERS = tuple(range(1, 11)) + tuple(range(12, 21, 2))


def _scale_free_fit(k: np.ndarray, n_breaks: int = 10) -> tuple[float, float]:
    """(R^2, slope) of log10 p(k) on log10 k, as WGCNA's ``scaleFreeFitIndex``."""
    k = np.asarray(k, dtype=float)
    lo, hi = float(k.min()), float(k.max())
    if hi <= lo:
        return 0.0, 0.0
    edges = np.linspace(lo, hi, n_breaks + 1)
    mids = (edges[:-1] + edges[1:]) / 2
    which = np.clip(np.searchsorted(edges, k, side="left") - 1, 0, n_breaks - 1)
    counts = np.bincount(which, minlength=n_breaks).astype(float)
    sums = np.bincount(which, weights=k, minlength=n_breaks)
    mean_k = np.where(counts > 0, sums / np.maximum(counts, 1), mids)
    mean_k = np.where(mean_k == 0, mids, mean_k)
    x = np.log10(mean_k)
    y = np.log10(counts / len(k) + 1e-9)
    if np.ptp(x) == 0:
        return 0.0, 0.0
    slope, intercept = np.polyfit(x, y, 1)
    resid = y - (slope * x + intercept)
    total = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - float((resid ** 2).sum()) / total if total > 1e-12 else 0.0
    return r2, float(slope)


def pick_soft_power(corr_abs: np.ndarray, r2_cut: float = 0.85) -> tuple[int, dict]:
    """WGCNA ``pickSoftThreshold``: smallest power whose network is scale-free enough.

    The criterion is ``-sign(slope) * R^2 > r2_cut``. When no power reaches it WGCNA
    leaves the choice to the user; here the power with the best signed R^2 is used
    and the fit is returned so the caller can record that the cut was missed.
    """
    fits = {}
    for power in WGCNA_POWERS:
        adj = corr_abs ** power
        np.fill_diagonal(adj, 0.0)
        r2, slope = _scale_free_fit(adj.sum(axis=1))
        fits[power] = -np.sign(slope) * r2
    for power in WGCNA_POWERS:
        if fits[power] > r2_cut:
            return power, {"signed_r2": fits[power], "cut_reached": True}
    best = max(fits, key=lambda p: fits[p])
    return best, {"signed_r2": fits[best], "cut_reached": False}


def wgcna_tom(expr: pd.DataFrame, tf_names: list[str] | None = None,
              power: int | None = None, max_genes: int = 6000, **_) -> pd.DataFrame:
    """WGCNA topological overlap (Zhang & Horvath 2005; Langfelder & Horvath 2008).

    Unsigned adjacency ``a_ij = |cor_ij|^beta`` with ``beta`` from
    :func:`pick_soft_power` unless ``power`` is given, then
    ``TOM_ij = (sum_u a_iu a_uj + a_ij) / (min(k_i, k_j) + 1 - a_ij)``. TOM is
    symmetric, so it ranks TF-gene pairs but does not orient them; the sign of the
    correlation is the ``direction``. The power actually used, how it was chosen and the
    scale-free fit are returned in ``attrs`` and recorded in the run's provenance. Builds a
    genes x genes matrix, so it refuses more than ``max_genes`` genes instead of exhausting
    memory.
    """
    if expr.shape[1] > int(max_genes):
        raise _error(f"wgcna_tom builds a genes x genes matrix: {expr.shape[1]} genes exceed "
                     f"max_genes={int(max_genes)}; keep fewer genes or raise max_genes")
    regulators, rows = _tf_index(expr, tf_names)
    X = expr.to_numpy(dtype=np.float64)
    centred = X - X.mean(axis=0)
    norm = np.linalg.norm(centred, axis=0)
    z = np.divide(centred, norm, out=np.zeros_like(centred), where=norm > 0)
    corr = np.clip(z.T @ z, -1.0, 1.0)
    cabs = np.abs(corr)
    if power is not None:
        beta, info = int(power), {"power": int(power), "power_selection": "given"}
    else:                                                      # WGCNA drops constant genes first
        live = norm > 0
        beta, fit = pick_soft_power(cabs[np.ix_(live, live)])
        info = {"power": int(beta), "signed_scale_free_r2": round(float(fit["signed_r2"]), 4),
                "power_selection": "scale_free_cut_reached" if fit["cut_reached"]
                else "fallback_best_signed_r2 (the 0.85 cut was not reached)"}
    adj = cabs ** beta
    np.fill_diagonal(adj, 0.0)
    k = adj.sum(axis=1)
    shared = adj[rows] @ adj                                   # sum_u a_iu a_uj, a_ii = 0
    denom = np.minimum(k[rows][:, None], k[None, :]) + 1.0 - adj[rows]
    tom = (shared + adj[rows]) / np.where(denom > 0, denom, 1.0)
    out = _long(regulators, expr.columns, tom, np.sign(corr[rows]))
    out.attrs.update(info)                                     # run_method records it in the provenance
    return out


# --------------------------------------------------------------------------- #
# proportionality
# --------------------------------------------------------------------------- #

def propr_rho(expr: pd.DataFrame, tf_names: list[str] | None = None, **_) -> pd.DataFrame:
    """Proportionality rho_p = 1 - var(x - y) / (var(x) + var(y)) = 2 cov / (var x + var y).

    Defined by Erb & Notredame (2016) and implemented in ``propr`` (Quinn et al., *Sci Rep*
    7:16252, 2017); proportionality measures were among the best performers in the scRNA-seq
    benchmark of association measures (Skinnider et al., *Nat Methods* 16, 2019). ``propr``
    applies it to centred log-ratio transformed counts; this adapter applies it to the
    normalised log expression every method here receives, so it differs from ``propr`` in that
    transform. Unlike correlation it also penalises genes with different variances.
    """
    regulators, rows = _tf_index(expr, tf_names)
    X = expr.to_numpy(dtype=np.float64)
    centred = X - X.mean(axis=0)
    cov = centred[:, rows].T @ centred / max(X.shape[0] - 1, 1)
    var = (centred ** 2).sum(axis=0) / max(X.shape[0] - 1, 1)
    denom = var[rows][:, None] + var[None, :]
    rho = np.divide(2.0 * cov, denom, out=np.zeros_like(cov), where=denom > 0)
    return _long(regulators, expr.columns, np.abs(rho), np.sign(rho))


# --------------------------------------------------------------------------- #
# ARACNe-style
# --------------------------------------------------------------------------- #

def aracne_dpi(expr: pd.DataFrame, tf_names: list[str] | None = None,
               n_bins: int = 8, pvalue: float = 1e-8, tolerance: float = 0.0,
               **_) -> pd.DataFrame:
    """ARACNe-style network: MI, a significance threshold, then the DPI.

    Margolin et al., *BMC Bioinformatics* 7(Suppl 1):S7 (2006). Differences from the
    ARACNe-AP / ARACNe3 programs (method ``aracne3``): MI is the equal-frequency
    plug-in estimator the other MI methods use, not adaptive partitioning or a
    kernel estimate; the MI threshold comes from the asymptotic G-test
    (``2 n MI ~ chi2((b_i - 1)(b_j - 1))``, ``b`` the occupied bins) at ``pvalue``,
    not from permutations; there is no bootstrap consolidation.

    The DPI is applied over TF-TF-target triangles: an edge ``TF i -> j`` is removed
    when another TF ``k`` has both ``I(i,k)`` and ``I(k,j)`` significant and no smaller
    than ``I(i,j) / (1 - tolerance)`` (ties go, as in ARACNe's ``<=``). Removed and non-significant edges get weight 0
    and tie at the bottom, as in an ARACNe output that simply omits them.
    """
    from scipy.stats import chi2

    regulators, rows = _tf_index(expr, tf_names)
    X = expr.to_numpy(dtype=np.float64)
    n_cells = X.shape[0]
    mi = _mi_matrix(X, rows, int(n_bins))
    codes = _equal_frequency_bins(X, int(n_bins))
    occupied = np.array([len(np.unique(codes[:, j])) for j in range(X.shape[1])])
    dof = (occupied[rows][:, None] - 1) * (occupied[None, :] - 1)
    with np.errstate(invalid="ignore"):
        pvals = np.where(dof > 0, chi2.sf(2.0 * n_cells * mi, np.maximum(dof, 1)), 1.0)
    significant = pvals < float(pvalue)
    mi_sig = np.where(significant, mi, 0.0)
    mi_sig[np.arange(len(rows)), rows] = 0.0                    # a TF is not its own target

    tf_mi = mi_sig[:, rows]                                     # TF x TF
    np.fill_diagonal(tf_mi, 0.0)
    keep = significant.copy()
    keep[np.arange(len(rows)), rows] = False
    for i in range(len(rows)):
        via = np.minimum(tf_mi[i][:, None], mi_sig).max(axis=0)   # best indirect path i -> k -> j
        keep[i] &= mi_sig[i] > via * (1.0 - float(tolerance))
    return _long(regulators, expr.columns, np.where(keep, mi, 0.0))


# --------------------------------------------------------------------------- #
# MRNET
# --------------------------------------------------------------------------- #

def mrnet(expr: pd.DataFrame, tf_names: list[str] | None = None,
          n_bins: int = 8, **_) -> pd.DataFrame:
    """MRNET (Meyer, Kontos, Lafitte & Bontempi 2007, *EURASIP J Bioinform Syst Biol* 2007:79879;
    implemented in the R package ``minet``).

    For every target, regulators are chosen one at a time by maximum relevance minus
    mean redundancy with those already chosen, ``s_i = I(i;t) - mean_k I(i;k)``; the
    selection stops when the best score is not positive. The weight is the score at
    which a regulator was chosen; regulators never chosen get 0 and tie at the
    bottom. ``minet`` symmetrises the result by taking the larger of the two
    directions; this adapter keeps the target-wise direction, which is what a
    TF -> gene ranking needs.
    """
    regulators, rows = _tf_index(expr, tf_names)
    X = expr.to_numpy(dtype=np.float64)
    n_tf, n_genes = len(rows), X.shape[1]
    mi = _mi_matrix(X, rows, int(n_bins))                       # TF x gene
    redundancy = mi[:, rows].copy()                             # TF x TF
    np.fill_diagonal(redundancy, 0.0)
    where_tf = {int(g): i for i, g in enumerate(rows)}
    block = np.zeros((n_tf, n_genes))
    for j in range(n_genes):
        relevance = mi[:, j]
        free = np.ones(n_tf, dtype=bool)
        if j in where_tf:
            free[where_tf[j]] = False
        redundant_sum = np.zeros(n_tf)
        chosen = 0
        for _ in range(n_tf):
            score = relevance - (redundant_sum / chosen if chosen else 0.0)
            score = np.where(free, score, -np.inf)
            best = int(np.argmax(score))
            if not np.isfinite(score[best]) or score[best] <= 0.0:
                break
            block[best, j] = score[best]
            free[best] = False
            redundant_sum += redundancy[:, best]
            chosen += 1
    return _long(regulators, expr.columns, block)
