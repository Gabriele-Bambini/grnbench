"""Ranking metrics for GRN edge lists against a gold standard.

Conventions (BEELINE-style, stated explicitly so reports are auditable):
- edges arrive ranked by descending weight;
- gold = set of (source, target) pairs counted as positives;
- prevalence pi = (#gold positives in the fixed support) / (#candidates) is
  the random precision baseline, not the finite-sample expectation of AP;
- AUPRC-ratio = AP / pi;  EPR = precision@k / pi with k = #gold hits.

Gene identifiers reach grnbench from two independent sources - the expression
matrix and the ground-truth network - and they disagree on case more often than
anyone expects. BEELINE ships mESC expression in MGI case (``Pou5f1``) and every
mouse network in upper case (``POU5F1``); exact string matching turned 64,557
reachable gold edges into zero hits, an AUPRC of NaN and not one line of output.
Matching is therefore case-insensitive by default, and the folding is *reported*
(see :func:`identifier_overlap`) rather than quietly assumed.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from grnbench.eval.universe import canonical_pairs, align_ranking, align_rankings, top_pairs

GoldLike = "set[tuple[str, str]] | frozenset[tuple[str, str]] | GoldIndex"


class GoldIndex:
    """Gold pairs prepared once for matching against many edge lists.

    Null models score the same gold set tens of times; normalising the pairs
    inside ``label_edges`` would redo that work on every draw, which on a
    million-pair ChIP-seq network is the whole runtime. Build this once and pass
    it down instead.
    """

    __slots__ = ("pairs", "case_sensitive", "lookup")

    def __init__(self, pairs, *, case_sensitive: bool = False) -> None:
        raw = list(pairs)
        checked = canonical_pairs(raw, case_sensitive=case_sensitive, name='gold', allow_self_loops=True)
        self.pairs = frozenset((str(s), str(t)) for s, t in raw)
        self.case_sensitive = bool(case_sensitive)
        self.lookup = frozenset(map(tuple, checked.to_numpy()))

    def __len__(self) -> int:
        return len(self.pairs)

    @classmethod
    def coerce(cls, gold: GoldLike, *, case_sensitive: bool = False) -> "GoldIndex":
        if isinstance(gold, cls):
            return gold
        return cls(gold, case_sensitive=case_sensitive)

    def labels(self, edges: pd.DataFrame) -> np.ndarray:
        if self.case_sensitive:
            keys = zip(edges["source"], edges["target"])
        else:
            keys = zip(edges["source"].astype(str).str.upper(),
                       edges["target"].astype(str).str.upper())
        look = self.lookup
        return np.fromiter((k in look for k in keys), dtype=np.int8,
                           count=len(edges))


def label_edges(edges: pd.DataFrame, gold: GoldLike,
                *, case_sensitive: bool = False) -> np.ndarray:
    return GoldIndex.coerce(gold, case_sensitive=case_sensitive).labels(edges)


def identifier_overlap(edges: pd.DataFrame, gold: GoldLike,
                       n_examples: int = 5) -> dict:
    """Do the edge list and the gold standard even speak the same gene names?

    Returned verdicts:

    ``ok``
        the two agree on case; nothing to report.
    ``case-mismatch``
        folding case recovers hits that exact matching loses - the numbers
        below are computed with folding, and the report says so.
    ``no-overlap``
        the identifiers do not match even case-insensitively, so every metric
        against this gold is vacuous. Different symbol namespaces (Ensembl IDs
        vs symbols), the wrong organism, or the wrong gold file.
    """
    # deliberately independent of how ``gold`` was configured: the whole point
    # is to compare exact against folded matching, so a case-sensitive
    # GoldIndex must not be able to hide the very mismatch being looked for
    pairs = gold.pairs if isinstance(gold, GoldIndex) else frozenset(gold)
    folded_pairs = frozenset((str(s).upper(), str(t).upper()) for s, t in pairs)
    src = edges["source"].astype(str)
    tgt = edges["target"].astype(str)
    exact = int(np.fromiter((k in pairs for k in zip(src, tgt)),
                            dtype=np.int8, count=len(edges)).sum())
    folded = int(np.fromiter((k in folded_pairs
                              for k in zip(src.str.upper(), tgt.str.upper())),
                             dtype=np.int8, count=len(edges)).sum())

    edge_genes = set(src) | set(tgt)
    gold_genes = {g for pair in pairs for g in pair}
    genes_exact = len(edge_genes & gold_genes)
    genes_folded = len({g.upper() for g in edge_genes} & {g.upper() for g in gold_genes})

    if folded == 0 and genes_folded == 0:
        verdict = "no-overlap"
    elif folded > exact:
        verdict = "case-mismatch"
    else:
        verdict = "ok"
    out = {
        "verdict": verdict,
        "matched_pairs_case_sensitive": exact,
        "matched_pairs_case_insensitive": folded,
        "gold_genes": len(gold_genes),
        "edge_genes": len(edge_genes),
        "shared_genes_case_sensitive": genes_exact,
        "shared_genes_case_insensitive": genes_folded,
    }
    if verdict != "ok":
        out["example_edge_genes"] = sorted(edge_genes)[:n_examples]
        out["example_gold_genes"] = sorted(gold_genes)[:n_examples]
    return out


def _tie_blocks(weights: np.ndarray) -> np.ndarray:
    """End positions (1-indexed, inclusive) of each run of equal weight.

    ``weights`` must already be in ranking order. Equal scores form one block:
    nothing inside a block is evidence that one edge outranks another.
    """
    w = np.asarray(weights, dtype=float)
    if w.size == 0:
        return np.empty(0, dtype=int)
    # NaN != NaN would split a run of NaNs into singletons; treat them as equal
    same = (w[1:] == w[:-1]) | (np.isnan(w[1:]) & np.isnan(w[:-1]))
    return np.flatnonzero(np.r_[~same, True]) + 1


def average_precision(y: np.ndarray, weights: np.ndarray | None = None) -> float:
    """Average precision.

    With ``weights``, precision and recall are read only at the boundaries
    between distinct scores, so the result no longer depends on the arbitrary
    order of edges that scored identically. Without them the old row-order
    behaviour is kept, which is correct exactly when every score is distinct.

    This matters more than it sounds: on a ranking with four score levels, six
    orderings of the *same* scores produced AUPRC between 0.319 and 0.478 and
    EPR between 3.0 and 5.0, while the tie-safe AUROC stayed at 0.9167 for all
    six - proof that the ranking was identical and only the tie-breaking moved.
    """
    y = np.asarray(y)
    n_pos = int(y.sum())
    if n_pos == 0:
        return float("nan")
    if weights is None:
        cdf = np.cumsum(y)
        precision_at_hits = cdf[y == 1] / np.arange(1, len(y) + 1)[y == 1]
        return float(precision_at_hits.mean())
    ends = _tie_blocks(weights)
    cdf = np.cumsum(y)
    hits = cdf[ends - 1].astype(float)          # hits at each block boundary
    counts = ends.astype(float)
    precision = hits / counts
    recall = hits / n_pos
    return float(np.sum(np.diff(np.r_[0.0, recall]) * precision))


def auc_roc(weights: np.ndarray, y: np.ndarray) -> float:
    """Rank-based AUROC (Mann-Whitney U), tie-safe via average ranks."""
    n_pos = int(y.sum())
    n_neg = len(y) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    ranks = rankdata(weights, method="average")
    u = ranks[y == 1].sum() - n_pos * (n_pos + 1) / 2.0
    return float(u / (n_pos * n_neg))


def precision_at_k(y: np.ndarray, k: int,
                   weights: np.ndarray | None = None) -> float:
    """Precision in the top ``k``.

    With ``weights``, this is the *expected* precision under random tie-breaking
    when the cut at ``k`` falls inside a block of equal scores. Taking the first
    ``k`` rows instead credits a method for an ordering it never asserted - the
    scores say those edges are indistinguishable.
    """
    y = np.asarray(y)
    if k <= 0 or len(y) == 0:
        return float("nan")
    k = min(int(k), len(y))
    if weights is None:
        return float(y[:k].mean())
    lo, hi, m = _boundary_block(weights, k)
    if m == 0:                                   # cut lands on a block edge
        return float(y[:k].mean())
    g = int(y[lo:hi].sum())                      # positives inside the block
    expected = float(y[:lo].sum()) + g * (k - lo) / m
    return expected / k


def _boundary_block(weights: np.ndarray, k: int) -> tuple[int, int, int]:
    """Locate the tied block straddling the cut at ``k``.

    Returns ``(lo, hi, size)``: ``lo`` edges sit entirely above the block, which
    spans ``y[lo:hi]``. ``size == 0`` means the cut coincides with a block
    boundary, so no tie-breaking is involved and the plain top-``k`` is exact.
    """
    ends = _tie_blocks(weights)
    if ends.size == 0 or k in set(ends.tolist()):
        return k, k, 0
    idx = int(np.searchsorted(ends, k, side="left"))
    hi = int(ends[idx])
    lo = int(ends[idx - 1]) if idx > 0 else 0
    return lo, hi, hi - lo


def precision_at_k_bounds(y: np.ndarray, k: int,
                          weights: np.ndarray) -> tuple[float, float]:
    """Worst and best precision@k achievable by reordering tied edges only.

    A wide interval means the headline number is a coin flip, not a measurement.
    """
    y = np.asarray(y)
    if k <= 0 or len(y) == 0:
        return float("nan"), float("nan")
    k = min(int(k), len(y))
    lo, hi, m = _boundary_block(weights, k)
    if m == 0:
        p = float(y[:k].sum()) / k
        return p, p
    base = float(y[:lo].sum())
    g = int(y[lo:hi].sum())
    taken = k - lo
    worst = base + max(0, g - (m - taken))     # tied positives pushed below the cut
    best = base + min(g, taken)                # tied positives pulled above it
    return worst / k, best / k


def pr_curve(y: np.ndarray, max_points: int = 200,
             weights: np.ndarray | None = None) -> dict:
    """Precision/recall along the ranking, thinned to ``max_points`` for plotting.

    Precision is evaluated at every rank — or, when ``weights`` are given, only
    where the score changes, matching the tie-aware :func:`average_precision`.
    The returned series is subsampled on a uniform grid (the last point is
    always kept) so an HTML report can draw the curve without shipping a million
    coordinates.
    """
    y = np.asarray(y, dtype=float)
    n_pos = float(y.sum())
    if len(y) == 0 or n_pos == 0:
        return {"recall": [], "precision": [], "baseline": 0.0, "n_pos": 0}
    cdf = np.cumsum(y)
    ranks = np.arange(1, len(y) + 1, dtype=float)
    if weights is not None:
        # read the curve only where the score actually changes, so it agrees
        # with the tie-aware AUPRC instead of drawing the sawtooth that the
        # arbitrary within-block order happens to trace
        ends = _tie_blocks(weights)
        cdf, ranks = cdf[ends - 1], ranks[ends - 1]
    precision = cdf / ranks
    recall = cdf / n_pos
    # thin over the series actually built: collapsing ties makes it shorter than
    # the edge list, and indexing it by len(y) would run off the end
    n_points = len(precision)
    if n_points > max_points:
        idx = np.unique(np.linspace(0, n_points - 1, max_points).astype(int))
    else:
        idx = np.arange(n_points)
    return {
        "recall": [round(float(v), 6) for v in recall[idx]],
        "precision": [round(float(v), 6) for v in precision[idx]],
        "baseline": round(float(n_pos / len(y)), 6),
        "n_pos": int(n_pos),
    }


def top_k_overlap(edge_sets: dict[str, pd.DataFrame], k: int = 1000,
                  *, case_sensitive: bool = False, universe=None) -> dict:
    """Pairwise Jaccard of the top-k edges of each method (methods agreeing?).

    Two methods scoring identically against a gold standard while sharing no
    edges are not two confirmations of the same biology; the report says so.
    That claim only holds if the comparison is real: edge lists compared here
    can come from different runs and different tools, so they are folded to a
    common case exactly as gold matching is. Left case-sensitive, two networks
    that differ only in ``Pou5f1`` vs ``POU5F1`` scored Jaccard 0.000 under a
    heading asserting the two methods had found nothing in common.
    """
    aligned = align_rankings(edge_sets, universe, case_sensitive=case_sensitive)
    tops = {name: top_pairs(df, k, case_sensitive=case_sensitive)
            for name, df in aligned.items()}
    names = list(tops)
    matrix = {}
    for a in names:
        row = {}
        for b in names:
            union = tops[a] | tops[b]
            row[b] = round(len(tops[a] & tops[b]) / len(union), 4) if union else float("nan")
        matrix[a] = row
    return {"k": int(k), "methods": names, "jaccard": matrix,
            "tie_policy": "include-complete-boundary-tie",
            "top_k_sizes": {n: len(s) for n, s in tops.items()}}


def metric_from_labels(metric: str, weights: np.ndarray, y: np.ndarray,
                       k: int) -> float:
    """One scalar metric from a ranking-ordered score vector and its labels.

    ``weights`` must be in ranking order (descending) and ``y`` aligned to it.
    This is the piece of :func:`evaluate_ranking` a null model needs thousands
    of times, without re-validating identifiers on every call.
    """
    y = np.asarray(y)
    n_pos = int(y.sum())
    if n_pos == 0 or len(y) == 0:
        return float("nan")
    prevalence = n_pos / len(y)
    if metric == "auroc":
        return auc_roc(weights, y)
    if metric in ("epr", "precision_at_k"):
        p = precision_at_k(y, k, weights)
        return p / prevalence if metric == "epr" else p
    if metric in ("auprc", "auprc_ratio"):
        ap = average_precision(y, weights)
        return ap / prevalence if metric == "auprc_ratio" else ap
    raise ValueError(f"unknown scalar metric {metric!r}")


def evaluate_ranking(edges: pd.DataFrame,
                     gold: GoldLike,
                     top_k: int | None = None,
                     with_curve: bool = True, *, universe=None) -> dict:
    """Full metric block for one ranked edge list.

    ``gold`` may be a bare set of pairs or a prebuilt :class:`GoldIndex`; pass
    the latter when scoring the same gold repeatedly.

    ``with_curve=False`` skips the PR curve — null models call this thousands of
    times and only ever read a scalar.
    """
    gi = GoldIndex.coerce(gold)
    input_hits = int(gi.labels(edges).sum())
    edges = align_ranking(edges, universe, case_sensitive=gi.case_sensitive)
    if top_k is not None and (isinstance(top_k, bool) or int(top_k) != top_k or top_k <= 0):
        raise ValueError('top_k must be a positive integer')
    y = gi.labels(edges)
    n_pos = int(y.sum())
    n_edges = len(y)
    prevalence = n_pos / n_edges
    k = n_pos if top_k is None else min(int(top_k), n_edges)

    w = edges["weight"].to_numpy(dtype=float)
    ap = average_precision(y, w)
    p_at_k = precision_at_k(y, k, w)
    out = {
        "n_edges": n_edges,
        **edges.attrs['evaluation_support'],
        "n_gold_in_list": input_hits,
        "n_gold_in_universe": n_pos,
        "n_gold_total": len(gi),
        "gold_outside_universe": len(gi) - n_pos,
        "gold_reachability": n_pos / len(gi) if len(gi) else float('nan'),
        "negative_label_policy": "unlabelled candidates treated as negatives; not certified biological negatives",
        "ap_baseline_note": "prevalence is a conventional baseline, not exact finite-sample random AP",
        "prevalence": prevalence,
        "auprc": ap,
        "auprc_ratio": (ap / prevalence) if prevalence > 0 else float("nan"),
        "epr_k": k,
        "precision_at_k": p_at_k,
        "epr": None,
        "auroc": auc_roc(w, y),
    }
    out["epr"] = (p_at_k / prevalence) if prevalence > 0 and \
        not np.isnan(p_at_k) else float("nan")
    # How much of the headline number is decided by edges the method scored
    # identically? A wide interval means the top-k cut fell inside a tied block
    # and EPR is a coin flip, not a measurement.
    lo_p, hi_p = precision_at_k_bounds(y, k, w)
    out["tie_block_at_k"] = int(_boundary_block(w, k)[2])
    out["epr_range"] = ([lo_p / prevalence, hi_p / prevalence]
                        if prevalence > 0 else [float("nan"), float("nan")])
    if with_curve:
        out["pr_curve"] = pr_curve(y, weights=w)
    return out
