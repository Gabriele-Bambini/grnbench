"""Reference inference on independent experimental-unit scores.

Paired Student t tests concern mean score differences, with Holm correction over
both directions of every pair. Bonferroni t intervals cover all pair differences
simultaneously under independent, normally distributed unit differences. These
are established methods, not a new estimator or a remedy for missing GRN labels.
"""
from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import t


def holm_adjust(pvalues) -> list[float]:
    """Holm stepdown adjusted p-values, returned in input order."""
    p = np.asarray(pvalues, dtype=float)
    if p.ndim != 1 or not np.isfinite(p).all() or ((p < 0) | (p > 1)).any():
        raise ValueError("p-values must be a finite vector in [0, 1]")
    order = np.argsort(p, kind="stable")
    adjusted = np.empty(len(p))
    adjusted[order] = np.minimum(1., np.maximum.accumulate(p[order] * np.arange(len(p), 0, -1)))
    return adjusted.tolist()


def scores_from_long(table: pd.DataFrame) -> pd.DataFrame:
    """Require exactly one observed score for every unit/method combination."""
    if not {"unit", "method", "score"}.issubset(table.columns):
        raise ValueError("scores need unit, method, score columns")
    table = table[["unit", "method", "score"]].copy()
    if table[["unit", "method"]].isna().any().any():
        raise ValueError("unit and method identifiers cannot be missing")
    table[["unit", "method"]] = table[["unit", "method"]].astype(str)
    if table.duplicated(["unit", "method"]).any():
        raise ValueError("duplicate unit/method observations; aggregate only by a justified design")
    frame = table.pivot(index="unit", columns="method", values="score")
    if frame.isna().any().any():
        raise ValueError("missing unit/method scores; incomplete pairs cannot be silently dropped")
    return frame


def compare_units(scores: pd.DataFrame, alpha: float = .05, margin: float = 0.) -> dict:
    """Compare means across independent rows; never treat edges/cells as units by default.

    A zero-variance sample supplies no estimate of sampling uncertainty: retain
    an unbounded interval and abstain, including for a constant nonzero gap.
    ``margin`` is in the input score's units and must be chosen before analysis.
    """
    if not np.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError("alpha must be between zero and one")
    if not np.isfinite(margin) or margin < 0:
        raise ValueError("margin must be finite and nonnegative")
    if not isinstance(scores, pd.DataFrame) or scores.shape[0] < 3 or scores.shape[1] < 2:
        raise ValueError("need at least three independent units and two methods")
    if not scores.index.is_unique or scores.index.hasnans or not scores.columns.is_unique or scores.columns.hasnans:
        raise ValueError("unit and method identifiers must be unique and not missing")
    names = [str(x) for x in scores.columns]
    if len(set(names)) != len(names) or any(not n.strip() for n in names):
        raise ValueError("method names must be nonempty and distinct strings")
    values = scores.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("all unit scores must be finite; missingness needs an explicit design")
    n, n_methods = values.shape
    n_pairs = n_methods * (n_methods - 1) // 2
    critical = float(t.ppf(1 - alpha / (2 * n_pairs), n - 1))
    if not np.isfinite(critical):
        raise ValueError("alpha is too small for a finite interval at this family size")
    pairs, raw_p = [], []
    for i, j in combinations(range(n_methods), 2):
        differences = values[:, i] - values[:, j]
        if not np.isfinite(differences).all():
            raise ValueError("score differences overflow; rescale the input scores")
        mean = float(differences.mean())
        se = float(differences.std(ddof=1) / np.sqrt(n))
        if not np.isfinite(mean) or not np.isfinite(se):
            raise ValueError("score moments overflow; rescale the input scores")
        if se == 0 or np.ptp(differences) == 0:
            ci, pa, pb, status = [None, None], 1., 1., "zero_variance_unresolved"
        else:
            ci = [mean - critical * se, mean + critical * se]
            pa = float(t.sf((mean - margin) / se, n - 1))
            pb = float(t.sf((-mean - margin) / se, n - 1))
            status = "estimated"
        pairs.append({"a": names[i], "b": names[j], "difference": mean,
                      "standard_error": se, "simultaneous_ci": ci,
                      "p_a_gt_b": pa, "p_b_gt_a": pb, "status": status})
        raw_p.extend([pa, pb])
    adjusted = holm_adjust(raw_p)
    for i, rec in enumerate(pairs):
        pa, pb = adjusted[2 * i:2 * i + 2]
        rec.update(p_a_gt_b_adjusted=pa, p_b_gt_a_adjusted=pb,
                   winner=rec["a"] if pa < alpha else rec["b"] if pb < alpha else None)
    return {
        "procedure": "paired_t_holm", "interval_procedure": "bonferroni_paired_t",
        "estimand": "mean difference of scores across the declared independent experimental units",
        "alpha": float(alpha), "margin": float(margin), "n_units": n,
        "methods": names, "n_pairs": n_pairs, "n_directed_hypotheses": 2 * n_pairs,
        "assumptions": ["Rows are independent experimental units, paired across all methods.",
                        "Student t inference is exact for normally distributed unit differences; otherwise approximate.",
                        "Methods, score, comparison family and practical margin were specified before evaluation.",
                        "This procedure does not correct incomplete labels or establish direct causal regulation."],
        "unresolved_is_equivalent": False, "pairs": pairs,
    }


def render_comparison(result: dict) -> str:
    lines = ["# Independent-unit score comparison", "",
             f"Units: {result['n_units']}; alpha: {result['alpha']}; margin: {result['margin']}.",
             "Paired t tests with Holm adjustment over both directions of all pairs.",
             "Intervals use simultaneous Bonferroni coverage; unresolved does not mean equivalent.", "",
             "| A | B | mean A−B | simultaneous interval | supported winner |",
             "|---|---|---:|---|---|"]
    for p in result["pairs"]:
        lo, hi = p["simultaneous_ci"]
        interval = "unbounded (zero sample variance)" if lo is None else f"[{lo:.6g}, {hi:.6g}]"
        safe = lambda s: str(s).replace("|", "\\|").replace("\n", " ")
        lines.append(f"| {safe(p['a'])} | {safe(p['b'])} | {p['difference']:.6g} | {interval} | {safe(p['winner'] or 'unresolved')} |")
    lines.extend(["", "## Assumptions", ""] + ["- " + a for a in result["assumptions"]])
    return "\n".join(lines) + "\n"
