"""Auto-recommendation of a preprocessing strategy from dataset statistics.

Heuristics only — every recommendation carries its reason and is recorded in
the run manifest. The user can always override with an explicit strategy name.
"""
from __future__ import annotations


def recommend_strategy(metrics: dict,
                       group_key: str | None = None,
                       n_groups: int | None = None,
                       median_cells_per_group: float | None = None) -> list[dict]:
    """Return ranked recommendations: ``[{strategy, rank, reason}, ...]``."""
    recs: list[dict] = []
    sparsity = metrics.get("sparsity", 0.0)
    lib_cv = metrics.get("library_cv", 0.0)
    median_cells = metrics.get("median_counts_per_cell", 0.0)

    if group_key is not None and n_groups is not None and median_cells_per_group is not None:
        if median_cells_per_group < 30 or sparsity >= 0.95:
            recs.append({
                "strategy": "pseudobulk",
                "reason": (
                    f"groups average {median_cells_per_group:.0f} cells with "
                    f"{sparsity:.0%} zero fraction: per-cell correlations are "
                    "noise-dominated, pseudobulk stabilizes them"
                ),
            })
    if lib_cv > 2.0 or median_cells < 300:
        recs.append({
            "strategy": "pearson_residuals",
            "reason": (
                f"library-size CV={lib_cv:.2f}, median counts/cell="
                f"{median_cells:.0f}: Pearson residuals stabilize variance "
                "across depth differences"
            ),
        })
    recs.append({
        "strategy": "cpm_log1p",
        "reason": "safe default for medium-depth data without strong imbalance",
    })
    recs.append({
        "strategy": "raw_counts",
        "reason": "for methods that normalize internally (e.g. some MI estimators)",
    })

    seen: set[str] = set()
    ranked: list[dict] = []
    for i, r in enumerate(recs):
        if r["strategy"] not in seen:
            r["rank"] = i + 1
            ranked.append(r)
            seen.add(r["strategy"])
    return ranked
