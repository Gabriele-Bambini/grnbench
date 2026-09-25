"""Evaluation subpackage (Parte 3): metrics, gold standards, nulls, stability."""
from __future__ import annotations

from grnbench.eval.gold import GoldStandard, coverage_report, load_gold
from grnbench.eval.metrics import (
    GoldIndex,
    auc_roc,
    average_precision,
    evaluate_ranking,
    identifier_overlap,
    label_edges,
    pr_curve,
    precision_at_k,
    top_k_overlap,
)
from grnbench.eval.nulls import degree_preserving_nulls, null_metrics, random_null
from grnbench.eval.stability import subsample_stability
from grnbench.eval.universe import align_ranking, align_rankings, load_universe
from grnbench.eval.significance import gold_thinning_sensitivity
from grnbench.eval.nulls import regulator_stratified_nulls

__all__ = [
    'align_ranking', 'align_rankings', 'load_universe',
    'gold_thinning_sensitivity', 'regulator_stratified_nulls',
    "GoldIndex",
    "GoldStandard",
    "coverage_report",
    "load_gold",
    "label_edges",
    "identifier_overlap",
    "average_precision",
    "auc_roc",
    "precision_at_k",
    "pr_curve",
    "top_k_overlap",
    "evaluate_ranking",
    "random_null",
    "degree_preserving_nulls",
    "null_metrics",
    "subsample_stability",
]
