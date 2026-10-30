"""Canonical pre-declared comparison family for real benchmark reports.

Declared before runs (2026-09-25): does the motif prior add to each base model,
does expression add to the prior alone, and do the newer models beat the Pearson default?
"""
from __future__ import annotations

DECLARED_PAIRS: list[tuple[str, str]] = [
    ("pearson_motif", "pearson_signed"),
    ("grnboost2_motif", "grnboost2"),
    ("regdiffusion_motif", "regdiffusion"),
    ("grnboost2_cistarget", "grnboost2"),
    ("pearson_motif", "motif_prior"),
    ("grnboost2", "pearson_signed"),
    ("regdiffusion", "pearson_signed"),
]
