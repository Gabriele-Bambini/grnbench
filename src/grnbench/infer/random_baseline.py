"""Negative control (D-037): every candidate TF-target pair gets a seeded uniform random weight.

It scores the same regulator x target block as the correlation baselines (regulators resolved from
the TF list, every gene a target, self-loops dropped by the registry), so its AUPRC ratio and EPR on
the same universe are an empirical null. Their null level depends on the number of evaluation
positives: close to 1 with hundreds, higher and noisier with tens, so the control checks for a gross
universe or evaluation error and is not a test against 1. It is a control, never a member of a
comparison family. The default seed is ``grnbench.seeds.RANDOM_BASELINE_SEED``: the benchmark runners
(real_bench, candidate_study) pass no seed, so they draw with it; ``grnbench infer`` and ``grnbench run``
pass ``--seed`` as for any model. The draw is ``DRAW_RULE``.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from grnbench.seeds import RANDOM_BASELINE_SEED
from grnbench.tfs import resolve_tf_names

# What the code below does, written into the run record of every run that scores the control (D-037).
DRAW_RULE = ("uniform [0, 1) weights from numpy.random.default_rng(seed).random(n_regulators * n_genes): one draw "
             "per regulator x gene cell in regulator-major order (regulators and genes in expression-column order), "
             "self-pairs included in the draw and then dropped by the registry")


def random_baseline(expr: pd.DataFrame, tf_names: list[str] | None = None,
                    seed: int = RANDOM_BASELINE_SEED, **_) -> pd.DataFrame:
    regulators = resolve_tf_names(expr.columns, tf_names)
    src, tgt = np.meshgrid(np.asarray(regulators, dtype=object), np.asarray(expr.columns, dtype=object),
                           indexing="ij")
    weight = np.random.default_rng(seed).random(src.size)
    return pd.DataFrame({"source": src.ravel(), "target": tgt.ravel(), "weight": weight})
