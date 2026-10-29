"""P6 of EXP-20261016-LM-1: time RegDiffusion and GRNBoost2 on random data of RPE1's shape.

No real data. Counts are negative binomial (same shape as the RPE1 control matrix:
3,000 cells x 8,749 genes), log1p(CPM)-normalised like the engine's baseline.
RegDiffusion: n_steps steps are timed and scaled to the 1,000-step default (per-step cost is
constant). GRNBoost2: the full regulator list against a subset of targets, scaled to all
genes (one independent regression per target).

usage: python kd_p6_timing.py regdiffusion <device> <steps>
       python kd_p6_timing.py grnboost2 <n_regulators> <n_targets> <workers>
"""
import json
import sys
import time

import numpy as np

CELLS, GENES, SEED = 3000, 8749, 0


def data():
    rng = np.random.default_rng(SEED)
    mu = rng.lognormal(mean=-1.0, sigma=1.5, size=GENES)
    counts = rng.negative_binomial(n=2, p=2 / (2 + mu[None, :]), size=(CELLS, GENES)).astype(np.float32)
    lib = counts.sum(axis=1, keepdims=True)
    lib[lib == 0] = 1
    x = np.log1p(counts / lib * 1e6).astype(np.float32)
    return x[:, x.std(axis=0) > 0]


def regdiffusion(device, steps):
    import regdiffusion as rd
    import torch

    x = data()
    t0 = time.time()
    trainer = rd.RegDiffusionTrainer(x, n_steps=int(steps), device=device, seed=0)
    trainer.train()
    adj = np.asarray(trainer.get_adj(), dtype=float)
    dt = time.time() - t0
    return {"model": "regdiffusion", "device": device, "torch": torch.__version__,
            "genes": int(x.shape[1]), "cells": CELLS, "steps": int(steps), "seconds": round(dt, 1),
            "seconds_per_1000_steps": round(dt / int(steps) * 1000, 1), "adj_shape": list(adj.shape)}


def grnboost2(n_reg, n_targets, workers):
    """Through grnbench's own runner (the code path the experiment uses); threads = workers."""
    import pandas as pd
    from grnbench.infer.tree_methods import grnboost2_runner

    x = data()
    genes = [f"G{i}" for i in range(x.shape[1])]
    rng = np.random.default_rng(1)
    regs = list(rng.choice(genes, size=int(n_reg), replace=False))
    others = [g for g in genes if g not in set(regs)]
    keep = regs + list(rng.choice(others, size=int(n_targets), replace=False))
    df = pd.DataFrame(x, columns=genes)[keep]
    t0 = time.time()
    net = grnboost2_runner(df, tf_names=regs, seed=42)
    dt = time.time() - t0
    n_fit = len(keep)
    return {"model": "grnboost2", "regulators": int(n_reg), "targets_fitted": n_fit, "workers": int(workers),
            "cells": CELLS, "seconds": round(dt, 1), "edges": int(len(net)),
            "seconds_scaled_to_all_genes": round(dt / n_fit * x.shape[1], 1)}


if __name__ == "__main__":
    what = sys.argv[1]
    out = regdiffusion(sys.argv[2], sys.argv[3]) if what == "regdiffusion" else grnboost2(*sys.argv[2:5])
    print(json.dumps(out), flush=True)
