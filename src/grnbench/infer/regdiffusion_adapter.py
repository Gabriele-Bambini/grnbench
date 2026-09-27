"""RegDiffusion adapter (Zhu & Slonim, J Comput Biol 31:1087, 2024; Apache-2.0).

A diffusion model whose learned gene x gene adjacency is the network. Input is
the log-normalized cells x genes matrix (what :func:`prepare_adata` hands every
``expr`` method); runs on CPU in about a minute for ~2,000 genes and 1,000 steps.
``adj[i, j]`` is the edge i -> j (RegDiffusion's own ``extract_node_sources``
reads column j); weight = |adj|, direction = its sign.
"""
from __future__ import annotations

import importlib.util

import numpy as np
import pandas as pd

REGDIFFUSION_UNAVAILABLE = (
    None if importlib.util.find_spec("regdiffusion") is not None
    else "regdiffusion not installed - pip install regdiffusion (needs torch)"
)


def regdiffusion_runner(expr: pd.DataFrame, tf_names: list[str] | None = None, *,
                        n_steps: int = 1000, seed: int = 0, device: str = "cpu",
                        **_) -> pd.DataFrame:
    import regdiffusion as rd

    from grnbench.infer.classical import _long, _regulators

    # RegDiffusion refuses constant genes, and a cell subsample can make a gene constant
    # (issue #13). They carry no co-variation, so they get no edges.
    expr = expr.loc[:, expr.std(axis=0).to_numpy() > 0]
    regulators = _regulators(expr, tf_names)
    X = expr.to_numpy(dtype=np.float32)
    trainer = rd.RegDiffusionTrainer(X, n_steps=int(n_steps), device=device, seed=int(seed))
    trainer.train()
    adj = np.asarray(trainer.get_adj(), dtype=float)
    rows = [expr.columns.get_loc(g) for g in regulators]
    block = adj[rows, :]
    return _long(regulators, expr.columns, np.abs(block), np.sign(block))
