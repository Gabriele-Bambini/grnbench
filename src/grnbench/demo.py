"""Synthetic demo dataset with a known regulatory network.

Raw UMI-like counts with the problems real droplet data have, so every
cleaning step has something to do, and a planted TF -> target network, so
every model can be scored against a true gold standard:

* ``n_tfs`` transcription factors, each regulating ``targets_per_tf`` genes
  (signed log-linear effects on the target's mean);
* three cell types with different TF activity programs;
* two batches with gene-wise multiplicative shifts;
* damaged cells (high mitochondrial fraction), empty-ish droplets (tiny
  libraries) and doublets (sums of two cells).

``write_demo`` writes ``counts.h5ad``, ``README.md`` (read by the context
reader), ``tfs.txt`` and ``gold.tsv`` (the planted network).
"""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

README = """# grnbench demo dataset (synthetic)

Synthetic single-cell RNA-seq, simulated to look like 10x Genomics Chromium v3
data: raw UMI counts, {n_cells} cells in two sequencing batches (column
`batch`), three cell types (column `cell_type`).

It contains damaged cells (high mitochondrial fraction), near-empty droplets
and doublets on purpose. The gene regulatory network used to simulate it is in
`gold.tsv`; the transcription factors are listed in `tfs.txt`.
"""


def make_demo(n_cells: int = 1500, n_genes: int = 2000, n_tfs: int = 30,
              targets_per_tf: int = 15, seed: int = 0):
    """Return ``(adata, tfs, gold)``; ``gold`` has columns source, target, sign."""
    rng = np.random.default_rng(seed)
    n_mito = 13
    tfs = [f"TF{i:02d}" for i in range(1, n_tfs + 1)]
    mito = [f"MT-GENE{i:02d}" for i in range(1, n_mito + 1)]
    others = [f"GENE{i:04d}" for i in range(1, n_genes - n_tfs - n_mito + 1)]
    genes = tfs + others + mito
    g_index = {g: i for i, g in enumerate(genes)}

    # planted network: each TF -> distinct-ish targets among the "others"
    edges = []
    for tf in tfs:
        for tgt in rng.choice(others[: n_tfs * targets_per_tf * 2], targets_per_tf, replace=False):
            edges.append((tf, str(tgt), int(rng.choice([-1, 1], p=[0.3, 0.7]))))
    gold = pd.DataFrame(edges, columns=["source", "target", "sign"]).drop_duplicates(
        ["source", "target"])

    # cell types with distinct TF programs; TF activity = log expression level
    cell_type = rng.choice(["alpha", "beta", "gamma"], n_cells, p=[0.45, 0.35, 0.20])
    program = {ct: rng.normal(0, 0.8, n_tfs) for ct in ("alpha", "beta", "gamma")}
    tf_log = np.stack([program[ct] for ct in cell_type]) + rng.normal(0, 0.6, (n_cells, n_tfs))
    z = (tf_log - tf_log.mean(0)) / tf_log.std(0)

    base = rng.normal(-1.0, 1.2, len(genes))
    base[: n_tfs] = rng.normal(0.8, 0.3, n_tfs)
    base[-n_mito:] = rng.normal(2.5, 0.3, n_mito)
    log_mu = np.tile(base, (n_cells, 1))
    log_mu[:, : n_tfs] += tf_log - tf_log.mean(0)
    for tf, tgt, sign in gold.itertuples(index=False):
        log_mu[:, g_index[tgt]] += sign * 0.7 * z[:, tfs.index(tf)]

    batch = rng.choice(["run1", "run2"], n_cells)
    shift = rng.normal(0, 0.25, len(genes))
    log_mu[batch == "run2"] += shift

    lib = np.exp(rng.normal(np.log(3000), 0.35, n_cells))
    lib[batch == "run2"] *= 1.3
    kind = rng.choice(["ok", "damaged", "empty", "doublet"], n_cells, p=[0.89, 0.04, 0.04, 0.03])
    lib[kind == "empty"] *= 0.03
    log_mu[kind == "damaged", -n_mito:] += np.log(8.0)
    lib[kind == "damaged"] *= 0.5

    mu = np.exp(log_mu)
    mu = mu / mu.sum(1, keepdims=True) * lib[:, None]
    theta = 5.0  # negative binomial dispersion
    X = rng.poisson(rng.gamma(theta, mu / theta)).astype(np.float32)
    dbl = np.flatnonzero(kind == "doublet")
    partners = rng.integers(0, n_cells, len(dbl))
    X[dbl] = X[dbl] + X[partners]

    obs = pd.DataFrame({"batch": batch, "cell_type": cell_type, "simulated_state": kind},
                       index=[f"cell{i:05d}" for i in range(n_cells)])
    adata = ad.AnnData(X=sparse.csr_matrix(X), obs=obs, var=pd.DataFrame(index=genes))
    return adata, tfs, gold


def write_demo(out_dir, **kwargs) -> dict[str, Path]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    adata, tfs, gold = make_demo(**kwargs)
    paths = {"counts": out / "counts.h5ad", "readme": out / "README.md",
             "tfs": out / "tfs.txt", "gold": out / "gold.tsv"}
    adata.write_h5ad(paths["counts"])
    paths["readme"].write_text(README.format(n_cells=adata.n_obs), encoding="utf-8")
    paths["tfs"].write_text("\n".join(tfs) + "\n", encoding="utf-8")
    gold.to_csv(paths["gold"], sep="\t", index=False)
    return paths
