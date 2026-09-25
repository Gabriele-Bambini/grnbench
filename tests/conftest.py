"""Shared fixtures: synthetic count datasets with known properties."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import anndata as ad
from scipy import sparse


@pytest.fixture()
def synthetic_adata() -> ad.AnnData:
    """400 cells x 600 genes raw counts.

    - 12 MT- genes (2% of space) boosted 30x in 100 cells -> elevated but
      non-failing mito fraction.
    - log-normal library sizes (CV ~1) -> healthy gates.
    - three groups in obs['group']: A(200), B(195), C(5 tiny).
    """
    rng = np.random.default_rng(0)
    n_cells, n_genes = 400, 600
    lib = np.exp(rng.normal(7.5, 0.8, n_cells)).astype(np.float32)
    base = rng.negative_binomial(2, 0.2, size=(n_cells, n_genes)).astype(np.float32)
    X = np.round(base * lib[:, None])  # integer-valued counts (float32 exact)
    genes = [f"GENE{i:03d}" for i in range(n_genes)]
    genes[:12] = [f"MT-GENE{j}" for j in range(12)]
    X[:100, :12] *= 30.0  # still integer
    obs = pd.DataFrame({"group": ["A"] * 200 + ["B"] * 195 + ["C"] * 5})
    Xs = sparse.csr_matrix(X.astype(np.float32))
    adata = ad.AnnData(X=Xs.copy(), obs=obs, var=pd.DataFrame(index=genes))
    adata.layers["counts"] = Xs
    return adata


@pytest.fixture()
def extreme_mito_adata() -> ad.AnnData:
    """50 cells where mitochondrial genes dominate -> mito gate must FAIL."""
    rng = np.random.default_rng(7)
    X = rng.poisson(50.0, size=(50, 9)).astype(np.float32)  # MT- genes huge
    normal = rng.poisson(0.5, size=(50, 1)).astype(np.float32)
    X = np.hstack([X * 20.0, normal])
    return ad.AnnData(X=X,
                      var=pd.DataFrame(index=[f"MT-G{i}" for i in range(9)] + ["ACTB"]))


@pytest.fixture()
def cascade_adata() -> ad.AnnData:
    """200 cells x 20 genes with a planted linear cascade TF_A -> G_B -> G_C.

    Remaining genes are independent Poisson noise; expression is integer
    counts so the ingest contract holds.
    """
    rng = np.random.default_rng(11)
    n_cells = 200
    A = rng.poisson(8.0, n_cells).astype(np.float64)
    B = np.round(np.maximum(0.0, 0.9 * A - 2.0 + rng.normal(0, 1.5, n_cells)))
    C = np.round(np.maximum(0.0, 0.9 * B - 2.0 + rng.normal(0, 1.5, n_cells)))
    genes = {"TF_A": A, "G_B": B, "G_C": C}
    for i in range(17):
        genes[f"N{i:02d}"] = rng.poisson(5.0, n_cells).astype(np.float64)
    cols = list(genes)
    X = np.column_stack([genes[g] for g in cols]).astype(np.float32)
    return ad.AnnData(X=X.astype(np.float32),
                      var=pd.DataFrame(index=cols))
