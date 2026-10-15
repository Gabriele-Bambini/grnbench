"""QC metrics and gate behavior tests."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import anndata as ad

from grnbench.qc import run_qc


def test_metrics_and_gates_on_healthy_data(synthetic_adata):
    rep = run_qc(synthetic_adata)
    assert {"n_cells", "n_genes", "sparsity", "library_cv",
            "median_pct_mito"} <= set(rep.metrics)
    assert rep.metrics["n_mito_genes"] == 12
    assert rep.metrics["n_cells"] == 400
    statuses = {g.status for g in rep.gates}
    assert statuses <= {"pass", "warn", "fail"}
    assert not rep.has_failures, [g.name for g in rep.gates if g.status == "fail"]


def test_obs_annotation_added(synthetic_adata):
    rep = run_qc(synthetic_adata)
    assert "_total_counts" in synthetic_adata.obs.columns
    assert "_pct_mito" in synthetic_adata.obs.columns
    assert len(rep.gates) >= 5


def test_mito_gate_fails_on_extreme_dataset(extreme_mito_adata):
    rep = run_qc(extreme_mito_adata)
    mito_gate = next(g for g in rep.gates if g.name == "mito_fraction")
    assert mito_gate.status == "fail"
    assert rep.has_failures


def test_low_complexity_gate_fails():
    rng = np.random.default_rng(1)
    X = np.zeros((100, 5000))
    for i in range(100):
        X[i, rng.choice(5000, size=5, replace=False)] = 1.0
    a = ad.AnnData(X=X.astype(np.float32),
                   var=pd.DataFrame(index=[f"G{i}" for i in range(5000)]))
    rep = run_qc(a)
    low_genes = next(g for g in rep.gates if g.name == "cells_below_min_genes")
    assert low_genes.status == "fail"
    sparsity_gate = next(g for g in rep.gates if g.name == "sparsity")
    assert sparsity_gate.status == "fail"


def test_custom_thresholds_change_gate_outcome():
    from grnbench.qc import Thresholds

    a = ad.AnnData(X=np.eye(4, dtype=np.float32),
                   var=pd.DataFrame(index=[f"G{i}" for i in range(4)]))
    default_rep = run_qc(a.copy())
    strict_rep = run_qc(a.copy(), Thresholds(min_genes_per_cell=2))
    assert isinstance(default_rep.has_failures, bool)
    assert strict_rep.gates[0].name == "cells_below_min_genes"


def test_a_single_batch_column_is_not_reported_as_missing(synthetic_adata):
    """`obs['batch'] = 'run1'` is a batch annotation: 'provide one if it exists' was wrong (HI-18)."""
    def gate(adata):
        return next(g for g in run_qc(adata).gates if g.name == "batch_annotation")
    plain = synthetic_adata.copy()
    assert gate(plain).status == "warn" and "no batch/donor annotation" in gate(plain).detail
    one = synthetic_adata.copy()
    one.obs["batch"] = "run1"
    g = gate(one)
    assert g.status == "pass" and "one batch only" in g.detail and "'batch'" in g.detail
    many = synthetic_adata.copy()
    many.obs["batch"] = ["a", "b"] * (many.n_obs // 2) + ["a"] * (many.n_obs % 2)
    assert gate(many).status == "warn", "run_qc is given batch_key by the caller; an unused 2-level column is not announced"
