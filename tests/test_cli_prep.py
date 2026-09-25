"""Tests for the grnbench-prep CLI wiring (grnbench.cli.main).

Nothing previously exercised this module directly: every other preprocessing
test calls ``run_prep`` itself, so an argument-marshalling bug in ``main()``
(the ``x or default`` treatment of an explicit ``0``, below) had no test that
could ever catch it.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
import pytest

from grnbench.cli import main


def _write_adata(path, n_cells=100, n_genes=400):
    rng = np.random.default_rng(0)
    X = rng.poisson(3.0, size=(n_cells, n_genes)).astype(np.float32)
    ad.AnnData(X=sp.csr_matrix(X),
              var=pd.DataFrame(index=[f"G{i}" for i in range(n_genes)])
              ).write_h5ad(path)


def test_explicit_zero_threshold_is_honoured_not_replaced_by_default(tmp_path):
    """``args.min_genes_per_cell or Thresholds.min_genes_per_cell`` used to read
    an explicit ``0`` as falsy and silently apply the default 200 instead."""
    src = tmp_path / "in.h5ad"
    _write_adata(src)
    out = tmp_path / "run"
    rc = main([str(src), "--out", str(out),
              "--min-genes-per-cell", "0", "--min-cells-per-gene", "0"])
    assert rc == 0
    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["args"]["thresholds"]["min_genes_per_cell"] == 0
    assert manifest["args"]["thresholds"]["min_cells_per_gene"] == 0


def test_nonzero_threshold_still_overrides_default(tmp_path):
    src = tmp_path / "in.h5ad"
    _write_adata(src)
    out = tmp_path / "run"
    rc = main([str(src), "--out", str(out), "--min-genes-per-cell", "5"])
    assert rc == 0
    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["args"]["thresholds"]["min_genes_per_cell"] == 5
    assert manifest["args"]["thresholds"]["min_cells_per_gene"] == 3  # untouched default


def test_missing_data_or_out_is_a_usage_error(capsys):
    assert main(["--out", "x"]) == 1
    assert main([]) == 1


def test_list_strategies_exits_zero_and_lists_names(capsys):
    assert main(["--list-strategies"]) == 0
    out = capsys.readouterr().out
    assert "cpm_log1p" in out and "pseudobulk" in out


def test_qc_gate_failure_exits_two_and_force_overrides(tmp_path):
    # heavy mitochondrial fraction fails the mito_fraction gate on its own,
    # without also filtering every cell out via min_genes_per_cell - so
    # --force has a non-empty result to actually produce.
    rng = np.random.default_rng(0)
    n_cells, n_genes = 50, 400
    genes = [f"G{i}" for i in range(n_genes - 4)] + ["MT-A", "MT-B", "MT-C", "MT-D"]
    X = rng.poisson(3.0, size=(n_cells, n_genes)).astype(np.float32)
    X[:, -4:] = rng.poisson(300.0, size=(n_cells, 4))  # median ~50% mito -> fails mito_fail=25%
    src = tmp_path / "in.h5ad"
    ad.AnnData(X=sp.csr_matrix(X), var=pd.DataFrame(index=genes)).write_h5ad(src)

    out1 = tmp_path / "run_fail"
    assert main([str(src), "--out", str(out1)]) == 2
    assert not (out1 / "manifest.json").exists()
    out2 = tmp_path / "run_forced"
    assert main([str(src), "--out", str(out2), "--force"]) == 0
    manifest = json.loads((out2 / "manifest.json").read_text())
    assert manifest["forced_past_gates"] is True


def test_max_doublet_rate_flag_reaches_thresholds(tmp_path):
    src = tmp_path / "in.h5ad"
    _write_adata(src)
    out = tmp_path / "run"
    assert main([str(src), "--out", str(out), "--max-doublet-rate", "0.5"]) == 0
    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["args"]["thresholds"]["doublet_rate_max"] == 0.5
    assert manifest["args"]["thresholds"]["min_genes_per_cell"] == 200  # untouched
