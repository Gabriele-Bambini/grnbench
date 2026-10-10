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
    diagnostics = json.loads((out / "diagnostics.json").read_text())
    assert manifest["args"]["thresholds"]["min_genes_per_cell"] == 0
    assert manifest["args"]["thresholds"]["min_cells_per_gene"] == 0
    assert diagnostics["prep"]["default"]["all_cells"]["outcome"] == "ok"
    choices = diagnostics["prep"]["default"]["all_cells"]["choices"]
    assert choices["requested_strategy"] == "auto"
    assert choices["strategy"] == manifest["strategy_chosen"]


def test_nonzero_threshold_still_overrides_default(tmp_path):
    src = tmp_path / "in.h5ad"
    _write_adata(src)
    out = tmp_path / "run"
    rc = main([str(src), "--out", str(out), "--min-genes-per-cell", "5"])
    assert rc == 0
    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["args"]["thresholds"]["min_genes_per_cell"] == 5
    assert manifest["args"]["thresholds"]["min_cells_per_gene"] == 3  # untouched default


def test_subset_run_uses_its_real_diagnostic_context(tmp_path):
    src = tmp_path / "in.h5ad"
    rng = np.random.default_rng(1)
    X = rng.poisson(3.0, size=(100, 400)).astype(np.float32)
    obs = pd.DataFrame({"cell_type": ["alpha"] * 60 + ["beta"] * 40})
    var = pd.DataFrame(index=[f"G{i}" for i in range(400)])
    ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var).write_h5ad(src)
    out = tmp_path / "subset"

    assert main([str(src), "--out", str(out), "--subset-key", "cell_type",
                 "--subset-value", "alpha"]) == 0
    diagnostics = json.loads((out / "diagnostics.json").read_text())
    assert set(diagnostics["prep"]["default"]) == {"cell_type=alpha"}


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
    failed = json.loads((out1 / "diagnostics.json").read_text())
    assert failed["outcome"] == "failed"
    assert failed["prep"]["default"]["all_cells"]["failure"]["category"] == "qc_gate"
    # The CLI explicitly recommends this retry.  A lone failure diagnostic is
    # safe to replace and must not make the recommended command impossible.
    assert main([str(src), "--out", str(out1), "--force"]) == 0
    manifest = json.loads((out1 / "manifest.json").read_text())
    assert manifest["forced_past_gates"] is True
    retried = json.loads((out1 / "diagnostics.json").read_text())
    assert retried["outcome"] == "ok"


def test_failed_subset_prep_uses_requested_context(tmp_path):
    src = tmp_path / "in.h5ad"
    _write_adata(src)
    out = tmp_path / "missing_subset"
    assert main([str(src), "--out", str(out), "--subset-key", "cell_type",
                 "--subset-value", "missing"]) == 2
    diagnostic = json.loads((out / "diagnostics.json").read_text())
    assert set(diagnostic["prep"]["default"]) == {"cell_type=missing"}


def test_nonempty_output_is_refused_without_adding_a_failure_diagnostic(tmp_path):
    src = tmp_path / "in.h5ad"
    _write_adata(src)
    out = tmp_path / "old_run"
    out.mkdir()
    old = out / "manifest.json"
    old.write_text('{"old": true}\n', encoding="utf-8")
    assert main([str(src), "--out", str(out)]) == 1
    assert old.read_text(encoding="utf-8") == '{"old": true}\n'
    assert not (out / "diagnostics.json").exists()


def test_real_missing_file_error_does_not_leak_quoted_windows_path(tmp_path):
    private = r"Z:\restricted\secret\absent.h5ad"
    out = tmp_path / "missing"
    assert main([private, "--out", str(out)]) == 1
    text = (out / "diagnostics.json").read_text(encoding="utf-8")
    assert private not in text
    assert "restricted" not in text
    assert "<path>" in text


def test_diagnostic_bug_warns_but_does_not_fail_successful_prep(tmp_path, monkeypatch):
    import grnbench.diagnostics as diagnostics

    src = tmp_path / "in.h5ad"
    _write_adata(src)
    out = tmp_path / "run"

    def broken_writer(*args, **kwargs):
        raise RuntimeError("diagnostic bug")

    monkeypatch.setattr(diagnostics, "write_diagnostics", broken_writer)
    with pytest.warns(RuntimeWarning, match="portable diagnostics unavailable"):
        assert main([str(src), "--out", str(out),
                     "--min-genes-per-cell", "0", "--min-cells-per-gene", "0"]) == 0
    assert (out / "manifest.json").exists()
    assert (out / "report.md").exists()
    assert not (out / "diagnostics.json").exists()


def test_max_doublet_rate_flag_reaches_thresholds(tmp_path):
    src = tmp_path / "in.h5ad"
    _write_adata(src)
    out = tmp_path / "run"
    assert main([str(src), "--out", str(out), "--max-doublet-rate", "0.5"]) == 0
    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["args"]["thresholds"]["doublet_rate_max"] == 0.5
    assert manifest["args"]["thresholds"]["min_genes_per_cell"] == 200  # untouched
