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
    report = (out / "report.md").read_text(encoding="utf-8")
    assert "## Pipeline diagnostics" in report
    assert "default" in report and "all_cells" in report
    assert "100 cells x 400 genes" in report


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


@pytest.mark.parametrize("extra, message", [
    (["--subset-value", "alpha"], "--subset-value needs --subset-key"),
    (["--subset-key", "cell_type"], "--subset-key needs --subset-value"),
    (["--n-hvg", "-5"], "--n-hvg must be at least 1"),
    (["--n-hvg", "0"], "--n-hvg must be at least 1"),
    (["--subsample-cells", "0"], "--subsample-cells must be at least 1"),
    (["--subsample-mode", "random"], "--subsample-mode needs --subsample-cells"),
    (["--subsample-cells", "10", "--subsample-mode", "stratified"],
     "--subsample-mode stratified needs"),
    (["--cells-per-metacell", "0"], "--cells-per-metacell must be at least 1"),
    (["--min-genes-per-cell", "-3"], "cannot be negative"),
    (["--max-doublet-rate", "2"], "fraction between 0 and 1"),
    (["--mito-warn", "50", "--mito-fail", "10"], "above --mito-fail"),
    (["--mito-warn", "40"], "above --mito-fail (25)"),      # inverted against the default failure level
], ids=lambda v: v if isinstance(v, str) else " ".join(v))
def test_options_that_would_be_ignored_or_fail_deep_inside_scanpy_are_refused_first(tmp_path, capsys, extra, message):
    """HI-02/HI-13: `--subset-value` alone processed the whole dataset and recorded a subset it never applied."""
    src = tmp_path / "in.h5ad"
    _write_adata(src)
    out = tmp_path / "run"
    assert main([str(src), "--out", str(out)] + extra) == 1
    assert message in capsys.readouterr().err
    assert not out.exists(), "a refused command leaves no output folder behind"


def test_valid_option_combinations_are_still_accepted(tmp_path):
    src = tmp_path / "in.h5ad"
    _write_adata(src)
    assert main([str(src), "--out", str(tmp_path / "ok"), "--min-genes-per-cell", "0", "--mito-warn", "5",
                 "--mito-fail", "5", "--max-doublet-rate", "0", "--n-hvg", "100"]) == 0


def test_hvg_request_above_gene_count_says_all_genes_are_kept(tmp_path, capsys):
    """#164: --n-hvg 99999 was recorded, but the effective all-gene fallback was silent in the CLI."""
    src = tmp_path / "in.h5ad"
    _write_adata(src, n_genes=400)
    out = tmp_path / "run"

    assert main([str(src), "--out", str(out), "--n-hvg", "99999"]) == 0

    terminal = capsys.readouterr().out
    manifest = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    warning = next(w for w in manifest["warnings"] if w.startswith("Feature selection:"))
    assert manifest["args"]["n_hvg"] == 99999  # requested value remains provenance
    assert "requested 99,999 HVGs" in warning
    assert "kept all 400 genes" in warning
    assert f"  ! {warning}" in terminal
    assert warning in (out / "report.md").read_text(encoding="utf-8")


def test_explicit_random_subsample_mode_reaches_the_manifest(tmp_path):
    src = tmp_path / "in.h5ad"
    _write_adata(src, n_cells=100)
    out = tmp_path / "run"
    assert main([str(src), "--out", str(out), "--subsample-cells", "50",
                 "--subsample-mode", "random"]) == 0
    manifest = json.loads((out / "manifest.json").read_text())
    diagnostics = json.loads((out / "diagnostics.json").read_text())
    step = [s for s in manifest["steps"] if s["step"] == "subsample_cells"][0]
    assert manifest["args"]["subsample_mode"] == "random"
    assert step["mode"] == "random" and step["stratification"] is None
    assert diagnostics["prep"]["default"]["all_cells"]["choices"]["subsample_mode"] == "random"


def test_a_gene_filter_that_removes_every_gene_says_so_in_our_words(tmp_path, capsys):
    """`--min-cells-per-gene 100000` used to end in sklearn's 'Found array with 0 feature(s)'."""
    src = tmp_path / "in.h5ad"
    _write_adata(src)
    assert main([str(src), "--out", str(tmp_path / "run"), "--min-cells-per-gene", "100000"]) == 2
    err = capsys.readouterr().err
    assert "every gene was filtered out" in err and "--min-cells-per-gene" in err and "sklearn" not in err


def test_the_terminal_summary_repeats_the_manifest_warnings(tmp_path, capsys):
    """Low TF coverage was only in report.md/manifest.json; the person reading the terminal never saw it (HI-09)."""
    src = tmp_path / "in.h5ad"
    _write_adata(src)
    tfs = tmp_path / "tfs.txt"
    tfs.write_text("G1\n" + "".join(f"NOT_A_GENE_{i}\n" for i in range(40)), encoding="utf-8")
    out = tmp_path / "run"
    assert main([str(src), "--out", str(out), "--tf-list", str(tfs)]) == 0
    text = capsys.readouterr().out
    manifest = json.loads((out / "manifest.json").read_text())
    coverage = [w for w in manifest["warnings"] if "TF coverage low" in w]
    assert coverage and "\nWarnings:\n" in text and f"  ! {coverage[0]}" in text
    assert text.index("Warnings:") < text.index("Recommendations:")


def test_an_output_folder_named_like_a_placeholder_keeps_the_failure_message(tmp_path, monkeypatch):
    """`--out path` is passed to the redaction as a sensitive value; '<path>' used to nest into
    '<<path>>' until the whole message was withheld (#177 re-review)."""
    from grnbench.diagnostics import read_diagnostics

    monkeypatch.chdir(tmp_path)
    assert main(["missing.h5ad", "--out", "path"]) == 1
    failure = read_diagnostics(tmp_path / "path" / "diagnostics.json")["prep"]["default"]["all_cells"]["failure"]
    assert failure["category"] == "prep_error"
    assert "<path>" in failure["message"] and "<<" not in failure["message"]
    assert failure["message"] != "<message withheld>"