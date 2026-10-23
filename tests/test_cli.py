"""CLI tests: the umbrella command, the one-shot `run`, and eval-with-stability."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from grnbench.cli_eval import main as eval_main
from grnbench.cli_infer import main as infer_main
from grnbench.cli_main import main as grnbench_main
from grnbench.cli_run import main as run_main


@pytest.fixture()
def dataset(tmp_path, synthetic_adata):
    """A small on-disk dataset plus a TF list that really matches its genes."""
    src = tmp_path / "data.h5ad"
    synthetic_adata.write_h5ad(src)
    tfs = tmp_path / "tfs.txt"
    tfs.write_text("GENE100\nGENE200\nGENE300\n", encoding="utf-8")
    return src, tfs


def _gold_from_edges(edges_path: Path, out: Path, n: int = 20) -> Path:
    """Build a gold standard out of real inferred pairs, spread over the ranking."""
    edges = pd.read_csv(edges_path, sep="\t")
    picked = edges.iloc[:: max(len(edges) // n, 1)].head(n)
    out.write_text(
        "regulator\ttarget\n"
        + "".join(f"{s}\t{t}\n" for s, t in zip(picked["source"], picked["target"])),
        encoding="utf-8")
    return out


# --------------------------------------------------------------------------- #
# umbrella command
# --------------------------------------------------------------------------- #

def test_umbrella_help_and_version(capsys):
    assert grnbench_main([]) == 0
    assert "grnbench <command>" in capsys.readouterr().out
    assert grnbench_main(["version"]) == 0
    import grnbench
    assert grnbench.__version__ in capsys.readouterr().out


def test_umbrella_lists_methods_and_strategies(capsys):
    assert grnbench_main(["methods"]) == 0
    out = capsys.readouterr().out
    assert "spearman_signed" in out and "builtin" in out
    assert grnbench_main(["strategies"]) == 0
    out = capsys.readouterr().out
    assert "cpm_log1p" in out and "scran_pooled" in out


def test_umbrella_rejects_unknown_command(capsys):
    assert grnbench_main(["frobnicate"]) == 1
    assert "unknown command" in capsys.readouterr().err


def test_umbrella_forwards_to_stage_cli(capsys):
    assert grnbench_main(["infer", "--methods", "list"]) == 0
    out = capsys.readouterr().out
    assert "spearman_signed" in out and "aracne3" in out


def test_methods_list_reports_unavailable_reasons(capsys):
    assert infer_main(["--methods", "list"]) == 0
    out = capsys.readouterr().out
    aracne = next(l for l in out.splitlines() if l.startswith("aracne3"))
    # availability is data: either it runs here, or the exact missing piece is named
    assert "available" in aracne.lower()
    if "UNAVAILABLE" in aracne:
        assert "java" in aracne or "ARACNE3_JAR" in aracne


# --------------------------------------------------------------------------- #
# grnbench run
# --------------------------------------------------------------------------- #

def test_run_without_gold_stops_after_inference(tmp_path, dataset, capsys):
    src, tfs = dataset
    out = tmp_path / "run1"
    rc = run_main(["--data", str(src), "--out", str(out), "--tf-list", str(tfs),
                   "--n-hvg", "120", "--methods", "spearman_signed", "--seed", "3"])
    assert rc == 0
    assert (out / "prep" / "processed.h5ad").exists()
    assert (out / "prep" / "report.md").exists()
    assert (out / "infer" / "spearman_signed.edges.tsv.gz").exists()
    assert (out / "infer" / "inference_manifest.json").exists()
    assert not (out / "eval").exists(), "eval must be skipped without --gold"

    manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "ok (eval skipped)"
    diagnostics = json.loads((out / "diagnostics.json").read_text(encoding="utf-8"))
    assert diagnostics["prep"]["default"]["all_cells"]["outcome"] == "ok"
    assert diagnostics["models"]["spearman_signed"]["all_cells"]["outcome"] == "ok"
    eval_stage = next(s for s in manifest["stages"] if s["stage"] == "eval")
    assert "no --gold given" in eval_stage["skipped"]
    assert "no --gold" in capsys.readouterr().out


def test_run_full_loop_with_gold_and_stability(tmp_path, dataset):
    src, tfs = dataset
    out1 = tmp_path / "stage"
    assert run_main(["--data", str(src), "--out", str(out1), "--tf-list", str(tfs),
                     "--n-hvg", "120", "--methods", "spearman_signed"]) == 0
    gold = _gold_from_edges(out1 / "infer" / "spearman_signed.edges.tsv.gz",
                            tmp_path / "gold.tsv")

    out = tmp_path / "run2"
    rc = run_main(["--data", str(src), "--out", str(out), "--tf-list", str(tfs),
                   "--n-hvg", "120", "--methods", "spearman_signed",
                   "--gold", str(gold), "--n-nulls", "10",
                   "--stability", "--stability-n", "3", "--stability-top-k", "50",
                   "--threads", "2"])
    assert rc == 0

    for artifact in ("prep/manifest.json", "prep/report.md",
                     "infer/inference_manifest.json",
                     "eval/eval_report.json", "eval/eval_report.md",
                     "eval/eval_report.html", "run_manifest.json"):
        assert (out / artifact).exists(), f"missing {artifact}"

    report = json.loads((out / "eval" / "eval_report.json").read_text(encoding="utf-8"))
    m = report["methods"]["spearman_signed"]
    assert m["n_gold_in_list"] > 0, "gold built from real edges must be reachable"
    assert m["auprc"] > 0 and m["null"]["n_nulls_used"] >= 1
    assert m["stability"]["n_subsamples"] == 3
    assert report["stability_config"]["threads"] == 2  # --threads reaches the stability re-runs
    assert 0.0 <= m["stability"]["jaccard_mean"] <= 1.0
    assert report["overlap"]["k"] == 1000
    assert "pr_curve" in m and m["pr_curve"]["recall"]

    html = (out / "eval" / "eval_report.html").read_text(encoding="utf-8")
    assert "<svg" in html and "spearman_signed" in html
    # the prep manifest is threaded through, so QC gates appear in the HTML
    assert "mito_fraction" in html

    run_manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    assert run_manifest["status"] == "ok"
    assert [s["stage"] for s in run_manifest["stages"]] == ["prep", "infer", "eval"]


def test_run_aborts_when_qc_gates_fail(tmp_path, extreme_mito_adata, capsys):
    src = tmp_path / "bad.h5ad"
    extreme_mito_adata.write_h5ad(src)
    out = tmp_path / "run_bad"
    rc = run_main(["--data", str(src), "--out", str(out)])
    assert rc == 2, "QC gate failure must propagate as exit code 2"
    assert not (out / "infer").exists(), "nothing downstream may run after an abort"
    manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "prep failed"
    diagnostics = json.loads((out / "diagnostics.json").read_text(encoding="utf-8"))
    assert diagnostics["outcome"] == "failed"
    assert diagnostics["status"] == "prep failed"
    assert diagnostics["prep"]["default"]["all_cells"]["failure"]["category"] == "qc_gate"
    assert "ABORT" in capsys.readouterr().err


def test_run_force_bypasses_qc_but_empty_preprocessing_result_still_fails(tmp_path, extreme_mito_adata):
    """``--force`` bypasses the mito QC gate, but the 10-gene fixture then has
    every cell filtered out by the (force-independent) min_genes_per_cell
    filter. Preprocessing itself now refuses that empty result loudly instead
    of writing a 0-cell processed.h5ad for inference to fail on later."""
    src = tmp_path / "bad.h5ad"
    extreme_mito_adata.write_h5ad(src)
    out = tmp_path / "run_forced"
    rc = run_main(["--data", str(src), "--out", str(out),
                   "--methods", "spearman_signed", "--force"])
    assert rc != 0, '--force cannot make a zero-gene matrix usable for inference'
    assert not (out / "prep" / "manifest.json").exists(), \
        "preprocessing raised before writing outputs; no partial manifest"
    run_manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    assert run_manifest["status"] == "prep failed"
    assert run_manifest["stages"][-1] == {"stage": "prep", "exit_code": rc, "out_dir": str(out / "prep")}


# --------------------------------------------------------------------------- #
# grnbench-eval standalone
# --------------------------------------------------------------------------- #

def test_eval_standalone_with_stability_and_plugin(tmp_path, dataset):
    src, tfs = dataset
    stage = tmp_path / "stage"
    assert run_main(["--data", str(src), "--out", str(stage), "--tf-list", str(tfs),
                     "--n-hvg", "120", "--methods", "spearman_signed"]) == 0
    edges = stage / "infer" / "spearman_signed.edges.tsv.gz"
    gold = _gold_from_edges(edges, tmp_path / "gold.tsv")

    out = tmp_path / "ev"
    rc = eval_main([str(edges), "--gold", str(gold), "--out", str(out),
                    "--n-nulls", "10", "--overlap-k", "200",
                    "--stability-data", str(stage / "prep" / "processed.h5ad"),
                    "--stability-n", "2", "--stability-top-k", "30",
                    "--stability-tf-list", str(tfs),
                    "--prep-manifest", str(stage / "prep" / "manifest.json")])
    assert rc == 0
    report = json.loads((out / "eval_report.json").read_text(encoding="utf-8"))
    assert report["stability_config"]["n_subsamples"] == 2
    assert report["methods"]["spearman_signed"]["stability"]["n_subsamples"] == 2
    md = (out / "eval_report.md").read_text(encoding="utf-8")
    assert "stability" in md and "Cross-method overlap" not in md  # single method


def test_eval_skips_unscorable_files_but_keeps_going(tmp_path, dataset, capsys):
    src, tfs = dataset
    stage = tmp_path / "stage"
    assert run_main(["--data", str(src), "--out", str(stage), "--tf-list", str(tfs),
                     "--n-hvg", "120", "--methods", "spearman_signed"]) == 0
    edges = stage / "infer" / "spearman_signed.edges.tsv.gz"
    gold = _gold_from_edges(edges, tmp_path / "gold.tsv")
    junk = tmp_path / "junk.edges.tsv"
    junk.write_text("a\tb\n1\t2\n", encoding="utf-8")

    out = tmp_path / "ev"
    rc = eval_main([str(junk), str(edges), "--gold", str(gold), "--out", str(out),
                    "--skip-nulls", "--no-html"])
    assert rc == 1
    assert "[SKIP] junk.edges.tsv" in capsys.readouterr().err
    report = json.loads((out / "eval_report.json").read_text(encoding="utf-8"))
    assert list(report["methods"]) == ["spearman_signed"]
    assert report['status'] == 'partial-failure'
    assert report['failures'][0]['stage'] == 'score'
    assert not (out / "eval_report.html").exists()


def test_misspelled_method_is_a_user_error_not_an_unavailability(tmp_path, dataset, capsys):
    """A typo must not be reported like a missing dependency.

    Regression: `--methods geni3` printed "[SKIP] geni3: unknown" - the same
    channel used for "arboreto not installed" - then crashed writing the
    manifest into a directory nothing had created.
    """
    src, _ = dataset
    rc = infer_main([str(src), "--out", str(tmp_path / "out"), "--methods", "geni3"])
    err = capsys.readouterr().err
    assert rc == 1, "an unknown method name must fail, not skip"
    assert "[UNKNOWN METHOD]" in err
    assert "genie3" in err, "the close match should be suggested"
    assert "registered methods:" in err
    assert not (tmp_path / "out" / "inference_manifest.json").exists()


def test_no_runnable_method_fails_cleanly_with_a_manifest(tmp_path, dataset, capsys):
    """Every requested method existing but unavailable is a distinct outcome."""
    src, _ = dataset
    table = {m["name"]: m for m in __import__(
        "grnbench.infer", fromlist=["list_methods"]).list_methods()}
    unavailable = [n for n, m in table.items() if not m["available"]]
    if not unavailable:
        pytest.skip("every registered method is available in this environment")
    rc = infer_main([str(src), "--out", str(tmp_path / "out2"),
                     "--methods", unavailable[0]])
    err = capsys.readouterr().err
    assert rc == 2, "a known-but-unavailable method is not the same as a typo"
    assert "[SKIP]" in err and "no requested method can run" in err


def test_eval_never_silently_drops_a_colliding_method_name(tmp_path, capsys):
    """Scoring the same method from two runs must score BOTH edge lists.

    Regression: the name was `path.name.split(".")[0]`, so
    runs/a/genie3.edges.tsv.gz and runs/b/genie3.edges.tsv.gz both derived
    "genie3" and the second silently overwrote the first - two files in, one
    method out, no warning.
    """
    from grnbench.infer import save_edges

    edges_a = pd.DataFrame([("TF1", "G1", 3.0), ("TF1", "G2", 2.0)],
                           columns=["source", "target", "weight"])
    edges_b = pd.DataFrame([("TF1", "G2", 3.0), ("TF1", "G1", 1.0)],
                           columns=["source", "target", "weight"])
    a, b = tmp_path / "runA", tmp_path / "runB"
    save_edges(edges_a, a, "genie3")
    save_edges(edges_b, b, "genie3")
    gold = tmp_path / "gold.tsv"
    gold.write_text("regulator\ttarget\nTF1\tG1\n", encoding="utf-8")

    out = tmp_path / "ev"
    rc = eval_main([str(a / "genie3.edges.tsv.gz"), str(b / "genie3.edges.tsv.gz"),
                    "--gold", str(gold), "--out", str(out),
                    "--skip-nulls", "--no-html"])
    assert rc == 0
    report = json.loads((out / "eval_report.json").read_text(encoding="utf-8"))
    assert len(report["methods"]) == 2, f"an edge list was dropped: {list(report['methods'])}"
    assert "genie3" in report["methods"]
    assert "runB/genie3" in report["methods"], list(report["methods"])
    assert report["name_collisions"] == [
        {"file": str(b / "genie3.edges.tsv.gz"), "derived": "genie3",
         "renamed_to": "runB/genie3"}]
    assert "[NAME COLLISION]" in capsys.readouterr().err
    # both are real, distinct rankings - the overlap matrix must see two methods
    assert len(report["overlap"]["methods"]) == 2


def test_eval_fails_loudly_on_a_bad_gold_spec(tmp_path, capsys):
    rc = eval_main(["--gold", "not_a_real_standard", "--out", str(tmp_path / "x")])
    assert rc == 1
    assert "[GOLD ERROR]" in capsys.readouterr().err


def test_eval_reports_a_bad_plugin_instead_of_crashing(tmp_path, capsys):
    rc = eval_main(["--gold", "whatever", "--out", str(tmp_path / "x"),
                    "--plugin", "./no_such_plugin.json"])
    assert rc == 1
    assert "[PLUGIN ERROR]" in capsys.readouterr().err


@pytest.mark.parametrize("module, required", [
    ("grnbench.cli_infer", []),
    ("grnbench.cli_run", ["--data", "d.h5ad", "--out", "out"]),
    ("grnbench.cli_eval", ["--gold", "gold.tsv"]),
])
@pytest.mark.parametrize("value", ["0", "-2"])
def test_threads_below_one_is_an_argument_error(module, required, value, capsys):
    """0 or a negative --threads used to reach sklearn/dask and fail there."""
    import importlib

    parser = importlib.import_module(module)._build_parser()
    with pytest.raises(SystemExit) as exit_:
        parser.parse_args([*required, "--threads", value])
    assert exit_.value.code == 2
    assert "--threads" in capsys.readouterr().err


def test_eval_stability_forwards_threads_to_subsample_stability(tmp_path, dataset, monkeypatch):
    """stability_config records --threads; this pins that the re-runs also receive it."""
    import grnbench.eval.stability as stability
    from grnbench.cli_eval import _build_parser, _run_stability

    seen = {}
    monkeypatch.setattr(stability, "subsample_stability",
                        lambda adata, name, **kwargs: seen.update(kwargs) or {})
    src, _ = dataset
    args = _build_parser().parse_args(["--gold", "gold.tsv", "--stability-data", str(src),
                                       "--threads", "3"])
    report = {"methods": {"spearman_signed": {}}}

    _run_stability(args, report, ["spearman_signed"])

    assert seen["threads"] == 3
# --------------------------------------------------------------------------- #
# defects found by the end-to-end test of 2026-10-09 (issue #74)
# --------------------------------------------------------------------------- #

def test_demo_never_overwrites_files_that_are_already_there(tmp_path, capsys):
    mine = tmp_path / "tfs.txt"
    mine.write_text("mine")
    assert grnbench_main(["demo", "--out", str(tmp_path)]) == 1
    assert mine.read_text() == "mine" and "must be empty" in capsys.readouterr().err
    assert grnbench_main(["demo", "--out", str(tmp_path / "fresh")]) == 0
    assert (tmp_path / "fresh" / "counts.h5ad").exists()


def test_demo_refuses_a_folder_holding_only_a_failed_diagnostic(tmp_path, capsys):
    """LD's review of #87: demo shared the retry exception of run outputs, wrote its files next to a failed
    prep's diagnostics.json and left that stale diagnostic in the demo folder."""
    from grnbench.cli import main as prep_main
    from grnbench.report import require_empty_output
    bad = tmp_path / "bad_mtx"
    bad.mkdir()
    (bad / "barcodes.tsv").touch()
    out = tmp_path / "prep_out"
    assert prep_main([str(bad), "--out", str(out)]) != 0
    assert [p.name for p in out.iterdir()] == ["diagnostics.json"]
    capsys.readouterr()
    assert grnbench_main(["demo", "--out", str(out)]) == 1
    assert "must be empty" in capsys.readouterr().err
    assert [p.name for p in out.iterdir()] == ["diagnostics.json"], "demo wrote into a non-empty folder"
    assert require_empty_output(out) == out, "a retry of the failed run still replaces its lone diagnostic"


def test_kb_refuses_unknown_options_and_stray_arguments(tmp_path, capsys, monkeypatch):
    """LD's review of #87: `kb definitely-not-an-option` and `kb --write-docs PATH extra` exited 0."""
    from grnbench import kb
    target = tmp_path / "kb.md"
    for bad, word in ((["definitely-not-an-option"], "definitely-not-an-option"), (["--bogus"], "--bogus"),
                      (["--write-docs", str(target), "extra"], "extra")):
        assert grnbench_main(["kb", *bad]) == 2, bad
        err = capsys.readouterr().err
        assert "usage: grnbench kb" in err and f"unrecognized arguments: {word}" in err, err
    assert not target.exists()
    assert grnbench_main(["kb", "--write-docs", str(target)]) == 0
    assert target.read_text(encoding="utf-8") == kb.render_markdown()
    monkeypatch.chdir(tmp_path)
    assert grnbench_main(["kb", "--write-docs"]) == 0, "PATH stays optional"
    assert (tmp_path / "docs" / "knowledge" / "preprocessing_kb.md").read_text(encoding="utf-8") == kb.render_markdown()
    capsys.readouterr()
    for unwritable in (str(tmp_path / "docs"), ""):    # a folder, an empty path: an error, not a traceback
        assert grnbench_main(["kb", "--write-docs", unwritable]) == 1, unwritable
        assert "grnbench kb: cannot write" in capsys.readouterr().err
    assert grnbench_main(["kb"]) == 0 and kb.option_table()[0]["option"] in capsys.readouterr().out


def test_builtin_commands_answer_help_and_refuse_stray_arguments(capsys, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    for cmd in ("methods", "strategies", "kb", "env", "demo", "ui", "version"):
        assert grnbench_main([cmd, "--help"]) == 0, cmd
        assert cmd in capsys.readouterr().out
    assert list(tmp_path.iterdir()) == [], "--help must not create anything (it used to write a folder named --help)"
    assert grnbench_main(["methods", "--bogus"]) == 2 and "takes no arguments" in capsys.readouterr().err
    assert grnbench_main(["demo", "--out"]) == 2 and "expected one argument" in capsys.readouterr().err


def test_env_reports_the_optional_tools_and_fails_only_for_an_r_you_named(capsys, monkeypatch, tmp_path):
    monkeypatch.setenv("GRNBENCH_RSCRIPT", "none")
    assert grnbench_main(["env", "--json"]) == 0, "R turned off is an answer, not a failure"
    rows = json.loads(capsys.readouterr().out)
    assert rows[0]["tool"] == "R" and rows[0]["found"] is False and "GRNBENCH_RSCRIPT=none" in rows[0]["install"]
    assert any(r["tool"] == "r:SoupX" and "ambient_rna=soupx" in r["used_for"] for r in rows)
    assert grnbench_main(["env"]) == 0
    out = capsys.readouterr().out
    assert out.splitlines()[0].split() == ["tool", "found", "version", "source"] and "r:scDblFinder" in out
    assert grnbench_main(["env", "--rscript", "no/such/Rscript"]) == 1, "the R you named does not work"
    assert "no Rscript at" in capsys.readouterr().out
    assert grnbench_main(["env", "--rscript", ""]) == 1 and "--rscript is empty" in capsys.readouterr().out
    monkeypatch.setenv("GRNBENCH_RSCRIPT", str(tmp_path / "Rscript"))
    assert grnbench_main(["env", "--json"]) == 1
    assert json.loads(capsys.readouterr().out)[0]["found"] is False, "--json still prints the rows"
    assert grnbench_main(["env", "--bogus"]) == 2 and "unrecognized arguments" in capsys.readouterr().err


def test_env_with_a_fake_r_lists_its_packages_and_asks_r_every_time(capsys, monkeypatch, fake_rscript):
    monkeypatch.setenv("GRNBENCH_RSCRIPT", str(fake_rscript.path))
    assert grnbench_main(["env"]) == 0
    out = capsys.readouterr().out
    soupx = next(line for line in out.splitlines() if line.startswith("r:SoupX "))
    assert soupx.split()[1:3] == ["yes", "1.6.2"]
    assert 'install: BiocManager::install("scDblFinder"' in out
    assert len(fake_rscript.calls()) == 1
    assert grnbench_main(["env"]) == 0
    assert len(fake_rscript.calls()) == 2, "a diagnostic does not trust the probe cache"


def test_ui_listens_on_this_computer_only_unless_asked(monkeypatch):
    pytest.importorskip("streamlit")
    seen = []
    monkeypatch.setattr("subprocess.call", lambda cmd, *a, **k: seen.append(cmd) or 0)
    assert grnbench_main(["ui"]) == 0 and grnbench_main(["ui", "--server.address", "0.0.0.0"]) == 0
    assert seen[0][-2:] != ["--server.address", "0.0.0.0"] and seen[0][seen[0].index("--server.address") + 1] == "localhost"
    assert seen[1].count("--server.address") == 1 and seen[1][seen[1].index("--server.address") + 1] == "0.0.0.0"


def test_python_dash_m_grnbench_works():
    import subprocess
    import sys
    out = subprocess.run([sys.executable, "-m", "grnbench", "version"], capture_output=True, text=True)
    assert out.returncode == 0 and "grnbench" in out.stdout


def test_an_unavailable_method_is_skipped_loudly_and_the_others_still_run_and_are_scored(tmp_path, dataset, capsys):
    """Naming a method that cannot run here (aracne3 needs java) used to abort `run` before evaluation."""
    src, tfs = dataset
    stage = tmp_path / "stage"
    assert run_main(["--data", str(src), "--out", str(stage), "--tf-list", str(tfs), "--n-hvg", "120",
                     "--methods", "spearman_signed"]) == 0
    gold = _gold_from_edges(stage / "infer" / "spearman_signed.edges.tsv.gz", tmp_path / "gold.tsv")
    out = tmp_path / "run"
    rc = run_main(["--data", str(src), "--out", str(out), "--tf-list", str(tfs), "--n-hvg", "120",
                   "--methods", "spearman_signed,aracne3", "--gold", str(gold), "--n-nulls", "10"])
    assert rc == 0
    assert "[SKIP] aracne3" in capsys.readouterr().err
    manifest = json.loads((out / "infer" / "inference_manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "partial" and manifest["errors"] == {} and "aracne3" in manifest["methods_skipped"]
    assert (out / "eval" / "eval_report.json").exists(), "the methods that ran are evaluated"


def test_nothing_runnable_is_still_an_error(tmp_path, dataset):
    src, tfs = dataset
    assert infer_main([str(src), "--out", str(tmp_path / "o"), "--methods", "aracne3"]) == 2


def test_infer_refuses_a_matrix_with_nan_and_says_why(tmp_path, dataset, capsys):
    """A NaN in X used to produce all-zero pearson/spearman networks and exit 0."""
    import anndata as ad
    import numpy as np
    src, tfs = dataset
    adata = ad.read_h5ad(src)
    X = np.asarray(adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X, dtype=float)
    X[0, 0] = np.nan
    adata.X = X
    bad = tmp_path / "with_nan.h5ad"
    adata.write_h5ad(bad)
    out = tmp_path / "o"
    assert infer_main([str(bad), "--out", str(out), "--tf-list", str(tfs), "--methods", "pearson_signed"]) == 1
    assert "non-finite" in capsys.readouterr().err
    manifest = json.loads((out / "inference_manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "failed" and "non-finite" in manifest["errors"]["setup"]
    assert not list(out.glob("*.edges.tsv.gz")), "no edge file may be written for a matrix the methods cannot use"


def test_infer_without_a_tf_list_says_every_gene_is_a_regulator(tmp_path, dataset, capsys):
    """Omitting --tf-list makes the problem genes x genes; nothing used to say so (HI-17/M2)."""
    import anndata as ad
    src, tfs = dataset
    n = ad.read_h5ad(src).n_vars
    out = tmp_path / "all_genes"
    assert infer_main([str(src), "--out", str(out), "--methods", "spearman_signed"]) == 0
    err = capsys.readouterr().err
    assert "[NO TF LIST]" in err and f"all {n:,} genes" in err and f"{n * (n - 1):,} candidate edges" in err
    manifest = json.loads((out / "inference_manifest.json").read_text(encoding="utf-8"))
    assert manifest["warnings"] and "no --tf-list" in manifest["warnings"][0]
    quiet = tmp_path / "with_tfs"
    assert infer_main([str(src), "--out", str(quiet), "--tf-list", str(tfs), "--methods", "spearman_signed"]) == 0
    assert "NO TF LIST" not in capsys.readouterr().err
    assert "warnings" not in json.loads((quiet / "inference_manifest.json").read_text(encoding="utf-8"))


def test_a_gold_standard_that_none_of_the_candidates_can_reach_is_a_failure_not_a_vacuous_success(
        tmp_path, dataset, capsys):
    """Empty, foreign or swapped reference files produced an all-NaN table, 'complete' and exit 0."""
    src, tfs = dataset
    stage = tmp_path / "stage"
    assert run_main(["--data", str(src), "--out", str(stage), "--tf-list", str(tfs), "--n-hvg", "120",
                     "--methods", "spearman_signed"]) == 0
    edges = stage / "infer" / "spearman_signed.edges.tsv.gz"
    real = _gold_from_edges(edges, tmp_path / "gold.tsv")

    foreign = tmp_path / "foreign.tsv"          # other gene names altogether
    foreign.write_text("regulator\ttarget\nSTAT1\tIRF1\nSTAT3\tMYC\n", encoding="utf-8")
    rc = eval_main([str(edges), "--gold", str(foreign), "--out", str(tmp_path / "e1"), "--skip-nulls"])
    err = capsys.readouterr().err
    assert rc == 1 and "[GOLD EMPTY]" in err and "vacuous" in err and "swapped" not in err
    report = json.loads((tmp_path / "e1" / "eval_report.json").read_text(encoding="utf-8"))
    assert report["status"] == "partial-failure"
    assert [f["stage"] for f in report["failures"]] == ["gold_overlap"]

    swapped = tmp_path / "swapped.tsv"          # the real pairs with the columns exchanged
    rows = [line.split("\t") for line in real.read_text(encoding="utf-8").splitlines()[1:]]
    swapped.write_text("regulator\ttarget\n" + "".join(f"{b}\t{a}\n" for a, b, *_ in rows), encoding="utf-8")
    rc = eval_main([str(edges), "--gold", str(swapped), "--out", str(tmp_path / "e2"), "--skip-nulls"])
    assert rc == 1 and "columns look swapped" in capsys.readouterr().err

    assert eval_main([str(edges), "--gold", str(real), "--out", str(tmp_path / "e3"), "--skip-nulls"]) == 0
    assert "[GOLD EMPTY]" not in capsys.readouterr().err


def test_run_checks_the_gold_standard_before_it_does_any_work(tmp_path, dataset, capsys):
    """A typo in --gold used to be found after the preprocessing and the inference were done (HI-15)."""
    src, tfs = dataset
    out = tmp_path / "run"
    assert run_main(["--data", str(src), "--out", str(out), "--tf-list", str(tfs), "--methods", "spearman_signed",
                     "--gold", str(tmp_path / "no_such_gold.tsv")]) == 1
    assert "[GOLD ERROR]" in capsys.readouterr().err
    assert not out.exists(), "nothing may have been written"


def test_prep_on_relative_mtx_folder_writes_diagnostics_without_warnings(tmp_path, monkeypatch):
    import warnings
    from grnbench.cli import main as prep_main
    monkeypatch.chdir(tmp_path)
    rel_mtx = Path("rel_mtx_dir")
    rel_mtx.mkdir()
    (rel_mtx / "barcodes.tsv").touch()
    out = tmp_path / "prep_out"

    with warnings.catch_warnings(record=True) as recorded_warnings:
        warnings.simplefilter("always")
        rc = prep_main([str(rel_mtx), "--out", str(out)])
    assert rc != 0
    assert (out / "diagnostics.json").exists()
    diag = json.loads((out / "diagnostics.json").read_text(encoding="utf-8"))
    assert diag["outcome"] == "failed"
    assert not [w for w in recorded_warnings if issubclass(w.category, RuntimeWarning) and "portable diagnostics" in str(w.message)]
