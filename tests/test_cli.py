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
                   "--stability", "--stability-n", "3", "--stability-top-k", "50"])
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
