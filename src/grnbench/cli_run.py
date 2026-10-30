"""grnbench run — the whole loop in one command.

    grnbench run --data d.h5ad --out runs/2026-09-05 --tf-list tfs.txt \\
                 --methods genie3,grnboost2,my-model --plugin my_model.json \\
                 --gold dorotheaABC --n-hvg 3000

Stages, each writing into its own subdirectory of ``--out``:

    prep/   processed.h5ad, qc_metrics.json, manifest.json, report.md
    infer/  <method>.edges.tsv.gz, inference_manifest.json
    eval/   eval_report.{json,md,html}

Every stage is the *same code* the standalone CLIs run, so a one-command run and
three separate runs produce identical artifacts. Stages are skipped explicitly
(never silently): no ``--gold`` means no evaluation, and the reason is printed
and recorded in ``run_manifest.json``.
"""
from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

from grnbench.report import thread_count


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="grnbench run",
        description="Preprocess, infer and evaluate in one command.")
    p.add_argument("--data", required=True,
                   help="input dataset (.h5ad/.csv/.tsv/.loom/10x dir/10x .h5)")
    p.add_argument("--out", required=True, help="run directory (created)")
    p.add_argument("--tf-list", default=None, help="TF symbol list (.txt/.csv/.tsv)")

    from grnbench.cli import N_HVG_HELP, n_hvg_value

    prep = p.add_argument_group("preprocessing")
    prep.add_argument("--strategy", default="auto")
    prep.add_argument("--n-hvg", type=n_hvg_value, default=None, metavar="all|N", help=N_HVG_HELP)
    prep.add_argument("--group-key", default=None)
    prep.add_argument("--label-key", default=None)
    prep.add_argument("--label-mapping", default=None)
    prep.add_argument("--gene-mapping", default=None,
                      help="two-column id -> symbol table for preprocessing")
    prep.add_argument("--stratify-key", default=None)
    prep.add_argument("--batch-key", default=None)
    prep.add_argument("--min-group-cells", type=int, default=None)
    prep.add_argument("--subsample-cells", type=int, default=None)
    prep.add_argument("--subsample-mode", choices=["random", "stratified"], default=None)
    prep.add_argument("--orientation", default="auto",
                      choices=["auto", "genes_x_cells", "cells_x_genes"],
                      help="CSV/TSV layout; 'auto' refuses to guess "
                           "when the axes are ambiguous")
    prep.add_argument("--no-downcast", action="store_true")
    prep.add_argument("--drop-doublets", action="store_true")
    prep.add_argument("--max-doublet-rate", type=float, default=None)
    prep.add_argument("--min-genes-per-cell", type=int, default=None)
    prep.add_argument("--min-cells-per-gene", type=int, default=None)
    prep.add_argument("--force", action="store_true",
                      help="proceed past QC gate failures (recorded in the manifest)")

    inf = p.add_argument_group("inference")
    inf.add_argument("--methods", default="available",
                     help="comma-separated names, 'available' (default set) or 'all'")
    inf.add_argument("--plugin", action="append", default=None, metavar="SPEC",
                     help="custom model spec file or 'module:Attribute' (repeatable)")
    inf.add_argument("--threads", type=thread_count, default=4,
                     help="CPU budget (>= 1) for random_forest, extra_trees, GENIE3, GRNBoost2 and ARACNe3; default 4, not every core; also used by the stability re-runs")

    ev = p.add_argument_group("evaluation")
    ev.add_argument("--gold", default=None,
                    help="gold pairs file or dorothea<levels> spec; omit to skip eval")
    ev.add_argument("--top-k", type=int, default=None)
    ev.add_argument("--n-nulls", type=int, default=50)
    ev.add_argument("--null-mode", choices=["degree", "random", "regulator"], default="random")
    ev.add_argument("--skip-nulls", action="store_true")
    ev.add_argument("--overlap-k", type=int, default=1000)
    ev.add_argument("--stability", action="store_true",
                    help="re-run each method on cell subsamples (slow, worth it)")
    ev.add_argument("--stability-n", type=int, default=5)
    ev.add_argument("--stability-frac", type=float, default=0.8)
    ev.add_argument("--stability-top-k", type=int, default=100)
    ev.add_argument("--no-html", action="store_true")

    p.add_argument("--seed", type=int, default=42)
    return p


def _opt(flag: str, value) -> list[str]:
    return [] if value is None else [flag, str(value)]


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    if args.gold and not Path(args.gold).exists() and not str(args.gold).lower().startswith("dorothea"):
        # the evaluation is the last stage: a typo here used to surface after hours of preprocessing and inference
        print(f"[GOLD ERROR] gold standard '{args.gold}': not a file and not a dorothea* spec "
              "(checked before any work starts)", file=sys.stderr)
        return 1
    from grnbench.report import require_empty_output
    try:
        require_empty_output(args.out)
    except ValueError as e:
        print(str(e), file=sys.stderr)
        return 1
    out = Path(args.out)
    prep_dir, infer_dir, eval_dir = out / "prep", out / "infer", out / "eval"
    started = datetime.now(timezone.utc).isoformat(timespec="seconds")
    stages: list[dict] = []

    # -- 1. preprocess ------------------------------------------------------- #
    print("=" * 62)
    print("[1/3] preprocessing")
    print("=" * 62)
    from grnbench.cli import main as prep_main
    prep_argv = [
        args.data, "--out", str(prep_dir), "--strategy", args.strategy,
        "--seed", str(args.seed),
        *_opt("--tf-list", args.tf_list), *_opt("--n-hvg", args.n_hvg),
        *_opt("--group-key", args.group_key), *_opt("--label-key", args.label_key),
        *_opt("--label-mapping", args.label_mapping),
        *_opt("--gene-mapping", args.gene_mapping),
        *_opt("--stratify-key", args.stratify_key),
        *_opt("--batch-key", args.batch_key),
        *_opt("--min-group-cells", args.min_group_cells),
        *_opt("--subsample-cells", args.subsample_cells),
        *_opt("--subsample-mode", args.subsample_mode),
        *_opt("--orientation", args.orientation),
        *_opt("--min-genes-per-cell", args.min_genes_per_cell),
        *_opt("--min-cells-per-gene", args.min_cells_per_gene),
        *_opt("--max-doublet-rate", args.max_doublet_rate),
        *(["--no-downcast"] if args.no_downcast else []),
        *(["--drop-doublets"] if args.drop_doublets else []),
        *(["--force"] if args.force else []),
    ]
    rc = prep_main(prep_argv)
    stages.append({"stage": "prep", "exit_code": rc, "out_dir": str(prep_dir)})
    if rc != 0:
        _write_manifest(out, args, started, stages, "prep failed")
        print("\n[ABORT] preprocessing failed; nothing downstream was run.",
              file=sys.stderr)
        return rc

    processed = prep_dir / "processed.h5ad"

    # -- 2. infer ------------------------------------------------------------ #
    print("\n" + "=" * 62)
    print("[2/3] inference")
    print("=" * 62)
    from grnbench.cli_infer import main as infer_main
    infer_argv = [
        str(processed), "--out", str(infer_dir), "--methods", args.methods,
        "--seed", str(args.seed), "--threads", str(args.threads),
        # deliberately NOT --assume-normalized: prepare_adata inspects the matrix,
        # so cpm_log1p/pearson_residuals output passes through untouched while a
        # raw_counts strategy still gets normalized before tree/correlation methods
        *_opt("--tf-list", args.tf_list),
    ]
    for spec in args.plugin or []:
        infer_argv += ["--plugin", spec]
    rc = infer_main(infer_argv)
    stages.append({"stage": "infer", "exit_code": rc, "out_dir": str(infer_dir)})
    if rc != 0:
        _write_manifest(out, args, started, stages, "inference failed")
        print("\n[ABORT] inference failed; evaluation skipped.", file=sys.stderr)
        return rc

    edge_files = sorted(str(p) for p in infer_dir.glob("*.edges.tsv.gz"))
    if not edge_files:
        stages.append({"stage": "eval", "skipped": "no edge lists were produced"})
        _write_manifest(out, args, started, stages, "no edges")
        print("\n[ABORT] no method produced an edge list.", file=sys.stderr)
        return 1

    # -- 3. evaluate --------------------------------------------------------- #
    print("\n" + "=" * 62)
    print("[3/3] evaluation")
    print("=" * 62)
    if not args.gold:
        reason = "no --gold given: inference artifacts written, scoring skipped"
        stages.append({"stage": "eval", "skipped": reason})
        print(reason)
        _write_manifest(out, args, started, stages, "ok (eval skipped)")
        _summary(out, prep_dir, infer_dir, None, edge_files)
        return 0

    from grnbench.cli_eval import main as eval_main
    from grnbench.diagnostics import run_document, try_write_diagnostics
    diagnostic_path = try_write_diagnostics(
        out / "diagnostics.json", run_document, out, outcome="partial", stages=stages)
    eval_argv = [
        *edge_files, "--gold", args.gold, "--out", str(eval_dir),
        "--universe", str(infer_dir / 'candidate_universe.tsv.gz'),
        "--n-nulls", str(args.n_nulls), "--null-mode", args.null_mode,
        "--overlap-k", str(args.overlap_k), "--seed", str(args.seed),
        "--prep-manifest", str(prep_dir / "manifest.json"),
        "--defer-report-output",
        *(["--diagnostics", str(diagnostic_path)] if diagnostic_path else []),
        *_opt("--top-k", args.top_k),
        *(["--skip-nulls"] if args.skip_nulls else []),
        *(["--no-html"] if args.no_html else []),
    ]
    for spec in args.plugin or []:
        eval_argv += ["--plugin", spec]
    if args.stability:
        eval_argv += [
            "--stability-data", str(processed),
            "--stability-n", str(args.stability_n),
            "--stability-frac", str(args.stability_frac),
            "--stability-top-k", str(args.stability_top_k),
            "--threads", str(args.threads),
            *_opt("--stability-tf-list", args.tf_list),
        ]
    rc = eval_main(eval_argv)
    stages.append({"stage": "eval", "exit_code": rc, "out_dir": str(eval_dir)})
    _write_manifest(out, args, started, stages, "ok" if rc == 0 else "eval failed")
    outputs = _refresh_evaluation_reports(out, prep_dir, eval_dir, html=not args.no_html)
    if outputs:
        print((eval_dir / "eval_report.md").read_text(encoding="utf-8"))
        print("\noutputs: " + " , ".join(str(path) for path in outputs))
    if rc == 0:
        _summary(out, prep_dir, infer_dir, eval_dir, edge_files)
    return rc


def _write_manifest(out: Path, args, started: str, stages: list[dict],
                    status: str) -> None:
    from grnbench.report import environment_block, write_json
    out.mkdir(parents=True, exist_ok=True)
    write_json(out / "run_manifest.json", {
        "grnbench_run": True,
        "started_utc": started,
        "finished_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "status": status,
        "args": vars(args),
        "stages": stages,
        "environment": environment_block(),
    })
    from grnbench.diagnostics import run_document, try_write_diagnostics
    try_write_diagnostics(out / "diagnostics.json", run_document,
                          out, outcome=status, stages=stages)


def _summary(out: Path, prep_dir: Path, infer_dir: Path,
             eval_dir: Path | None, edge_files: list[str]) -> None:
    print("\n" + "=" * 62)
    print(f"run complete -> {out}")
    print("=" * 62)
    print(f"  prep    {prep_dir / 'report.md'}")
    print(f"  infer   {len(edge_files)} edge list(s) in {infer_dir}")
    if eval_dir is not None:
        html = eval_dir / "eval_report.html"
        print(f"  eval    {eval_dir / 'eval_report.md'}")
        if html.exists():
            print(f"  report  {html}")
    print(f"  manifest {out / 'run_manifest.json'}")


def _refresh_evaluation_reports(out: Path, prep_dir: Path, eval_dir: Path,
                                *, html: bool) -> list[Path]:
    """Regenerate static reports after the final root diagnostic exists.

    ``grnbench-eval`` finishes before the one-shot run knows its terminal
    outcome.  The score JSON is unchanged; only Markdown/HTML are rendered
    again from that JSON and the final portable diagnostic.
    """
    import json

    report_path = eval_dir / "eval_report.json"
    if not report_path.exists():
        return []
    from grnbench.cli_eval import write_evaluation_reports
    from grnbench.diagnostics import read_diagnostics
    from grnbench.report import file_fingerprint

    report = json.loads(report_path.read_text(encoding="utf-8"))
    prep_path = prep_dir / "manifest.json"
    prep = json.loads(prep_path.read_text(encoding="utf-8")) if prep_path.exists() else None
    diagnostic_path = out / "diagnostics.json"
    diagnostics = read_diagnostics(diagnostic_path)
    if diagnostics is not None:
        diagnostic_fingerprint = file_fingerprint(diagnostic_path)
        diagnostic_fingerprint["path"] = diagnostic_path.name
        report.setdefault("provenance", {}).setdefault("additional_inputs", {})[
            "diagnostics"] = diagnostic_fingerprint
    return write_evaluation_reports(report, eval_dir, prep_manifest=prep,
                                    diagnostics=diagnostics, html=html)


if __name__ == "__main__":
    raise SystemExit(main())
