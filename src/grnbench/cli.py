"""Command-line interface.

    grnbench-prep DATA --out DIR [--strategy auto] [--tf-list FILE]
                      [--n-hvg N] [--group-key COL] [--force]

Exit codes: 0 success - 2 QC gate failure (use --force to override) - 1 error.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from grnbench.hygiene import load_label_mapping
from grnbench.pipeline import PrepError, run_prep
from grnbench.qc import Thresholds
from grnbench.report import require_empty_output
from grnbench.strategies import list_strategies
from grnbench.tfs import load_tf_list


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="grnbench-prep",
        description="Bulletproof preprocessing for GRN-bound scRNA-seq datasets.",
    )
    p.add_argument("data", nargs="?", help="input dataset (.h5ad/.csv/.tsv/.loom/10x dir)")
    p.add_argument("--out", help="output directory for the run")
    p.add_argument("--strategy", default="auto",
                   help="auto (default) or one of the registered strategies")
    p.add_argument("--tf-list", default=None,
                   help="TF symbol list (.txt/.csv/.tsv); enables TF-union reporting")
    p.add_argument("--n-hvg", type=int, default=None,
                   help="select N HVGs unioned with TFs (default: keep all filtered genes)")
    p.add_argument("--group-key", default=None,
                   help="obs column used by pseudobulk and recommendations; cells with "
                        "a missing value are DROPPED (logged)")
    p.add_argument("--label-key", default=None,
                   help="obs column with class labels: hygiene + imbalance audit; cells "
                        "with a missing value are KEPT and excluded from stratification")
    p.add_argument("--label-mapping", default=None,
                   help="two-column file mapping raw -> clean labels")
    p.add_argument("--stratify-key", default=None,
                   help="obs column for stratified subsampling (default: group/label key)")
    p.add_argument("--batch-key", default=None,
                   help="batch/donor column (auto-detected if omitted); required by "
                        "--strategy cpm_log1p_combat")
    p.add_argument("--min-group-cells", type=int, default=None,
                   help="pseudobulk: drop groups smaller than this (loudly)")
    p.add_argument("--no-downcast", action="store_true",
                   help="keep float64 matrices as-is")
    p.add_argument("--min-genes-per-cell", type=int, default=None)
    p.add_argument("--min-cells-per-gene", type=int, default=None)
    p.add_argument("--max-genes-per-cell", type=int, default=None,
                   help="drop cells with more detected genes (doublet-like outliers)")
    p.add_argument("--max-counts-per-cell", type=float, default=None,
                   help="drop cells with a larger library")
    p.add_argument("--mito-warn", type=float, default=None,
                   help="dataset gate: warn when the median %% mitochondrial exceeds this (default 10)")
    p.add_argument("--mito-fail", type=float, default=None,
                   help="dataset gate: refuse when the median %% mitochondrial exceeds this (default 25)")
    p.add_argument("--max-pct-mito", type=float, default=None,
                   help="drop cells whose mitochondrial fraction (%%) is higher")
    p.add_argument("--max-doublet-rate", type=float, default=None,
                   help="--drop-doublets refuses to drop when the call rate exceeds this "
                        "fraction (default 0.20): a failed threshold, not a doublet population")
    p.add_argument("--subsample-cells", type=int, default=None)
    p.add_argument("--orientation", default="auto",
                   choices=["auto", "genes_x_cells", "cells_x_genes"],
                   help="CSV/TSV layout. 'auto' detects it and refuses to "
                        "guess when ambiguous; declare it when told to "
                        "(genes_x_cells = genes as rows, the GEO layout)")
    p.add_argument("--drop-doublets", action="store_true",
                   help="call doublets with scrublet (per batch when a batch column "
                        "exists) and drop them; scores kept in obs; never on by default")
    p.add_argument("--force", action="store_true",
                   help="proceed past QC gate failures (recorded in manifest)")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--subset-key", default=None,
                   help="obs column selecting one cellular context (with --subset-value)")
    p.add_argument("--subset-value", default=None, help="value of --subset-key to keep")
    p.add_argument("--cell-calling", default=None, choices=["knee"],
                   help="call cells in an unfiltered droplet matrix before QC")
    p.add_argument("--exclude-qc-genes", action="store_true",
                   help="drop mitochondrial/ribosomal genes from the features after QC")
    p.add_argument("--per-batch-bounds", default=None, metavar="JSON",
                   help="JSON file {batch: {min_genes_per_cell, max_genes_per_cell, "
                        "max_counts_per_cell, max_pct_mito_cell}} for per-batch cell QC")
    p.add_argument("--cells-per-metacell", type=int, default=None,
                   help="metacell size for --strategy metacells (default 20)")
    p.add_argument("--list-strategies", action="store_true")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)

    if args.list_strategies:
        for s in list_strategies():
            status = "available" if s["available"] else f"UNAVAILABLE ({s['unavailable_reason']})"
            print(f"{s['name']:20s} {status:45s} {s['description']}")
        return 0

    if not args.data or not args.out:
        print("error: both DATA and --out are required (or use --list-strategies)",
              file=sys.stderr)
        return 1

    # `x or default` would treat an explicit 0 as falsy and silently revert to
    # the default (200 / 3) - exactly the case a user passing 0 to disable a
    # filter needs honoured, not ignored. Only `None` means "not given".
    overrides = {k: v for k, v in {
        "min_genes_per_cell": args.min_genes_per_cell,
        "min_cells_per_gene": args.min_cells_per_gene,
        "doublet_rate_max": args.max_doublet_rate,
        "max_genes_per_cell": args.max_genes_per_cell,
        "max_counts_per_cell": args.max_counts_per_cell,
        "max_pct_mito_cell": args.max_pct_mito,
        "mito_warn": args.mito_warn,
        "mito_fail": args.mito_fail,
    }.items() if v is not None}
    thresholds = Thresholds(**overrides) if overrides else None

    output_claimed = False

    def write_failure(exc, category=None):
        """Write only into the empty directory created for this failed attempt."""
        if not output_claimed:
            return
        from grnbench.diagnostics import context_label, failed_prep_document, try_write_diagnostics
        try_write_diagnostics(
            Path(args.out) / "diagnostics.json",
            failed_prep_document,
            exc,
            choices=vars(args),
            context=context_label(subset_key=args.subset_key, subset_value=args.subset_value),
            category=category,
            sensitive_values=(args.data, args.out, args.tf_list, args.label_mapping,
                              args.per_batch_bounds),
        )

    try:
        # Claim the output before any work.  This makes it explicit whether a
        # later failure belongs to this attempt; a pre-existing run must never
        # gain a diagnostic from a command that was refused at the boundary.
        require_empty_output(args.out)
        output_claimed = True
        res = run_prep(
            args.data, args.out,
            strategy=args.strategy,
            tf_list=(load_tf_list(args.tf_list) if args.tf_list else None),
            n_hvg=args.n_hvg,
            group_key=args.group_key,
            label_key=args.label_key,
            label_mapping=(load_label_mapping(args.label_mapping)
                           if args.label_mapping else None),
            stratify_key=args.stratify_key,
            batch_key=args.batch_key,
            min_group_cells=args.min_group_cells,
            downcast=not args.no_downcast,
            thresholds=thresholds,
            subsample_cells=args.subsample_cells,
            orientation=args.orientation,
            drop_doublets=args.drop_doublets,
            force=args.force,
            seed=args.seed,
            subset_key=args.subset_key,
            subset_value=args.subset_value,
            exclude_qc_genes=args.exclude_qc_genes,
            cell_calling=args.cell_calling,
            per_batch_bounds=(json.loads(Path(args.per_batch_bounds).read_text(encoding="utf-8"))
                              if args.per_batch_bounds else None),
            cells_per_metacell=args.cells_per_metacell,
        )
    except PrepError as e:
        gate = str(e).startswith("QC gate failure")
        write_failure(e, "qc_gate" if gate else "prep_refused")
        print(f"[{'QC GATE FAILURE' if gate else 'PREP REFUSED'}] {e}", file=sys.stderr)
        if gate:
            print("(re-run with --force to proceed anyway)", file=sys.stderr)
        return 2
    except Exception as e:  # noqa: BLE001 - CLI boundary
        write_failure(f"{type(e).__name__}: {e}", "prep_error")
        print(f"[ERROR] {type(e).__name__}: {e}", file=sys.stderr)
        return 1

    print(f"\nGRN-Bench prep complete - strategy: {res.strategy}")
    print("\nQC gates:")
    for g in res.manifest["qc"]["gates"]:
        print(f"  [{g['status'].upper():4s}] {g['name']:24s} {g['detail']}")
    print("\nRecommendations:")
    for r in res.manifest["recommendations"]:
        marker = "*" if r["strategy"] == res.strategy else " "
        print(f" {marker} {r['rank']}. {r['strategy']:18s} {r['reason']}")
    print(f"\nOutputs:\n  {res.processed_path}\n  {res.manifest_path}\n"
          f"  {res.diagnostics_path}\n  {res.report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
