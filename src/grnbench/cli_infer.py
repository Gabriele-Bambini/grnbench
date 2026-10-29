"""grnbench-infer CLI: run GRN inference methods on a (preprocessed) dataset.

    grnbench-infer DATA --out DIR [--methods m1,m2] [--tf-list FILE]
                        [--assume-normalized] [--seed 42]

Default method set (``--methods available``) = every *available* registered method that
is in the default set; ``--methods all`` adds the methods registered later (marked
``[extra]`` by ``--methods list``). Methods with missing external dependencies are
skipped with an explicit notice, never silently.
Writes ``<method>.edges.tsv.gz`` per method plus ``inference_manifest.json``.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from grnbench.infer import (
    MethodError,
    list_methods,
    load_edges,
    run_method,
    save_edges,
)
from grnbench.report import (environment_block, file_fingerprint, thread_count, write_json,
                             sequence_identity, require_empty_output)


def default_method_names(methods_table: list[dict], *, everything: bool = False) -> list[str]:
    """Names ``--methods available`` (default set) or ``--methods all`` runs, sorted."""
    return sorted(m["name"] for m in methods_table
                  if m["available"] and (everything or m.get("default", True)))


def _write_state(out_dir: Path, manifest: dict) -> None:
    """Keep the detailed provenance and portable diagnostic in sync."""
    write_json(out_dir / "inference_manifest.json", manifest)
    from grnbench.diagnostics import model_document, try_write_diagnostics
    try_write_diagnostics(out_dir / "diagnostics.json", model_document, manifest)


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="grnbench-infer",
                                description="Run GRN inference methods (canonical edge lists).")
    p.add_argument("data", nargs="?", help="input dataset (.h5ad; counts or preprocessed)")
    p.add_argument("--out", help="output directory")
    p.add_argument("--methods", default="available",
                   help="comma-separated names, 'available' (default set), 'all' (every available "
                        "method), or 'list'")
    p.add_argument("--tf-list", default=None, help="regulator whitelist (.txt/.csv/.tsv)")
    p.add_argument("--assume-normalized", action="store_true",
                   help="skip internal cpm_log1p even if X looks like counts")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--threads", type=thread_count, default=4,
                   help="CPU budget (>= 1) for random_forest, extra_trees, GENIE3, GRNBoost2 and ARACNe3; default 4, not every core")
    p.add_argument("--plugin", action="append", default=None, metavar="SPEC",
                   help="custom model: spec file (.json/.yaml) describing a command/"
                        "Docker adapter, or 'package.module:Attribute' for a Python "
                        "one (repeatable)")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)

    plugins: list[dict] = []
    if args.plugin:
        from grnbench.plugins import register_plugins
        try:
            plugins = register_plugins(args.plugin)
        except MethodError as e:
            print(f"[PLUGIN ERROR] {e}", file=sys.stderr)
            return 1
        for p in plugins:
            print(f"[PLUGIN] registered '{p['name']}' ({p['origin']})")

    methods_table = list_methods()
    if args.methods == "list":
        for m in methods_table:
            status = "available" if m["available"] else f"UNAVAILABLE ({m['unavailable_reason']})"
            extra = "" if m.get("default", True) else " [extra]"
            print(f"{m['name']:16s} {m['origin']:12s} {status:50s} {m['description']}{extra}")
        return 0
    if not args.data or not args.out:
        print("error: both DATA and --out are required (or use --methods list)",
              file=sys.stderr)
        return 1

    known = {m["name"] for m in methods_table}
    available = {m["name"] for m in methods_table if m["available"]}
    if args.methods in ("available", "all"):
        chosen = default_method_names(methods_table, everything=args.methods == "all")
    else:
        chosen = [m.strip() for m in args.methods.split(",") if m.strip()]

    # A misspelled name is a user error, NOT an unavailable method: reporting
    # "geni3: unknown" next to "genie3: arboreto not installed" makes a typo
    # look like a missing dependency. Refuse the run and name the alternatives.
    unknown = [m for m in chosen if m not in known]
    if unknown:
        import difflib
        for m in unknown:
            close = difflib.get_close_matches(m, sorted(known), n=3, cutoff=0.6)
            hint = f" Did you mean: {', '.join(close)}?" if close else ""
            print(f"[UNKNOWN METHOD] '{m}' is not registered.{hint}", file=sys.stderr)
        print(f"registered methods: {', '.join(sorted(known))}\n"
              "(add a custom model with --plugin; see `grnbench-infer --methods list`)",
              file=sys.stderr)
        return 1

    skipped = [m for m in chosen if m not in available]
    chosen = [m for m in chosen if m in available]
    for m in skipped:
        reason = next((x["unavailable_reason"] for x in methods_table if x["name"] == m),
                      "unknown")
        print(f"[SKIP] {m}: {reason}", file=sys.stderr)
    tf_list = None
    manifest = {'args': vars(args), 'runs': {}, 'errors': {}, 'status': 'running',
                'environment': environment_block()}
    out_dir = None
    try:
        out_dir = require_empty_output(args.out)
        if args.tf_list:
            from grnbench.tfs import load_tf_list
            tf_list = load_tf_list(args.tf_list)
        import anndata as ad
        from grnbench.infer.base import prepare_expression
        from grnbench.tfs import resolve_tf_names
        import pandas as pd

        adata = ad.read_h5ad(args.data)
        _, norm_note = prepare_expression(adata, assume_normalized=args.assume_normalized)
        effective_tfs = resolve_tf_names(adata.var_names, tf_list)
        notes: list[str] = []
        if tf_list is None:      # the universe is genes x genes: say so before any model runs (D-003)
            n = len(effective_tfs)
            notes.append(f"no --tf-list: all {n:,} genes are treated as regulators, {n * (n - 1):,} candidate edges for "
                         "each method (time and memory grow with the square of the gene count). "
                         "Pass --tf-list to restrict the regulators.")
            print(f"[NO TF LIST] {notes[-1]}", file=sys.stderr)
        manifest.update({
            "input": file_fingerprint(Path(args.data)),
            "effective_genes": sequence_identity(adata.var_names),
            "effective_tfs": sequence_identity(effective_tfs),
            "tf_input": file_fingerprint(args.tf_list) if args.tf_list else None,
            "normalization": norm_note,
            "tf_list_size": len(tf_list) if tf_list else None,
            "seed": args.seed,
            "plugins": plugins,
            "methods_skipped": {
                m: next((x["unavailable_reason"] for x in methods_table
                         if x["name"] == m), "unknown")
                for m in skipped
            },
            "environment": environment_block(),
            "runs": {},
        })
        if notes:
            manifest["warnings"] = notes
        # created up front: the manifest must be writable even if every method
        # is skipped, otherwise write_json below dies with a FileNotFoundError
        universe = pd.DataFrame(((s, str(t)) for s in effective_tfs for t in adata.var_names
                                 if s != str(t)), columns=['source', 'target'])
        if universe.empty:
            raise ValueError('No candidate edges: effective genes and TFs must define a nonempty universe')
        universe_path = out_dir / 'candidate_universe.tsv.gz'
        universe.to_csv(universe_path, sep='\t', index=False, compression='gzip')
        manifest['candidate_universe'] = file_fingerprint(universe_path)
        manifest['candidate_universe']['n_pairs'] = len(universe)
        # a skipped method is not an error: its reason stays in methods_skipped and the others still run
        if not chosen:
            manifest['errors']['selection'] = 'No requested method is available'
            print('error: no requested method can run in this environment', file=sys.stderr)
        _write_state(out_dir, manifest)
        for name in chosen:
            print(f"running {name} ...")
            try:
                edges, prov = run_method(name, adata, tf_list=effective_tfs,
                                         assume_normalized=args.assume_normalized,
                                         seed=args.seed, threads=args.threads)
                path = save_edges(edges, out_dir, name)
                prov["edges_file"] = str(path)
                prov['edges_fingerprint'] = file_fingerprint(path)
                manifest["runs"][name] = prov
                print(f"  -> {len(edges)} edges -> {path}")
            except Exception as e:
                manifest['errors'][name] = f'{type(e).__name__}: {e}'
                print(f'[METHOD ERROR] {name}: {e}', file=sys.stderr)
            _write_state(out_dir, manifest)
        manifest['status'] = ('partial_failure' if manifest['errors'] and manifest['runs'] else
                              'failed' if manifest['errors'] else
                              'partial' if skipped else 'complete')
        _write_state(out_dir, manifest)
    except Exception as e:  # noqa: BLE001 - CLI boundary
        if out_dir is not None:
            manifest['status'] = 'failed'
            manifest['errors']['setup'] = f'{type(e).__name__}: {e}'
            _write_state(out_dir, manifest)
        print(f"[ERROR] {type(e).__name__}: {e}", file=sys.stderr)
        return 1

    print(f"\ninference complete: {len(manifest['runs'])} method(s)")
    return (2 if not chosen else 1) if manifest['errors'] else 0


if __name__ == "__main__":
    raise SystemExit(main())
