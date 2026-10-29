"""grnbench-eval CLI: score one or more ranked edge lists against a gold standard.

    grnbench-eval --gold GOLD.tsv --out DIR m1.edges.tsv.gz m2.edges.tsv.gz

Writes ``eval_report.json`` + ``eval_report.md`` + ``eval_report.html`` comparing
all methods, with conditional score-permutation p-values, cross-method top-k overlap
and — when ``--stability-data`` points at the processed matrix — subsample
stability re-run per method.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

from grnbench.eval.gold import coverage_report, load_gold
from grnbench.eval.metrics import (
    GoldIndex,
    evaluate_ranking,
    identifier_overlap,
    top_k_overlap,
)
from grnbench.eval.nulls import null_metrics
from grnbench.eval.universe import load_universe, align_rankings
from grnbench.report import thread_count


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="grnbench-eval",
                                description="Benchmark ranked GRN edge lists vs a gold standard.")
    p.add_argument("edges", nargs="*", help="edge list files (.edges.tsv.gz / .tsv / .csv)")
    p.add_argument("--gold", required=True,
                   help="path to gold pairs file, or a dorothea<levels> spec like dorotheaABC")
    p.add_argument("--out", default="eval_out", help="output directory")
    p.add_argument('--universe', help='fixed source,target CSV/TSV pair file; missing predictions form a bottom tie')
    p.add_argument("--top-k", type=int, default=None,
                   help="k for precision@k / EPR (default: reachable gold count, frozen before thinning)")
    p.add_argument("--n-nulls", type=int, default=50)
    p.add_argument("--null-mode", choices=["degree", "random", "regulator"], default="random")
    p.add_argument("--skip-nulls", action="store_true", help="faster: no null p-values")
    p.add_argument("--overlap-k", type=int, default=1000,
                   help="top-k used for the cross-method Jaccard matrix")
    p.add_argument("--seed", type=int, default=42)

    g = p.add_argument_group("descriptive gold-thinning sensitivity (needs >= 2 edge lists)")
    g.add_argument('--gold-thinning-sensitivity', action='store_true', help='descriptive paired annotation-removal scenarios; no significance inference')
    g.add_argument("--thinning-metric", dest='significance_metric', default="epr",
                   choices=["epr", "auprc", "auprc_ratio", "precision_at_k", "auroc"])
    g.add_argument("--thinning-draws", dest='significance_boot', type=int, default=200,
                   help="number of gold-thinning scenarios (default 200)")
    g.add_argument("--thinning-frac", dest='significance_frac', type=float, default=0.5,
                   help="fraction of reachable gold pairs kept per scenario (default 0.5)")
    p.add_argument("--case-sensitive", action="store_true",
                   help="require gene symbols to match case exactly. Off by "
                        "default: expression matrices and ground-truth networks "
                        "routinely disagree on case (BEELINE's mESC data is MGI "
                        "case, its mouse networks are upper case) and exact "
                        "matching turns that into a silent AUPRC of NaN.")

    s = p.add_argument_group("subsample stability (re-runs each method)")
    s.add_argument("--stability-data", default=None,
                   help="processed .h5ad to re-run methods on; enables stability")
    s.add_argument("--stability-n", type=int, default=5, help="number of subsamples")
    s.add_argument('--algorithm-seed', type=int, default=42, help='fixed algorithm seed for stability reruns')
    s.add_argument('--vary-algorithm-seed', action='store_true', help='also vary algorithm seed; measures combined perturbations')
    s.add_argument("--stability-frac", type=float, default=0.8,
                   help="fraction of cells kept per subsample")
    s.add_argument("--stability-top-k", type=int, default=100,
                   help="top-k edge set compared across subsamples")
    s.add_argument("--threads", type=thread_count, default=4,
                   help="CPU budget (>= 1) for the stability re-runs of random_forest, extra_trees, GENIE3, "
                        "GRNBoost2 and ARACNe3; default 4, not every core")
    s.add_argument("--stability-methods", default=None,
                   help="comma-separated subset (default: every scored method that "
                        "is registered and available)")
    s.add_argument("--stability-tf-list", default=None,
                   help="regulator whitelist for the re-runs (match the original run)")
    s.add_argument("--assume-normalized", action="store_true",
                   help="stability re-runs: skip internal cpm_log1p")

    r = p.add_argument_group("reporting")
    r.add_argument("--prep-manifest", default=None,
                   help="manifest.json from grnbench-prep, embedded in the HTML report")
    r.add_argument("--diagnostics", default=None,
                   help="schema-v1 diagnostics.json rendered consistently in Markdown and HTML")
    r.add_argument("--plugin", action="append", default=None, metavar="SPEC",
                   help="custom model: spec file (.json/.yaml) or 'module:Attribute' "
                        "(repeatable; needed to re-run plugin methods for stability)")
    r.add_argument("--no-html", action="store_true", help="skip eval_report.html")
    r.add_argument("--defer-report-output", action="store_true", help=argparse.SUPPRESS)
    return p


def _method_name(path: Path, taken: set[str]) -> tuple[str, str | None]:
    """Derive a method name from an edge-list filename, keeping it unique.

    Scoring the same method from two runs (``runs/a/genie3.edges.tsv.gz`` and
    ``runs/b/genie3.edges.tsv.gz``) is an obvious thing to want, and both derive
    "genie3". Colliding names are qualified with the parent directory rather than
    letting the second file silently overwrite the first; the caller reports the
    rename. Returns ``(name, collided_with)``.
    """
    base = path.name.split(".")[0] or path.stem or "method"
    if base not in taken:
        return base, None
    parent = path.parent.name
    candidate = f"{parent}/{base}" if parent else base
    n = 2
    while candidate in taken:
        candidate = f"{parent}/{base}#{n}" if parent else f"{base}#{n}"
        n += 1
    return candidate, base


def _render_markdown(report: dict, diagnostics: dict | None = None) -> str:
    lines = ["# GRN-Bench evaluation report", "",
             f"- Gold standard: **{report['gold']}** ({report['gold_pairs']} pairs)",
             f"- Null model: {report['null_mode']} ({report['n_nulls']} draws)", "",
             "| method | AUPRC | AUPRC-ratio | EPR | AUROC | null p (AUPRC) | stability |",
             "|---|---|---|---|---|---|---|"]
    for name, m in report["methods"].items():
        null_cell = "-" if m.get("null") is None else "%.3f" % m["null"]["p_value"]
        stab = m.get("stability")
        stab_cell = "-" if not stab else "%.3f" % stab["jaccard_mean"]
        # when k lands inside a block of equal scores the EPR is an expectation
        # over orderings the method never asserted; show the interval it spans
        epr_cell = f"{m['epr']:.2f}"
        rng = m.get("epr_range") or []
        if (m.get("tie_block_at_k") or 0) > 1 and len(rng) == 2                 and abs(rng[1] - rng[0]) > 1e-9:
            epr_cell += f" ({rng[0]:.2f}-{rng[1]:.2f})"
        lines.append(
            f"| {name} | {m['auprc']:.4f} | {m['auprc_ratio']:.2f} | "
            f"{epr_cell} | {m['auroc']:.4f} | {null_cell} | {stab_cell} |"
        )
    lines.append("")
    lines += ['Scores use explicit U when supplied; otherwise they are conditional on the shared output support. '
              'Missing predictions form a tied bottom block. Unlabelled pairs are not certified biological negatives. '
              'AP/prevalence is a conventional ratio, not normalization by exact finite-sample random AP. '
              'Null p-values refer to the stated score-exchangeability hypothesis and are unadjusted across methods.', '']
    lines += ['| method | support origin | missing predictions | gold in U / total |', '|---|---|---|---|']
    for name, m in report['methods'].items():
        lines.append(f"| {name} | {m.get('support_origin', 'unspecified')} | {m.get('missing_predictions', 0)} | {m.get('n_gold_in_universe', m['n_gold_in_list'])} / {m.get('n_gold_total', report['gold_pairs'])} |")
    lines.append('')
    # If k fell inside a block of identically-scored edges, the single EPR is an
    # expectation over orderings the method never asserted. Say so here too, or
    # the Markdown and the HTML disagree about how solid the number is.
    tied = {n: m for n, m in report["methods"].items()
            if (m.get("tie_block_at_k") or 0) > 1
            and len(m.get("epr_range") or []) == 2
            and abs(m["epr_range"][1] - m["epr_range"][0]) > 1e-9}
    if tied:
        lines += ["### Tied scores at the top-k cut", "",
                  "| method | edges tied at the cut | EPR if tie-broken worst / best |",
                  "|---|---|---|"]
        for name, m in tied.items():
            lo, hi = m["epr_range"]
            lines.append(f"| {name} | {m['tie_block_at_k']:,} | "
                         f"{lo:.2f} / {hi:.2f} |")
        lines += ["", "The EPR above is the expectation under a random tie-break. "
                      "A wide interval means the method did not actually rank "
                      "those edges, so the headline number is a summary, not a "
                      "measurement.", ""]
    cov = report.get("coverage", {})
    if cov:
        lines += ["## Gold reachability", "",
                  "| method | gold reachable in candidates | TF overlap |",
                  "|---|---|---|"]
        for name, c in cov.items():
            lines.append(f"| {name} | {c['gold_reachable']}/{c['gold_pairs']} | "
                         f"{c['tf_coverage_in_data']} |")
        lines.append("")
    ids = report.get("identifiers") or {}
    flagged = {n: d for n, d in ids.items() if d.get("verdict") != "ok"}
    if flagged:
        lines += ["## Gene-identifier matching", "",
                  "| method | verdict | pairs matched (exact / folded) | "
                  "shared genes (exact / folded) |", "|---|---|---|---|"]
        for name, d in flagged.items():
            lines.append(
                f"| {name} | **{d['verdict']}** | "
                f"{d['matched_pairs_case_sensitive']:,} / "
                f"{d['matched_pairs_case_insensitive']:,} | "
                f"{d['shared_genes_case_sensitive']:,} / "
                f"{d['shared_genes_case_insensitive']:,} |")
        lines += ["",
                  "`case-mismatch`: the edge list and the gold agree on genes but "
                  "not on case; scores above use case-folded matching.  "
                  "`no-overlap`: they share no gene identifiers at all, so every "
                  "metric against this gold is vacuous.", ""]
    sig = report.get("gold_thinning_sensitivity")
    if sig:
        from grnbench.eval.significance import render_significance_md
        lines.append(render_significance_md(sig))

    ov = report.get("overlap") or {}
    if len(ov.get("methods") or []) >= 2:
        names = ov["methods"]
        lines += [f"## Cross-method overlap (top-{ov['k']} Jaccard)", "",
                  "| | " + " | ".join(names) + " |",
                  "|---" * (len(names) + 1) + "|"]
        for a in names:
            row = [("-" if a == b else f"{ov['jaccard'][a][b]:.3f}") for b in names]
            lines.append(f"| **{a}** | " + " | ".join(row) + " |")
        lines.append("")
    from grnbench.diagnostics import render_diagnostics_markdown
    lines.append(render_diagnostics_markdown(diagnostics))
    return "\n".join(lines)


def write_evaluation_reports(report: dict, out_dir: Path, *, prep_manifest: dict | None = None,
                             diagnostics: dict | None = None, html: bool = True) -> list[Path]:
    """Write every evaluation report from one result and one validated diagnostic structure."""
    from grnbench.report import write_json
    write_json(out_dir / "eval_report.json", report)
    md_path = out_dir / "eval_report.md"
    md_path.write_text(_render_markdown(report, diagnostics), encoding="utf-8")
    outputs = [out_dir / "eval_report.json", md_path]
    if html:
        from grnbench.html_report import render_html
        html_path = out_dir / "eval_report.html"
        html_path.write_text(render_html(report, prep_manifest, diagnostics), encoding="utf-8")
        outputs.append(html_path)
    return outputs


def _run_stability(args, report: dict, names: list[str]) -> None:
    """Re-run each scored method on cell subsamples; annotate ``report`` in place."""
    import anndata as ad

    from grnbench.eval.stability import subsample_stability
    from grnbench.infer import list_methods

    available = {m["name"] for m in list_methods() if m["available"]}
    wanted = ([m.strip() for m in args.stability_methods.split(",") if m.strip()]
              if args.stability_methods else names)

    tf_list = None
    if args.stability_tf_list:
        from grnbench.tfs import load_tf_list
        tf_list = load_tf_list(args.stability_tf_list)

    adata = ad.read_h5ad(args.stability_data)
    report["stability_config"] = {
        "data": str(args.stability_data),
        "n_subsamples": int(args.stability_n),
        "frac": float(args.stability_frac),
        "top_k": int(args.stability_top_k),
        "tf_list_size": len(tf_list) if tf_list else None,
        "seed": int(args.seed),
        "threads": int(args.threads),
    }
    for name in wanted:
        if name not in report["methods"]:
            record = {'stage': 'stability', 'method': name, 'error': 'not among the scored edge lists'}
            report.setdefault('failures' if args.stability_methods is not None else 'stability_skips', []).append(record)
            print(f"[STABILITY SKIP] {name}: not among the scored edge lists",
                  file=sys.stderr)
            continue
        if name not in available:
            record = {'stage': 'stability', 'method': name, 'error': 'not a registered/available method'}
            report.setdefault('failures' if args.stability_methods is not None else 'stability_skips', []).append(record)
            print(f"[STABILITY SKIP] {name}: not a registered/available method here "
                  "(pass --plugin for custom models)", file=sys.stderr)
            continue
        print(f"stability: re-running {name} on {args.stability_n} subsamples ...")
        try:
            report["methods"][name]["stability"] = subsample_stability(
                adata, name, tf_list=tf_list,
                frac=float(args.stability_frac),
                n_subsamples=int(args.stability_n),
                top_k=int(args.stability_top_k),
                assume_normalized=args.assume_normalized,
                seed=int(args.seed),
                algorithm_seed=args.algorithm_seed,
                vary_algorithm_seed=args.vary_algorithm_seed,
                threads=int(args.threads),
            )
        except Exception as e:  # noqa: BLE001 - one method must not sink the report
            report.setdefault('failures', []).append({'stage': 'stability', 'method': name, 'error': str(e)})
            print(f"[STABILITY FAIL] {name}: {type(e).__name__}: {e}", file=sys.stderr)


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    from grnbench.diagnostics import read_diagnostics
    try:
        diagnostics = read_diagnostics(args.diagnostics)
    except ValueError as exc:
        print(f"[REPORT FAIL] {exc}", file=sys.stderr)
        return 2
    from grnbench.report import (environment_block, write_json, file_fingerprint,
                                 sequence_identity, require_empty_output)
    try:
        out_dir = require_empty_output(args.out)
    except (ValueError, OSError) as exc:
        print(f'[OUTPUT ERROR] {exc}', file=sys.stderr)
        return 1
    failures = []
    if args.plugin:
        from grnbench.infer import MethodError
        from grnbench.plugins import register_plugins
        try:
            plugins = register_plugins(args.plugin)
        except MethodError as e:
            print(f"[PLUGIN ERROR] {e}", file=sys.stderr)
            return 1
    else:
        plugins = []

    try:
        gold = load_gold(args.gold)
    except Exception as e:  # noqa: BLE001 - CLI boundary
        print(f"[GOLD ERROR] {e}", file=sys.stderr)
        return 1
    if not args.edges:
        print("error: provide at least one edge-list file", file=sys.stderr)
        return 1

    from grnbench.infer import load_edges

    # one normalisation of the gold for the whole run: nulls re-score it
    # n_nulls times per method and a ChIP-seq network has ~10^6 pairs
    try:
        gold_index = GoldIndex(gold.pairs, case_sensitive=args.case_sensitive)
        universe = load_universe(args.universe, case_sensitive=args.case_sensitive) if args.universe else None
        provenance = {
            'arguments': vars(args).copy(),
            'argv': list(argv if argv is not None else sys.argv[1:]),
            'gold_input': file_fingerprint(args.gold) if Path(args.gold).is_file()
                          else {'spec': args.gold, 'origin': gold.origin},
            'universe_input': file_fingerprint(args.universe) if args.universe else None,
            'edge_inputs': {},
            'additional_inputs': {},
            'identifier_matching': 'case-sensitive' if args.case_sensitive else 'uppercase-injective',
        }
        for key in ('stability_data', 'stability_tf_list', 'prep_manifest'):
            value = getattr(args, key)
            if value:
                provenance['additional_inputs'][key] = file_fingerprint(value)
        if args.diagnostics:
            from grnbench.diagnostics import diagnostics_file
            diagnostic_source = diagnostics_file(args.diagnostics)
            diagnostic_fingerprint = file_fingerprint(diagnostic_source)
            diagnostic_fingerprint['path'] = diagnostic_source.name
            provenance['additional_inputs']['diagnostics'] = diagnostic_fingerprint
            provenance['arguments']['diagnostics'] = diagnostic_source.name
            for index, token in enumerate(provenance['argv']):
                if token == '--diagnostics' and index + 1 < len(provenance['argv']):
                    provenance['argv'][index + 1] = diagnostic_source.name
                elif token.startswith('--diagnostics='):
                    provenance['argv'][index] = f'--diagnostics={diagnostic_source.name}'
        provenance['plugin_files'] = [file_fingerprint(spec) for spec in (args.plugin or []) if Path(spec).is_file()]
    except Exception as e:
        print(f'[INPUT ERROR] {e}', file=sys.stderr)
        return 1

    methods_metrics: dict = {}
    coverage: dict = {}
    identifiers: dict = {}
    edge_frames: dict = {}
    collisions: list[dict] = []
    for path_str in args.edges:
        path = Path(path_str)
        name, collided = _method_name(path, set(provenance['edge_inputs']))
        if collided is not None:
            collisions.append({"file": str(path), "derived": collided,
                               "renamed_to": name})
            print(f"[NAME COLLISION] {path.name} also derives '{collided}'; "
                  f"scoring it as '{name}'", file=sys.stderr)
        try:
            provenance['edge_inputs'][name] = file_fingerprint(path)
            edges = load_edges(path)
            methods_metrics[name] = evaluate_ranking(edges, gold_index,
                                                     top_k=args.top_k, universe=universe)
            methods_metrics[name]["file"] = str(path)
        except Exception as e:  # noqa: BLE001
            failures.append({'stage': 'score', 'method': name, 'error': str(e)})
            print(f"[SKIP] {path.name}: {e}", file=sys.stderr)
            continue
        edge_frames[name] = edges
        candidates = universe if universe is not None else edges
        coverage[name] = coverage_report({name: candidates}, gold,
                                         case_sensitive=args.case_sensitive)[name]

        # a gene-name mismatch produces NaNs, not errors; say so out loud
        diag = identifier_overlap(candidates, gold_index)
        diag['scope'] = 'candidate-universe' if universe is not None else 'output-list'
        identifiers[name] = diag
        if diag["verdict"] == "no-overlap":
            print(f"[GOLD MISMATCH] {name}: shares NO gene identifier with gold "
                  f"'{gold.name}' - every metric below is vacuous.\n"
                  f"    edge list uses: {diag['example_edge_genes']}\n"
                  f"    gold uses:      {diag['example_gold_genes']}\n"
                  "    different symbol namespaces, or the wrong organism/gold file?",
                  file=sys.stderr)
        elif diag["verdict"] == "case-mismatch":
            print(f"[GOLD CASE] {name}: gene symbols differ in case between the "
                  f"edge list and gold '{gold.name}' "
                  f"({diag['matched_pairs_case_sensitive']:,} pairs match exactly, "
                  f"{diag['matched_pairs_case_insensitive']:,} case-insensitively). "
                  "Scored case-insensitively; pass --case-sensitive to refuse.",
                  file=sys.stderr)
        if coverage[name]["gold_reachable"] == 0:
            # No positive among the candidate pairs: AUPRC, EPR and every null are undefined. A run that
            # looks complete but scores nothing must not exit 0 (a pipeline would read it as a result).
            fold = (lambda x: str(x)) if args.case_sensitive else (lambda x: str(x).upper())
            seen = {(fold(s_), fold(t_)) for s_, t_ in zip(candidates["source"], candidates["target"])}
            swapped = any((fold(t_), fold(s_)) in seen for s_, t_ in gold.pairs)
            reason = ("no reference edge is among this method's candidate pairs, so every metric is vacuous"
                      + ("; the reference's source and target columns look swapped" if swapped else ""))
            failures.append({"stage": "gold_overlap", "method": name, "error": reason})
            print(f"[GOLD EMPTY] {name}: {reason}", file=sys.stderr)
        if not args.skip_nulls:
            try:
                methods_metrics[name]["null"] = null_metrics(
                    edges, gold_index, n_nulls=args.n_nulls, mode=args.null_mode,
                    seed=args.seed, universe=universe, top_k=args.top_k)
            except Exception as e:  # noqa: BLE001
                methods_metrics[name]["null"] = None
                failures.append({'stage': 'null', 'method': name, 'error': str(e)})
                print(f"[WARN] nulls failed for {name}: {e}", file=sys.stderr)

    if not methods_metrics:
        print("error: no edge list could be scored", file=sys.stderr)
        return 1

    try:
        aligned_frames = align_rankings(edge_frames, universe, case_sensitive=args.case_sensitive)
    except ValueError as e:
        print(f'[SUPPORT ERROR] {e}', file=sys.stderr)
        return 1
    if universe is None and len(edge_frames) > 1:
        for metrics in methods_metrics.values():
            metrics['support_origin'] = 'shared-output-support'

    # Snapshot the actual canonical reference and candidate set, including
    # remotely resolved gold specs. Pair identities do not depend on row order.
    import pandas as pd
    def pair_identity(pairs):
        result = sequence_identity(json.dumps(pair, separators=(',', ':')) for pair in pairs)
        result.pop('values')  # exact values are retained once in the hashed TSV
        result['encoding'] = 'sorted pairs encoded as compact JSON arrays'
        return result
    gold_pairs = sorted(gold_index.lookup)
    candidate_pairs = sorted(zip(next(iter(aligned_frames.values())).source,
                                 next(iter(aligned_frames.values())).target))
    gold_snapshot = out_dir / 'effective_gold.tsv'
    universe_snapshot = out_dir / 'candidate_universe.tsv.gz'
    pd.DataFrame(gold_pairs, columns=['source', 'target']).to_csv(gold_snapshot, sep='\t', index=False)
    pd.DataFrame(candidate_pairs, columns=['source', 'target']).to_csv(universe_snapshot, sep='\t', index=False,
                                                                   compression='gzip')
    provenance['effective_gold'] = {'identity': pair_identity(gold_pairs),
                                    'snapshot': file_fingerprint(gold_snapshot)}
    provenance['candidate_universe'] = {'identity': pair_identity(candidate_pairs),
                                        'snapshot': file_fingerprint(universe_snapshot)}
    provenance['reachable_gold'] = pair_identity(sorted(set(gold_pairs) & set(candidate_pairs)))

    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "gold": gold.name,
        "gold_origin": gold.origin,
        "gold_pairs": len(gold),
        "universe": str(args.universe) if args.universe else None,
        "evaluation_scope": 'explicit-universe' if universe is not None else 'conditional-on-shared-output-support',
        "failures": failures,
        "null_mode": args.null_mode if not args.skip_nulls else "skipped",
        "n_nulls": 0 if args.skip_nulls else args.n_nulls,
        "case_sensitive_symbols": bool(args.case_sensitive),
        "methods": methods_metrics,
        "coverage": coverage,
        "identifiers": identifiers,
        "overlap": top_k_overlap(edge_frames, k=args.overlap_k, universe=universe,
                                 case_sensitive=args.case_sensitive),
        "name_collisions": collisions,
        "plugins": plugins,
        "environment": environment_block(),
        "provenance": provenance,
    }

    if args.gold_thinning_sensitivity:
        if len(edge_frames) < 2:
            failures.append({'stage': 'gold_thinning', 'error': 'needs at least two edge lists'})
            print("[THINNING SKIP] needs at least two scorable edge lists",
                  file=sys.stderr)
        else:
            from grnbench.eval.significance import gold_thinning_sensitivity
            try:
                report["gold_thinning_sensitivity"] = gold_thinning_sensitivity(
                    edge_frames, gold_index,
                    metric=args.significance_metric,
                    n_draws=int(args.significance_boot),
                    frac=float(args.significance_frac),
                    top_k=args.top_k, seed=int(args.seed), universe=universe)
                print('gold-thinning sensitivity: descriptive scenarios computed')
            except Exception as e:  # noqa: BLE001 - never sink the report
                failures.append({'stage': 'gold_thinning', 'error': str(e)})
                print(f"[THINNING FAIL] {type(e).__name__}: {e}",
                      file=sys.stderr)

    if args.stability_data:
        try:
            _run_stability(args, report, list(methods_metrics))
        except Exception as e:  # noqa: BLE001
            failures.append({'stage': 'stability', 'error': str(e)})
            print(f"[WARN] stability stage failed: {type(e).__name__}: {e}",
                  file=sys.stderr)

    prep_manifest = None
    if args.prep_manifest:
        try:
            prep_manifest = json.loads(Path(args.prep_manifest).read_text(encoding="utf-8"))
        except Exception as e:  # noqa: BLE001
            print(f"[WARN] could not read --prep-manifest: {e}", file=sys.stderr)

    report['status'] = 'partial-failure' if failures else 'success'
    outputs = write_evaluation_reports(report, out_dir, prep_manifest=prep_manifest,
                                       diagnostics=diagnostics, html=not args.no_html)
    md = (out_dir / "eval_report.md").read_text(encoding="utf-8")

    if not args.defer_report_output:
        print(md)
        print("\noutputs: " + " , ".join(str(o) for o in outputs))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
