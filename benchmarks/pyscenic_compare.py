#!/usr/bin/env python3
"""Compare grnbench cisTarget pruning with the official pySCENIC 0.12.1 CLI.

The comparison deliberately reuses one expression matrix and one *full* GRNBoost2
adjacency table for both implementations.  Large intermediate files stay in the
experiment directory outside git; ``report`` writes only portable JSON/Markdown.

Commands are resumable and terminal failures remain visible:

``prepare``
    Recreate the registered GRNBoost2 adjacency, verify its first 20,000 rows
    against the kept real-benchmark artifact, and export common inputs.
``run``
    Run the unmodified ``pyscenic ctx`` executable and pySCENIC's own module
    generator in its isolated Python environment.
``report``
    Require every registered unit to be terminal, then compare modules,
    motif-supported pairs and descriptive reference-recovery scores.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

TERMINAL = {"completed", "permanent_error"}
REGULON_SUFFIX = re.compile(r"\s*\([+-]\)\s*$")


def _json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
                   encoding="utf-8")
    tmp.replace(path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _scrub(text: str) -> str:
    home = str(Path.home())
    return str(text).replace(home, "~") if home else str(text)


def read_units(path: str | Path) -> list[dict[str, str]]:
    units = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {"task", "arm", "context", "species"}
    missing = required - set(units.columns)
    if missing:
        raise ValueError(f"units file missing columns: {sorted(missing)}")
    if units.empty or units.duplicated(["task", "arm", "context"]).any():
        raise ValueError("units must be non-empty and unique by task/arm/context")
    if not set(units["species"]) <= {"human", "mouse"}:
        raise ValueError("species must be human or mouse")
    return units.to_dict("records")


def unit_id(unit: dict[str, str]) -> str:
    return f"{unit['task']}__{unit['arm']}__{unit['context']}"


def pair_union(modules: list[dict]) -> set[tuple[str, str]]:
    """Union TF-target membership, excluding the TF self-pair pySCENIC adds."""
    pairs = set()
    for module in modules:
        tf = str(module["tf"])
        for gene in module.get("genes", []):
            gene = str(gene)
            if tf != gene:
                pairs.add((tf, gene))
    return pairs


def load_regulon_pairs(path: str | Path) -> set[tuple[str, str]]:
    """Read the documented pySCENIC JSON ``{regulon_name: [targets...]}``."""
    payload = _json(Path(path))
    if not isinstance(payload, dict):
        raise ValueError("pySCENIC regulon JSON must be an object")
    pairs = set()
    for name, targets in payload.items():
        tf = REGULON_SUFFIX.sub("", str(name))
        if not isinstance(targets, list):
            raise ValueError(f"regulon {name!r} targets must be a list")
        for target in targets:
            target = str(target)
            if target != tf:
                pairs.add((tf, target))
    return pairs


def set_metrics(candidate: set, reference: set) -> dict:
    intersection = len(candidate & reference)
    union = len(candidate | reference)
    return {
        "candidate": len(candidate),
        "reference": len(reference),
        "intersection": intersection,
        "union": union,
        "exact": candidate == reference,
        "jaccard": intersection / union if union else 1.0,
        "precision": intersection / len(candidate) if candidate else (1.0 if not reference else 0.0),
        "recall": intersection / len(reference) if reference else (1.0 if not candidate else 0.0),
    }


def rerank_by_support(edges: pd.DataFrame, support: set[tuple[str, str]]) -> pd.DataFrame:
    """Put supported edges first and preserve base order within both blocks."""
    out = edges[["source", "target", "weight"]].copy()
    out = out.sort_values("weight", ascending=False, kind="stable").reset_index(drop=True)
    base_score = 1.0 - np.arange(len(out), dtype=float) / max(len(out), 1)
    mask = np.fromiter(((str(s), str(t)) in support
                        for s, t in zip(out["source"], out["target"])),
                       dtype=bool, count=len(out))
    out["weight"] = base_score + mask.astype(float)
    return out.sort_values("weight", ascending=False, kind="stable").reset_index(drop=True)


def _load_edges(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    rename = {"TF": "source", "target": "target", "importance": "weight"}
    frame = frame.rename(columns=rename)
    required = ["source", "target", "weight"]
    if not set(required) <= set(frame.columns):
        raise ValueError(f"{path.name}: missing edge columns")
    return frame[required].sort_values("weight", ascending=False, kind="stable").reset_index(drop=True)


def _same_kept_edges(full: pd.DataFrame, kept_path: Path) -> None:
    kept = pd.read_csv(kept_path, sep="\t")[["source", "target", "weight"]]
    top = full.head(len(kept)).reset_index(drop=True)
    if list(zip(top.source.astype(str), top.target.astype(str))) != \
            list(zip(kept.source.astype(str), kept.target.astype(str))):
        raise RuntimeError("recreated GRNBoost2 adjacency does not match kept artifact pairs")
    if not np.allclose(top.weight.to_numpy(float), kept.weight.to_numpy(float),
                       rtol=1e-12, atol=1e-14):
        raise RuntimeError("recreated GRNBoost2 adjacency does not match kept artifact weights")


def _grnbench_support(edges: pd.DataFrame, expr: pd.DataFrame, species: str,
                      cistarget_dir: Path) -> tuple[list[dict], set[tuple[str, str]]]:
    from grnbench.infer.cistarget import (load_annotation, load_rankings, modules_from_edges,
                                          prune, resource_files)
    modules = modules_from_edges(edges, expr)
    dbs, annotation = resource_files(species, cistarget_dir)
    if not dbs or annotation is None:
        raise FileNotFoundError(f"no complete cisTarget resources for {species}")
    tf_motifs = load_annotation(str(annotation))
    genes = tuple(map(str, expr.columns))
    support: dict[tuple[str, str], float] = {}
    for database in dbs:
        motifs, matched, ranks, n_total = load_rankings(str(database), genes)
        for pair, score in prune(modules, motifs, matched, ranks, n_total, tf_motifs).items():
            support[pair] = max(support.get(pair, -np.inf), score)
    portable = [{"tf": str(tf), "genes": sorted(map(str, targets))}
                for tf, targets in modules]
    return portable, set(support)


def prepare(args) -> int:
    import anndata as ad
    from grnbench.autoprep import resolve_strategy
    from grnbench.infer import prepare_expression, run_method

    source = Path(args.source_run).expanduser().resolve()
    out = Path(args.out).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    records = []
    for unit in read_units(args.units):
        uid = unit_id(unit)
        target = out / "units" / uid
        marker = target / "PREPARED.json"
        if marker.exists():
            records.append(_json(marker))
            print(f"skip prepared {uid}")
            continue
        artifact = source / "artifacts" / unit["task"]
        context = artifact / unit["arm"] / unit["context"]
        required = [context / "processed.h5ad", context / "grnboost2.tsv.gz",
                    artifact / "reference_selection.tsv", artifact / "reference_evaluation.tsv",
                    artifact / "decisions.json"]
        missing = [str(p) for p in required if not p.exists()]
        if missing:
            raise FileNotFoundError(f"{uid}: missing source artifacts: {missing}")

        decisions = _json(artifact / "decisions.json")
        choices = decisions["chosen_plans"][unit["arm"]]["choices"]
        strategy = resolve_strategy(choices)[0]
        assume_normalized = strategy != "raw_counts"
        selection = pd.read_csv(artifact / "reference_selection.tsv", sep="\t")
        evaluation = pd.read_csv(artifact / "reference_evaluation.tsv", sep="\t")
        tfs = sorted(set(selection.source.astype(str)) | set(evaluation.source.astype(str)))
        adata = ad.read_h5ad(context / "processed.h5ad")
        edges, provenance = run_method("grnboost2", adata, tf_list=tfs,
                                       assume_normalized=assume_normalized)
        _same_kept_edges(edges, context / "grnboost2.tsv.gz")
        expr, normalization = prepare_expression(adata, assume_normalized=assume_normalized)

        target.mkdir(parents=True, exist_ok=True)
        expression_path = target / "expression.tsv"
        adjacency_path = target / "adjacencies.tsv"
        reference_path = target / "reference_evaluation.tsv"
        expr.to_csv(expression_path, sep="\t", index=True, index_label="CellID")
        edges.rename(columns={"source": "TF", "weight": "importance"}) \
            [["TF", "target", "importance"]].to_csv(adjacency_path, sep="\t", index=False)
        evaluation[["source", "target"]].to_csv(reference_path, sep="\t", index=False)
        modules, support = _grnbench_support(edges, expr, unit["species"],
                                             Path(args.cistarget_dir).expanduser())
        _write_json(target / "grnbench_modules.json", modules)
        pd.DataFrame(sorted(support), columns=["source", "target"]).to_csv(
            target / "grnbench_support.tsv", sep="\t", index=False)
        record = {
            **unit,
            "id": uid,
            "status": "prepared",
            "strategy": strategy,
            "normalization": normalization,
            "n_cells": int(expr.shape[0]),
            "n_genes": int(expr.shape[1]),
            "n_tfs": len(tfs),
            "n_adjacencies": len(edges),
            "n_grnbench_modules": len(modules),
            "n_grnbench_supported_pairs": len(support),
            "grnboost2_provenance": provenance,
            "files": {p.name: _sha256(p) for p in
                      (expression_path, adjacency_path, reference_path,
                       target / "grnbench_modules.json", target / "grnbench_support.tsv")},
        }
        _write_json(marker, record)
        records.append(record)
        print(f"prepared {uid}: {len(edges)} adjacencies, {len(support)} supported pairs")
    _write_json(out / "PREPARED.json", {"status": "prepared", "units": records})
    return 0


def _verify_prepared(unit_dir: Path, record: dict) -> None:
    for name, expected in record["files"].items():
        path = unit_dir / name
        if not path.is_file() or _sha256(path) != expected:
            raise RuntimeError(f"prepared input changed or missing: {name}")


def run(args) -> int:
    prepared = Path(args.prepared).expanduser().resolve()
    root = Path(args.cistarget_dir).expanduser().resolve()
    pyscenic = Path(args.pyscenic).expanduser().resolve()
    if not pyscenic.is_file():
        raise FileNotFoundError(f"pySCENIC executable not found: {_scrub(pyscenic)}")
    py = pyscenic.with_name("python")
    if not py.is_file():
        raise FileNotFoundError(f"isolated Python not found beside pySCENIC: {_scrub(py)}")
    exporter = Path(__file__).with_name("pyscenic_reference_export.py").resolve()
    prepared_index = _json(prepared / "PREPARED.json")
    for record in prepared_index["units"]:
        uid = record["id"]
        unit_dir = prepared / "units" / uid
        run_path = unit_dir / "run.json"
        if run_path.exists() and _json(run_path).get("status") in TERMINAL:
            print(f"skip terminal {uid}")
            continue
        _verify_prepared(unit_dir, record)
        from grnbench.infer.cistarget import resource_files
        dbs, annotation = resource_files(record["species"], root)
        if not dbs or annotation is None:
            raise FileNotFoundError(f"no complete cisTarget resources for {record['species']}")
        log_path = unit_dir / "pyscenic.log"
        regulons_path = unit_dir / "pyscenic_regulons.json"
        modules_path = unit_dir / "pyscenic_modules.json"
        command = [str(pyscenic), "ctx", str(unit_dir / "adjacencies.tsv"),
                   *map(str, dbs), "--annotations_fname", str(annotation),
                   "--expression_mtx_fname", str(unit_dir / "expression.tsv"),
                   "--mode", "custom_multiprocessing", "--num_workers", str(args.workers),
                   "-o", str(regulons_path)]
        export_command = [str(py), str(exporter), "--adjacencies",
                          str(unit_dir / "adjacencies.tsv"), "--expression",
                          str(unit_dir / "expression.tsv"), "--out", str(modules_path)]
        started = time.time()
        payload = {"id": uid, "status": "running", "attempt": 1,
                   "command": [Path(x).name if i == 0 else _scrub(x)
                               for i, x in enumerate(command)]}
        _write_json(run_path, payload)
        try:
            with log_path.open("w", encoding="utf-8") as log:
                first = subprocess.run(export_command, stdout=log, stderr=subprocess.STDOUT,
                                       text=True, check=False)
                second = (subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                                         text=True, check=False)
                          if first.returncode == 0 else None)
            code = second.returncode if second is not None else first.returncode
            if code != 0 or not modules_path.is_file() or not regulons_path.is_file():
                tail = log_path.read_text(encoding="utf-8", errors="replace")[-2000:]
                payload.update(status="permanent_error", exit_code=int(code),
                               error=_scrub(tail))
            else:
                modules_payload = _json(modules_path)
                regulons_payload = _json(regulons_path)
                if not isinstance(modules_payload, list) or not modules_payload:
                    payload.update(status="permanent_error", exit_code=0,
                                   error="pySCENIC returned no activating modules")
                elif not isinstance(regulons_payload, dict) or not regulons_payload:
                    payload.update(status="permanent_error", exit_code=0,
                                   error="pySCENIC returned no regulons")
                else:
                    payload.update(status="completed", exit_code=0,
                                   files={modules_path.name: _sha256(modules_path),
                                          regulons_path.name: _sha256(regulons_path)})
        except Exception as exc:  # a terminal result, never silently omitted
            payload.update(status="permanent_error", exit_code=None,
                           error=_scrub(f"{type(exc).__name__}: {exc}"))
        payload["seconds"] = round(time.time() - started, 3)
        _write_json(run_path, payload)
        print(f"{uid}: {payload['status']} ({payload['seconds']}s)")
    terminal = []
    for record in prepared_index["units"]:
        path = prepared / "units" / record["id"] / "run.json"
        terminal.append(path.is_file() and _json(path).get("status") in TERMINAL)
    _write_json(prepared / "RUN_COMPLETE.json",
                {"complete": all(terminal), "expected": len(terminal),
                 "terminal": sum(terminal)})
    return 0 if all(terminal) else 2


def _pairs_tsv(path: Path) -> set[tuple[str, str]]:
    frame = pd.read_csv(path, sep="\t", dtype=str)
    return set(zip(frame.source, frame.target))


def _finite(value):
    value = float(value)
    return value if np.isfinite(value) else None


def _write_portable_outputs(output: Path, summary: dict) -> None:
    """Generate figures and LaTeX from COMPARISON.json values only."""
    completed = [row for row in summary["units"] if row["status"] == "completed"]
    viz = output / "viz"
    latex = output / "latex"
    viz.mkdir(parents=True, exist_ok=True)
    latex.mkdir(parents=True, exist_ok=True)
    index = ["# EXP-20261018-LD-1 figures", ""]
    if completed:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        labels = [row["species"] for row in completed]
        x = np.arange(len(labels))
        fig, ax = plt.subplots(figsize=(6.4, 3.8))
        ax.bar(x - 0.18, [row["modules"]["jaccard"] for row in completed], 0.36,
               label="module pairs")
        ax.bar(x + 0.18, [row["support"]["jaccard"] for row in completed], 0.36,
               label="motif-supported pairs")
        ax.set_xticks(x, labels)
        ax.set_ylim(0, 1)
        ax.set_ylabel("Jaccard overlap")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(viz / "implementation_overlap.png", dpi=160)
        plt.close(fig)

        methods = ("grnboost2", "grnbench_cistarget", "pyscenic_ctx")
        fig, ax = plt.subplots(figsize=(7.2, 4.0))
        width = 0.24
        for i, method in enumerate(methods):
            values = [row["scores"][method]["auprc_ratio"] for row in completed]
            ax.bar(x + (i - 1) * width, values, width, label=method)
        ax.axhline(1, color="black", linewidth=0.8, linestyle="--")
        ax.set_xticks(x, labels)
        ax.set_ylabel("AUPRC ratio (1 = random)")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(viz / "reference_recovery.png", dpi=160)
        plt.close(fig)
        index += ["- `implementation_overlap.png`: overlap of module and final supported-pair sets.",
                  "- `reference_recovery.png`: descriptive recovery on the registered CollecTRI "
                  "evaluation half.", ""]
    else:
        index += ["No figures produced: every registered unit ended in a permanent pySCENIC error.", ""]
    index += ["## Not produced", "",
              "- QC and embeddings: the experiment compares implementations on identical already-kept "
              "inputs; it does not compare preprocessing arms.",
              "- Attributions and stability: neither implementation is a learned end-to-end model here, "
              "and one deterministic adjacency is deliberately held fixed.", ""]
    (viz / "index.md").write_text("\n".join(index), encoding="utf-8")

    def tex_number(value) -> str:
        return "NA" if value is None else f"{float(value):.4f}"

    rows = []
    for row in completed:
        scores = row["scores"]
        rows.append(
            f"{row['species']} & {row['modules']['jaccard']:.4f} & "
            f"{row['support']['jaccard']:.4f} & "
            f"{tex_number(scores['grnboost2']['auprc_ratio'])} & "
            f"{tex_number(scores['grnbench_cistarget']['auprc_ratio'])} & "
            f"{tex_number(scores['pyscenic_ctx']['auprc_ratio'])} \\\\"
        )
    body = "\n".join(rows) if rows else "No completed units."
    tex = r"""\documentclass{article}
\usepackage[margin=1in]{geometry}
\usepackage{booktabs}
\title{EXP-20261018-LD-1: pySCENIC reference comparison}
\begin{document}
\maketitle
This is an exploratory, descriptive implementation comparison on identical expression
matrices and GRNBoost2 adjacencies. No inferential test is performed.
\begin{center}
\begin{tabular}{lrrrrr}
\toprule
Species & module Jaccard & support Jaccard & base AUPRC-ratio & grnbench & pySCENIC \\
\midrule
""" + body + r"""
\bottomrule
\end{tabular}
\end{center}
\end{document}
"""
    (latex / "report.tex").write_text(tex, encoding="utf-8")


def report(args) -> int:
    from grnbench.eval.metrics import evaluate_ranking

    prepared = Path(args.prepared).expanduser().resolve()
    output = Path(args.out).expanduser().resolve()
    expected = read_units(args.units)
    index = _json(prepared / "PREPARED.json")
    by_id = {record["id"]: record for record in index["units"]}
    expected_ids = [unit_id(unit) for unit in expected]
    if set(by_id) != set(expected_ids):
        raise RuntimeError("prepared units do not exactly match the registered unit manifest")

    results = []
    lines = ["# pySCENIC 0.12.1 reference comparison", "",
             "Exploratory, descriptive method-validation on the exact same expression matrices "
             "and GRNBoost2 adjacencies. No p-values are computed (n = 2 datasets).", ""]
    for uid in expected_ids:
        unit_dir = prepared / "units" / uid
        run_path = unit_dir / "run.json"
        if not run_path.exists() or _json(run_path).get("status") not in TERMINAL:
            raise RuntimeError(f"registered unit is not terminal: {uid}")
        run_record = _json(run_path)
        base = {"id": uid, "task": by_id[uid]["task"], "species": by_id[uid]["species"],
                "status": run_record["status"], "seconds": run_record.get("seconds")}
        if run_record["status"] == "permanent_error":
            base["error"] = run_record.get("error", "unknown error")
            results.append(base)
            lines += [f"## {uid}", "", "**Permanent error.**", "", "```text",
                      str(base["error"])[-1000:], "```", ""]
            continue

        grn_modules = pair_union(_json(unit_dir / "grnbench_modules.json"))
        py_modules = pair_union(_json(unit_dir / "pyscenic_modules.json"))
        grn_support = _pairs_tsv(unit_dir / "grnbench_support.tsv")
        py_support = load_regulon_pairs(unit_dir / "pyscenic_regulons.json")
        edges = _load_edges(unit_dir / "adjacencies.tsv")
        universe = set(zip(edges.source.astype(str), edges.target.astype(str)))
        gold = _pairs_tsv(unit_dir / "reference_evaluation.tsv")
        rankings = {
            "grnboost2": edges,
            "grnbench_cistarget": rerank_by_support(edges, grn_support),
            "pyscenic_ctx": rerank_by_support(edges, py_support),
        }
        scores = {}
        explicit_universe = edges[["source", "target"]]
        for name, ranking in rankings.items():
            metric = evaluate_ranking(ranking, gold, with_curve=False,
                                      universe=explicit_universe)
            scores[name] = {k: _finite(metric[k]) for k in ("auprc_ratio", "epr", "auprc")}
        base.update(
            modules=set_metrics(grn_modules, py_modules),
            support=set_metrics(grn_support, py_support),
            support_outside_universe={"grnbench": len(grn_support - universe),
                                      "pyscenic": len(py_support - universe)},
            scores=scores,
        )
        results.append(base)
        lines += [f"## {uid}", "",
                  f"- Module pairs: grnbench {len(grn_modules)}, pySCENIC {len(py_modules)}, "
                  f"Jaccard {base['modules']['jaccard']:.6f}, exact {base['modules']['exact']}.",
                  f"- Motif-supported pairs: grnbench {len(grn_support)}, pySCENIC {len(py_support)}, "
                  f"Jaccard {base['support']['jaccard']:.6f}, exact {base['support']['exact']}.",
                  f"- AUPRC ratio: GRNBoost2 {scores['grnboost2']['auprc_ratio']}; "
                  f"grnbench {scores['grnbench_cistarget']['auprc_ratio']}; "
                  f"pySCENIC {scores['pyscenic_ctx']['auprc_ratio']}.", ""]

    terminal = len(results)
    completed = sum(r["status"] == "completed" for r in results)
    errors = terminal - completed
    summary = {"experiment": "EXP-20261018-LD-1", "expected": len(expected_ids),
               "terminal": terminal, "completed": completed, "permanent_errors": errors,
               "all_terminal": terminal == len(expected_ids), "units": results,
               "interpretation": "exploratory descriptive method validation; no confirmatory claim"}
    output.mkdir(parents=True, exist_ok=True)
    _write_json(output / "COMPARISON.json", summary)
    (output / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")
    _write_json(output / "COMPLETE.json",
                {"complete": terminal == len(expected_ids), "expected": len(expected_ids),
                 "terminal": terminal, "completed": completed, "permanent_errors": errors})
    _write_portable_outputs(output, summary)
    print(f"terminal={terminal}/{len(expected_ids)} completed={completed} permanent_errors={errors}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--source-run", required=True)
    p.add_argument("--units", required=True)
    p.add_argument("--cistarget-dir", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=prepare)
    p = sub.add_parser("run")
    p.add_argument("--prepared", required=True)
    p.add_argument("--cistarget-dir", required=True)
    p.add_argument("--pyscenic", required=True)
    p.add_argument("--workers", type=int, default=8)
    p.set_defaults(func=run)
    p = sub.add_parser("report")
    p.add_argument("--prepared", required=True)
    p.add_argument("--units", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=report)
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
