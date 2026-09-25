"""Engine mass test on real datasets (no planted network).

What can be measured without a true network:

* **detection** - protocol and organism found by the engine vs the curated
  metadata (``truth.json`` written by fetch_cellxgene.py, or given per dataset);
* **robustness** - every dataset either runs end to end or is refused with a
  clear reason (e.g. BEELINE inputs are log-normalized, not counts);
* **reference recovery** - for human/mouse data, the CollecTRI network split in
  selection/evaluation halves exactly as in the synthetic bench: pipelines are
  chosen on one half and scored on the other.

    python benchmarks/real_bench.py plan --data-root D --out R/real_tasks.json
    python benchmarks/real_bench.py run  --tasks R/real_tasks.json --refs D/references --out R --jobs 4
    python benchmarks/engine_bench.py report --out R      # same report format

Model comparison (``--models``): after screening, every listed inference method is
scored on the *baseline* and *auto* pipelines (evaluation half only), AUPRC ratio and
EPR, dataset = unit. ``X_motif`` / ``X_cistarget`` reuse the edges of ``X`` computed on
the same context, then add the cisTarget motif prior (GRNBENCH_CISTARGET_DIR).

    python benchmarks/real_bench.py run --tasks T --refs D/references --out R         --models pearson_signed,spearman_signed,clr,grnboost2,motif_prior,pearson_motif,grnboost2_motif
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import tempfile
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from engine_bench import BASELINE, EVAL_METHODS, keep_awake, score_candidates  # noqa: E402

MIN_REFERENCE_EDGES = 50
# Regulators scored: the TFs with the most reference targets present in the data. Evaluation
# only has power for TFs with many known targets, and mutual-information methods scale with
# TFs x genes (800 TFs made one dataset take ~40 min). A bench setting, reported as such.
MAX_TFS = 150


def make_tasks(data_root: Path) -> list[dict]:
    root = Path(data_root)
    tasks = []
    for truth in sorted(root.glob("cellxgene/**/truth.json")):
        rel = truth.parent.relative_to(root / "cellxgene").as_posix().replace("/", "__")
        tasks.append({"id": f"cxg__{rel}", "path": str(truth.parent / "data.h5ad"),
                      "truth": json.loads(truth.read_text(encoding="utf-8")), "kind": "cellxgene"})
    pbmc = next(iter(sorted(root.glob("10x_pbmc3k/filtered/**/matrix.mtx*"))), None)
    if pbmc is not None:
        tasks.append({"id": "pbmc3k_filtered", "path": str(pbmc.parent), "kind": "10x",
                      "truth": {"technology": "umi", "organism": "human"}})
    raw = next(iter(sorted(root.glob("10x_pbmc3k/raw/**/matrix.mtx*"))), None)
    if raw is not None:
        tasks.append({"id": "pbmc3k_raw", "path": str(raw.parent), "kind": "10x_raw",
                      "truth": {"technology": "umi", "organism": "human", "raw_droplets": True}})
    paul = root / "scanpy_paul15" / "paul15.h5ad"
    if paul.exists():
        tasks.append({"id": "paul15", "path": str(paul), "kind": "h5ad",
                      "truth": {"technology": "umi", "organism": "mouse", "time_course": True}})
    for expr in sorted(root.glob("beeline/**/ExpressionData.csv"))[:7]:
        tasks.append({"id": f"beeline__{expr.parent.name}", "path": str(expr), "kind": "beeline",
                      "truth": {"expect_refusal": "not raw counts (BEELINE inputs are normalized)"}})
    return tasks


def _reference(refs: Path | None, organism: str | None, genes) -> set[tuple[str, str]]:
    if refs is None or organism not in ("human", "mouse"):
        return set()
    path = Path(refs) / f"collectri_{organism}.tsv"
    if not path.exists():
        return set()
    ref = pd.read_csv(path, sep="\t")
    upper = {str(g).upper(): str(g) for g in genes}
    pairs = set()
    for s, t in zip(ref["source"].astype(str), ref["target"].astype(str)):
        if s.upper() in upper and t.upper() in upper and s.upper() != t.upper():
            pairs.add((upper[s.upper()], upper[t.upper()]))
    degree: dict[str, int] = {}
    for s, _ in pairs:
        degree[s] = degree.get(s, 0) + 1
    keep = set(sorted(degree, key=lambda k: (-degree[k], k))[:MAX_TFS])
    return {(s, t) for s, t in pairs if s in keep}


COMBINERS = ("_motif", "_cistarget")

# Pre-declared comparison family for the report (declared 2026-09-25 before the run):
# does the motif prior add to each base model, does expression add to the prior alone,
# and do the newer models beat the Pearson default?
DECLARED_PAIRS = [("pearson_motif", "pearson_signed"), ("grnboost2_motif", "grnboost2"),
                  ("regdiffusion_motif", "regdiffusion"), ("grnboost2_cistarget", "grnboost2"),
                  ("pearson_motif", "motif_prior"), ("grnboost2", "pearson_signed"),
                  ("regdiffusion", "pearson_signed")]


BASE_ALIAS = {"pearson": "pearson_signed", "spearman": "spearman_signed"}


def _split_combined(method: str) -> tuple[str, str | None]:
    """'pearson_motif' -> ('pearson_signed', '_motif'); plain methods -> (method, None)."""
    for suffix in COMBINERS:
        if method.endswith(suffix) and method != "motif_prior":
            base = method[: -len(suffix)]
            return BASE_ALIAS.get(base, base), suffix
    return method, None


ARTIFACT_TOP_EDGES = 20000     # per model/arm/context edge dump size when --keep-artifacts


def _write_edges_artifact(edges: pd.DataFrame, path: Path) -> None:
    top = edges.sort_values("weight", ascending=False, kind="stable").head(ARTIFACT_TOP_EDGES)
    path.parent.mkdir(parents=True, exist_ok=True)
    top[["source", "target", "weight"]].to_csv(path, sep="\t", index=False, compression="gzip")


def score_models(records_by_cand, plans, arms: dict, methods: list[str], tfs, eval_pairs,
                 artifacts_dir: Path | None = None) -> dict:
    """AUPRC ratio + EPR of every method on the arms' candidates, mean over contexts.

    Returns ``{f"m_{method}_{arm}_{metric}": value}`` plus per-method seconds.

    ``artifacts_dir``, when given, writes per arm/context the top-scored edges
    of every model (``<arm>/ctx<i>/<model>.tsv.gz``) and the processed AnnData
    actually scored (``<arm>/ctx<i>/processed.h5ad`` + ``context.json``). This
    disables the baseline==auto dedup so every arm gets its own files; it does
    not change the returned scores.
    """
    import anndata as ad
    from engine_bench import MAX_SCORE_CELLS
    from grnbench.eval.metrics import GoldIndex, evaluate_ranking
    from grnbench.infer import prepare_expression, run_method
    from grnbench.infer.cistarget import cistarget_rerank, motif_rank_product
    from grnbench.tfs import resolve_tf_names

    gi = GoldIndex.coerce(set(eval_pairs))
    out, done = {}, {}
    for arm, cand in arms.items():
        if cand is None or cand not in records_by_cand:
            continue
        if cand in done and artifacts_dir is None:   # auto == baseline: same numbers, no rerun
            for k, v in done[cand].items():
                out[k.replace("__ARM__", arm)] = v
            continue
        assume = plans[cand].value("strategy") != "raw_counts"
        vals = {m: {"auprc_ratio": [], "epr": [], "seconds": 0.0} for m in methods}
        for i, rec in enumerate(records_by_cand[cand]):
            if not rec.get("processed"):
                continue
            adata = ad.read_h5ad(rec["processed"])
            if adata.n_obs > MAX_SCORE_CELLS:
                keep = np.sort(np.random.default_rng(0).choice(adata.n_obs, MAX_SCORE_CELLS, replace=False))
                adata = adata[keep].copy()
            genes = pd.Index(adata.var_names).astype(str)
            sources = resolve_tf_names(genes, tfs)
            universe = pd.DataFrame([(s_, g) for s_ in sources for g in genes if s_ != g],
                                    columns=["source", "target"])
            expr, _ = prepare_expression(adata, assume_normalized=assume)
            ctx_dir = None
            if artifacts_dir is not None:
                ctx_dir = Path(artifacts_dir) / arm / f"ctx{i}"
                ctx_dir.mkdir(parents=True, exist_ok=True)
                adata.write_h5ad(ctx_dir / "processed.h5ad", compression="gzip")
                (ctx_dir / "context.json").write_text(json.dumps(
                    {"context": rec.get("context"), "n_obs": int(adata.n_obs),
                     "n_vars": int(adata.n_vars)}, indent=1), encoding="utf-8")
            cache: dict[str, pd.DataFrame] = {}
            for m in methods:
                t0 = time.time()
                try:
                    base, combiner = _split_combined(m)
                    if base not in cache:
                        cache[base], _ = run_method(base, adata, tf_list=tfs, assume_normalized=assume)
                    edges = cache[base]
                    if combiner == "_motif":
                        edges = motif_rank_product(edges, expr)
                    elif combiner == "_cistarget":
                        edges = cistarget_rerank(edges, expr)
                    if ctx_dir is not None:
                        _write_edges_artifact(edges, ctx_dir / f"{m}.tsv.gz")
                    r = evaluate_ranking(edges, gi, with_curve=False, universe=universe)
                    vals[m]["auprc_ratio"].append(float(r.get("auprc_ratio", np.nan)))
                    vals[m]["epr"].append(float(r.get("epr", np.nan)))
                except Exception as exc:  # noqa: BLE001 - a failing method is recorded, not fatal
                    vals[m]["auprc_ratio"].append(np.nan)
                    vals[m]["epr"].append(np.nan)
                    out[f"m_{m}_{arm}_error"] = f"{type(exc).__name__}: {exc}"[:300]
                vals[m]["seconds"] += time.time() - t0
        block = {}
        for m, v in vals.items():
            for metric in ("auprc_ratio", "epr"):
                x = [y for y in v[metric] if np.isfinite(y)]
                block[f"m_{m}___ARM___{metric}"] = float(np.mean(x)) if x else np.nan
            block[f"m_{m}___ARM___seconds"] = round(v["seconds"], 1)
        done[cand] = block
        for k, v in block.items():
            out[k.replace("__ARM__", arm)] = v
    return out


def run_real_task(task: dict, refs: Path | None, workroot: Path, keep: bool = False,
                  models: list[str] | None = None, artifacts_dir: Path | None = None) -> dict:
    from grnbench.autoprep import diagnose, execute_plan, plan_preprocessing
    from grnbench.screen import rank_candidates, screen_candidates

    t0 = time.time()
    row = {"task": task["id"], "scenario": f"real_{task['kind']}", "seed": 0, "error": None,
           "truth": json.dumps(task["truth"])}
    work = Path(workroot) / task["id"]
    shutil.rmtree(work, ignore_errors=True)
    try:
        orientation = "genes_x_cells" if task["kind"] == "beeline" else "auto"
        try:
            diag = diagnose(task["path"], orientation=orientation)
        except Exception as exc:  # noqa: BLE001
            if task["truth"].get("expect_refusal"):
                row["refused_as_expected"] = True
                row["refusal"] = f"{type(exc).__name__}: {exc}"[:500]
                row["seconds"] = round(time.time() - t0, 1)
                return row
            raise
        if task["truth"].get("expect_refusal"):
            row["refused_as_expected"] = False
        truth = task["truth"]
        row["n_cells_in"], row["n_genes_in"] = int(diag.adata.n_obs), int(diag.adata.n_vars)
        pairs = _reference(refs, truth.get("organism"), diag.adata.var_names)
        tfs = sorted({s for s, _ in pairs}) or None
        plan = plan_preprocessing(diag, tf_list=tfs)
        prof = plan.profile
        for key in ("technology", "organism"):
            if key in truth:
                row[f"detect_{key}"] = str(prof.get(key)).lower() == str(truth[key]).lower()
                row[f"got_{key}"] = prof.get(key)
        if truth.get("raw_droplets"):
            row["detect_cell_calling_best"] = plan.selection["results"]["cell_calling"]["best_overall"] in {
                "emptydrops", "knee"}
        if truth.get("time_course"):
            row["detect_time_course"] = bool(prof.get("time_course"))
        row["rules_choices"] = json.dumps(plan.choices)
        row["n_warnings"] = len(plan.warnings)
        row["reference_edges"] = len(pairs)

        if len(pairs) >= MIN_REFERENCE_EDGES:
            res = screen_candidates(diag, task["path"], work / "screen", tf_list=tfs,
                                    reference_pairs=pairs, max_candidates=4,
                                    extra=[{"name": "baseline", "choices": BASELINE}],
                                    proxy_methods=EVAL_METHODS)
            scores = score_candidates(res["records"], res["plans"], tfs, res["reference"]["evaluation"])
            noref = res["rows"].drop(columns=[c for c in res["rows"] if c.startswith("reference_")])
            noref_rank, _ = rank_candidates(noref)
            arms = {"baseline": "baseline", "rules": "recommended", "auto": res["best"],
                    "auto_noref": str(noref_rank.index[0]) if len(noref_rank) else None}
            finite = {k: v["eval_mean"] for k, v in scores.items() if np.isfinite(v["eval_mean"])}
            arms["oracle"] = max(finite, key=finite.get) if finite else None
            for arm, cand in arms.items():
                row[f"{arm}_candidate"] = cand
                for k, v in scores.get(cand, {}).items():
                    row[f"{arm}_{k}"] = v
                row[f"{arm}_failed_contexts"] = sum(1 for r in res["records"].get(cand, []) if r.get("error"))
            for m in EVAL_METHODS:
                best_m = res["best_per_method"].get(m)
                row[f"permodel_{m}_candidate"] = best_m
                row[f"permodel_{m}_auto"] = scores.get(best_m, {}).get(f"eval_{m}", np.nan)
                row[f"permodel_{m}_baseline"] = scores.get("baseline", {}).get(f"eval_{m}", np.nan)
                row[f"permodel_{m}_rules"] = scores.get("recommended", {}).get(f"eval_{m}", np.nan)
                vals = [v.get(f"eval_{m}", np.nan) for v in scores.values()]
                vals = [v for v in vals if np.isfinite(v)]
                row[f"permodel_{m}_oracle"] = max(vals) if vals else np.nan
            task_artifacts = Path(artifacts_dir) / task["id"] if artifacts_dir is not None else None
            if task_artifacts is not None:
                task_artifacts.mkdir(parents=True, exist_ok=True)
                for half, pairs_ in (("selection", res["reference"]["selection"]),
                                     ("evaluation", res["reference"]["evaluation"])):
                    pd.DataFrame(sorted(pairs_), columns=["source", "target"]).to_csv(
                        task_artifacts / f"reference_{half}.tsv", sep="\t", index=False)
                decisions = {"rules_choices": plan.choices, "candidates": dict(arms), "chosen_plans": {}}
                for arm, cand in arms.items():
                    cand_plan = res["plans"].get(cand)
                    if cand_plan is None:
                        continue
                    decisions["chosen_plans"][arm] = {
                        "candidate": cand, "selection": cand_plan.selection,
                        "choices": cand_plan.choices, "warnings": cand_plan.warnings}
                (task_artifacts / "decisions.json").write_text(
                    json.dumps(decisions, indent=1, default=str), encoding="utf-8")
            if models:
                row.update(score_models(res["records"], res["plans"],
                                        {"baseline": "baseline", "auto": res["best"]},
                                        models, tfs, res["reference"]["evaluation"],
                                        artifacts_dir=task_artifacts))
        else:
            recs = execute_plan(plan, task["path"], work / "clean", tf_list=tfs)
            row["contexts_ok"] = sum(1 for r in recs if r["processed"])
            row["contexts_failed"] = sum(1 for r in recs if r["error"])
            if row["contexts_failed"]:
                row["context_errors"] = " | ".join(r["error"] for r in recs if r["error"])[:1000]
    except Exception as exc:  # noqa: BLE001 - robustness is part of what is measured
        row["error"] = f"{type(exc).__name__}: {exc}"
        row["traceback"] = traceback.format_exc()[-4000:]
    row["seconds"] = round(time.time() - t0, 1)
    if not keep:
        shutil.rmtree(work, ignore_errors=True)
    return row


def _run_and_save(args):
    task, refs, out, workroot, models, keep_artifacts = args
    path = Path(out) / "rows" / f"{task['id']}.json"
    if path.exists():
        return task["id"], "skipped"
    artifacts_dir = Path(out) / "artifacts" if keep_artifacts else None
    row = run_real_task(task, refs, workroot, models=models, artifacts_dir=artifacts_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(row, default=float), encoding="utf-8")
    return task["id"], row.get("error") or f"{row['seconds']}s"


def _model_matrix(df: pd.DataFrame, arm: str, metric: str) -> pd.DataFrame:
    cols = {c[2:-len(f"_{arm}_{metric}")]: c for c in df
            if c.startswith("m_") and c.endswith(f"_{arm}_{metric}")}
    mat = pd.DataFrame({m: df[c] for m, c in cols.items()})
    mat.index = df["task"].str.replace("cxg__", "", regex=False).str[:45]
    return mat


def models_section(df: pd.DataFrame) -> list[str]:
    """Per-model table + paired comparisons over datasets (compare_units, Holm over the
    declared family). Datasets missing a score for either member of a pair are dropped
    from that pair only, and n is printed per pair."""
    from grnbench.eval.comparison import compare_units, holm_adjust

    if not any(c.startswith("m_") for c in df):
        return []
    lines = ["## Models (CollecTRI evaluation half; unit = dataset)", ""]
    for metric, label in (("auprc_ratio", "AUPRC ratio"), ("epr", "EPR")):
        auto, base = _model_matrix(df, "auto", metric), _model_matrix(df, "baseline", metric)
        order = auto.mean().sort_values(ascending=False).index
        lines += [f"### {label} (1 = random)", "",
                  "| model | n | mean auto | median auto | mean baseline | auto > baseline | s per dataset (auto) |",
                  "|---|---|---|---|---|---|---|"]
        secs = _model_matrix(df, "auto", "seconds")
        for m in order:
            a, b = auto[m], base.get(m, pd.Series(np.nan, index=auto.index))
            k = a.notna() & b.notna()
            lines.append(f"| {m} | {int(a.notna().sum())} | {a.mean():.3f} | {a.median():.3f} | {b[k].mean():.3f} | "
                         f"{int(((a - b)[k] > 1e-9).sum())}/{int(k.sum())} | {secs[m].mean():.0f} |")
        lines.append("")
        for arm, mat in (("auto", auto), ("baseline", base)):
            recs, raw = [], []
            for a_name, b_name in DECLARED_PAIRS:
                if a_name not in mat or b_name not in mat:
                    continue
                pair = mat[[a_name, b_name]].dropna()
                if len(pair) < 3:
                    continue
                res = compare_units(pair)["pairs"][0]
                recs.append((a_name, b_name, len(pair), res))
                raw += [res["p_a_gt_b"], res["p_b_gt_a"]]
            if not recs:
                continue
            adj = holm_adjust(raw)
            lines += [f"#### Declared comparisons, {label}, {arm} pipeline (paired t, Holm over "
                      f"{len(raw)} directed hypotheses)", "",
                      "| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |",
                      "|---|---|---|---|---|---|---|---|"]
            for i, (a_name, b_name, n, res) in enumerate(recs):
                pa, pb = adj[2 * i], adj[2 * i + 1]
                lo, hi = res["simultaneous_ci"]
                ci = "unbounded" if lo is None else f"[{lo:+.3f}, {hi:+.3f}]"
                win = a_name if pa < .05 else b_name if pb < .05 else "unresolved"
                lines.append(f"| {a_name} | {b_name} | {n} | {res['difference']:+.3f} | {ci} | "
                             f"{pa:.3g} | {pb:.3g} | {win} |")
            lines.append("")
    auto = _model_matrix(df, "auto", "auprc_ratio")
    if len(auto.columns):
        lines += ["### Per dataset, AUPRC ratio (auto pipeline)", "",
                  "| dataset | " + " | ".join(auto.columns) + " |",
                  "|---|" + "---|" * len(auto.columns)]
        for idx, r in auto.iterrows():
            lines.append(f"| {idx} | " + " | ".join("-" if pd.isna(v) else f"{v:.2f}" for v in r) + " |")
        lines.append("")
    errs = [c for c in df if c.startswith("m_") and c.endswith("_error")]
    if errs:
        lines += ["### Model errors", ""]
        for c in errs:
            for t, e in zip(df["task"], df[c]):
                if isinstance(e, str):
                    lines.append(f"- {t} / {c[2:-6]}: {e[:200]}")
        lines.append("")
    return lines


def report(out: Path) -> str:
    """REPORT_REAL.md: one row per dataset plus detection / robustness / recovery summaries."""
    rows = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((Path(out) / "rows").glob("*.json"))]
    df = pd.DataFrame(rows)
    df.to_csv(Path(out) / "all_rows_real.csv", index=False)
    n = len(df)
    err = df["error"].notna() if "error" in df else pd.Series(False, index=df.index)
    refused = (df["refused_as_expected"].map(lambda v: v is True) if "refused_as_expected" in df
               else pd.Series(False, index=df.index))
    lines = [f"# Engine mass test on real datasets - {n} datasets", "",
             f"- ran end to end: {int((~err & ~refused).sum())}",
             f"- refused as expected (not raw counts): {int(refused.sum())}",
             f"- errors: {int(err.sum())}", ""]
    det = [c for c in df if c.startswith("detect_")]
    if det:
        lines += ["## Detection vs curated metadata", "", "| check | n | correct |", "|---|---|---|"]
        for c in sorted(det):
            v = df[c].dropna().astype(bool)
            lines.append(f"| {c[7:]} | {len(v)} | {v.mean():.0%} |")
        lines.append("")
    m_cols = [c for c in df if c.startswith("permodel_") and c.endswith("_auto")]
    scored = df[df.get("reference_edges", pd.Series(0, index=df.index)).fillna(0) >= MIN_REFERENCE_EDGES]
    if len(scored) and m_cols:
        lines += ["## Reference recovery (CollecTRI evaluation half, AUPRC ratio; 1 = random)", "",
                  "| model | n | auto > baseline | auto = baseline | mean auto | mean baseline | mean rules |",
                  "|---|---|---|---|---|---|---|"]
        for m in EVAL_METHODS:
            a, b = scored.get(f"permodel_{m}_auto"), scored.get(f"permodel_{m}_baseline")
            if a is None:
                continue
            keep = a.notna() & b.notna()
            d = (a - b)[keep]
            lines.append(f"| {m} | {int(keep.sum())} | {(d > 1e-9).mean():.0%} | {(d.abs() <= 1e-9).mean():.0%} | "
                         f"{a[keep].mean():.2f} | {b[keep].mean():.2f} | "
                         f"{scored.loc[keep, f'permodel_{m}_rules'].mean():.2f} |")
        lines.append("")
    lines += models_section(df)
    lines += ["## Per dataset", "",
              "| dataset | truth | detected | ref edges | outcome | Pearson auto/base | CLR auto/base | s |",
              "|---|---|---|---|---|---|---|---|"]
    for _, r in df.iterrows():
        truth = json.loads(r["truth"]) if isinstance(r.get("truth"), str) else {}
        outcome = ("ERROR: " + str(r["error"])[:60]) if isinstance(r.get("error"), str) else             ("refused: not counts" if r.get("refused_as_expected") is True else
             (f"ok ({int(r['contexts_ok'])} contexts)" if pd.notna(r.get("contexts_ok")) else "ok, scored"))
        def pair(m):
            a, b = r.get(f"permodel_{m}_auto"), r.get(f"permodel_{m}_baseline")
            return "-" if a is None or pd.isna(a) else f"{a:.2f}/{b:.2f}"
        lines.append(f"| {r['task'][:55]} | {truth.get('assay', truth.get('technology', ''))} {truth.get('organism', '')} | "
                     f"{r.get('got_technology', '')} {r.get('got_organism', '')} | {r.get('reference_edges', '')} | "
                     f"{outcome} | {pair('pearson_signed')} | {pair('clr')} | {r.get('seconds', '')} |")
    text = "\n".join(lines) + "\n"
    (Path(out) / "REPORT_REAL.md").write_text(text, encoding="utf-8")
    return text


def final_report(out: Path, tasks_json: str | None = None) -> None:
    """Last pipeline step: figures (grnbench viz) and the LaTeX report (grnbench report), as the
    thesis pipeline emitted its final report. A failure here is logged, never fatal."""
    try:
        from grnbench.texreport import build_report
        from grnbench.viz import build_run

        print(build_run(out, tasks_json=tasks_json) / "index.md")
        print(build_report(str(out)))
    except Exception as exc:  # noqa: BLE001 - the numbers above are already written
        print(f"final report failed ({type(exc).__name__}: {exc}); REPORT_REAL.md is complete",
              file=sys.stderr)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("plan")
    a.add_argument("--data-root", required=True)
    a.add_argument("--out", required=True)
    b = sub.add_parser("run")
    b.add_argument("--tasks", required=True)
    b.add_argument("--refs", default=None)
    b.add_argument("--out", required=True)
    b.add_argument("--jobs", type=int, default=1)
    b.add_argument("--index", type=int, default=None)
    b.add_argument("--chunk", type=int, default=1)
    b.add_argument("--workroot", default=None)
    b.add_argument("--models", default=None,
                   help="comma-separated inference methods scored on the baseline and auto pipelines")
    b.add_argument("--keep-artifacts", action="store_true",
                   help="write edges, processed h5ad and decisions per task under OUT/artifacts/ "
                        "(default off; scores are unaffected either way)")
    c = sub.add_parser("report")
    c.add_argument("--out", required=True)
    c.add_argument("--tasks-for-viz", default=None, help="tasks.json with raw data paths (QC/embedding figures)")
    c.add_argument("--no-final-report", action="store_true", help="REPORT_REAL.md only (no figures, no LaTeX)")
    args = p.parse_args(argv)
    if args.cmd == "report":
        print(report(Path(args.out)))
        if not args.no_final_report:
            final_report(Path(args.out), args.tasks_for_viz)
        return 0
    if args.cmd == "plan":
        tasks = make_tasks(Path(args.data_root))
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(tasks, indent=0), encoding="utf-8")
        print(f"{len(tasks)} real-data tasks -> {args.out}")
        return 0
    tasks = json.loads(Path(args.tasks).read_text(encoding="utf-8"))
    if args.index is not None:
        tasks = tasks[args.index * args.chunk:(args.index + 1) * args.chunk]
    keep_awake()
    workroot = Path(args.workroot or tempfile.mkdtemp(prefix="gb_real_"))
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    models = [m.strip() for m in args.models.split(",") if m.strip()] if args.models else None
    jobs = [(t, args.refs, args.out, workroot, models, args.keep_artifacts) for t in tasks]
    if args.jobs > 1:
        from multiprocessing import get_context
        with get_context("spawn").Pool(args.jobs, maxtasksperchild=1) as pool:
            for tid, status in pool.imap_unordered(_run_and_save, jobs):
                print(f"{tid}: {status}", flush=True)
    else:
        for j in jobs:
            print("%s: %s" % _run_and_save(j), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
