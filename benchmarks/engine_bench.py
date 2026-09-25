"""Mass test of the cleaning engine on datasets with a known regulatory network.

For every dataset (simulated scenario, or a real dataset with a reference):

1. diagnose -> profile -> knowledge-base selection (what the engine detects);
2. screen the candidate pipelines - the rules' recommendation, one-change
   alternatives and a fixed **baseline** (standard scanpy-like cleaning) - on
   the *selection half* of the network;
3. score every candidate's cleaned data with fast inference methods on the
   *evaluation half* (never seen by the selection);
4. measure cleaning accuracy against ``obs['simulated_state']``.

Arms compared (paired per dataset): ``baseline``, ``rules`` (recommendation, no
screening), ``auto`` (screening with the selection half), ``auto_noref``
(screening without any reference) and ``oracle`` (best candidate on the
evaluation half - an upper bound for the candidate set).

    python benchmarks/engine_bench.py plan --presets-reps 10 --random 100 --out R/tasks.json
    python benchmarks/engine_bench.py run  --tasks R/tasks.json --out R --jobs 6
    python benchmarks/engine_bench.py run  --tasks R/tasks.json --out R --index $SLURM_ARRAY_TASK_ID --chunk 5
    python benchmarks/engine_bench.py report --out R
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

EVAL_METHODS = ["pearson_signed", "spearman_signed", "clr"]
BASELINE = {"cell_qc_thresholds": "fixed", "doublets": "none", "context_scope": "pooled",
            "normalization": "cpm_log1p", "batch": "none", "sparsity": "single_cells",
            "hvg": "tf_plus_1000", "gene_filter": "min_cells", "subsample": "none",
            "confounders": "none"}
BAD_STATES = ("damaged", "empty", "doublet", "raw_empty")
MAX_SCORE_CELLS = 3000


# -- task list ----------------------------------------------------------------------------- #

def make_tasks(presets_reps: int, n_random: int, presets=None) -> list[dict]:
    from grnbench.simulate import PRESETS
    names = presets or list(PRESETS)
    tasks = [{"id": f"{n}__s{s}", "kind": "preset", "name": n, "seed": s}
             for n in names for s in range(presets_reps)]
    tasks += [{"id": f"random__s{s}", "kind": "random", "name": "random", "seed": s}
              for s in range(n_random)]
    return tasks


# -- one task ----------------------------------------------------------------------------- #

def _scenario(task):
    from grnbench.simulate import PRESETS, random_scenario
    if task["kind"] == "random":
        return random_scenario(task["seed"])
    return PRESETS[task["name"]].with_seed(task["seed"])


def check_detection(diag, plan, expected: dict) -> dict:
    prof, sel, ch = plan.profile, plan.selection["results"], plan.choices
    got = {
        "technology": prof.get("technology"), "organism": prof.get("organism"),
        "batch_detected": prof.get("batch_key") is not None,
        "doublets": ch.get("doublets"), "doublets_removed": bool(plan.value("drop_doublets")),
        "normalization": ch.get("normalization"), "batch": ch.get("batch"),
        "identifier_mapping_best": sel["identifier_mapping"]["best_overall"],
        "cell_calling_best": sel["cell_calling"]["best_overall"],
        "context_scope": ch.get("context_scope"),
        "conflict": "technology" if "technology" in diag.context.conflicts() else None,
    }
    out = {}
    for key, want in expected.items():
        have = got.get(key)
        out[f"detect_{key}"] = bool(have in want) if isinstance(want, set) else bool(have == want)
    return out


def qc_accuracy(records, obs_all: pd.DataFrame, strategy: str) -> dict:
    import anndata as ad
    if strategy in ("metacells", "pseudobulk") or "simulated_state" not in obs_all:
        return {}
    kept, scope = set(), pd.Series(False, index=obs_all.index)
    for rec in records:
        if not rec.get("processed"):
            continue
        kept |= set(ad.read_h5ad(rec["processed"], backed="r").obs_names)
        if rec.get("subset_key"):
            scope |= obs_all[rec["subset_key"]].astype(str) == str(rec["subset_value"])
        else:
            scope[:] = True
    state = obs_all.loc[scope, "simulated_state"]
    is_kept = state.index.isin(list(kept))
    bad = state.isin(BAD_STATES).to_numpy()
    out = {"qc_ok_retained": float(is_kept[~bad].mean()) if (~bad).any() else np.nan}
    if bad.any():
        out["qc_bad_removed"] = float(1 - is_kept[bad].mean())
    for s in BAD_STATES:
        m = (state == s).to_numpy()
        if m.any():
            out[f"qc_removed_{s}"] = float(1 - is_kept[m].mean())
    return out


def score_candidates(records_by_cand, plans, tfs, eval_pairs) -> dict[str, dict]:
    import anndata as ad
    from grnbench.eval.metrics import GoldIndex, evaluate_ranking
    from grnbench.infer import run_method
    from grnbench.tfs import resolve_tf_names

    gi = GoldIndex.coerce(set(eval_pairs))
    scores = {}
    for name, recs in records_by_cand.items():
        strategy = plans[name].value("strategy")
        per_method = {m: [] for m in EVAL_METHODS}
        for rec in recs:
            if not rec.get("processed"):
                continue
            adata = ad.read_h5ad(rec["processed"])
            if adata.n_obs > MAX_SCORE_CELLS:
                keep = np.sort(np.random.default_rng(0).choice(adata.n_obs, MAX_SCORE_CELLS, replace=False))
                adata = adata[keep].copy()
            genes = pd.Index(adata.var_names).astype(str)
            sources = resolve_tf_names(genes, tfs)
            universe = pd.DataFrame([(s, g) for s in sources for g in genes if s != g],
                                    columns=["source", "target"])
            for m in EVAL_METHODS:
                try:
                    edges, _ = run_method(m, adata, tf_list=tfs, assume_normalized=strategy != "raw_counts")
                    r = evaluate_ranking(edges, gi, with_curve=False, universe=universe)
                    per_method[m].append(float(r.get("auprc_ratio", np.nan)))
                except Exception:  # noqa: BLE001 - a failing method is recorded as NaN
                    per_method[m].append(np.nan)
        row = {f"eval_{m}": float(np.nanmean(v)) if v and np.isfinite(v).any() else np.nan
               for m, v in per_method.items()}
        vals = [v for v in row.values() if np.isfinite(v)]
        row["eval_mean"] = float(np.mean(vals)) if vals else np.nan
        row["n_contexts"] = sum(1 for r in recs if r.get("processed"))
        scores[name] = row
    return scores


def run_task(task: dict, workroot: Path, keep: bool = False) -> dict:
    from grnbench.autoprep import diagnose, plan_preprocessing
    from grnbench.screen import rank_candidates, screen_candidates
    from grnbench.simulate import write_scenario

    t0 = time.time()
    row = {"task": task["id"], "scenario": task["name"], "seed": task["seed"], "error": None}
    work = Path(workroot) / task["id"]
    shutil.rmtree(work, ignore_errors=True)
    try:
        scn = _scenario(task)
        row.update({f"scn_{k}": v for k, v in scn.__dict__.items() if k != "name"})
        paths = write_scenario(work / "data", scn)
        tfs = paths["tfs"].read_text(encoding="utf-8").split()
        gold = pd.read_csv(paths["gold"], sep="\t")
        pairs = set(zip(gold["source"].astype(str), gold["target"].astype(str)))
        diag = diagnose(paths["counts"])
        plan = plan_preprocessing(diag, tf_list=tfs)
        row.update(check_detection(diag, plan, scn.expected()))
        row["rules_choices"] = json.dumps(plan.choices)
        row["n_warnings"] = len(plan.warnings)

        res = screen_candidates(diag, paths["counts"], work / "screen", tf_list=tfs,
                                reference_pairs=pairs, max_candidates=4, seed=task["seed"],
                                extra=[{"name": "baseline", "choices": BASELINE}],
                                proxy_methods=EVAL_METHODS)
        for name, recs in res["records"].items():
            for rec, ctx in zip(recs, res["plans"][name].contexts):
                rec["subset_key"], rec["subset_value"] = ctx.get("subset_key"), ctx.get("subset_value")
        errs = [r["error"] for recs in res["records"].values() for r in recs if r.get("error")]
        row["n_failed_contexts"] = len(errs)
        row["n_gate_refusals"] = sum("QC gate failure" in e for e in errs)
        if errs:
            row["first_context_error"] = errs[0][:300]
        scores = score_candidates(res["records"], res["plans"], tfs, res["reference"]["evaluation"])
        noref = res["rows"].drop(columns=[c for c in res["rows"] if c.startswith("reference_")])
        noref_rank, _ = rank_candidates(noref)
        arms = {"baseline": "baseline", "rules": "recommended", "auto": res["best"],
                "auto_noref": str(noref_rank.index[0]) if len(noref_rank) else None}
        finite = {k: v["eval_mean"] for k, v in scores.items() if np.isfinite(v["eval_mean"])}
        arms["oracle"] = max(finite, key=finite.get) if finite else None
        obs_all = diag.adata.obs
        for arm, cand in arms.items():
            row[f"{arm}_candidate"] = cand
            if cand is None or cand not in scores:
                continue
            for k, v in scores[cand].items():
                row[f"{arm}_{k}"] = v
            for k, v in qc_accuracy(res["records"][cand], obs_all,
                                    res["plans"][cand].value("strategy")).items():
                row[f"{arm}_{k}"] = v
            row[f"{arm}_failed_contexts"] = sum(1 for r in res["records"][cand] if r.get("error"))
        # Per model: the pipeline chosen on the selection half for THAT model,
        # scored with that model on the evaluation half.
        for m in EVAL_METHODS:
            best_m = res["best_per_method"].get(m)
            row[f"permodel_{m}_candidate"] = best_m
            row[f"permodel_{m}_auto"] = scores.get(best_m, {}).get(f"eval_{m}", np.nan)
            row[f"permodel_{m}_baseline"] = scores.get("baseline", {}).get(f"eval_{m}", np.nan)
            row[f"permodel_{m}_rules"] = scores.get("recommended", {}).get(f"eval_{m}", np.nan)
            vals = {k: v.get(f"eval_{m}", np.nan) for k, v in scores.items()}
            vals = {k: v for k, v in vals.items() if np.isfinite(v)}
            row[f"permodel_{m}_oracle"] = max(vals.values()) if vals else np.nan
        row["candidates"] = json.dumps({k: {kk: vv for kk, vv in v.items()} for k, v in scores.items()},
                                       default=float)
    except Exception as exc:  # noqa: BLE001 - robustness is part of what is measured
        row["error"] = f"{type(exc).__name__}: {exc}"
        row["traceback"] = traceback.format_exc()[-4000:]
    row["seconds"] = round(time.time() - t0, 1)
    if not keep:
        shutil.rmtree(work, ignore_errors=True)
    return row


def _run_and_save(args):
    task, out, workroot, keep = args
    path = Path(out) / "rows" / f"{task['id']}.json"
    if path.exists():
        return task["id"], "skipped"
    row = run_task(task, workroot, keep)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(row, default=float), encoding="utf-8")
    return task["id"], row.get("error") or f"{row['seconds']}s"


def keep_awake() -> None:
    """Ask Windows not to idle-sleep while this process runs (released when it exits).

    A per-process power request, like a video player's - no system setting changes.
    """
    if sys.platform != "win32":
        return
    import ctypes
    ES_CONTINUOUS, ES_SYSTEM_REQUIRED = 0x80000000, 0x00000001
    ctypes.windll.kernel32.SetThreadExecutionState(ES_CONTINUOUS | ES_SYSTEM_REQUIRED)


# -- report ------------------------------------------------------------------------------------ #

def report(out: Path) -> str:
    from scipy.stats import wilcoxon

    rows = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((out / "rows").glob("*.json"))]
    df = pd.DataFrame(rows)
    df.to_csv(out / "all_rows.csv", index=False)
    ok = df[df["error"].isna()] if "error" in df else df
    lines = [f"# Engine mass test - {len(df)} datasets, {len(ok)} completed, "
             f"{len(df) - len(ok)} errors", ""]
    arms = [a for a in ["baseline", "rules", "auto", "auto_noref", "oracle"] if f"{a}_eval_mean" in ok]

    def paired(a: str, b: str, frame: pd.DataFrame) -> dict:
        x, y = frame[f"{a}_eval_mean"], frame[f"{b}_eval_mean"]
        m = x.notna() & y.notna()
        x, y = x[m], y[m]
        if len(x) == 0:
            return {"n": 0}
        d = x - y
        res = {"n": int(len(x)), "win": float((d > 0).mean()), "tie": float((d == 0).mean()),
               "median_rel": float(np.median(d / y.abs().clip(lower=1e-9))),
               "mean_a": float(x.mean()), "mean_b": float(y.mean())}
        if (d != 0).sum() >= 6:
            res["p_wilcoxon"] = float(wilcoxon(x, y).pvalue)
        return res

    lines += ["## Overall (paired over datasets = independent units)", "",
              "| comparison | n | win rate | median relative gain | mean A | mean B | Wilcoxon p |",
              "|---|---|---|---|---|---|---|"]
    for a, b in [("auto", "baseline"), ("rules", "baseline"), ("auto_noref", "baseline"),
                 ("auto", "rules"), ("oracle", "auto")]:
        if a in arms and b in arms:
            r = paired(a, b, ok)
            lines.append(f"| {a} vs {b} | {r['n']} | {r.get('win', np.nan):.0%} | "
                         f"{r.get('median_rel', np.nan):+.1%} | {r.get('mean_a', np.nan):.2f} | "
                         f"{r.get('mean_b', np.nan):.2f} | {r.get('p_wilcoxon', np.nan):.2g} |")
    lines += ["", "## Per model: pipeline chosen for that model vs baseline (evaluation half)", "",
              "| model | n | auto wins | median relative gain | mean auto | mean baseline | "
              "mean rules | mean oracle | Wilcoxon p |", "|---|---|---|---|---|---|---|---|---|"]
    for m in EVAL_METHODS:
        a_col, b_col = f"permodel_{m}_auto", f"permodel_{m}_baseline"
        if a_col not in ok:
            continue
        x, y = ok[a_col], ok[b_col]
        keep = x.notna() & y.notna()
        x, y = x[keep], y[keep]
        if not len(x):
            continue
        d = x - y
        pval = wilcoxon(x, y).pvalue if (d != 0).sum() >= 6 else np.nan
        lines.append(f"| {m} | {len(x)} | {(d > 0).mean():.0%} | "
                     f"{np.median(d / y.abs().clip(lower=1e-9)):+.1%} | {x.mean():.2f} | {y.mean():.2f} | "
                     f"{ok.loc[keep, f'permodel_{m}_rules'].mean():.2f} | "
                     f"{ok.loc[keep, f'permodel_{m}_oracle'].mean():.2f} | {pval:.2g} |")
    lines += ["", "## Per scenario and model: auto-for-model vs baseline (mean AUPRC ratio)", "",
              "| scenario | n | " + " | ".join(f"{m} auto / base" for m in EVAL_METHODS) + " |",
              "|---|---|" + "---|" * len(EVAL_METHODS)]
    for scn, g in ok.groupby("scenario"):
        cells = []
        for m in EVAL_METHODS:
            if f"permodel_{m}_auto" in g:
                cells.append(f"{g[f'permodel_{m}_auto'].mean():.1f} / {g[f'permodel_{m}_baseline'].mean():.1f}")
            else:
                cells.append("-")
        lines.append(f"| {scn} | {len(g)} | " + " | ".join(cells) + " |")
    lines += ["", "## Per scenario: mean AUPRC ratio on the evaluation half (mean over models)", "",
              "| scenario | n paired/total | " + " | ".join(arms) + " | auto wins vs baseline |",
              "|---|---|" + "---|" * (len(arms) + 1)]
    for scn, g in ok.groupby("scenario"):
        r = paired("auto", "baseline", g) if {"auto", "baseline"} <= set(arms) else {}
        both = g[[f"{a}_eval_mean" for a in arms]].dropna()      # paired: every arm present
        lines.append(f"| {scn} | {len(both)}/{len(g)} | "
                     + " | ".join(f"{both[f'{a}_eval_mean'].mean():.2f}" for a in arms)
                     + f" | {r.get('win', np.nan):.0%} |")
    det = [c for c in ok if c.startswith("detect_")]
    if det:
        lines += ["", "## Detection accuracy (share of datasets where the engine got it right)", "",
                  "| check | n | accuracy | scenarios where it failed |", "|---|---|---|---|"]
        for c in sorted(det):
            v = ok[c].dropna().astype(bool)
            fails = sorted(set(ok.loc[ok[c] == False, "scenario"]))  # noqa: E712
            lines.append(f"| {c[7:]} | {len(v)} | {v.mean():.0%} | {', '.join(fails) or '-'} |")
    qc = [c for c in ok if c.startswith(("auto_qc_", "baseline_qc_"))]
    if qc:
        lines += ["", "## Cleaning accuracy (mean over datasets)", "", "| measure | baseline | auto |",
                  "|---|---|---|"]
        for m in sorted({c.split("_qc_", 1)[1] for c in qc}):
            b, a = ok.get(f"baseline_qc_{m}"), ok.get(f"auto_qc_{m}")
            lines.append(f"| {m} | {b.mean() if b is not None else np.nan:.1%} | "
                         f"{a.mean() if a is not None else np.nan:.1%} |")
    if len(df) - len(ok):
        lines += ["", "## Errors", ""]
        for _, r in df[df["error"].notna()].iterrows():
            lines.append(f"- `{r['task']}`: {r['error'][:300]}")
    lines += ["", f"Median runtime per dataset: {df['seconds'].median():.0f} s."]
    text = "\n".join(lines) + "\n"
    (out / "REPORT.md").write_text(text, encoding="utf-8")
    return text


# -- cli ----------------------------------------------------------------------------------------- #

def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("plan")
    a.add_argument("--presets-reps", type=int, default=10)
    a.add_argument("--random", type=int, default=100)
    a.add_argument("--presets", default=None, help="comma-separated subset of presets")
    a.add_argument("--out", required=True)
    b = sub.add_parser("run")
    b.add_argument("--tasks", required=True)
    b.add_argument("--out", required=True)
    b.add_argument("--jobs", type=int, default=1)
    b.add_argument("--index", type=int, default=None, help="SLURM array index (with --chunk)")
    b.add_argument("--chunk", type=int, default=1)
    b.add_argument("--workroot", default=None)
    b.add_argument("--keep", action="store_true")
    c = sub.add_parser("report")
    c.add_argument("--out", required=True)
    args = p.parse_args(argv)

    if args.cmd == "plan":
        tasks = make_tasks(args.presets_reps, args.random,
                           args.presets.split(",") if args.presets else None)
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(tasks, indent=0), encoding="utf-8")
        print(f"{len(tasks)} tasks -> {args.out}")
        return 0
    if args.cmd == "report":
        print(report(Path(args.out)))
        return 0

    tasks = json.loads(Path(args.tasks).read_text(encoding="utf-8"))
    if args.index is not None:
        tasks = tasks[args.index * args.chunk:(args.index + 1) * args.chunk]
    keep_awake()
    workroot = Path(args.workroot or tempfile.mkdtemp(prefix="gb_engine_"))
    jobs = [(t, args.out, workroot, args.keep) for t in tasks]
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    if args.jobs > 1:
        from multiprocessing import get_context
        with get_context("spawn").Pool(args.jobs, maxtasksperchild=1) as pool:
            for tid, status in pool.imap_unordered(_run_and_save, jobs):
                print(f"{tid}: {status}", flush=True)
    else:
        for j in jobs:
            tid, status = _run_and_save(j)
            print(f"{tid}: {status}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
