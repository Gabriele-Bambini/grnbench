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
import hashlib
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
DEFAULT_BENCH_MODELS = ["pearson_signed", "clr", "grnboost2", "regdiffusion"]
BENCH_METHOD_PARAMS = {"grnboost2": {"n_estimators": 100},
                       "regdiffusion": {"n_steps": 100}}
BASELINE = {"cell_qc_thresholds": "fixed", "doublets": "none", "context_scope": "pooled",
            "normalization": "cpm_log1p", "batch": "none", "sparsity": "single_cells",
            "hvg": "tf_plus_1000", "gene_filter": "min_cells", "subsample": "none",
            "confounders": "none"}
BAD_STATES = ("damaged", "empty", "doublet", "raw_empty")
MAX_SCORE_CELLS = 3000
TASK_SCHEMA_VERSION = 1
TERMINAL_STATUSES = {"completed", "expected_refusal", "permanent_error"}
RETRYABLE_STATUSES = {"pending", "running", "transient_error"}


# -- task list ----------------------------------------------------------------------------- #

def make_tasks(presets_reps: int, n_random: int, presets=None) -> list[dict]:
    from grnbench.simulate import PRESETS
    names = presets or list(PRESETS)
    tasks = [{"id": f"{n}__s{s}", "kind": "preset", "name": n, "seed": s}
             for n in names for s in range(presets_reps)]
    tasks += [{"id": f"random__s{s}", "kind": "random", "name": "random", "seed": s}
              for s in range(n_random)]
    return tasks


STRESS_FACTORS = {
    "depth": [
        ("mild", {"median_library": 1500.0}, 1500.0),
        ("moderate", {"median_library": 800.0}, 800.0),
        ("severe", {"median_library": 400.0}, 400.0),
        ("extreme", {"median_library": 200.0}, 200.0),
    ],
    "cells": [
        ("mild", {"n_cells": 1500}, 1500),
        ("moderate", {"n_cells": 800}, 800),
        ("severe", {"n_cells": 400}, 400),
        ("extreme", {"n_cells": 200}, 200),
    ],
    "doublets": [
        ("mild", {"doublet_frac": 0.03}, 0.03),
        ("moderate", {"doublet_frac": 0.06}, 0.06),
        ("severe", {"doublet_frac": 0.10}, 0.10),
        ("extreme", {"doublet_frac": 0.15}, 0.15),
    ],
    "ambient": [
        ("mild", {"ambient_frac": 0.05}, 0.05),
        ("moderate", {"ambient_frac": 0.10}, 0.10),
        ("severe", {"ambient_frac": 0.20}, 0.20),
        ("extreme", {"ambient_frac": 0.30}, 0.30),
    ],
    "damaged": [
        ("mild", {"damaged_frac": 0.03}, 0.03),
        ("moderate", {"damaged_frac": 0.06}, 0.06),
        ("severe", {"damaged_frac": 0.10}, 0.10),
        ("extreme", {"damaged_frac": 0.20}, 0.20),
    ],
    "batch_confounding": [
        ("mild", {"n_batches": 3, "batch_shift_sd": 0.4, "batch_confounding": 0.25}, 0.25),
        ("moderate", {"n_batches": 3, "batch_shift_sd": 0.4, "batch_confounding": 0.50}, 0.50),
        ("severe", {"n_batches": 3, "batch_shift_sd": 0.4, "batch_confounding": 0.75}, 0.75),
        ("extreme", {"n_batches": 3, "batch_shift_sd": 0.4, "batch_confounding": 1.00}, 1.00),
    ],
    "batch_depth": [
        ("mild", {"n_batches": 2, "batch_depth_ratio": 1.25}, 1.25),
        ("moderate", {"n_batches": 2, "batch_depth_ratio": 1.50}, 1.50),
        ("severe", {"n_batches": 2, "batch_depth_ratio": 2.00}, 2.00),
        ("extreme", {"n_batches": 2, "batch_depth_ratio": 3.00}, 3.00),
    ],
    "type_imbalance": [
        ("mild", {"type_balance_ratio": 2.0}, 2.0),
        ("moderate", {"type_balance_ratio": 5.0}, 5.0),
        ("severe", {"type_balance_ratio": 10.0}, 10.0),
        ("extreme", {"type_balance_ratio": 20.0}, 20.0),
    ],
}


def make_stress_tasks(n_seeds: int = 5,
                      shared_fractions=(1.0, 0.75, 0.5, 0.25)) -> list[dict]:
    """Controlled 740-task grid: clean + 8 one-factor curves + 4 mixed curves.

    Every biological regime has 37 scenarios and every scenario has ``n_seeds``
    independent replicates.  The explicit task metadata drives the generated
    breaking-point table; scenario names are never parsed to recover design.
    """
    tasks = []
    for shared in shared_fractions:
        sf = f"sf{int(round(100 * shared)):03d}"
        base = {"n_cells": 3000, "median_library": 3000.0,
                "type_balance_ratio": 1.0, "shared_edge_fraction": float(shared)}
        designs = [("clean", 0, "baseline", 0.0, {})]
        for stressor, levels in STRESS_FACTORS.items():
            for severity, (label, changes, value) in enumerate(levels, start=1):
                designs.append((stressor, severity, label, value, changes))
        for severity in range(1, 5):
            changes = {}
            values = {}
            for stressor, levels in STRESS_FACTORS.items():
                _label, change, value = levels[severity - 1]
                changes.update(change)
                values[stressor] = value
            designs.append(("mixed", severity, f"mixed_{severity}", severity,
                            {**changes, "_mixed_values": values}))
        for stressor, severity, label, value, raw_changes in designs:
            changes = dict(raw_changes)
            mixed_values = changes.pop("_mixed_values", None)
            scenario = {**base, **changes}
            for seed in range(n_seeds):
                task = {
                    "id": f"{sf}__{stressor}__l{severity}__s{seed}",
                    "kind": "stress",
                    "name": f"{sf}_{stressor}_{label}",
                    "seed": seed,
                    "shared_edge_fraction": float(shared),
                    "stressor": stressor,
                    "severity": severity,
                    "severity_label": label,
                    "severity_value": value,
                    "scenario": scenario,
                }
                if mixed_values is not None:
                    task["mixed_values"] = mixed_values
                tasks.append(task)
    return tasks


def _canonical_tasks(tasks: list[dict]) -> bytes:
    return json.dumps(tasks, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()


def task_digest(tasks: list[dict]) -> str:
    return hashlib.sha256(_canonical_tasks(tasks)).hexdigest()


def write_task_manifest(path: str | Path, tasks: list[dict], *, task_type: str) -> dict:
    ids = [str(t.get("id")) for t in tasks]
    if len(ids) != len(set(ids)):
        raise ValueError("task ids must be unique")
    payload = {"schema_version": TASK_SCHEMA_VERSION, "task_type": task_type,
               "task_count": len(tasks), "tasks_sha256": task_digest(tasks), "tasks": tasks}
    _atomic_json(Path(path), payload)
    return payload


def load_task_manifest(path: str | Path) -> tuple[list[dict], dict]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if isinstance(payload, list):                 # legacy manifests remain runnable
        tasks = payload
        meta = {"schema_version": 0, "task_type": "legacy", "task_count": len(tasks),
                "tasks_sha256": task_digest(tasks)}
    elif isinstance(payload, dict) and isinstance(payload.get("tasks"), list):
        tasks = payload["tasks"]
        meta = {k: v for k, v in payload.items() if k != "tasks"}
        if int(meta.get("task_count", -1)) != len(tasks):
            raise ValueError(f"{path}: task_count does not match tasks")
        actual = task_digest(tasks)
        if meta.get("tasks_sha256") != actual:
            raise ValueError(f"{path}: task manifest hash mismatch")
    else:
        raise ValueError(f"{path}: expected a task list or versioned task manifest")
    ids = [str(t.get("id")) for t in tasks]
    if len(ids) != len(set(ids)) or any(i in ("", "None") for i in ids):
        raise ValueError(f"{path}: task ids must be present and unique")
    return tasks, meta


def ensure_run_config(out: str | Path, manifest_meta: dict, *, runner: str,
                      scientific_config: dict) -> dict:
    """Pin score-affecting run options on first invocation and reject drift."""
    payload = {"schema_version": 1, "runner": runner,
               "tasks_sha256": manifest_meta["tasks_sha256"],
               "scientific_config": scientific_config}
    payload["config_sha256"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    path = Path(out) / "RUN_CONFIG.json"
    if path.exists():
        previous = json.loads(path.read_text(encoding="utf-8"))
        if previous != payload:
            raise ValueError(
                f"{path}: run configuration differs from the existing immutable configuration; "
                "use a new output directory"
            )
    else:
        _atomic_json(path, payload)
    return payload


# -- one task ----------------------------------------------------------------------------- #

def _scenario(task):
    from grnbench.simulate import PRESETS, Scenario, random_scenario
    if task["kind"] == "stress":
        return Scenario(name=task["name"], seed=int(task["seed"]), **task["scenario"])
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


def _truth_groups(rec: dict, eval_pairs, eval_pairs_by_context: dict[str, set] | None):
    if not eval_pairs_by_context:
        return {"all": set(eval_pairs)}
    value = rec.get("subset_value")
    if value is not None and str(value) in eval_pairs_by_context:
        return {str(value): eval_pairs_by_context[str(value)]}
    return eval_pairs_by_context


def score_candidates(records_by_cand, plans, tfs, eval_pairs,
                     eval_pairs_by_context: dict[str, set] | None = None,
                     methods: list[str] | None = None,
                     method_params: dict[str, dict] | None = None) -> dict[str, dict]:
    import anndata as ad
    from grnbench.eval.metrics import GoldIndex, evaluate_ranking
    from grnbench.infer import run_method
    from grnbench.tfs import resolve_tf_names

    methods = methods or EVAL_METHODS
    method_params = method_params or {}
    scores = {}
    for name, recs in records_by_cand.items():
        strategy = plans[name].value("strategy")
        per_method = {m: {"auprc_ratio": [], "epr": []} for m in methods}
        per_method_errors = {m: [] for m in methods}
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
            groups = _truth_groups(rec, eval_pairs, eval_pairs_by_context)
            for m in methods:
                try:
                    edges, _ = run_method(m, adata, tf_list=tfs,
                                          assume_normalized=strategy != "raw_counts",
                                          **method_params.get(m, {}))
                    for pairs in groups.values():
                        gi = GoldIndex.coerce(set(pairs))
                        r = evaluate_ranking(edges, gi, with_curve=False, universe=universe)
                        for metric in ("auprc_ratio", "epr"):
                            per_method[m][metric].append(float(r.get(metric, np.nan)))
                except Exception as exc:  # noqa: BLE001 - a failing method is recorded as NaN
                    per_method[m]["auprc_ratio"].append(np.nan)
                    per_method[m]["epr"].append(np.nan)
                    per_method_errors[m].append(f"{type(exc).__name__}: {exc}"[:300])
        row = {}
        for method, metrics in per_method.items():
            for metric, values in metrics.items():
                finite = [v for v in values if np.isfinite(v)]
                key = f"eval_{method}" if metric == "auprc_ratio" else f"eval_{method}_{metric}"
                row[key] = float(np.mean(finite)) if finite else np.nan
            if per_method_errors[method]:
                row[f"eval_{method}_error"] = " | ".join(sorted(set(per_method_errors[method])))[:1000]
        vals = [row[f"eval_{m}"] for m in methods if np.isfinite(row[f"eval_{m}"])]
        row["eval_mean"] = float(np.mean(vals)) if vals else np.nan
        row["n_contexts"] = sum(len(_truth_groups(r, eval_pairs, eval_pairs_by_context))
                                for r in recs if r.get("processed"))
        scores[name] = row
    return scores


def score_selected_models(records_by_cand, plans, arms: dict[str, str | None], tfs,
                          eval_pairs, eval_pairs_by_context, methods: list[str]) -> dict:
    """Score requested models only on distinct candidates used by an arm.

    Screening remains deliberately cheap and chooses one common ``auto``
    pipeline from the selection half.  Expensive models are then evaluated on
    the held-out half for baseline/rules/auto without multiplying their cost by
    every rejected candidate.
    """
    chosen = {cand: records_by_cand[cand] for cand in set(arms.values())
              if cand is not None and cand in records_by_cand}
    scores = score_candidates(chosen, plans, tfs, eval_pairs, eval_pairs_by_context,
                              methods=methods, method_params=BENCH_METHOD_PARAMS)
    out = {}
    for arm, cand in arms.items():
        if cand not in scores:
            continue
        for method in methods:
            out[f"m_{method}_{arm}_auprc_ratio"] = scores[cand].get(f"eval_{method}", np.nan)
            out[f"m_{method}_{arm}_epr"] = scores[cand].get(f"eval_{method}_epr", np.nan)
            if scores[cand].get(f"eval_{method}_error"):
                out[f"m_{method}_{arm}_error"] = scores[cand][f"eval_{method}_error"]
    return out


def run_task(task: dict, workroot: Path, keep: bool = False,
             models: list[str] | None = None, attempt: int = 1) -> dict:
    from grnbench.autoprep import diagnose, plan_preprocessing
    from grnbench.screen import rank_candidates, screen_candidates
    from grnbench.simulate import write_scenario

    t0 = time.time()
    row = {"task": task["id"], "scenario": task["name"], "seed": task["seed"], "error": None,
           "status": "running", "attempt": int(attempt)}
    for key in ("shared_edge_fraction", "stressor", "severity", "severity_label",
                "severity_value", "mixed_values"):
        if key in task:
            row[key] = task[key]
    work = Path(workroot) / task["id"]
    shutil.rmtree(work, ignore_errors=True)
    try:
        scn = _scenario(task)
        row.update({f"scn_{k}": v for k, v in scn.__dict__.items() if k != "name"})
        paths = write_scenario(work / "data", scn)
        tfs = paths["tfs"].read_text(encoding="utf-8").split()
        gold = pd.read_csv(paths["gold"], sep="\t")
        pairs = set(zip(gold["source"].astype(str), gold["target"].astype(str)))
        pairs_by_context = {}
        for context, context_path in paths.get("gold_by_context", {}).items():
            frame = pd.read_csv(context_path, sep="\t")
            pairs_by_context[str(context)] = set(zip(frame["source"].astype(str),
                                                      frame["target"].astype(str)))
        diag = diagnose(paths["counts"])
        plan = plan_preprocessing(diag, tf_list=tfs)
        row.update(check_detection(diag, plan, scn.expected()))
        row["rules_choices"] = json.dumps(plan.choices)
        row["n_warnings"] = len(plan.warnings)

        res = screen_candidates(diag, paths["counts"], work / "screen", tf_list=tfs,
                                reference_pairs=pairs, max_candidates=4, seed=task["seed"],
                                extra=[{"name": "baseline", "choices": BASELINE}],
                                proxy_methods=EVAL_METHODS,
                                reference_pairs_by_context=pairs_by_context)
        for name, recs in res["records"].items():
            for rec, ctx in zip(recs, res["plans"][name].contexts):
                rec["subset_key"], rec["subset_value"] = ctx.get("subset_key"), ctx.get("subset_value")
        errs = [r["error"] for recs in res["records"].values() for r in recs if r.get("error")]
        row["n_failed_contexts"] = len(errs)
        row["n_gate_refusals"] = sum("QC gate failure" in e for e in errs)
        if errs:
            row["first_context_error"] = errs[0][:300]
        scores = score_candidates(res["records"], res["plans"], tfs,
                                  res["reference"]["evaluation"],
                                  res["reference"].get("evaluation_by_context"))
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
        if models:
            row.update(score_selected_models(
                res["records"], res["plans"],
                {k: arms[k] for k in ("baseline", "rules", "auto")},
                tfs, res["reference"]["evaluation"],
                res["reference"].get("evaluation_by_context"), models,
            ))
        row["candidates"] = json.dumps({k: {kk: vv for kk, vv in v.items()} for k, v in scores.items()},
                                       default=float)
        row["status"] = "completed"
    except Exception as exc:  # noqa: BLE001 - robustness is part of what is measured
        row["error"] = f"{type(exc).__name__}: {exc}"
        row["traceback"] = traceback.format_exc()[-4000:]
        row["status"] = "transient_error"
    row["seconds"] = round(time.time() - t0, 1)
    if not keep:
        shutil.rmtree(work, ignore_errors=True)
    return row


def _atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, default=float)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def _normalise_status(row: dict) -> str:
    status = row.get("status")
    if status:
        return str(status)
    if row.get("refused_as_expected") is True:
        return "expected_refusal"
    return "completed" if not row.get("error") else "transient_error"


def _run_and_save(args):
    task, out, workroot, keep, models, max_attempts, retry_permanent = args
    path = Path(out) / "rows" / f"{task['id']}.json"
    previous = None
    if path.exists():
        try:
            previous = json.loads(path.read_text(encoding="utf-8"))
            status = _normalise_status(previous)
            if status in TERMINAL_STATUSES and not (status == "permanent_error" and retry_permanent):
                return task["id"], status, "skipped terminal result"
        except (OSError, ValueError):
            previous = {"attempt": 0, "error": "unreadable previous row"}
    attempt = int((previous or {}).get("attempt", 0)) + 1
    row = run_task(task, workroot, keep, models=models, attempt=attempt)
    if row["status"] == "transient_error" and attempt >= max_attempts:
        row["status"] = "permanent_error"
    prior_errors = list((previous or {}).get("prior_errors", []))
    if previous and previous.get("error"):
        prior_errors.append({"attempt": previous.get("attempt"), "error": previous["error"]})
    if prior_errors:
        row["prior_errors"] = prior_errors
    _atomic_json(path, row)
    if keep:
        source = Path(workroot) / task["id"]
        target = Path(out) / "artifacts" / task["id"]
        if source.exists():
            shutil.rmtree(target, ignore_errors=True)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(source), str(target))
    return task["id"], row["status"], row.get("error") or f"{row['seconds']}s"


def pool_results(fn, jobs, n_workers: int):
    """Yield ``fn(job)`` for every job from spawned workers, one task per worker process.

    ProcessPoolExecutor, not multiprocessing.Pool: when the OOM killer takes a worker, Pool waits
    for it forever (EXP-20261002-LM-2 job 679448 sat ~20 h until its time limit), while the
    executor raises ``BrokenProcessPool`` at once. Rows already written are kept, so the caller
    reports the run as retryable and a resubmission resumes it.
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed
    from multiprocessing import get_context

    with ProcessPoolExecutor(n_workers, mp_context=get_context("spawn"), max_tasks_per_child=1) as pool:
        for fut in as_completed([pool.submit(fn, j) for j in jobs]):
            yield fut.result()


def keep_awake() -> None:
    """Ask Windows not to idle-sleep while this process runs (released when it exits).

    A per-process power request, like a video player's - no system setting changes.
    """
    if sys.platform != "win32":
        return
    import ctypes
    ES_CONTINUOUS, ES_SYSTEM_REQUIRED = 0x80000000, 0x00000001
    ctypes.windll.kernel32.SetThreadExecutionState(ES_CONTINUOUS | ES_SYSTEM_REQUIRED)


def is_sentinel_task(task: dict) -> bool:
    """Pre-declared small artifact set: clean, extreme depth and extreme mixed, seed 0."""
    return (task.get("kind") == "stress" and int(task.get("seed", -1)) == 0
            and (task.get("stressor") == "clean"
                 or (int(task.get("severity", -1)) == 4
                     and task.get("stressor") in {"depth", "mixed"})))


# -- report ------------------------------------------------------------------------------------ #

def audit_completeness(out: Path, tasks: list[dict], manifest_meta: dict) -> dict:
    planned = [str(t["id"]) for t in tasks]
    rows, invalid_files = [], []
    for path in sorted((out / "rows").glob("*.json")):
        try:
            row = json.loads(path.read_text(encoding="utf-8"))
            row["_row_file"] = path.name
            rows.append(row)
        except (OSError, ValueError) as exc:
            invalid_files.append({"file": path.name, "error": f"{type(exc).__name__}: {exc}"})
    by_id: dict[str, list[dict]] = {}
    for row in rows:
        by_id.setdefault(str(row.get("task")), []).append(row)
    missing = [task for task in planned if task not in by_id]
    duplicate_rows = sorted(task for task, group in by_id.items() if len(group) > 1)
    extra = sorted(task for task in by_id if task not in set(planned))
    status_by_task = {task: _normalise_status(by_id[task][0])
                      for task in planned if task in by_id and len(by_id[task]) == 1}
    incomplete = sorted(task for task, status in status_by_task.items()
                        if status not in TERMINAL_STATUSES)
    counts = {status: sum(s == status for s in status_by_task.values())
              for status in sorted(TERMINAL_STATUSES | RETRYABLE_STATUSES)}
    config_errors = []
    config_path = out / "RUN_CONFIG.json"
    if int(manifest_meta.get("schema_version", 0)) >= 1:
        try:
            config = json.loads(config_path.read_text(encoding="utf-8"))
            if config.get("tasks_sha256") != manifest_meta.get("tasks_sha256"):
                config_errors.append("RUN_CONFIG tasks_sha256 differs from task manifest")
            claimed = config.get("config_sha256")
            unhashed = {k: v for k, v in config.items() if k != "config_sha256"}
            actual = hashlib.sha256(
                json.dumps(unhashed, sort_keys=True, separators=(",", ":")).encode()
            ).hexdigest()
            if claimed != actual:
                config_errors.append("RUN_CONFIG config_sha256 mismatch")
        except (OSError, ValueError) as exc:
            config_errors.append(f"RUN_CONFIG unavailable: {type(exc).__name__}: {exc}")
    complete = not (missing or incomplete or duplicate_rows or extra or invalid_files or config_errors)
    return {
        "schema_version": 1,
        "complete": complete,
        "task_count": len(planned),
        "tasks_sha256": manifest_meta.get("tasks_sha256", task_digest(tasks)),
        "status_counts": counts,
        "missing": missing,
        "incomplete": incomplete,
        "duplicate_rows": duplicate_rows,
        "extra_rows": extra,
        "invalid_row_files": invalid_files,
        "run_config_errors": config_errors,
        "unexpected_errors": [task for task, status in status_by_task.items()
                              if status == "permanent_error"],
    }


def enforce_completeness(out: Path, tasks_path: str | Path, *, partial: bool = False) -> dict:
    tasks, meta = load_task_manifest(tasks_path)
    audit = audit_completeness(out, tasks, meta)
    complete_path, partial_path = out / "COMPLETE.json", out / "PARTIAL.json"
    if audit["complete"]:
        _atomic_json(complete_path, audit)
        partial_path.unlink(missing_ok=True)
    elif partial:
        _atomic_json(partial_path, audit)
        complete_path.unlink(missing_ok=True)
    else:
        raise RuntimeError(
            "run is incomplete: "
            f"missing={len(audit['missing'])}, retryable={len(audit['incomplete'])}, "
            f"duplicates={len(audit['duplicate_rows'])}, extra={len(audit['extra_rows'])}, "
            f"invalid={len(audit['invalid_row_files'])}, config={len(audit['run_config_errors'])}; "
            "rerun tasks or pass --partial for a watermarked report"
        )
    return audit


def _mean_ci(values) -> tuple[float, float, float]:
    from scipy.stats import t
    x = np.asarray([v for v in values if np.isfinite(v)], dtype=float)
    if len(x) == 0:
        return np.nan, np.nan, np.nan
    mean = float(x.mean())
    if len(x) < 2 or float(x.std(ddof=1)) == 0:
        return mean, mean, mean
    half = float(t.ppf(0.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x)))
    return mean, mean - half, mean + half


def _write_stress_figures(levels: pd.DataFrame, out: Path) -> Path:
    """Plot every controlled stress curve; values come only from ``stress_levels.csv``."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    viz = out / "viz"
    viz.mkdir(exist_ok=True)
    made = []
    regimes = sorted(levels["shared_edge_fraction"].dropna().unique(), reverse=True)
    for stressor, stress in levels.groupby("stressor", sort=True):
        fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), sharex=True, sharey=True)
        for ax, shared in zip(axes.flat, regimes):
            panel = stress[stress["shared_edge_fraction"] == shared]
            for model, curve in panel.groupby("model", sort=True):
                curve = curve.sort_values("severity")
                x = np.r_[0, curve["severity"].to_numpy(dtype=float)]
                y = np.r_[curve["clean_mean_auprc_ratio"].iloc[0],
                          curve["mean_auprc_ratio"].to_numpy(dtype=float)]
                lo = np.r_[curve["clean_ci95_low"].iloc[0],
                           curve["ci95_low"].to_numpy(dtype=float)]
                hi = np.r_[curve["clean_ci95_high"].iloc[0],
                           curve["ci95_high"].to_numpy(dtype=float)]
                ax.plot(x, y, marker="o", label=model)
                ax.fill_between(x, lo, hi, alpha=0.12)
            labels = ["clean", *curve["severity_label"].astype(str).tolist()]
            ax.set_xticks(range(len(labels)), labels, rotation=25, ha="right")
            ax.axhline(1.0, color="black", linestyle="--", linewidth=0.8)
            ax.set_title(f"shared edges = {shared:.0%}")
            ax.set_ylabel("AUPRC ratio (1 = random)")
            ax.grid(alpha=0.2)
        handles, labels = axes.flat[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", ncol=max(1, len(labels)))
        fig.suptitle(f"Controlled stress: {stressor} (mean and 95% t interval; n = 5 seeds)", y=0.995)
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        name = "".join(c if c.isalnum() or c in "-_" else "_" for c in str(stressor))
        path = viz / f"stress_{name}.png"
        fig.savefig(path, dpi=120)
        plt.close(fig)
        made.append(path.name)
    lines = ["# Synthetic stress-curve figures", "",
             "Generated from `stress_levels.csv`; each point is the mean across five simulator seeds, "
             "and bands are 95% Student-t intervals. AUPRC ratio 1 is random.", ""]
    lines.extend(f"- [{name}]({name})" for name in made)
    lines.append("")
    (viz / "index.md").write_text("\n".join(lines), encoding="utf-8")
    return viz


def stress_breaking_points(df: pd.DataFrame, out: Path) -> list[str]:
    """Generate descriptive breakpoints for every model/regime/stressor.

    These are diagnostic thresholds, not confirmatory hypothesis tests: five
    simulator seeds estimate a 95% t interval; a later registered experiment
    must confirm any claim on real datasets.
    """
    if "stressor" not in df:
        return []
    model_cols = {c[len("m_"):-len("_auto_auprc_ratio")]: c for c in df
                  if c.startswith("m_") and c.endswith("_auto_auprc_ratio")}
    if not model_cols:
        return []
    records, summaries = [], []
    for shared, regime in df.groupby("shared_edge_fraction", dropna=False):
        clean = regime[regime["stressor"] == "clean"]
        for model, col in model_cols.items():
            clean_mean, clean_lo, clean_hi = _mean_ci(clean[col].dropna())
            for stressor, curve in regime[regime["stressor"] != "clean"].groupby("stressor"):
                levels = []
                for severity, level in curve.groupby("severity"):
                    mean, lo, hi = _mean_ci(level[col].dropna())
                    finite_rate = float(level[col].notna().mean())
                    rec = {"shared_edge_fraction": shared, "model": model,
                           "stressor": stressor, "severity": int(severity),
                           "severity_label": str(level["severity_label"].iloc[0]),
                           "severity_value": level["severity_value"].iloc[0],
                           "n": int(level[col].notna().sum()), "mean_auprc_ratio": mean,
                           "ci95_low": lo, "ci95_high": hi,
                           "clean_mean_auprc_ratio": clean_mean,
                           "clean_ci95_low": clean_lo, "clean_ci95_high": clean_hi,
                           "relative_to_clean": (mean / clean_mean
                                                 if np.isfinite(mean) and clean_mean > 0 else np.nan),
                           "finite_rate": finite_rate,
                           "above_random_95ci": bool(np.isfinite(lo) and lo > 1.0)}
                    records.append(rec)
                    levels.append(rec)
                statistical = ("clean" if not (np.isfinite(clean_lo) and clean_lo > 1.0) else
                               next((r["severity_label"] for r in levels
                                     if not r["above_random_95ci"]), "not reached"))
                relative = next((r["severity_label"] for r in levels
                                 if np.isfinite(r["relative_to_clean"])
                                 and r["relative_to_clean"] <= 0.5), "not reached")
                operational = next((r["severity_label"] for r in levels
                                    if r["finite_rate"] < 1.0), "not reached")
                summaries.append({"shared_edge_fraction": shared, "model": model,
                                  "stressor": stressor, "statistical_break": statistical,
                                  "half_clean_break": relative,
                                  "operational_break": operational})
    pd.DataFrame(records).to_csv(out / "stress_levels.csv", index=False)
    pd.DataFrame(summaries).to_csv(out / "breaking_points.csv", index=False)
    _write_stress_figures(pd.DataFrame(records), out)
    lines = ["", "## Controlled stress-test breaking points", "",
             "Exploratory diagnostics: `statistical` is the first level whose 95% t interval no longer "
             "lies wholly above random (AUPRC ratio 1); `half-clean` is the first level at or below "
             "50% of the clean mean; `operational` is the first level with any missing model score.", "",
             "| shared edges | model | stressor | statistical | half-clean | operational |",
             "|---|---|---|---|---|---|"]
    for row in summaries:
        lines.append(f"| {row['shared_edge_fraction']:.0%} | {row['model']} | {row['stressor']} | "
                     f"{row['statistical_break']} | {row['half_clean_break']} | "
                     f"{row['operational_break']} |")
    return lines + [""]


def synthetic_models_section(df: pd.DataFrame) -> list[str]:
    models = sorted({c[len("m_"):-len("_auto_auprc_ratio")] for c in df
                     if c.startswith("m_") and c.endswith("_auto_auprc_ratio")})
    if not models:
        return []
    lines = ["", "## Requested inference models", "",
             "The same held-out context truths are used for all arms. The common `auto` pipeline is "
             "chosen by the fast proxy methods on the selection half; expensive models never see "
             "the evaluation half during pipeline selection.", ""]
    for metric, label in (("auprc_ratio", "AUPRC ratio"), ("epr", "EPR")):
        lines += [f"### {label} (1 = random)", "",
                  "| model | n auto | mean baseline | mean rules | mean auto | auto > baseline |",
                  "|---|---|---|---|---|---|"]
        for model in models:
            auto = df.get(f"m_{model}_auto_{metric}", pd.Series(np.nan, index=df.index))
            base = df.get(f"m_{model}_baseline_{metric}", pd.Series(np.nan, index=df.index))
            rules = df.get(f"m_{model}_rules_{metric}", pd.Series(np.nan, index=df.index))
            paired = auto.notna() & base.notna()
            lines.append(f"| {model} | {int(auto.notna().sum())} | {base.mean():.3f} | "
                         f"{rules.mean():.3f} | {auto.mean():.3f} | "
                         f"{int((auto[paired] > base[paired]).sum())}/{int(paired.sum())} |")
        lines.append("")
    errors = [c for c in df if c.startswith("m_") and c.endswith("_error")]
    if errors:
        lines += ["### Model errors", ""]
        for col in errors:
            for task, error in zip(df["task"], df[col]):
                if isinstance(error, str):
                    lines.append(f"- `{task}` / `{col[2:-6]}`: {error[:300]}")
        lines.append("")
    return lines


def report(out: Path, tasks_path: str | Path | None = None, *, partial: bool = False) -> str:
    from scipy.stats import wilcoxon

    tasks_path = Path(tasks_path or out / "tasks.json")
    audit = enforce_completeness(out, tasks_path, partial=partial)
    rows = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((out / "rows").glob("*.json"))]
    df = pd.DataFrame(rows)
    df.to_csv(out / "all_rows.csv", index=False)
    status = df.apply(lambda r: _normalise_status(r.to_dict()), axis=1)
    ok = df[status == "completed"]
    lines = [f"# {'PARTIAL - ' if partial and not audit['complete'] else ''}Engine mass test - "
             f"{audit['task_count']} planned datasets", "",
             f"- completed: {audit['status_counts'].get('completed', 0)}",
             f"- expected refusals: {audit['status_counts'].get('expected_refusal', 0)}",
             f"- permanent errors: {audit['status_counts'].get('permanent_error', 0)}",
             f"- missing or retryable: {len(audit['missing']) + len(audit['incomplete'])}", ""]
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
    lines += synthetic_models_section(ok)
    lines += stress_breaking_points(df, out)
    if audit["unexpected_errors"]:
        lines += ["", "## Errors", ""]
        for _, r in df[df["task"].isin(audit["unexpected_errors"])].iterrows():
            lines.append(f"- `{r['task']}`: {r['error'][:300]}")
    lines += ["", f"Median runtime per completed dataset: {ok['seconds'].median():.0f} s."]
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
    a.add_argument("--stress-grid", action="store_true",
                   help="plan the controlled 4-regime x 37-scenario stress grid")
    a.add_argument("--stress-seeds", type=int, default=5)
    a.add_argument("--out", required=True)
    b = sub.add_parser("run")
    b.add_argument("--tasks", required=True)
    b.add_argument("--out", required=True)
    b.add_argument("--jobs", type=int, default=1)
    b.add_argument("--index", type=int, default=None, help="SLURM array index (with --chunk)")
    b.add_argument("--chunk", type=int, default=1)
    b.add_argument("--workroot", default=None)
    b.add_argument("--keep", action="store_true", help="keep work artifacts for every task")
    b.add_argument("--keep-artifacts", action="store_true",
                   help="keep the 12 pre-declared stress sentinels (clean, extreme depth and "
                        "extreme mixed at seed 0 in each network regime)")
    b.add_argument("--models", default=None,
                   help="comma-separated models evaluated on baseline/rules/auto; "
                        "use 'default' for pearson, CLR, GRNBoost2 and RegDiffusion")
    b.add_argument("--max-attempts", type=int, default=2,
                   help="row failures become permanent after this many invocations")
    b.add_argument("--retry-permanent", action="store_true")
    c = sub.add_parser("report")
    c.add_argument("--out", required=True)
    c.add_argument("--tasks", default=None)
    c.add_argument("--partial", action="store_true",
                   help="write a visibly watermarked report even when tasks are missing/retryable")
    args = p.parse_args(argv)

    if args.cmd == "plan":
        tasks = (make_stress_tasks(args.stress_seeds) if args.stress_grid else
                 make_tasks(args.presets_reps, args.random,
                            args.presets.split(",") if args.presets else None))
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        manifest = write_task_manifest(args.out, tasks,
                                       task_type="synthetic-stress" if args.stress_grid else "synthetic")
        print(f"{len(tasks)} tasks -> {args.out} sha256={manifest['tasks_sha256']}")
        return 0
    if args.cmd == "report":
        print(report(Path(args.out), args.tasks, partial=args.partial))
        return 0

    tasks, manifest_meta = load_task_manifest(args.tasks)
    if args.index is not None:
        tasks = tasks[args.index * args.chunk:(args.index + 1) * args.chunk]
    keep_awake()
    workroot = Path(args.workroot or tempfile.mkdtemp(prefix="gb_engine_"))
    if args.max_attempts < 1:
        p.error("--max-attempts must be >= 1")
    models = (DEFAULT_BENCH_MODELS if args.models == "default" else
              [m.strip() for m in args.models.split(",") if m.strip()] if args.models else [])
    ensure_run_config(args.out, manifest_meta, runner="engine_bench",
                      scientific_config={"models": models, "proxy_methods": EVAL_METHODS,
                                         "method_params": {m: BENCH_METHOD_PARAMS.get(m, {})
                                                           for m in models},
                                         "max_score_cells": MAX_SCORE_CELLS,
                                         "context_truth_scoring": "macro_by_context",
                                         "artifact_policy": ("all" if args.keep else
                                                             "sentinels" if args.keep_artifacts
                                                             else "none")})
    jobs = [(t, args.out, workroot,
             args.keep or (args.keep_artifacts and is_sentinel_task(t)),
             models, args.max_attempts, args.retry_permanent)
            for t in tasks]
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    retryable = False
    if args.jobs > 1:
        from concurrent.futures.process import BrokenProcessPool
        try:
            for tid, status, detail in pool_results(_run_and_save, jobs, args.jobs):
                retryable |= status in RETRYABLE_STATUSES
                print(f"{tid}: {status}: {detail}", flush=True)
        except BrokenProcessPool as exc:
            retryable = True
            print(f"worker died ({exc}); finished rows are saved, resubmit to resume", flush=True)
    else:
        for j in jobs:
            tid, status, detail = _run_and_save(j)
            retryable |= status in RETRYABLE_STATUSES
            print(f"{tid}: {status}: {detail}", flush=True)
    return 2 if retryable else 0


if __name__ == "__main__":
    sys.exit(main())
