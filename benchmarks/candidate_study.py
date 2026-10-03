"""Candidate study: every cleaning pipeline scored on real data under controlled degradations.

Exploratory instrument, not a verdict. For each real CELLxGENE dataset and each degradation
condition (clean, thinned depth, ambient RNA, doublets, damaged cells, mixed) the engine's
candidates, the fixed baseline and one-change variants of the baseline are all run and scored on
BOTH halves of the CollecTRI reference (selection half = what screening sees, evaluation half =
never seen), so selection rules and "where does cleaning help" can be analysed offline from one
table (``candidates.csv``).

    python benchmarks/candidate_study.py plan   --data-root D --out R/tasks.json [--conditions clean,mixed]
    python benchmarks/candidate_study.py run    --tasks R/tasks.json --refs D/references --out R --jobs 4
    python benchmarks/candidate_study.py report --out R

Every condition goes through a degraded COPY of the dataset (also ``clean``), so all share one
code path; the copy keeps the README/text files the context detector reads.
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
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import real_bench  # noqa: E402
from engine_bench import BASELINE, EVAL_METHODS, keep_awake, score_candidates  # noqa: E402

CONDITIONS = ("clean", "depth_30", "depth_10", "ambient_20", "doublets_20", "damaged_15", "mixed")
DOUBLET_FRAC, AMBIENT_FRAC, DAMAGED_FRAC, DAMAGED_KEEP = 0.2, 0.2, 0.15, 0.3
# Baseline with one change each: (decision, option).
VARIANTS = [("cell_qc_thresholds", "mad"), ("cell_qc_thresholds", "mad_per_batch"),
            ("doublets", "scrublet"), ("gene_filter", "min_cells_drop_qc_genes"),
            ("context_scope", "per_cell_type"), ("context_scope", "selected_type"),
            ("normalization", "pearson_residuals"), ("normalization", "scran"),
            ("batch", "center_per_batch"), ("batch", "combat"), ("batch", "harmony"),
            ("sparsity", "metacells"),
            ("hvg", "tf_plus_500"), ("hvg", "tf_plus_2000")]
ENGINE_SOURCES = ("engine", "baseline")
EPS = 1e-9


# -- degradations (integer UMI counts, CSR) ------------------------------------------------- #

def _csr(adata) -> sp.csr_matrix:
    X = sp.csr_matrix(adata.X)
    X.data = np.rint(X.data).astype(np.float32)
    return X


def _thin(X: sp.csr_matrix, keep: float, rng, rows=None, genes=None) -> sp.csr_matrix:
    """Binomial thinning of every count (or only rows x genes selected by the boolean masks)."""
    X = X.copy()
    sel = np.ones(X.nnz, dtype=bool)
    if rows is not None:
        sel &= np.repeat(rows, np.diff(X.indptr))
    if genes is not None:
        sel &= genes[X.indices]
    X.data[sel] = rng.binomial(X.data[sel].astype(np.int64), keep)
    X.eliminate_zeros()
    return X


def mito_mask(adata) -> np.ndarray:
    """Mitochondrial genes as the engine detects them (qc: MITO_PREFIXES on upper-cased names),
    from var_names and, when present, the symbol column (Ensembl var_names)."""
    from grnbench.ingest import MITO_PREFIXES, SYMBOL_COLUMNS

    labels = [pd.Index(adata.var_names).astype(str)]
    labels += [pd.Index(adata.var[c].astype(str)) for c in SYMBOL_COLUMNS if c in adata.var.columns]
    mask = np.zeros(adata.n_vars, dtype=bool)
    for names in labels:
        up = names.str.upper()
        for pref in MITO_PREFIXES:
            mask |= np.asarray(up.str.startswith(pref))
    return mask


def _ambient(X, frac, rng):
    profile = np.asarray(X.sum(axis=0)).ravel().astype(np.float64)
    profile /= profile.sum()
    lib = np.asarray(X.sum(axis=1)).ravel()
    # Poisson(frac * lib_i * profile_g) per gene == Poisson(frac * lib_i) total, then multinomial
    k = rng.poisson(frac * lib)
    cells = np.repeat(np.arange(X.shape[0]), k)
    genes = rng.choice(X.shape[1], size=int(k.sum()), p=profile)
    extra = sp.coo_matrix((np.ones(len(cells), dtype=np.float32), (cells, genes)), shape=X.shape).tocsr()
    return (X + extra).tocsr()


def _doublets(obs_all, X, frac, rng):
    n = X.shape[0]
    k = min(n - 1, int(round(frac * n)))
    first = rng.choice(n, k, replace=False)
    second = (first + rng.integers(1, n, size=k)) % n          # any other cell
    if "cell_type" in obs_all:                                  # prefer a different cell type
        ct = obs_all["cell_type"].astype(str).to_numpy()
        tries = rng.integers(0, n, size=(k, 10))
        for j in range(tries.shape[1]):
            cand = tries[:, j]
            bad = (ct[second] == ct[first]) & (cand != first) & (ct[cand] != ct[first])
            second = np.where(bad, cand, second)
    new_X = (X[first] + X[second]).tocsr()
    obs = obs_all.iloc[first].copy()
    obs.index = [f"{obs_all.index[i]}-dbl{j}" for j, i in enumerate(first)]
    return sp.vstack([X, new_X], format="csr"), pd.concat([obs_all, obs])


def degrade(adata, condition: str, seed: int = 0):
    """Return a degraded copy of ``adata`` (counts in X). ``clean`` is the identity."""
    import anndata as ad

    if condition not in CONDITIONS:
        raise ValueError(f"unknown condition {condition!r}; choose from {CONDITIONS}")
    rng = np.random.default_rng(seed)
    X, obs = _csr(adata), adata.obs.copy()
    steps = {"clean": [], "depth_30": ["depth_30"], "depth_10": ["depth_10"],
             "ambient_20": ["ambient"], "doublets_20": ["doublets"], "damaged_15": ["damaged"],
             "mixed": ["doublets", "ambient", "depth_30"]}[condition]
    for step in steps:
        if step.startswith("depth"):
            X = _thin(X, int(step[6:]) / 100, rng)
        elif step == "ambient":
            X = _ambient(X, AMBIENT_FRAC, rng)
        elif step == "doublets":
            X, obs = _doublets(obs, X, DOUBLET_FRAC, rng)
        elif step == "damaged":
            hit = np.zeros(X.shape[0], dtype=bool)
            hit[rng.choice(X.shape[0], int(round(DAMAGED_FRAC * X.shape[0])), replace=False)] = True
            X = _thin(X, DAMAGED_KEEP, rng, rows=hit, genes=~mito_mask(adata))
    return ad.AnnData(X, obs=obs, var=adata.var.copy())


def write_degraded_copy(path, condition: str, seed: int, dest: Path) -> Path:
    """Load ``path`` as the pipeline does, degrade, write ``dest/<name>.h5ad`` and copy the text
    files the context detector reads, so context detection is unchanged."""
    from grnbench.context import discover_text_files
    from grnbench.ingest import load_any

    dest.mkdir(parents=True, exist_ok=True)
    adata, _ = load_any(path)
    out = dest / (Path(path).stem + ".h5ad")
    degrade(adata, condition, seed).write_h5ad(out)
    for f in discover_text_files(path):
        shutil.copy2(f, dest / Path(f).name)
    return out


# -- tasks ------------------------------------------------------------------------------- #

def make_tasks(data_root: Path, conditions=CONDITIONS, max_datasets: int | None = None) -> list[dict]:
    datasets = [t for t in real_bench.make_tasks(Path(data_root)) if t["kind"] == "cellxgene"]
    datasets = datasets[:max_datasets] if max_datasets else datasets
    return [{"id": f"{d['id']}__{c}", "dataset": d["id"], "condition": c, "path": d["path"],
             "truth": d["truth"], "seed": 0} for d in datasets for c in conditions]


# -- candidates -------------------------------------------------------------------------- #

def _choices_key(choices: dict) -> str:
    return json.dumps(choices, sort_keys=True)


def variant_candidates():
    """``([{name, choices, source}], skipped)`` for BASELINE with one change each."""
    from grnbench.kb import DECISION_INDEX

    out, skipped = [], []
    for decision, option in VARIANTS:
        ok, why = DECISION_INDEX[decision].option(option).availability()
        if not ok:
            skipped.append({"variant": f"{decision}={option}", "reason": why})
            continue
        out.append({"name": f"base+{decision}={option}", "choices": {**BASELINE, decision: option},
                    "source": "variant"})
    return out, skipped


def _mean(x):
    x = [v for v in x if v is not None and np.isfinite(v)]
    return float(np.mean(x)) if x else float("nan")


def candidate_row(name, recs, table: pd.DataFrame) -> dict:
    """Screening-table means over the contexts of one candidate."""
    rows = table[table["candidate"] == name]
    ok = rows[rows["error"].isna()]
    errs = [r["error"] for r in recs if r.get("error")]
    out = {"n_contexts": int(len(ok)), "failed_contexts": len(errs),
           "error": errs[0][:300] if errs else None}
    for col in ("edge_stability", "depth_leakage", "batch_leakage"):
        out[col] = _mean(ok[col]) if col in ok else float("nan")
    for m in EVAL_METHODS:
        col = f"reference_auprc_ratio__{m}"
        out[f"sel_{m}"] = _mean(ok[col]) if col in ok else float("nan")
    out["sel_mean"] = _mean([out[f"sel_{m}"] for m in EVAL_METHODS])
    out["screen_seconds"] = round(float(ok["screen_seconds"].sum()), 1) if "screen_seconds" in ok else None
    return out


def candidate_genes(recs) -> set:
    """Union over a candidate's processed contexts of the genes it kept."""
    import anndata as ad

    out = set()
    for rec in recs:
        if rec.get("processed"):
            out |= set(map(str, ad.read_h5ad(rec["processed"], backed="r").var_names))
    return out


# -- one task ---------------------------------------------------------------------------- #

def run_study_task(task: dict, refs: Path | None, workroot: Path) -> dict:
    from grnbench.autoprep import diagnose, plan_preprocessing
    from grnbench.profile import build_profile
    from grnbench.screen import rank_candidates, screen_candidates
    from grnbench.selector import select

    t0 = time.time()
    row = {"task": task["id"], "dataset": task["dataset"], "condition": task["condition"],
           "seed": task["seed"], "error": None, "rules": "recommended"}
    work = Path(workroot) / task["id"]
    shutil.rmtree(work, ignore_errors=True)
    try:
        # task lists from real_bench.make_tasks hold portable "~/..." paths since #29
        src = Path(task["path"]).expanduser()
        path = write_degraded_copy(src, task["condition"], task["seed"], work / "data")
        diag = diagnose(path)
        row["n_cells_in"], row["n_genes_in"] = int(diag.adata.n_obs), int(diag.adata.n_vars)
        pairs = real_bench._reference(refs, task["truth"].get("organism"), diag.adata.var_names)
        tfs = sorted({s for s, _ in pairs}) or None
        row["reference_edges"] = len(pairs)
        if len(pairs) < real_bench.MIN_REFERENCE_EDGES:
            row["skipped"] = f"reference_edges {len(pairs)} < {real_bench.MIN_REFERENCE_EDGES}"
            row["seconds"] = round(time.time() - t0, 1)
            return row

        prof = build_profile(diag.adata, diag.qc.metrics, diag.context,
                             id_verdict=diag.id_verdict, tf_list=tfs)
        selection = select(prof)
        engine = selection.candidates(max_candidates=4, always=("context_scope",))
        specs = [{"name": c["name"], "choices": c["choices"], "source": "engine"} for c in engine]
        variants, row["skipped_variants"] = variant_candidates()
        specs += [{"name": "baseline", "choices": {**selection.chosen(), **BASELINE}, "source": "baseline"}]
        specs += [{**v, "choices": {**selection.chosen(), **v["choices"]}} for v in variants]

        seen, run_specs, info = {}, [], {}
        for s in specs:
            key = _choices_key(s["choices"])
            info[s["name"]] = {"name": s["name"], "source": s["source"], "choices": s["choices"],
                               "alias_of": seen.get(key)}
            if key in seen:
                continue
            if s["source"] == "variant":        # a variant that cannot even be planned is not fatal
                try:
                    plan_preprocessing(diag, tf_list=tfs, selection=selection, choices=s["choices"])
                except Exception as exc:  # noqa: BLE001
                    info[s["name"]]["error"] = f"{type(exc).__name__}: {exc}"[:300]
                    continue
            seen[key] = s["name"]
            if s["source"] != "engine":
                run_specs.append(s)

        res = screen_candidates(diag, path, work / "screen", tf_list=tfs, reference_pairs=pairs,
                                max_candidates=4, seed=task["seed"], selection=selection,
                                extra=[{"name": s["name"], "choices": s["choices"]} for s in run_specs],
                                proxy_methods=EVAL_METHODS)
        evaluation = res["reference"]["evaluation"]
        genes = {n: candidate_genes(r) for n, r in res["records"].items()}
        genes = {n: g for n, g in genes.items() if g}
        ref_cache = {}
        done = {}
        for name in sorted(res["records"], key=lambda n: n != "baseline"):     # baseline first
            t1 = time.time()
            sc = score_candidates({name: res["records"][name]}, {name: res["plans"][name]}, tfs, evaluation,
                                  pair_ref="baseline", pair_genes=genes, ref_cache=ref_cache)[name]
            done[name] = {**candidate_row(name, res["records"][name], res["rows"]), **sc,
                          "n_genes": len(genes[name]) if name in genes else None}
            done[name]["eval_seconds"] = round(time.time() - t1, 1)
            done[name]["seconds"] = round(done[name]["eval_seconds"] + (done[name]["screen_seconds"] or 0), 1)
        cands = []
        for s in specs:
            c = dict(info[s["name"]])
            target = c["alias_of"] or s["name"]
            if target in done:
                c.update(done[target])
            elif "error" not in c:
                c["error"] = "not run"
            c["choices"] = s["choices"]
            cands.append(c)
        row["candidates"] = cands

        subset = {c["alias_of"] or c["name"] for c in cands if c["source"] in ENGINE_SOURCES}
        sub = res["rows"][res["rows"]["candidate"].isin(subset)]
        ranking, _ = rank_candidates(sub)
        ok = len(ranking) and np.isfinite(ranking["score"].iloc[0])
        row["engine_auto"] = str(ranking.index[0]) if ok else None
        # The pre-D-021 reference-free rule, pinned so rows stay comparable across code versions:
        # rank_candidates now returns the baseline without a reference (D-021).
        row["engine_auto_noref"] = rule_choice(row, "engine", _noref_shipped)
    except Exception as exc:  # noqa: BLE001 - robustness is part of what is measured
        row["error"] = f"{type(exc).__name__}: {exc}"
        row["traceback"] = traceback.format_exc()[-4000:]
    row["seconds"] = round(time.time() - t0, 1)
    shutil.rmtree(work, ignore_errors=True)
    return row


def _run_and_save(args):
    task, refs, out, workroot = args
    path = Path(out) / "rows" / f"{task['id']}.json"
    if path.exists():
        return task["id"], "skipped"
    row = run_study_task(task, Path(refs) if refs else None, workroot)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(row, default=float), encoding="utf-8")
    os.replace(tmp, path)
    return task["id"], row.get("error") or row.get("skipped") or f"{row['seconds']}s"


# -- report ------------------------------------------------------------------------------ #

def _lookup(row: dict, name):
    return next((c for c in row.get("candidates", []) if c["name"] == name), None) if name else None


def _val(row: dict, name, field: str = "eval_mean") -> float:
    c = _lookup(row, name)
    v = c.get(field) if c else None
    return float(v) if v is not None and np.isfinite(v) else float("nan")


def task_arms(row: dict) -> dict:
    """eval_mean of each arm for one task (NaN when not available)."""
    cands = row.get("candidates", [])
    allv = [_val(row, c["name"]) for c in cands]
    engv = [_val(row, c["name"]) for c in cands if c["source"] in ENGINE_SOURCES]
    return {"baseline": _val(row, "baseline"), "engine_auto": _val(row, row.get("engine_auto")),
            "engine_auto_noref": _val(row, row.get("engine_auto_noref")),
            "rules": _val(row, row.get("rules")),
            "oracle_engine_baseline": float(np.nanmax(engv)) if np.isfinite(engv).any() else float("nan"),
            "oracle_all": float(np.nanmax(allv)) if np.isfinite(allv).any() else float("nan")}


# Selection rules compared offline: (label, candidate pool, score from screening fields only).
# The evaluation half never enters a score. "engine" = the engine's candidates + baseline, "all" =
# every candidate the study ran (variants included).
def _noref_shipped(c):
    return c["edge_stability"] - c["depth_leakage"] - max(np.nan_to_num(c["batch_leakage"]), 0.0)


# 4th item: "pair" = judged as pair_eval_mean - pair_base_eval_mean (genes shared with the baseline).
RULES = [
    ("with reference: selection half, engine pool", "engine", lambda c: c["sel_mean"]),
    ("with reference: selection half, all candidates", "all", lambda c: c["sel_mean"]),
    ("with reference: selection half, all candidates, judged on shared genes", "all",
     lambda c: c["sel_mean"], "pair"),
    ("no reference: shipped (stability - depth - batch), engine pool", "engine", _noref_shipped),
    ("no reference: stability - depth leakage, engine pool", "engine",
     lambda c: c["edge_stability"] - c["depth_leakage"]),
    ("no reference: stability - depth leakage, all candidates", "all",
     lambda c: c["edge_stability"] - c["depth_leakage"]),
    ("no reference: stability only, all candidates", "all", lambda c: c["edge_stability"]),
]


def rule_choice(row: dict, pool: str, score, field: str = "eval_mean") -> str | None:
    """Name of the candidate a rule picks: highest finite score among candidates that ran."""
    field = "pair_eval_mean" if field == "pair" else field
    best, best_s = None, -np.inf
    for c in row.get("candidates", []):
        if c.get("alias_of") or not np.isfinite(_val(row, c["name"], field)):
            continue
        if pool == "engine" and c["source"] not in ENGINE_SOURCES:
            continue
        try:
            s = float(score(c))
        except (KeyError, TypeError):
            continue
        if np.isfinite(s) and s > best_s:
            best, best_s = c["name"], s
    return best


def rules_section(arms_by_cond: dict, conds: list[str]) -> list[str]:
    failed = {cond: sum(1 for r in arms_by_cond[cond][0] if not np.isfinite(_val(r, "baseline")))
              for cond in conds}
    lines = ["", "## Selection rules: mean (chosen - baseline) eval, wins/losses vs baseline", "",
             "Tasks where the baseline itself failed (QC gate refusal or error) are left out of the "
             "mean; `+k` counts those where the rule still produced a scored network.", "",
             "Baseline failed: " + ", ".join(f"{c} {failed[c]}/{len(arms_by_cond[c][0])}" for c in conds), "",
             "| rule | " + " | ".join(conds) + " |", "|---|" + "---|" * len(conds)]
    for label, pool, score, *fld in RULES:
        field = fld[0] if fld else "eval_mean"
        cells = []
        for cond in conds:
            d, rescued = [], 0
            for r in arms_by_cond[cond][0]:
                pick = rule_choice(r, pool, score, field)
                base = _val(r, pick, "pair_base_eval_mean") if field == "pair" else _val(r, "baseline")
                if not np.isfinite(base):
                    rescued += int(pick is not None)
                    continue
                d.append(_val(r, pick, "pair_eval_mean" if field == "pair" else field) - base
                         if pick else np.nan)
            d = np.array([x for x in d if np.isfinite(x)])
            cell = (f"{d.mean():+.3f} ({int((d > EPS).sum())}/{int((d < -EPS).sum())})" if len(d) else "-")
            cells.append(cell + (f" +{rescued}" if rescued else ""))
        lines.append(f"| {label} | " + " | ".join(cells) + " |")
    return lines


def candidates_table(rows: list[dict]) -> pd.DataFrame:
    recs = []
    for r in rows:
        for c in r.get("candidates", []):
            recs.append({"task": r["task"], "dataset": r["dataset"], "condition": r["condition"],
                         "seed": r["seed"], "n_cells_in": r.get("n_cells_in"),
                         "n_genes_in": r.get("n_genes_in"), "reference_edges": r.get("reference_edges"),
                         "engine_auto": r.get("engine_auto"), "engine_auto_noref": r.get("engine_auto_noref"),
                         **{k: (json.dumps(v, sort_keys=True) if k == "choices" else v) for k, v in c.items()}})
    return pd.DataFrame(recs)


def variant_deltas(ok_rows: list[dict], name: str) -> tuple[list[float], list[float]]:
    """Per-task (variant - baseline) eval_mean on each arm's own universe, and on the genes shared."""
    own = [_val(r, name) - _val(r, "baseline") for r in ok_rows]
    shared = [_val(r, name, "pair_eval_mean") - _val(r, name, "pair_base_eval_mean") for r in ok_rows]
    return [x for x in own if np.isfinite(x)], [x for x in shared if np.isfinite(x)]


def summary(rows: list[dict]) -> str:
    conds = [c for c in CONDITIONS if any(r["condition"] == c for r in rows)]
    lines = ["# Candidate study - SUMMARY", "",
             "EXPLORATORY. Descriptive means only; no p-values, no claims. Unit = task "
             "(dataset x condition); eval = AUPRC ratio on the CollecTRI evaluation half, "
             "mean of the evaluation methods (1 = random).", "",
             "## Per condition", "",
             "| condition | ok | skipped | error | baseline | engine_auto | engine_auto_noref | rules | "
             "oracle engine+baseline | oracle all | noref worse than baseline |",
             "|---|---|---|---|---|---|---|---|---|---|---|"]
    arms_by_cond = {}
    for cond in conds:
        rs = [r for r in rows if r["condition"] == cond]
        ok = [r for r in rs if not r.get("error") and not r.get("skipped")]
        arms = pd.DataFrame([task_arms(r) for r in ok]) if ok else pd.DataFrame()
        arms_by_cond[cond] = (ok, arms)
        f = lambda k: f"{arms[k].mean():.3f}" if len(arms) and arms[k].notna().any() else "-"  # noqa: E731
        both = arms[["engine_auto_noref", "baseline"]].dropna() if len(arms) else pd.DataFrame()
        worse = f"{int((both['engine_auto_noref'] < both['baseline'] - EPS).sum())}/{len(both)}"
        lines.append(f"| {cond} | {len(ok)} | {sum(1 for r in rs if r.get('skipped'))} | "
                     f"{sum(1 for r in rs if r.get('error'))} | {f('baseline')} | {f('engine_auto')} | "
                     f"{f('engine_auto_noref')} | {f('rules')} | {f('oracle_engine_baseline')} | "
                     f"{f('oracle_all')} | {worse} |")
    names = sorted({c["name"] for r in rows for c in r.get("candidates", []) if c["source"] == "variant"})
    lines += ["", "## Variants: mean (eval_mean - baseline eval_mean), n tasks in brackets", ""]
    lines += ["| variant | " + " | ".join(conds) + " |", "|---|" + "---|" * len(conds)]
    shared_lines = ["", "## Variants vs baseline on shared genes: mean (pair_eval_mean - pair_base_eval_mean), "
                    "both on TF x the genes the variant and the baseline share, n tasks in brackets", "",
                    "| variant | " + " | ".join(conds) + " |", "|---|" + "---|" * len(conds)]
    for name in names:
        cells, shared = [], []
        for cond in conds:
            own, pair = variant_deltas(arms_by_cond[cond][0], name)
            cells.append(f"{np.mean(own):+.3f} [{len(own)}]" if own else "-")
            shared.append(f"{np.mean(pair):+.3f} [{len(pair)}]" if pair else "-")
        lines.append(f"| {name} | " + " | ".join(cells) + " |")
        shared_lines.append(f"| {name} | " + " | ".join(shared) + " |")
    lines += shared_lines
    lines += rules_section(arms_by_cond, conds)
    skipped = sorted({s["variant"] for r in rows for s in r.get("skipped_variants", [])})
    if skipped:
        lines += ["", "Variants skipped (option unavailable): " + ", ".join(skipped)]
    errs = [r for r in rows if r.get("error")]
    if errs:
        lines += ["", "## Task errors", ""] + [f"- {r['task']}: {str(r['error'])[:200]}" for r in errs]
    return "\n".join(lines) + "\n"


# Confirmatory family of EXP-20261003-LM-1, fixed before its run. Units = datasets; endpoint =
# eval_mean (AUPRC ratio on the CollecTRI evaluation half, mean of EVAL_METHODS); both arms are
# scored on the dataset's own task. "Gene-set preserving" = same hvg and gene_filter choices as
# the baseline, so both arms are judged on the same TF x gene universe.
def _preserves_gene_set(c) -> bool:
    ch = c.get("choices") or {}
    return all(ch.get(k) == BASELINE[k] for k in ("hvg", "gene_filter"))


CONFIRM_ARMS = {
    "baseline": lambda r: "baseline",
    "noref_shipped": lambda r: rule_choice(r, "engine", _noref_shipped),
    "ref_geneset_preserving": lambda r: _rule_pool(r, _preserves_gene_set),
}
CONFIRM_FAMILY = [("baseline", "noref_shipped"), ("ref_geneset_preserving", "baseline")]
CONFIRM_PRIMARY_CONDITION = "clean"


def _rule_pool(row: dict, keep) -> str | None:
    best, best_s = None, -np.inf
    for c in row.get("candidates", []):
        if c.get("alias_of") or not keep(c) or not np.isfinite(_val(row, c["name"])):
            continue
        s = c.get("sel_mean")
        if s is not None and np.isfinite(s) and s > best_s:
            best, best_s = c["name"], s
    return best


def confirm(rows: list[dict], condition: str) -> tuple[list[dict], str]:
    """Paired t over datasets for CONFIRM_FAMILY in one condition, Holm over both directions."""
    from grnbench.eval.comparison import compare_units, holm_adjust

    ok = [r for r in rows if r["condition"] == condition and not r.get("error") and not r.get("skipped")]
    recs, raw = [], []
    for a, b in CONFIRM_FAMILY:
        pairs = []
        for r in ok:
            va, vb = _val(r, CONFIRM_ARMS[a](r)), _val(r, CONFIRM_ARMS[b](r))
            if np.isfinite(va) and np.isfinite(vb):
                pairs.append((va, vb))
        rec = {"a": a, "b": b, "n": len(pairs), "excluded": len(ok) - len(pairs)}
        if len(pairs) >= 3:
            res = compare_units(pd.DataFrame(pairs, columns=[a, b]))["pairs"][0]
            rec.update(diff=float(res["difference"]), ci=res["simultaneous_ci"])
            # compare_units abstains (None) on a zero-variance sample: no evidence either way
            raw += [1.0 if res["p_a_gt_b"] is None else res["p_a_gt_b"],
                    1.0 if res["p_b_gt_a"] is None else res["p_b_gt_a"]]
        recs.append(rec)
    if raw and len(raw) == 2 * len(recs):
        adj = holm_adjust(raw)
        for i, r in enumerate(recs):
            r["p_a_gt_b"], r["p_b_gt_a"] = float(adj[2 * i]), float(adj[2 * i + 1])
            r["supported"] = r["a"] if r["p_a_gt_b"] < .05 else r["b"] if r["p_b_gt_a"] < .05 else "unresolved"
    lines = [f"### {condition}", "", "| A | B | n | excluded | mean A-B | 95% interval | p(A>B) Holm | "
             "p(B>A) Holm | supported |", "|---|---|---|---|---|---|---|---|---|"]
    for r in recs:
        ci = r.get("ci")
        ci_txt = (f"[{ci[0]:+.3f}, {ci[1]:+.3f}]"
                  if ci and all(v is not None and np.isfinite(v) for v in ci) else "unbounded")
        lines.append(f"| {r['a']} | {r['b']} | {r['n']} | {r['excluded']} | "
                     f"{r.get('diff', float('nan')):+.3f} | {ci_txt} | {r.get('p_a_gt_b', float('nan')):.3g} | "
                     f"{r.get('p_b_gt_a', float('nan')):.3g} | {r.get('supported', '-')} |")
    return recs, "\n".join(lines)


def confirm_report(out: Path) -> str:
    rows = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((Path(out) / "rows").glob("*.json"))]
    conds = [c for c in CONDITIONS if any(r["condition"] == c for r in rows)]
    text = ["# Confirmatory comparisons (CONFIRM_FAMILY)", "",
            f"Primary condition: `{CONFIRM_PRIMARY_CONDITION}`; other conditions are secondary. "
            "Paired t over datasets, one-sided per direction (each direction its own hypothesis), Holm over the directed hypotheses of "
            "each condition. Arms: " + ", ".join(CONFIRM_ARMS) + ".", ""]
    allrecs = {}
    for cond in [CONFIRM_PRIMARY_CONDITION] + [c for c in conds if c != CONFIRM_PRIMARY_CONDITION]:
        if cond not in conds:
            continue
        recs, table = confirm(rows, cond)
        allrecs[cond] = recs
        text += [table, ""]
    (Path(out) / "CONFIRM.json").write_text(json.dumps(allrecs, indent=1, default=float), encoding="utf-8")
    body = "\n".join(text)
    (Path(out) / "CONFIRM.md").write_text(body, encoding="utf-8")
    return body


def report(out: Path) -> str:
    rows = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((Path(out) / "rows").glob("*.json"))]
    candidates_table(rows).to_csv(Path(out) / "candidates.csv", index=False)
    text = summary(rows)
    (Path(out) / "SUMMARY.md").write_text(text, encoding="utf-8")
    return text


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("plan")
    a.add_argument("--data-root", required=True)
    a.add_argument("--out", required=True)
    a.add_argument("--conditions", default=",".join(CONDITIONS))
    a.add_argument("--max-datasets", type=int, default=None)
    b = sub.add_parser("run")
    b.add_argument("--tasks", required=True)
    b.add_argument("--refs", default=None)
    b.add_argument("--out", required=True)
    b.add_argument("--jobs", type=int, default=1)
    b.add_argument("--index", type=int, default=None)
    b.add_argument("--chunk", type=int, default=1)
    b.add_argument("--workroot", default=None)
    b.add_argument("--max-attempts", type=int, default=1,
                   help="accepted because hpc/engine_bench.sbatch passes it; rows are not retried here "
                        "(a failed task is a row with its error; a missing one reruns on resubmission)")
    c = sub.add_parser("report")
    c.add_argument("--out", required=True)
    d = sub.add_parser("confirm", help="CONFIRM_FAMILY tests (EXP-20261003-LM-1) -> CONFIRM.md")
    d.add_argument("--out", required=True)
    args = p.parse_args(argv)
    if args.cmd == "report":
        print(report(Path(args.out)))
        return 0
    if args.cmd == "confirm":
        print(confirm_report(Path(args.out)))
        return 0
    if args.cmd == "plan":
        conds = [x.strip() for x in args.conditions.split(",") if x.strip()]
        bad = [x for x in conds if x not in CONDITIONS]
        if bad:
            p.error(f"unknown conditions {bad}; choose from {list(CONDITIONS)}")
        tasks = make_tasks(Path(args.data_root), conds, args.max_datasets)
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(tasks, indent=0), encoding="utf-8")
        print(f"{len(tasks)} candidate-study tasks -> {args.out}")
        return 0
    tasks = json.loads(Path(args.tasks).read_text(encoding="utf-8"))
    if args.index is not None:
        tasks = tasks[args.index * args.chunk:(args.index + 1) * args.chunk]
    keep_awake()
    workroot = Path(args.workroot or tempfile.mkdtemp(prefix="gb_cs_"))
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    jobs = [(t, args.refs, args.out, workroot) for t in tasks]
    if args.jobs > 1:
        # ProcessPoolExecutor, not multiprocessing.Pool: when the OOM killer takes a worker, Pool
        # waits for it forever (EXP-20261002-LM-2 job 679448 idled ~20 h to its time limit);
        # the executor raises BrokenProcessPool, so the job ends and a resubmission resumes.
        from concurrent.futures import ProcessPoolExecutor, as_completed
        from concurrent.futures.process import BrokenProcessPool
        from multiprocessing import get_context
        with ProcessPoolExecutor(args.jobs, mp_context=get_context("spawn"),
                                 max_tasks_per_child=1) as pool:
            futures = [pool.submit(_run_and_save, j) for j in jobs]
            try:
                for fut in as_completed(futures):
                    print("%s: %s" % fut.result(), flush=True)
            except BrokenProcessPool as exc:
                print(f"worker died ({exc}); finished rows are saved, resubmit to resume", flush=True)
                return 3
    else:
        for j in jobs:
            print("%s: %s" % _run_and_save(j), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
