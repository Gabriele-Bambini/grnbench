"""Screen candidate preprocessing pipelines on the data itself.

Each candidate plan (the recommended one and one-change alternatives from the
selector) is executed and scored with cheap, comparable measurements:

* **retention** - cells, genes and TFs that survive;
* **depth leakage** - median |Spearman rho| between genes and log library size
  after processing (lower = normalization removed more depth signal);
* **batch leakage** - silhouette of batch labels on the processed PCA
  (lower = batches mixed; only with a batch column);
* **edge stability** - Jaccard overlap of the top-k edges inferred from two
  disjoint halves of the cells (higher = less noise-driven);
* **reference recovery** - AUPRC ratio of a fast proxy method against the
  *selection half* of a reference network.

Anti-circularity: the reference is split in two halves by edges within each
TF (:func:`split_reference`). Screening only ever sees the selection half; the
evaluation half is written next to the results and reserved for benchmarking
the models, so the pipeline is never chosen on the edges it is judged on.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

SCREEN_METHOD = "pearson_signed"   # seconds per dataset; a proxy, not a verdict
STABILITY_TOP_K = 500
STABILITY_HALF_CELLS = 150         # equal-size halves: stability must not reward cell count
MAX_SCREEN_CELLS = 3000


def split_reference(pairs, seed: int = 0) -> tuple[set, set]:
    """Split reference edges into (selection, evaluation) halves, balanced per TF."""
    rng = np.random.default_rng(seed)
    by_tf: dict[str, list] = {}
    for s, t in sorted(pairs):
        by_tf.setdefault(str(s), []).append((str(s), str(t)))
    select, evaluate = set(), set()
    for edges in by_tf.values():
        order = rng.permutation(len(edges))
        for rank, i in enumerate(order):
            (select if rank % 2 == 0 else evaluate).add(edges[i])
    return select, evaluate


def _dense(X) -> np.ndarray:
    return np.asarray(X.toarray() if hasattr(X, "toarray") else X, dtype=np.float32)


def _library(adata) -> np.ndarray:
    src = adata.layers["counts"] if "counts" in adata.layers else adata.X
    return np.asarray(src.sum(axis=1)).ravel()


def depth_leakage(adata, n_genes: int = 300) -> float:
    """Median |Spearman rho| between the most variable genes and log library size."""
    from scipy.stats import rankdata
    X = _dense(adata.X)
    lib = np.log1p(_library(adata))
    if X.shape[0] < 10 or np.std(lib) == 0:
        return float("nan")
    var = X.var(axis=0)
    cols = np.argsort(-var)[: min(n_genes, X.shape[1])]
    cols = cols[var[cols] > 0]
    r_lib = rankdata(lib)
    r_lib = (r_lib - r_lib.mean()) / r_lib.std()
    rhos = []
    for c in cols:
        r = rankdata(X[:, c])
        sd = r.std()
        if sd > 0:
            rhos.append(abs(float(np.mean((r - r.mean()) / sd * r_lib))))
    return float(np.median(rhos)) if rhos else float("nan")


def batch_leakage(adata, batch_key: str | None, seed: int = 0) -> float | None:
    from sklearn.decomposition import PCA
    from sklearn.metrics import silhouette_score

    if not batch_key or batch_key not in adata.obs.columns:
        return None
    labels = adata.obs[batch_key].astype(str).to_numpy()
    if len(set(labels)) < 2 or adata.n_obs < 20:
        return None
    X = _dense(adata.X)
    n = min(20, X.shape[0] - 1, X.shape[1] - 1)
    pcs = PCA(n_components=n, random_state=seed).fit_transform(X - X.mean(axis=0))
    return float(silhouette_score(pcs, labels, random_state=seed))


def _top_pairs(edges: pd.DataFrame, k: int) -> set:
    top = edges.reindex(edges["weight"].abs().sort_values(ascending=False).index).head(k)
    return set(zip(top["source"].astype(str).str.upper(), top["target"].astype(str).str.upper()))


def edge_stability(adata, tf_list, method: str = SCREEN_METHOD, k: int = STABILITY_TOP_K,
                   seed: int = 0) -> float:
    from grnbench.infer import run_method

    if adata.n_obs < 40:
        return float("nan")
    # Fixed-size halves: with all cells, a pooled pipeline (more cells) always
    # looks more stable than a per-cell-type one, whatever its quality.
    half = min(adata.n_obs // 2, STABILITY_HALF_CELLS)
    order = np.random.default_rng(seed).permutation(adata.n_obs)
    halves = [adata[np.sort(order[:half])].copy(), adata[np.sort(order[half:2 * half])].copy()]
    tops = []
    for h in halves:
        edges, _ = run_method(method, h, tf_list=tf_list, assume_normalized=True)
        tops.append(_top_pairs(edges, k))
    union = tops[0] | tops[1]
    return float(len(tops[0] & tops[1]) / len(union)) if union else float("nan")


def screen_processed(path, tf_list=None, reference_pairs=None, batch_key: str | None = None,
                     method: str = SCREEN_METHOD, seed: int = 0, proxy_methods=None,
                     reference_groups: dict[str, set] | None = None) -> dict:
    """All screening measurements for one processed.h5ad."""
    import anndata as ad
    from grnbench.eval.metrics import evaluate_ranking
    from grnbench.infer import run_method
    from grnbench.tfs import resolve_tf_names

    t0 = time.time()
    adata = ad.read_h5ad(path)
    X = adata.X.data if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    n_bad = int((~np.isfinite(X)).sum())
    if n_bad:
        raise ValueError(f"processed matrix has {n_bad} non-finite value(s) (NaN/inf); cannot screen")
    if adata.n_obs > MAX_SCREEN_CELLS:
        idx = np.sort(np.random.default_rng(seed).choice(adata.n_obs, MAX_SCREEN_CELLS, replace=False))
        adata = adata[idx].copy()
    tfs = resolve_tf_names(adata.var_names, tf_list) if tf_list else None
    row = {"n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars),
           "n_tfs": len(tfs) if tfs else None,
           "depth_leakage": depth_leakage(adata),
           "batch_leakage": batch_leakage(adata, batch_key, seed),
           "edge_stability": edge_stability(adata, tfs, method, seed=seed)}
    groups = reference_groups or ({"all": set(reference_pairs)} if reference_pairs else {})
    if groups:
        genes = pd.Index(adata.var_names).astype(str)
        sources = tfs if tfs else list(genes)
        universe = pd.DataFrame([(s, g) for s in sources for g in genes if s != g],
                                columns=["source", "target"])
        ratios = []
        for m_name in (proxy_methods or [method]):
            try:
                edges, _ = run_method(m_name, adata, tf_list=tfs, assume_normalized=True)
                per_context = []
                reachable = []
                for pairs in groups.values():
                    m = evaluate_ranking(edges, set(pairs), with_curve=False, universe=universe)
                    per_context.append(float(m.get("auprc_ratio", np.nan)))
                    reachable.append(int(m.get("n_gold_in_universe", 0)))
                finite_contexts = [v for v in per_context if np.isfinite(v)]
                value = float(np.mean(finite_contexts)) if finite_contexts else float("nan")
            except Exception:  # noqa: BLE001 - a proxy that cannot run scores NaN
                reachable, value = [], float("nan")
            row[f"reference_auprc_ratio__{m_name}"] = value
            ratios.append(value)
            if reachable:
                # This is a sum over context-specific evaluations: one shared edge can therefore
                # contribute once per context. It is not the size of a deduplicated pooled universe.
                row.setdefault("reference_context_edge_instances_in_universe", int(sum(reachable)))
        finite = [v for v in ratios if np.isfinite(v)]
        row["reference_auprc_ratio"] = float(np.mean(finite)) if finite else float("nan")
        row["reference_contexts"] = int(len(groups))
    row["screen_seconds"] = round(time.time() - t0, 2)
    return row


# (metric, higher_is_better, relative change that counts as a difference)
_COMPARE = [("n_cells", True, 0.05), ("n_tfs", True, 0.05), ("depth_leakage", False, 0.10),
            ("batch_leakage", False, 0.10), ("edge_stability", True, 0.05),
            ("reference_auprc_ratio", True, 0.05)]


def pros_cons(rows: pd.DataFrame, baseline: str = "recommended") -> dict[str, dict]:
    """Per candidate: which measurements are better/worse than the baseline."""
    out = {}
    if baseline not in set(rows["candidate"]):
        return out
    base = rows[rows["candidate"] == baseline].mean(numeric_only=True)
    for name, grp in rows.groupby("candidate", sort=False):
        cur = grp.mean(numeric_only=True)
        pros, cons = [], []
        for metric, higher, tol in _COMPARE:
            b, c = base.get(metric), cur.get(metric)
            if b is None or c is None or pd.isna(b) or pd.isna(c) or name == baseline:
                continue
            rel = (c - b) / abs(b) if b else (c - b)
            if abs(rel) < tol:
                continue
            (pros if (rel > 0) == higher else cons).append(f"{metric}: {b:.3g} -> {c:.3g}")
        out[name] = {"pros": pros, "cons": cons}
    return out


def rank_candidates(rows: pd.DataFrame, method: str | None = None) -> tuple[pd.DataFrame, str]:
    """Order candidates by a stated rule; returns ``(ranking, rule)``.

    With a reference: AUPRC ratio on the selection half - of ``method`` when
    given (the best cleaning depends on the model: per-cell-type contexts help
    correlation methods and starve mutual-information ones), else the mean over
    the proxy methods; edge stability breaks ties. Without a reference: edge
    stability (on equal-size cell halves) minus depth and batch leakage. Failed
    contexts rank a candidate last - a pipeline that cannot run is not better.
    """
    ok = rows[rows["error"].isna()] if "error" in rows else rows
    failed = set(rows.loc[rows["error"].notna(), "candidate"]) if "error" in rows else set()
    if ok.empty:
        return (pd.DataFrame({"score": -np.inf}, index=pd.Index(sorted(failed), name="candidate")),
                "every candidate pipeline failed on at least one context")
    agg = ok.groupby("candidate", sort=False).mean(numeric_only=True)
    for needed in ("edge_stability", "depth_leakage"):
        if needed not in agg:
            agg[needed] = np.nan
    col = f"reference_auprc_ratio__{method}" if method else "reference_auprc_ratio"
    if col in agg and agg[col].notna().any():
        agg["score"] = agg[col] + 1e-3 * agg["edge_stability"].fillna(0)
        rule = (f"AUPRC ratio of {method or 'the proxy methods (mean)'} on the selection half of "
                "the reference (edge stability breaks ties)")
    else:
        agg["score"] = (agg["edge_stability"].fillna(0) - agg["depth_leakage"].fillna(0)
                        - agg.get("batch_leakage", pd.Series(0, index=agg.index)).fillna(0).clip(lower=0))
        rule = "edge stability - depth leakage - batch leakage (no reference given)"
    agg.loc[agg.index.isin(failed), "score"] = -np.inf
    return agg.sort_values("score", ascending=False), rule


def screen_candidates(diag, input_path, out_dir, *, tf_list=None, reference_pairs=None,
                      reference_pairs_by_context: dict[str, set] | None = None,
                      selection=None, max_candidates: int = 4, seed: int = 0,
                      method: str = SCREEN_METHOD, extra: list[dict] | None = None,
                      proxy_methods: list[str] | None = None,
                      always_vary: tuple = ("context_scope",)) -> dict:
    """Execute and screen the recommended plan and its one-change alternatives.

    ``extra`` adds user-proposed candidates: ``[{"name": ..., "choices": {decision: option}}]``
    (choices not given keep the recommended option).

    Returns ``{"rows": DataFrame, "pros_cons": {...}, "plans": {name: PrepPlan},
    "records": {name: [context records]}, "reference": {...}}`` and writes
    ``screening.csv`` and ``reference_split.json`` into ``out_dir``.
    """
    from grnbench.autoprep import execute_plan, plan_preprocessing
    from grnbench.profile import build_profile
    from grnbench.selector import select

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    if selection is None:
        prof = build_profile(diag.adata, diag.qc.metrics, diag.context,
                             id_verdict=diag.id_verdict, tf_list=tf_list)
        selection = select(prof)
    if reference_pairs_by_context:
        context_union = set().union(*(set(v) for v in reference_pairs_by_context.values()))
        reference_pairs = set(reference_pairs or set()) | context_union
    select_half = eval_half = None
    select_by_context = eval_by_context = None
    if reference_pairs:
        select_half, eval_half = split_reference(reference_pairs, seed)
        if reference_pairs_by_context:
            select_by_context = {str(k): set(v) & select_half
                                 for k, v in reference_pairs_by_context.items()}
            eval_by_context = {str(k): set(v) & eval_half
                               for k, v in reference_pairs_by_context.items()}
        (out / "reference_split.json").write_text(json.dumps(
            {"selection_half": sorted(map(list, select_half)),
             "evaluation_half": sorted(map(list, eval_half)),
             "selection_by_context": ({k: sorted(map(list, v))
                                        for k, v in select_by_context.items()}
                                       if select_by_context else None),
             "evaluation_by_context": ({k: sorted(map(list, v))
                                         for k, v in eval_by_context.items()}
                                        if eval_by_context else None),
             "note": "screening used only the selection half; benchmark models on the evaluation half"},
            indent=1), encoding="utf-8")
    rows, plans, records = [], {}, {}
    candidates = selection.candidates(max_candidates=max_candidates, always=always_vary)
    for e in extra or []:
        candidates.append({"name": e["name"], "choices": {**selection.chosen(), **e["choices"]}})
    for cand in candidates:
        name = cand["name"]
        plan = plan_preprocessing(diag, tf_list=tf_list, selection=selection, choices=cand["choices"])
        plans[name] = plan
        from grnbench.autoprep import safe_folder_name
        safe = safe_folder_name(name)
        recs = execute_plan(plan, input_path, out / safe, tf_list=tf_list, seed=seed)
        records[name] = recs
        for rec, context in zip(recs, plan.contexts):
            rec["subset_key"] = context.get("subset_key")
            rec["subset_value"] = context.get("subset_value")
            row = {"candidate": name, "context": rec["context"], "error": rec["error"],
                   "strategy": plan.value("strategy")}
            if rec["processed"]:
                groups = None
                if select_by_context:
                    context_value = context.get("subset_value")
                    if context_value is not None and str(context_value) in select_by_context:
                        groups = {str(context_value): select_by_context[str(context_value)]}
                    else:
                        # A pooled prediction is judged against every context truth separately,
                        # then macro-averaged.  Scoring it against the union would reward edges
                        # that are correct in only one context as if they were universally true.
                        groups = select_by_context
                try:
                    row.update(screen_processed(rec["processed"], tf_list, select_half,
                                                plan.value("batch_key"), method, seed,
                                                proxy_methods=proxy_methods,
                                                reference_groups=groups))
                except Exception as exc:  # noqa: BLE001 - one bad context fails its candidate only
                    row["error"] = rec["error"] = f"screening failed: {type(exc).__name__}: {exc}"
            rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(out / "screening.csv", index=False)
    ranking, rule = rank_candidates(table) if len(table) else (pd.DataFrame(), "")
    best = (str(ranking.index[0]) if len(ranking) and np.isfinite(ranking["score"].iloc[0])
            else None)
    best_per_method = {}
    if reference_pairs and len(table):
        for m_name in (proxy_methods or [method]):
            r_m, _ = rank_candidates(table, m_name)
            if len(r_m) and np.isfinite(r_m["score"].iloc[0]):
                best_per_method[m_name] = str(r_m.index[0])
    (out / "screening_summary.json").write_text(json.dumps(
        {"best": best, "best_per_method": best_per_method, "rule": rule,
         "ranking": ranking.reset_index().to_dict("records"),
         "choices": {k: p.choices for k, p in plans.items()}}, indent=1, default=str), encoding="utf-8")
    return {"rows": table, "pros_cons": pros_cons(table[table["error"].isna()]) if len(table) else {},
            "ranking": ranking, "rule": rule, "best": best, "best_per_method": best_per_method,
            "plans": plans, "records": records,
            "reference": {"selection": select_half, "evaluation": eval_half,
                          "selection_by_context": select_by_context,
                          "evaluation_by_context": eval_by_context}}
