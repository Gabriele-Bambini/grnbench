"""Interface logic, kept free of Streamlit so it can be tested and scripted.

The app (``grnbench.ui.app``) only lays these functions out on screen:
analyse a dataset, describe the knowledge-base decisions, draw the QC
figures, run the models on every cleaned context, benchmark them against the
evaluation half of the reference, and package the run folder.
"""
from __future__ import annotations

import io
import json
import time
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd


def read_pairs(path) -> set[tuple[str, str]]:
    """Reference/gold edges from a table with source/target (or TF/target) columns."""
    from grnbench.eval.gold import load_gold
    return {(str(s), str(t)) for s, t in load_gold(str(path)).pairs}


def analyze(data_path, *, text_files=(), urls=(), notes=None, tf_list=None, choices=None):
    """Diagnosis + knowledge-base plan for one dataset: ``{"diag", "plan"}``."""
    from grnbench.autoprep import diagnose, plan_preprocessing
    diag = diagnose(data_path, text_files=text_files, urls=urls, notes=notes)
    plan = plan_preprocessing(diag, tf_list=tf_list, choices=choices)
    return {"diag": diag, "plan": plan}


def facts_table(ctx) -> pd.DataFrame:
    rows = [{"fact": f.key, "value": f.value, "source": f.source, "weight": f.weight,
             "evidence": f.snippet} for f in ctx.facts]
    return pd.DataFrame(rows, columns=["fact", "value", "source", "weight", "evidence"])


def gates_table(qc) -> pd.DataFrame:
    return pd.DataFrame([g.to_dict() for g in qc.gates])


def decision_rows(plan, decision_id: str) -> pd.DataFrame:
    """Every option of one decision: score, availability, evidence, reasons, references."""
    from grnbench.kb import DECISION_INDEX, REFERENCES
    res = plan.selection["results"][decision_id]
    d = DECISION_INDEX[decision_id]
    rows = []
    for so in res["options"]:
        opt = d.option(so["option"])
        refs = sorted({r for x in so["reasons"] for r in x["refs"]} | set(opt.refs))
        rows.append({
            "option": so["option"], "name": so["name"], "score": round(so["score"], 2),
            "chosen": so["option"] == plan.choices.get(decision_id),
            "available": so["available"], "why unavailable": so["unavailable_reason"],
            "evidence": opt.evidence,
            "reasons": " | ".join(f"{x['delta']:+g} {x['text']}" for x in so["reasons"]),
            "pros": "; ".join(opt.pros), "cons": "; ".join(opt.cons),
            "references": "; ".join(REFERENCES[r].citation.split(" - ")[0] for r in refs),
        })
    return pd.DataFrame(rows)


def qc_figures(diag, plan) -> dict:
    """Plotly figures of the per-cell QC distributions with the plan's thresholds."""
    import plotly.express as px

    obs = diag.adata.obs
    color = plan.value("batch_key") if plan.value("batch_key") in obs.columns else None
    df = pd.DataFrame({"detected genes": obs["_n_genes"].to_numpy(),
                       "log10 counts per cell": np.log10(obs["_total_counts"].to_numpy() + 1),
                       "% mitochondrial": obs["_pct_mito"].to_numpy()})
    if color:
        df[color] = obs[color].astype(str).to_numpy()
    figs = {}
    lines = {"detected genes": [("min_genes_per_cell", "min"), ("max_genes_per_cell", "max")],
             "log10 counts per cell": [("max_counts_per_cell", "max")],
             "% mitochondrial": [("max_pct_mito_cell", "max")]}
    for col, keys in lines.items():
        fig = px.histogram(df, x=col, color=color, nbins=80, barmode="overlay", opacity=0.6,
                           height=300)
        for key, label in keys:
            v = plan.value(key)
            if v is not None:
                x = np.log10(v + 1) if col.startswith("log10") else v
                fig.add_vline(x=x, line_dash="dash", line_color="crimson",
                              annotation_text=f"{label} {v:g}", annotation_position="top")
        fig.update_layout(margin=dict(l=10, r=10, t=30, b=10))
        figs[col] = fig
    ctype = plan.profile.get("cell_type_key")
    if ctype and color and ctype != color:
        tab = pd.crosstab(obs[ctype].astype(str), obs[color].astype(str))
        figs["composition"] = px.imshow(tab, text_auto=True, aspect="auto", height=300,
                                        labels=dict(x=color, y=ctype, color="cells"))
    elif ctype:
        counts = obs[ctype].astype(str).value_counts().rename_axis(ctype).reset_index(name="cells")
        figs["composition"] = px.bar(counts, x=ctype, y="cells", height=300)
    return figs


def retention_table(records) -> pd.DataFrame:
    """Cells/genes before and after every cleaning step, per context."""
    rows = []
    for rec in records:
        if not rec.get("manifest"):
            rows.append({"context": rec["context"], "step": "FAILED", "detail": rec.get("error")})
            continue
        for st in json.loads(Path(rec["manifest"]).read_text(encoding="utf-8"))["steps"]:
            rows.append({"context": rec["context"], "step": st["step"],
                         "cells": st.get("n_cells_after"), "genes": st.get("n_genes_after"),
                         "detail": st.get("detail", "")})
    return pd.DataFrame(rows)


def register_user_model(py_path, attr: str) -> list[dict]:
    """Register a user's model (``infer(adata, tf_list, **params)``) as a method."""
    from grnbench.plugins import register_plugins
    return register_plugins([f"{py_path}:{attr}"], replace=True)


def run_models(records, methods, out_dir, *, tf_list=None, strategy: str = "cpm_log1p",
               progress=None) -> pd.DataFrame:
    """Run every method on every cleaned context; write edge lists; one row per run."""
    import anndata as ad
    from grnbench.infer import run_method, save_edges

    out = Path(out_dir)
    rows = []
    todo = [(r, m) for r in records if r.get("processed") for m in methods]
    for i, (rec, method) in enumerate(todo):
        if progress:
            progress(i / max(len(todo), 1), f"{method} on {rec['context']}")
        adata = ad.read_h5ad(rec["processed"])
        ctx_dir = out / Path(rec["out_dir"]).name
        t0 = time.time()
        row = {"context": rec["context"], "method": method, "edges": None, "error": None,
               "processed": rec["processed"]}
        try:
            edges, _ = run_method(method, adata, tf_list=tf_list,
                                  assume_normalized=strategy != "raw_counts")
            row["edges"] = str(save_edges(edges, ctx_dir, method))
            row["n_edges"] = int(len(edges))
        except Exception as exc:  # noqa: BLE001 - a failing model is a result, not a crash
            row["error"] = f"{type(exc).__name__}: {exc}"
        row["seconds"] = round(time.time() - t0, 2)
        rows.append(row)
    if progress:
        progress(1.0, "done")
    return pd.DataFrame(rows)


def benchmark(runs: pd.DataFrame, reference_pairs, *, tf_list=None) -> pd.DataFrame:
    """Score every edge list on its context's full candidate universe."""
    import anndata as ad
    from grnbench.eval.metrics import GoldIndex, evaluate_ranking
    from grnbench.infer import load_edges
    from grnbench.tfs import resolve_tf_names

    gi = GoldIndex.coerce(set(reference_pairs))
    rows = []
    universes: dict[str, pd.DataFrame] = {}
    for run in runs.itertuples(index=False):
        if not run.edges:
            continue
        edges = load_edges(run.edges)
        key = str(run.processed)
        if key not in universes:
            genes = pd.Index(ad.read_h5ad(key, backed="r").var_names).astype(str)
            sources = resolve_tf_names(genes, tf_list) if tf_list else list(genes)
            universes[key] = pd.DataFrame([(s, g) for s in sources for g in genes if s != g],
                                          columns=["source", "target"])
        m = evaluate_ranking(edges, gi, with_curve=False, universe=universes[key])
        rows.append({"context": run.context, "method": run.method,
                     "AUPRC": m.get("auprc"), "AUPRC ratio": m.get("auprc_ratio"),
                     "EPR": m.get("epr"), "AUROC": m.get("auroc"),
                     "reference edges in universe": m.get("n_gold_in_universe")})
    return pd.DataFrame(rows)


def zip_folder(folder) -> bytes:
    buf = io.BytesIO()
    root = Path(folder)
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for f in sorted(root.rglob("*")):
            if f.is_file() and f.name.lower() != "desktop.ini":
                zf.write(f, f.relative_to(root))
    return buf.getvalue()


def new_workdir(base) -> Path:
    stamp = time.strftime("%Y%m%d-%H%M%S")
    path = Path(base) / f"run-{stamp}"
    path.mkdir(parents=True, exist_ok=True)
    return path


def summarize_profile(profile: dict) -> list[tuple[str, str]]:
    """Human-readable profile lines for the metric cards."""
    def fmt(v, spec=""):
        if v is None or (isinstance(v, float) and np.isnan(v)):
            return "-"
        return format(v, spec) if spec else str(v)
    return [
        ("Cells", fmt(profile["n_cells"], ",")), ("Genes", fmt(profile["n_genes"], ",")),
        ("Organism", fmt(profile["organism"])), ("Protocol", fmt(profile["platform"] or profile["technology"])),
        ("Zeros", fmt(profile["sparsity"], ".0%")), ("Median genes/cell", fmt(profile["median_genes_per_cell"], ",.0f")),
        ("Median counts/cell", fmt(profile["median_counts_per_cell"], ",.0f")),
        ("Median % mito", fmt(profile["median_pct_mito"], ".1f")),
        ("Batch column", fmt(profile["batch_key"])), ("Batches", fmt(profile["n_batches"])),
        ("Batch x cell type (Cramer's V)", fmt(profile["batch_celltype_cramers_v"], ".2f")),
        ("Batch silhouette", fmt(profile["batch_silhouette"], ".2f")),
        ("Cell-type column", fmt(profile["cell_type_key"])), ("Cell types", fmt(profile["n_cell_types"])),
        ("TFs detected", f"{profile['tf_detected']}/{profile['tf_requested']}" if profile["tf_requested"] else "-"),
    ]
