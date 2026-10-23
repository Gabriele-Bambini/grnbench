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


def analyze(data_path, *, text_files=(), urls=(), notes=None, tf_list=None, choices=None,
            orientation="auto"):
    """Diagnosis + knowledge-base plan for one dataset: ``{"diag", "plan"}``."""
    from grnbench.autoprep import diagnose, plan_preprocessing
    diag = diagnose(data_path, orientation=orientation, text_files=text_files, urls=urls, notes=notes)
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
        m = json.loads(Path(rec["manifest"]).read_text(encoding="utf-8"))
        qm = (m.get("qc") or {}).get("metrics") or {}
        n, g = qm.get("n_cells"), qm.get("n_genes")           # the size as loaded
        if n is not None and g is not None:
            rows.append({"context": rec["context"], "step": "loaded", "cells": n, "genes": g,
                         "detail": "as read from the file"})
        for st in m["steps"]:
            n, g = st.get("n_cells_after", n), st.get("n_genes_after", g)   # a step that keeps the matrix keeps the size
            rows.append({"context": rec["context"], "step": st["step"], "cells": n, "genes": g,
                         "detail": st.get("detail", "")})
    return pd.DataFrame(rows)


def runs_view(runs: pd.DataFrame) -> pd.DataFrame:
    """What to show of ``run_models``: no file paths, and the error column only when something failed."""
    view = runs.drop(columns=["processed", "edges"], errors="ignore")
    if "error" in view:
        view = view.drop(columns="error") if view["error"].isna().all() else view.fillna({"error": ""})
    return view


# --------------------------------------------------------------------------- #
# plain-language summaries for the screens
# --------------------------------------------------------------------------- #
GATE_LABELS = {
    "cells_below_min_genes": "Cells with very few detected genes", "mito_fraction": "Mitochondrial RNA share",
    "library_size_cv": "Evenness of sequencing depth", "sparsity": "Zeros in the count matrix",
    "small_dataset": "Number of cells", "small_gene_space": "Number of genes", "zero_count_cells": "Empty cells",
    "class_imbalance": "Balance of cell labels", "group_imbalance": "Balance of groups",
    "batch_annotation": "Batch information", "gene_id_type": "Gene identifiers",
    "duplicate_barcodes": "Cell barcodes", "doublet_heuristic": "Possible doublets",
}


def plural(n: int, word: str) -> str:
    """``1 model``, ``4 models``, ``0 models``."""
    return f"{n} {word}" + ("" if n == 1 else "s")


# What a warning on each check means for the person, in one sentence (the engine's own detail says what was measured).
GATE_MEANING = {
    "cells_below_min_genes": "Many cells have too few detected genes to be useful; the cleaning removes them.",
    "mito_fraction": "A high share of mitochondrial RNA usually marks damaged cells; the cleaning can remove the worst ones.",
    "library_size_cv": "Cells were sequenced to very different depths, which can create false correlations; the cleaning corrects for depth.",
    "sparsity": "Most entries are zero (dropout), which makes correlation-based models less reliable.",
    "small_dataset": "Few cells were measured, so every network estimate will be noisy.",
    "small_gene_space": "Fewer genes were measured than usual, so some regulators or targets may be missing from the data. Not blocking.",
    "zero_count_cells": "Cells with no counts at all are dropped before normalisation.",
    "class_imbalance": "Some groups have very few cells, so a network for them is unreliable.",
    "group_imbalance": "Some groups have very few cells, so a network for them is unreliable.",
    "batch_annotation": "Without batch information, differences between runs can look like regulation.",
    "gene_id_type": "Gene names that are not symbols will not match a TF list.",
    "duplicate_barcodes": "Some cell barcodes appeared twice (pooled libraries?); they were made unique.",
    "doublet_heuristic": "Some cells look like two cells captured together; the cleaning can remove them.",
}


def quality_checks(qc) -> pd.DataFrame:
    """The QC gates with names a biologist reads: ``check``, ``status`` (pass / warn / fail), ``detail``, ``meaning``."""
    rows = [{"check": GATE_LABELS.get(g.name, g.name.replace("_", " ").capitalize()), "status": g.status,
             "detail": g.detail, "meaning": GATE_MEANING.get(g.name, "")} for g in qc.gates]
    return pd.DataFrame(rows, columns=["check", "status", "detail", "meaning"])


def quality_verdict(qc) -> dict:
    """One sentence about the whole QC report: ``{"level": ok|warn|fail, "headline", "counts", "attention"}``."""
    checks = quality_checks(qc)
    counts = {k: int((checks["status"] == k).sum()) for k in ("pass", "warn", "fail")}
    n = len(checks)
    if counts["fail"]:
        level, headline = "fail", (f"{counts['fail']} of {n} quality checks failed: the data may not be usable as it is"
                                   + (f", and {counts['warn']} more need a look." if counts["warn"] else "."))
    elif counts["warn"]:
        level, headline = "warn", f"{counts['pass']} of {n} quality checks passed; {counts['warn']} to keep an eye on."
    else:
        level, headline = "ok", f"All {n} quality checks passed." if n else "No quality check was run."
    attention = checks[checks["status"] != "pass"].sort_values("status", key=lambda s: s.map({"fail": 0, "warn": 1}))
    return {"level": level, "headline": headline, "counts": counts, "attention": attention.to_dict("records")}


def plain_warning(text: str) -> str:
    """A plan warning about an option that cannot run, in words (``text`` unchanged when it is some other warning)."""
    import re
    m = re.match(r"^(\w+): best option '(.+?)' is unavailable \((.+)\); using '(.+?)'\.$", str(text))
    if not m:
        return str(text)
    from grnbench.kb import DECISIONS
    decision, best, why, used = m.groups()
    d = next((d for d in DECISIONS if d.id == decision), None)
    if d is None:
        return str(text)
    try:
        used = d.option(used).name
    except Exception:  # noqa: BLE001 - an id the knowledge base does not know is shown as it is
        pass
    return f"**{d.question}** The highest-scoring option, {best}, cannot run here ({why}); the engine uses \u201c{used}\u201d instead."


def clean_reason(text: str) -> str:
    """Drop the knowledge base's ``[decision = option]`` prefix from a reason shown to a person."""
    import re
    return re.sub(r"^\s*\[[^\]]*\]\s*", "", str(text or "")).strip()


def _plain_value(v) -> str:
    if v is None:
        return "off"
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    return str(v)


def decision_caption(dec, plan=None) -> str:
    """The line under a plan value: its reason and source, or after a user override what the engine chose and why,
    the engine's source when it is not the data, and the option the user chose when a choice set the value."""
    if getattr(dec, "auto_reason", None) is None:
        return f"{clean_reason(dec.reason)} (source: {dec.source})"
    why = clean_reason(dec.auto_reason) + (f"; source: {dec.auto_source}" if dec.auto_source != "data" else "")
    text = f"auto: {_plain_value(dec.auto_value)} ({why}) → you: {_plain_value(dec.value)}"
    chosen = plan.chosen_options(dec.key) if plan is not None else []
    return text + (f" (chosen: {', '.join(chosen)})" if chosen else "")


def pipeline_overview(plan) -> list[dict]:
    """One row per decision of the plan, in the knowledge base's order: stage, question, chosen option, doubt, impact."""
    from grnbench.kb import DECISIONS, STAGES
    rows = []
    for stage in STAGES:
        for d in DECISIONS:
            if d.stage != stage.id or d.id not in plan.selection["results"]:
                continue
            res = plan.selection["results"][d.id]
            rows.append({"id": d.id, "stage": stage.name, "question": d.question, "impact": d.impact,
                         "choice": d.option(plan.choices[d.id]).name, "uncertain": bool(res["uncertain"])})
    return rows


def reproduce_commands(plan, data_path, *, tf_list=None, seed: int) -> list[str]:
    """One ``grnbench prep`` command per context of ``plan``, with the folders ``execute_plan`` writes and ``seed``."""
    from grnbench.autoprep import unique_folder_names
    folders = unique_folder_names([str(c["name"]) for c in plan.contexts])
    cmds = []
    for ctx, name in zip(plan.contexts, folders):
        cli = plan.cli_args(data_path, "OUT/" + name, "tf_list.txt" if tf_list else None, context=ctx,
                            bounds_path="per_batch_bounds.json", seed=seed)
        cmds.append(shell_command(["grnbench", "prep", *cli]))
    return cmds


def shell_command(argv) -> str:
    """``argv`` as one line a POSIX shell (Git Bash too) reads back exactly: labels carry spaces, ``|``, ``(``, ``$``."""
    import shlex
    return shlex.join(str(a) for a in argv)


def clean_path(text) -> Path | None:
    """A path typed or pasted into a box: spaces and the quotes of "Copy as path" removed, ``~`` expanded."""
    t = str(text or "").strip().strip("\"'").strip()
    return Path(t).expanduser() if t else None


def demo_summary(paths: dict) -> str:
    """One line about a written demo dataset: size, regulators, known links."""
    import anndata as ad
    adata = ad.read_h5ad(paths["counts"], backed="r")
    n_cells, n_genes = adata.shape
    adata.file.close()
    tfs = len([line for line in Path(paths["tfs"]).read_text().splitlines() if line.strip()])
    links = max(len(Path(paths["gold"]).read_text().splitlines()) - 1, 0)
    return f"{n_cells:,} cells x {n_genes:,} genes, {tfs} regulators, {links} known links"


def run_summary(records, runs, bench=None) -> dict:
    """Headline numbers of a finished run: cells and genes kept, models run, the highest score (not a ranking)."""
    tbl = retention_table(records) if records else pd.DataFrame()
    ok = tbl[tbl["step"] != "FAILED"] if len(tbl) else tbl
    out = {"contexts": 0, "cells_loaded": None, "cells_used": None, "genes_used": None,
           "models_ok": 0, "models_total": 0, "top": None, "reference_edges": None}
    if len(ok):
        last = ok.groupby("context", sort=False).tail(1)
        out["contexts"] = int(len(last))
        out["cells_used"] = int(last["cells"].sum())
        out["genes_used"] = int(last["genes"].max())
        loaded = ok[ok["step"] == "loaded"]
        out["cells_loaded"] = int(loaded["cells"].max()) if len(loaded) else None
    if runs is not None and len(runs):
        out["models_total"] = int(len(runs))
        out["models_ok"] = int(runs["edges"].notna().sum()) if "edges" in runs else 0
    if bench is not None and len(bench) and "AUPRC ratio" in bench:
        scored = bench.dropna(subset=["AUPRC ratio"])
        if len(scored):
            row = scored.loc[scored["AUPRC ratio"].idxmax()]
            out["top"] = {"method": str(row["method"]), "context": str(row["context"]), "ratio": float(row["AUPRC ratio"])}
            if "reference edges in universe" in scored:
                out["reference_edges"] = int(scored["reference edges in universe"].max())
    return out


def register_user_model(py_path, attr: str) -> list[dict]:
    """Register a user's model (``infer(adata, tf_list, **params)``) as a method."""
    from grnbench.plugins import register_plugins
    return register_plugins([f"{py_path}:{attr}"], replace=True)


def _check_edges(edges: pd.DataFrame, genes) -> None:
    """Fail a model's run, with a message, when its edge list cannot be scored or drawn: a missing or
    non-finite score, a duplicate pair or self-loop, or a name that is not a gene of the data. Checked
    once here so the explorer pages can trust what they read."""
    from grnbench.eval.universe import validate_edges

    clean = validate_edges(edges)
    known = set(pd.Index(genes).astype(str).str.upper())
    unknown = sorted((set(clean["source"]) | set(clean["target"])) - known)
    if unknown:
        raise ValueError(f"{len(unknown)} edge endpoint(s) are not genes of the data, e.g. {unknown[:3]}")


RUN_COLUMNS = ["context", "method", "edges", "error", "processed", "n_edges", "seconds"]


def _run_models_here(records, methods, out_dir, *, tf_list=None, strategy: str = "cpm_log1p",
                     progress=None) -> pd.DataFrame:
    """Run every method on every cleaned context; write edge lists; one row per run."""
    import anndata as ad
    from grnbench.infer import run_method, save_edges

    out = Path(out_dir)
    rows = []
    todo = [(r, m) for r in records if r.get("processed") for m in methods]
    for i, (rec, method) in enumerate(todo):
        if progress:
            progress(i / max(len(todo), 1), f"{method} on {rec['context']} (job {i + 1} of {len(todo)})")
        adata = ad.read_h5ad(rec["processed"])
        ctx_dir = out / Path(rec["out_dir"]).name
        t0 = time.time()
        row = {"context": rec["context"], "method": method, "edges": None, "error": None,
               "processed": rec["processed"]}
        try:
            edges, prov = run_method(method, adata, tf_list=tf_list,
                                     assume_normalized=strategy != "raw_counts")
            _check_edges(edges, adata.var_names)
            row["edges"] = str(save_edges(edges, ctx_dir, method))
            row["n_edges"] = int(len(edges))
            row["seconds"] = round(time.time() - t0, 2)
            (ctx_dir / f"{method}.provenance.json").write_text(
                json.dumps({**prov, "context": rec["context"], "seconds": row["seconds"]}, indent=1, default=str),
                encoding="utf-8")
        except Exception as exc:  # noqa: BLE001 - a failing model is a result, not a crash
            row["error"] = f"{type(exc).__name__}: {exc}"
        row["seconds"] = round(time.time() - t0, 2)
        rows.append(row)
    if progress:
        progress(1.0, "done")
    return pd.DataFrame(rows, columns=RUN_COLUMNS)


def _models_child(conn, plugin_refs, args, kwargs) -> None:
    import multiprocessing
    import os
    import threading

    parent = multiprocessing.parent_process()

    def stop_if_orphaned():                 # the server was stopped or crashed: do not keep computing for nobody
        while parent is not None and parent.is_alive():
            time.sleep(1.0)
        os._exit(1)
    threading.Thread(target=stop_if_orphaned, daemon=True).start()
    try:
        from grnbench.plugins import register_plugins
        register_plugins(plugin_refs, replace=True)     # the user's own models, registered in the page's process
        conn.send(("done", _run_models_here(*args, progress=lambda f, m: conn.send(("progress", f, m)), **kwargs)))
    except BaseException as exc:  # noqa: BLE001 - handed to the page
        conn.send(("fail", f"{type(exc).__name__}: {exc}"))


def run_models(records, methods, out_dir, *, tf_list=None, strategy: str = "cpm_log1p",
               progress=None) -> pd.DataFrame:
    """``_run_models_here`` in a child process that this function kills if the page is stopped.

    Streamlit's Stop (and a rerun, and a closed tab) is only honoured at the next ``st.*`` call, which a model
    computing in the page's own thread never makes. So the page polls the child and calls ``progress`` (a
    ``st.*`` call) every half second: that is where Stop raises, and ``finally`` then ends the child.
    """
    import multiprocessing as mp
    from grnbench.plugins import _registered_refs

    ctx = mp.get_context("spawn")           # never fork a server that has threads (dask, tornado)
    recv, send = ctx.Pipe(duplex=False)
    child = ctx.Process(target=_models_child,       # not a daemon: a model may start its own workers (joblib, dask)
                        args=(send, list(dict.fromkeys(_registered_refs.values())), (records, methods, out_dir),
                              dict(tf_list=tf_list, strategy=strategy)))
    child.start()
    last = (0.0, "starting the model process")
    try:
        while True:
            if recv.poll(0.5):
                kind, *rest = recv.recv()
                if kind == "done":
                    return rest[0]
                if kind == "fail":
                    raise RuntimeError(rest[0])
                last = tuple(rest)
            elif not child.is_alive():
                raise RuntimeError(f"the model process ended without a result (exit code {child.exitcode})")
            if progress:
                progress(*last)
    finally:
        if child.is_alive():
            child.kill()
        child.join()


def screening_scope(screen) -> dict:
    """What the winner of a screening kept: ``{"best", "kept", "left_out"}`` (context names).

    The engine ranks candidates by the mean of their own contexts, so a candidate scored on fewer cell types
    than the recommended pipeline (``context_scope=selected_type`` keeps only the largest, the easiest) can win
    by leaving the hard ones out. The interface cannot change that ranking; it can say what was left out, so
    the cleaning is never narrowed silently. ``left_out`` is empty unless the winner's contexts are a strict
    subset of the recommended pipeline's.
    """
    plans = (screen or {}).get("plans") or {}
    best = (screen or {}).get("best")
    if best not in plans or "recommended" not in plans:
        return {"best": best, "kept": [], "left_out": []}
    names = lambda plan: [c["name"] for c in plan.contexts]      # noqa: E731
    kept, rec = names(plans[best]), names(plans["recommended"])
    left_out = [n for n in rec if n not in kept] if set(kept) < set(rec) else []
    return {"best": best, "kept": kept, "left_out": left_out}


def reference_for_scoring(screen, reference, uploaded=None, uploaded_name: str = "file"):
    """Which reference network a score may be computed on: ``(pairs or None, label, warning or None)``.

    One definition for every page that scores or flags edges (7, 9 and the quick run). After a screening
    that split the reference, the cleaning was chosen on the selection half, so only the evaluation half
    is fair; with no split the full reference is. An upload replaces both, with a warning when it shares
    edges with the selection half (scores on those edges are optimistic).
    """
    split = (screen or {}).get("reference") or {}
    selection, evaluation = split.get("selection"), split.get("evaluation")
    if uploaded is not None:
        if not uploaded:       # an empty file is still the person's choice: never swap it for the stored reference
            return None, "", f"{uploaded_name} has no edges (only a header row?), so nothing can be scored"
        pairs = set(uploaded)
        shared = len(_ref_upper(pairs) & _ref_upper(selection)) if selection else 0
        warning = (f"{shared} of the uploaded edges are in the half of the reference the engine used to choose "
                   "the cleaning pipeline, so scores on them are optimistic; upload edges it never saw, or "
                   "clear this upload to use the held-out half.") if shared else None
        return pairs, f"uploaded: {uploaded_name}", warning
    if selection and not evaluation:
        return None, "", ("the reference was too small to hold out an evaluation half after choosing the "
                          "cleaning pipeline on it, so no fair score is possible")
    if evaluation:
        return set(evaluation), "evaluation half (never used to choose the pipeline)", None
    if reference:
        return set(reference), "full reference (the pipeline was not chosen on it)", None
    return None, "", None


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
    out = pd.DataFrame(rows)
    if len(out) and not (out["reference edges in universe"] > 0).any():
        raise ValueError("no reference edge can be scored on this data: both its genes must be genes of the data "
                         "(species? gene ID type?) and its regulator must be in the TF list, if you gave one")
    return out


def zip_folder(folder) -> bytes:
    buf = io.BytesIO()
    root = Path(folder)
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for f in sorted(root.rglob("*")):
            if f.is_file() and f.name.lower() != "desktop.ini":
                zf.write(f, f.relative_to(root))
    return buf.getvalue()


def output_folder(base) -> Path:
    """The sidebar's Output folder, expanded, created and checked; ``ValueError`` says what to fix in one sentence."""
    text = str(base or "").strip()
    if not text:
        raise ValueError("The output folder is empty: type the folder where the demo data, uploads and runs go.")
    path = Path(text).expanduser()
    try:
        path.mkdir(parents=True, exist_ok=True)
        import tempfile
        tempfile.TemporaryFile(dir=path).close()      # a folder can exist and still be read-only
    except OSError as exc:
        raise ValueError(f"Cannot write to `{path}`: {exc.strerror or exc}. Choose another output folder.") from exc
    return path


def fresh_folder(path) -> Path:
    """``path`` if it is free or empty, else ``path-2``, ``path-3``...: a second run never meets the first one's files.

    The engine refuses to write into a folder that already holds results, and nothing is ever deleted for the person.
    """
    path = Path(path)
    candidate, n = path, 1
    while candidate.exists() and any(candidate.iterdir()):
        n += 1
        candidate = path.with_name(f"{path.name}-{n}")
    return candidate


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


# --------------------------------------------------------------------------- #
# looking at the data
# --------------------------------------------------------------------------- #
MAX_EMBED_CELLS = 3000
OKABE_ITO = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442", "#999999"]


def _colors(n: int) -> list[str]:
    import plotly.express as px
    return OKABE_ITO + list(px.colors.qualitative.Alphabet) if n > len(OKABE_ITO) else OKABE_ITO


def _normalised(adata):
    """The normalised AnnData: raw counts are CPM/1e4 + log1p'd on a copy; a sparse matrix stays sparse."""
    from grnbench.infer import prepare_adata
    return prepare_adata(adata)[0]


def _dense_columns(X, cols) -> np.ndarray:
    """float32 array of the columns ``cols`` of ``X``. Only such a bounded block is ever made dense: the whole
    matrix of a 247,914 x 8,749 Perturb-seq dataset would take 8.7 GB."""
    block = X[:, cols]
    return (block.toarray() if hasattr(block, "toarray") else np.asarray(block)).astype(np.float32)


def _gene_moments(X) -> tuple[np.ndarray, np.ndarray]:
    """Mean and variance of every gene (column), E[x^2] - E[x]^2; a sparse matrix is summed in float64, never made dense."""
    if hasattr(X, "toarray"):
        n = X.shape[0]
        mean = np.asarray(X.sum(axis=0, dtype=np.float64)).ravel() / n
        sq = np.asarray(X.multiply(X).sum(axis=0, dtype=np.float64)).ravel() / n
    else:
        X = np.asarray(X, dtype=np.float32)
        mean, sq = X.mean(axis=0), (X ** 2).mean(axis=0)
    return mean, sq - mean ** 2


def subsample_cells(adata, max_cells: int = MAX_EMBED_CELLS, seed: int = 0):
    if adata.n_obs <= max_cells:
        return adata
    rng = np.random.default_rng(seed)
    return adata[np.sort(rng.choice(adata.n_obs, max_cells, replace=False))].copy()


def embedding_frame(adata, method: str = "PCA", *, max_cells: int = MAX_EMBED_CELLS,
                    n_top: int = 2000, seed: int = 0) -> pd.DataFrame:
    """2-D coordinates (``x``, ``y``) of up to ``max_cells`` cells next to their annotations.

    PCA of the ``n_top`` most variable genes of the normalised data; ``UMAP`` is computed on the
    first PCs and needs umap-learn. ``out.attrs`` records the method and how many cells exist.
    """
    from sklearn.decomposition import PCA

    if method not in ("PCA", "UMAP"):
        raise ValueError(f"unknown embedding '{method}' (PCA or UMAP)")
    sub = subsample_cells(adata, max_cells, seed)
    X = _normalised(sub).X
    if X.shape[0] < 3 or X.shape[1] < 2:
        raise ValueError(f"too small to embed: {X.shape[0]} cells x {X.shape[1]} genes")
    top = np.argsort(-_gene_moments(X)[1], kind="stable")[: min(n_top, X.shape[1])]
    D = _dense_columns(X, top)          # the sampled cells x the most variable genes only
    D -= D.mean(axis=0)
    pcs = PCA(n_components=max(2, min(20, D.shape[0] - 1, D.shape[1])), random_state=seed).fit_transform(D)
    if method == "UMAP":
        try:
            import umap
        except ImportError as exc:
            raise RuntimeError("UMAP needs umap-learn (pip install umap-learn); PCA works without it") from exc
        xy = umap.UMAP(random_state=seed, n_neighbors=min(15, D.shape[0] - 1)).fit_transform(pcs)
    else:
        xy = pcs[:, :2]
    out = sub.obs.copy()
    out = out.loc[:, ~out.columns.duplicated()]
    out.insert(0, "y", xy[:, 1])
    out.insert(0, "x", xy[:, 0])
    out["detected genes"] = cell_detected_genes(sub)
    out.attrs.update(method=method, n_cells_total=int(adata.n_obs))
    return out


def color_choices(frame: pd.DataFrame) -> list[str]:
    """Annotations worth colouring an embedding by: few-level categories first, then numbers."""
    cats, nums = [], []
    for col in frame.columns:
        if col in ("x", "y"):
            continue
        s = frame[col]
        if pd.api.types.is_numeric_dtype(s) and not pd.api.types.is_bool_dtype(s):
            if s.notna().any():
                nums.append(col)
        elif 2 <= s.astype(str).nunique() <= 30:
            cats.append(col)
    return cats + nums


def embedding_figure(frame: pd.DataFrame, color: str | None = None, title: str = "", height: int = 430):
    import plotly.express as px

    d = frame.copy()
    method = frame.attrs.get("method", "PCA")
    kw = dict(x="x", y="y", height=height, opacity=0.75, hover_name=d.index.astype(str),
              labels={"x": f"{method} 1", "y": f"{method} 2"})
    if color and color in d:
        s = d[color]
        continuous = pd.api.types.is_numeric_dtype(s) and not pd.api.types.is_bool_dtype(s) and s.nunique() > 12
        if continuous:
            fig = px.scatter(d, color=color, color_continuous_scale="Viridis", **kw)
        else:
            d[color] = s.astype(str)
            fig = px.scatter(d, color=color, color_discrete_sequence=_colors(d[color].nunique()), **kw)
    else:
        fig = px.scatter(d, **kw)
    fig.update_traces(marker=dict(size=5))
    fig.update_layout(margin=dict(l=10, r=10, t=40, b=10), title=title)
    return fig


def _counts_matrix(adata):
    """The matrix whose zeros mean "no count": the counts layer when there is one, else X if non-negative."""
    mat = adata.layers["counts"] if "counts" in adata.layers else adata.X
    data = mat.data if hasattr(mat, "data") and hasattr(mat, "toarray") else np.asarray(mat)
    if data.size and np.nanmin(data) < 0:
        return None, "X has negative values (centred or scaled) and there is no counts layer"
    return mat, "counts layer" if "counts" in adata.layers else "X"


def gene_detection(adata) -> tuple[np.ndarray, str]:
    """Fraction of cells with a count above zero for every gene, and where it was read from.

    A centred or Pearson-residual matrix has no zeros that mean "not seen", so detection is read from
    the counts layer the pipeline keeps; without one and with negative values it is NaN, not a guess.
    """
    mat, source = _counts_matrix(adata)
    if mat is None:
        return np.full(adata.n_vars, np.nan), source
    return np.asarray((mat > 0).mean(axis=0)).ravel(), source


def cell_detected_genes(adata) -> np.ndarray:
    """Number of genes with a count above zero in every cell (NaN when it cannot be read, see above)."""
    mat, _ = _counts_matrix(adata)
    if mat is None:
        return np.full(adata.n_obs, np.nan)
    return np.asarray((mat > 0).sum(axis=1)).ravel().astype(float)


def gene_summary(adata, tfs=None) -> pd.DataFrame:
    """Per gene: detection rate, mean and variance of the normalised expression, TF or not.

    ``out.attrs["detection_source"]`` says where the detection rate was read from.
    """
    norm = _normalised(adata)
    mean, variance = _gene_moments(norm.X)
    detection, source = gene_detection(adata)
    out = pd.DataFrame({"gene": norm.var_names.astype(str), "detection": detection,
                        "mean": mean, "variance": variance})
    if tfs is not None:
        wanted = {str(t).upper() for t in tfs}
        out["TF"] = out["gene"].str.upper().isin(wanted)
    out = out.sort_values("detection", ascending=False, kind="stable", na_position="last").reset_index(drop=True)
    out.attrs["detection_source"] = source
    return out


def expression_frame(adata, genes, group_key: str | None = None, *, max_cells: int = 5000,
                     seed: int = 0) -> pd.DataFrame:
    """Long table (cell, gene, expression[, group]) of the normalised expression of ``genes``."""
    sub = subsample_cells(adata, max_cells, seed)
    names = pd.Index(sub.var_names.astype(str))
    found = [g for g in genes if g in names]
    if not found:
        return pd.DataFrame(columns=["cell", "gene", "expression"])
    X = _dense_columns(_normalised(sub).X, [names.get_loc(g) for g in found])     # only the chosen genes
    parts = []
    for j, g in enumerate(found):
        part = pd.DataFrame({"cell": sub.obs_names.astype(str), "gene": g, "expression": X[:, j]})
        if group_key and group_key in sub.obs:
            part[group_key] = sub.obs[group_key].astype(str).to_numpy()
        parts.append(part)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=["cell", "gene", "expression"])


def expression_figure(frame: pd.DataFrame, group_key: str | None = None, height: int = 420):
    import plotly.express as px

    if frame.empty:
        return px.scatter(height=height)
    if group_key and group_key in frame:
        n_genes = frame["gene"].nunique()
        fig = px.violin(frame, x=group_key, y="expression", color=group_key, box=True, points=False,
                        facet_col="gene", facet_col_wrap=min(n_genes, 3), height=height,
                        color_discrete_sequence=_colors(frame[group_key].nunique()))
    else:
        fig = px.violin(frame, x="gene", y="expression", color="gene", box=True, points=False, height=height,
                        color_discrete_sequence=_colors(frame["gene"].nunique()))
    fig.update_layout(margin=dict(l=10, r=10, t=30, b=10), showlegend=False)
    return fig


def load_stage(diag, records, stage: str):
    """The AnnData behind a "what to look at" choice: the data as loaded, or one cleaned context."""
    import anndata as ad
    if stage == "As loaded (raw counts)":
        return diag.adata
    for rec in records or []:
        if rec.get("processed") and stage == f"Cleaned: {rec['context']}":
            return ad.read_h5ad(rec["processed"])
    raise KeyError(stage)


def stage_choices(records) -> list[str]:
    return ["As loaded (raw counts)"] + [f"Cleaned: {r['context']}" for r in records or [] if r.get("processed")]


# --------------------------------------------------------------------------- #
# looking at the results
# --------------------------------------------------------------------------- #
def _ref_upper(reference) -> set[tuple[str, str]]:
    return {(str(s).upper(), str(t).upper()) for s, t in (reference or [])}


def edges_table(edges: pd.DataFrame, reference=None, *, tf: str | None = None, gene: str | None = None,
                top: int = 200) -> pd.DataFrame:
    """Ranked edges (best first) with their rank, optionally filtered, flagged when in the reference.

    Tied edges share the best rank of their block (competition ranking), so the table does not
    depend on the order of the rows in the file.
    """
    e = edges.sort_values("weight", ascending=False, kind="stable").reset_index(drop=True)
    e.insert(0, "rank", e["weight"].rank(method="min", ascending=False).astype(int))
    if tf:
        e = e[e["source"].astype(str).str.upper() == tf.strip().upper()]
    if gene:
        g = gene.strip().upper()
        e = e[e["target"].astype(str).str.upper().str.contains(g, regex=False)
              | e["source"].astype(str).str.upper().str.contains(g, regex=False)]
    if reference:
        ref = _ref_upper(reference)
        e = e.assign(**{"in reference": [(str(s).upper(), str(t).upper()) in ref
                                         for s, t in zip(e["source"], e["target"])]})
    return e.head(top).reset_index(drop=True)


def network_figure(edges: pd.DataFrame, reference=None, *, top_n: int = 60, tf: str | None = None,
                   seed: int = 0, height: int = 560):
    """The ``top_n`` best edges (of one TF if given) as a network: TFs squares, genes circles, an arrow from
    regulator to target, green edges are in the reference. A -> B and B -> A are two edges, drawn side by side."""
    import networkx as nx
    import plotly.graph_objects as go

    e = edges_table(edges, reference, tf=tf, top=top_n)
    fig = go.Figure()
    fig.update_layout(height=height, margin=dict(l=10, r=10, t=10, b=10), showlegend=True,
                      xaxis=dict(visible=False), yaxis=dict(visible=False), plot_bgcolor="rgba(0,0,0,0)")   # the page's own background
    if e.empty:
        fig.add_annotation(text=(f"No edges for TF '{tf}' in this model: type a TF name exactly (capitals do not matter)"
                                 if tf else "no edges to draw"), showarrow=False)
        return fig
    G = nx.DiGraph()        # source -> target: an undirected graph merged A -> B with B -> A and lost the direction
    ref = _ref_upper(reference)
    for r in e.itertuples(index=False):
        G.add_edge(str(r.source), str(r.target), weight=float(r.weight),
                   hit=(str(r.source).upper(), str(r.target).upper()) in ref)
    pos = nx.spring_layout(G.to_undirected(as_view=True), seed=seed, k=1.6 / max(len(G) ** 0.5, 1))
    sources = set(e["source"].astype(str))
    size = {n: min(14 + 3 * G.degree(n), 30) if n in sources else 9 for n in G}     # a hub must not swallow its neighbours

    def ends(u, v):
        """Where u -> v is drawn: when v -> u is an edge too, each of the two is shifted to its own right."""
        (x0, y0), (x1, y1) = pos[u], pos[v]
        if not G.has_edge(v, u):
            return x0, y0, x1, y1
        d = max(float(np.hypot(x1 - x0, y1 - y0)), 1e-9)
        dx, dy = 0.025 * (y1 - y0) / d, -0.025 * (x1 - x0) / d
        return x0 + dx, y0 + dy, x1 + dx, y1 + dy

    for hit, color, label in ((True, "#009E73", "in reference"), (False, "#B9B9B9", "not in reference")):
        xs, ys, width = [], [], 2 if hit else 1
        # Plotly drops an arrow, line and head, shorter than standoff + startstandoff + head backoff (1.3 x arrowsize x
        # arrowwidth for arrowhead 2) after cutting its tail to the empty label box (<= 2.83 px). The arrow therefore
        # runs the whole edge, ending at the rim of the target and starting `slack` px inside the rim of the source:
        # its head is drawn whenever the two markers do not overlap, at any plot width. An arrow drawn from the
        # middle of the edge lost its head between two close TF squares.
        slack = 1.3 * 1.4 * width + 3
        for u, v, d in G.edges(data=True):
            if d["hit"] == hit:
                x0, y0, x1, y1 = ends(u, v)
                xs += [x0, x1, None]
                ys += [y0, y1, None]
                fig.add_annotation(x=x1, y=y1, ax=x0, ay=y0, xref="x", yref="y", axref="x", ayref="y", text="",
                                   showarrow=True, arrowhead=2, arrowsize=1.4, arrowwidth=width, arrowcolor=color,
                                   standoff=size[v] / 2, startstandoff=max(0.0, size[u] / 2 - slack))
        if xs:
            fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", line=dict(color=color, width=width),
                                     name=label if reference else "edge", hoverinfo="skip"))
    # print a name only where it does not land on another name or node (pixel estimate; hover shows every name)
    px = {n: ((pos[n][0] + 1) * 450, (1 - pos[n][1]) * height / 2) for n in G}
    boxes = [(x - size[n] / 2, y - size[n] / 2, x + size[n] / 2, y + size[n] / 2) for n, (x, y) in px.items()]
    shown = set()
    for n in sorted(G, key=lambda n: (-G.degree(n), n)):
        (x, y), w = px[n], 6.2 * len(n) + 4
        box = (x - w / 2, y - size[n] / 2 - 16, x + w / 2, y - size[n] / 2 - 2)
        if not any(box[0] < b[2] and b[0] < box[2] and box[1] < b[3] and b[1] < box[3] for b in boxes):
            shown.add(n)
            boxes.append(box)
    for is_tf, symbol, color, label in ((True, "square", "#0072B2", "TF"), (False, "circle", "#E69F00", "gene")):
        nodes = [n for n in G if (n in sources) == is_tf]
        if nodes:
            fig.add_trace(go.Scatter(x=[pos[n][0] for n in nodes], y=[pos[n][1] for n in nodes],
                                     mode="markers+text",
                                     text=[n if n in shown else "" for n in nodes], textposition="top center",
                                     textfont=dict(size=10),
                                     marker=dict(symbol=symbol, size=[size[n] for n in nodes],
                                                 color=color, line=dict(width=1, color="white")),
                                     name=label, hovertext=[f"{n}: {G.degree(n)} edge(s)" for n in nodes], hoverinfo="text"))
    return fig


def add_random_line(fig, y: float = 1.0, *, horizontal_bars: bool = False):
    """Dotted line at ``y`` (a random ranking) named in the legend: a label on the plot sits on the last bar."""
    import plotly.graph_objects as go
    (fig.add_vline if horizontal_bars else fig.add_hline)(**{"x" if horizontal_bars else "y": y},
                                                          line_dash="dot", line_color="gray")
    fig.add_trace(go.Scatter(x=[None], y=[None], mode="lines", line=dict(color="gray", dash="dot"),
                             name=f"random ({y:g})"))


RANDOM_LEVEL = {"AUPRC ratio": 1.0, "EPR": 1.0, "AUROC": 0.5}


def score_figure(bench: pd.DataFrame, metric: str = "AUPRC ratio", height: int | None = None):
    """Models ranked by one score, best on top, with the random level marked (a bar chart, horizontal)."""
    import plotly.express as px
    from grnbench.ui import theme
    df = bench.dropna(subset=[metric]).sort_values(metric, ascending=True, kind="stable")
    several = df["context"].nunique() > 1
    pattern = "kind" if "kind" in df and df["kind"].nunique() > 1 else None
    fig = px.bar(df, x=metric, y="method", color="context" if several else None, pattern_shape=pattern,
                 orientation="h", barmode="group", text_auto=".2f" if metric == "AUROC" else ".1f",
                 height=height or max(220, 46 * len(df) + 90), labels={"method": ""},
                 color_discrete_sequence=None if several else [theme.ACCENT])
    fig.update_traces(textposition="outside", cliponaxis=False)
    # the legend sits under the bars: at the top right it covered the longest bar's value on a narrow screen
    fig.update_layout(margin=dict(l=10, r=30, t=10, b=10), legend=dict(orientation="h", yanchor="top", y=-0.18, x=0),
                      yaxis=dict(categoryorder="array", categoryarray=list(df["method"].drop_duplicates())))
    if metric in RANDOM_LEVEL:
        add_random_line(fig, RANDOM_LEVEL[metric], horizontal_bars=True)
    return fig


def _top_set(edges: pd.DataFrame, k: int) -> set:
    """Top-``k`` pairs with the repository's tie rule (``eval.universe.top_pairs``: a block of equal
    scores at the cut is kept whole, names folded to upper case). Only the head is validated."""
    from grnbench.eval.universe import top_pairs

    w = edges["weight"].to_numpy(dtype=float)
    if not np.isfinite(w).all():
        raise ValueError("non-finite scores")
    n = len(w)
    kk = min(int(k), n)
    thr = np.partition(w, n - kk)[n - kk]
    head = edges.loc[w >= thr, ["source", "target", "weight"]]
    return top_pairs(head, len(head))


def overlap_matrix(edges_by_method: dict, k: int = 200) -> pd.DataFrame:
    """Jaccard overlap of the top-``k`` edges of every pair of methods (tie blocks kept whole).

    ``attrs["sizes"]`` holds the realised set sizes (larger than ``k`` when scores tie at the cut) and
    ``attrs["skipped"]`` the methods that could not be compared, with the reason.
    """
    tops, skipped = {}, {}
    for m, e in edges_by_method.items():
        if e is None or not len(e):
            continue
        try:
            tops[m] = _top_set(e, k)
        except ValueError as exc:
            skipped[m] = str(exc)
    names = list(tops)
    mat = pd.DataFrame(np.eye(len(names)), index=names, columns=names)
    for a in names:
        for b in names:
            if a != b:
                union = len(tops[a] | tops[b])
                mat.loc[a, b] = len(tops[a] & tops[b]) / union if union else 0.0
    mat.attrs.update(sizes={m: len(t) for m, t in tops.items()}, skipped=skipped)
    return mat


def overlap_figure(mat: pd.DataFrame, k: int):
    import plotly.express as px
    fig = px.imshow(mat, text_auto=".2f", zmin=0, zmax=1, color_continuous_scale="Blues", aspect="auto",
                    height=120 + 48 * max(len(mat), 3), labels=dict(color="Jaccard"))
    fig.update_layout(margin=dict(l=10, r=10, t=30, b=10), title=f"Share of the top-{k} edges two models agree on")
    return fig


def _universe(genes, tf_list) -> dict:
    """The candidate pairs, every regulator x every other gene, described instead of listed (8,749 genes without a
    TF list are 76 million pairs): ``sources`` and ``targets`` folded to upper case as the scoring folds them, ``n``."""
    from grnbench.tfs import resolve_tf_names
    genes = pd.Index(genes).astype(str)
    if genes.str.upper().has_duplicates:
        raise ValueError("universe: duplicate gene names, or names that differ only in case")
    sources = set(pd.Index(resolve_tf_names(genes, tf_list) if tf_list else genes).str.upper())
    targets = set(genes.str.upper())
    return {"sources": sources, "targets": targets, "n": len(sources) * len(targets) - len(sources & targets)}


def _pr_curve_on_universe(edges: pd.DataFrame, gi, universe: dict, max_points: int = 200) -> dict:
    """``evaluate_ranking(edges, gi, universe=<every pair listed>)["pr_curve"]``, without listing the pairs.

    On an explicit universe the pairs a model did not score share one tied block at the bottom
    (``eval.universe.align_ranking``), and ``metrics.pr_curve`` reads the curve only at the end of each tied
    block, so that block is a single point: every reference pair of the universe at rank ``n``.
    """
    from grnbench.eval.metrics import _tie_blocks
    from grnbench.eval.universe import validate_edges

    src, tgt, n = universe["sources"], universe["targets"], universe["n"]
    if n == 0:
        raise ValueError("empty candidate universe")
    clean = validate_edges(edges, allow_empty=True)
    if not (clean["source"].isin(src) & clean["target"].isin(tgt)).all():
        raise ValueError("prediction contains pairs outside the explicit universe")
    clean = clean.sort_values("weight", ascending=False, kind="stable")
    n_pos = sum(1 for s, t in gi.lookup if s != t and s in src and t in tgt)
    if n_pos == 0:
        return {"recall": [], "precision": [], "baseline": 0.0, "n_pos": 0}
    ends = _tie_blocks(clean["weight"].to_numpy(dtype=float))
    hits, ranks = np.cumsum(gi.labels(clean))[ends - 1].astype(float), ends.astype(float)
    if n > len(clean):                  # the pairs the model did not score
        hits, ranks = np.r_[hits, n_pos], np.r_[ranks, n]
    precision, recall = hits / ranks, hits / n_pos
    idx = (np.unique(np.linspace(0, len(precision) - 1, max_points).astype(int)) if len(precision) > max_points
           else np.arange(len(precision)))
    return {"recall": [round(float(v), 6) for v in recall[idx]],
            "precision": [round(float(v), 6) for v in precision[idx]],
            "baseline": round(float(n_pos / n), 6), "n_pos": int(n_pos)}


def pr_figure(edges_by_method: dict, reference, genes, tf_list=None, height: int = 430):
    """Precision-recall curve of every method on the full TF x gene universe, with the random line."""
    import plotly.graph_objects as go
    from grnbench.eval.metrics import GoldIndex

    gi = GoldIndex.coerce(set(reference))
    universe = _universe(genes, tf_list)
    fig = go.Figure()
    baseline = None
    for m, edges in edges_by_method.items():
        if edges is None or not len(edges):
            continue
        curve = _pr_curve_on_universe(edges, gi, universe)
        if curve.get("recall"):
            fig.add_trace(go.Scatter(x=curve["recall"], y=curve["precision"], mode="lines", name=m,
                                     hovertemplate=f"{m}<br>recall %{{x:.3f}}<br>precision %{{y:.3f}}"))
            baseline = curve.get("baseline", baseline)
    if baseline:
        fig.add_hline(y=baseline, line_dash="dot", line_color="grey", annotation_text="random", annotation_position="top left")
    fig.update_layout(height=height, xaxis_title="recall", yaxis_title="precision", yaxis_range=[0, 1.02],
                      margin=dict(l=10, r=10, t=30, b=10), title="Which model ranks reference edges highest?")
    return fig


def _expected_hits(weights: np.ndarray, is_ref: np.ndarray, k: int) -> float:
    """Reference hits in the top ``k`` when ties at the cut are broken at random (expected value)."""
    n = len(weights)
    if n == 0:
        return float("nan")
    order = np.argsort(-weights, kind="stable")
    w, y = weights[order], is_ref[order].astype(float)
    kk = min(int(k), n)
    thr = w[kk - 1]
    above, block = w > thr, w == thr
    need = kk - int(above.sum())
    return float(y[above].sum() + need * y[block].mean())


def _fold(names) -> pd.Index:
    return pd.Index(names).astype(str).str.upper()


def _reference_in_data(reference, genes) -> pd.DataFrame:
    ref = pd.DataFrame(list(reference), columns=["source", "target"]).astype(str)
    ref = ref.assign(source=_fold(ref["source"]), target=_fold(ref["target"]))
    ref = ref[ref["source"] != ref["target"]].drop_duplicates()
    if genes is not None:
        present = set(_fold(genes))
        ref = ref[ref["source"].isin(present) & ref["target"].isin(present)]
    return ref


def tf_recovery_table(edges_by_method: dict, reference, genes=None, *, k_per_tf: int = 20,
                      max_tfs: int = 30) -> pd.DataFrame:
    """Fraction of each TF's reference targets that a method puts in the TF's top-``k_per_tf`` targets.

    With ``genes`` only reference TFs and targets that are in the data count, in the denominator too;
    a TF without any stored edge for a method gets NaN, not 0. Ties at the cut give the expected number
    of hits under random tie-breaking. Rows are the TFs ranked by at least one method, largest in-data
    regulon first; ``reference targets`` is that regulon size. A random ranking recovers about
    ``k / (genes - 1)``.
    """
    ref = _reference_in_data(reference, genes)
    targets = ref.groupby("source")["target"].apply(set)
    folded = {}
    for m, edges in edges_by_method.items():
        e = edges.assign(source=_fold(edges["source"]), target=_fold(edges["target"]))
        folded[m] = {tf: g for tf, g in e.groupby("source")}
    ranked = {tf for per in folded.values() for tf in per}
    tfs = [tf for tf in targets.apply(len).sort_values(ascending=False, kind="stable").index if tf in ranked][:max_tfs]
    out = pd.DataFrame(index=tfs, columns=list(edges_by_method), dtype=float)
    for m, per in folded.items():
        for tf in tfs:
            e = per.get(tf)
            if e is None or e.empty:
                continue
            hits = _expected_hits(e["weight"].to_numpy(dtype=float), e["target"].isin(targets[tf]).to_numpy(), k_per_tf)
            out.loc[tf, m] = hits / max(len(targets[tf]), 1)
    out.insert(0, "reference targets", [len(targets[t]) for t in tfs])
    return out


def tf_recovery_figure(tbl: pd.DataFrame, k_per_tf: int):
    import plotly.express as px
    heat = tbl.drop(columns="reference targets").astype(float)
    heat.index = [f"{tf} (n={n})" for tf, n in zip(tbl.index, tbl["reference targets"])]
    fig = px.imshow(heat, aspect="auto", zmin=0, zmax=1, color_continuous_scale="Viridis",
                    height=180 + 22 * len(heat), labels=dict(color="recovered at top-k"))
    fig.update_layout(margin=dict(l=10, r=10, t=40, b=10),
                      title=f"Fraction of each TF's reference targets (in the data) in its top-{k_per_tf} predictions")
    return fig


def tf_detail(edges: pd.DataFrame, reference, tf: str, k: int = 20, genes=None) -> tuple[pd.DataFrame, list[str], int]:
    """One TF: its top-``k`` predictions flagged in/not in the reference, the in-data reference targets it
    missed, and how many of its reference targets are not in the data at all (when ``genes`` is given)."""
    top = edges_table(edges, reference, tf=tf, top=k)
    top = top.rename(columns={"rank": "overall rank"})
    top.insert(0, "rank for this TF", np.arange(1, len(top) + 1))
    ref = _reference_in_data(reference, None)
    mine = set(ref.loc[ref["source"] == tf.strip().upper(), "target"])
    present = mine & set(_fold(genes)) if genes is not None else mine
    seen = set(top["target"].astype(str).str.upper())
    return top, sorted(t for t in present if t not in seen), len(mine) - len(present)


# --------------------------------------------------------------------------- #
# one-click run
# --------------------------------------------------------------------------- #
def baseline_plan(diag, tf_list=None):
    """The fixed baseline pipeline (D-021): what the engine returns when it has no reference."""
    from grnbench.autoprep import plan_preprocessing
    from grnbench.profile import Profile
    from grnbench.screen import BASELINE_CHOICES
    from grnbench.selector import select

    plan = plan_preprocessing(diag, tf_list=tf_list)
    selection = select(Profile(**plan.profile))
    return plan_preprocessing(diag, tf_list=tf_list, selection=selection,
                              choices={**selection.chosen(), **BASELINE_CHOICES})


def quick_run(data_path, base, methods, *, tf_list=None, reference=None, text_files=(), urls=(), notes=None,
              screen: bool = False, progress=None, orientation: str = "auto") -> dict:
    """Read, clean, run the models and (with a reference) score them, in one call.

    Cleaning is the fixed baseline pipeline (D-021) unless ``screen`` asks the engine to compare
    candidate pipelines, which only changes the choice when a reference network is given.
    """
    from grnbench.autoprep import DEFAULT_SEED, execute_plan

    def tick(frac, msg):
        if progress:
            progress(min(max(frac, 0.0), 1.0), msg)

    work = new_workdir(base)
    tick(0.02, "Reading the data and its context")
    res = analyze(data_path, text_files=text_files, urls=urls, notes=notes, tf_list=tf_list,
                  orientation=orientation)
    diag, plan = res["diag"], res["plan"]
    screening = None
    if screen:
        from grnbench.profile import Profile
        from grnbench.screen import screen_candidates
        from grnbench.selector import select
        tick(0.1, "Comparing candidate pipelines")
        screening = screen_candidates(diag, data_path, work / "screening", tf_list=tf_list,
                                      reference_pairs=reference, selection=select(Profile(**plan.profile)))
        plan = screening["plans"][screening["best"]] if screening["best"] in screening["plans"] else plan
    else:
        plan = baseline_plan(diag, tf_list)
    tick(0.3, "Cleaning")
    records = execute_plan(plan, data_path, work / "clean", tf_list=tf_list, seed=DEFAULT_SEED)
    if not any(r.get("processed") for r in records):
        raise RuntimeError("Nothing was run: the data could not be cleaned. "
                           + "; ".join(f"{r['context']}: {r['error']}" for r in records))
    runs = run_models(records, methods, work / "models", tf_list=tf_list, strategy=plan.value("strategy"),
                      progress=lambda f, msg: tick(0.4 + 0.5 * f, msg))
    eval_ref, scored_on, warning = reference_for_scoring(screening, reference)
    if reference is not None and not reference:
        warning = "the reference file has no edges (only a header row?), so nothing was scored"
    bench = None
    if eval_ref:
        try:
            bench = benchmark(runs, eval_ref, tf_list=tf_list)
        except ValueError as exc:      # the models' results stay; only the score is withheld, with the reason
            warning = str(exc)
    tick(1.0, "Done")
    return {"workdir": work, "diag": diag, "plan": plan, "records": records, "runs": runs, "bench": bench,
            "screen": screening, "reference_used": eval_ref, "reference_warning": warning,
            "scored_on": f"{scored_on} · {len(eval_ref):,} links" if bench is not None and eval_ref else None}
