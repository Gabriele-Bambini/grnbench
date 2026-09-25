"""grnbench interface: clean single-cell data for GRN inference, run models, compare.

Launch with ``grnbench ui`` (or ``streamlit run src/grnbench/ui/app.py``).
Every screen is a thin layer over the library: what the interface does can be
reproduced from the command line with the plan it saves (``plan.json``) and
the ``grnbench prep`` commands it shows.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

from grnbench.ui import logic

st.set_page_config(page_title="grnbench", page_icon=":dna:", layout="wide")
S = st.session_state

STEPS = ["1 · Data", "2 · Profile", "3 · Pipeline", "4 · Compare pipelines",
         "5 · Clean", "6 · Models", "7 · Benchmark"]
DOWNSTREAM = ["screen", "records", "runs", "bench"]


def need(key: str, message: str) -> None:
    if S.get(key) is None:
        st.info(message)
        show_status()
        st.stop()


def reset_downstream(*keys: str) -> None:
    for k in keys:
        S.pop(k, None)


def save_upload(upload, folder: Path) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / upload.name
    path.write_bytes(upload.getbuffer())
    return path


# -- sidebar ------------------------------------------------------------------------ #
st.sidebar.title("grnbench")
st.sidebar.caption("Context-aware data cleaning for gene regulatory network inference")
step = st.sidebar.radio("Workflow", STEPS, key="step")
base = st.sidebar.text_input("Output folder", value=S.get("base", str(Path.home() / "grnbench_runs")))
S["base"] = base
status_box = st.sidebar.container()   # filled after the page ran, so it shows this run's result


def show_status() -> None:
    if S.get("workdir"):
        status_box.caption(f"Run folder: `{S['workdir']}`")
    for label, key in [("dataset analysed", "plan"), ("pipelines compared", "screen"),
                       ("data cleaned", "records"), ("models run", "runs"), ("benchmark", "bench")]:
        status_box.write(("✅ " if S.get(key) is not None else "▫️ ") + label)


# -- 1 data --------------------------------------------------------------------------- #
def page_data() -> None:
    st.header("1 · Data and context")
    st.write("Load raw counts in any common format. README files, download pages and notes are "
             "read to understand the protocol, organism and design.")
    source = st.radio("Dataset", ["Demo dataset", "File or folder path", "Upload a file"], horizontal=True)
    data_path = tf_path = gold_path = None
    if source == "Demo dataset":
        st.caption("Synthetic 10x-like UMI counts with a known regulatory network, two batches, "
                   "three cell types, damaged cells, empty droplets and doublets.")
        data_path = S.get("demo_paths", {}).get("counts")
        if st.button("Create demo dataset"):
            from grnbench.demo import write_demo
            S["demo_paths"] = write_demo(Path(base) / "demo")
            st.rerun()
        if data_path:
            st.success(f"Demo ready: `{data_path}`")
            tf_path, gold_path = S["demo_paths"]["tfs"], S["demo_paths"]["gold"]
    elif source == "File or folder path":
        p = st.text_input("Path to .h5ad / 10x .h5 / 10x MTX folder / .loom / .csv / .tsv")
        data_path = Path(p) if p else None
        if data_path and not data_path.exists():
            st.error("Path not found.")
            data_path = None
    else:
        up = st.file_uploader("Counts file", type=["h5ad", "h5", "loom", "csv", "tsv", "txt", "gz"])
        if up:
            data_path = save_upload(up, Path(base) / "uploads")

    st.subheader("Context (optional)")
    c1, c2 = st.columns(2)
    readmes = c1.file_uploader("README / metadata files", accept_multiple_files=True,
                               type=["md", "txt", "json", "yaml", "yml", "rst"])
    urls = c2.text_area("Download pages or accessions (one per line: GEO GSE..., Zenodo, CELLxGENE ...)")
    notes = st.text_area("Notes about the experiment (protocol, organism, design...)")

    st.subheader("Regulators and reference network (optional)")
    c1, c2 = st.columns(2)
    tf_up = c1.file_uploader("TF list (one symbol per line)", type=["txt", "tsv", "csv"])
    ref_up = c2.file_uploader("Reference network (source/target columns)", type=["tsv", "csv", "txt", "gz"])
    use_dorothea = c2.checkbox("No reference file: download DoRothEA (needs internet + decoupler)")

    if st.button("Analyse dataset", type="primary", disabled=data_path is None):
        work = logic.new_workdir(base)
        text_files = [save_upload(f, work / "context") for f in readmes or []]
        if tf_up:
            tf_path = save_upload(tf_up, work / "inputs")
        if ref_up:
            gold_path = save_upload(ref_up, work / "inputs")
        from grnbench.tfs import load_tf_list
        tfs = load_tf_list(tf_path) if tf_path else None
        reference = None
        with st.spinner("Reading data, context and QC..."):
            if gold_path:
                reference = logic.read_pairs(gold_path)
            elif use_dorothea:
                try:
                    reference = logic.read_pairs("dorotheaABC")
                except Exception as exc:  # noqa: BLE001
                    st.warning(f"DoRothEA not available: {exc}")
            res = logic.analyze(data_path, text_files=text_files,
                                urls=[u for u in urls.splitlines() if u.strip()],
                                notes=notes or None, tf_list=tfs)
        reset_downstream(*DOWNSTREAM)
        S.update(workdir=work, data_path=str(data_path), tfs=tfs, reference=reference,
                 diag=res["diag"], plan=res["plan"], user_choices={})
        st.success("Analysis done - see 2 · Profile and 3 · Pipeline.")


# -- 2 profile ------------------------------------------------------------------------- #
def page_profile() -> None:
    need("plan", "Analyse a dataset first (1 · Data).")
    diag, plan = S["diag"], S["plan"]
    st.header("2 · Dataset profile")
    cards = logic.summarize_profile(plan.profile)
    for i in range(0, len(cards), 5):
        cols = st.columns(5)
        for col, (label, value) in zip(cols, cards[i:i + 5]):
            col.metric(label, value)
    if plan.profile.get("notes"):
        for n in plan.profile["notes"]:
            st.warning(n)

    st.subheader("What the sources say")
    ctx = diag.context
    st.caption("Read: " + ", ".join(ctx.sources_read))
    for f in ctx.sources_failed:
        st.warning(f"Could not read {f['source']}: {f['reason']}")
    for key, votes in ctx.conflicts().items():
        st.warning(f"Sources disagree on **{key}**: {votes} → using `{ctx.get(key)}`")
    st.dataframe(logic.facts_table(ctx), use_container_width=True, hide_index=True)

    st.subheader("Quality control")
    figs = logic.qc_figures(diag, plan)
    cols = st.columns(2)
    for i, (name, fig) in enumerate(figs.items()):
        cols[i % 2].plotly_chart(fig, use_container_width=True, key=f"qc_{name}")
    gates = logic.gates_table(diag.qc)
    st.dataframe(gates.style.map(lambda v: {"fail": "background-color:#f8d7da",
                                            "warn": "background-color:#fff3cd"}.get(v, ""),
                                 subset=["status"]), use_container_width=True, hide_index=True)


# -- 3 pipeline ------------------------------------------------------------------------ #
def page_pipeline() -> None:
    need("plan", "Analyse a dataset first (1 · Data).")
    from grnbench.autoprep import plan_preprocessing
    from grnbench.kb import DECISIONS, EVIDENCE_LEVELS, KB_VERSION, STAGES

    plan, diag = S["plan"], S["diag"]
    st.header("3 · Cleaning pipeline")
    st.caption(f"Chosen by rules from the knowledge base v{KB_VERSION} applied to this dataset's profile. "
               "Every choice shows its reasons, evidence level and references; change any of them.")
    with st.expander("Evidence levels"):
        for k, v in EVIDENCE_LEVELS.items():
            st.write(f"**{k}** - {v}")
    for w in plan.warnings:
        st.warning(w)

    user = S.setdefault("user_choices", {})
    for stage in STAGES:
        decisions = [d for d in DECISIONS if d.stage == stage.id]
        if not decisions:
            continue
        st.markdown(f"#### {stage.order}. {stage.name}")
        for d in decisions:
            res = plan.selection["results"][d.id]
            chosen = plan.choices[d.id]
            flag = " · ⚠ uncertain" if res["uncertain"] else ""
            with st.expander(f"{d.question} → **{chosen}**{flag}  ({d.impact} impact)"):
                st.write(d.why_critical)
                rows = logic.decision_rows(plan, d.id)
                available = rows.loc[rows["available"], "option"].tolist()
                if chosen not in available:
                    available = [chosen] + available
                pick = st.selectbox("Option", available, index=available.index(chosen),
                                    key=f"pick_{d.id}")
                if pick != chosen:
                    user[d.id] = pick
                st.dataframe(rows, use_container_width=True, hide_index=True)
    if user and st.button(f"Apply {len(user)} change(s) and rebuild the plan", type="primary"):
        S["plan"] = plan_preprocessing(diag, tf_list=S.get("tfs"), choices=dict(user))
        reset_downstream(*DOWNSTREAM)
        st.rerun()

    st.subheader("Parameters")
    editable = ["min_genes_per_cell", "max_genes_per_cell", "max_counts_per_cell",
                "max_pct_mito_cell", "min_cells_per_gene", "subsample_cells", "cells_per_metacell"]
    for key in editable:
        if key not in plan.decisions:
            continue
        dec = plan.decisions[key]
        c1, c2 = st.columns([1, 3])
        cur = dec.value
        new = c1.number_input(dec.label, value=float(cur) if cur is not None else 0.0,
                              key=f"num_{key}", help="0 = off" if cur is None else None)
        c2.caption(f"{dec.reason}  _(source: {dec.source})_")
        new_v = None if (cur is None and new == 0) else (int(new) if isinstance(cur, int) or cur is None else new)
        if new_v != cur:
            plan.set(key, new_v, "set by the user in the interface")
            reset_downstream(*DOWNSTREAM)
    st.write("**Cellular contexts:** " + ", ".join(c["name"] for c in plan.contexts))

    st.subheader("Reproduce from a terminal")
    ctx0 = plan.contexts[0]
    cli = plan.cli_args(S["data_path"], "OUT/" + str(ctx0["name"]).replace("=", "-"),
                        "tf_list.txt" if S.get("tfs") else None, context=ctx0,
                        bounds_path="per_batch_bounds.json")
    st.code("grnbench prep " + " ".join(f'"{a}"' if " " in a else a for a in cli), language="bash")
    import json
    st.download_button("Download plan.json", json.dumps(plan.to_dict(), indent=2, default=str),
                       file_name="plan.json")


# -- 4 compare pipelines ----------------------------------------------------------------- #
def page_compare() -> None:
    need("plan", "Analyse a dataset first (1 · Data).")
    from grnbench.kb import DECISIONS
    from grnbench.screen import screen_candidates
    from grnbench.selector import select
    from grnbench.profile import Profile

    plan = S["plan"]
    st.header("4 · Compare candidate pipelines on the data")
    st.write("The recommended pipeline and one-change alternatives are run and measured: retention, "
             "depth and batch leakage, edge stability across cell halves and - with a reference "
             "network - recovery of its **selection half**. The other half is kept for the model "
             "benchmark, so the pipeline is never chosen on the edges it is judged on.")
    selection = select(Profile(**plan.profile))
    auto = selection.candidates(max_candidates=4)
    st.write("Automatic alternatives (uncertain decisions): " +
             (", ".join(c["name"] for c in auto[1:]) or "none - every decision was clear-cut"))

    st.subheader("Add your own alternatives")
    extra = S.setdefault("extra_candidates", [])
    c1, c2, c3 = st.columns([2, 2, 1])
    dec = c1.selectbox("Decision", [d.id for d in DECISIONS])
    opts = [o.id for o in next(d for d in DECISIONS if d.id == dec).options if o.availability()[0]]
    opt = c2.selectbox("Option", opts)
    if c3.button("Add"):
        extra.append({"name": f"{dec}={opt}", "choices": {dec: opt}})
    if extra:
        st.write("Custom: " + ", ".join(e["name"] for e in extra))
        if st.button("Clear custom alternatives"):
            extra.clear()
            st.rerun()

    if st.button("Run and compare", type="primary"):
        with st.spinner("Running every candidate pipeline..."):
            S["screen"] = screen_candidates(S["diag"], S["data_path"], Path(S["workdir"]) / "screening",
                                            tf_list=S.get("tfs"), reference_pairs=S.get("reference"),
                                            selection=selection, extra=list(extra))
    res = S.get("screen")
    if res is None:
        return
    import plotly.express as px
    rows = res["rows"]
    st.success(f"Best: **{res['best']}** - rule: {res['rule']}")
    st.dataframe(res["ranking"].round(3), use_container_width=True)
    metrics = [m for m in ["reference_auprc_ratio", "edge_stability", "depth_leakage", "batch_leakage",
                           "n_cells"] if m in rows and rows[m].notna().any()]
    cols = st.columns(2)
    for i, m in enumerate(metrics):
        fig = px.bar(rows, x="candidate", y=m, color="context", barmode="group", height=320)
        cols[i % 2].plotly_chart(fig, use_container_width=True, key=f"scr_{m}")
    st.subheader("Pros and cons against the recommended pipeline")
    for name, pc in res["pros_cons"].items():
        if name == "recommended":
            continue
        st.markdown(f"**{name}** - pros: {', '.join(pc['pros']) or '-'}; cons: {', '.join(pc['cons']) or '-'}")
    adopt = st.selectbox("Pipeline to use for cleaning", list(res["plans"]),
                         index=list(res["plans"]).index(res["best"]) if res["best"] in res["plans"] else 0)
    if st.button("Use this pipeline"):
        S["plan"] = res["plans"][adopt]
        reset_downstream("records", "runs", "bench")
        st.success(f"Pipeline '{adopt}' adopted - go to 5 · Clean.")


# -- 5 clean --------------------------------------------------------------------------- #
def page_clean() -> None:
    need("plan", "Analyse a dataset first (1 · Data).")
    from grnbench.autoprep import execute_plan

    plan = S["plan"]
    st.header("5 · Clean the data")
    st.write(f"Strategy **{plan.value('strategy')}**, contexts: "
             + ", ".join(c["name"] for c in plan.contexts))
    if st.button("Run cleaning", type="primary"):
        with st.spinner("Cleaning..."):
            S["records"] = execute_plan(plan, S["data_path"], Path(S["workdir"]) / "clean",
                                        tf_list=S.get("tfs"))
        reset_downstream("runs", "bench")
    recs = S.get("records")
    if not recs:
        return
    for r in recs:
        (st.error if r["error"] else st.success)(f"{r['context']}: {r['error'] or r['processed']}")
    st.dataframe(logic.retention_table(recs), use_container_width=True, hide_index=True)
    with st.expander("Exact commands"):
        for r in recs:
            st.code(" ".join(f'"{a}"' if " " in a else a for a in r["cli"]), language="bash")


# -- 6 models --------------------------------------------------------------------------- #
def page_models() -> None:
    need("records", "Clean the data first (5 · Clean).")
    from grnbench.infer import list_methods

    plan = S["plan"]
    st.header("6 · Models")
    rec = {m["name"]: m for m in plan.methods}
    table = pd.DataFrame([{"method": m["name"], "available": m["available"], "origin": m["origin"],
                           "recommended": rec.get(m["name"], {}).get("recommended", False),
                           "note": rec.get(m["name"], {}).get("reason", "") or m.get("unavailable_reason") or ""}
                          for m in list_methods()])
    st.dataframe(table, use_container_width=True, hide_index=True)

    st.subheader("Your model")
    st.caption("A .py file with a class exposing `infer(adata, tf_list, **params)` (or a function "
               "`f(adata, tf_list)`) returning columns source, target, weight. It gets exactly the "
               "same cleaned data as the built-in models.")
    c1, c2, c3 = st.columns([2, 1, 1])
    up = c1.file_uploader("Model file", type=["py"])
    attr = c2.text_input("Class or function name")
    if c3.button("Register", disabled=not (up and attr)):
        path = save_upload(up, Path(S["workdir"]) / "user_models")
        try:
            reg = logic.register_user_model(path, attr)
            S.setdefault("user_methods", set()).update(r["name"] for r in reg)
            st.success("Registered: " + ", ".join(r["name"] for r in reg))
            st.rerun()
        except Exception as exc:  # noqa: BLE001
            st.error(f"Could not register the model: {exc}")

    avail = table.loc[table["available"], "method"].tolist()
    default = [m for m in avail if rec.get(m, {}).get("recommended")] + sorted(S.get("user_methods", set()))
    methods = st.multiselect("Models to run", avail, default=[m for m in default if m in avail])
    if st.button("Run models", type="primary", disabled=not methods):
        bar = st.progress(0.0)
        S["runs"] = logic.run_models(S["records"], methods, Path(S["workdir"]) / "models",
                                     tf_list=S.get("tfs"), strategy=plan.value("strategy"),
                                     progress=lambda f, msg: bar.progress(f, text=msg))
        reset_downstream("bench")
    runs = S.get("runs")
    if runs is not None:
        st.dataframe(runs.drop(columns=["processed"]), use_container_width=True, hide_index=True)


# -- 7 benchmark ------------------------------------------------------------------------- #
def page_benchmark() -> None:
    need("runs", "Run the models first (6 · Models).")
    st.header("7 · Benchmark")
    screen = S.get("screen")
    if screen and screen["reference"]["evaluation"]:
        ref, which = screen["reference"]["evaluation"], "evaluation half (never used to choose the pipeline)"
    else:
        ref, which = S.get("reference"), "full reference (the pipeline was not chosen on it)"
    up = st.file_uploader("Or score against another reference network", type=["tsv", "csv", "txt", "gz"])
    if up:
        ref, which = logic.read_pairs(save_upload(up, Path(S["workdir"]) / "inputs")), f"uploaded: {up.name}"
    if not ref:
        st.info("No reference network: load one in 1 · Data or upload it here.")
        show_status()
        st.stop()
    st.caption(f"Reference: {which} - {len(ref)} edges.")
    if st.button("Score models", type="primary"):
        S["bench"] = logic.benchmark(S["runs"], ref, tf_list=S.get("tfs"))
    bench = S.get("bench")
    if bench is None or bench.empty:
        return
    import plotly.express as px
    user = S.get("user_methods", set())
    bench = bench.assign(kind=["your model" if m in user else "standard" for m in bench["method"]])
    st.dataframe(bench.round(3), use_container_width=True, hide_index=True)
    for metric in ["AUPRC ratio", "EPR", "AUROC"]:
        fig = px.bar(bench, x="method", y=metric, color="context", pattern_shape="kind",
                     barmode="group", height=340)
        if metric in ("AUPRC ratio", "EPR"):
            fig.add_hline(y=1, line_dash="dot", annotation_text="random")
        st.plotly_chart(fig, use_container_width=True, key=f"bench_{metric}")
    st.caption("Contexts from one dataset are not independent replicates: differences here describe "
               "this dataset; claims across datasets need `grnbench compare` over independent units.")
    st.download_button("Download scores (CSV)", bench.to_csv(index=False), file_name="benchmark.csv")
    st.download_button("Download the whole run (zip)", logic.zip_folder(S["workdir"]),
                       file_name=f"{Path(S['workdir']).name}.zip")


{STEPS[0]: page_data, STEPS[1]: page_profile, STEPS[2]: page_pipeline, STEPS[3]: page_compare,
 STEPS[4]: page_clean, STEPS[5]: page_models, STEPS[6]: page_benchmark}[step]()
show_status()
