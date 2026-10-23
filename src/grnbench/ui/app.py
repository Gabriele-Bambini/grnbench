"""grnbench interface: clean single-cell data for GRN inference, run models, compare.

Launch with ``grnbench ui`` (or ``streamlit run src/grnbench/ui/app.py``).
Every screen is a thin layer over the library: a cleaning can be reproduced
from the command line with the ``grnbench prep`` commands the interface shows
(the Clean screen shows one command per context). They use the work-folder
paths, so they reproduce on the same machine or after the paths are rewritten.
The plan is written as ``plan.json`` in the cleaning output folder (and offered
as a download on the pipeline screen) before any context runs: it records the
choices, not that they succeeded; the per-context outcome and commands are in
``contexts.json``. No grnbench command or screen reads ``plan.json`` back (the
library can: ``PrepPlan.load``).
"""
from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

from grnbench import __version__
from grnbench.ui import catalog, logic, steps, style
from grnbench.ui.steps import (BENCHMARK, CLEAN, COMPARE, DATA, EXPLORE_DATA, EXPLORE_RESULTS, MODELS, PIPELINE,
                               PROFILE, QUICK)

st.set_page_config(
    page_title="grnbench", page_icon=":material/hub:", layout="wide", initial_sidebar_state="auto",
    menu_items={"About": f"**grnbench {__version__}** - context-aware data cleaning for gene regulatory network "
                         "inference. Every number on these pages comes from the library; the commands shown "
                         "next to the results reproduce them."})
S = st.session_state
# Streamlit forgets the value of every widget that is not drawn in a run, so changing page used to wipe what
# had been typed or ticked. Writing a key back before its widget is drawn keeps it (no streamlit version needed).
KEPT_INPUTS = ("quick_source", "quick_pathbox", "quick_orientation", "quick_preset", "quick_screen",
               "data_source", "data_pathbox", "data_orientation", "data_urls", "data_notes", "data_dorothea",
               "models_preset")
for _key in KEPT_INPUTS:
    if _key in S:
        S[_key] = S[_key]
S.setdefault("step", QUICK)

DOWNSTREAM = ["screen", "records", "runs", "bench", "bench_on", "reference_warning"]
TYPED = "set by the user in the interface"      # the reason of a value typed in a box of step 3


class PageDone(Exception):
    """Raised by a screen that has nothing more to draw (``st.stop()`` would also skip the sidebar menu)."""


def goto(step: str) -> None:
    """Button callback: show another screen on the next run."""
    S["step"] = step


def empty_state(icon: str, title: str, text: str, primary: tuple[str, str], secondary: tuple[str, str] | None = None):
    """A screen that has nothing to show yet says why and offers the way forward, then stops the page."""
    with st.container(border=True, key="empty-state"):
        st.markdown(f"## :material/{icon}:")
        st.markdown(f"### {title}")
        st.markdown(text)
        _, c1, c2, _ = st.columns([1, 1.4, 1.4, 1] if secondary else [1.2, 1.6, 0.01, 1.2])
        c1.button(primary[0], on_click=goto, args=(primary[1],), type="primary", key="empty-primary",
                  width="stretch")
        if secondary:
            c2.button(secondary[0], on_click=goto, args=(secondary[1],), key="empty-secondary", width="stretch")
    raise PageDone


def reset_session() -> None:
    """Forget this session's data, plan and results (the output folder stays); files on disk are not touched."""
    for key in [k for k in S.keys() if k != "base"]:
        del S[key]
    S["step"] = QUICK


def need_data() -> None:
    if S.get("plan") is None:
        empty_state("upload_file", "No dataset yet",
                    "Load your data first, or try the built-in demo dataset: it takes one click.",
                    ("Load my data", DATA), ("Try the demo", QUICK))


def need_cleaned() -> None:
    if S.get("records") is None:
        empty_state("cleaning_services", "The data has not been cleaned yet",
                    "Models run on cleaned data. Clean it here, or let the quick run do everything at once.",
                    ("Clean the data", CLEAN), ("Quick run", QUICK))


def need_models() -> None:
    if S.get("runs") is None:
        empty_state("hub", "No model has been run yet",
                    "Run at least one model to see scores and networks here.",
                    ("Run models", MODELS), ("Quick run", QUICK))


def footer(back: str | None = None, forward: tuple[str, str, bool] | None = None, alt: tuple[str, str] | None = None) -> None:
    """Bottom of a workflow screen: where you came from, an alternative, where to go next (``(label, step, ready)``)."""
    st.markdown('<div class="gb-foot"></div>', unsafe_allow_html=True)
    left, mid, right = st.columns([1.2, 1.5, 1.5])
    if back:
        left.button(f"Back: {back}", on_click=goto, args=(back,), icon=":material/arrow_back:", key="foot-back",
                    width="stretch")
    if alt:
        mid.button(alt[0], on_click=goto, args=(alt[1],), key="foot-alt", width="stretch")
    if forward:
        label, target, ready = forward
        right.button(label, on_click=goto, args=(target,), icon=":material/arrow_forward:", icon_position="right",
                     type="primary" if ready else "secondary", key="foot-next", width="stretch")


def reset_downstream(*keys: str) -> None:
    for k in keys:
        S.pop(k, None)
    for k in [k for k in S.keys() if k.startswith("_") and k.endswith("_cache")]:
        S.pop(k, None)       # explorer caches: whatever was drawn from the old data or edges is stale


def cached(name: str, key, make, what: str = "Working..."):
    """Compute ``make()`` once per ``key`` in this session (the explorer tabs redraw on every widget change).

    The spinner is drawn in the sidebar: the heavy blocks run inside tab panes that are not on screen."""
    cache = S.setdefault(f"_{name}_cache", {})
    if key not in cache:
        with status_box, st.spinner(what):
            cache[key] = make()
    return cache[key]


@contextmanager
def guarded(what: str):
    """One heavy tab failing must not blank the page: say what failed and let the other tabs work."""
    try:
        yield
    except Exception as exc:  # noqa: BLE001 - shown to the person, not hidden
        st.error(f"Could not {what}: {type(exc).__name__}: {exc}")


def save_upload(upload, folder: Path) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / upload.name
    path.write_bytes(upload.getbuffer())
    return path


def read_upload(what: str, parse, path):
    """``parse(path)`` for an uploaded TF list or reference network; a bad file gets one red line, not a traceback."""
    try:
        return parse(path)
    except Exception as exc:  # noqa: BLE001 - UnicodeDecodeError, EmptyDataError, TFListError, missing columns, a corrupt .gz ...
        st.error(f"Could not read {Path(path).name} as {what}: {type(exc).__name__}: {exc}")
        raise PageDone from None



SOURCES = ["Demo dataset", "Upload a file", "File or folder path"]
PRESET_HELP = {
    "Fast (seconds)": "Correlation and information-theory models that finish in seconds. The right start for a first look.",
    "Fast, including extras": "The same, plus the extra fast models that have not been scored on real data yet.",
    "Recommended by the engine": "The engine's own pick for this dataset. It can include models that take many minutes.",
    "Everything available": "Every model that can run on this computer. This can take hours on a real dataset.",
    "None": "Tick the models you want yourself in the list below.",
}


def pick_dataset(key: str):
    """Dataset chooser shared by the quick run and step 1: ``(data_path, tf_path, gold_path)``."""
    source = st.radio("Dataset", SOURCES, horizontal=True, key=f"{key}_source", label_visibility="collapsed")
    data_path = tf_path = gold_path = None
    if source == "Demo dataset":
        data_path = S.get("demo_paths", {}).get("counts")
        if data_path:
            tf_path, gold_path = S["demo_paths"]["tfs"], S["demo_paths"]["gold"]
            st.success(f"Demo dataset ready: {S.get('demo_info', 'synthetic counts')}.", icon=":material/check_circle:")
        else:
            st.markdown("A small synthetic dataset with a known regulatory network: the quickest way to see "
                        "everything work. It has two batches, three cell types and deliberately damaged cells, "
                        "empty droplets and doublets, so you can watch the cleaning do its job.")
            if st.button("Create the demo dataset", key=f"{key}_demo", type="primary", icon=":material/science:"):
                from grnbench.demo import write_demo
                with st.spinner("Writing the demo dataset..."):
                    S["demo_paths"] = write_demo(Path(base) / "demo")
                    S["demo_info"] = logic.demo_summary(S["demo_paths"])
                st.rerun()
    elif source == "Upload a file":
        up = st.file_uploader("Counts file", type=["h5ad", "h5", "loom", "csv", "tsv", "txt", "gz"], key=f"{key}_upload",
                              help="h5ad, 10x .h5, loom, csv, tsv or txt (also gzipped), up to 200 MB. For a larger "
                                   "file, or a 10x folder, choose 'File or folder path'.")
        if up:
            data_path = save_upload(up, Path(base) / "uploads")
    else:
        text = st.text_input("Path to the data", key=f"{key}_pathbox", placeholder="~/data/my_dataset.h5ad",
                             help="A .h5ad, 10x .h5, .loom, .csv or .tsv file, or a 10x folder (matrix.mtx, barcodes, "
                                  "features), on the computer that runs grnbench.")
        data_path = logic.clean_path(text)
        if data_path and not data_path.exists():
            st.error("There is nothing at that path: check the spelling, and that it is on the computer running "
                     "grnbench.", icon=":material/folder_off:")
            data_path = None
    if source != "Demo dataset":
        with st.expander("Is it a CSV / TSV table? Say whether rows are cells or genes", icon=":material/table_rows:"):
            st.selectbox("CSV / TSV layout", ["auto", "genes_x_cells", "cells_x_genes"], key=f"{key}_orientation",
                         help="Only for .csv/.tsv tables. auto detects it and refuses to guess; pick genes_x_cells "
                              "(rows are genes, the usual GEO layout) or cells_x_genes when it asks you to.")
    return data_path, tf_path, gold_path


def model_picker(key: str, recommended=(), *, all_regulators: bool = False) -> list[str]:
    """A preset, a one-line summary and, in an expander, the full table; returns the available models ticked."""
    from grnbench.infer import list_methods

    table = catalog.model_table(list_methods(), recommended, S.get("user_methods", set()),
                                all_regulators=all_regulators)
    # Always open on the fast preset (like the Quick run): the engine's own recommendation includes tree
    # ensembles that take tens of minutes per context, and this page is busy until every job ends.
    preset = st.radio("Start from", catalog.PRESETS, horizontal=True, key=f"{key}_preset", index=0)
    st.caption(PRESET_HELP[preset])
    ticked = set(catalog.preset_selection(table, preset, recommended))
    if preset != "None":     # your own models ride along with every preset except "None"
        ticked |= set(S.get("user_methods", set())) & set(table.loc[table["available"], "method"])
    shown = table.assign(run=table["method"].isin(ticked))
    with st.expander("Choose individual models", icon=":material/tune:"):
        tiers = st.multiselect("Show", list(catalog.TIERS), default=list(catalog.TIERS), key=f"{key}_tiers",
                               help="default: the engine's default set (what `--methods available` runs); extra: "
                               "registered but off by default, never scored on real data; yours: a model you registered")
        shown = shown[shown["tier"].isin(tiers)]
        edited = st.data_editor(
            shown, hide_index=True, key=f"{key}_editor_{preset}_{'-'.join(tiers)}",
            disabled=[c for c in shown.columns if c != "run"], height=min(60 + 35 * len(shown), 700),
            column_order=["run", "method", "family", "tier", "speed", "available", "why not available", "what it is"],
            column_config={"run": st.column_config.CheckboxColumn("run", help="tick the models to run"),
                           "what it is": st.column_config.TextColumn(width="large")})
        st.caption("Speeds are indicative for about 1,000 genes x 1,500 cells and grow with the number of genes. "
                   "Models marked *extra* were never scored on real data.")
    chosen = edited.loc[edited["run"] & edited["available"], "method"].tolist()
    chosen += [m for m in ticked if m not in set(shown["method"]) and m not in chosen]
    words = {"seconds": "fast", "about a minute": "about a minute each", "minutes or more": "slow (minutes or more)",
             "unknown speed": "of unknown speed"}
    speeds = table.set_index("method").loc[chosen, "speed"].value_counts() if chosen else pd.Series(dtype=int)
    slow = [m for m in chosen if table.set_index("method").loc[m, "speed"] == catalog.SPEED_LABEL[2]]
    if slow:
        st.warning(f"{', '.join(slow)}: minutes or more for each context (tens of minutes to hours on thousands of "
                   "cells), and this page is busy until every job ends. Untick them for a first look.")
    st.caption(f"{logic.plural(len(chosen), 'model')} selected"
               + (": " + ", ".join(f"{n} {words[k]}" for k, n in speeds.items()) if len(speeds) else "") + ".")
    return chosen


# -- sidebar ------------------------------------------------------------------------ #
GLOSSARY = {
    "TF (transcription factor)": "a gene whose protein switches other genes on or off: the possible regulators of a network.",
    "Network model": "a method that ranks TF -> gene links by how strongly the data suggest the TF regulates the gene.",
    "Reference network": "TF -> gene links known from the literature or from experiments, used to score the models.",
    "Context": "the cells one network is built from: all cells pooled, or one cell type at a time.",
    "Selection half / evaluation half": "the reference is split in two: one half helps choose the cleaning pipeline, "
                                        "the other half is only used to score the models, so they are never graded "
                                        "on links the pipeline was picked on.",
    "AUPRC ratio": "area under the precision-recall curve divided by what a random ranking gets. 1 = random; "
                   "10 = ten times better than random at putting true links first.",
    "EPR": "early precision ratio: the share of true links among the top k predictions (k = number of reference "
           "links found), divided by the random value.",
    "AUROC": "area under the ROC curve. 0.5 = random, 1 = perfect.",
    "Batch": "cells processed separately (run, donor, lane...) whose technical differences can look like biology.",
    "Doublet": "a droplet that captured two cells: it looks like one cell with mixed identity.",
}


def environment_panel() -> None:
    """Optional tools of the engine (``grnbench env``). Read-only, and R starts only when the button is pressed."""
    st.caption("The cleaning engine's optional tools: R and its packages, Python and command-line extras. "
               "grnbench installs nothing: a missing tool shows the command to run yourself.")
    if st.button("Check R", key="env-check", icon=":material/search:", width="stretch",
                 help="Starts R to list its packages (about a second). Nothing starts R until you press it."):
        from grnbench import rbridge
        with st.spinner("Asking R which packages it has..."):
            S["env_report"] = rbridge.environment_report(fresh=True)    # pressing the button means: ask R now
    rows = S.get("env_report")
    if rows is None:
        st.caption("Not checked in this session yet.")
        return
    no_r = not rows[0]["found"]                     # the first row is R itself
    for row in rows:
        mark = ":green[:material/check_circle:]" if row["found"] else ":red[:material/cancel:]"
        detail = " · ".join(x for x in (row["version"], row["source"]) if x)
        st.markdown(f"{mark} **{row['tool']}**" + (f" {detail}" if detail else ""))
        needed = "needed for: " + (", ".join(row["used_for"]) or "-")
        if row["found"]:
            if row["install"]:                      # R found, but other candidates were skipped: say why
                st.code(row["install"], language=None)   # paths verbatim, not read as markdown
        elif no_r and row["tool"].startswith("r:"):
            st.caption(f"{needed}; install R first (above)")
        elif row["install"]:
            st.caption(needed)
            st.code(row["install"], language="r" if row["tool"].startswith("r:") else "bash")


with st.sidebar:
    st.markdown(style.brand_html(), unsafe_allow_html=True)
    nav_slot = st.container()          # filled after the page ran, so the ticks show this run's progress
    st.divider()
    with st.popover("Glossary", icon=":material/menu_book:", width="stretch"):
        for term, meaning in GLOSSARY.items():
            st.markdown(f"**{term}**: {meaning}")
    with st.popover("Environment", icon=":material/memory:", width="stretch"):
        environment_panel()
    S.setdefault("base", str(Path.home() / "grnbench_runs"))
    with st.expander("Settings", icon=":material/settings:"):
        base = st.text_input("Output folder", key="base",
                             help="Demo data, uploads and every run are written here. Your own files are never changed.")
        if S.get("workdir"):
            st.caption(f"Current run: `{Path(str(S['workdir'])).name}`")
        st.button("Start over", on_click=reset_session,
                  icon=":material/restart_alt:", key="start-over", width="stretch",
                  help="Forget the data, plan and results of this session. Files already written stay on disk.")
    status_box = st.container()        # spinners of the explorer tabs appear here
base_error = None
try:
    base = str(logic.output_folder(base))     # expanded, created and checked: every later use can trust it
except ValueError as exc:
    base_error = str(exc)


def step_done(key: str) -> bool:
    """A step counts as done when it produced something: cleaned data, edges, scores."""
    value = S.get(key)
    if value is None:
        return False
    if key == "records":
        return any(r.get("processed") for r in value)
    if key == "runs":
        return "edges" in value and bool(value["edges"].notna().any())
    if hasattr(value, "empty"):
        return not value.empty
    return True


def mark_seen(step: str) -> None:
    """Review screens (profile, pipeline) count as done once their content has been on screen."""
    S.setdefault("seen", set()).add(step)


def is_done(step: str) -> bool:
    return {QUICK: lambda: bool(S.get("quick_done")), DATA: lambda: step_done("plan"),
            PROFILE: lambda: PROFILE in S.get("seen", ()), PIPELINE: lambda: PIPELINE in S.get("seen", ()),
            COMPARE: lambda: step_done("screen"), CLEAN: lambda: step_done("records"),
            MODELS: lambda: step_done("runs"), BENCHMARK: lambda: step_done("bench")}.get(step, lambda: False)()


def render_nav() -> None:
    """The sidebar menu: three groups, a tick on what is finished, the current screen highlighted."""
    with nav_slot:
        group = None
        for i, item in enumerate(steps.STEPS):
            if item.group != group:
                group = item.group
                st.markdown(style.group_html(group), unsafe_allow_html=True)
            state = "cur" if item.key == S["step"] else ("done" if is_done(item.key) else "todo")
            if item.group == steps.WORKFLOW:
                icon = "check_circle" if state == "done" else ("arrow_circle_right" if state == "cur" else "radio_button_unchecked")
            else:
                icon = item.icon
            label = item.key + ("  :gray[optional]" if item.optional else "")
            st.button(label, key=f"nav-{state}-{i}", icon=f":material/{icon}:", on_click=goto, args=(item.key,),
                      width="stretch")


# -- 1 data --------------------------------------------------------------------------- #
def page_data() -> None:
    with st.container(border=True, key="card-data"):
        st.markdown(style.card_title_html(1, "Dataset", "Where are your raw counts?"), unsafe_allow_html=True)
        data_path, tf_path, gold_path = pick_dataset("data")
    with st.expander("Add context: README files, download pages, notes (optional)", icon=":material/description:"):
        st.caption("Text the engine reads to recognise the protocol, the organism and the design. It is combined with "
                   "what the count values themselves say, and the two are compared on the Profile page.")
        c1, c2 = st.columns(2)
        readmes = c1.file_uploader("README or metadata files", accept_multiple_files=True, key="data_readmes",
                                   type=["md", "txt", "json", "yaml", "yml", "rst"])
        urls = c2.text_area("Download pages or accessions", key="data_urls", placeholder="GSE123456",
                            help="One per line: GEO GSE..., Zenodo, CELLxGENE ...")
        notes = st.text_area("Notes about the experiment", key="data_notes",
                             placeholder="10x Chromium v3, human PBMC, two donors, sequenced in two runs...")
    with st.container(border=True, key="card-regulators"):
        demo_has_both = tf_path is not None and gold_path is not None
        st.markdown(style.card_title_html(2, "Regulators and reference network",
                                          "The demo dataset brings both: leave the boxes empty to use them." if demo_has_both
                                          else "Both are optional."), unsafe_allow_html=True)
        c1, c2 = st.columns(2)
        tf_up = c1.file_uploader("TF list", type=["txt", "tsv", "csv"], key="data_tf",
                                 help="One gene symbol per line (.txt), or a .csv/.tsv with a symbol column.")
        c1.caption("The genes allowed to act as regulators.")
        ref_up = c2.file_uploader("Reference network", type=["tsv", "csv", "txt", "gz"], key="data_ref",
                                  help="Known regulator -> target links, in two columns named source and target.")
        c2.caption("Known regulator -> target links, used to score models and to compare pipelines.")
        use_dorothea = c2.checkbox("No file? Download DoRothEA", key="data_dorothea", disabled=bool(gold_path or ref_up),
                                   help="A curated human reference network. Needs internet and the decoupler package. "
                                        "Not used when a reference network is already given.")
    if st.button("Analyse dataset", type="secondary" if S.get("plan") is not None else "primary", disabled=data_path is None,
                 icon=":material/search_insights:"):
        work = logic.new_workdir(base)
        text_files = [save_upload(f, work / "context") for f in readmes or []]
        if tf_up:
            tf_path = save_upload(tf_up, work / "inputs")
        if ref_up:
            gold_path = save_upload(ref_up, work / "inputs")
        from grnbench.tfs import load_tf_list
        tfs = read_upload("a TF list (one gene symbol per line)", load_tf_list, tf_path) if tf_path else None
        reference = None
        with st.spinner("Reading data, context and QC..."):
            if gold_path:
                reference = read_upload("a reference network (a source and a target column)", logic.read_pairs, gold_path)
            elif use_dorothea:
                try:
                    reference = logic.read_pairs("dorotheaABC")
                except Exception as exc:  # noqa: BLE001
                    st.warning(f"DoRothEA not available: {exc}")
            res = logic.analyze(data_path, text_files=text_files,
                                urls=[u for u in urls.splitlines() if u.strip()],
                                notes=notes or None, tf_list=tfs, orientation=S.get("data_orientation", "auto"))
        reset_downstream(*DOWNSTREAM)
        S.update(workdir=work, data_path=str(data_path), tfs=tfs, reference=reference,
                 diag=res["diag"], plan=res["plan"], user_choices={})
        st.success("Analysis done. Next: check what the engine found about your data.", icon=":material/task_alt:")
    elif data_path is None:
        st.caption("Choose a dataset in card 1 to enable the analysis.")
    if S.get("plan") is not None:
        facts = logic.summarize_profile(S["plan"].profile)[:4]
        st.markdown(f"**Currently analysed:** `{Path(str(S['data_path'])).name}`")
        cols = st.columns(4)
        for col, (label, value) in zip(cols, facts):
            col.metric(label, value, border=True)
    footer(forward=(f"Continue: {PROFILE}", PROFILE, S.get("plan") is not None))


# -- 2 profile ------------------------------------------------------------------------- #
QC_NOTES = {
    "detected genes": "Genes detected per cell. Cells left of the dashed line have too few and are removed as empty or damaged.",
    "log10 counts per cell": "Total counts per cell (log10). Cells right of the dashed line, if any, are removed as likely multiplets.",
    "% mitochondrial": "Share of counts from mitochondrial genes. High values mean a damaged cell; cells right of the dashed line are removed.",
    "composition": "Cells per cell type and batch. If a cell type sits in one batch only, batch and biology cannot be told apart.",
}
STATUS_ICON = {"pass": ":green[:material/check_circle:]", "warn": ":orange[:material/warning:]",
               "fail": ":red[:material/error:]"}


def page_profile() -> None:
    need_data()
    mark_seen(PROFILE)
    diag, plan = S["diag"], S["plan"]
    verdict = logic.quality_verdict(diag.qc)
    {"ok": st.success, "warn": st.warning, "fail": st.error}[verdict["level"]](
        verdict["headline"], icon={"ok": ":material/check_circle:", "warn": ":material/warning:",
                                   "fail": ":material/error:"}[verdict["level"]])
    for item in verdict["attention"]:
        st.markdown(f"{STATUS_ICON[item['status']]} **{item['check']}**: {item['detail']}"
                    + (f"  \n:gray[{item['meaning']}]" if item["meaning"] else ""))
    for n in plan.profile.get("notes") or []:
        st.warning(n)

    cards = dict(logic.summarize_profile(plan.profile))
    key_facts = ["Cells", "Genes", "Organism", "Protocol", "Batches"]
    cols = st.columns(len(key_facts))
    for col, label in zip(cols, key_facts):
        col.metric(label, cards[label], border=True)
    with st.expander("All measurements", icon=":material/straighten:"):
        rest = [(k, v) for k, v in cards.items() if k not in key_facts]
        for i in range(0, len(rest), 4):
            cols = st.columns(4)
            for col, (label, value) in zip(cols, rest[i:i + 4]):
                col.metric(label, value)

    ctx = diag.context
    for f in ctx.sources_failed:
        st.warning(f"Could not read {f['source']}: {f['reason']}")
    for key, votes in ctx.conflicts().items():
        st.warning(f"Sources disagree on **{key}**: {votes} → using `{ctx.get(key)}`")

    st.subheader("Quality control", anchor=False)
    figs = logic.qc_figures(diag, plan)
    cols = st.columns(2)
    for i, (name, fig) in enumerate(figs.items()):
        cols[i % 2].plotly_chart(fig, key=f"qc_{name}")
        cols[i % 2].caption(QC_NOTES.get(name, ""))
    checks = logic.quality_checks(diag.qc)
    with st.expander(f"All {len(checks)} quality checks", icon=":material/fact_check:"):
        for _, row in checks.iterrows():
            st.markdown(f"{STATUS_ICON[row['status']]} **{row['check']}**: {row['detail']}"
                        + (f"  \n:gray[{row['meaning']}]" if row["meaning"] and row["status"] != "pass" else ""))
    with st.expander("What the sources say about this dataset", icon=":material/menu_book:"):
        st.caption("Read: " + ", ".join(ctx.sources_read))
        st.dataframe(logic.facts_table(ctx), hide_index=True)
    footer(back=DATA, forward=(f"Continue: {PIPELINE}", PIPELINE, True))


# -- 3 pipeline ------------------------------------------------------------------------ #
def option_cards(rows: pd.DataFrame) -> None:
    """Every option of a decision as a short block: score, evidence, why it scored so, pros, cons, references."""
    for _, r in rows.iterrows():
        mark = ":green[:material/check_circle:] " if r["chosen"] else ""
        head = f"{mark}**{r['name']}** · score {r['score']:g} · evidence: {r['evidence']}"
        if not r["available"]:
            head += f" · :gray[not available: {r['why unavailable']}]"
        st.markdown(head)
        bits = [f"{label}: {r[col]}" for label, col in (("Why", "reasons"), ("Pros", "pros"), ("Cons", "cons"),
                                                         ("References", "references")) if r[col]]
        if bits:
            st.caption("  \n".join(bits))


def decision_editor(plan, d, user: dict, *, expanded: bool = False) -> None:
    """One decision of the plan: its question and answer, the reasons, and a box to choose another option."""
    res = plan.selection["results"][d.id]
    chosen = plan.choices[d.id]
    flag = "  :orange[unsure]" if res["uncertain"] else ""
    with st.expander(f"{d.question} → **{d.option(chosen).name}**{flag}", expanded=expanded):
        st.write(d.why_critical)
        rows = logic.decision_rows(plan, d.id)
        available = rows.loc[rows["available"], "option"].tolist()
        if chosen not in available:
            available = [chosen] + available
        pick = st.selectbox("Option", available, index=available.index(chosen), key=f"pick_{d.id}",
                            format_func=lambda o: d.option(o).name)
        if pick != chosen:
            user[d.id] = pick
        option_cards(rows)


def box_value(key: str, value):
    """What the box of a threshold shows for a plan value: 0 where the plan has none (off); all but the
    mitochondrial percentage are counts."""
    return float(value or 0.0) if key == "max_pct_mito_cell" else int(value or 0)


def box_off(dec) -> bool:
    """0 = off in the box of a decision where the plan or the engine has no bound (a 0 bound drops every cell)."""
    engine = dec.auto_value if dec.auto_reason is not None else dec.value
    return dec.value is None or engine is None


def from_box(dec, number):
    """The plan value a number in the box of a decision stands for (the inverse of box_value)."""
    if box_off(dec) and number == 0:
        return None
    return float(number) if dec.key == "max_pct_mito_cell" else int(number)


def page_pipeline() -> None:
    need_data()
    mark_seen(PIPELINE)
    from grnbench.autoprep import plan_preprocessing
    from grnbench.kb import DECISIONS, EVIDENCE_LEVELS, KB_VERSION, STAGES

    plan, diag = S["plan"], S["diag"]
    st.caption(f"Chosen by rules from the knowledge base v{KB_VERSION}, applied to this dataset's profile. "
               "Every choice shows its reasons, evidence level and references.")
    if plan.warnings:
        said = [logic.plain_warning(w) for w in plan.warnings]
        st.warning(said[0] if len(said) == 1 else "\n".join(f"- {w}" for w in said))
        notes = logic.option_notes(plan)
        if notes:
            st.caption("  \n".join(notes))

    overview = logic.pipeline_overview(plan)
    unsure = [r for r in overview if r["uncertain"]]
    st.subheader("At a glance", anchor=False)
    if unsure:
        st.markdown(f":orange[:material/help:] **The engine is unsure about {len(unsure)} of {len(overview)} decisions** "
                    "(below, under *Worth a look*). The others are clear-cut for this data.")
    else:
        st.markdown(f":green[:material/check_circle:] **All {len(overview)} decisions are clear-cut for this data.**")
    glance = pd.DataFrame([{"stage": r["stage"], "decision": r["question"], "chosen": r["choice"],
                            "confidence": "unsure" if r["uncertain"] else "clear", "impact": r["impact"]}
                           for r in overview])
    st.dataframe(glance.style.map(lambda v: "color: #B45309; font-weight: 600" if v == "unsure" else "",
                                  subset=["confidence"]), hide_index=True, height=min(60 + 35 * len(overview), 560))

    user = {}      # pending = dropdowns that differ from the plan in THIS run; never kept in session_state
    by_id = {d.id: d for d in DECISIONS}
    if unsure:
        st.subheader("Worth a look", anchor=False)
        st.caption("Close calls: more than one option scored about the same. Read the reasons and change the choice "
                   "if you know something the data cannot tell.")
        for r in unsure:
            decision_editor(plan, by_id[r["id"]], user, expanded=True)
    st.subheader("All decisions", anchor=False)
    with st.expander("What do the evidence levels mean?", icon=":material/help:"):
        for k, v in EVIDENCE_LEVELS.items():
            st.write(f"**{k}**: {v}")
    for stage in STAGES:
        decisions = [d for d in DECISIONS if d.stage == stage.id and d.id in plan.selection["results"]
                     and not plan.selection["results"][d.id]["uncertain"]]
        if not decisions:
            continue
        st.markdown(f"##### {stage.name}")
        for d in decisions:
            decision_editor(plan, d, user)
    if user and st.button(f"Apply {len(user)} change(s) and rebuild the plan", type="primary", icon=":material/refresh:"):
        for k in [k for k in S.keys() if k.startswith("num_")]:     # a number sent with this click: the boxes, drawn
            dec = plan.decisions.get(k.removeprefix("num_"))           # below, have not written it into the plan yet
            if dec is not None and from_box(dec, S[k]) != dec.value:
                plan.set(dec.key, from_box(dec, S[k]), TYPED)
        new_plan = plan_preprocessing(diag, tf_list=S.get("tfs"), choices={**plan.choices, **user})
        for key, dec in plan.decisions.items():      # a value typed in a box stays typed
            if dec.reason == TYPED and key in new_plan.decisions:
                new_plan.set(key, dec.value, TYPED)
        S["plan"] = new_plan
        reset_downstream(*DOWNSTREAM)
        # A box keeps its value through the rebuild (on the server and in the browser) and would write the old plan's
        # value into the new one. Writing the new value before the boxes below are drawn makes both show it.
        for k in [k for k in S.keys() if k.startswith("num_")]:
            dec = new_plan.decisions.get(k.removeprefix("num_"))
            if dec is None:
                del S[k]
            else:
                S[k] = box_value(dec.key, dec.value)
        st.rerun()

    st.subheader("Thresholds and sizes", anchor=False)
    editable = ["min_genes_per_cell", "max_genes_per_cell", "max_counts_per_cell",
                "max_pct_mito_cell", "min_cells_per_gene", "subsample_cells", "cells_per_metacell"]
    shown = [k for k in editable if k in plan.decisions]
    for i in range(0, len(shown), 2):
        cols = st.columns(2)
        for col, key in zip(cols, shown[i:i + 2]):
            dec = plan.decisions[key]
            cur = dec.value
            as_int = key != "max_pct_mito_cell"      # every other threshold is a count
            new = col.number_input(dec.label, value=box_value(key, cur),
                                   step=1 if as_int else 0.5, key=f"num_{key}", help="0 = off" if box_off(dec) else None)
            new_v = from_box(dec, new)
            if new_v != cur:
                plan.set(key, new_v, TYPED)
                reset_downstream(*DOWNSTREAM)
            col.caption(logic.decision_caption(dec, plan))      # after set: an edit shows what the engine chose at once
    st.markdown("**Cellular contexts:** " + ", ".join(f"`{c['name']}`" for c in plan.contexts))

    with st.expander("Reproduce from a terminal", icon=":material/terminal:"):
        from grnbench.autoprep import DEFAULT_SEED
        st.caption("One command per cellular context, with the seed *Run cleaning* uses.")
        for cmd in logic.reproduce_commands(plan, S["data_path"], tf_list=S.get("tfs"), seed=DEFAULT_SEED):
            st.code(cmd, language="bash")
        import json
        st.download_button("Download plan.json", json.dumps(plan.to_dict(), indent=2, default=str),
                           file_name="plan.json", icon=":material/download:")
    footer(back=PROFILE, forward=(f"Continue: {CLEAN}", CLEAN, True), alt=("Compare pipelines first (optional)", COMPARE))


def warn_if_contexts_left_out(screen) -> None:
    """Say when the winning pipeline cleans fewer cell types than the recommended one (see logic.screening_scope)."""
    scope = logic.screening_scope(screen)
    if scope["left_out"]:
        st.warning(f"**{scope['best']}** keeps only {', '.join(scope['kept'])}; it leaves out "
                   f"{', '.join(scope['left_out'])}. The engine scores each candidate on the mean of its own contexts, "
                   "so a candidate that drops hard cell types can score higher without being better on the same "
                   "cells. Compare it with the others on the context they share before adopting it.")


# -- 4 compare pipelines ----------------------------------------------------------------- #
METRIC_NOTES = {
    "reference_auprc_ratio": ("Recovery of the reference", "How well the pipeline's network recovers the selection half of "
                              "the reference network (1 = random). Higher is better."),
    "edge_stability": ("Edge stability", "Share of the strongest links found again in an independent half of the cells. "
                       "Higher is better."),
    "depth_leakage": ("Depth leakage", "How strongly the most variable genes still follow sequencing depth. Lower is better."),
    "batch_leakage": ("Batch leakage", "How well the batches can still be told apart after cleaning (0 = fully mixed). "
                      "Lower is better."),
    "n_cells": ("Cells kept", "Cells that survive cleaning. More is better, unless the extra cells are poor."),
}


def page_compare() -> None:
    need_data()
    from grnbench.kb import DECISIONS
    from grnbench.screen import screen_candidates
    from grnbench.selector import select
    from grnbench.profile import Profile

    plan = S["plan"]
    selection = select(Profile(**plan.profile))
    rec = selection.chosen()
    differs = [f"{k} = {v}" for k, v in plan.choices.items() if v != rec.get(k)]
    differs += [f"{d.label} = {d.value}" for d in plan.decisions.values()
                if d.reason == TYPED]
    if differs:
        st.warning("The current plan differs from the engine's recommendation (" + "; ".join(differs) + "). "
                   "The candidates below start from the recommendation, not from the current plan, so none of "
                   "them has these differences, and 'Use this pipeline' replaces the current plan with the one "
                   "you pick. To keep the current plan, do not adopt a candidate.")
    auto = selection.candidates(max_candidates=4, always=("context_scope",))     # exactly what screen_candidates runs
    extra = S.setdefault("extra_candidates", [])
    with st.container(border=True, key="card-candidates"):
        st.markdown(style.card_title_html(1, "What will be compared", "Each candidate is cleaned and measured on your data."),
                    unsafe_allow_html=True)
        st.markdown("Automatic alternatives (the uncertain decisions, plus the cellular context, which is always "
                    "tried): " + (", ".join(f"`{c['name']}`" for c in auto[1:]) or "none - every decision was clear-cut") +
                    ". The fixed `baseline` pipeline is always screened too; without a reference network it is "
                    "the default choice.")
        with st.expander("Add your own alternative", icon=":material/add_circle:"):
            c1, c2, c3 = st.columns([2, 2, 1], vertical_alignment="bottom")
            dec = c1.selectbox("Decision", [d.id for d in DECISIONS])
            point = next(d for d in DECISIONS if d.id == dec)
            opts = [o.id for o in point.options if point.availability(o.id)[0]]
            opt = c2.selectbox("Option", opts)
            if c3.button("Add", width="stretch"):
                name = f"{dec}={opt}"
                if opt == rec.get(dec):
                    st.warning(f"`{name}` is already the recommended pipeline's choice: nothing to compare.")
                elif name in {c["name"] for c in auto + extra}:
                    st.warning(f"`{name}` is already in the comparison.")
                else:
                    extra.append({"name": name, "choices": {dec: opt}})
        if extra:
            st.markdown("Custom: " + ", ".join(f"`{e['name']}`" for e in extra))
            if st.button("Clear custom alternatives"):
                extra.clear()
                st.rerun()

    if st.button("Run and compare", type="primary", icon=":material/compare_arrows:"):
        with st.spinner("Running every candidate pipeline: about a minute for a small dataset, longer for a large one..."):
            S["screen"] = screen_candidates(S["diag"], S["data_path"], logic.fresh_folder(Path(S["workdir"]) / "screening"),
                                            tf_list=S.get("tfs"), reference_pairs=S.get("reference"),
                                            selection=selection, extra=list(extra))
    res = S.get("screen")
    if res is None:
        footer(back=PIPELINE, forward=(f"Skip: {CLEAN}", CLEAN, False))
        return
    import plotly.express as px
    rows = res["rows"]
    st.subheader("Results", anchor=False)
    st.success(f"Best: **{res['best']}** - rule: {res['rule']}", icon=":material/emoji_events:")
    warn_if_contexts_left_out(res)
    ranking = res["ranking"].round(3)
    sizes = {c: st.column_config.NumberColumn(c, format="%.0f", help="Averaged over the contexts of the candidate.")
             for c in ("n_cells", "n_genes", "n_tfs") if c in ranking}
    st.dataframe(ranking, column_config=sizes)
    metrics = [m for m in METRIC_NOTES if m in rows and rows[m].notna().any()]
    if metrics:
        tabs = st.tabs([METRIC_NOTES[m][0] for m in metrics])
        for tab, m in zip(tabs, metrics):
            with tab:
                st.caption(METRIC_NOTES[m][1])
                st.plotly_chart(px.bar(rows, x="candidate", y=m, color="context", barmode="group", height=330),
                                key=f"scr_{m}")
    st.subheader("Pros and cons against the recommended pipeline", anchor=False)
    pc = [{"candidate": n, "better": "; ".join(v["pros"]) or "-", "worse": "; ".join(v["cons"]) or "-"}
          for n, v in res["pros_cons"].items() if n != "recommended"]
    if pc:
        st.dataframe(pd.DataFrame(pc), hide_index=True)
    with st.container(border=True, key="card-adopt"):
        st.markdown(style.card_title_html(2, "Choose the pipeline for cleaning"), unsafe_allow_html=True)
        st.caption("Nothing changes until you press Use this pipeline. Continue keeps the pipeline you have now.")
        c1, c2 = st.columns([3, 1], vertical_alignment="bottom")
        adopt = c1.selectbox("Pipeline to use for cleaning", list(res["plans"]),
                             index=list(res["plans"]).index(res["best"]) if res["best"] in res["plans"] else 0)
        if c2.button("Use this pipeline", type="primary", width="stretch"):
            S["plan"] = res["plans"][adopt]
            S["user_choices"] = {}   # step 3's edits and pending picks belonged to the plan just replaced
            reset_downstream("records", "runs", "bench")
            st.success(f"Pipeline '{adopt}' adopted - go to 5 · Clean."
                       + (" It replaced the previous plan." if differs else ""))
    footer(back=PIPELINE, forward=(f"Continue: {CLEAN}", CLEAN, True))


# -- 5 clean --------------------------------------------------------------------------- #
def page_clean() -> None:
    need_data()
    from grnbench.autoprep import DEFAULT_SEED, execute_plan

    plan = S["plan"]
    with st.container(border=True, key="card-clean"):
        st.markdown(f"**Strategy:** `{plan.value('strategy')}` · **Contexts:** "
                    + ", ".join(f"`{c['name']}`" for c in plan.contexts))
        st.caption("A context is the set of cells one network is built from. Each one is cleaned separately.")
        if st.button("Run cleaning", type="primary", icon=":material/cleaning_services:"):
            with st.spinner("Cleaning..."):
                S["records"] = execute_plan(plan, S["data_path"], logic.fresh_folder(Path(S["workdir"]) / "clean"),
                                            tf_list=S.get("tfs"), seed=DEFAULT_SEED)    # the seed step 3 shows
            reset_downstream("runs", "bench", "bench_on")
    recs = S.get("records")
    if not recs:
        footer(back=PIPELINE)
        return
    tbl = logic.retention_table(recs)
    st.subheader("What the cleaning did", anchor=False)
    for r in recs:
        rows = tbl[(tbl["context"] == r["context"]) & (tbl["step"] != "FAILED")]
        if r["error"] or rows.empty:
            st.error(f"{r['context']}: {r['error']}")
            continue
        c = st.columns([1.4, 1, 1, 1])
        c[0].markdown(f"**{r['context']}**")
        c[1].metric("Cells kept", f"{rows['cells'].iloc[-1]:,.0f}", f"of {rows['cells'].iloc[0]:,.0f}", delta_color="off",
                    delta_arrow="off", border=True)
        c[2].metric("Genes kept", f"{rows['genes'].iloc[-1]:,.0f}", f"of {rows['genes'].iloc[0]:,.0f}", delta_color="off",
                    delta_arrow="off", border=True)
        c[3].metric("Share of cells kept", f"{rows['cells'].iloc[-1] / max(rows['cells'].iloc[0], 1):.0%}", border=True)
    with st.expander("Step by step: what each step removed", icon=":material/format_list_numbered:"):
        st.dataframe(tbl, hide_index=True)
    with st.expander("Exact commands", icon=":material/terminal:"):
        for r in recs:
            st.code(logic.shell_command(r["cli"]), language="bash")
    footer(back=PIPELINE, forward=(f"Continue: {MODELS}", MODELS, True), alt=("Look at the cleaned data", EXPLORE_DATA))


# -- 6 models --------------------------------------------------------------------------- #
def page_models() -> None:
    need_cleaned()

    plan = S["plan"]
    with st.container(border=True, key="card-models"):
        st.markdown(style.card_title_html(None, "Choose the models", "Every model gets exactly the same cleaned data."),
                    unsafe_allow_html=True)
        recommended = [m["name"] for m in plan.methods if m.get("recommended")]
        methods = model_picker("models", recommended, all_regulators=S.get("tfs") is None)
    with st.expander("Add your own model", icon=":material/extension:"):
        st.caption("A .py file with a class exposing `infer(adata, tf_list, **params)` (or a function "
                   "`f(adata, tf_list)`) returning columns source, target, weight. It gets exactly the "
                   "same cleaned data as the built-in models.")
        c1, c2, c3 = st.columns([2, 1, 1], vertical_alignment="bottom")
        up = c1.file_uploader("Model file", type=["py"])
        attr = c2.text_input("Class or function name")
        if c3.button("Register", disabled=not (up and attr)):
            path = save_upload(up, Path(S["workdir"]) / "user_models")
            try:
                reg = logic.register_user_model(path, attr)
                S.setdefault("user_methods", set()).update(r["name"] for r in reg)
                st.success("Registered: " + ", ".join(r["name"] for r in reg))
            except Exception as exc:  # noqa: BLE001
                st.error(f"Could not register the model: {exc}")
    if st.button("Run models", type="primary", disabled=not methods, icon=":material/play_arrow:"):
        with st.status("Running the models...", expanded=True) as status:
            bar = st.progress(0.0)
            S["runs"] = logic.run_models(S["records"], methods, logic.fresh_folder(Path(S["workdir"]) / "models"),
                                         tf_list=S.get("tfs"), strategy=plan.value("strategy"),
                                         progress=lambda f, msg: bar.progress(f, text=msg))
            status.update(label="Models finished", state="complete", expanded=False)
        reset_downstream("bench")
    runs = S.get("runs")
    if runs is not None:
        st.subheader("Results", anchor=False)
        for _, row in (runs[runs["edges"].isna()].iterrows() if "edges" in runs else []):
            st.warning(f"{row['method']} did not finish on {row['context']}: {row['error']}")
        st.dataframe(logic.runs_view(runs), hide_index=True, column_config={
            "n_edges": st.column_config.NumberColumn("links", help="Regulator -> gene links the model returned.", format="%d"),
            "seconds": st.column_config.NumberColumn("seconds", format="%.1f")})
    footer(back=CLEAN, forward=(f"Continue: {BENCHMARK}", BENCHMARK, runs is not None),
           alt=("Explore the networks", EXPLORE_RESULTS) if runs is not None else None)


# -- 7 benchmark ------------------------------------------------------------------------- #
def page_benchmark() -> None:
    need_models()
    with st.container(border=True, key="card-reference"):
        st.markdown(style.card_title_html(None, "Reference network", "Known regulator -> target links the models are scored against."),
                    unsafe_allow_html=True)
        up = st.file_uploader("Or score against another reference network", type=["tsv", "csv", "txt", "gz"])
        ref, which, warn = logic.reference_for_scoring(
            S.get("screen"), S.get("reference"),
            read_upload("a reference network (a source and a target column)", logic.read_pairs, save_upload(up, Path(S["workdir"]) / "inputs")) if up else None, up.name if up else "")
        if warn:
            st.warning(warn)
        if not ref:
            st.info("No reference network yet: load one in step 1 (Data) or upload it above.", icon=":material/info:")
            footer(back=MODELS, forward=("Look at the networks instead", EXPLORE_RESULTS, False))
            raise PageDone
        st.markdown(f"**Reference:** {which} · {len(ref):,} links")
    if st.button("Score models", type="primary", icon=":material/leaderboard:"):
        S.pop("bench", None)
        with guarded("score the models"):
            S["bench"] = logic.benchmark(S["runs"], ref, tf_list=S.get("tfs"))
            S["bench_on"] = f"{which} · {len(ref):,} links"
    bench = S.get("bench")
    if bench is None or bench.empty:
        footer(back=MODELS)
        return
    user = S.get("user_methods", set())
    bench = bench.assign(kind=["your model" if m in user else "built-in" for m in bench["method"]])
    summ = logic.run_summary(S["records"], S["runs"], bench)
    if summ["top"]:
        st.metric("Highest AUPRC ratio here", f"{summ['top']['method']}", f"{summ['top']['ratio']:.1f}x random", delta_color="off",
                  border=True, help="The highest AUPRC ratio on this dataset. One dataset is not a ranking of the models.")
    metric = st.radio("Show", ["AUPRC ratio", "EPR", "AUROC"], horizontal=True, key="bench_metric")
    st.plotly_chart(logic.score_figure(bench, metric), key=f"bench_{metric}")
    show_score_table(bench)
    c1, c2, _ = st.columns([1, 1.2, 2])
    c1.download_button("Download scores (CSV)", bench.to_csv(index=False), file_name="benchmark.csv",
                       icon=":material/download:", width="stretch")
    c2.download_button("Download the whole run (zip)", logic.zip_folder(S["workdir"]),
                       file_name=f"{Path(S['workdir']).name}.zip", icon=":material/folder_zip:", width="stretch")
    footer(back=MODELS, forward=("Explore the networks", EXPLORE_RESULTS, True))


# -- shared result widgets ------------------------------------------------------------------- #
SCORE_HELP = {
    "AUPRC ratio": "Area under the precision-recall curve divided by what a random ranking gets. 1 = random.",
    "EPR": "Early precision ratio: share of true links among the top k predictions (k = number of reference links "
           "found), divided by the random value. 1 = random.",
    "AUROC": "Area under the ROC curve. 0.5 = random, 1 = perfect.",
    "AUPRC": "Area under the precision-recall curve (raw, before dividing by the random level).",
}


def show_score_table(bench) -> None:
    """The numbers behind the chart, with a tooltip on every score column and how to read them."""
    st.dataframe(bench.sort_values("AUPRC ratio", ascending=False), hide_index=True,
                 column_order=[c for c in ["method", "context", "AUPRC ratio", "EPR", "AUROC", "AUPRC",
                                           "reference edges in universe"] if c in bench],
                 column_config={col: st.column_config.NumberColumn(col, help=SCORE_HELP[col],
                                                                   format="%.3f" if col in ("AUPRC", "AUROC") else "%.1f")
                                for col in SCORE_HELP if col in bench})
    with st.expander("How to read these scores", icon=":material/help:"):
        for name, text in SCORE_HELP.items():
            st.markdown(f"**{name}**: {text}")
        st.markdown("Contexts from one dataset are not independent replicates: these numbers describe this dataset. "
                    "They are not a ranking of the models; a claim across datasets needs `grnbench compare` over "
                    "independent datasets.")


def show_scores(bench, key: str) -> None:
    """Models ordered by score with the random level marked, then the numbers."""
    st.plotly_chart(logic.score_figure(bench), key=f"{key}_scores")
    st.caption("AUPRC ratio: how many times better than a random ranking a model puts true links first (1 = random). "
               "Higher is better." + (f" Scored on: {S['bench_on']}." if S.get("bench_on") else ""))
    show_score_table(bench)


def show_run_results() -> None:
    """What a quick run produced: headline numbers, scores, cleaning and model tables, where to go next."""
    summ = logic.run_summary(S["records"], S["runs"], S.get("bench"))
    st.subheader("Results", anchor=False)
    failed = [r for r in S["records"] if r["error"]]
    for r in failed:
        st.error(f"Cleaning failed for {r['context']}: {r['error']}")
    st.success(f"Analysis complete: {summ['models_ok']} of {summ['models_total']} model run(s) finished on "
               f"{summ['cells_used']:,} cells.", icon=":material/task_alt:")
    if S.get("screen"):
        st.info(f"Cleaning pipeline chosen by the engine: **{S['screen']['best']}** (contexts: "
                + ", ".join(c["name"] for c in S["plan"].contexts) + ").", icon=":material/tune:")
        warn_if_contexts_left_out(S["screen"])
    c = st.columns(4)
    c[0].metric("Cells used", f"{summ['cells_used']:,}", f"of {summ['cells_loaded']:,} loaded" if summ["cells_loaded"] else None,
                delta_color="off", delta_arrow="off", border=True,
                help="Cells left after quality control: empty droplets, damaged cells and doublets are removed.")
    c[1].metric("Genes", f"{summ['genes_used']:,}", border=True)
    c[2].metric("Models run", f"{summ['models_ok']} of {summ['models_total']}", border=True)
    top = summ["top"]
    c[3].metric("Highest AUPRC ratio here", f"{top['method']}" if top else "not scored", border=True,
                delta=f"{top['ratio']:.1f}x random" if top else None, delta_color="off",
                help="The highest AUPRC ratio on this dataset. One dataset is not a ranking of the models."
                if top else "Give a reference network to score the models.")
    for _, row in S["runs"][S["runs"]["edges"].isna()].iterrows() if "edges" in S["runs"] else []:
        st.warning(f"{row['method']} did not finish on {row['context']}: {row['error']}")
    bench = S.get("bench")
    t_scores, t_clean, t_models = st.tabs(["Scores", "Cleaning", "Models"])
    with t_scores:
        if S.get("reference_warning"):
            st.warning(S["reference_warning"])
        if bench is not None and not bench.empty:
            show_scores(bench, "quick")
        elif not S.get("reference_warning"):
            st.info("No reference network, so the models were run but not scored. Look at their links in "
                    "Explore results, or add a reference in step 7.", icon=":material/info:")
    with t_clean:
        st.dataframe(logic.retention_table(S["records"]), hide_index=True)
    with t_models:
        st.dataframe(logic.runs_view(S["runs"]), hide_index=True)
    st.markdown("**What next?**")
    n1, n2, n3, _ = st.columns([1.3, 1.3, 1.3, 1])
    n1.button("Explore the networks", on_click=goto, args=(EXPLORE_RESULTS,), type="primary",
              icon=":material/share:", width="stretch", key="next-results")
    n2.button("Explore your data", on_click=goto, args=(EXPLORE_DATA,), icon=":material/scatter_plot:",
              width="stretch", key="next-data")
    n3.button("Step-by-step workflow", on_click=goto, args=(PROFILE,), icon=":material/route:",
              width="stretch", key="next-flow", help="Review the profile, the pipeline and every setting.")
    with st.expander("Where are my files?", icon=":material/folder_open:"):
        st.code(str(S["workdir"]), language=None)
        st.download_button("Download the whole run (zip)", logic.zip_folder(S["workdir"]),
                           file_name=f"{Path(S['workdir']).name}.zip", icon=":material/download:")


# -- 0 quick run ------------------------------------------------------------------------- #
def page_quick() -> None:
    done = bool(S.get("quick_done")) and S.get("runs") is not None
    if done:
        show_run_results()
        st.markdown('<div class="gb-foot"></div>', unsafe_allow_html=True)
    # after a run the form folds away: the results are what the person came for
    form = (st.expander("Change the data or the models, then run again", icon=":material/replay:")
            if done else st.container())
    with form:
        quick_form()


def quick_form() -> None:
    with st.container(border=True, key="card-data"):
        st.markdown(style.card_title_html(1, "Your data", "Try the demo to see how it works, or use your own counts."),
                    unsafe_allow_html=True)
        data_path, tf_path, gold_path = pick_dataset("quick")
    with st.container(border=True, key="card-regulators"):
        st.markdown(style.card_title_html(2, "Regulators and reference network",
                                          "Both are optional. Each one makes the result better."),
                    unsafe_allow_html=True)
        c1, c2 = st.columns(2)
        tf_up = c1.file_uploader("TF list", type=["txt", "tsv", "csv"], key="quick_tf",
                                 help="One gene symbol per line (.txt), or a .csv/.tsv with a symbol column.")
        c1.caption("The genes allowed to act as regulators. Without it every gene is a candidate regulator, "
                   "which is much slower.")
        ref_up = c2.file_uploader("Reference network", type=["tsv", "csv", "txt", "gz"], key="quick_ref",
                                  help="Known regulator -> target links, in two columns named source and target.")
        c2.caption("Known regulator -> target links. They let grnbench score each model.")
        if tf_path is not None:
            st.caption("The demo comes with both; uploading a file here replaces the demo's.")
    no_tfs = tf_path is None and tf_up is None
    with st.container(border=True, key="card-models"):
        st.markdown(style.card_title_html(3, "Models", "Start from a preset, or open the list to choose each model."),
                    unsafe_allow_html=True)
        if no_tfs and data_path is not None:
            st.warning("No TF list: every gene is a candidate regulator, so each model scores about genes x genes "
                       "pairs (about 1,000,000 edges at 1,000 genes) and lasso, tigress and community take minutes "
                       "to hours, not seconds. Upload a TF list above, or leave those models unticked.")
        methods = model_picker("quick", all_regulators=no_tfs)
    with st.expander("Advanced options", icon=":material/settings_suggest:"):
        screen = st.checkbox("Let the engine compare candidate cleaning pipelines first (slower)", key="quick_screen",
                             help="It only changes the choice when you give a reference network.")
    ready = data_path is not None and bool(methods)
    st.markdown('<div class="gb-foot"></div>', unsafe_allow_html=True)
    left, right = st.columns([3.2, 1])
    if not data_path:
        left.markdown("**Choose a dataset to start.**")
        left.caption("Create the demo dataset, upload a file or give a path in card 1.")
    elif not methods:
        left.markdown("**Choose at least one model.**")
        left.caption("Pick a preset or tick models in card 3.")
    else:
        dataset = "the demo dataset" if tf_path is not None else Path(str(data_path)).name
        left.markdown(f"**Ready to run:** {dataset} · "
                      f"{'your TF list' if (tf_up or tf_path) else 'every gene as a regulator'} · "
                      f"{'a reference network' if (ref_up or gold_path) else 'no reference, so models will not be scored'}"
                      f" · {logic.plural(len(methods), 'model')}")
    if right.button("Run", type="primary", disabled=not ready, key="quick_run", icon=":material/play_arrow:",
                    width="stretch"):
        from grnbench.tfs import load_tf_list
        inputs = Path(base) / "uploads"
        if tf_up:
            tf_path = save_upload(tf_up, inputs)
        if ref_up:
            gold_path = save_upload(ref_up, inputs)
        tfs = read_upload("a TF list (one gene symbol per line)", load_tf_list, tf_path) if tf_path else None
        reference = read_upload("a reference network (a source and a target column)", logic.read_pairs, gold_path) if gold_path else None
        with st.status("Running the analysis...", expanded=True) as status:
            bar = st.progress(0.0)
            try:
                out = logic.quick_run(data_path, base, methods, tf_list=tfs, reference=reference, screen=screen,
                                      orientation=S.get("quick_orientation", "auto"),
                                      progress=lambda f, msg: bar.progress(f, text=msg))
            except Exception as exc:  # noqa: BLE001 - say what stopped it, not a stack trace
                status.update(label="The run stopped", state="error", expanded=True)
                st.error(f"The run stopped: {type(exc).__name__}: {exc}")
                raise PageDone from None
            status.update(label="Analysis complete", state="complete", expanded=False)
        reset_downstream(*DOWNSTREAM)
        S.update(workdir=out["workdir"], data_path=str(data_path), tfs=tfs, reference=reference,
                 diag=out["diag"], plan=out["plan"], user_choices={}, records=out["records"], runs=out["runs"],
                 reference_warning=out["reference_warning"])
        if out["bench"] is not None and not out["bench"].empty:
            S["bench"] = out["bench"]
            S["bench_on"] = out["scored_on"]
        if out["screen"] is not None:
            S["screen"] = out["screen"]
        S["quick_done"] = True
        st.rerun()          # draw the results above the (now folded) form


# -- 8 explore data ---------------------------------------------------------------------- #
def stage_data(stage: str):
    return cached("stage", (str(S.get("workdir")), stage), lambda: logic.load_stage(S["diag"], S.get("records"), stage))


def embedding_cached(stage: str, method: str, max_cells: int):
    return cached("embed", (str(S.get("workdir")), stage, method, max_cells),
                  lambda: logic.embedding_frame(stage_data(stage), method, max_cells=max_cells),
                  what=f"Computing the {method} map" + (" (UMAP takes 20-30 s the first time after the app starts)"
                                                       if method == "UMAP" else "") + "...")


def page_explore_data() -> None:
    need_data()
    st.caption("Exploration of this dataset, not a result: nothing here is a registered experiment or a claim.")
    stages = logic.stage_choices(S.get("records"))
    stage = st.radio("Look at", stages, horizontal=True, key="explore_stage", index=min(1, len(stages) - 1))
    adata = stage_data(stage)
    c = st.columns(4)
    c[0].metric("Cells", f"{adata.n_obs:,}", border=True)
    c[1].metric("Genes", f"{adata.n_vars:,}", border=True)
    detection, _ = logic.gene_detection(adata)
    c[2].metric("Zeros", "n/a" if np.isnan(detection).all() else f"{1 - np.nanmean(detection):.0%}", border=True,
                help="share of count entries that are zero; n/a when only a centred or scaled matrix is left")
    tfs = S.get("tfs")
    if tfs:
        c[3].metric("TFs present", f"{len({str(g).upper() for g in adata.var_names} & {t.upper() for t in tfs})}/{len(tfs)}",
                    border=True)
    if len(stages) == 1:
        st.info("This is the data as loaded. Clean it to compare before and after.", icon=":material/info:")
        st.button("Clean the data", on_click=goto, args=(CLEAN,), key="explore-clean", icon=":material/cleaning_services:")

    tab_map, tab_genes, tab_obs = st.tabs(["Map of the cells", "Genes and TFs", "Cell annotations"])
    with tab_map:
        st.caption("Each dot is a cell; cells with similar expression sit close together. Colour them by any "
                   "annotation to see whether the groups stay apart.")
        c1, c2, c3 = st.columns([1, 1, 2])
        method = c1.radio("Map", ["PCA", "UMAP"], horizontal=True, key="explore_method",
                          help="PCA is instant; UMAP takes longer and needs umap-learn")
        max_cells = c2.select_slider("Cells drawn", [1000, 2000, 3000, 5000], value=3000, key="explore_n")
        try:
            frame = embedding_cached(stage, method, max_cells)
        except Exception as exc:  # noqa: BLE001
            st.error(f"Could not draw the map: {exc}")
            frame = None
        if frame is not None:
            options = logic.color_choices(frame)
            color = c3.selectbox("Colour by", options or [None], key="explore_color")
            if frame.attrs["n_cells_total"] > len(frame):
                st.caption(f"{len(frame):,} of {frame.attrs['n_cells_total']:,} cells drawn (random sample).")
            if len(stages) > 1 and stage != stages[0]:
                try:
                    raw = embedding_cached(stages[0], method, max_cells)
                except Exception as exc:  # noqa: BLE001
                    st.warning(f"Could not draw the raw map: {exc}")
                    raw = None
                left, right = st.columns(2)
                if raw is not None:
                    left.plotly_chart(logic.embedding_figure(raw, color if color in raw else None,
                                                             title=f"Before: {stages[0]}"),
                                      key="explore_map_raw")
                right.plotly_chart(logic.embedding_figure(frame, color, title=f"After: {stage}"),
                                   key="explore_map")
                st.caption("Two separate maps: each is drawn from its own random sample of cells and has its own "
                           "axes, so a dot on one is not the same cell or place on the other. Compare their shapes: "
                           "damaged cells and empty droplets that stretch the raw map should be gone, and real groups "
                           "should stay apart.")
            else:
                st.plotly_chart(logic.embedding_figure(frame, color, title=f"{stage}: {method}"),
                                key="explore_map")
    with tab_genes:
        tf_key = None if tfs is None else tuple(sorted({str(t).upper() for t in tfs}))
        summary = cached("summary", (str(S.get("workdir")), stage, tf_key),
                         lambda: logic.gene_summary(adata, tfs), what="Summarising the genes...")
        source = summary.attrs.get("detection_source", "")
        if summary["detection"].isna().all():
            st.warning(f"Detection cannot be shown for this stage: {source}.")
        c1, c2 = st.columns([2, 1])
        query = c1.text_input("Find a gene", key="explore_find")
        only_tf = c2.checkbox("Only TFs", value=False, disabled="TF" not in summary, key="explore_only_tf")
        view = summary
        if query:
            view = view[view["gene"].str.upper().str.contains(query.strip().upper(), regex=False)]
        if only_tf and "TF" in view:
            view = view[view["TF"]]
        st.dataframe(view.head(500).round(4), hide_index=True)
        st.caption(f"detection = fraction of cells with a count above zero (read from the {source}); mean and "
                   "variance are of the normalised expression. Showing up to 500 rows.")
        default_genes = (view.loc[view["TF"], "gene"].head(3).tolist() if "TF" in view and view["TF"].any()
                         else view["gene"].head(3).tolist())
        genes = st.multiselect("Plot the expression of", summary["gene"].tolist(), default=default_genes,
                               max_selections=6, key=f"explore_genes_{stage}")
        group_opts = [None] + [col for col in adata.obs.columns
                               if 2 <= adata.obs[col].astype(str).nunique() <= 12]
        group = st.selectbox("Split by", group_opts, key="explore_group")
        if genes:
            st.plotly_chart(logic.expression_figure(logic.expression_frame(adata, genes, group), group),
                            key="explore_expr")
    with tab_obs:
        st.dataframe(adata.obs.head(300))
        st.download_button("Download the cell annotations (CSV)", adata.obs.to_csv(), file_name="cell_annotations.csv")


# -- 9 explore results ------------------------------------------------------------------- #
def reference_key(reference) -> tuple:
    return (len(reference), hash(frozenset(reference))) if reference else (0, 0)


def page_explore_results() -> None:
    need_models()
    from grnbench.infer import load_edges

    st.info("Exploration of this one dataset, not a result and not a ranking of the models: what agrees "
            "or is recovered here describes these runs. A claim across datasets needs `grnbench compare` over "
            "independent datasets.")
    runs = S["runs"]
    ok = runs[runs["edges"].notna()] if "edges" in runs else runs.iloc[0:0]
    if ok.empty:
        st.warning("No model produced edges. See the errors in 6 · Models.")
        raise PageDone
    contexts = sorted(ok["context"].unique())
    with st.container(border=True, key="card-view"):
        c1, c2 = st.columns([1, 3]) if len(contexts) > 1 else (None, st.container())
        context = c1.selectbox("Cellular context", contexts, key="res_context") if c1 is not None else contexts[0]
        sel = ok[ok["context"] == context]
        methods = c2.multiselect("Models", sel["method"].tolist(), default=sel["method"].tolist(), key="res_methods")
    if not methods:
        st.info("Pick at least one model.")
        raise PageDone
    paths = sel.set_index("method")["edges"]
    edges = {}
    for m in methods:
        try:
            edges[m] = cached("edges", paths[m], lambda p=paths[m]: load_edges(p))
        except Exception as exc:  # noqa: BLE001
            st.warning(f"{m}: its edges could not be read ({type(exc).__name__}: {exc}).")
    methods = [m for m in methods if m in edges]
    if not methods:
        raise PageDone
    with st.expander("Reference network to compare with (optional)", icon=":material/rule:"):
        up = st.file_uploader("Reference network", type=["tsv", "csv", "txt", "gz"], key="res_ref")
        reference, which, warn = logic.reference_for_scoring(
            S.get("screen"), S.get("reference"),
            read_upload("a reference network (a source and a target column)", logic.read_pairs, save_upload(up, Path(S["workdir"]) / "inputs")) if up else None, up.name if up else "")
        if warn:
            st.warning(warn)
        if reference:
            st.caption(f"Reference: {which} - {len(reference):,} links.")
        else:
            st.caption("No reference network: links are shown without a right/wrong flag.")
    import anndata as ad
    processed = sel["processed"].iloc[0]
    genes = cached("genes", processed, lambda: pd.Index(ad.read_h5ad(processed, backed="r").var_names).astype(str))
    key = (tuple(str(paths[m]) for m in methods), reference_key(reference))

    tab_net, tab_edges, tab_agree, tab_pr, tab_tf = st.tabs(
        ["Network", "Top links", "Do the models agree?", "Precision-recall", "Per TF"])
    with tab_edges, guarded("list the edges"):
        c1, c2, c3, c4 = st.columns([2, 1, 1, 1])
        m = c1.selectbox("Model", methods, key="res_edges_model")
        tf = c2.text_input("TF", key="res_edges_tf")
        gene = c3.text_input("Gene contains", key="res_edges_gene")
        top = c4.number_input("Rows", 10, 5000, 200, step=50, key="res_edges_top")
        tbl = logic.edges_table(edges[m], reference, tf=tf or None, gene=gene or None, top=int(top))
        st.dataframe(tbl, hide_index=True)
        st.caption("Edges with the same score share a rank.")
        st.download_button("Download this table (CSV)", tbl.to_csv(index=False), file_name=f"{m}_edges.csv")
    with tab_net, guarded("draw the network"):
        c1, c2, c3 = st.columns([2, 1, 1])
        m = c1.selectbox("Model", methods, key="res_net_model")
        tf = c2.text_input("Only this TF (optional)", key="res_net_tf")
        n = c3.slider("Edges drawn", 10, 150, 50, key="res_net_n")
        st.plotly_chart(logic.network_figure(edges[m], reference, top_n=n, tf=tf or None),
                        key="res_net")
        st.caption("Squares are TFs, circles are genes; arrows point from the regulator to its target, and green "
                   "edges are in the reference network. "
                   "Only the best-ranked edges are drawn. A name is printed only where it fits: hover any node for its name.")
    with tab_agree, guarded("compare the models"):
        k = st.slider("Compare the top", 50, 2000, 200, step=50, key="res_agree_k")
        if len(methods) < 2:
            st.info("Pick two or more models.")
        else:
            mat = cached("overlap", (key[0], k), lambda: logic.overlap_matrix(edges, k))
            for m, why in mat.attrs["skipped"].items():
                st.warning(f"{m} left out of the comparison: {why}.")
            if len(mat) < 2:
                st.info("Fewer than two models have scores that can be compared.")
            else:
                st.plotly_chart(logic.overlap_figure(mat, k), key="res_agree")
                sizes = mat.attrs["sizes"]
                st.caption("1.00 = the same edges, 0 = none in common. Models built on the same signal agree "
                           "more; agreement is not accuracy. Edges tied at the cut are all kept, so a set can be "
                           "larger than the number asked: "
                           + ", ".join(f"{m} {n:,}" for m, n in sizes.items()) + ".")
    with tab_pr, guarded("draw the precision-recall curves"):
        if not reference:
            st.info("Needs a reference network (1 · Data, or upload one above).")
        else:
            tfs = S.get("tfs")
            fig = cached("pr", (*key, hash(tuple(sorted(map(str, tfs)))) if tfs else None),
                         lambda: logic.pr_figure(edges, reference, genes, tfs),
                         what="Scoring the models against the reference (about 30 s for 4 models x 1M edges)...")
            st.plotly_chart(fig, key="res_pr")
            st.caption("Precision against recall over the full TF x gene universe; the dotted line is what a "
                       "random ranking would give.")
    with tab_tf, guarded("compare the models TF by TF"):
        if not reference:
            st.info("Needs a reference network (1 · Data, or upload one above).")
        else:
            k = st.slider("Top predictions per TF", 5, 100, 20, key="res_tf_k")
            tbl = cached("tfrec", (*key, k), lambda: logic.tf_recovery_table(edges, reference, genes, k_per_tf=k),
                         what="Counting the recovered edges per TF...")
            if tbl.empty:
                st.info("None of the reference's regulators has stored edges for these models in this data.")
            else:
                st.plotly_chart(logic.tf_recovery_figure(tbl, k), key="res_tf")
                st.caption(f"Only reference regulators and targets that are genes of this data count. A random "
                           f"ranking would recover about {k / max(len(genes) - 1, 1):.1%}. Blank = the model stored "
                           "no edge for that TF. A TF has few reference targets, so one TF is an anecdote, not a claim.")
                c1, c2 = st.columns(2)
                tf = c1.selectbox("Look at one TF", tbl.index.tolist(), key="res_tf_pick")
                m = c2.selectbox("Model", methods, key="res_tf_model")
                top, missed, n_out = logic.tf_detail(edges[m], reference, tf, k, genes)
                st.dataframe(top, hide_index=True)
                st.write(f"**Reference targets of {tf} not in its top {k}** ({len(missed)}): "
                         + (", ".join(missed[:60]) + (" ..." if len(missed) > 60 else "") if missed else "none"))
                if n_out:
                    st.caption(f"{n_out} more reference target(s) of {tf} are not genes of this data and are left out.")


PAGES = {QUICK: page_quick, DATA: page_data, PROFILE: page_profile, PIPELINE: page_pipeline,
         COMPARE: page_compare, CLEAN: page_clean, MODELS: page_models, BENCHMARK: page_benchmark,
         EXPLORE_DATA: page_explore_data, EXPLORE_RESULTS: page_explore_results}
if S["step"] not in PAGES:
    S["step"] = QUICK
st.markdown(style.header_html(S["step"]), unsafe_allow_html=True)
try:
    if base_error:
        st.error(base_error + " Change it under Settings in the sidebar.", icon=":material/folder_off:")
    else:
        PAGES[S["step"]]()
except PageDone:
    pass
except Exception:        # a crash in one screen must not take the menu with it
    render_nav()
    raise
render_nav()
