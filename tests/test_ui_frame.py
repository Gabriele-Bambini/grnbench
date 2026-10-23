"""Interface frame: the list of screens, the product theme, the sidebar menu, headers, empty states, footers and the
plain-language summaries the screens are built from. Helpers are called directly; pages are rendered with
Streamlit's AppTest (no browser)."""
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from grnbench import cli_main
from grnbench.autoprep import execute_plan
from grnbench.demo import write_demo
from grnbench.ui import logic, steps, style, theme

ROOT = Path(__file__).resolve().parents[1]
APP = ROOT / "src" / "grnbench" / "ui" / "app.py"


def _at(timeout=180):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    return AppTest.from_file(str(APP), default_timeout=timeout)


def _goto(at, step):
    at.session_state["step"] = step
    at.run()
    assert not at.exception, (step, [e.value for e in at.exception])


def _text(elements):
    return " ".join(str(e.value) for e in elements)


def _expander_labels(block):
    """AppTest reports an expander that has an icon as a ``status`` block, one without as an ``expander``."""
    return [e.label for e in list(block.expander) + list(block.status)]


@pytest.fixture(scope="module")
def demo(tmp_path_factory):
    d = tmp_path_factory.mktemp("uiframe")
    paths = write_demo(d / "demo", n_cells=900, n_genes=700, n_tfs=12, targets_per_tf=10)
    tfs = paths["tfs"].read_text().split()
    res = logic.analyze(paths["counts"], tf_list=tfs, choices={"context_scope": "pooled"})
    return d, paths, tfs, res


@pytest.fixture(scope="module")
def finished(demo):
    """A cleaned dataset with two models run and scored: what the screens after step 5 are drawn from."""
    d, paths, tfs, res = demo
    records = execute_plan(res["plan"], paths["counts"], d / "clean", tf_list=tfs)
    runs = logic.run_models(records, ["pearson_signed", "spearman_signed"], d / "models", tf_list=tfs)
    gold = pd.read_csv(paths["gold"], sep="\t")
    bench = logic.benchmark(runs, set(zip(gold.source, gold.target)), tf_list=tfs)
    return records, runs, bench


def _with_state(demo, **extra):
    d, paths, tfs, res = demo
    at = _at()
    state = dict(workdir=d, data_path=str(paths["counts"]), tfs=tfs, reference=None, diag=res["diag"],
                 plan=res["plan"].__class__.from_dict(res["plan"].to_dict()))
    state.update(extra)
    for k, v in state.items():
        at.session_state[k] = v
    at.run()
    assert not at.exception, [e.value for e in at.exception]
    return at


# --------------------------------------------------------------------------- #
# the screens
# --------------------------------------------------------------------------- #

def test_the_screens_are_one_ordered_list_with_a_headline_each():
    assert steps.LABELS[0] == steps.QUICK and len(set(steps.LABELS)) == len(steps.LABELS)
    order = [steps.GROUPS.index(s.group) for s in steps.STEPS]
    assert order == sorted(order), "screens are listed group by group"
    for s in steps.STEPS:
        assert s.title and s.lede and s.icon and s.group in steps.GROUPS, s.key
        assert s.lede.endswith("."), f"{s.key}: the lede is a sentence"
    assert steps.WORKFLOW_KEYS == [steps.DATA, steps.PROFILE, steps.PIPELINE, steps.COMPARE, steps.CLEAN,
                                   steps.MODELS, steps.BENCHMARK]
    assert steps.position(steps.DATA) == (1, 7) and steps.position(steps.BENCHMARK) == (7, 7)
    assert steps.position(steps.QUICK) is None and steps.position(steps.EXPLORE_RESULTS) is None
    assert [s.key for s in steps.STEPS if s.optional] == [steps.COMPARE]
    assert [s.key for s in steps.STEPS if s.group == steps.EXPLORE] == [steps.EXPLORE_DATA, steps.EXPLORE_RESULTS]


def test_every_screen_header_names_the_screen_and_its_place_in_the_workflow():
    for s in steps.STEPS:
        html = style.header_html(s.key)
        assert s.title in html and s.lede in html and "<style>" in html, s.key
        pos = steps.position(s.key)
        assert (f"Step {pos[0]} of {pos[1]} · {s.short}" in html) if pos else (s.group in html), s.key
        assert ('class="gb-opt"' in html) == s.optional, s.key
    assert "__ACCENT__" not in style.stylesheet() and theme.ACCENT in style.stylesheet()
    assert "__ACCENT__" not in style.brand_html() and "grnbench" in style.brand_html()


def test_card_titles_number_their_badge_and_indent_their_subtitle():
    html = style.card_title_html(2, "Models", "Pick some")
    assert 'class="gb-num">2<' in html and "Models" in html and "Pick some" in html and "margin-left:2.15rem" in html
    plain = style.card_title_html(None, "Models")
    assert "gb-num" not in plain and "gb-card-sub" not in plain
    assert "margin-left:0" in style.card_title_html(None, "Models", "Sub")


def test_the_stylesheet_hooks_match_the_keys_the_app_gives_its_widgets():
    """The CSS selects widgets by the ``st-key-*`` class Streamlit derives from their key: renaming a key in the
    app without the stylesheet would silently lose the look, so the two are checked against each other."""
    app = APP.read_text(encoding="utf-8")
    assert 'key=f"nav-{state}-{i}"' in app
    for state in ("cur", "done", "todo"):
        assert f"st-key-nav-{state}-" in style.CSS
    for key in ("explore_stage", "explore_method", "bench_metric", "empty-state"):
        assert f"st-key-{key}" in style.CSS and f'key="{key}"' in app, key
    for suffix in ("_source", "_preset"):
        assert f'[class*="{suffix}"]' in style.CSS and f'{suffix}"' in app, suffix


# --------------------------------------------------------------------------- #
# the product theme, handed to Streamlit by `grnbench ui`
# --------------------------------------------------------------------------- #

def _flat(flags):
    return dict(zip(flags[0::2], flags[1::2]))


def test_theme_flags_carry_the_accent_the_chart_colours_and_a_quiet_toolbar():
    flat = _flat(theme.streamlit_flags())
    assert re.fullmatch(r"#[0-9A-Fa-f]{6}", theme.ACCENT) and re.fullmatch(r"#[0-9A-Fa-f]{6}", theme.ACCENT_ON_DARK)
    assert flat["--theme.light.primaryColor"] == theme.ACCENT and flat["--theme.dark.primaryColor"] == theme.ACCENT
    assert "--theme.primaryColor" not in flat, "a colour without a variant would pin one look for every browser"
    assert json.loads(flat["--theme.chartCategoricalColors"]) == theme.CATEGORICAL
    assert json.loads(flat["--theme.chartSequentialColors"]) == theme.SEQUENTIAL
    assert len(json.loads(flat["--theme.headingFontSizes"])) == 3
    assert flat["--browser.gatherUsageStats"] == "false" and flat["--client.toolbarMode"] == "minimal"
    names = theme.streamlit_flags()[0::2]
    assert len(names) == len(set(names)) and all(isinstance(v, str) for v in theme.streamlit_flags()[1::2])


def test_a_theme_option_the_user_gives_wins_in_either_spelling():
    spaced = theme.streamlit_flags(["--theme.baseRadius", "1rem", "--server.port", "9000"])
    assert "--theme.baseRadius" not in spaced and "1rem" not in spaced and "--theme.buttonRadius" in spaced
    joined = theme.streamlit_flags(["--theme.baseRadius=1rem"])
    assert "--theme.baseRadius" not in joined and "--theme.buttonRadius" in joined
    assert "--browser.gatherUsageStats" not in theme.streamlit_flags(["--browser.gatherUsageStats=true"])
    assert theme.streamlit_flags(["--server.port", "9000"]) == theme.streamlit_flags(), "unrelated options change nothing"


def test_a_colour_the_user_names_replaces_ours_in_the_variants_it_covers():
    both = theme.streamlit_flags(["--theme.primaryColor", "#123456"])
    assert not any("primaryColor" in f for f in both[0::2]) and "--theme.light.greenColor" in both
    dark_only = theme.streamlit_flags(["--theme.dark.primaryColor=#123456"])
    assert "--theme.light.primaryColor" in dark_only and "--theme.dark.primaryColor" not in dark_only
    own = theme.streamlit_flags(["--theme.base", "dark"])
    assert not [f for f in own[0::2] if ".light." in f or ".dark." in f], "their own base: none of our colours"
    assert "--theme.baseRadius" in own


def test_the_stylesheet_follows_the_light_or_dark_variant_of_the_theme():
    css = style.stylesheet()
    assert "light-dark(" in css and theme.ACCENT_ON_DARK in css and "__" not in css
    assert f"--gb-accent-text: {theme.ACCENT};" in css, "browsers without light-dark() keep the plain accent"


def test_every_theme_option_is_one_streamlit_run_accepts():
    """Checked against the Streamlit that is installed: a name it does not know would stop ``grnbench ui``."""
    pytest.importorskip("streamlit")
    # Streamlit's help text has non-ASCII characters: on a Windows console (cp1252) the child dies writing
    # them to the pipe unless it is told to use UTF-8.
    env = {**os.environ, "PYTHONIOENCODING": "utf-8", "PYTHONUTF8": "1"}
    out = subprocess.run([sys.executable, "-m", "streamlit", "run", "--help"], capture_output=True, text=True,
                         encoding="utf-8", timeout=120, env=env)
    assert out.returncode == 0, out.stderr
    for name in theme.streamlit_flags()[0::2]:
        assert name in out.stdout, f"{name} is not an option of `streamlit run`"


def test_grnbench_ui_hands_the_theme_to_streamlit(monkeypatch):
    pytest.importorskip("streamlit")
    calls = []
    monkeypatch.setattr(subprocess, "call", lambda cmd, *a, **k: calls.append(cmd) or 0)
    assert cli_main.main(["ui", "--server.port", "8765", "--theme.primaryColor", "#123456"]) == 0
    cmd = calls[0]
    assert cmd[:4] == [sys.executable, "-m", "streamlit", "run"]
    # the same file, whichever way the path is spelled (a symlinked home on the cluster, #124)
    assert Path(cmd[4]).resolve() == APP.resolve()
    assert "--theme.baseRadius" in cmd and cmd.count("--theme.primaryColor") == 1 and cmd.count("#123456") == 1
    assert cmd[-4:] == ["--server.port", "8765", "--theme.primaryColor", "#123456"], "the user's own options come last"


# --------------------------------------------------------------------------- #
# the pages: headlines, menu, empty states, footers
# --------------------------------------------------------------------------- #

def test_every_screen_draws_its_own_headline():
    at = _at(timeout=60)
    at.run()
    for s in steps.STEPS:
        _goto(at, s.key)
        assert s.title in _text(at.main.markdown), s.key


def test_the_menu_groups_the_screens_marks_the_current_one_and_navigates():
    at = _at(timeout=60)
    at.run()
    nav = [b.key for b in at.sidebar.button if b.key.startswith("nav-")]
    assert len(nav) == len(steps.STEPS) and nav[0].startswith("nav-cur-")
    assert all(k.startswith("nav-todo-") for k in nav[1:]), "nothing is ticked before anything ran"
    groups = [re.sub(r"<[^>]+>", "", m.value) for m in at.sidebar.markdown if "gb-group" in m.value]
    assert groups == list(steps.GROUPS)
    assert any("optional" in b.label for b in at.sidebar.button if b.key.startswith("nav-")), "the optional step says so"
    at.button(key=f"nav-todo-{steps.LABELS.index(steps.DATA)}").click().run()
    assert not at.exception and at.session_state["step"] == steps.DATA
    assert at.button(key=f"nav-cur-{steps.LABELS.index(steps.DATA)}") is not None
    assert steps.DATA in _text(at.main.markdown) or steps.BY_KEY[steps.DATA].title in _text(at.main.markdown)


@pytest.mark.parametrize("step, title, target", [
    (steps.PROFILE, "No dataset yet", steps.DATA),
    (steps.PIPELINE, "No dataset yet", steps.DATA),
    (steps.EXPLORE_DATA, "No dataset yet", steps.DATA),
])
def test_a_screen_that_needs_data_says_so_and_offers_the_way_forward(step, title, target):
    at = _at(timeout=60)
    at.run()
    _goto(at, step)
    assert title in _text(at.main.markdown)
    assert {b.key for b in at.main.button} >= {"empty-primary", "empty-secondary"}
    assert len([b for b in at.sidebar.button if b.key.startswith("nav-")]) == len(steps.STEPS), "the menu is still drawn"
    at.button(key="empty-primary").click().run()
    assert not at.exception and at.session_state["step"] == target


def test_cleaned_data_and_models_are_asked_for_in_that_order(demo):
    at = _with_state(demo)
    _goto(at, steps.MODELS)
    assert "has not been cleaned" in _text(at.main.markdown)
    at.button(key="empty-primary").click().run()
    assert at.session_state["step"] == steps.CLEAN
    _goto(at, steps.BENCHMARK)
    assert "has not been cleaned" in _text(at.main.markdown) or "No model has been run" in _text(at.main.markdown)


def test_the_footer_buttons_walk_the_workflow_forward_and_back(demo):
    at = _with_state(demo)
    _goto(at, steps.DATA)
    nxt = at.button(key="foot-next")
    assert nxt.label.startswith("Continue") and steps.PROFILE in nxt.label
    nxt.click().run()
    assert not at.exception and at.session_state["step"] == steps.PROFILE
    back = at.button(key="foot-back")
    assert steps.DATA in back.label
    back.click().run()
    assert at.session_state["step"] == steps.DATA


def test_the_data_page_says_when_the_demo_brings_both_the_regulators_and_the_reference(demo, tmp_path):
    d, paths, tfs, res = demo
    at = _at(timeout=60)
    at.session_state["base"] = str(tmp_path)
    at.session_state["demo_paths"] = paths
    at.session_state["demo_info"] = "900 cells"
    at.run()
    _goto(at, steps.DATA)
    next(r for r in at.radio if r.label == "Dataset").set_value("Demo dataset").run()
    assert "The demo dataset brings both" in _text(at.main.markdown)
    assert at.checkbox(key="data_dorothea").disabled, "a reference is already given: DoRothEA would be ignored"
    next(r for r in at.radio if r.label == "Dataset").set_value("File or folder path").run()
    assert "Both are optional" in _text(at.main.markdown) and not at.checkbox(key="data_dorothea").disabled


def test_a_finished_quick_run_shows_the_results_first_and_folds_the_form(demo, finished):
    records, runs, bench = finished
    at = _with_state(demo, records=records, runs=runs, bench=bench, quick_done=True)
    _goto(at, steps.QUICK)
    assert "Analysis complete" in _text(at.main.success)
    assert [m.label for m in at.main.metric] == ["Cells used", "Genes", "Models run", "Highest AUPRC ratio here"]
    assert any("run again" in label for label in _expander_labels(at.main)), "the form is folded away after a run"
    assert {b.key for b in at.main.button} >= {"next-results", "next-data", "next-flow"}
    at.button(key="next-results").click().run()
    assert at.session_state["step"] == steps.EXPLORE_RESULTS and not at.exception


# --------------------------------------------------------------------------- #
# the plain-language summaries
# --------------------------------------------------------------------------- #

def _qc(*gates):
    return SimpleNamespace(gates=[SimpleNamespace(name=n, status=s, detail=d) for n, s, d in gates])


def test_plural_counts_things_in_words():
    assert [logic.plural(n, "model") for n in (0, 1, 2)] == ["0 models", "1 model", "2 models"]


def test_quality_checks_use_names_a_biologist_reads():
    t = logic.quality_checks(_qc(("mito_fraction", "warn", "x"), ("some_new_gate", "pass", "y")))
    assert list(t["check"]) == ["Mitochondrial RNA share", "Some new gate"]
    assert list(t.columns) == ["check", "status", "detail", "meaning"] and list(t["status"]) == ["warn", "pass"]
    assert "mitochondrial" in t.loc[0, "meaning"] and t.loc[1, "meaning"] == "", "an unknown check has no invented meaning"
    assert logic.quality_checks(_qc()).empty
    assert set(logic.GATE_MEANING) <= set(logic.GATE_LABELS), "every explained check has a plain name too"


def test_the_quality_verdict_is_one_honest_sentence():
    ok = logic.quality_verdict(_qc(("sparsity", "pass", ""), ("small_dataset", "pass", "")))
    assert ok["level"] == "ok" and ok["headline"] == "All 2 quality checks passed." and ok["attention"] == []
    warn = logic.quality_verdict(_qc(("sparsity", "pass", ""), ("mito_fraction", "warn", "high"), ("small_dataset", "pass", "")))
    assert warn["level"] == "warn" and warn["headline"] == "2 of 3 quality checks passed; 1 to keep an eye on."
    assert warn["counts"] == {"pass": 2, "warn": 1, "fail": 0}
    bad = logic.quality_verdict(_qc(("mito_fraction", "warn", "w"), ("zero_count_cells", "fail", "f"), ("sparsity", "pass", "")))
    assert bad["level"] == "fail" and "1 of 3 quality checks failed" in bad["headline"] and "1 more need a look" in bad["headline"]
    assert [a["status"] for a in bad["attention"]] == ["fail", "warn"], "what failed comes first"
    assert logic.quality_verdict(_qc())["headline"] == "No quality check was run."


def test_a_reason_is_shown_without_the_knowledge_base_prefix():
    assert logic.clean_reason("[normalization = cpm] Counts differ in depth.") == "Counts differ in depth."
    assert logic.clean_reason("  [x]   plain  ") == "plain" and logic.clean_reason(None) == ""
    assert logic.clean_reason("no prefix [kept]") == "no prefix [kept]"


def test_a_value_the_user_set_shows_what_the_engine_chose_and_why():
    from grnbench.autoprep import Decision, PrepPlan
    d = Decision("max_pct_mito_cell", "Max % mito", 24.8, "median 17.9% mito; median + 3 MADs = 24.8%", "data")
    assert logic.decision_caption(d) == "median 17.9% mito; median + 3 MADs = 24.8% (source: data)"
    d.value, d.reason, d.source = 12.5, "set by the user in the interface", "user"
    d.auto_value, d.auto_reason, d.auto_source = 24.8, "median 17.9% mito; median + 3 MADs = 24.8%", "data"
    assert logic.decision_caption(d) == "auto: 24.8 (median 17.9% mito; median + 3 MADs = 24.8%) → you: 12.5"
    off = Decision("max_genes_per_cell", "Max genes", 3000, "set by the user in the interface", "user",
                   auto_value=None, auto_reason="no spread in detected genes: no upper bound", auto_source="data")
    assert logic.decision_caption(off) == "auto: off (no spread in detected genes: no upper bound) → you: 3000"
    plan = PrepPlan(choices={"cell_qc_thresholds": "fixed", "doublets": "none"},
                    selection={"results": {"cell_qc_thresholds": {"chosen": "mad"}, "doublets": {"chosen": "scrublet"}}})
    fixed = Decision("min_genes_per_cell", "Min genes", 200, "[cell_qc_thresholds = fixed] fixed bounds", "user",
                     auto_value=397, auto_reason="median 597 genes/cell; lower 3-MAD bound", auto_source="data")
    assert logic.decision_caption(fixed, plan) == \
        "auto: 397 (median 597 genes/cell; lower 3-MAD bound) → you: 200 (chosen: fixed)"
    same = PrepPlan(choices={"cell_qc_thresholds": "mad"}, selection={"results": {"cell_qc_thresholds": {"chosen": "mad"}}})
    assert logic.decision_caption(d, same) == logic.decision_caption(d)        # a typed value names no option
    doublets = Decision("drop_doublets", "Remove doublets", False, "[doublets = none] none", "user",
                        auto_value=True, auto_reason="[doublets = scrublet] droplet data", auto_source="knowledge base")
    assert logic.decision_caption(doublets, plan) == \
        "auto: True (droplet data; source: knowledge base) → you: False (chosen: none)"


def test_a_pasted_path_loses_its_quotes_and_spaces_and_expands_the_home_folder():
    assert logic.clean_path('  "C:\\data\\set.h5ad" ') == Path("C:\\data\\set.h5ad")
    assert logic.clean_path("'/data/x.h5ad'") == Path("/data/x.h5ad")
    assert logic.clean_path("~/x.h5ad") == Path.home() / "x.h5ad"
    assert logic.clean_path("") is None and logic.clean_path(None) is None and logic.clean_path("   ") is None


def test_the_pipeline_overview_has_one_row_per_decision_in_the_knowledge_base_order(demo):
    plan = demo[3]["plan"]
    rows = logic.pipeline_overview(plan)
    assert rows and {r["id"] for r in rows} == set(plan.selection["results"])
    assert len({r["id"] for r in rows}) == len(rows)
    for r in rows:
        assert {"id", "stage", "question", "impact", "choice", "uncertain"} <= set(r)
        assert isinstance(r["choice"], str) and r["choice"] and isinstance(r["uncertain"], bool)
    stages = [r["stage"] for r in rows]
    assert stages == sorted(stages, key=stages.index), "rows of one stage stay together"


def test_the_demo_summary_counts_cells_genes_regulators_and_known_links(demo):
    paths = demo[1]
    gold_rows = len(pd.read_csv(paths["gold"], sep="\t"))
    assert logic.demo_summary(paths) == f"900 cells x 700 genes, 12 regulators, {gold_rows} known links"


def test_the_cleaning_table_starts_from_the_size_as_loaded(demo, finished):
    records, _, _ = finished
    t = logic.retention_table(records)
    for _, one in t.groupby("context", sort=False):
        assert one.iloc[0]["step"] == "loaded" and one.iloc[0]["detail"] == "as read from the file"
        assert one["cells"].is_monotonic_decreasing, "cleaning never adds cells"
    assert t.iloc[0]["cells"] == 900 and t.iloc[0]["genes"] == 700


def test_the_run_summary_gives_the_headline_numbers_and_never_a_ranking_claim(finished):
    records, runs, bench = finished
    s = logic.run_summary(records, runs, bench)
    assert s["cells_loaded"] == 900 and 0 < s["cells_used"] <= 900 and s["genes_used"] <= 700
    assert s["contexts"] == len(records) and s["models_ok"] == 2 and s["models_total"] == 2
    assert s["top"]["ratio"] == pytest.approx(bench["AUPRC ratio"].max()) and s["top"]["method"] in {"pearson_signed", "spearman_signed"}
    assert s["reference_edges"] and s["reference_edges"] > 0
    empty = logic.run_summary([], None)
    assert empty["models_total"] == 0 and empty["top"] is None and empty["cells_used"] is None
    unscored = logic.run_summary(records, runs, None)
    assert unscored["top"] is None and unscored["models_ok"] == 2


def test_the_score_chart_puts_the_best_model_on_top_and_marks_the_random_level():
    bench = pd.DataFrame({"method": ["a", "b", "c"], "context": ["all_cells"] * 3, "AUPRC ratio": [3.0, 9.0, 6.0],
                          "AUROC": [0.55, 0.9, 0.7], "EPR": [2.0, 8.0, 5.0]})
    fig = logic.score_figure(bench)
    assert list(fig.layout.yaxis.categoryarray) == ["a", "c", "b"], "best last = on top of a horizontal chart"
    assert any(t.name == "random (1)" for t in fig.data), "the random level is named in the legend"
    assert fig.data[0].marker.color == theme.ACCENT, "one context: the product accent, not a palette"
    auroc = logic.score_figure(bench, "AUROC")
    assert any(t.name == "random (0.5)" for t in auroc.data)
    assert list(logic.score_figure(bench.assign(**{"AUPRC ratio": [3.0, None, 6.0]})).layout.yaxis.categoryarray) == ["a", "c"], \
        "a model without a score is left out, not drawn as zero"
    two = bench.assign(context=["x", "y", "x"])
    assert len({t.name for t in logic.score_figure(two).data if t.name in ("x", "y")}) == 2, "several contexts are told apart"


def test_the_network_figure_takes_the_page_background_so_names_stay_readable_in_dark_mode():
    edges = pd.DataFrame({"source": ["T1", "T1", "T2"], "target": ["G1", "G2", "G1"], "weight": [3.0, 2.0, 1.0]})
    assert logic.network_figure(edges).layout.plot_bgcolor == "rgba(0,0,0,0)"
    assert logic.network_figure(edges.iloc[0:0]).layout.plot_bgcolor == "rgba(0,0,0,0)"


def test_the_eyebrow_repeats_the_menu_name_so_the_title_and_the_menu_match_up():
    assert [steps.BY_KEY[k].short for k in steps.WORKFLOW_KEYS] == [
        "Data", "Profile", "Pipeline", "Compare pipelines", "Clean", "Models", "Benchmark"]
    assert steps.BY_KEY[steps.EXPLORE_DATA].short == "Data" and steps.BY_KEY[steps.EXPLORE_RESULTS].short == "Results"
    assert "Explore · Results" in style.header_html(steps.EXPLORE_RESULTS)
    assert "Step 7 of 7 · Benchmark" in style.header_html(steps.BENCHMARK)
    assert "Start here<" in style.header_html(steps.QUICK)


def test_a_warning_about_an_option_that_cannot_run_is_put_in_words():
    raw = "ambient_rna: best option 'DecontX' is unavailable (missing: r:celda); using 'none'."
    said = logic.plain_warning(raw)
    assert "DecontX" in said and "r:celda" in said and "ambient_rna" not in said and "cannot run here" in said
    assert "?" in said, "it starts from the question the decision answers"
    other = "Sources disagree on organism: {'human': 2} - using 'human'."
    assert logic.plain_warning(other) == other, "any other warning is shown as the engine wrote it"
    assert logic.plain_warning("nosuchdecision: best option 'X' is unavailable (missing: y); using 'z'.").startswith("nosuch")


def test_the_pipeline_page_shows_warnings_and_option_names_in_words(demo):
    d, paths, tfs, res = demo
    plan = res["plan"].__class__.from_dict(res["plan"].to_dict())
    plan.warnings = ["ambient_rna: best option 'DecontX' is unavailable (missing: r:celda); using 'none'."]
    at = _with_state(demo, plan=plan)
    _goto(at, steps.PIPELINE)
    assert "cannot run here" in _text(at.main.warning) and "ambient_rna:" not in _text(at.main.warning)
    from grnbench.kb import DECISIONS
    pick = at.selectbox(key="pick_normalization")
    shown = pick.format_func(pick.value)
    assert shown == next(d for d in DECISIONS if d.id == "normalization").option(pick.value).name != pick.value, \
        "options show their names, not their ids"


def test_one_primary_button_at_a_time_on_the_data_page(demo):
    fresh = _at(timeout=60)
    fresh.run()
    _goto(fresh, steps.DATA)
    assert next(b for b in fresh.main.button if b.label == "Analyse dataset").proto.type == "primary", "nothing analysed yet"
    at = _with_state(demo)
    _goto(at, steps.DATA)
    analyse = next(b for b in at.main.button if b.label == "Analyse dataset")
    assert analyse.proto.type == "secondary", "the analysis is done: Continue is the primary action now"
    assert at.button(key="foot-next").proto.type == "primary"


def test_the_quick_run_says_what_it_scored_on_and_how_to_read_the_chart(demo, finished):
    records, runs, bench = finished
    at = _with_state(demo, records=records, runs=runs, bench=bench, quick_done=True,
                     bench_on="full reference (the pipeline was not chosen on it) · 444 links")
    _goto(at, steps.QUICK)
    captions = _text(at.main.caption)
    assert "Scored on: full reference" in captions and "times better than a random ranking" in captions
    assert at.main.metric[0].delta.startswith("of "), "cells used says how many were loaded"


def test_cleaning_again_clears_what_was_scored_on_the_old_cleaning(demo, finished):
    records, runs, bench = finished
    at = _with_state(demo, records=records, runs=runs, bench=bench, bench_on="x · 1 links")
    _goto(at, steps.CLEAN)
    next(b for b in at.main.button if b.label == "Run cleaning").click().run()
    assert not at.exception and "bench_on" not in at.session_state and "bench" not in at.session_state
