"""Interface regressions found by driving the whole tool as a newcomer (issue #70): each test failed on the
interface before its fix. Pages are rendered with Streamlit's AppTest (no browser); the logic is called directly."""
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from grnbench.autoprep import execute_plan
from grnbench.demo import write_demo
from grnbench.ui import logic

ROOT = Path(__file__).resolve().parents[1]
APP = ROOT / "src" / "grnbench" / "ui" / "app.py"


def _at(timeout=180):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    return AppTest.from_file(str(APP), default_timeout=timeout)


@pytest.fixture(scope="module")
def demo(tmp_path_factory):
    d = tmp_path_factory.mktemp("uireg")
    paths = write_demo(d / "demo", n_cells=900, n_genes=700, n_tfs=12, targets_per_tf=10)
    tfs = paths["tfs"].read_text().split()
    res = logic.analyze(paths["counts"], tf_list=tfs, choices={"context_scope": "pooled"})
    return d, paths, tfs, res


def _with_state(demo, **extra):
    d, paths, tfs, res = demo
    at = _at()
    for k, v in dict(workdir=d, data_path=str(paths["counts"]), tfs=tfs, reference=None, diag=res["diag"],
                     plan=res["plan"].__class__.from_dict(res["plan"].to_dict()), **extra).items():
        at.session_state[k] = v
    at.run()
    assert not at.exception, [e.value for e in at.exception]
    return at


def _goto(at, step):
    at.session_state["step"] = step
    at.run()
    assert not at.exception, (step, [e.value for e in at.exception])


# --------------------------------------------------------------------------- #
# the dataset choosers
# --------------------------------------------------------------------------- #

def test_step1_file_path_source_analyses_the_dataset(demo, monkeypatch):
    """The path box of step 1 shared its key with the session value of the same name, so Analyse always crashed."""
    _, paths, _, _ = demo
    at = _at()
    at.run()
    _goto(at, "1 · Data")
    next(r for r in at.radio if r.label == "Dataset").set_value("File or folder path").run()
    next(t for t in at.text_input if t.label.startswith("Path to")).set_value(str(paths["counts"])).run()
    next(b for b in at.button if b.label == "Analyse dataset").click().run()
    assert not at.exception, [e.value for e in at.exception]
    assert any("Analysis done" in s.value for s in at.success)
    assert at.session_state["plan"] is not None and at.session_state["data_path"] == str(paths["counts"])


def test_a_table_whose_layout_cannot_be_guessed_can_be_read_by_naming_the_layout(tmp_path):
    """The engine refuses to guess whether rows are cells or genes; the interface had no way to say."""
    rng = np.random.default_rng(0)
    genes = [f"G{i:03d}" for i in range(300)]
    cells = [f"c{i:03d}" for i in range(200)]
    path = tmp_path / "counts_small.csv"
    pd.DataFrame(rng.poisson(2.0, size=(300, 200)), index=genes, columns=cells).to_csv(path)
    with pytest.raises(Exception, match="orientation"):
        logic.analyze(path)
    res = logic.analyze(path, orientation="genes_x_cells")
    assert res["plan"].value("orientation") == "genes_x_cells"
    at = _at()
    at.run()
    next(r for r in at.radio if r.label == "Dataset").set_value("File or folder path").run()
    assert any(s.label == "CSV / TSV layout" for s in at.selectbox), "the layout is a visible control"


# --------------------------------------------------------------------------- #
# a run that did nothing must not look finished
# --------------------------------------------------------------------------- #

def test_quick_run_refuses_when_no_context_could_be_cleaned(tmp_path):
    paths = write_demo(tmp_path / "demo", n_cells=900, n_genes=700, n_tfs=12, targets_per_tf=10)
    with pytest.raises(RuntimeError, match="Nothing was run"):
        logic.quick_run(paths["counts"], tmp_path / "runs", ["pearson_signed"], tf_list=["NOTAGENE1", "NOTAGENE2"])


def test_run_models_with_nothing_to_run_keeps_its_columns(tmp_path):
    runs = logic.run_models([{"context": "all_cells", "processed": None, "error": "x", "out_dir": str(tmp_path)}],
                            ["pearson_signed"], tmp_path / "models")
    assert list(runs.columns) == logic.RUN_COLUMNS and runs.empty


def test_the_sidebar_ticks_only_what_produced_something(demo):
    at = _with_state(demo, records=[{"context": "all_cells", "processed": None, "error": "boom"}],
                     runs=pd.DataFrame(columns=logic.RUN_COLUMNS))
    _goto(at, "Explore data")
    done = {b.label.split("  ")[0] for b in at.sidebar.button if b.key.startswith("nav-done-")}
    assert "1 · Data" in done, "the dataset was analysed"
    assert not done & {"5 · Clean", "6 · Models", "7 · Benchmark"}, done


# --------------------------------------------------------------------------- #
# bad uploads and bad references say what is wrong
# --------------------------------------------------------------------------- #

def test_a_corrupt_tf_list_is_one_red_line_not_a_traceback(demo):
    _, paths, _, _ = demo
    at = _at()
    at.run()
    next(r for r in at.radio if r.label == "Dataset").set_value("File or folder path").run()
    next(t for t in at.text_input if t.label.startswith("Path to")).set_value(str(paths["counts"])).run()
    at.file_uploader(key="quick_tf").set_value(("corrupt.csv", b"\xc3\x28\xff\xfe\x00\x81", "text/csv")).run()
    next(b for b in at.button if b.label == "Run").click().run()
    assert not at.exception, [e.value for e in at.exception]
    assert any("Could not read corrupt.csv" in e.value for e in at.error), [e.value for e in at.error]


def test_a_reference_that_shares_nothing_with_the_data_is_explained(demo):
    d, paths, tfs, _ = demo
    runs_dir = d / "ui13"
    recs = execute_plan(demo[3]["plan"], paths["counts"], runs_dir / "clean", tf_list=tfs)
    runs = logic.run_models(recs, ["pearson_signed"], runs_dir / "models", tf_list=tfs)
    gold = logic.read_pairs(paths["gold"])
    assert logic.benchmark(runs, gold, tf_list=tfs)["reference edges in universe"].gt(0).all()
    with pytest.raises(ValueError, match="no reference edge"):
        logic.benchmark(runs, {("NOPE1", "NOPE2")}, tf_list=tfs)
    out = logic.quick_run(paths["counts"], runs_dir / "q_nomatch", ["pearson_signed"], tf_list=tfs,
                          reference={("NOPE1", "NOPE2")})
    assert out["bench"] is None and out["reference_warning"]
    assert out["runs"]["error"].isna().all() and out["runs"]["n_edges"].gt(0).all(), "the models' results are kept"
    out = logic.quick_run(paths["counts"], runs_dir / "q_empty", ["pearson_signed"], tf_list=tfs, reference=set())
    assert out["bench"] is None and "no edges" in out["reference_warning"]
    pairs, label, warn = logic.reference_for_scoring(None, gold, set(), "header_only.tsv")
    assert pairs is None and "no edges" in warn and "header_only.tsv" in warn, "an empty upload is not swapped"


# --------------------------------------------------------------------------- #
# step 3: what is pending is what is on the page
# --------------------------------------------------------------------------- #

def _apply_labels(at):
    return [b.label for b in at.button if b.label.startswith("Apply")]


def test_leaving_step3_does_not_leave_a_hidden_pending_change(demo):
    at = _with_state(demo, user_choices={})
    _goto(at, "3 · Pipeline")
    assert _apply_labels(at) == []
    at.selectbox(key="pick_normalization").set_value("pearson_residuals").run()
    assert _apply_labels(at) == ["Apply 1 change(s) and rebuild the plan"]
    _goto(at, "4 · Compare pipelines")
    _goto(at, "3 · Pipeline")
    shown = at.selectbox(key="pick_normalization").value
    assert shown == "pearson_residuals" or _apply_labels(at) == [], (shown, _apply_labels(at))


def test_picking_the_original_option_again_withdraws_the_change(demo):
    at = _with_state(demo, user_choices={})
    _goto(at, "3 · Pipeline")
    at.selectbox(key="pick_normalization").set_value("pearson_residuals").run()
    at.selectbox(key="pick_normalization").set_value("cpm_log1p").run()
    assert _apply_labels(at) == []


def test_a_second_apply_keeps_the_first_and_an_adopted_pipeline(demo):
    from grnbench.autoprep import plan_preprocessing
    at = _with_state(demo, user_choices={})
    _goto(at, "3 · Pipeline")
    at.selectbox(key="pick_normalization").set_value("pearson_residuals").run()
    [b for b in at.button if b.label.startswith("Apply")][0].click().run()
    assert _apply_labels(at) == [], "nothing stays pending after Apply"
    at.selectbox(key="pick_batch").set_value("center_per_batch").run()
    [b for b in at.button if b.label.startswith("Apply")][0].click().run()
    ch = at.session_state["plan"].choices
    assert ch["batch"] == "center_per_batch" and ch["normalization"] == "pearson_residuals"
    adopted = plan_preprocessing(demo[3]["diag"], tf_list=demo[2], choices={"batch": "center_per_batch"})
    at = _with_state(demo, user_choices={})
    at.session_state["plan"] = adopted               # what step 4 does with "Use this pipeline"
    _goto(at, "3 · Pipeline")
    at.selectbox(key="pick_normalization").set_value("pearson_residuals").run()
    [b for b in at.button if b.label.startswith("Apply")][0].click().run()
    ch = at.session_state["plan"].choices
    assert ch["normalization"] == "pearson_residuals" and ch["batch"] == "center_per_batch"


def test_step4_says_when_the_plan_was_edited_in_step3(demo):
    at = _with_state(demo, user_choices={})
    _goto(at, "3 · Pipeline")
    at.selectbox(key="pick_normalization").set_value("pearson_residuals").run()
    [b for b in at.button if b.label.startswith("Apply")][0].click().run()
    _goto(at, "4 · Compare pipelines")
    text = " ".join(w.value for w in at.warning)
    assert "pearson_residuals" in text and "recommendation" in text and "replaces" in text, text
    d, paths, tfs, _ = demo        # the shared plan overrides the scope; the engine's own recommendation does not
    at = _with_state(demo, user_choices={})
    at.session_state["plan"] = logic.analyze(paths["counts"], tf_list=tfs)["plan"]
    _goto(at, "4 · Compare pipelines")
    assert not [w for w in at.warning if "recommendation" in w.value], "quiet when the plan is the recommendation"


# --------------------------------------------------------------------------- #
# the model picker and the pages around it
# --------------------------------------------------------------------------- #

def test_typed_inputs_survive_a_visit_to_another_page():
    at = _at(60)
    at.run()
    at.radio(key="quick_source").set_value("File or folder path").run()
    at.text_input(key="quick_pathbox").set_value("some/data.h5ad").run()
    at.radio(key="quick_preset").set_value("Fast, including extras").run()
    for other in ("1 · Data", "2 · Profile"):          # 2 · Profile ends with st.stop() (no data yet)
        _goto(at, other)
        _goto(at, "Quick run")
        assert at.radio(key="quick_source").value == "File or folder path"
        assert at.text_input(key="quick_pathbox").value == "some/data.h5ad"
        assert at.radio(key="quick_preset").value == "Fast, including extras"


def _register_ridge(at):
    src = ROOT / "examples" / "plugins" / "ridge_model.py"
    next(w for w in at.file_uploader if w.label == "Model file").set_value((src.name, src.read_bytes(), "text/x-python"))
    next(w for w in at.text_input if w.label == "Class or function name").set_value("RidgeGRN")
    at.run()
    next(b for b in at.button if b.label == "Register").click().run()
    assert not at.exception and not at.error


def test_registering_a_model_confirms_it_and_none_selects_nothing(demo):
    d, paths, tfs, res = demo
    recs = execute_plan(res["plan"], paths["counts"], d / "ui18" / "clean", tf_list=tfs)
    at = _with_state(demo, records=recs, user_choices={})
    _goto(at, "6 · Models")
    _register_ridge(at)
    assert any("Registered: ridge" in s.value for s in at.success), "the confirmation must stay on screen"
    next(r for r in at.radio if r.label == "Start from").set_value("None").run()
    captions = [c.value for c in at.caption if "selected" in c.value]
    assert captions and captions[0].startswith("0 models selected"), captions
    assert next(b for b in at.button if b.label == "Run models").disabled


# --------------------------------------------------------------------------- #
# the output folder and the Stop button
# --------------------------------------------------------------------------- #

def test_output_folder_is_checked_once_with_a_clear_message(tmp_path):
    assert logic.output_folder(tmp_path / "new" / "deep").is_dir()
    with pytest.raises(ValueError, match="empty"):
        logic.output_folder("   ")
    blocker = tmp_path / "a_file"
    blocker.write_text("x")
    with pytest.raises(ValueError, match="Cannot write"):
        logic.output_folder(blocker / "inside")


def test_stop_kills_the_model_process(demo, tmp_path):
    """Streamlit's Stop only raises at the next st.* call; the polling progress call is that call, and the
    model runs in a child process that is killed when it raises, so nothing keeps computing invisibly."""
    import multiprocessing

    d, paths, tfs, res = demo
    slow = tmp_path / "slow_model.py"
    slow.write_text("import time\nimport pandas as pd\n\n\ndef infer(adata, tf_list=None, **params):\n"
                    "    time.sleep(600)\n    return pd.DataFrame(columns=['source', 'target', 'weight'])\n")
    logic.register_user_model(slow, "infer")
    recs = execute_plan(res["plan"], paths["counts"], d / "stop" / "clean", tf_list=tfs)

    class Stopped(BaseException):          # what Streamlit raises inside the script thread
        pass

    seen = []

    def progress(frac, msg):
        seen.append(msg)
        if "infer on" in msg:
            raise Stopped()

    t0 = time.time()
    with pytest.raises(Stopped):
        logic.run_models(recs, ["infer"], d / "stop" / "models", tf_list=tfs, progress=progress)
    assert time.time() - t0 < 120, "the call returned as soon as Stop was raised"
    assert multiprocessing.active_children() == [], "the child process is gone"
    assert any("infer on" in m for m in seen)


def test_step4_refuses_no_op_and_repeated_alternatives(demo):
    at = _with_state(demo, user_choices={})
    _goto(at, "4 · Compare pipelines")
    add = next(b for b in at.button if b.label == "Add")
    add.click().run()            # the preselected decision/option pair is the recommended choice itself
    assert any("nothing to compare" in w.value for w in at.warning), [w.value for w in at.warning]
    assert "extra_candidates" not in at.session_state or not at.session_state["extra_candidates"]
    next(s for s in at.selectbox if s.label == "Decision").set_value("normalization").run()
    next(s for s in at.selectbox if s.label == "Option").set_value("pearson_residuals").run()
    next(b for b in at.button if b.label == "Add").click().run()
    assert [c["name"] for c in at.session_state["extra_candidates"]] == ["normalization=pearson_residuals"]
    next(b for b in at.button if b.label == "Add").click().run()
    assert any("already in the comparison" in w.value for w in at.warning)
    assert len(at.session_state["extra_candidates"]) == 1
    assert any("cellular context" in m.value for m in at.markdown), "the preview names what is always tried"


# --------------------------------------------------------------------------- #
# a screening winner that scored on fewer cell types than the rest
# --------------------------------------------------------------------------- #

def _screen_that_dropped_cell_types(res):
    """A screening result in which context_scope=selected_type won on alpha alone (the demo's three cell types)."""
    plan_cls = res["plan"].__class__

    def plan_with(*names):
        plan = plan_cls.from_dict(res["plan"].to_dict())
        plan.contexts = [{"name": n} for n in names]
        return plan
    cts = [f"cell_type={c}" for c in ("alpha", "beta", "gamma")]
    best = "context_scope=selected_type"
    rows = pd.DataFrame({"candidate": ["recommended"] * 3 + [best], "context": cts + cts[:1],
                         "n_cells": [300, 200, 100, 300], "reference_auprc_ratio": [6.0, 5.0, 4.0, 6.0]})
    ranking = pd.DataFrame({"candidate": [best, "recommended"], "score": [6.0, 5.0]})
    return dict(best=best, rule="highest mean score", rows=rows, ranking=ranking, pros_cons={},
                plans={"recommended": plan_with(*cts), best: plan_with(cts[0])}), cts


def test_a_screening_winner_that_drops_cell_types_says_so(demo):
    screen, cts = _screen_that_dropped_cell_types(demo[3])
    scope = logic.screening_scope(screen)
    assert scope["kept"] == cts[:1] and scope["left_out"] == cts[1:]
    pooled = {**screen, "best": "baseline", "plans": {**screen["plans"], "baseline": screen["plans"]["recommended"].__class__
                                                      .from_dict(demo[3]["plan"].to_dict())}}
    pooled["plans"]["baseline"].contexts = [{"name": "all_cells"}]
    assert logic.screening_scope(pooled)["left_out"] == [], "a pooled baseline does not narrow anything"
    assert logic.screening_scope(None)["left_out"] == [] and logic.screening_scope({})["kept"] == []
    at = _with_state(demo, screen=screen)
    _goto(at, "4 · Compare pipelines")
    text = " ".join(w.value for w in at.warning)
    assert "keeps only cell_type=alpha" in text and "cell_type=beta" in text and "cell_type=gamma" in text, text


# --------------------------------------------------------------------------- #
# a second dataset in the same session
# --------------------------------------------------------------------------- #

def _metric(at, label):
    return next(m.value for m in at.metric if m.label == label)


def _caches(at):
    return {k: v for k, v in at.session_state.items() if k.startswith("_") and k.endswith("_cache")}


def test_a_second_dataset_never_shows_what_was_cached_for_the_first(demo, tmp_path):
    """B1 (PR #71 review): the explorer's caches outlived the dataset they were drawn from, so a second dataset
    analysed in the same session could be shown with the first one's cells, map and annotations. Look at one
    dataset, analyse another in the same session: the explorer shows the second and keeps nothing of the first."""
    d, _, _, res = demo
    other = write_demo(tmp_path / "other", n_cells=500, n_genes=400, n_tfs=8, targets_per_tf=8)
    at = _with_state(demo, base=str(tmp_path / "runs"))
    _goto(at, "Explore data")
    first = res["diag"].adata
    assert (_metric(at, "Cells"), _metric(at, "Genes")) == (f"{first.n_obs:,}", f"{first.n_vars:,}")
    assert {"_stage_cache", "_embed_cache", "_summary_cache"} <= set(_caches(at)), "the first dataset was drawn"
    _goto(at, "1 · Data")
    next(r for r in at.radio if r.label == "Dataset").set_value("File or folder path").run()
    next(t for t in at.text_input if t.label.startswith("Path to")).set_value(str(other["counts"])).run()
    next(b for b in at.button if b.label == "Analyse dataset").click().run()
    assert not at.exception and not at.error, [e.value for e in [*at.exception, *at.error]]
    second = at.session_state["diag"].adata
    assert second.shape != first.shape and str(at.session_state["workdir"]) != str(d)
    _goto(at, "Explore data")
    assert not at.error, [e.value for e in at.error]
    assert (_metric(at, "Cells"), _metric(at, "Genes")) == (f"{second.n_obs:,}", f"{second.n_vars:,}")
    maps = list(at.session_state["_embed_cache"].values())
    assert maps and all(f.attrs["n_cells_total"] == second.n_obs for f in maps)
    summaries = list(at.session_state["_summary_cache"].values())
    assert summaries and all(len(s) == second.n_vars for s in summaries)
    stale = [(name, key) for name, cache in _caches(at).items() for key in cache
             if isinstance(key, tuple) and key[0] == str(d)]
    assert stale == [], "nothing cached for the first dataset is kept"


# --------------------------------------------------------------------------- #
# a second run in the same run folder
# --------------------------------------------------------------------------- #

def test_a_second_run_gets_its_own_folder_instead_of_the_engine_refusing_to_write(tmp_path):
    """After a quick run, 'Run cleaning' in step 5 failed with 'Output directory must be empty' (the files of the
    first run were in the way). The interface never deletes: it writes the next run beside the first."""
    first = logic.fresh_folder(tmp_path / "clean")
    assert first == tmp_path / "clean"
    first.mkdir()
    assert logic.fresh_folder(tmp_path / "clean") == tmp_path / "clean"      # an empty folder is still free
    (first / "all_cells").mkdir()
    (first / "all_cells" / "processed.h5ad").write_bytes(b"x")
    second = logic.fresh_folder(tmp_path / "clean")
    assert second == tmp_path / "clean-2"
    second.mkdir()
    (second / "f").write_bytes(b"y")
    assert logic.fresh_folder(tmp_path / "clean") == tmp_path / "clean-3"
    assert (first / "all_cells" / "processed.h5ad").read_bytes() == b"x", "nothing of the first run was touched"


def test_cleaning_again_after_a_quick_run_works(demo):
    d, paths, tfs, res = demo
    at = _with_state(demo)
    at.session_state["workdir"] = str(d / "again")
    _goto(at, "5 · Clean")
    next(b for b in at.button if b.label == "Run cleaning").click().run()
    assert not at.exception and not at.error, [e.value for e in at.error]
    again = next(b for b in at.button if b.label == "Run cleaning")
    again.click().run()
    assert not at.exception and not at.error, [e.value for e in at.error]
    assert (d / "again" / "clean").is_dir() and (d / "again" / "clean-2").is_dir()
