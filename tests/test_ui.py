"""Interface: logic end to end on the demo, and the Streamlit app renders."""
from pathlib import Path

import pandas as pd
import pytest

from grnbench.autoprep import execute_plan
from grnbench.demo import write_demo
from grnbench.ui import logic

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def analysed(tmp_path_factory):
    d = tmp_path_factory.mktemp("ui")
    paths = write_demo(d / "demo", n_cells=900, n_genes=700, n_tfs=12, targets_per_tf=10)
    tfs = paths["tfs"].read_text().split()
    res = logic.analyze(paths["counts"], tf_list=tfs, choices={"context_scope": "pooled"})
    return d, paths, tfs, res


def test_profile_tables_and_figures(analysed):
    _, _, _, res = analysed
    diag, plan = res["diag"], res["plan"]
    facts = logic.facts_table(diag.context)
    assert {"fact", "value", "source", "evidence"} <= set(facts.columns) and len(facts)
    rows = logic.decision_rows(plan, "doublets")
    assert rows["chosen"].sum() == 1 and {"score", "available", "references"} <= set(rows.columns)
    figs = logic.qc_figures(diag, plan)
    assert {"detected genes", "log10 counts per cell", "% mitochondrial", "composition"} <= set(figs)
    assert dict(logic.summarize_profile(plan.profile))["Cells"] == "900"


def test_models_including_a_user_model_are_run_and_benchmarked(analysed):
    d, paths, tfs, res = analysed
    recs = execute_plan(res["plan"], paths["counts"], d / "clean", tf_list=tfs)
    reg = logic.register_user_model(ROOT / "examples" / "plugins" / "ridge_model.py", "RidgeGRN")
    assert reg[0]["name"] == "ridge"
    runs = logic.run_models(recs, ["pearson_signed", "ridge", "no_such_method"], d / "models", tf_list=tfs)
    assert runs.loc[runs["method"] == "no_such_method", "error"].notna().all()
    assert runs.loc[runs["method"] != "no_such_method", "edges"].notna().all()
    retention = logic.retention_table(recs)
    assert {"filter_cells", "strategy"} <= set(retention["step"])
    gold = pd.read_csv(paths["gold"], sep="\t")
    bench = logic.benchmark(runs, set(zip(gold.source, gold.target)), tf_list=tfs)
    assert set(bench["method"]) == {"pearson_signed", "ridge"}
    assert (bench["AUPRC ratio"] > 1).all()      # planted network is recoverable
    assert logic.zip_folder(d / "clean")[:2] == b"PK"


def test_app_renders_every_page_without_data():
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(str(ROOT / "src" / "grnbench" / "ui" / "app.py"), default_timeout=60)
    at.run()
    assert not at.exception
    for step in at.sidebar.radio[0].options:
        at.sidebar.radio[0].set_value(step).run()
        assert not at.exception, step
