"""Interface: logic end to end on the demo, and the Streamlit app renders."""
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

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


def test_reproduce_commands_give_one_shell_safe_line_per_context_with_its_seed():
    import shlex
    from grnbench.autoprep import PrepPlan, unique_folder_names
    # real labels carry spaces, ':' and '|' (AGENTS.md section 13); "CD4 T" and "CD4:T" share one safe folder name
    names = ["CD4 T", "CD4:T", "T|reg", "NK(bright)", 'x"y', "a$b", "-control"]
    plan = PrepPlan(contexts=[{"name": n, "subset_key": "cell_type", "subset_value": n} for n in names])
    cmds = logic.reproduce_commands(plan, "data.h5ad", tf_list=["TF1"], seed=0)
    assert len(cmds) == len(names) and "'--subset-value=CD4 T'" in cmds[0]
    for cmd, ctx, folder in zip(cmds, plan.contexts, unique_folder_names(names)):
        argv = shlex.split(cmd)
        assert argv == ["grnbench", "prep", *plan.cli_args("data.h5ad", "OUT/" + folder, "tf_list.txt", context=ctx,
                                                           bounds_path="per_batch_bounds.json", seed=0)]
        assert f"--subset-value={ctx['subset_value']}" in argv and argv[-2:] == ["--seed", "0"]


def test_shell_command_pastes_back_to_the_same_argv():
    import shlex
    argv = ["grnbench", "prep", r"C:\data\my counts.h5ad", "--subset-value", "T|reg", "--x", "a$b`c'd\"e"]
    assert shlex.split(logic.shell_command(argv)) == argv


def test_analyze_and_reproduction_commands_thread_an_explicit_gene_mapping(tmp_path):
    ids = [f"ENSG{i:011d}" for i in range(250)]
    src = tmp_path / "ensembl.h5ad"
    rng = np.random.default_rng(8)
    counts = np.vstack([
        rng.poisson(depth, len(ids)) for depth in np.linspace(1.0, 5.0, 20)
    ]).astype(np.float32)
    ad.AnnData(
        sp.csr_matrix(counts),
        var=pd.DataFrame(index=ids),
    ).write_h5ad(src)
    mapping = tmp_path / "gene map.tsv"
    mapping.write_text(
        "id\tsymbol\n" + "".join(f"{gene_id}\tGENE{i}\n" for i, gene_id in enumerate(ids)),
        encoding="utf-8",
    )

    res = logic.analyze(src, gene_mapping_path=mapping)
    commands = logic.reproduce_commands(
        res["plan"], src, gene_mapping_path=mapping, seed=0)

    assert list(res["diag"].adata.var_names[:3]) == ["GENE0", "GENE1", "GENE2"]
    assert commands and "--gene-mapping '" in commands[0]
    assert str(mapping) in commands[0]


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
    from grnbench.ui import steps
    for step in steps.LABELS:
        at.session_state["step"] = step
        at.run()
        assert not at.exception, step
        assert not at.error, step
        assert step in [b.label.split("  ")[0] for b in at.sidebar.button], "the menu lists every screen"


def test_environment_popover_starts_r_only_when_asked(fake_rscript, monkeypatch):
    """The sidebar's Environment popover must never start R on page load, on any page."""
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    from grnbench.ui import steps
    monkeypatch.setenv("GRNBENCH_RSCRIPT", str(fake_rscript.path))
    at = AppTest.from_file(str(ROOT / "src" / "grnbench" / "ui" / "app.py"), default_timeout=60)
    at.run()
    for step in steps.LABELS:
        at.session_state["step"] = step
        at.run()
        assert not at.exception, step
    assert fake_rscript.calls() == [], "R was started without the button"
    at.button(key="env-check").click().run()
    assert not at.exception
    assert [c["script"] for c in fake_rscript.calls()] == ["_probe.R"]
    rows = at.session_state["env_report"]
    assert rows[0]["tool"] == "R" and rows[0]["found"] and rows[0]["version"] == "4.5.3"
    codes = [c.value for c in at.code]
    assert any('BiocManager::install("scDblFinder")' in c for c in codes), "a missing package shows its install command"
    at.run()
    assert len(fake_rscript.calls()) == 1, "the stored report is shown again without starting R"
    at.button(key="env-check").click().run()
    assert len(fake_rscript.calls()) == 2, "pressing Check R again asks R again, not the cache"
