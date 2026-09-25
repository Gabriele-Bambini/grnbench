"""Plugin contract tests: Python objects, entry points, external commands."""
from __future__ import annotations

import importlib.util
import json
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from grnbench import plugins
from grnbench.infer import base as infer_base
from grnbench.infer import list_methods, run_method
from grnbench.plugins import (
    PluginError,
    canonical_edge_frame,
    command_method_spec,
    load_plugin,
    python_method_spec,
    register_plugins,
)

HAS_PYARROW = importlib.util.find_spec("pyarrow") is not None


@pytest.fixture(autouse=True)
def clean_registry():
    """Every test starts from the built-in registry and leaves it unchanged."""
    saved = dict(infer_base.REGISTRY)
    saved_flag = plugins._loaded_entry_points
    saved_refs = dict(plugins._registered_refs)
    yield
    infer_base.REGISTRY.clear()
    infer_base.REGISTRY.update(saved)
    plugins._loaded_entry_points = saved_flag
    plugins._registered_refs.clear()
    plugins._registered_refs.update(saved_refs)


# --------------------------------------------------------------------------- #
# python plugins
# --------------------------------------------------------------------------- #

class DemoModel:
    """First gene regulates everything (deterministic, contract-shaped)."""

    name = "demo-model"
    description = "toy plugin used in tests"

    def infer(self, adata, tf_list=None, **params):
        genes = [str(g) for g in adata.var_names]
        regs = [g for g in (tf_list or genes[:1]) if g in genes]
        rows = [(r, t, 1.0 / (1 + i))
                for r in regs for i, t in enumerate(genes)]
        return pd.DataFrame(rows, columns=["source", "target", "weight"])


def demo_function(adata, tf_list=None, **params):
    genes = [str(g) for g in adata.var_names]
    return pd.DataFrame({"source": [genes[0]] * len(genes),
                         "target": genes,
                         "weight": np.linspace(1.0, 0.1, len(genes))})


def test_class_plugin_runs_through_the_normal_pipeline(cascade_adata):
    spec = python_method_spec(DemoModel)
    assert spec.name == "demo-model" and spec.takes == "adata"
    infer_base.register_method(spec)

    edges, prov = run_method("demo-model", cascade_adata, tf_list=["TF_A"])
    assert list(edges.columns[:3]) == ["source", "target", "weight"]
    assert (edges["source"] == "TF_A").all()
    assert not (edges["source"] == edges["target"]).any(), "self-loops must be dropped"
    assert edges["weight"].is_monotonic_decreasing
    assert prov["origin"] == "python"
    assert prov["normalization"].startswith("raw counts detected")


def test_function_plugin_and_explicit_name(cascade_adata):
    spec = python_method_spec(demo_function, name="fn-model")
    infer_base.register_method(spec)
    edges, prov = run_method("fn-model", cascade_adata)
    assert len(edges) > 0 and prov["method"] == "fn-model"


def test_plugin_receives_normalized_adata_with_metadata(cascade_adata):
    seen = {}

    def spy(adata, tf_list=None, **params):
        seen["max"] = float(np.asarray(adata.X).max())
        seen["n_obs"] = adata.n_obs
        seen["has_counts_layer"] = "counts" in adata.layers
        return demo_function(adata, tf_list)

    infer_base.register_method(python_method_spec(spy, name="spy"))
    run_method("spy", cascade_adata)
    # cpm_log1p was applied on a copy -> log-scale values, raw counts still reachable
    assert seen["max"] < 12.0
    assert seen["n_obs"] == cascade_adata.n_obs
    assert seen["has_counts_layer"]


def test_plugin_cannot_silently_shadow_a_builtin():
    spec = python_method_spec(demo_function, name="spearman_signed")
    with pytest.raises(infer_base.MethodError, match="already registered"):
        infer_base.register_method(spec)
    infer_base.register_method(spec, replace=True)  # explicit override is allowed


def test_bad_plugin_shapes_raise_plugin_error():
    class NoInfer:
        name = "x"

    with pytest.raises(PluginError, match="no .infer"):
        python_method_spec(NoInfer)
    with pytest.raises(PluginError, match="neither callable"):
        python_method_spec(object())
    with pytest.raises(PluginError, match="not a module path"):
        plugins.import_object("no_colon_here")
    with pytest.raises(PluginError, match="cannot import"):
        plugins.import_object("definitely_not_a_module_xyz:Thing")


def test_plugin_returning_garbage_is_rejected(cascade_adata):
    infer_base.register_method(python_method_spec(
        lambda adata, tf_list=None, **p: "not a dataframe", name="junk"))
    with pytest.raises(infer_base.MethodError, match="must return a DataFrame"):
        run_method("junk", cascade_adata)


# --------------------------------------------------------------------------- #
# entry points
# --------------------------------------------------------------------------- #

def test_entry_point_discovery(monkeypatch):
    class FakeEP:
        name = "ep-model"
        value = "tests.test_plugins:DemoModel"

        def load(self):
            return DemoModel

    class BrokenEP:
        name = "broken-model"
        value = "nope:Nope"

        def load(self):
            raise ImportError("no such module")

    monkeypatch.setattr("importlib.metadata.entry_points",
                        lambda **kw: [FakeEP(), BrokenEP()])
    plugins._loaded_entry_points = False
    plugins.ensure_plugins_loaded()

    table = {m["name"]: m for m in list_methods()}
    assert table["ep-model"]["available"] and table["ep-model"]["origin"] == "entry_point"
    # a plugin that fails to import stays visible with its reason
    assert not table["broken-model"]["available"]
    assert "no such module" in table["broken-model"]["unavailable_reason"]


# --------------------------------------------------------------------------- #
# external command plugins
# --------------------------------------------------------------------------- #

_WORKER = textwrap.dedent("""
    import argparse, sys
    import pandas as pd

    ap = argparse.ArgumentParser()
    ap.add_argument("--expr", required=True)
    ap.add_argument("--tfs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--params", required=True)
    ap.add_argument("--fmt", default="tsv")
    a = ap.parse_args()

    expr = pd.read_parquet(a.expr) if a.fmt == "parquet" else pd.read_csv(a.expr, sep="\\t", index_col=0)
    tfs = [t for t in pd.read_csv(a.tfs)["tf"].astype(str).tolist() if t in expr.columns]
    regs = tfs or list(expr.columns[:1])
    rows = [(r, t, float(abs(expr[r].corr(expr[t])))) for r in regs for t in expr.columns]
    out = pd.DataFrame(rows, columns=["regulator", "target", "MI"])
    out.to_parquet(a.out) if a.fmt == "parquet" else out.to_csv(a.out, sep="\\t", index=False)
""")


def _worker_spec(tmp_path: Path, fmt: str) -> Path:
    script = tmp_path / "worker.py"
    script.write_text(_WORKER, encoding="utf-8")
    spec = {
        "name": f"cmd-model-{fmt}",
        "description": "external command plugin under test",
        "format": fmt,
        "timeout_s": 300,
        "command": [sys.executable, str(script),
                    "--expr", "{expr}", "--tfs", "{tfs}", "--out", "{out}",
                    "--params", "{params}", "--fmt", fmt],
    }
    p = tmp_path / f"plugin_{fmt}.json"
    p.write_text(json.dumps(spec), encoding="utf-8")
    return p


def test_command_plugin_tsv_roundtrip(tmp_path, cascade_adata):
    registered = register_plugins([str(_worker_spec(tmp_path, "tsv"))])
    assert registered == [{"name": "cmd-model-tsv", "origin": "command",
                           "requested_as": str(_worker_spec(tmp_path, "tsv"))}]

    edges, prov = run_method("cmd-model-tsv", cascade_adata, tf_list=["TF_A"])
    assert list(edges.columns[:3]) == ["source", "target", "weight"]
    assert (edges["source"] == "TF_A").all()
    assert prov["origin"] == "command"
    # the worker returns |corr|, so the planted cascade must lead the ranking
    assert set(edges.head(2)["target"]) & {"G_B", "G_C"}


@pytest.mark.skipif(not HAS_PYARROW, reason="pyarrow not installed")
def test_command_plugin_parquet_roundtrip(tmp_path, cascade_adata):
    register_plugins([str(_worker_spec(tmp_path, "parquet"))])
    edges, _ = run_method("cmd-model-parquet", cascade_adata, tf_list=["TF_A"])
    assert len(edges) > 0 and edges["weight"].is_monotonic_decreasing


def test_command_plugin_failure_is_loud(tmp_path, cascade_adata):
    script = tmp_path / "boom.py"
    script.write_text("import sys; sys.stderr.write('kaboom\\n'); sys.exit(3)",
                      encoding="utf-8")
    spec = tmp_path / "boom.json"
    spec.write_text(json.dumps({
        "name": "boom", "format": "tsv",
        "command": [sys.executable, str(script)],
    }), encoding="utf-8")
    register_plugins([str(spec)])
    with pytest.raises(PluginError, match="rc=3"):
        run_method("boom", cascade_adata)


def test_command_plugin_writing_nothing_is_loud(tmp_path, cascade_adata):
    script = tmp_path / "quiet.py"
    script.write_text("pass", encoding="utf-8")
    spec = tmp_path / "quiet.json"
    spec.write_text(json.dumps({
        "name": "quiet", "format": "tsv",
        "command": [sys.executable, str(script)],
    }), encoding="utf-8")
    register_plugins([str(spec)])
    with pytest.raises(PluginError, match="wrote no edges"):
        run_method("quiet", cascade_adata)


def test_missing_executable_is_unavailable_not_a_crash():
    spec = command_method_spec({
        "name": "ghost", "command": ["definitely-not-on-path-xyz", "--go"],
        "format": "tsv",
    })
    infer_base.register_method(spec)
    entry = next(m for m in list_methods() if m["name"] == "ghost")
    assert not entry["available"]
    assert "not on PATH" in entry["unavailable_reason"]
    with pytest.raises(infer_base.MethodError, match="unavailable"):
        run_method("ghost", None)


def test_docker_style_placeholders_are_substituted(tmp_path, cascade_adata):
    """`{workdir}` + `{expr_name}` is the documented Docker mount pattern."""
    script = tmp_path / "echo_args.py"
    script.write_text(textwrap.dedent("""
        import sys, pathlib, pandas as pd
        args = sys.argv[1:]
        workdir = pathlib.Path(args[0])
        assert args[1] == "expr.tsv", args
        assert (workdir / args[1]).exists()
        pd.DataFrame({"source": ["A"], "target": ["B"], "weight": [1.0]}).to_csv(
            workdir / "edges.tsv", sep="\\t", index=False)
    """), encoding="utf-8")
    spec = tmp_path / "docker_like.json"
    spec.write_text(json.dumps({
        "name": "docker-like", "format": "tsv",
        "command": [sys.executable, str(script), "{workdir}", "{expr_name}"],
    }), encoding="utf-8")
    register_plugins([str(spec)])
    edges, _ = run_method("docker-like", cascade_adata)
    assert edges.iloc[0]["source"] == "A"


# --------------------------------------------------------------------------- #
# spec files + edge-frame canonicalization
# --------------------------------------------------------------------------- #

def test_spec_file_with_python_key(tmp_path, cascade_adata):
    spec = tmp_path / "py.json"
    spec.write_text(json.dumps({
        "name": "via-file", "python": "test_plugins:DemoModel",
    }), encoding="utf-8")
    register_plugins([str(spec)])
    edges, prov = run_method("via-file", cascade_adata)
    assert prov["origin"] == "python" and len(edges) > 0


def test_yaml_spec_files_work_when_pyyaml_is_available(tmp_path, cascade_adata):
    """YAML specs are documented, so pyyaml must be installable via an extra."""
    pytest.importorskip("yaml")
    spec = tmp_path / "plug.yaml"
    spec.write_text(
        "name: yaml-model\n"
        "description: loaded from YAML\n"
        "python: test_plugins:DemoModel\n",
        encoding="utf-8")
    register_plugins([str(spec)])
    edges, prov = run_method("yaml-model", cascade_adata, tf_list=["TF_A"])
    assert len(edges) > 0 and prov["origin"] == "python"


def test_spec_file_can_hold_several_methods(tmp_path):
    spec = tmp_path / "many.json"
    spec.write_text(json.dumps({"methods": [
        {"name": "m1", "command": ["echo"], "format": "tsv"},
        {"name": "m2", "command": ["echo"], "format": "tsv"},
    ]}), encoding="utf-8")
    names = [s.name for s in plugins.load_plugin_file(spec)]
    assert names == ["m1", "m2"]


def test_bad_spec_files_fail_loudly(tmp_path):
    bad = tmp_path / "bad.json"
    bad.write_text("{not json", encoding="utf-8")
    with pytest.raises(PluginError, match="invalid JSON"):
        plugins.load_plugin_file(bad)

    nameless = tmp_path / "nameless.json"
    nameless.write_text(json.dumps({"command": ["echo"]}), encoding="utf-8")
    with pytest.raises(PluginError, match="needs a 'name'"):
        plugins.load_plugin_file(nameless)

    cmdless = tmp_path / "cmdless.json"
    cmdless.write_text(json.dumps({"name": "x"}), encoding="utf-8")
    with pytest.raises(PluginError, match="needs a 'command'"):
        plugins.load_plugin_file(cmdless)

    with pytest.raises(PluginError, match="not an existing spec file"):
        load_plugin("./nope.json")


def test_registering_the_same_plugin_twice_is_a_no_op(tmp_path):
    """`grnbench run` passes --plugin to both the infer and the eval stage."""
    spec = str(_worker_spec(tmp_path, "tsv"))
    first = register_plugins([spec])
    second = register_plugins([spec])          # must not raise
    assert first == second
    assert len([m for m in list_methods() if m["name"] == "cmd-model-tsv"]) == 1


def test_a_different_plugin_may_not_steal_a_taken_name(tmp_path):
    a = tmp_path / "a.json"
    a.write_text(json.dumps({"name": "clash", "command": ["echo"],
                             "format": "tsv"}), encoding="utf-8")
    b = tmp_path / "b.json"
    b.write_text(json.dumps({"name": "clash", "command": ["echo"],
                             "format": "tsv"}), encoding="utf-8")
    register_plugins([str(a)])
    with pytest.raises(PluginError, match="already registered"):
        register_plugins([str(b)])


def test_file_path_object_reference(tmp_path, cascade_adata):
    """A one-off model needs no packaging: point --plugin straight at the file."""
    mod = tmp_path / "one_off.py"
    mod.write_text(textwrap.dedent("""
        import numpy as np, pandas as pd
        class OneOff:
            name = "one-off"
            def infer(self, adata, tf_list=None, **p):
                genes = [str(g) for g in adata.var_names]
                return pd.DataFrame({"source": [genes[0]] * len(genes),
                                     "target": genes,
                                     "weight": np.linspace(1, 0.1, len(genes))})
    """), encoding="utf-8")
    register_plugins([f"{mod}:OneOff"])
    edges, prov = run_method("one-off", cascade_adata)
    assert len(edges) > 0 and prov["origin"] == "python"


def test_split_object_ref_tolerates_windows_drive_letters():
    assert plugins.split_object_ref("pkg.mod:Cls") == ("pkg.mod", "Cls")
    assert plugins.split_object_ref(r"C:\models\x.py:Cls") == (r"C:\models\x.py", "Cls")
    # a bare Windows path has no attribute: must not be read as head:attr
    assert plugins.split_object_ref(r"C:\models\x.py") == (None, None)
    assert plugins.split_object_ref("plain.json") == (None, None)


@pytest.mark.skipif(importlib.util.find_spec("sklearn") is None,
                    reason="scikit-learn not installed")
def test_shipped_ridge_example_plugin_works(cascade_adata):
    spec_file = Path(__file__).resolve().parent.parent / "examples" / "plugins" / "ridge.json"
    assert spec_file.exists(), "the documented example plugin must ship with the package"
    register_plugins([str(spec_file)])
    edges, prov = run_method("ridge", cascade_adata, tf_list=["TF_A"])
    assert list(edges.columns[:3]) == ["source", "target", "weight"]
    assert "direction" in edges.columns
    assert (edges["source"] == "TF_A").all()
    assert edges["weight"].is_monotonic_decreasing
    assert prov["origin"] == "python"
    # the planted cascade TF_A -> G_B must outrank the pure-noise genes
    assert "G_B" in set(edges.head(3)["target"])


def test_canonical_edge_frame_alias_matching_and_errors():
    df = pd.DataFrame({"Regulator": ["A"], "Target": ["B"], "MI": [0.5],
                       "Sign": [-1]})
    out = canonical_edge_frame(df)
    assert list(out.columns) == ["source", "target", "weight", "direction"]
    assert out.iloc[0]["direction"] == -1

    positional = pd.DataFrame({"x": ["A"], "y": ["B"], "z": [0.5]})
    assert canonical_edge_frame(positional).iloc[0]["weight"] == 0.5

    # an exact canonical name beats an alias, whatever the column order:
    # resolving "source" by position would silently pick the wrong regulator
    ambiguous = pd.DataFrame({"regulator": ["ALIAS"], "source": ["EXACT"],
                              "target": ["B"], "score": [1.0], "weight": [2.0]})
    resolved = canonical_edge_frame(ambiguous)
    assert resolved.iloc[0]["source"] == "EXACT"
    assert resolved.iloc[0]["weight"] == 2.0

    with pytest.raises(PluginError, match=">= 3 columns"):
        canonical_edge_frame(pd.DataFrame({"a": [1], "b": [2]}))
    with pytest.raises(PluginError, match="non-numeric weight"):
        canonical_edge_frame(pd.DataFrame({"source": ["A"], "target": ["B"],
                                           "weight": ["not-a-number"]}))
