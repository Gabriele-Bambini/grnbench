"""The optional R bridge: finding an Rscript, probing it once, running R scripts, reporting the environment.

Every test except the real-R ones runs a fake Rscript (tests/conftest.py), so no R is needed; the real-R
tests skip when no R is found.
"""
from __future__ import annotations

import importlib.resources
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

from grnbench import kb, rbridge
from grnbench.rbridge import RBridgeError, RRuntime

WINDOWS = os.name == "nt"


@pytest.fixture()
def bare(tmp_path, monkeypatch):
    """No R anywhere: no GRNBENCH_RSCRIPT or GRNBENCH_R_ENV, no conda, no R_HOME, nothing on PATH, a fresh cache."""
    for name in ("GRNBENCH_RSCRIPT", "GRNBENCH_R_ENV", "R_HOME", "CONDA_EXE", "R_LIBS", "R_LIBS_SITE", "R_LIBS_USER"):
        monkeypatch.delenv(name, raising=False)
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("USERPROFILE", str(home))
    system = [str(Path(os.environ.get("SystemRoot", r"C:\Windows")) / "System32")] if WINDOWS else []
    monkeypatch.setenv("PATH", os.pathsep.join(system))
    monkeypatch.setenv("GRNBENCH_CACHE_DIR", str(tmp_path / "cache"))
    if WINDOWS:     # the fakes are .cmd launchers; real installs have Rscript.exe
        monkeypatch.setattr(rbridge, "_EXE", ".cmd")
    rbridge.clear_memo()
    yield home
    rbridge.clear_memo()


def _conda_env(fake, root: Path, name: str, *, rscript: bool = True) -> Path:
    """A conda-like env folder (conda-meta, R library) with the fake as its Rscript."""
    prefix = root / "envs" / name
    if rscript:
        fake.launcher(prefix / ("Scripts" if WINDOWS else "bin"))
    (prefix / "lib" / "R" / "library").mkdir(parents=True)
    (prefix / "conda-meta").mkdir()
    return prefix


def _list_envs(home: Path, *prefixes: Path) -> None:
    (home / ".conda").mkdir(exist_ok=True)
    (home / ".conda" / "environments.txt").write_text("".join(f"{p}\n" for p in prefixes), encoding="utf-8")


def _script(path: Path, posix: str, windows: str) -> Path:
    """A tiny executable script (``.cmd`` on Windows)."""
    path = path.with_suffix(".cmd") if WINDOWS else path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(windows + "\r\n" if WINDOWS else "#!/bin/sh\n" + posix + "\n", encoding="utf-8")
    path.chmod(0o755)
    return path


def _on_path(monkeypatch, folder: Path) -> None:
    monkeypatch.setenv("PATH", str(folder) + os.pathsep + os.environ["PATH"])


# --------------------------------------------------------------------------- #
# detection
# --------------------------------------------------------------------------- #

def test_without_any_r_nothing_is_found_and_nothing_is_started(bare, fake_rscript):
    assert rbridge.find_rscript() is None
    assert fake_rscript.calls() == []


def test_grnbench_rscript_none_turns_r_off_without_starting_anything(bare, fake_rscript, monkeypatch):
    _on_path(monkeypatch, fake_rscript.path.parent)
    monkeypatch.setenv("GRNBENCH_RSCRIPT", "")             # empty means unset
    assert rbridge.find_rscript() is not None               # R on PATH is found ...
    rbridge.clear_memo()
    before = len(fake_rscript.calls())
    for value in ("none", "NONE", " None "):
        monkeypatch.setenv("GRNBENCH_RSCRIPT", value)      # ... until it is turned off
        assert rbridge.find_rscript() is None
    assert len(fake_rscript.calls()) == before


def test_detection_order_explicit_env_conda_r_home_path(bare, fake_rscript, tmp_path, monkeypatch):
    on_path = fake_rscript.launcher(tmp_path / "on_path")
    _on_path(monkeypatch, on_path.parent)
    r_home = tmp_path / "r_home"
    fake_rscript.launcher(r_home / "bin")
    monkeypatch.setenv("R_HOME", str(r_home))
    prefix = _conda_env(fake_rscript, tmp_path, "grnbench-r-test")
    _list_envs(bare, _conda_env(fake_rscript, tmp_path, "something-else"), prefix)
    monkeypatch.setenv("GRNBENCH_R_ENV", "grnbench-r-test")
    monkeypatch.setenv("GRNBENCH_RSCRIPT", str(fake_rscript.launcher(tmp_path / "env_var")))
    explicit = fake_rscript.launcher(tmp_path / "explicit")

    assert rbridge.find_rscript(str(explicit)).source == "--rscript"
    assert rbridge.find_rscript().source == "GRNBENCH_RSCRIPT"
    monkeypatch.delenv("GRNBENCH_RSCRIPT")
    rt = rbridge.find_rscript()
    assert rt.source == "conda:grnbench-r-test" and Path(rt.rscript).parent.parent == prefix
    assert rt.library == str(prefix / "lib" / "R" / "library") and rt.conda_prefix == str(prefix)
    assert rt.path_prepend[0] == str(prefix) and rt.version == "4.5.3"
    if WINDOWS:
        assert [Path(p).relative_to(prefix).as_posix() for p in rt.path_prepend[1:]] == [
            "Library/mingw-w64/bin", "Library/usr/bin", "Library/bin", "Scripts", "bin"]
    else:
        assert rt.path_prepend == (str(prefix / "bin"),)
    monkeypatch.delenv("GRNBENCH_R_ENV")                  # the default env grnbench-r does not exist: skipped
    assert rbridge.find_rscript().source == "R_HOME"
    monkeypatch.delenv("R_HOME")
    rt = rbridge.find_rscript()
    assert rt.source == "PATH" and Path(rt.rscript) == on_path and rt.conda_prefix is None
    assert rt.library == str(fake_rscript.library), "an R that is not a conda env reports its own library"


def test_a_named_env_that_is_missing_is_an_error_not_a_switch(bare, fake_rscript, monkeypatch):
    _on_path(monkeypatch, fake_rscript.path.parent)
    monkeypatch.setenv("GRNBENCH_R_ENV", "no-such-env")
    with pytest.raises(RBridgeError, match="GRNBENCH_R_ENV=no-such-env: no conda env called no-such-env") as info:
        rbridge.find_rscript()
    assert "environments.txt" in str(info.value) and "conda" in str(info.value)
    assert fake_rscript.calls() == [], "the R on PATH was not used instead"


def test_an_existing_env_without_a_working_r_is_an_error_even_with_r_on_path(bare, fake_rscript, tmp_path,
                                                                                monkeypatch):
    _on_path(monkeypatch, fake_rscript.path.parent)
    empty = _conda_env(fake_rscript, tmp_path, "grnbench-r", rscript=False)       # the default name
    _list_envs(bare, empty)
    with pytest.raises(RBridgeError, match="conda:grnbench-r: no Rscript at"):
        rbridge.find_rscript()
    broken = tmp_path / "envs2" / "grnbench-r"
    _script(broken / ("Scripts" if WINDOWS else "bin") / "Rscript", "exit 1", "@exit /b 1")
    (broken / "conda-meta").mkdir()
    _list_envs(bare, broken)
    with pytest.raises(RBridgeError, match="conda:grnbench-r: .* does not run: exited with code 1"):
        rbridge.find_rscript()
    assert fake_rscript.calls() == []


def test_an_explicit_or_env_rscript_that_does_not_work_is_an_error(bare, fake_rscript, tmp_path, monkeypatch):
    _on_path(monkeypatch, fake_rscript.path.parent)
    with pytest.raises(RBridgeError, match="no Rscript at"):
        rbridge.find_rscript(str(tmp_path / "nowhere" / "Rscript"))
    with pytest.raises(RBridgeError, match="--rscript is empty"):
        rbridge.find_rscript("")
    broken = _script(tmp_path / "broken", "exit 1", "@exit /b 1")
    with pytest.raises(RBridgeError, match="code 1"):
        rbridge.find_rscript(str(broken))
    monkeypatch.setenv("GRNBENCH_RSCRIPT", str(broken))
    with pytest.raises(RBridgeError, match="GRNBENCH_RSCRIPT"):
        rbridge.find_rscript()                          # never a silent fall-through to the R on PATH


def test_relative_rscript_paths_become_absolute(bare, fake_rscript, tmp_path, monkeypatch):
    monkeypatch.chdir(fake_rscript.path.parent.parent)
    relative = os.path.join("bin", fake_rscript.path.name)
    rt = rbridge.find_rscript(relative)
    assert Path(rt.rscript).is_absolute() and Path(rt.rscript) == fake_rscript.path
    monkeypatch.setenv("GRNBENCH_RSCRIPT", relative)
    assert Path(rbridge.find_rscript().rscript) == fake_rscript.path


def test_a_candidate_that_fails_is_skipped_and_the_reason_kept(bare, fake_rscript, tmp_path, monkeypatch):
    monkeypatch.setenv("R_HOME", str(tmp_path / "empty_r_home"))
    _on_path(monkeypatch, fake_rscript.path.parent)
    rt, notes = rbridge._find()
    assert rt.source == "PATH"
    assert any(n.startswith("R_HOME: no Rscript at") for n in notes)
    assert any(n.startswith("conda:grnbench-r: no conda env called grnbench-r") for n in notes)
    r_row = rbridge.environment_report()[0]
    assert r_row["tool"] == "R" and r_row["found"] and "# skipped: R_HOME: no Rscript at" in r_row["install"]
    hint = rbridge.requirement_status("r:scDblFinder").install_hint
    assert "# skipped: R_HOME: no Rscript at" in hint, "the hint says another R was passed over"


def _fake_conda(tmp_path: Path, monkeypatch, body: str) -> Path:
    """CONDA_EXE pointing at a Python script; every call adds a line to the returned file."""
    asked = tmp_path / "conda_calls.txt"
    code = tmp_path / "conda_fake.py"
    code.write_text(f"import json, sys\nopen({str(asked)!r}, 'a').write('x\\n')\n{body}\n", encoding="utf-8")
    if WINDOWS:
        conda = tmp_path / "conda.cmd"
        conda.write_text(f'@"{sys.executable}" "{code}" %*\r\n', encoding="utf-8")
    else:
        conda = tmp_path / "conda"
        conda.write_text(f"#!{sys.executable}\n" + code.read_text(encoding="utf-8"), encoding="utf-8")
        conda.chmod(0o755)
    monkeypatch.setenv("CONDA_EXE", str(conda))
    return asked


def test_conda_env_list_is_asked_only_when_environments_txt_has_no_match(bare, fake_rscript, tmp_path, monkeypatch):
    prefix = _conda_env(fake_rscript, tmp_path, "grnbench-r")
    asked = _fake_conda(tmp_path, monkeypatch, f"print(json.dumps({{'envs': [{str(prefix)!r}]}}))")
    assert rbridge.find_rscript().source == "conda:grnbench-r"
    rbridge.clear_memo()
    assert rbridge.find_rscript().source == "conda:grnbench-r"
    assert asked.read_text().count("x") == 1, "conda env list runs at most once per cache"
    assert rbridge.find_rscript(fresh=True).source == "conda:grnbench-r"
    assert asked.read_text().count("x") == 2, "fresh=True asks conda again"


def test_a_failed_conda_env_list_is_never_cached_and_says_why(bare, fake_rscript, tmp_path, monkeypatch):
    asked = _fake_conda(tmp_path, monkeypatch, "print('CondaError: broken', file=sys.stderr); sys.exit(1)")
    _on_path(monkeypatch, fake_rscript.path.parent)
    for n in (1, 2):
        rt, notes = rbridge._find()
        assert rt.source == "PATH"
        assert any("env list --json` exited with code 1: CondaError: broken" in note for note in notes), notes
        assert asked.read_text().count("x") == n
    monkeypatch.setenv("GRNBENCH_R_ENV", "grnbench-r")
    with pytest.raises(RBridgeError, match="CondaError: broken"):
        rbridge.find_rscript()


# --------------------------------------------------------------------------- #
# probe cost and validity
# --------------------------------------------------------------------------- #

def test_the_probe_runs_once_per_cache_and_again_after_an_install(bare, fake_rscript, monkeypatch):
    monkeypatch.setenv("GRNBENCH_RSCRIPT", str(fake_rscript.path))
    rt = rbridge.find_rscript()
    rbridge.probe(rt, rbridge.R_REQUIREMENTS)
    assert rbridge.find_rscript() == rt
    rbridge.clear_memo()                                # a new worker process: the disk cache answers
    assert rbridge.find_rscript() == rt
    assert len(fake_rscript.calls()) == 1
    stamp = time.time() + 10                            # installing a package changes the library folder
    os.utime(fake_rscript.library, (stamp, stamp))
    rbridge.find_rscript()
    assert len(fake_rscript.calls()) == 2
    for f in Path(os.environ["GRNBENCH_CACHE_DIR"]).glob("*.json"):
        f.write_text("{ not json", encoding="utf-8")   # a corrupt cache is a miss, never an error
    rbridge.clear_memo()
    assert rbridge.find_rscript() == rt
    assert len(fake_rscript.calls()) == 3
    assert rbridge.find_rscript(fresh=True) == rt and len(fake_rscript.calls()) == 4, "fresh=True asks R again"
    rbridge.probe(rt, rbridge.R_REQUIREMENTS, fresh=True)
    assert len(fake_rscript.calls()) == 5


def test_touching_any_library_folder_invalidates_the_probe(bare, fake_rscript, tmp_path, monkeypatch):
    second = tmp_path / "site-library"
    second.mkdir()
    monkeypatch.setenv("FAKE_R_LIBRARY", os.pathsep.join([str(fake_rscript.library), str(second)]))
    monkeypatch.setenv("GRNBENCH_RSCRIPT", str(fake_rscript.path))
    rbridge.find_rscript()
    rbridge.find_rscript()
    assert len(fake_rscript.calls()) == 1
    stamp = time.time() + 10
    os.utime(second, (stamp, stamp))
    rbridge.find_rscript()
    assert len(fake_rscript.calls()) == 2


def test_probe_cache_depends_on_the_callers_r_library_settings_outside_conda(bare, fake_rscript, tmp_path,
                                                                             monkeypatch):
    monkeypatch.setenv("GRNBENCH_RSCRIPT", str(fake_rscript.path))
    rbridge.find_rscript()
    monkeypatch.setenv("R_LIBS_USER", str(tmp_path))      # another R library for a system R: ask again
    rbridge.find_rscript()
    assert len(fake_rscript.calls()) == 2


def test_probe_reports_versions_and_missing_packages(bare, fake_rscript, monkeypatch):
    monkeypatch.setenv("GRNBENCH_RSCRIPT", str(fake_rscript.path))
    info = rbridge.probe(rbridge.find_rscript(), ["SoupX", "celda"])
    assert info["version"] == "4.5.3" and info["packages"] == {"SoupX": "1.6.2", "celda": None}
    assert info["library"] == [str(fake_rscript.library)] and info["built"]["SoupX"] == "4.5.3"


def test_a_package_built_for_another_r_is_not_usable(bare, fake_rscript, monkeypatch):
    monkeypatch.setenv("GRNBENCH_RSCRIPT", str(fake_rscript.path))
    monkeypatch.setenv("FAKE_R_BUILT", json.dumps({"SoupX": "4.4.1"}))
    st = rbridge.requirement_status("r:SoupX")
    assert not st.ok and st.reason == "built for R 4.4" and "reinstall" in st.install_hint
    row = next(r for r in rbridge.environment_report() if r["tool"] == "r:SoupX")
    assert not row["found"] and row["version"] == "1.6.2" and "built for R 4.4" in row["install"]
    assert rbridge.requirement_status("r:jsonlite").ok, "built for 4.5.x on R 4.5.3 is fine"


def test_an_unwritable_cache_folder_is_not_an_error(bare, fake_rscript, tmp_path, monkeypatch):
    blocker = tmp_path / "a_file"
    blocker.write_text("x", encoding="utf-8")
    monkeypatch.setenv("GRNBENCH_CACHE_DIR", str(blocker / "cache"))
    monkeypatch.setenv("GRNBENCH_RSCRIPT", str(fake_rscript.path))
    assert rbridge.find_rscript().version == "4.5.3"


def test_r_requirements_come_from_the_knowledge_base():
    from_kb = {r.split(":", 1)[1] for d in kb.DECISIONS for o in d.options for r in o.requires if r.startswith("r:")}
    assert rbridge.R_REQUIREMENTS == tuple(sorted(from_kb | {"jsonlite"}))
    assert {"SoupX", "scDblFinder", "jsonlite"} <= set(rbridge.R_REQUIREMENTS)


# --------------------------------------------------------------------------- #
# requirement status and install hints
# --------------------------------------------------------------------------- #

def test_requirement_status_for_r_python_and_command_line_tools(bare, fake_rscript, monkeypatch):
    monkeypatch.setenv("GRNBENCH_RSCRIPT", str(fake_rscript.path))
    rt = rbridge.find_rscript()
    missing = rbridge.requirement_status("r:scDblFinder", rt)
    assert not missing.ok and missing.reason == "missing: r:scDblFinder"
    assert 'BiocManager::install("scDblFinder")' in missing.install_hint and "4.5.3" in missing.install_hint
    assert "GRNBENCH_RSCRIPT" in missing.install_hint
    assert rbridge.requirement_status("r:SoupX", rt).ok
    assert rbridge.requirement_status("r:SoupX").ok, "the default finds the R itself"
    no_r = rbridge.requirement_status("r:SoupX", None)
    assert no_r.reason == "missing: r:SoupX" and "docs/R_BRIDGE.md" in no_r.install_hint

    assert rbridge.requirement_status("py:json").ok
    py = rbridge.requirement_status("py:no_such_module_xyz")
    assert py.reason == "missing: py:no_such_module_xyz" and py.install_hint == "pip install no_such_module_xyz"
    assert "scvi-tools" in rbridge.requirement_status("py:scvi").install_hint or rbridge.requirement_status("py:scvi").ok

    _on_path(monkeypatch, fake_rscript.path.parent)
    assert rbridge.requirement_status("cli:Rscript").ok
    cli = rbridge.requirement_status("cli:no-such-tool-xyz")
    assert cli.reason == "missing: cli:no-such-tool-xyz" and "PATH" in cli.install_hint


def test_requirement_status_says_when_r_is_turned_off(bare, monkeypatch):
    monkeypatch.setenv("GRNBENCH_RSCRIPT", "none")
    st = rbridge.requirement_status("r:SoupX")
    assert not st.ok and st.reason == "missing: r:SoupX" and "GRNBENCH_RSCRIPT=none" in st.install_hint


def test_install_hints_name_the_right_repository_and_the_env_library():
    plain = RRuntime("Rscript", "4.5.3", "PATH", None)
    conda = RRuntime("Rscript", "4.5.3", "conda:grnbench-r", "/env/lib/R/library", ("/env/bin",), "/env")
    assert rbridge._r_install_hint("scDblFinder", plain).startswith('BiocManager::install("scDblFinder")  #')
    assert rbridge._r_install_hint("SoupX", plain).startswith('install.packages("SoupX")  #')
    assert rbridge._r_install_hint("sctransform", conda).startswith('install.packages("sctransform", lib = .Library)')
    assert rbridge._r_install_hint("celda", conda).startswith('BiocManager::install("celda", lib = .Library)')
    assert rbridge._r_install_hint("DoubletFinder", plain).startswith(
        'remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")')
    assert "not on CRAN or Bioconductor" in rbridge._r_install_hint("someRandomPkg", plain)


# --------------------------------------------------------------------------- #
# subprocess environment
# --------------------------------------------------------------------------- #

def test_subprocess_env_pins_the_library_only_for_a_conda_r():
    base = {"PATH": "orig", "R_LIBS": "a", "R_LIBS_SITE": "b", "R_HOME": "c", "R_LIBS_USER": "d", "KEEP": "1"}
    rt = RRuntime("Rscript", "4.5.3", "conda:grnbench-r", "/env/lib/R/library", ("/env", "/env/bin"), "/env")
    env = rt.subprocess_env(base)
    assert env["PATH"] == os.pathsep.join(["/env", "/env/bin", "orig"])
    assert env["R_LIBS_USER"] == "/env/lib/R/library" and env["KEEP"] == "1"
    assert not {"R_LIBS", "R_LIBS_SITE", "R_HOME"} & set(env)
    assert base["R_HOME"] == "c", "the caller's mapping is not changed"
    own = RRuntime("Rscript", "4.5.3", "PATH", "/usr/lib/R/library").subprocess_env(base)
    assert own["R_LIBS_USER"] == "d" and own["R_LIBS"] == "a" and own["R_LIBS_SITE"] == "b", "R's own settings stay"
    assert "R_HOME" not in own and own["PATH"] == "orig"
    kept = RRuntime("Rscript", "4.5.3", "R_HOME", None).subprocess_env(base)
    assert kept["R_HOME"] == "c"
    assert "PATH" in RRuntime("Rscript", "4.5.3", "PATH", None).subprocess_env()


@pytest.mark.skipif(not WINDOWS, reason="environment names are case-insensitive on Windows only")
def test_subprocess_env_merges_path_spellings_on_windows():
    env = RRuntime("Rscript", "4.5.3", "conda:x", None, ("C:\\env",), "C:\\env").subprocess_env(
        {"Path": "C:\\orig", "r_home": "x"})
    assert [k for k in env if k.upper() == "PATH"] == ["PATH"] and env["PATH"] == "C:\\env;C:\\orig"
    assert not [k for k in env if k.upper() == "R_HOME"]


def test_the_child_gets_vanilla_the_path_prefix_and_the_env_library(bare, fake_rscript, tmp_path, monkeypatch):
    prefix = _conda_env(fake_rscript, tmp_path, "grnbench-r")
    _list_envs(bare, prefix)
    monkeypatch.setenv("R_LIBS_USER", str(tmp_path / "personal"))
    rt = rbridge.find_rscript()
    call = fake_rscript.calls()[-1]
    assert call["vanilla"] and call["script"] == "_probe.R"
    assert call["path"][:len(rt.path_prepend)] == list(rt.path_prepend)
    assert call["r_libs_user"] == str(prefix / "lib" / "R" / "library")
    monkeypatch.setenv("GRNBENCH_RSCRIPT", str(fake_rscript.path))      # not a conda R: its own settings
    rbridge.find_rscript()
    assert fake_rscript.calls()[-1]["r_libs_user"] == str(tmp_path / "personal")


# --------------------------------------------------------------------------- #
# run_r / run_external
# --------------------------------------------------------------------------- #

@pytest.fixture()
def method(tmp_path) -> Path:
    script = tmp_path / "method.R"
    script.write_text("x <- runif(5)\n", encoding="utf-8")
    return script


def test_run_r_returns_outputs_and_records_seed_and_inputs(bare, fake_rscript, method, tmp_path, monkeypatch):
    monkeypatch.setenv("GRNBENCH_RSCRIPT", str(fake_rscript.path))
    counts = tmp_path / "counts.h5ad"
    counts.write_text("x", encoding="utf-8")
    work = tmp_path / "work"
    import numpy as np
    res = rbridge.run_r(method, inputs={"counts": counts}, seed=7, workdir=work, expect={"table": "out.tsv"},
                        params={"k": 3, "n": np.int64(4), "x": np.float32(0.5), "flag": np.bool_(True),
                                "v": np.arange(3), "where": tmp_path})
    assert res.outputs == {"table": work.resolve() / "out.tsv"} and res.outputs["table"].read_text() == "table seed=7\n"
    sent = json.loads((work / "params.json").read_text(encoding="utf-8"))
    assert sent["seed"] == 7 and sent["expect"] == {"table": "out.tsv"}
    assert sent["params"] == {"k": 3, "n": 4, "x": 0.5, "flag": True, "v": [0, 1, 2], "where": str(tmp_path)}
    assert sent["inputs"] == {"counts": str(counts.resolve())} and Path(sent["workdir"]) == work.resolve()
    assert res.session["seed"] == 7 and res.log == work.resolve() / "r.log" and res.log.exists()
    assert res.runtime.source == "GRNBENCH_RSCRIPT" and res.seconds > 0 and res.skipped == ()
    assert [c["script"] for c in fake_rscript.calls()] == ["_probe.R", "_wrapper.R"]
    call = fake_rscript.calls()[-1]       # only the wrapper on the command line; the rest is in params.json
    assert call["argv"] == ["--vanilla", str(rbridge.WRAPPER)] and Path(call["cwd"]) == work.resolve()
    assert sent["script"] == str(method.resolve())


def test_run_r_refuses_nan_in_params(bare, fake_rscript, method, tmp_path):
    rt = rbridge.find_rscript(str(fake_rscript.path))
    for bad in (float("nan"), float("inf"), object()):
        with pytest.raises(RBridgeError, match="params must be JSON values"):
            rbridge.run_r(method, params={"x": bad}, seed=1, workdir=tmp_path / "w", rt=rt)
    assert not (tmp_path / "w").exists()


def test_run_r_records_the_r_candidates_it_skipped(bare, fake_rscript, method, tmp_path, monkeypatch):
    monkeypatch.setenv("R_HOME", str(tmp_path / "empty_r_home"))
    _on_path(monkeypatch, fake_rscript.path.parent)
    res = rbridge.run_r(method, seed=1, workdir=tmp_path / "w", expect={"t": "o.tsv"})
    assert any(s.startswith("R_HOME: no Rscript at") for s in res.skipped)
    assert any(s.startswith("R_HOME: no Rscript at") for s in res.manifest_block()["skipped"])


def test_run_r_timeout_kills_the_whole_tree_and_names_the_log(bare, fake_rscript, method, tmp_path, monkeypatch):
    psutil = pytest.importorskip("psutil")
    monkeypatch.setenv("FAKE_R_MODE", "sleep")
    child_pid = tmp_path / "child.pid"
    monkeypatch.setenv("FAKE_R_CHILD_PID", str(child_pid))
    rt = rbridge.find_rscript(str(fake_rscript.path))
    started = time.monotonic()
    with pytest.raises(RBridgeError, match="R timed out after 3 s; log: ") as info:
        rbridge.run_r(method, seed=1, workdir=tmp_path / "w", timeout_s=3, rt=rt, expect={"t": "o.tsv"})
    assert time.monotonic() - started < 20, "the process tree is killed, not waited for"
    assert str((tmp_path / "w" / "r.log").resolve()) in str(info.value)
    assert "may still be running" not in str(info.value)
    pid = int(child_pid.read_text())

    def gone() -> bool:
        try:
            return psutil.Process(pid).status() == psutil.STATUS_ZOMBIE
        except psutil.NoSuchProcess:
            return True

    deadline = time.monotonic() + 10
    while not gone():
        assert time.monotonic() < deadline, f"R's own child (pid {pid}) outlived the timeout"
        time.sleep(0.2)


def test_an_unconfirmed_kill_says_the_process_may_still_run(bare, fake_rscript, method, tmp_path, monkeypatch):
    monkeypatch.setenv("FAKE_R_MODE", "sleep")
    real = rbridge._kill_tree
    monkeypatch.setattr(rbridge, "_kill_tree", lambda proc: real(proc) and False)
    rt = rbridge.find_rscript(str(fake_rscript.path))
    with pytest.raises(RBridgeError) as info:
        rbridge.run_r(method, seed=1, workdir=tmp_path / "w", timeout_s=1, rt=rt)
    assert str(info.value).endswith(")") and "; process may still be running (pid " in str(info.value)


def test_an_interrupt_kills_r_and_propagates(bare, fake_rscript, method, tmp_path, monkeypatch):
    monkeypatch.setenv("FAKE_R_MODE", "sleep")
    rt = rbridge.find_rscript(str(fake_rscript.path))
    killed = []
    real_kill = rbridge._kill_tree
    monkeypatch.setattr(rbridge, "_kill_tree", lambda proc: killed.append(proc) or real_kill(proc))
    real_wait = subprocess.Popen.wait
    state = {"raised": False}

    def wait(self, timeout=None):
        if not state["raised"]:
            state["raised"] = True
            raise KeyboardInterrupt
        return real_wait(self, timeout)

    monkeypatch.setattr(subprocess.Popen, "wait", wait)
    with pytest.raises(KeyboardInterrupt):
        rbridge.run_r(method, seed=1, workdir=tmp_path / "w", timeout_s=60, rt=rt)
    assert len(killed) == 1 and killed[0].poll() is not None


def test_run_r_failure_shows_the_exit_code_and_the_end_of_the_log(bare, fake_rscript, method, tmp_path, monkeypatch):
    monkeypatch.setenv("FAKE_R_MODE", "fail")
    rt = rbridge.find_rscript(str(fake_rscript.path))
    with pytest.raises(RBridgeError) as info:
        rbridge.run_r(method, seed=1, workdir=tmp_path / "w", rt=rt)
    msg = str(info.value)
    assert "code 3" in msg and "Error: boom" in msg and "log line 29" in msg and "log line 0" not in msg


def test_run_r_refuses_missing_outputs_and_a_bad_session(bare, fake_rscript, method, tmp_path, monkeypatch):
    rt = rbridge.find_rscript(str(fake_rscript.path))
    monkeypatch.setenv("FAKE_R_MODE", "no_output")
    with pytest.raises(RBridgeError, match="did not write.*table"):
        rbridge.run_r(method, seed=1, workdir=tmp_path / "a", rt=rt, expect={"table": "out.tsv"})
    monkeypatch.setenv("FAKE_R_MODE", "bad_session")
    with pytest.raises(RBridgeError, match="session.json"):
        rbridge.run_r(method, seed=1, workdir=tmp_path / "b", rt=rt, expect={"table": "out.tsv"})
    monkeypatch.setenv("FAKE_R_MODE", "ok")
    rbridge.run_r(method, seed=1, workdir=tmp_path / "c", rt=rt, expect={"table": "out.tsv"})
    monkeypatch.setenv("FAKE_R_MODE", "no_output")
    with pytest.raises(RBridgeError, match="did not write"):    # a stale output of an earlier run never counts
        rbridge.run_r(method, seed=1, workdir=tmp_path / "c", rt=rt, expect={"table": "out.tsv"})


def test_expected_outputs_stay_inside_the_work_folder(bare, fake_rscript, method, tmp_path):
    rt = rbridge.find_rscript(str(fake_rscript.path))
    work = tmp_path / "w"
    (work / "folder").mkdir(parents=True)
    victim = tmp_path / "precious.txt"
    victim.write_text("keep", encoding="utf-8")
    for bad, why in ((str(victim), "relative to the work folder"), ("../precious.txt", "outside the work folder"),
                     ("sub/../../precious.txt", "outside the work folder"), (".", "is the work folder itself"),
                     ("", "is the work folder itself"),
                     ("params.json", "bridge's own files"), ("r.log", "bridge's own files"),
                     ("session.json", "bridge's own files"), ("folder", "is a folder")):
        with pytest.raises(RBridgeError, match=why):
            rbridge.run_r(method, seed=1, workdir=work, rt=rt, expect={"out": bad})
    assert victim.read_text(encoding="utf-8") == "keep" and fake_rscript.calls()[-1]["script"] == "_probe.R"
    assert rbridge.run_r(method, seed=1, workdir=work, rt=rt, expect={"out": "sub/../ok.tsv"}).outputs["out"] == \
        work.resolve() / "ok.tsv"


def test_a_log_that_cannot_be_written_is_reported_as_such(bare, fake_rscript, method, tmp_path):
    rt = rbridge.find_rscript(str(fake_rscript.path))
    (tmp_path / "w" / "r.log").mkdir(parents=True)
    with pytest.raises(RBridgeError, match="cannot write the log"):
        rbridge.run_r(method, seed=1, workdir=tmp_path / "w", rt=rt)


def test_run_r_seed_must_be_a_plain_int_in_r_range(bare, fake_rscript, method, tmp_path):
    rt = rbridge.find_rscript(str(fake_rscript.path))
    for bad, error in ((True, TypeError), ("1", TypeError), (1.0, TypeError), (-1, ValueError), (2**31, ValueError)):
        with pytest.raises(error):
            rbridge.run_r(method, seed=bad, workdir=tmp_path / "w", rt=rt)
    assert not (tmp_path / "w").exists(), "a refused seed starts nothing and writes nothing"
    assert rbridge.run_r(method, seed=2**31 - 1, workdir=tmp_path / "w", rt=rt).session["seed"] == 2**31 - 1


def test_run_r_without_r_says_how_to_get_it(bare, method, tmp_path):
    with pytest.raises(RBridgeError, match="docs/R_BRIDGE.md"):
        rbridge.run_r(method, seed=1, workdir=tmp_path / "w")


def test_manifest_block_replaces_the_home_folder(bare, fake_rscript, method, tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("USERPROFILE", str(tmp_path))
    rt = rbridge.find_rscript(str(fake_rscript.path))
    block = rbridge.run_r(method, seed=5, workdir=tmp_path / "w", rt=rt, expect={"t": "o.tsv"}).manifest_block()
    assert set(block) == {"version", "source", "rscript", "library", "packages", "seed", "rng_kind",
                          "rng_kind_after", "set_seed_calls", "rng_changed", "seconds", "log", "skipped"}
    assert block["seed"] == 5 and block["version"] == "4.5.3" and block["packages"] == {"jsonlite": "2.0.0"}
    assert block["rng_changed"] is False and block["set_seed_calls"] == [] and block["skipped"] == []
    assert block["rscript"].startswith("~/") and block["library"].startswith("~/") and block["log"] == "~/w/r.log"
    text = json.dumps(block)
    assert str(tmp_path) not in text and str(tmp_path).replace("\\", "/") not in text


def test_error_messages_write_the_home_folder_as_tilde(bare, monkeypatch, tmp_path):
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("USERPROFILE", str(tmp_path))
    with pytest.raises(RBridgeError) as info:
        rbridge.find_rscript(str(tmp_path / "nowhere" / "Rscript"))
    assert str(tmp_path) not in str(info.value) and "~/nowhere/Rscript" in str(info.value)
    assert str(RBridgeError(f"see {tmp_path}{os.sep}x.log")) == f"see ~{os.sep}x.log"
    quoted = str(RBridgeError(f"no such file: {str(tmp_path / 'x.log')!r}"))      # repr doubles backslashes
    assert str(tmp_path) not in quoted and str(tmp_path).replace("\\", "\\\\") not in quoted and "~" in quoted


def test_a_refused_expect_path_does_not_print_the_home_folder(bare, fake_rscript, method, tmp_path):
    rt = rbridge.find_rscript(str(fake_rscript.path))
    with pytest.raises(RBridgeError, match="relative to the work folder") as info:
        rbridge.run_r(method, seed=1, workdir=tmp_path / "w", rt=rt, expect={"out": str(bare / "victim.txt")})
    assert str(bare) not in str(info.value) and str(bare).replace("\\", "\\\\") not in str(info.value)


def test_grnbench_env_exits_1_when_looking_for_r_is_an_error(bare, fake_rscript, tmp_path, monkeypatch, capsys):
    from grnbench.cli_main import main as grnbench_main
    assert grnbench_main(["env", "--json"]) == 0, "no R anywhere is an answer, not a failure"
    capsys.readouterr()
    monkeypatch.setenv("GRNBENCH_R_ENV", "no-such-env")
    assert grnbench_main(["env", "--json"]) == 1, "a named env that does not exist"
    assert "no conda env called no-such-env" in json.loads(capsys.readouterr().out)[0]["install"]
    monkeypatch.delenv("GRNBENCH_R_ENV")
    _list_envs(bare, _conda_env(fake_rscript, tmp_path, "grnbench-r", rscript=False))
    assert grnbench_main(["env", "--json"]) == 1, "the default env exists but has no R"
    assert "no Rscript at" in json.loads(capsys.readouterr().out)[0]["install"]
    assert fake_rscript.calls() == [], "no other R was tried instead"


def test_grnbench_env_starts_a_broken_r_once_and_never_crashes(bare, tmp_path, monkeypatch, capsys):
    from grnbench.cli_main import main as grnbench_main
    starts = tmp_path / "starts.txt"
    broken = _script(tmp_path / "broken", f'echo x >> "{starts}"; exit 1', f'@echo x>> "{starts}"\r\n@exit /b 1')
    assert grnbench_main(["env", "--rscript", str(broken)]) == 1
    assert len(starts.read_text().splitlines()) == 1, "a failing R is asked once, not again for the exit code"
    capsys.readouterr()

    def broken_kb():
        raise ZeroDivisionError("kb is broken")

    monkeypatch.setattr(rbridge, "_kb_uses", broken_kb)
    assert grnbench_main(["env", "--json"]) == 1, "an internal error is a failure, not a traceback"
    assert "internal error (ZeroDivisionError)" in json.loads(capsys.readouterr().out)[0]["install"]


def test_run_external_logs_and_fails_loudly(tmp_path, monkeypatch):
    ok = rbridge.run_external([sys.executable, "-c", "print('hello')"], workdir=tmp_path / "a", timeout_s=60,
                              version_cmd=[sys.executable, "--version"])
    assert ok.returncode == 0 and "hello" in ok.log.read_text() and ok.version.startswith("Python 3")
    again = rbridge.run_external([sys.executable, "-c", "print('again')"], workdir=tmp_path / "a", timeout_s=60)
    assert again.log != ok.log and "hello" in ok.log.read_text(), "every call keeps its own log"
    assert again.log.name.startswith(Path(sys.executable).stem + "-")
    with pytest.raises(RBridgeError, match="code 4") as info:
        rbridge.run_external([sys.executable, "-c", "import sys; print('bad input'); sys.exit(4)"],
                             workdir=tmp_path / "b", timeout_s=60)
    assert "bad input" in str(info.value)
    started = time.monotonic()
    with pytest.raises(RBridgeError, match="timed out after 1 s; log: "):
        rbridge.run_external([sys.executable, "-c", "import time; time.sleep(60)"], workdir=tmp_path / "c",
                             timeout_s=1)
    assert time.monotonic() - started < 15
    with pytest.raises(RBridgeError, match="cannot start"):
        rbridge.run_external([str(tmp_path / "no_such_tool")], workdir=tmp_path / "d", timeout_s=5)
    tool = _script(tmp_path / "bin" / "tool", "echo from-tool", "@echo from-tool")
    monkeypatch.chdir(tmp_path)
    rel = os.path.join("bin", tool.name)
    res = rbridge.run_external([rel], workdir=tmp_path / "e", timeout_s=60)       # runs in e, found from here
    assert "from-tool" in res.log.read_text()


# --------------------------------------------------------------------------- #
# environment report and packaged R files
# --------------------------------------------------------------------------- #

def test_environment_report_rows(bare, fake_rscript, monkeypatch):
    monkeypatch.setenv("GRNBENCH_RSCRIPT", str(fake_rscript.path))
    rows = rbridge.environment_report(fresh=True)
    assert all(set(r) == {"tool", "found", "version", "source", "used_for", "install"} for r in rows)
    by_tool = {r["tool"]: r for r in rows}
    assert rows[0]["tool"] == "R" and rows[0]["found"] and rows[0]["version"] == "4.5.3"
    assert "GRNBENCH_RSCRIPT" in rows[0]["source"]
    assert {f"r:{p}" for p in rbridge.R_REQUIREMENTS} <= set(by_tool)
    assert by_tool["r:SoupX"]["found"] and by_tool["r:SoupX"]["version"] == "1.6.2"
    assert "ambient_rna=soupx" in by_tool["r:SoupX"]["used_for"]
    assert not by_tool["r:scDblFinder"]["found"] and "BiocManager::install" in by_tool["r:scDblFinder"]["install"]
    assert "py:scanpy" in by_tool and by_tool["py:scanpy"]["found"] and by_tool["py:scanpy"]["version"]
    assert "cli:cellbender" in by_tool
    tools = [r["tool"] for r in rows]
    assert tools.index("py:scanpy") < tools.index("cli:cellbender")
    assert len(fake_rscript.calls()) == 1, "one probe answers every R row"
    rbridge.environment_report(fresh=True)
    assert len(fake_rscript.calls()) == 2, "fresh=True asks R again"


def test_environment_report_never_raises(bare, tmp_path, monkeypatch):
    rows = rbridge.environment_report(str(tmp_path / "missing" / "Rscript"))
    assert rows[0]["tool"] == "R" and not rows[0]["found"] and "no Rscript at" in rows[0]["install"]
    monkeypatch.setenv("GRNBENCH_RSCRIPT", "none")
    rows = rbridge.environment_report()
    assert not rows[0]["found"] and "GRNBENCH_RSCRIPT=none" in rows[0]["install"]
    assert not next(r for r in rows if r["tool"] == "r:SoupX")["found"]

    def broken():
        raise ZeroDivisionError("kb is broken")

    monkeypatch.setattr(rbridge, "_kb_uses", broken)
    rows = rbridge.environment_report()
    assert rows[0]["install"] == "internal error (ZeroDivisionError): kb is broken"


def test_the_r_files_ship_inside_the_package():
    for name in ("_wrapper.R", "_probe.R"):
        assert (importlib.resources.files("grnbench") / "r" / name).is_file(), name


# --------------------------------------------------------------------------- #
# the knowledge base sees R packages through the bridge, only when it is on (step 2c)
# --------------------------------------------------------------------------- #

def test_the_knowledge_base_asks_r_only_when_the_bridge_is_on(bare, fake_rscript, monkeypatch):
    monkeypatch.setenv("GRNBENCH_RSCRIPT", str(fake_rscript.path))     # an R with SoupX, without scDblFinder
    soupx, scdbl = kb.DECISION_INDEX["ambient_rna"], kb.DECISION_INDEX["doublets"]
    monkeypatch.delenv("GRNBENCH_R_BRIDGE", raising=False)
    assert not kb.r_bridge_on()
    assert soupx.availability("soupx") == (False, "missing: r:SoupX"), "off: the registered text"
    kb.option_table()
    assert fake_rscript.calls() == [], "off: planning and the option table never start R"
    for value in ("1", "on", "TRUE", " yes "):
        monkeypatch.setenv("GRNBENCH_R_BRIDGE", value)
        assert kb.r_bridge_on(), value
    assert soupx.option("soupx").availability() == (True, ""), "on: the package is found in R"
    ok, reason = soupx.availability("soupx")
    assert not ok and reason == "no grnbench executor for ambient_rna=soupx yet", "found, but nothing runs it yet"
    assert scdbl.availability("scdblfinder") == (False, "missing: r:scDblFinder")
    assert [c["script"] for c in fake_rscript.calls()] == ["_probe.R"], "one probe answers every R requirement"
    monkeypatch.setenv("GRNBENCH_R_BRIDGE", "0")
    assert soupx.availability("soupx") == (False, "missing: r:SoupX")


def test_the_bridge_on_with_no_usable_r_reads_as_missing_never_as_an_error(bare, fake_rscript, tmp_path,
                                                                           monkeypatch):
    monkeypatch.setenv("GRNBENCH_R_BRIDGE", "1")
    monkeypatch.setenv("GRNBENCH_RSCRIPT", "none")
    assert kb.DECISION_INDEX["ambient_rna"].availability("soupx") == (False, "missing: r:SoupX")
    monkeypatch.setenv("GRNBENCH_RSCRIPT", str(tmp_path / "no" / "Rscript"))     # a named R that is not there
    assert kb.DECISION_INDEX["ambient_rna"].availability("soupx") == (False, "missing: r:SoupX")
    assert fake_rscript.calls() == []


def test_the_bridge_on_starts_a_broken_r_once_per_process_not_once_per_package(bare, tmp_path, monkeypatch):
    starts = tmp_path / "starts.txt"
    broken = _script(tmp_path / "broken", f'echo x >> "{starts}"; exit 1', f'@echo x>> "{starts}"\r\n@exit /b 1')
    monkeypatch.setenv("GRNBENCH_R_BRIDGE", "1")
    monkeypatch.setenv("GRNBENCH_RSCRIPT", str(broken))
    kb.option_table()
    kb.option_table()
    assert kb.DECISION_INDEX["ambient_rna"].availability("soupx") == (False, "missing: r:SoupX")
    assert len(starts.read_text().splitlines()) == 1, "a failing R is asked once, then remembered as unusable"
    rbridge.clear_memo()                                   # a new session (or Check R) asks again
    kb.option_table()
    assert len(starts.read_text().splitlines()) == 2


def test_no_text_names_rpy2_any_more():
    from grnbench import strategies
    texts = [kb._REQ_HINT["r"]] + [s.unavailable or "" for s in strategies.REGISTRY.values()]
    root = Path(__file__).resolve().parents[1]
    texts += [line for line in (root / "README.md").read_text(encoding="utf-8").splitlines() if "scran" in line]
    assert not [t for t in texts if "rpy2" in t], "D-032: R runs through the bridge, never rpy2"


# --------------------------------------------------------------------------- #
# real R (skips without one; on a machine with the grnbench-r env it must run)
# --------------------------------------------------------------------------- #

@pytest.fixture()
def real_r():
    try:
        rt = rbridge.find_rscript()
    except RBridgeError as exc:
        pytest.skip(f"R is configured but does not run: {exc}")
    if rt is None:
        pytest.skip("no R found (GRNBENCH_RSCRIPT, the grnbench-r conda env, R_HOME, PATH)")
    return rt


def test_real_r_is_reproducible_with_a_fixed_seed(real_r, tmp_path):
    rt = real_r
    script = tmp_path / "runif.R"      # with the byte-order mark Windows editors add: R must still read it
    script.write_text('\ufeffx <- runif(5)\nwrite.table(data.frame(x = x), file.path(workdir, expect$table), '
                      'sep = "\\t", row.names = FALSE, quote = FALSE)\n', encoding="utf-8")

    def run(seed: int, name: str):
        return rbridge.run_r(script, seed=seed, workdir=tmp_path / name, rt=rt, expect={"table": "out.tsv"},
                             timeout_s=300)

    a, b, c = run(7, "a"), run(7, "b"), run(8, "c")
    first = a.outputs["table"].read_bytes()
    assert first == b.outputs["table"].read_bytes() and first != c.outputs["table"].read_bytes()
    assert len(first.splitlines()) == 6
    assert a.session["r_version"].startswith("R version") and a.session["seed"] == 7
    assert a.session["rng_kind"] == ["Mersenne-Twister", "Inversion", "Rejection"]
    assert "jsonlite" in a.session["packages"]
    real = lambda p: os.path.normcase(os.path.realpath(p))  # noqa: E731 - symlinked homes on the cluster
    assert real(a.session["lib_paths"][0]) == real(rt.library)
    if rbridge.probe(rt, ["SoupX"])["packages"]["SoupX"]:
        assert rbridge.requirement_status("r:SoupX", rt).ok
    print(f"real R: {rt.version} from {rt.source}, {a.seconds:.2f} s per run")


def test_real_r_records_a_method_that_changes_the_generator(real_r, tmp_path):
    script = tmp_path / "rng.R"
    script.write_text('RNGkind("Wichmann-Hill")\nset.seed(99)\n'
                      'writeLines(format(runif(2)), file.path(workdir, expect$out))\n', encoding="utf-8")
    block = rbridge.run_r(script, seed=1, workdir=tmp_path / "w", rt=real_r, expect={"out": "o.txt"},
                          timeout_s=300).manifest_block()
    assert block["rng_changed"] is True and block["rng_kind"][0] == "Mersenne-Twister"
    assert block["rng_kind_after"][0] == "Wichmann-Hill"
    assert block["set_seed_calls"] == [{"seed": 99, "kind": None, "normal.kind": None, "sample.kind": None}]
    plain = tmp_path / "plain.R"
    plain.write_text("writeLines(format(runif(2)), file.path(workdir, expect$out))\n", encoding="utf-8")
    block = rbridge.run_r(plain, seed=1, workdir=tmp_path / "p", rt=real_r, expect={"out": "o.txt"},
                          timeout_s=300).manifest_block()
    assert block["rng_changed"] is False and block["set_seed_calls"] == []


def test_real_r_works_in_a_folder_outside_the_windows_code_page(real_r, tmp_path):
    work = tmp_path / "run_ψ_ł"          # neither character is in cp1252
    script = tmp_path / "list.R"
    script.write_text('stopifnot(is.list(params$v), params$n == 3)\n'
                      'writeLines(c(workdir, inputs$data), file.path(workdir, expect$out), useBytes = FALSE)\n',
                      encoding="utf-8")
    data = tmp_path / "dati_ψ.txt"
    data.write_text("x", encoding="utf-8")
    res = rbridge.run_r(script, seed=1, workdir=work, rt=real_r, inputs={"data": data},
                        params={"v": [1, 2], "n": 3}, expect={"out": "o.txt"}, timeout_s=300)
    lines = res.outputs["out"].read_text(encoding="utf-8").splitlines()
    assert Path(lines[0]) == work.resolve() and Path(lines[1]) == data.resolve()


def test_real_r_error_reaches_the_message(real_r, tmp_path):
    script = tmp_path / "boom.R"
    script.write_text('stop("boom from the method")\n', encoding="utf-8")
    with pytest.raises(RBridgeError, match="boom from the method"):
        rbridge.run_r(script, seed=1, workdir=tmp_path / "w", rt=real_r, timeout_s=300)
