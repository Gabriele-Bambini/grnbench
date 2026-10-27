"""ARACNe3 adapter: the exact command, the input files, the output parser, availability.

The real program never runs here (it is a native binary, not installed on laptops or CI).
``subprocess.run`` is replaced by a stand-in that records the call, reads the files the adapter
wrote and answers the way ARACNe3 does at califano-lab/ARACNe3 commit 3d8791a
(``src/app/ARACNe3.cpp`` for the flags, ``src/app/io.cpp`` for the file formats).
"""
from __future__ import annotations

import contextlib
import errno
import json
import os
import stat
import subprocess
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from grnbench.infer import MethodError, list_methods, run_method
from grnbench.infer import aracne3

# the header writeConsolidatedNetwork() prints (src/app/io.cpp)
CONSOLIDATED_HEADER = "regulator.values\ttarget.values\tmi.values\tscc.values\tcount.values\tlog.p.values"
NETWORK_ROWS = [("TF_A", "G_B", 0.42, 0.9, 1, -7.1), ("TF_A", "G_C", 0.30, -0.6, 1, -7.1)]


class FakeAracne3:
    """Stands in for ``subprocess.run``: records argv and keyword arguments, keeps the text of the
    ``-e`` and ``-r`` files, and writes ``<-o>/consolidated-net_<--runid>.tsv`` like ARACNe3."""

    def __init__(self):
        self.calls: list[tuple[list[str], dict]] = []
        self.inputs: dict[str, str] = {}
        self.returncode, self.stdout, self.stderr = 0, "", ""
        self.write_network, self.time_out = True, False

    def __call__(self, argv, **kwargs):
        argv = list(argv)
        self.calls.append((argv, kwargs))
        if self.time_out:
            raise subprocess.TimeoutExpired(argv, kwargs.get("timeout"))
        value = {flag: argv[i + 1] for i, flag in enumerate(argv[:-1]) if flag.startswith("-")}
        for flag in ("-e", "-r"):
            self.inputs[flag] = Path(value[flag]).read_bytes().decode("utf-8")
        if self.returncode == 0 and self.write_network:
            out = Path(value["-o"])
            out.mkdir(exist_ok=True)          # the real program creates its output folder
            lines = [CONSOLIDATED_HEADER] + ["\t".join(map(str, row)) for row in NETWORK_ROWS]
            (out / f"consolidated-net_{value['--runid']}.tsv").write_text("\n".join(lines) + "\n")
        return subprocess.CompletedProcess(argv, self.returncode, self.stdout, self.stderr)


@pytest.fixture
def fake_bin(tmp_path, monkeypatch) -> Path:
    """An executable file standing in for the built ARACNe3_app_release, named by GRNBENCH_ARACNE3_BIN."""
    exe = tmp_path / "bin" / "ARACNe3_app_release"
    exe.parent.mkdir()
    exe.write_text("#!/bin/sh\nexit 99\n")        # never executed: subprocess.run is replaced
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv(aracne3.ENV_BIN, str(exe))
    return exe


@pytest.fixture
def fake_run(monkeypatch) -> FakeAracne3:
    fake = FakeAracne3()
    monkeypatch.setattr(aracne3.subprocess, "run", fake)
    return fake


@pytest.fixture
def fake_run_in_adapter(monkeypatch) -> FakeAracne3:
    """``fake_run`` for the adapter only: a CLI around it still calls the real ``subprocess.run``
    (``platform.platform()`` does on Windows, for the manifest's environment block)."""
    fake = FakeAracne3()
    monkeypatch.setattr(aracne3, "subprocess", SimpleNamespace(run=fake, TimeoutExpired=subprocess.TimeoutExpired))
    return fake


@pytest.fixture
def workdir(tmp_path, monkeypatch) -> Path:
    """The adapter's temporary folder, made a known path so the argv can be compared literally."""
    work = tmp_path / "work"
    work.mkdir()
    monkeypatch.setattr(aracne3, "tempfile",
                        SimpleNamespace(TemporaryDirectory=lambda **_: contextlib.nullcontext(str(work))))
    return work


@pytest.fixture
def expr() -> pd.DataFrame:
    """Six cells x three genes, cells as rows (the orientation the adapter contract hands over)."""
    values = np.random.default_rng(3).gamma(2.0, 1.5, size=(6, 3))
    values[0, 0] = 1.234567891e-7                   # rounds to 0 at six decimals
    return pd.DataFrame(values, index=[f"cell{i}" for i in range(6)], columns=["TF_A", "G_B", "G_C"])


# --------------------------------------------------------------------------- #
# the command
# --------------------------------------------------------------------------- #

def test_runner_calls_the_native_program_with_the_flags_aracne3_reads(fake_bin, workdir, fake_run, expr):
    aracne3.run_aracne3(expr, tf_names=["TF_A"])
    (argv, _), = fake_run.calls
    assert argv == [
        str(fake_bin),
        "-e", str(workdir / "expr.tsv"),
        "-r", str(workdir / "regulators.txt"),
        "-o", str(workdir / "out"),
        "--seed", "42",
        "--threads", "1",
        "--runid", "grnbench",
    ]


def test_runner_starts_aracne3_inside_its_own_temporary_folder(fake_bin, workdir, fake_run, expr):
    """ARACNe3 caches its null model in ./.ARACNe3_cached/ and skips the seeded draws when one exists,
    so a working folder shared between runs would make a run depend on the runs before it."""
    aracne3.run_aracne3(expr)
    (_, kwargs), = fake_run.calls
    assert kwargs["cwd"] == str(workdir)


def test_seed_aracne3_threads_and_extra_args_reach_the_command(fake_bin, workdir, fake_run, expr):
    edges = aracne3.run_aracne3(expr, seed=7, aracne3_threads=8,
                                extra_args=["-x", "30", "--alpha", "0.01", "--noMaxEnt"])
    (argv, _), = fake_run.calls
    assert argv[7:] == ["--seed", "7", "--threads", "8", "--runid", "grnbench",
                        "-x", "30", "--alpha", "0.01", "--noMaxEnt"]
    assert edges.attrs["threads"] == 8 and "threads_not_used" not in edges.attrs


# --------------------------------------------------------------------------- #
# threads: the engine's generic count never reaches ARACNe3
# --------------------------------------------------------------------------- #

def test_the_generic_threads_the_clis_pass_does_not_reach_aracne3(fake_bin, workdir, fake_run, cascade_adata):
    """grnbench infer and grnbench run pass threads=4 (hpc/grnbench.sbatch the job's CPU count) to every
    method. With more than one thread ARACNe3's null model does not repeat from the seed
    (EXP-20261018-GB-1), so aracne3 stays single-threaded, and its provenance says what ARACNe3 got."""
    edges, prov = run_method("aracne3", cascade_adata, tf_list=["TF_A"], seed=42, threads=4)
    (argv, _), = fake_run.calls
    assert argv[argv.index("--threads") + 1] == "1"
    assert prov["params"]["threads"] == "4"                  # what the caller asked of every method
    assert prov["method_info"]["threads"] == 1               # what ARACNe3 got
    assert prov["method_info"]["threads_not_used"].startswith("threads=4 ")
    assert "aracne3_threads" in prov["method_info"]["threads_not_used"]


def test_an_explicit_aracne3_threads_request_is_honoured_and_recorded(fake_bin, workdir, fake_run, cascade_adata):
    edges, prov = run_method("aracne3", cascade_adata, tf_list=["TF_A"], seed=42, threads=4, aracne3_threads=3)
    (argv, _), = fake_run.calls
    assert argv[argv.index("--threads") + 1] == "3"
    assert prov["params"]["aracne3_threads"] == "3"
    assert prov["method_info"]["threads"] == 3


def test_grnbench_infer_runs_aracne3_single_threaded_and_records_it(fake_bin, workdir, fake_run_in_adapter,
                                                                   cascade_adata, tmp_path):
    from grnbench.cli_infer import main as infer_main

    data, out = tmp_path / "data.h5ad", tmp_path / "infer"
    cascade_adata.write_h5ad(data)
    assert infer_main([str(data), "--out", str(out), "--methods", "aracne3"]) == 0   # --threads defaults to 4
    (argv, _), = fake_run_in_adapter.calls
    assert argv[argv.index("--threads") + 1] == "1"
    manifest = json.loads((out / "inference_manifest.json").read_text(encoding="utf-8"))
    assert manifest["args"]["threads"] == 4
    assert manifest["runs"]["aracne3"]["method_info"]["threads"] == 1


@pytest.mark.parametrize("bad", [0, -2])
def test_an_aracne3_thread_count_below_one_is_refused_before_running(bad, fake_bin, fake_run, expr):
    with pytest.raises(MethodError, match="aracne3_threads"):
        aracne3.run_aracne3(expr, aracne3_threads=bad)
    assert fake_run.calls == []


@pytest.mark.parametrize("kwargs", [{"n_threads": 8}, {"timeout": 10}, {"extra": ["-x", "3"]}])
def test_a_parameter_aracne3_does_not_take_is_refused_not_dropped(kwargs, fake_bin, fake_run, expr):
    name = next(iter(kwargs))
    with pytest.raises(MethodError, match=f"aracne3 does not take {name}.*aracne3_threads"):
        aracne3.run_aracne3(expr, **kwargs)
    assert fake_run.calls == []


@pytest.mark.parametrize("args, message", [
    (["-o", "elsewhere"], "cannot be passed"),
    (["--seed", "7"], "cannot be passed"),
    (["--threads", "8"], "cannot be passed"),
    (["--runid", "mine"], "cannot be passed"),
    (["--noConsolidate"], "cannot be passed"),
    (["--noAlpha"], "spells it --noalpha"),          # the README's spelling, which the code never reads
    (["-FWER"], "spells it --FWER"),                 # README example 5
    (["--bogus"], "does not read"),
    (["--alpha"], "needs a value"),
])
def test_extra_args_that_aracne3_would_misread_or_silently_ignore_are_refused(args, message, fake_bin,
                                                                              fake_run, expr):
    with pytest.raises(MethodError, match=message):
        aracne3.run_aracne3(expr, extra_args=args)
    assert fake_run.calls == []


# --------------------------------------------------------------------------- #
# the input files
# --------------------------------------------------------------------------- #

def test_expression_file_has_genes_as_rows_and_cells_as_columns(fake_bin, workdir, fake_run, expr):
    aracne3.run_aracne3(expr)
    text = fake_run.inputs["-e"]
    assert "\r" not in text and text.endswith("\n")
    header, *rows = text[:-1].split("\n")
    assert header.count("\t") == len(expr)           # ARACNe3 counts samples from the header's tabs
    assert header.split("\t")[1:] == list(expr.index)
    assert [row.split("\t")[0] for row in rows] == list(expr.columns)
    values = np.array([[float(v) for v in row.split("\t")[1:]] for row in rows])
    np.testing.assert_allclose(values, expr.to_numpy().T, rtol=1e-8)   # no ties made by rounding


def test_regulator_file_lists_the_tfs_that_are_genes_of_the_matrix_one_per_line(fake_bin, workdir,
                                                                                fake_run, expr):
    aracne3.run_aracne3(expr, tf_names=["TF_A", "NOT_A_GENE", "G_C"])
    assert fake_run.inputs["-r"] == "TF_A\nG_C\n"


def test_without_a_tf_list_every_gene_is_a_candidate_regulator(fake_bin, workdir, fake_run, expr):
    aracne3.run_aracne3(expr, tf_names=None)
    assert fake_run.inputs["-r"] == "TF_A\nG_B\nG_C\n"


def test_a_tf_list_with_no_gene_of_the_matrix_is_refused_before_running(fake_bin, fake_run, expr):
    with pytest.raises(MethodError, match="none of the 2 regulators"):
        aracne3.run_aracne3(expr, tf_names=["X1", "X2"])
    assert fake_run.calls == []


@pytest.mark.parametrize("shape", [(65_536, 2), (2, 65_536)])
def test_matrices_beyond_aracne3s_16_bit_counters_are_refused_before_running(shape, fake_bin, fake_run):
    """More cells than 65,535 make ARACNe3 exit with a row-length error; more genes wrap its 16-bit gene
    ids without a word. Both are refused with the reason before anything is written."""
    big = pd.DataFrame(np.zeros(shape), columns=[f"g{i}" for i in range(shape[1])])
    with pytest.raises(MethodError, match="65,535"):
        aracne3.run_aracne3(big)
    assert fake_run.calls == []


# --------------------------------------------------------------------------- #
# the output
# --------------------------------------------------------------------------- #

def test_parser_reads_the_consolidated_network_aracne3_writes():
    text = (CONSOLIDATED_HEADER + "\n"
            "TP53\tCDKN1A\t0.52\t0.61\t3\t-12.5\n"
            "MYC\tMDM2\t0.31\t-0.44\t1\t-3.2\n"
            "MYC\tNPM1\t0.12\t0\t1\t-3.2\n")
    net = aracne3.parse_aracne3_text(text)
    assert list(net.columns) == ["source", "target", "weight", "direction"]
    assert net.to_dict("records") == [
        {"source": "TP53", "target": "CDKN1A", "weight": 0.52, "direction": 1},   # weight = MI
        {"source": "MYC", "target": "MDM2", "weight": 0.31, "direction": -1},     # sign of Spearman
        {"source": "MYC", "target": "NPM1", "weight": 0.12, "direction": 0},
    ]


def test_parser_refuses_a_table_that_is_not_aracne3s_consolidated_network():
    with pytest.raises(MethodError, match="regulator.values"):
        aracne3.parse_aracne3_text("Regulator\tTarget\tScore\tSign\nTP53\tCDKN1A\t12.4\t+\n")


def test_a_consolidated_network_without_edges_is_an_error_that_says_how_to_relax_the_pruning():
    with pytest.raises(MethodError, match=r"kept no edges.*--alpha"):
        aracne3.parse_aracne3_text(CONSOLIDATED_HEADER + "\n")


def test_run_method_turns_the_consolidated_network_into_canonical_signed_edges(fake_bin, fake_run,
                                                                              cascade_adata):
    edges, prov = run_method("aracne3", cascade_adata, tf_list=["TF_A"])
    assert list(edges.columns) == ["source", "target", "weight", "direction"]
    assert edges[["source", "target"]].to_numpy().tolist() == [["TF_A", "G_B"], ["TF_A", "G_C"]]
    assert edges["weight"].tolist() == [0.42, 0.30]
    assert edges["direction"].tolist() == [1, -1]
    assert prov["regulators"] == "tf_list"


def test_a_failed_run_is_a_method_error_with_the_programs_last_output(fake_bin, fake_run, expr):
    fake_run.returncode = 1
    fake_run.stdout = "usage: ARACNe3_app_release -e path/to/matrix.txt -r path/to/regulators.txt -o path/to/output"
    with pytest.raises(MethodError, match=r"rc=1.*usage: ARACNe3_app_release"):
        aracne3.run_aracne3(expr)


def test_a_run_that_writes_no_consolidated_network_is_an_error(fake_bin, fake_run, expr):
    fake_run.write_network = False
    with pytest.raises(MethodError, match="consolidated-net_grnbench.tsv"):
        aracne3.run_aracne3(expr)


@pytest.mark.parametrize("returncode, write_network", [(1, True), (0, False)])
def test_the_programs_output_in_an_error_names_the_temporary_folder_as_tmp(returncode, write_network, fake_bin,
                                                                          workdir, fake_run, expr):
    """ARACNe3's messages name the files it was given, all inside the adapter's temporary folder, which
    can be under home; grnbench infer writes the error text into inference_manifest.json."""
    fake_run.returncode, fake_run.write_network = returncode, write_network
    fake_run.stderr = f'error: file open failed "{workdir / "expr.tsv"}".'
    with pytest.raises(MethodError) as err:
        aracne3.run_aracne3(expr)
    assert str(workdir) not in str(err.value)
    assert f'"<tmp>{os.sep}expr.tsv"' in str(err.value)


def test_the_programs_path_in_its_usage_line_is_written_as_its_name(fake_bin, workdir, fake_run, expr):
    """ARACNe3 prints ``usage: <argv[0]> -e ...`` (src/app/ARACNe3.cpp): argv[0] is the program's full path."""
    fake_run.returncode = 1
    fake_run.stdout = f"usage: {fake_bin} -e path/to/matrix.txt -r path/to/regulators.txt -o path/to/output/directory"
    with pytest.raises(MethodError, match="usage: ARACNe3_app_release -e path/to/matrix.txt") as err:
        aracne3.run_aracne3(expr)
    assert str(fake_bin.parent) not in str(err.value)


def test_a_failed_run_reaches_the_inference_manifest_without_the_temporary_folder(fake_bin, workdir,
                                                                                  fake_run_in_adapter,
                                                                                  cascade_adata, tmp_path):
    from grnbench.cli_infer import main as infer_main

    data, out = tmp_path / "data.h5ad", tmp_path / "infer"
    cascade_adata.write_h5ad(data)
    fake = fake_run_in_adapter
    fake.returncode = 1
    fake.stderr = f'Directory Created: "{workdir / "out"}".\nerror: file open failed "{workdir / "expr.tsv"}".'
    assert infer_main([str(data), "--out", str(out), "--methods", "aracne3"]) == 1
    text = (out / "inference_manifest.json").read_text(encoding="utf-8")
    error = json.loads(text)["errors"]["aracne3"]
    assert error.startswith("MethodError: aracne3 failed (rc=1)") and "<tmp>" in error
    assert str(workdir) not in text and json.dumps(str(workdir))[1:-1] not in text


def test_a_run_past_its_timeout_is_a_method_error(fake_bin, fake_run, expr):
    fake_run.time_out = True
    with pytest.raises(MethodError, match="timed out after 5 s"):
        aracne3.run_aracne3(expr, timeout_s=5)


class BusyFolder:
    """A temporary folder the file system refuses to delete, as the cluster's BeeGFS home did after a
    finished run in EXP-20261019-GB-1 (``OSError: [Errno 16] Device or resource busy: 'subnets_log'``)."""

    def __init__(self, path: Path):
        self.path = path

    def __enter__(self) -> str:
        return str(self.path)

    def __exit__(self, *exc_info):
        raise OSError(errno.EBUSY, os.strerror(errno.EBUSY), "subnets_log")


@pytest.fixture
def busy_folder(tmp_path, monkeypatch) -> Path:
    work = tmp_path / "grnbench_aracne3_busy"
    work.mkdir()
    monkeypatch.setattr(aracne3, "tempfile", SimpleNamespace(TemporaryDirectory=lambda **_: BusyFolder(work)))
    return work


def test_a_finished_network_is_kept_when_its_temporary_folder_cannot_be_deleted(fake_bin, busy_folder,
                                                                                fake_run, expr):
    with pytest.warns(RuntimeWarning, match="ARACNE3 FOLDER LEFT.*grnbench_aracne3_busy"):
        edges = aracne3.run_aracne3(expr, tf_names=["TF_A"])
    assert edges[["source", "target"]].to_numpy().tolist() == [["TF_A", "G_B"], ["TF_A", "G_C"]]
    assert edges.attrs["temporary_folder_left"] == (
        f"grnbench_aracne3_busy (OSError: [Errno {errno.EBUSY}] {os.strerror(errno.EBUSY)})")


def test_the_folder_left_behind_is_named_in_the_run_provenance(fake_bin, busy_folder, fake_run, cascade_adata):
    with pytest.warns(RuntimeWarning, match="ARACNE3 FOLDER LEFT"):
        edges, prov = run_method("aracne3", cascade_adata, tf_list=["TF_A"])
    assert len(edges) == 2
    assert prov["method_info"]["temporary_folder_left"].startswith("grnbench_aracne3_busy (OSError")


def test_a_failed_run_reports_its_own_error_when_its_folder_cannot_be_deleted_either(fake_bin, busy_folder,
                                                                                     fake_run, expr):
    fake_run.returncode, fake_run.stderr = 1, "gene GENE_X has no values"
    with pytest.raises(MethodError, match=r"rc=1.*GENE_X"):
        aracne3.run_aracne3(expr)


def test_an_os_error_of_the_run_is_raised_as_it_is_not_as_the_error_it_was_raised_from(fake_bin, workdir,
                                                                                         monkeypatch, expr):
    def refused(argv, **kwargs):
        try:
            raise ValueError("an inner detail of the library")
        except ValueError:
            raise PermissionError(errno.EACCES, "permission denied", argv[0])
    monkeypatch.setattr(aracne3.subprocess, "run", refused)
    with pytest.raises(PermissionError, match="permission denied"):
        aracne3.run_aracne3(expr)


# --------------------------------------------------------------------------- #
# availability
# --------------------------------------------------------------------------- #

def test_without_the_variable_or_the_program_on_path_aracne3_is_unavailable_and_says_what_to_set(monkeypatch):
    monkeypatch.delenv(aracne3.ENV_BIN, raising=False)
    monkeypatch.delenv("ARACNE3_JAR", raising=False)
    monkeypatch.setattr(aracne3.shutil, "which", lambda name: None)
    ok, why = aracne3.aracne3_availability()
    assert not ok
    assert aracne3.ENV_BIN in why and "ARACNe3_app_release" in why
    assert "java" not in why.lower()


def test_an_old_aracne3_jar_variable_is_named_as_no_longer_used(monkeypatch):
    monkeypatch.delenv(aracne3.ENV_BIN, raising=False)
    monkeypatch.setenv("ARACNE3_JAR", "aracne3.jar")
    monkeypatch.setattr(aracne3.shutil, "which", lambda name: None)
    ok, why = aracne3.aracne3_availability()
    assert not ok and "ARACNE3_JAR is no longer used" in why


def test_the_program_is_found_on_path_under_its_cmake_names(monkeypatch):
    monkeypatch.delenv(aracne3.ENV_BIN, raising=False)
    asked = []

    def which(name):
        asked.append(name)
        return f"/opt/aracne3/{name}" if name == "ARACNe3_app" else None

    monkeypatch.setattr(aracne3.shutil, "which", which)
    assert aracne3.aracne3_availability() == (True, "")
    assert asked == ["ARACNe3_app_release", "ARACNe3_app"]


def test_a_wrong_variable_is_reported_not_replaced_by_a_program_on_path(monkeypatch, tmp_path):
    monkeypatch.setenv(aracne3.ENV_BIN, str(tmp_path / "not-there"))
    monkeypatch.setattr(aracne3.shutil, "which", lambda name: f"/usr/local/bin/{name}")
    ok, why = aracne3.aracne3_availability()
    assert not ok and aracne3.ENV_BIN in why
    assert str(tmp_path) not in why          # the reason is written into run manifests: no personal paths


def test_the_variable_pointing_at_an_executable_makes_aracne3_available(fake_bin):
    assert aracne3.aracne3_availability() == (True, "")


def test_a_relative_variable_is_run_by_its_absolute_path_from_the_temporary_folder(fake_bin, workdir, fake_run,
                                                                                   expr, monkeypatch):
    """The program starts with cwd = its temporary folder, where a relative path would not resolve:
    available must mean runnable, so the path the run gets is absolute."""
    monkeypatch.chdir(fake_bin.parent.parent)
    monkeypatch.setenv(aracne3.ENV_BIN, os.path.join(".", "bin", fake_bin.name))
    assert aracne3.aracne3_availability() == (True, "")
    aracne3.run_aracne3(expr)
    (argv, kwargs), = fake_run.calls
    assert os.path.isabs(argv[0]) and Path(argv[0]) == fake_bin
    assert kwargs["cwd"] == str(workdir)


@pytest.mark.skipif(os.name == "nt", reason="Windows files have no execute bit")
def test_a_binary_without_execute_permission_is_unavailable(monkeypatch, tmp_path):
    exe = tmp_path / "ARACNe3_app_release"
    exe.write_text("not executable")
    exe.chmod(0o644)
    monkeypatch.setenv(aracne3.ENV_BIN, str(exe))
    ok, why = aracne3.aracne3_availability()
    assert not ok and "not executable" in why


def test_registry_lists_aracne3_as_unavailable_and_refuses_to_run_it_without_the_program(monkeypatch, tmp_path,
                                                                                         cascade_adata, expr):
    monkeypatch.setenv(aracne3.ENV_BIN, str(tmp_path / "missing"))
    entry = next(m for m in list_methods() if m["name"] == "aracne3")
    assert not entry["available"] and aracne3.ENV_BIN in entry["unavailable_reason"]
    with pytest.raises(MethodError, match="unavailable"):
        run_method("aracne3", cascade_adata)
    with pytest.raises(MethodError, match="aracne3 unavailable"):
        aracne3.run_aracne3(expr)
