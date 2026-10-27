"""benchmarks/aracne3_smoke.py (EXP-20261018-GB-1): the wrapper, the comparisons, the subset, the scrub.

The real ARACNe3 never runs here: ``subprocess.run`` is replaced by a stand-in that reads the
adapter's input files and answers the way ARACNe3 does at califano-lab/ARACNe3 commit 3d8791a
(stdout sample counts, ``log_<runid>.txt`` timings, ``./.ARACNe3_cached/``, the consolidated file).
The smoke script's wrapper sits between the adapter and that stand-in, as it sits between the
adapter and the real program on the HPC.
"""
from __future__ import annotations

import json
import math
import stat
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
import aracne3_smoke as smoke  # noqa: E402

from grnbench.infer import aracne3  # noqa: E402

HEADER = "regulator.values\ttarget.values\tmi.values\tscc.values\tcount.values\tlog.p.values"


class FakeProgram:
    """Answers like ARACNe3: counts samples from the header's tabs, writes its log, its null-model
    cache in the working folder and ``<-o>/consolidated-net_<--runid>.tsv``. ``shift_from`` changes
    the weights (and the null model) from that call on; ``empty_unless`` keeps no edges unless that
    flag is in argv."""

    def __init__(self, shift_from: int | None = None, empty_unless: str | None = None):
        self.n, self.shift_from, self.empty_unless = 0, shift_from, empty_unless

    def __call__(self, argv, **kwargs):
        self.n += 1
        argv = [str(a) for a in argv]
        value = {flag: argv[i + 1] for i, flag in enumerate(argv[:-1]) if flag.startswith("-")}
        lines = Path(value["-e"]).read_text(encoding="utf-8").splitlines()
        samples = lines[0].count("\t")
        genes = [line.split("\t")[0] for line in lines[1:]]
        regs = Path(value["-r"]).read_text(encoding="utf-8").split()
        shift = 0.01 if self.shift_from is not None and self.n >= self.shift_from else 0.0
        out = Path(value["-o"])
        out.mkdir(exist_ok=True)
        cache = Path(kwargs["cwd"]) / ".ARACNe3_cached"
        cache.mkdir(exist_ok=True)
        (cache / f"Nssamp-{math.ceil(0.632 * samples)}_Nnull-1000000.txt").write_text(f"{0.1 + shift}\n")
        (out / f"log_{value['--runid']}.txt").write_text(     # whole seconds, as ARACNe3 writes them
            " ".join(argv) + " \n\nGene expression matrix & regulators list read time: 0s\n"
            "\nMutual Information null model calculation time: 5s\nSUCCESS!\n")
        rows = []
        if self.empty_unless is None or self.empty_unless in argv:
            rows = [f"{r}\t{genes[(i + 1) % len(genes)]}\t{0.5 - 0.1 * i + shift}\t{(-1) ** i * 0.4}\t1\t-9.2"
                    for i, r in enumerate(regs[:3])]
        (out / f"consolidated-net_{value['--runid']}.tsv").write_text("\n".join([HEADER, *rows]) + "\n")
        stdout = (f"Beginning ARACNe3 instance.  See logs and progress reports in \"{out}/log_grnbench.txt\".\n"
                  f"\nTotal N Samples: {samples}\nSubsampled N Samples: {math.ceil(0.632 * samples)}\nSUCCESS!\n")
        return subprocess.CompletedProcess(argv, 0, stdout, "")


@pytest.fixture
def program(tmp_path, monkeypatch):
    exe = tmp_path / "home" / "tools" / "ARACNe3_app_release"
    exe.parent.mkdir(parents=True)
    exe.write_text("#!/bin/sh\nexit 99\n")        # never executed: subprocess.run is replaced
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv(aracne3.ENV_BIN, str(exe))

    def install(fake: FakeProgram) -> FakeProgram:
        monkeypatch.setattr(subprocess, "run", fake)
        return fake
    return install


@pytest.fixture
def scrub(tmp_path) -> smoke.Scrubber:
    return smoke.Scrubber(home=str(tmp_path / "home"), user="someone")


def _strings(obj):
    if isinstance(obj, str):
        yield obj
    elif isinstance(obj, dict):
        for v in obj.values():
            yield from _strings(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from _strings(v)


def _adata(n_cells: int = 12, n_genes: int = 6) -> ad.AnnData:
    rng = np.random.default_rng(1)
    a = ad.AnnData(rng.random((n_cells, n_genes)))
    a.obs_names = [f"c{i}" for i in range(n_cells)]
    a.var_names = ["TF1", "TF2", "TF3", "G4", "G5", "G6"][:n_genes]
    return a


def test_runs_parse_and_identical_repeats_compare_identical(program, scrub, tmp_path):
    fake = program(FakeProgram())
    out = tmp_path / "out"
    result = smoke.smoke(_adata(), ["TF1", "TF2", "TF3"], [1, 1, 4], 42, out, scrub)
    assert [r["label"] for r in result["runs"]] == ["t1_a", "t1_b", "t4_a"] and fake.n == 3
    first = result["runs"][0]
    assert first["status"] == "ok" and first["parsed"] and first["program_calls"] == 1
    assert first["aracne3"]["returncode"] == 0 and first["n_edges"] == 3
    # the header the adapter wrote has one tab per cell, and ARACNe3's own count is recorded
    assert first["aracne3"]["inputs"]["-e"]["first_line_tabs"] == 12 and first["total_samples"] == 12
    assert first["aracne3"]["inputs"]["-r"]["lines"] == 3
    assert first["log_times_s"]["Mutual Information null model calculation"] == 5.0
    assert first["consolidated"]["header"] == HEADER.split("\t") and first["consolidated"]["rows"] == 3
    assert first["consolidated"]["columns"] == {
        "regulators_with_edges": 3, "mi_min": pytest.approx(0.3), "mi_median": pytest.approx(0.4),
        "mi_max": pytest.approx(0.5), "mi_zero_rows": 0, "scc_positive": 2, "scc_negative": 1, "scc_zero": 0,
        "count_values": [1], "distinct_log_p_values": 1}
    argv = [a.replace("\\", "/") for a in first["aracne3"]["argv"]]
    assert argv[:3] == ["~/tools/ARACNe3_app_release", "-e", "<tmp>/expr.tsv"]
    assert argv[3:] == ["-r", "<tmp>/regulators.txt", "-o", "<tmp>/out", "--seed", "42",
                        "--threads", "1", "--runid", "grnbench"]
    assert (out / "runs" / "t1_a" / "aracne3_out" / "consolidated-net_grnbench.tsv").exists()
    assert [(c["a"], c["b"]) for c in result["comparisons"]] == [("t1_a", "t1_b"), ("t1_a", "t4_a")]
    assert all(c["identical_edges"] and c["same_consolidated_file"] and c["same_null_model"]
               for c in result["comparisons"])
    assert result["verdict"]["every_run_returned_a_network"]


def _threads_of(run: dict) -> str:
    argv = run["aracne3"]["argv"]
    return argv[argv.index("--threads") + 1]


def test_every_run_gets_the_thread_count_of_its_label(program, scrub, tmp_path):
    """The adapter does not pass the engine's generic ``threads`` on (review of #106), so the script asks
    with ``aracne3_threads``: every run, not only the first, gets the count it is labelled and compared with."""
    fake = program(FakeProgram())
    result = smoke.smoke(_adata(), ["TF1", "TF2", "TF3"], [1, 1, 4, 4], 42, tmp_path / "out", scrub)
    got = [(r["label"], r["threads"], _threads_of(r), r["provenance"]["method_info"]["threads"])
           for r in result["runs"]]
    assert got == [("t1_a", 1, "1", 1), ("t1_b", 1, "1", 1), ("t4_a", 4, "4", 4), ("t4_b", 4, "4", 4)]
    assert fake.n == 4


def test_the_relaxed_repeats_keep_the_thread_count_of_their_label(program, scrub, tmp_path):
    program(FakeProgram(empty_unless="--alpha"))
    result = smoke.smoke(_adata(), ["TF1", "TF2", "TF3"], [1, 4], 42, tmp_path / "out", scrub)
    assert [(r["label"], _threads_of(r)) for r in result["runs"]] == [
        ("t1_a", "1"), ("t4_a", "4"), ("t1_a_alpha0.2", "1"), ("t4_a_alpha0.2", "4")]


def test_a_changed_run_is_reported_as_different(program, scrub, tmp_path):
    program(FakeProgram(shift_from=3))
    result = smoke.smoke(_adata(), ["TF1", "TF2", "TF3"], [1, 1, 4, 4], 42, tmp_path / "out", scrub)
    by_pair = {(c["a"], c["b"]): c for c in result["comparisons"]}
    assert by_pair[("t1_a", "t1_b")]["identical_edges"]
    cross = by_pair[("t1_a", "t4_a")]
    assert not cross["identical_edges"] and not cross["same_consolidated_file"] and not cross["same_null_model"]
    assert cross["shared"] == 3 and cross["jaccard"] == 1.0 and cross["max_abs_weight_diff"] == pytest.approx(0.01)
    assert by_pair[("t4_a", "t4_b")]["identical_edges"]


def test_no_edges_repeats_every_run_relaxed(program, scrub, tmp_path):
    program(FakeProgram(empty_unless="--alpha"))
    result = smoke.smoke(_adata(), ["TF1", "TF2", "TF3"], [1, 4], 42, tmp_path / "out", scrub)
    assert result["relaxed"]
    assert [r["label"] for r in result["runs"]] == ["t1_a", "t4_a", "t1_a_alpha0.2", "t4_a_alpha0.2"]
    assert [r["status"] for r in result["runs"]] == ["error", "error", "ok", "ok"]
    assert smoke.NO_EDGES in result["runs"][0]["error"]
    assert result["runs"][2]["extra_args"] == smoke.RELAX
    assert result["runs"][2]["aracne3"]["argv"][-2:] == smoke.RELAX
    assert [(c["a"], c["b"]) for c in result["comparisons"]] == [("t1_a", "t4_a"),
                                                                ("t1_a_alpha0.2", "t4_a_alpha0.2")]
    assert not result["verdict"]["every_run_returned_a_network"]


def test_summary_names_neither_home_nor_user_and_the_scan_finds_a_leak(program, scrub, tmp_path):
    program(FakeProgram())
    out = tmp_path / "out"
    result = smoke.smoke(_adata(), ["TF1", "TF2", "TF3"], [1], 42, out, scrub)
    summary = {"experiment": smoke.EXP, "code_commit": "abc", "aracne3": {"source_commit": "3d8791a",
               "binary_sha256": "x"}, "environment": {"slurm_cpus_per_task": "4"},
               "subset": {"cells": 12, "cells_in_dataset": 12, "genes": 6, "regulators": 3, "subset_seed": 0},
               **result}
    assert smoke.write_summary(summary, out, scrub) == []
    written = json.loads((out / "summary.json").read_text(encoding="utf-8"))
    assert written["privacy_scan"] == "clean"
    assert not [s for s in _strings(written) if str(tmp_path) in s or "someone" in s]
    assert "| t1_a | 1 | ok | 0 |" in (out / "SUMMARY.md").read_text(encoding="utf-8")
    (out / "stray.txt").write_text(f"made by someone in {tmp_path / 'home'}\n")
    assert smoke.privacy_scan(out, scrub) == ["stray.txt"]


# the head of log_grnbench.txt that the real program wrote in EXP-20261018-GB-1 (run t1_a, job 705969)
REAL_LOG = """~/tools/ARACNe3/build/src/app/ARACNe3_app_release -e <tmp>/expr.tsv -r <tmp>/regulators.txt -o <tmp>/out --seed 42 --threads 1 --runid grnbench

---------Sun Oct 18 03:40:36 2026 CEST---------
Beginning ARACNe3 instance...

Gene expression matrix & regulators list read time: 0s

Mutual Information null model calculation time: 5s

Creating subnetwork(s) time: 0s
Total subnetworks generated: 1

Consolidating subnetwork(s) time: 0s
Subnetworks consolidated: 1

Writing final network...
SUCCESS!
"""


def test_log_facts_read_the_real_program_s_log(tmp_path):
    (tmp_path / "log_grnbench.txt").write_text(REAL_LOG, encoding="utf-8")
    facts = smoke.log_facts(tmp_path, "grnbench", "\nTotal N Samples: 500\nSubsampled N Samples: 317\n")
    assert facts["log_times_s"] == {"Gene expression matrix & regulators list read": 0.0,
                                    "Mutual Information null model calculation": 5.0,
                                    "Creating subnetwork(s)": 0.0, "Consolidating subnetwork(s)": 0.0}
    assert facts["log_says_success"] and facts["total_samples"] == 500 and facts["subsampled_samples"] == 317


def test_resummarize_rereads_the_run_folder_and_changes_nothing_else(program, scrub, tmp_path):
    program(FakeProgram())
    out = tmp_path / "out"
    result = smoke.smoke(_adata(), ["TF1", "TF2", "TF3"], [1], 42, out, scrub)
    summary = {"experiment": smoke.EXP, "code_commit": "abc", "aracne3": {"source_commit": "3d8791a",
               "binary_sha256": "x"}, "environment": {"slurm_cpus_per_task": "4"},
               "subset": {"cells": 12, "cells_in_dataset": 12, "genes": 6, "regulators": 3, "subset_seed": 0},
               **result}
    smoke.write_summary(summary, out, scrub)
    before = json.loads((out / "summary.json").read_text(encoding="utf-8"))
    stale = json.loads(json.dumps(before))
    stale["runs"][0]["log_times_s"], stale["runs"][0]["total_samples"] = {}, None   # as the old parser left it
    del stale["runs"][0]["consolidated"]["columns"]                                 # as the run wrote it
    (out / "summary.json").write_text(json.dumps(stale), encoding="utf-8")
    assert smoke.messages_stdout(out / "runs" / "t1_a" / "messages.txt").count("Total N Samples: 12") == 1
    assert smoke.resummarize(out, scrub) == []
    assert json.loads((out / "summary.json").read_text(encoding="utf-8")) == before


def test_scrubber_keeps_numbers_that_contain_the_user_digits():
    s = smoke.Scrubber(home="/home/1234567", user="1234567")
    s.tmp.append("/tmp/job/grnbench_aracne3_x")
    assert s("/home/1234567/tools/ARACNe3 -o /tmp/job/grnbench_aracne3_x/out") == "~/tools/ARACNe3 -o <tmp>/out"
    assert s("user 1234567 mi 0.1234567 big 21234567") == "user <user> mi 0.1234567 big 21234567"
    assert s.leaks("/home/1234567/x") and s.leaks("by 1234567") and not s.leaks("0.1234567")


def test_subset_is_fixed_by_its_seed_and_puts_the_most_variable_tfs_first():
    rng = np.random.default_rng(3)
    counts = rng.poisson(lam=np.linspace(0.2, 6.0, 20), size=(80, 20)).astype(float)
    counts[:, 0] = 0.0                                   # never detected: never chosen
    a = ad.AnnData(counts)
    a.obs_names = [f"c{i}" for i in range(80)]
    a.var_names = [f"TF{i}" if i < 8 else f"G{i}" for i in range(20)]
    sources = {f"tf{i}" for i in range(8)}               # CollecTRI symbols, matched case-insensitively
    sel, tfs, facts = smoke.select_subset(a, sources, 50, 10, 3, 0.1, 0)
    again, tfs_again, _ = smoke.select_subset(a, sources, 50, 10, 3, 0.1, 0)
    assert list(sel.obs_names) == list(again.obs_names) and tfs == tfs_again
    assert sel.shape == (50, 10) and len(tfs) == 3 and "TF0" not in sel.var_names
    assert set(tfs) <= set(sel.var_names) and all(t.startswith("TF") for t in tfs)
    x = sel.X.toarray() if hasattr(sel.X, "toarray") else sel.X
    var = dict(zip(sel.var_names, x.var(axis=0)))
    assert [var[t] for t in tfs] == sorted((var[t] for t in tfs), reverse=True)
    assert facts["cells"] == 50 and facts["regulators"] == 3 and facts["genes"] == 10
