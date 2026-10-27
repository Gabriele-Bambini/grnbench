"""ARACNe3 smoke run: the ``aracne3`` adapter against the real program (EXP-20261018-GB-1, issue #109).

Calls ``grnbench.infer.run_method("aracne3", ...)`` on a small fixed subset of PBMC3k (filtered)
once per entry of ``--threads`` (passed as ``aracne3_threads``, ARACNe3's ``--threads``), with seed
``--seed`` every time, and records what the real program did. Run it only where ARACNe3 is built,
inside a job (``hpc/aracne3_smoke.sbatch``), never on a laptop.

Subset, fixed before the run: ``--cells`` cells drawn without replacement (``--subset-seed``);
``cpm_log1p`` on all genes; of the genes detected in at least ``--min-detected`` of those cells,
the ``--tfs`` most variable CollecTRI regulators (sources of ``references/collectri_human.tsv``)
and the most variable other genes, ``--genes`` in all. The matrix goes in normalized
(``assume_normalized=True``), as the engine hands it to every method; the TFs are the ``tf_list``.

The adapter's ``subprocess.run`` is wrapped, not replaced: the real call runs unchanged, and its
argv, exit status, stdout, stderr, input files, output folder and null-model cache are kept before
the adapter deletes its temporary folder. If the first run keeps no edges (the adapter's "kept no
edges" error), every run is repeated once with ``extra_args=["--alpha", "0.2"]``, the relaxation
that error message names.

Writes into ``--out``: ``summary.json`` and ``SUMMARY.md``, the subset (``subset_cells.txt``,
``subset_genes.txt``, ``regulators.txt``) and per run ``runs/<label>/`` with ``edges.tsv`` (what
``run_method`` returned), ``messages.txt`` and ``aracne3_out/`` (ARACNe3's own files). Texts have
the home folder replaced by ``~``, the adapter's temporary folder by ``<tmp>`` and the user name by
``<user>``; network tables (gene names and numbers) are copied as written. Exit status: 0 every
run returned a network; 1 a run did not; 2 grnbench is not imported from this checkout; 3 a
written file still names the home folder or the user. The summary is written for 0, 1 and 3.
``--resummarize`` rewrites the summary of a finished run folder from its own files (analysis code
fixed after a run, AGENTS.md section 7); it runs nothing.

    python benchmarks/aracne3_smoke.py --data-root ~/grnbench_data --out ~/gb_runs/EXP-20261018-GB-1
    python benchmarks/aracne3_smoke.py --resummarize --out benchmarks/results/EXP-20261018-GB-1
"""
from __future__ import annotations

import argparse
import getpass
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EXP = "EXP-20261018-GB-1"
PBMC3K = "10x_pbmc3k/filtered/filtered_gene_bc_matrices/hg19"
COLLECTRI = "references/collectri_human.tsv"
RELAX = ["--alpha", "0.2"]
NO_EDGES = "kept no edges"
# "<what> time: <seconds>s" lines of ARACNe3's log_<runid>.txt (whole seconds, e.g.
# "Mutual Information null model calculation time: 5s"), sample counts on its stdout
# (src/app/ARACNe3.cpp at 3d8791a)
LOG_TIME = re.compile(r"^(?P<what>[A-Za-z][^:\n]*?) time: (?P<seconds>[0-9.eE+-]+)s?\s*$", re.M)
SAMPLES = re.compile(r"^(Total|Subsampled) N Samples: (\d+)\s*$", re.M)


class Scrubber:
    """Replaces the home folder (as given and resolved) with ``~``, the adapter's temporary folders
    with ``<tmp>`` and the user name with ``<user>``: none of them may reach a committed file."""

    def __init__(self, home: str | None = None, user: str | None = None):
        home = os.path.expanduser("~") if home is None else home
        self.homes = sorted({h for h in (home, os.path.realpath(home)) if h not in ("", "/", "~")},
                            key=len, reverse=True)
        if user is None:
            try:
                user = getpass.getuser()
            except Exception:  # noqa: BLE001 - no user name to hide is not an error
                user = ""
        # a whole word only, and not after a dot: the digits of a decimal are not the user
        self.user = re.compile(rf"(?<![\w.]){re.escape(user)}(?!\w)") if user else None
        self.tmp: list[str] = []

    def __call__(self, text: str) -> str:
        for t in sorted(set(self.tmp), key=len, reverse=True):
            text = text.replace(t, "<tmp>")
        for h in self.homes:
            text = text.replace(h, "~")
        return self.user.sub("<user>", text) if self.user else text

    def leaks(self, text: str) -> bool:
        """The home folder as written or as JSON escapes it (backslashes), or the user name."""
        homes = {v for h in self.homes for v in (h, json.dumps(h)[1:-1])}
        return any(h in text for h in homes) or bool(self.user and self.user.search(text))


def scrub_obj(obj, scrub: Scrubber):
    if isinstance(obj, str):
        return scrub(obj)
    if isinstance(obj, dict):
        return {k: scrub_obj(v, scrub) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [scrub_obj(v, scrub) for v in obj]
    return obj


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def is_network_table(path: Path) -> bool:
    """ARACNe3's network tables (consolidated and per subnetwork): gene names and numbers only."""
    return path.suffix == ".tsv" and (path.name.startswith("consolidated-net_") or path.parent.name == "subnets")


def copy_output(src: Path, dest: Path, scrub: Scrubber) -> None:
    """ARACNe3's output folder: network tables byte for byte, every other text scrubbed."""
    for f in sorted(p for p in src.rglob("*") if p.is_file()):
        to = dest / f.relative_to(src)
        to.parent.mkdir(parents=True, exist_ok=True)
        if is_network_table(f):
            shutil.copyfile(f, to)
        else:
            to.write_text(scrub(f.read_text(encoding="utf-8", errors="replace")), encoding="utf-8")


class Capture:
    """Stands between the adapter and ``subprocess.run``: the real call runs unchanged, and what it
    did is kept before the adapter's temporary folder is deleted."""

    def __init__(self, dest: Path, scrub: Scrubber):
        self.dest, self.scrub, self.calls = dest, scrub, []

    def run(self, cmd, **kwargs):
        work = Path(kwargs["cwd"])
        self.scrub.tmp += [str(work), os.path.realpath(work)]
        argv = [str(a) for a in cmd]
        inputs = {}
        for flag in ("-e", "-r"):
            if flag in argv[:-1]:
                path = Path(argv[argv.index(flag) + 1])
                if path.is_file():
                    with open(path, "rb") as fh:
                        lines = fh.read().splitlines()
                    inputs[flag] = {"sha256": sha256(path), "lines": len(lines),
                                    "first_line_tabs": lines[0].count(b"\t") if lines else 0}
        t0 = time.perf_counter()
        proc = subprocess.run(cmd, **kwargs)  # noqa: S603 - the adapter's own argv, unchanged
        call = {"argv": [self.scrub(a) for a in argv], "returncode": proc.returncode,
                "seconds": round(time.perf_counter() - t0, 3), "inputs": inputs,
                "stdout": self.scrub(proc.stdout or ""), "stderr": self.scrub(proc.stderr or "")}
        if (work / "out").is_dir():
            copy_output(work / "out", self.dest / "aracne3_out", self.scrub)
        cache = work / ".ARACNe3_cached"
        call["null_model_sha256"] = ({p.name: sha256(p) for p in sorted(cache.iterdir()) if p.is_file()}
                                     if cache.is_dir() else {})
        self.calls.append(call)
        return proc


def consolidated_info(out_dir: Path, run_id: str) -> dict:
    """The consolidated file as ARACNe3 wrote it: hash, header, rows, first rows and, when it has the
    six columns, what its columns hold (``mi.values`` is the MI on all cells, while the edges were
    chosen on the subsample)."""
    net = out_dir / f"consolidated-net_{run_id}.tsv"
    if not net.exists():
        return {"file": net.name, "present": False}
    lines = net.read_text(encoding="utf-8", errors="replace").splitlines()
    info = {"file": net.name, "present": True, "sha256": sha256(net),
            "header": lines[0].split("\t") if lines else [],
            "rows": sum(1 for line in lines[1:] if line.strip()),
            "first_rows": [line.split("\t") for line in lines[1:6]]}
    if info["rows"] and info["header"] == ["regulator.values", "target.values", "mi.values",
                                           "scc.values", "count.values", "log.p.values"]:
        t = pd.read_csv(net, sep="\t")
        mi, scc = t["mi.values"].astype(float), t["scc.values"].astype(float)
        info["columns"] = {
            "regulators_with_edges": int(t["regulator.values"].nunique()),
            "mi_min": float(mi.min()), "mi_median": float(mi.median()), "mi_max": float(mi.max()),
            "mi_zero_rows": int((mi == 0).sum()),
            "scc_positive": int((scc > 0).sum()), "scc_negative": int((scc < 0).sum()),
            "scc_zero": int((scc == 0).sum()),
            "count_values": sorted(int(c) for c in t["count.values"].unique()),
            "distinct_log_p_values": int(t["log.p.values"].nunique())}
    return info


def log_facts(out_dir: Path, run_id: str, stdout: str) -> dict:
    """Timings from ARACNe3's log and the sample counts it printed (its reading of the header)."""
    log = out_dir / f"log_{run_id}.txt"
    text = log.read_text(encoding="utf-8", errors="replace") if log.exists() else ""
    facts = {"log_times_s": {m["what"].strip(): float(m["seconds"]) for m in LOG_TIME.finditer(text)},
             "log_says_success": "SUCCESS!" in text}
    for kind, n in SAMPLES.findall(stdout):
        facts[f"{kind.lower()}_samples"] = int(n)
    return facts


def run_once(label: str, adata, tfs: list[str], threads: int, seed: int, extra: list[str],
             out: Path, scrub: Scrubber) -> tuple[dict, pd.DataFrame | None]:
    from grnbench.infer import MethodError, aracne3, run_method

    dest = out / "runs" / label
    shutil.rmtree(dest, ignore_errors=True)
    dest.mkdir(parents=True)
    cap = Capture(dest, scrub)
    # aracne3_threads, not threads: the adapter does not pass the engine's generic count on (review of #106)
    params = {"aracne3_threads": threads, "seed": seed, **({"extra_args": list(extra)} if extra else {})}
    row: dict = {"label": label, "threads": threads, "seed": seed, "extra_args": list(extra)}
    edges = None
    t0 = time.perf_counter()
    try:
        with mock.patch.object(aracne3, "subprocess",
                               SimpleNamespace(run=cap.run, TimeoutExpired=subprocess.TimeoutExpired)):
            edges, prov = run_method("aracne3", adata, tf_list=tfs, assume_normalized=True, **params)
        row.update(status="ok", error=None, n_edges=int(len(edges)), provenance=scrub_obj(prov, scrub),
                   top_edges=json.loads(edges.head(5).to_json(orient="records")))
        edges.to_csv(dest / "edges.tsv", sep="\t", index=False, lineterminator="\n")
    except MethodError as exc:
        row.update(status="error", error=scrub(str(exc)), n_edges=0)
    except Exception as exc:  # noqa: BLE001 - not the adapter's refusal: kept, and the run counts as failed
        row.update(status="crash", error=scrub(f"{type(exc).__name__}: {exc}"), n_edges=0)
    row["seconds"] = round(time.perf_counter() - t0, 3)
    row["program_calls"] = len(cap.calls)
    call = cap.calls[-1] if cap.calls else {}
    row["aracne3"] = {k: v for k, v in call.items() if k not in ("stdout", "stderr")}
    out_dir = dest / "aracne3_out"
    row["consolidated"] = consolidated_info(out_dir, aracne3.RUN_ID)
    row.update(log_facts(out_dir, aracne3.RUN_ID, call.get("stdout", "")))
    row["parsed"] = row["status"] == "ok"
    (dest / "messages.txt").write_text(
        f"exit status: {call.get('returncode')}\n--- stdout ---\n{call.get('stdout', '')}\n"
        f"--- stderr ---\n{call.get('stderr', '')}\n--- adapter ---\n{row['error'] or 'ok'}\n",
        encoding="utf-8")
    return row, edges


def compare(a: dict, b: dict, edges: dict[str, pd.DataFrame | None]) -> dict:
    """Two runs: same consolidated file (bytes), same null model, same edges and weights."""
    ca, cb = a["consolidated"], b["consolidated"]
    na, nb = a["aracne3"].get("null_model_sha256") or {}, b["aracne3"].get("null_model_sha256") or {}
    res = {"a": a["label"], "b": b["label"],
           "same_consolidated_file": bool(ca.get("sha256")) and ca.get("sha256") == cb.get("sha256"),
           "same_null_model": bool(na) and na == nb}
    ea, eb = edges.get(a["label"]), edges.get(b["label"])
    if ea is None or eb is None:
        return {**res, "comparable": False}
    wa = ea.set_index(["source", "target"])
    wb = eb.set_index(["source", "target"])
    shared = wa.index.intersection(wb.index)
    union = wa.index.union(wb.index)
    dw = (wa.loc[shared, "weight"] - wb.loc[shared, "weight"]).abs()
    flips = int((wa.loc[shared, "direction"] != wb.loc[shared, "direction"]).sum()) if "direction" in wa else 0
    biggest = float(dw.max()) if len(shared) else None
    return {**res, "comparable": True, "edges_a": int(len(wa)), "edges_b": int(len(wb)),
            "shared": int(len(shared)), "jaccard": round(len(shared) / len(union), 6) if len(union) else 1.0,
            "max_abs_weight_diff": biggest, "direction_flips": flips,
            "identical_edges": len(shared) == len(union) and not biggest and flips == 0}


def comparison_pairs(rows: list[dict]) -> list[tuple[dict, dict]]:
    """Each run against the first run with its thread count (repeats), and the first run of every
    other thread count against the very first run (threads)."""
    first: dict[int, dict] = {}
    pairs = []
    for r in rows:
        if r["threads"] in first:
            pairs.append((first[r["threads"]], r))
        else:
            if first:
                pairs.append((rows[0], r))
            first[r["threads"]] = r
    return pairs


def select_subset(adata, sources: set[str], n_cells: int, n_genes: int, n_tfs: int,
                  min_detected: float, subset_seed: int):
    """(normalized cells x genes AnnData, regulators most variable first, facts)."""
    from grnbench.strategies import cpm_log1p

    rng = np.random.default_rng(subset_seed)
    cells = np.sort(rng.choice(adata.n_obs, size=min(n_cells, adata.n_obs), replace=False))
    sub = adata[cells].copy()
    detected = np.asarray((sub.X > 0).sum(axis=0)).ravel() / sub.n_obs
    norm = cpm_log1p(sub)
    x = norm.X.toarray() if hasattr(norm.X, "toarray") else np.asarray(norm.X)
    var = x.var(axis=0)
    names = np.asarray(norm.var_names.astype(str))
    upper = {s.upper() for s in sources}
    is_tf = np.array([n.upper() in upper for n in names])
    order = np.argsort(-var, kind="stable")              # most variable first; ties by position
    keep = detected >= min_detected
    tf_idx = [int(i) for i in order if keep[i] and is_tf[i]][:n_tfs]
    other = [int(i) for i in order if keep[i] and not is_tf[i]][:max(n_genes - len(tf_idx), 0)]
    genes = sorted(tf_idx + other)
    facts = {"cells_in_dataset": int(adata.n_obs), "genes_in_dataset": int(adata.n_vars),
             "cells": int(len(cells)), "genes_detected": int(keep.sum()),
             "regulators_detected": int((keep & is_tf).sum()), "genes": int(len(genes)),
             "regulators": int(len(tf_idx))}
    return norm[:, genes].copy(), [str(names[i]) for i in tf_idx], facts


def input_files(data_root: Path, rels: list[str], registry: Path) -> list[dict]:
    """sha256 of each input and whether it is the file registered in experiments/datasets.tsv."""
    reg = pd.read_csv(registry, sep="\t", dtype=str) if registry.exists() else pd.DataFrame()
    known = dict(zip(reg.get("path", []), reg.get("sha256", [])))
    rows = []
    for rel in rels:
        p = data_root / rel
        digest = sha256(p)
        rows.append({"path": rel, "bytes": p.stat().st_size, "sha256": digest,
                     "registered_sha256": known.get(rel), "matches_registry": known.get(rel) == digest})
    return rows


def write_markdown(summary: dict, path: Path) -> None:
    s = summary
    lines = [f"# {s['experiment']}: ARACNe3 smoke run (written by benchmarks/aracne3_smoke.py)", "",
             f"- code: `{s['code_commit']}`; ARACNe3 source `{s['aracne3']['source_commit']}`, "
             f"program sha256 `{s['aracne3']['binary_sha256']}`",
             f"- subset: {s['subset']['cells']} of {s['subset']['cells_in_dataset']} cells, "
             f"{s['subset']['genes']} genes, {s['subset']['regulators']} regulators "
             f"(subset seed {s['subset']['subset_seed']}); CPUs given to the job: "
             f"{s['environment']['slurm_cpus_per_task']}",
             f"- relaxed rerun (`{' '.join(RELAX)}`): {'yes' if s['relaxed'] else 'no'}; "
             f"private-text scan: {s['privacy_scan']}", "",
             "| run | threads | status | exit | seconds | samples read | edges | consolidated rows | null model sha256 |",
             "|---|---|---|---|---|---|---|---|---|"]
    for r in s["runs"]:
        null = ",".join(v[:12] for v in (r["aracne3"].get("null_model_sha256") or {}).values()) or "-"
        lines.append(f"| {r['label']} | {r['threads']} | {r['status']} | {r['aracne3'].get('returncode')} | "
                     f"{r['seconds']} | {r.get('total_samples', '-')} | {r['n_edges']} | "
                     f"{r['consolidated'].get('rows', '-')} | {null} |")
    lines += ["", "| a | b | same file | same null model | edges a / b | shared | Jaccard | "
              "max abs weight diff | identical edges |", "|---|---|---|---|---|---|---|---|---|"]
    for c in s["comparisons"]:
        lines.append(f"| {c['a']} | {c['b']} | {c['same_consolidated_file']} | {c['same_null_model']} | "
                     f"{c.get('edges_a', '-')} / {c.get('edges_b', '-')} | {c.get('shared', '-')} | "
                     f"{c.get('jaccard', '-')} | {c.get('max_abs_weight_diff', '-')} | "
                     f"{c.get('identical_edges', '-')} |")
    lines += ["", "| run | regulators with edges | MI min / median / max | rows with MI 0 | "
              "Spearman + / - / 0 | count.values | distinct log.p.values |", "|---|---|---|---|---|---|---|"]
    for r in s["runs"]:
        c = r["consolidated"].get("columns")
        if c:
            lines.append(f"| {r['label']} | {c['regulators_with_edges']} | {c['mi_min']} / {c['mi_median']} / "
                         f"{c['mi_max']} | {c['mi_zero_rows']} | {c['scc_positive']} / {c['scc_negative']} / "
                         f"{c['scc_zero']} | {c['count_values']} | {c['distinct_log_p_values']} |")
    for r in s["runs"]:
        if r["status"] != "ok":
            lines += ["", f"{r['label']}: {r['error']}"]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def privacy_scan(out: Path, scrub: Scrubber) -> list[str]:
    return [f.relative_to(out).as_posix() for f in sorted(p for p in out.rglob("*") if p.is_file())
            if scrub.leaks(f.read_text(encoding="utf-8", errors="replace"))]


def program_facts(src: Path, scrub: Scrubber) -> dict:
    from grnbench.infer import aracne3

    binary, why = aracne3._find_binary()

    def git(*args: str) -> str:
        try:
            return subprocess.run(["git", "-C", str(src), *args], capture_output=True, text=True,
                                  timeout=60).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            return ""
    return {"binary": scrub(binary or ""), "unavailable_reason": why,
            "binary_sha256": sha256(Path(binary)) if binary else None,
            "source_commit": git("rev-parse", "HEAD") or None,
            "submodules": git("submodule", "status").splitlines()}


def run_labels(threads: list[int]) -> list[str]:
    seen: dict[int, int] = {}
    labels = []
    for t in threads:
        seen[t] = seen.get(t, 0) + 1
        labels.append(f"t{t}_{chr(96 + seen[t])}")
    return labels


def smoke(sel, tfs: list[str], threads: list[int], seed: int, out: Path, scrub: Scrubber) -> dict:
    """Every run (and the relaxed repeat when the first keeps no edges), then the comparisons."""
    labels = run_labels(threads)

    def series(extra: list[str], suffix: str) -> tuple[list[dict], dict]:
        rows, edges = [], {}
        for label, t in zip(labels, threads):
            row, e = run_once(label + suffix, sel, tfs, t, seed, extra, out, scrub)
            rows.append(row)
            edges[row["label"]] = e
            print(f"{row['label']}: {row['status']} exit {row['aracne3'].get('returncode')} "
                  f"{row['n_edges']} edges {row['seconds']} s", flush=True)
        return rows, edges

    rows, edges = series([], "")
    comparisons = [compare(a, b, edges) for a, b in comparison_pairs(rows)]
    relaxed = rows[0]["status"] != "ok" and NO_EDGES in (rows[0]["error"] or "")
    if relaxed:
        more, more_edges = series(RELAX, "_alpha0.2")
        comparisons += [compare(a, b, more_edges) for a, b in comparison_pairs(more)]
        rows += more
    return {"seed": seed, "relaxed": relaxed, "runs": rows, "comparisons": comparisons,
            "verdict": {"every_run_returned_a_network": all(r["status"] == "ok" for r in rows),
                        "identical_edges": {f"{c['a']}~{c['b']}": c.get("identical_edges")
                                            for c in comparisons}}}


def write_summary(summary: dict, out: Path, scrub: Scrubber) -> list[str]:
    """summary.json and SUMMARY.md, then the scan of every written file (its result goes in both)."""
    def write() -> None:
        clean = scrub_obj(summary, scrub)
        (out / "summary.json").write_text(json.dumps(clean, indent=2) + "\n", encoding="utf-8")
        write_markdown(clean, out / "SUMMARY.md")

    summary["privacy_scan"] = "pending"
    write()
    bad = privacy_scan(out, scrub)
    summary["privacy_scan"] = bad or "clean"
    write()
    return bad


def messages_stdout(path: Path) -> str:
    """ARACNe3's stdout as ``run_once`` kept it in ``messages.txt``."""
    text = path.read_text(encoding="utf-8") if path.exists() else ""
    start, end = text.find("--- stdout ---\n"), text.find("\n--- stderr ---")
    return text[start + len("--- stdout ---\n"):end] if start >= 0 and end > start else ""


def resummarize(out: Path, scrub: Scrubber) -> list[str]:
    """Rewrite ``summary.json`` and ``SUMMARY.md`` of a finished run from its own files: the log
    facts of every run are read again from ``runs/<label>/``; nothing runs, nothing else changes."""
    from grnbench.infer import aracne3

    summary = json.loads((out / "summary.json").read_text(encoding="utf-8"))
    for r in summary["runs"]:
        run = out / "runs" / r["label"]
        r["consolidated"] = consolidated_info(run / "aracne3_out", aracne3.RUN_ID)
        r.update(log_facts(run / "aracne3_out", aracne3.RUN_ID, messages_stdout(run / "messages.txt")))
    return write_summary(summary, out, scrub)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data-root", type=Path, default=Path("~/grnbench_data"))
    ap.add_argument("--registry", type=Path, default=ROOT / "experiments" / "datasets.tsv")
    ap.add_argument("--aracne3-src", type=Path, default=Path("~/tools/ARACNe3"))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--cells", type=int, default=500)
    ap.add_argument("--genes", type=int, default=300)
    ap.add_argument("--tfs", type=int, default=30)
    ap.add_argument("--min-detected", type=float, default=0.1)
    ap.add_argument("--subset-seed", type=int, default=0)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--threads", default="1,1,4,4", help="one run per entry, e.g. 1,1,4,4")
    ap.add_argument("--resummarize", action="store_true",
                    help="only rewrite summary.json and SUMMARY.md of the finished run in --out")
    args = ap.parse_args(argv)
    if args.resummarize:
        bad = resummarize(args.out.expanduser(), Scrubber())
        return 3 if bad else 0

    import anndata
    import scanpy

    import grnbench
    from grnbench.ingest import load_any

    scrub = Scrubber()
    here = Path(grnbench.__file__).resolve().parent
    if here != (ROOT / "src" / "grnbench").resolve():
        print(f"grnbench is imported from {scrub(str(here))}, not from this checkout: set PYTHONPATH=src",
              file=sys.stderr)
        return 2
    data_root, out = args.data_root.expanduser(), args.out.expanduser()
    out.mkdir(parents=True, exist_ok=True)
    adata, _ = load_any(data_root / PBMC3K)
    sources = set(pd.read_csv(data_root / COLLECTRI, sep="\t")["source"].astype(str))
    sel, tfs, facts = select_subset(adata, sources, args.cells, args.genes, args.tfs,
                                    args.min_detected, args.subset_seed)
    (out / "subset_cells.txt").write_text("".join(f"{c}\n" for c in sel.obs_names), encoding="utf-8")
    (out / "subset_genes.txt").write_text("".join(f"{g}\n" for g in sel.var_names), encoding="utf-8")
    (out / "regulators.txt").write_text("".join(f"{t}\n" for t in tfs), encoding="utf-8")

    threads = [int(t) for t in str(args.threads).split(",") if t.strip()]
    result = smoke(sel, tfs, threads, args.seed, out, scrub)
    commit_file = ROOT / "GIT_COMMIT"
    rels = [f"{PBMC3K}/{n}" for n in ("matrix.mtx", "genes.tsv", "barcodes.tsv")] + [COLLECTRI]
    summary = {
        "experiment": EXP,
        "code_commit": commit_file.read_text(encoding="utf-8").strip() if commit_file.exists() else "unknown",
        "grnbench_imported_from": here.relative_to(ROOT.resolve()).as_posix(),
        "aracne3": program_facts(args.aracne3_src.expanduser(), scrub),
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
                        "pandas": pd.__version__, "anndata": anndata.__version__,
                        "scanpy": scanpy.__version__, "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
                        "slurm_cpus_per_task": os.environ.get("SLURM_CPUS_PER_TASK"),
                        "omp_num_threads": os.environ.get("OMP_NUM_THREADS")},
        "inputs": input_files(data_root, rels, args.registry),
        "subset": {**facts, "subset_seed": args.subset_seed, "min_detected": args.min_detected,
                   "requested": {"cells": args.cells, "genes": args.genes, "tfs": args.tfs}},
        **result,
    }
    bad = write_summary(summary, out, scrub)
    if bad:
        print(f"home folder or user name still in: {', '.join(bad)}", file=sys.stderr)
        return 3
    return 0 if result["verdict"]["every_run_returned_a_network"] else 1


if __name__ == "__main__":
    sys.exit(main())
