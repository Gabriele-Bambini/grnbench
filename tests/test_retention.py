"""Retention ledger: every cleaning step accounts for the cells, genes and named TFs it removes.

Accounting and warnings only: which genes and cells are kept must not change.
"""
from __future__ import annotations

import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from grnbench.diagnostics import read_diagnostics, write_diagnostics
from grnbench.pipeline import run_prep

PREFIX = "TF protection:"


def _write(tmp_path, extra: dict[str, np.ndarray] | None = None, *, n_cells: int = 300,
           n_genes: int = 300, obs: dict | None = None, seed: int = 0, name: str = "in.h5ad"):
    """Poisson(3) counts (every cell passes the 200-gene floor) plus named extra columns."""
    rng = np.random.default_rng(seed)
    X = rng.poisson(3.0, size=(n_cells, n_genes)).astype("float32")
    names = [f"G{i:03d}" for i in range(n_genes)]
    cols = [X]
    for gene, column in (extra or {}).items():
        cols.append(np.asarray(column, dtype="float32").reshape(-1, 1))
        names.append(gene)
    frame = pd.DataFrame(obs or {}, index=[f"c{i}" for i in range(n_cells)])
    path = tmp_path / name
    ad.AnnData(X=np.hstack(cols), obs=frame, var=pd.DataFrame(index=names)).write_h5ad(path)
    return path


def _in_cells(n_cells: int, cells, value: float = 2.0) -> np.ndarray:
    col = np.zeros(n_cells, dtype="float32")
    col[list(cells)] = value
    return col


def _row(manifest: dict, step: str) -> dict:
    rows = [r for r in manifest["retention"] if r["step"] == step]
    assert len(rows) == 1, (step, [r["step"] for r in manifest["retention"]])
    return rows[0]


def _tf_warnings(manifest: dict) -> list[str]:
    return [w for w in manifest["warnings"] if w.startswith(PREFIX)]


def test_gene_filter_warns_about_a_detected_tf_and_keeps_the_genes_unchanged(tmp_path):
    n = 300
    rng = np.random.default_rng(3)
    path = _write(tmp_path, {
        "TFDET": _in_cells(n, [5]),                         # detected in one cell: < 3 cells
        "TFZERO": np.zeros(n),                              # never detected
        "TFKEEP": rng.poisson(3.0, n),
    })
    res = run_prep(path, tmp_path / "out", strategy="cpm_log1p", force=True, seed=1,
                   tf_list=["TFDET", "TFZERO", "TFKEEP", "NOT_IN_DATA"])
    on_disk = json.loads(res.manifest_path.read_text(encoding="utf-8"))
    row = _row(on_disk, "filter_genes")
    assert (row["genes_before"], row["genes_after"], row["genes_removed"]) == (303, 301, 2)
    assert (row["cells_before"], row["cells_after"], row["cells_removed"]) == (300, 300, 0)
    assert (row["tfs_before"], row["tfs_after"], row["tfs_removed"]) == (3, 1, 2)
    assert row["tfs_removed_names"] == ["TFDET", "TFZERO"]
    assert row["tfs_removed_detected"] == 1
    assert row["tfs_removed_detected_names"] == ["TFDET"]
    for key in ("genes_before", "tfs_removed", "tfs_removed_detected", "cells_removed"):
        assert type(row[key]) is int, key                  # never a numpy scalar written as a string

    warns = _tf_warnings(on_disk)
    assert len(warns) == 1, on_disk["warnings"]
    assert "'filter_genes'" in warns[0] and "1 TF" in warns[0] and "TFDET (1 cell)" in warns[0]
    assert "TFZERO" not in warns[0], "an undetected TF carries no signal and is not a loss"

    processed = ad.read_h5ad(res.processed_path)
    assert "TFDET" not in processed.var_names and "TFZERO" not in processed.var_names
    assert processed.n_vars == 301, "accounting only: the gene filter keeps exactly what it kept"


def test_steps_that_remove_nothing_record_zeros_in_step_order(tmp_path):
    path = _write(tmp_path, {"TFKEEP": np.random.default_rng(1).poisson(3.0, 300)})
    res = run_prep(path, tmp_path / "out", strategy="cpm_log1p", force=True, seed=1,
                   tf_list=["TFKEEP"], subsample_cells=120)
    manifest = res.manifest
    ledger_steps = [r["step"] for r in manifest["retention"]]
    assert ledger_steps == ["filter_cells", "filter_genes", "subsample_cells", "strategy"]
    assert ledger_steps == [s["step"] for s in manifest["steps"] if s["step"] in ledger_steps]
    for step in ("filter_cells", "filter_genes", "strategy"):
        row = _row(manifest, step)
        assert row["cells_removed"] == 0 and row["genes_removed"] == 0, row
        assert row["tfs_removed"] == 0 and row["tfs_removed_names"] == [], row
        assert row["tfs_before"] == row["tfs_after"] == 1
        assert row["aggregated"] is False
    sub = _row(manifest, "subsample_cells")
    assert (sub["cells_before"], sub["cells_after"], sub["cells_removed"]) == (300, 120, 180)
    assert _tf_warnings(manifest) == []


def test_without_a_tf_list_the_ledger_counts_cells_and_genes_only(tmp_path):
    path = _write(tmp_path, {"ZERO": np.zeros(300)})
    res = run_prep(path, tmp_path / "out", strategy="cpm_log1p", force=True, seed=1)
    row = _row(res.manifest, "filter_genes")
    assert row["genes_removed"] == 1
    assert not any(key.startswith("tfs_") for key in row), row
    assert _tf_warnings(res.manifest) == []


def test_detection_is_judged_in_the_context_and_both_gene_filters_are_checked(tmp_path):
    n = 300
    ct = np.array(["A"] * 150 + ["B"] * 150)
    rng = np.random.default_rng(5)
    # TFGLOBAL: 2 cells of A -> under the global 3-cell floor; detected in A, absent from B.
    # TFCTX: 2 cells of A + 5 of B -> passes the global floor, then fails the floor inside A.
    extra = {"TFGLOBAL": _in_cells(n, [0, 1]), "TFCTX": _in_cells(n, [2, 3, 200, 201, 202, 203, 204]),
             "TFKEEP": rng.poisson(3.0, n)}
    path = _write(tmp_path, extra, obs={"ct": ct})
    tfs = ["TFGLOBAL", "TFCTX", "TFKEEP"]

    in_a = run_prep(path, tmp_path / "a", strategy="cpm_log1p", force=True, seed=1, tf_list=tfs,
                    subset_key="ct", subset_value="A").manifest
    glob = _row(in_a, "filter_genes")
    assert glob["tfs_removed_names"] == ["TFGLOBAL"] and glob["tfs_removed_detected"] == 1
    ctx = _row(in_a, "filter_genes_in_context")
    assert ctx["tfs_removed_names"] == ["TFCTX"] and ctx["tfs_removed_detected"] == 1
    sub = _row(in_a, "subset_context")
    assert (sub["cells_before"], sub["cells_after"], sub["cells_removed"]) == (300, 150, 150)
    warns = _tf_warnings(in_a)
    assert len(warns) == 2, in_a["warnings"]
    # the count is the context's cells: TFCTX is non-zero in 7 cells overall, 2 of them in A
    assert "'filter_genes'" in warns[0] and "TFGLOBAL (2 cells)" in warns[0]
    assert "'filter_genes_in_context'" in warns[1] and "TFCTX (2 cells)" in warns[1]

    in_b = run_prep(path, tmp_path / "b", strategy="cpm_log1p", force=True, seed=1, tf_list=tfs,
                    subset_key="ct", subset_value="B").manifest
    glob_b = _row(in_b, "filter_genes")
    assert glob_b["tfs_removed_names"] == ["TFGLOBAL"], "the global filter removes it in every context"
    assert glob_b["tfs_removed_detected"] == 0, "never seen in B: removing it loses nothing there"
    assert _row(in_b, "filter_genes_in_context")["tfs_removed"] == 0
    assert _tf_warnings(in_b) == []


def test_qc_gene_exclusion_is_checked_like_any_gene_filter(tmp_path):
    rng = np.random.default_rng(2)
    path = _write(tmp_path, {"RPL22": rng.poisson(3.0, 300), "TFKEEP": rng.poisson(3.0, 300)})
    res = run_prep(path, tmp_path / "out", strategy="cpm_log1p", force=True, seed=1,
                   tf_list=["RPL22", "TFKEEP"], exclude_qc_genes=True)
    row = _row(res.manifest, "exclude_qc_genes")
    assert row["genes_removed"] == 1 and row["tfs_removed_names"] == ["RPL22"]
    warns = _tf_warnings(res.manifest)
    assert len(warns) == 1 and "'exclude_qc_genes'" in warns[0]
    n_cells = int((ad.read_h5ad(path)[:, "RPL22"].X > 0).sum())
    assert f"RPL22 ({n_cells} cells)" in warns[0]
    assert "RPL22" not in ad.read_h5ad(res.processed_path).var_names


def test_names_are_capped_sorted_and_counted(tmp_path):
    n = 300
    rng = np.random.default_rng(4)
    extra = {f"TFZ{i:02d}": np.zeros(n) for i in reversed(range(25))}      # undetected
    extra.update({f"TFD{i}": _in_cells(n, [i]) for i in reversed(range(7))})  # detected once
    extra["TFKEEP"] = rng.poisson(3.0, n)
    path = _write(tmp_path, extra)
    res = run_prep(path, tmp_path / "out", strategy="cpm_log1p", force=True, seed=1,
                   tf_list=list(extra))
    row = _row(res.manifest, "filter_genes")
    assert row["tfs_removed"] == 32 and row["tfs_removed_detected"] == 7
    expected = sorted(list(extra)[:-1], key=str.casefold)
    assert row["tfs_removed_names"] == expected[:20]
    assert row["tfs_removed_detected_names"] == [f"TFD{i}" for i in range(7)]
    warn = _tf_warnings(res.manifest)[0]
    assert "7 TF" in warn and "TFD5" not in warn and "and 2 more" in warn
    assert ", ".join(f"TFD{i} (1 cell)" for i in range(5)) in warn


def test_feature_selection_and_aggregation_rows(tmp_path):
    n = 300
    rng = np.random.default_rng(6)
    path = _write(tmp_path, {"TFKEEP": rng.poisson(3.0, n)}, obs={"grp": [f"g{i % 3}" for i in range(n)]})
    fs = run_prep(path, tmp_path / "fs", strategy="cpm_log1p", force=True, seed=1,
                  tf_list=["TFKEEP"], n_hvg=40).manifest
    row = _row(fs, "feature_selection")
    assert row["genes_before"] == 301 and row["genes_after"] == fs["feature_selection"]["n_features_kept"]
    assert row["tfs_removed"] == 0, "D-003: feature selection never drops a TF in the data"

    pb = run_prep(path, tmp_path / "pb", strategy="pseudobulk", group_key="grp", force=True,
                  seed=1, tf_list=["TFKEEP"]).manifest
    row = _row(pb, "strategy")
    assert row["aggregated"] is True and row["cells_after"] == 3
    assert row["cells_removed"] is None, "aggregating cells into groups is not removing them"


def test_portable_diagnostics_and_reports_carry_the_ledger(tmp_path):
    n = 300
    path = _write(tmp_path, {"TFDET": _in_cells(n, [5]),
                             "TFKEEP": np.random.default_rng(1).poisson(3.0, n)})
    res = run_prep(path, tmp_path / "out", strategy="cpm_log1p", force=True, seed=1,
                   tf_list=["TFDET", "TFKEEP"])
    doc = read_diagnostics(res.diagnostics_path)
    entry = doc["prep"]["default"]["all_cells"]
    rows = {r["stage"]: r for r in entry["retention"]}
    manifest_row = _row(res.manifest, "filter_genes")
    assert rows["filter_genes"]["genes_before"] == manifest_row["genes_before"]
    assert rows["filter_genes"]["tfs_removed_names"] == ["TFDET"]

    report = res.report_path.read_text(encoding="utf-8")
    assert "Retention ledger" in report
    line = next(x for x in report.splitlines() if x.startswith("| filter_genes |"))
    assert "302 -> 301 (-1)" in line and "TFDET" in line
    cells = next(x for x in report.splitlines() if x.startswith("| filter_cells |"))
    assert "| 300 -> 300 |" in cells and "(-0)" not in report, "a zero removal shows no '(-0)'"

    from grnbench.html_report import render_html
    from grnbench import texreport
    html = render_html({"gold": "x", "gold_pairs": 0, "methods": {}, "overlap": {}}, diagnostics=doc)
    assert "Retention ledger" in html and "TFDET" in html and "302 -&gt; 301" in html
    tex = texreport._diagnostics_tex(doc)
    assert "TFDET" in tex
    assert r"\begin{tabular}{llllp{" in tex, "both TF name columns wrap (p{}), they hold up to 20 names"

    # generated, not typed: change the number in the diagnostic, the report follows
    entry["retention"][[r["stage"] for r in entry["retention"]].index("filter_genes")]["genes_before"] = 999
    from grnbench.diagnostics import render_diagnostics_markdown
    assert "999 -> 301" in render_diagnostics_markdown(doc)


def test_a_long_gene_name_cut_inside_a_token_does_not_fail_the_run(tmp_path):
    """Review of #177, B1: the 80-character cut used to end inside an A/C/G/T run, the reader's
    second sanitising pass then saw a barcode, and run_prep failed after writing its outputs."""
    long = "X" * 60 + "-" + "ACGT" * 5 + "ZZZZZ"
    path = _write(tmp_path, {long: _in_cells(300, [4]),
                             "TFKEEP": np.random.default_rng(0).poisson(3.0, 300)})
    res = run_prep(path, tmp_path / "out", strategy="cpm_log1p", force=True, seed=1,
                   tf_list=[long, "TFKEEP"])
    assert res.report_path.exists()
    assert _row(res.manifest, "filter_genes")["tfs_removed_names"] == [long], "the manifest keeps the name"
    portable = read_diagnostics(res.diagnostics_path)["prep"]["default"]["all_cells"]["retention"]
    name = next(r for r in portable if r["stage"] == "filter_genes")["tfs_removed_names"][0]
    # cut at 80 -> '...-ACGTACGTACGTACGTAC…'; the 18-letter run before the '…' is then a barcode
    assert name == "X" * 60 + "-<barcode>…", name
    assert "Retention ledger" in res.report_path.read_text(encoding="utf-8")


@pytest.mark.parametrize("mutate,match", [
    (lambda r: r.update(unexpected=1), "unsupported fields"),
    (lambda r: r.update(genes_before=-1), "non-negative integer"),
    (lambda r: r.update(tfs_removed_names=["AAACCCAAGAAACCAT-1"]), "sanitised"),
    (lambda r: r.update(tfs_removed_names=[f"TF{i}" for i in range(21)]), "at most"),
])
def test_reader_validates_ledger_rows_strictly(tmp_path, mutate, match):
    from grnbench.diagnostics import new_document
    row = {"stage": "filter_genes", "cells_before": 10, "cells_after": 10, "cells_removed": 0,
           "genes_before": 5, "genes_after": 4, "genes_removed": 1, "aggregated": False,
           "tfs_before": 2, "tfs_after": 1, "tfs_removed": 1, "tfs_removed_detected": 0,
           "tfs_removed_names": ["TFZERO"], "tfs_removed_detected_names": []}
    doc = new_document(outcome="ok")
    doc["prep"] = {"default": {"all_cells": {
        "outcome": "ok", "failure": None, "choices": {}, "steps": [],
        "dimensions": {"before": None, "after": None}, "retention": [row]}}}
    write_diagnostics(tmp_path / "ok.json", doc)
    assert read_diagnostics(tmp_path / "ok.json")["prep"]["default"]["all_cells"]["retention"] == [row]
    mutate(row)
    (tmp_path / "bad.json").write_text(json.dumps(doc), encoding="utf-8")
    with pytest.raises(ValueError, match=match):
        read_diagnostics(tmp_path / "bad.json")
