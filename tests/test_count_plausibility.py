"""Issue #42: integer input that is not raw counts, and specific advice when ingest refuses.

``looks_like_counts`` proves integrality only, so rounded normalize_total output and a
binarized 0/1 matrix used to go through QC, doublets and every model as if they were
counts. The plausibility check warns loudly and records its signals (it never refuses
until its thresholds are calibrated on real data); refusals now say what the data look
like and which other matrices in the file might hold counts.
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.io
import scipy.sparse as sp

from grnbench.ingest import (
    IngestError,
    _validate_counts_layer,
    count_plausibility,
    load_any,
    plausibility_warning,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
import count_plausibility as cp_script  # noqa: E402


def _umi(n_cells=400, n_genes=300, seed=0, depth=3000.0, depth_sd=0.5):
    """Raw UMI counts: gamma-Poisson per gene with a lognormal depth per cell."""
    rng = np.random.default_rng(seed)
    mu = rng.lognormal(0.0, 1.2, n_genes)
    mu /= mu.sum()
    lib = rng.lognormal(np.log(depth), depth_sd, n_cells)
    return rng.poisson(rng.gamma(2.0, mu[None, :] * lib[:, None] / 2.0)).astype(np.float32)


def _full_length(seed=0):
    """Full-length-like read counts: variable depth, molecules amplified, almost no 1s."""
    rng = np.random.default_rng(seed)
    mol = _umi(seed=seed, depth=5000.0, depth_sd=0.6)
    return rng.poisson(mol * rng.gamma(4.0, 30.0 / 4.0, size=mol.shape)).astype(np.float32)


def _normalized(raw, target=1e4):
    return raw / raw.sum(axis=1, keepdims=True) * target


def _write(tmp_path, X, name="in.h5ad", **layers):
    a = ad.AnnData(X=sp.csr_matrix(np.asarray(X, dtype=np.float32)),
                   var=pd.DataFrame(index=[f"G{i}" for i in range(X.shape[1])]))
    for key, value in layers.items():
        a.layers[key] = sp.csr_matrix(np.asarray(value, dtype=np.float32))
    p = tmp_path / name
    p.parent.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(p)
    return p


# ------------------------------------------------------------ the signals #

def test_rounded_normalize_total_is_flagged():
    sig = count_plausibility(np.round(_normalized(_umi())))
    assert sig["cv_cell_totals"] < 0.05
    assert "near_identical_totals" in sig["flags"] and sig["warn"]


def test_rounded_cpm_has_no_ones_and_is_flagged():
    sig = count_plausibility(sp.csr_matrix(np.round(_normalized(_umi(), 1e6))))
    assert {"near_identical_totals", "no_ones"} <= set(sig["flags"]) and sig["warn"]
    assert "no non-zero value is 1" in plausibility_warning(sig)


def test_normalized_then_gene_subset_then_rounded_is_flagged_by_the_combination():
    """Totals are no longer identical once genes are dropped, but still tight, and the
    rounding skipped 1: neither signal warns alone, together they do."""
    raw = _umi()
    keep = np.argsort(-raw.sum(axis=0))[10:]
    sig = count_plausibility(np.round(_normalized(raw, 1e6)[:, keep]))
    assert set(sig["flags"]) == {"tight_totals", "no_ones"} and sig["warn"]
    assert "unusually tight" in plausibility_warning(sig)


def test_binarized_matrix_is_flagged():
    sig = count_plausibility(sp.csr_matrix((_umi() > 0).astype(np.float32)))
    assert "binary" in sig["flags"] and sig["warn"]
    assert "binarized" in plausibility_warning(sig)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_simulated_raw_umi_counts_are_not_flagged(seed):
    sig = count_plausibility(sp.csr_matrix(_umi(seed=seed)))
    assert sig["cv_cell_totals"] > 0.3 and 0 < sig["share_ones"] < 1
    assert sig["flags"] == [] and not sig["warn"]
    assert plausibility_warning(sig) is None


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_full_length_counts_with_few_ones_are_not_flagged(seed):
    """Few 1s is normal for full-length protocols: alone it must never warn."""
    sig = count_plausibility(_full_length(seed))
    assert sig["share_ones"] < 0.01 and sig["max_value"] > 1000
    assert "near_identical_totals" not in sig["flags"] and not sig["warn"]


def test_tiny_matrices_are_not_judged():
    sig = count_plausibility(np.ones((5, 4)))
    assert not sig["checked_totals"] and not sig["checked_values"]
    assert sig["flags"] == [] and not sig["warn"]


def test_signals_record_the_provisional_thresholds():
    sig = count_plausibility(_umi())
    assert "provisional" in sig["thresholds"]["status"]
    assert sig["thresholds"]["totals_cv_near_identical"] < sig["thresholds"]["totals_cv_tight"]


# ---------------------------------------------- warning, report, manifest #

def test_ingest_warns_loudly_and_records_but_does_not_refuse(tmp_path):
    raw = _umi()
    p = _write(tmp_path, np.round(_normalized(raw)))
    with pytest.warns(RuntimeWarning, match="COUNTS LOOK NON-RAW"):
        adata, rep = load_any(p)
    assert adata.shape == raw.shape and "counts" in adata.layers  # loaded, not refused
    assert rep.count_plausibility["warn"] is True
    assert "near_identical_totals" in rep.count_plausibility["flags"]
    assert len(rep.warnings) == 1 and "issue #42" in rep.warnings[0]
    assert rep.to_dict()["count_plausibility"]["flags"] == rep.count_plausibility["flags"]
    stored = adata.uns["grnbench_ingest"]
    assert "COUNTS LOOK NON-RAW" in stored["warnings"][0]


def test_raw_counts_load_without_warning(tmp_path):
    p = _write(tmp_path, _umi())
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        _, rep = load_any(p)
    assert rep.warnings == [] and rep.count_plausibility["flags"] == []
    assert rep.count_plausibility["cv_cell_totals"] > 0.3


def test_finding_reaches_manifest_and_report(tmp_path):
    from grnbench.pipeline import run_prep

    raw = _umi(n_cells=300, n_genes=400)
    src = _write(tmp_path, np.round(_normalized(raw)))
    with pytest.warns(RuntimeWarning, match="COUNTS LOOK NON-RAW"):
        res = run_prep(src, tmp_path / "run", force=True)
    manifest = json.loads(res.manifest_path.read_text(encoding="utf-8"))
    assert "near_identical_totals" in manifest["ingest"]["count_plausibility"]["flags"]
    assert manifest["ingest"]["count_plausibility"]["thresholds"]["totals_cv_near_identical"] > 0
    assert any(w.startswith("COUNTS LOOK NON-RAW") for w in manifest["warnings"])
    report = res.report_path.read_text(encoding="utf-8")
    assert "COUNTS LOOK NON-RAW" in report
    assert "Count plausibility:" in report and "near_identical_totals" in report


def test_engine_plan_carries_the_ingest_warning(tmp_path):
    from grnbench.autoprep import diagnose, plan_preprocessing

    src = _write(tmp_path, np.round(_normalized(_umi())))
    with pytest.warns(RuntimeWarning, match="COUNTS LOOK NON-RAW"):
        diag = diagnose(src)
    assert diag.ingest["warnings"]
    plan = plan_preprocessing(diag)
    assert any(w.startswith("COUNTS LOOK NON-RAW") for w in plan.warnings)


# ---------------------------------------------------------- refusal advice #

def test_refusal_names_log_normalized_data(tmp_path):
    p = _write(tmp_path, np.log1p(_normalized(_umi())))
    with pytest.raises(IngestError, match="log-normalized") as exc:
        load_any(p)
    msg = str(exc.value)
    assert msg.startswith("Input does not look like raw counts")
    assert "re-export the dataset with raw counts" in msg
    assert "assume_counts" not in msg  # the option it used to advertise does not exist


def test_refusal_names_cpm_like_data(tmp_path):
    p = _write(tmp_path, _normalized(_umi(), 1e6))
    with pytest.raises(IngestError, match="normalize_total or CPM"):
        load_any(p)


def test_refusal_names_scaled_data(tmp_path):
    logn = np.log1p(_normalized(_umi()))
    p = _write(tmp_path, (logn - logn.mean(axis=0)) / (logn.std(axis=0) + 1e-9))
    with pytest.raises(IngestError, match="scaled or centered"):
        load_any(p)


def test_refusal_of_other_non_integer_values_stays_generic(tmp_path):
    rng = np.random.default_rng(3)
    tpm = _normalized(_umi(), 1e6) * rng.uniform(0.2, 5.0, size=(1, 300))
    p = _write(tmp_path, tpm)
    with pytest.raises(IngestError, match="TPM/FPKM"):
        load_any(p)


def test_candidate_layer_is_suggested_but_never_used(tmp_path):
    raw = _umi()
    p = _write(tmp_path, np.log1p(_normalized(raw)), umi=raw, spliced=np.round(raw * 0.8),
               velocity=np.log1p(raw))
    with pytest.raises(IngestError) as exc:
        load_any(p)  # refused: the integer layers are not used in place of X
    msg = str(exc.value)
    assert "layers['umi']" in msg and "layers['spliced']" in msg
    assert "layers['velocity']" not in msg  # not integer: not a candidate
    assert "None is used automatically" in msg and "store it as layers['counts']" in msg


def test_raw_that_is_not_counts_is_reported(tmp_path):
    logn = np.log1p(_normalized(_umi()))
    a = ad.AnnData(X=sp.csr_matrix(logn.astype(np.float32)),
                   var=pd.DataFrame(index=[f"G{i}" for i in range(logn.shape[1])]))
    a.raw = a.copy()
    p = tmp_path / "raw_is_log.h5ad"
    a.write_h5ad(p)
    with pytest.raises(IngestError, match=r"\.raw is present but is not counts either"):
        load_any(p)


def test_broken_counts_layer_suggests_x_when_x_holds_counts(tmp_path):
    raw = _umi()
    p = _write(tmp_path, raw, counts=np.log1p(_normalized(raw)))
    with pytest.raises(IngestError, match=r"layers\['counts'\] does not look like raw counts") as exc:
        load_any(p)
    msg = str(exc.value)
    assert "log-normalized" in msg and "matrices in this file: X" in msg
    assert "Fix or drop layers['counts']" in msg


def test_counts_layer_validator_still_works_without_the_anndata():
    layer = np.array([[1.5, 2.0], [3.0, 4.7]])
    with pytest.raises(IngestError, match="raw counts"):
        _validate_counts_layer(np.ones_like(layer), layer)


# ---------------------------------------------------------- calibration script #

def _mtx_folder(folder: Path, X):
    folder.mkdir(parents=True)
    scipy.io.mmwrite(str(folder / "matrix.mtx"), sp.coo_matrix(np.asarray(X).T.astype(np.int64)))
    (folder / "genes.tsv").write_text(
        "".join(f"ENSG{i:011d}\tG{i}\n" for i in range(X.shape[1])), encoding="utf-8")
    (folder / "barcodes.tsv").write_text(
        "".join(f"CELL{i}-1\n" for i in range(X.shape[0])), encoding="utf-8")


def test_calibration_script_reports_signals_relative_paths_and_summary(tmp_path, capsys):
    root = tmp_path / "data"
    raw = _umi()
    _write(root, raw, name="cellxgene/umi/a/data.h5ad")
    _write(root, _full_length(), name="cellxgene/smart-seq2/b/data.h5ad")
    _write(root, np.round(_normalized(raw)), name="rounded/data.h5ad")
    _write(root, np.log1p(_normalized(raw)), name="beeline_like/data.h5ad")
    _mtx_folder(root / "10x" / "filtered" / "hg19", _umi(seed=4))
    out_tsv = tmp_path / "plausibility.tsv"

    assert cp_script.main([str(root), "--tsv", str(out_tsv)]) == 0
    out = capsys.readouterr().out
    assert str(tmp_path) not in out  # no personal path in what gets pasted into a record
    assert out.strip().splitlines()[-1] == (
        "files: 5, would warn: 1, ok: 3, refused: 1, errors: 0")
    rows = {line.split("\t")[0]: line.split("\t") for line in out.splitlines()
            if line and not line.startswith(("#", "path", "files:"))}
    assert set(rows) == {"cellxgene/umi/a/data.h5ad", "cellxgene/smart-seq2/b/data.h5ad",
                         "rounded/data.h5ad", "beeline_like/data.h5ad", "10x/filtered/hg19"}
    assert rows["rounded/data.h5ad"][1] == "warn"
    assert "near_identical_totals" in rows["rounded/data.h5ad"][9]
    assert rows["cellxgene/smart-seq2/b/data.h5ad"][1] == "ok"
    assert rows["beeline_like/data.h5ad"][1] == "refused"
    assert rows["10x/filtered/hg19"][1] == "ok"
    assert out_tsv.read_text(encoding="utf-8").count("\n") == 6  # header + 5 rows


def test_calibration_script_fails_when_nothing_was_found(tmp_path, capsys):
    (tmp_path / "empty").mkdir()
    assert cp_script.main([str(tmp_path / "empty")]) == 1
    assert "files: 0, would warn: 0" in capsys.readouterr().out
