"""Preprocessing gaps found by running the CLI on a deliberately messy dataset.

Each test here reproduces something the pipeline used to get wrong or leave
undone: auditing only one of two label columns, folding missing labels into a
literal "nan" stratum, carrying duplicated barcodes through untouched, flagging
doublets without any way to act on them, and detecting a batch column it then
never used.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import anndata as ad
import scipy.sparse as sp

from grnbench.pipeline import PrepError, run_prep


def _messy(n_cells=600, n_genes=500, seed=0, dup_barcodes=0):
    rng = np.random.default_rng(seed)
    X = rng.poisson(3.0, size=(n_cells, n_genes)).astype(np.float32)
    obs = pd.DataFrame(index=[f"cell{i}" for i in range(n_cells)])
    obs["donor"] = np.repeat([f"D{i}" for i in range(4)], n_cells // 4)
    obs["cell_type"] = ["A"] * (n_cells - 100) + ["B"] * 90 + ["  Rare  "] * 6 + ["?"] * 4
    if dup_barcodes:
        names = list(obs.index)
        for i in range(dup_barcodes):
            names[i + 1] = names[0]
        obs.index = pd.Index(names)
    return ad.AnnData(X=sp.csr_matrix(X), obs=obs,
                      var=pd.DataFrame(index=[f"G{i}" for i in range(n_genes)]))


# ------------------------------------------------ imbalance audits both keys #

def test_group_key_and_label_key_are_both_audited(tmp_path):
    """With --group-key donor (balanced) AND --label-key cell_type (500/90/6),
    the class-imbalance gate used to look only at donor and report "balanced"."""
    src = tmp_path / "in.h5ad"
    _messy().write_h5ad(src)
    res = run_prep(src, tmp_path / "run", group_key="donor", label_key="cell_type")
    gates = {g["name"]: g for g in res.manifest["qc"]["gates"]}
    assert gates["class_imbalance"]["status"] == "warn"
    assert "cell_type" in gates["class_imbalance"]["detail"]
    assert gates["group_imbalance"]["status"] == "pass"
    assert "donor" in gates["group_imbalance"]["detail"]


def test_group_key_alone_keeps_the_single_gate(tmp_path):
    src = tmp_path / "in.h5ad"
    _messy().write_h5ad(src)
    res = run_prep(src, tmp_path / "run", group_key="cell_type")
    names = [g["name"] for g in res.manifest["qc"]["gates"]]
    assert "class_imbalance" in names and "group_imbalance" not in names


# ----------------------------------------- missing labels and stratification #

def test_unlabelled_cells_are_excluded_from_stratified_subsample(tmp_path):
    """Missing labels ("?") became a literal "nan" stratum that competed for
    quota like a real class. They are now excluded and counted."""
    src = tmp_path / "in.h5ad"
    _messy().write_h5ad(src)
    res = run_prep(src, tmp_path / "run", label_key="cell_type",
                   stratify_key="cell_type", subsample_cells=120)
    step = next(s for s in res.manifest["steps"] if s["step"] == "subsample_cells")
    assert step["n_unlabelled_excluded"] == 4
    assert "nan" not in step["stratification"]["classes_dropped_examples"]
    assert step["stratification"]["n_classes"] == 3
    out = ad.read_h5ad(res.processed_path)
    assert out.n_obs == 120 and out.obs["cell_type"].notna().all()
    assert any("4 cells" in w and "cell_type" in w for w in res.manifest["warnings"])


def test_unlabelled_cells_are_kept_but_reported_without_subsampling(tmp_path):
    src = tmp_path / "in.h5ad"
    _messy().write_h5ad(src)
    res = run_prep(src, tmp_path / "run", label_key="cell_type")
    out = ad.read_h5ad(res.processed_path)
    assert int(out.obs["cell_type"].isna().sum()) == 4
    assert any("4 cells" in w and "no usable" in w for w in res.manifest["warnings"])


# ------------------------------------------------------- duplicate barcodes #

def test_duplicate_barcodes_are_made_unique_and_preserved(tmp_path):
    src = tmp_path / "in.h5ad"
    _messy(dup_barcodes=3).write_h5ad(src)
    res = run_prep(src, tmp_path / "run")
    step = next(s for s in res.manifest["steps"] if s["step"] == "dedup_barcodes")
    assert step["n_duplicates"] == 3
    out = ad.read_h5ad(res.processed_path)
    assert out.obs_names.is_unique
    assert "barcode_as_loaded" in out.obs.columns
    assert (out.obs["barcode_as_loaded"] == "cell0").sum() == 4
    gate = next(g for g in res.manifest["qc"]["gates"] if g["name"] == "duplicate_barcodes")
    assert gate["status"] == "warn" and "made unique" in gate["detail"]


# ----------------------------------------------------------------- doublets #

def _two_populations_with_doublets(n_per=250, n_genes=400, n_doublets=30, seed=3):
    """Two expression programs plus synthetic doublets made by summing one cell
    of each - the textbook scrublet-detectable case."""
    rng = np.random.default_rng(seed)
    base = rng.gamma(2.0, 1.0, size=n_genes)
    prog_a = base.copy(); prog_a[: n_genes // 2] *= 6
    prog_b = base.copy(); prog_b[n_genes // 2:] *= 6
    A = rng.poisson(prog_a, size=(n_per, n_genes))
    B = rng.poisson(prog_b, size=(n_per, n_genes))
    D = A[:n_doublets] + B[:n_doublets]
    X = np.vstack([A, B, D]).astype(np.float32)
    obs = pd.DataFrame({"truth": ["A"] * n_per + ["B"] * n_per + ["doublet"] * n_doublets},
                       index=[f"c{i}" for i in range(len(X))])
    return ad.AnnData(X=sp.csr_matrix(X), obs=obs,
                      var=pd.DataFrame(index=[f"G{i}" for i in range(n_genes)]))


def test_drop_doublets_uses_scrublet_and_logs_what_it_removed(tmp_path):
    src = tmp_path / "in.h5ad"
    _two_populations_with_doublets().write_h5ad(src)
    res = run_prep(src, tmp_path / "run", drop_doublets=True, seed=0)
    step = next(s for s in res.manifest["steps"] if s["step"] == "drop_doublets")
    assert step["method"] == "scrublet"
    out = ad.read_h5ad(res.processed_path)
    assert {"doublet_score", "predicted_doublet"} <= set(out.obs.columns)
    assert step["n_removed"] == 530 - out.n_obs
    # what went should be the synthetic doublets, not the singlets
    removed_doublets = 30 - int(out.obs["truth"].eq("doublet").sum())
    removed_singlets = 500 - int(out.obs["truth"].ne("doublet").sum())
    # scrublet's automatic threshold is conservative on a two-cluster toy (it
    # caught 15/30 here at threshold 0.27); a third is the floor that says the
    # calls are pointed at the doublets rather than at random cells
    assert removed_doublets >= 10, f"caught only {removed_doublets}/30 doublets"
    assert removed_singlets <= 25, f"dropped {removed_singlets} singlets"


def test_doublet_fallback_is_loud_and_recorded(tmp_path, monkeypatch):
    """If scrublet cannot run (e.g. scikit-image missing), the heuristic runs, a
    RuntimeWarning is raised, and the manifest records the fallback and its reason."""
    import scanpy as sc

    def broken(*args, **kwargs):
        raise ValueError("threshold is None and thus scrublet requires skimage, "
                         "but skimage is not installed.")

    monkeypatch.setattr(sc.pp, "scrublet", broken)
    src = tmp_path / "in.h5ad"
    _two_populations_with_doublets().write_h5ad(src)
    with pytest.warns(RuntimeWarning, match="DOUBLET FALLBACK"):
        res = run_prep(src, tmp_path / "run", drop_doublets=True, seed=0)
    step = next(s for s in res.manifest["steps"] if s["step"] == "drop_doublets")
    assert step["method"] == "heuristic" and step["fallback"] is True
    assert "scikit-image" in step["fallback_reason"]
    assert any("DOUBLET FALLBACK" in w for w in res.manifest["warnings"])


def test_scrublet_run_is_not_marked_as_fallback(tmp_path):
    pytest.importorskip("skimage")
    src = tmp_path / "in.h5ad"
    _two_populations_with_doublets().write_h5ad(src)
    res = run_prep(src, tmp_path / "run", drop_doublets=True, seed=0)
    step = next(s for s in res.manifest["steps"] if s["step"] == "drop_doublets")
    assert step["method"] == "scrublet" and step["fallback"] is False and step["fallback_reason"] is None


def test_doublets_are_never_dropped_by_default(tmp_path):
    src = tmp_path / "in.h5ad"
    _two_populations_with_doublets().write_h5ad(src)
    res = run_prep(src, tmp_path / "run")
    assert all(s["step"] != "drop_doublets" for s in res.manifest["steps"])
    assert ad.read_h5ad(res.processed_path).n_obs == 530


def test_drop_doublets_falls_back_to_heuristic_loudly(tmp_path, monkeypatch):
    import scanpy as sc

    def boom(*a, **k):
        raise RuntimeError("scrublet exploded")
    monkeypatch.setattr(sc.pp, "scrublet", boom)
    src = tmp_path / "in.h5ad"
    _two_populations_with_doublets().write_h5ad(src)
    res = run_prep(src, tmp_path / "run", drop_doublets=True)
    step = next(s for s in res.manifest["steps"] if s["step"] == "drop_doublets")
    assert step["method"] == "heuristic"
    assert any("scrublet" in w and "heuristic" in w for w in res.manifest["warnings"])


# ------------------------------------------------------------------- ComBat #

def _batched(n_per=200, n_genes=300, seed=5):
    rng = np.random.default_rng(seed)
    rate = rng.gamma(2.0, 2.0, size=n_genes)
    b1 = rng.poisson(rate, size=(n_per, n_genes))
    b2 = rng.poisson(rate * 3.0, size=(n_per, n_genes))  # batch 2: 3x depth, gene-wise shift
    shift = rng.choice([0.3, 1.0, 3.0], size=n_genes)
    b2 = rng.poisson(rate * shift, size=(n_per, n_genes))
    X = np.vstack([b1, b2]).astype(np.float32)
    obs = pd.DataFrame({"batch": ["b1"] * n_per + ["b2"] * n_per},
                       index=[f"c{i}" for i in range(2 * n_per)])
    return ad.AnnData(X=sp.csr_matrix(X), obs=obs,
                      var=pd.DataFrame(index=[f"G{i}" for i in range(n_genes)]))


def _batch_gap(adata):
    X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
    b = adata.obs["batch"].to_numpy()
    return float(np.abs(X[b == "b1"].mean(0) - X[b == "b2"].mean(0)).mean())


def test_combat_strategy_removes_the_batch_offset(tmp_path):
    src = tmp_path / "in.h5ad"
    _batched().write_h5ad(src)
    plain = run_prep(src, tmp_path / "plain", strategy="cpm_log1p", batch_key="batch")
    fixed = run_prep(src, tmp_path / "combat", strategy="cpm_log1p_combat", batch_key="batch")
    gap_plain = _batch_gap(ad.read_h5ad(plain.processed_path))
    gap_fixed = _batch_gap(ad.read_h5ad(fixed.processed_path))
    assert gap_fixed < 0.2 * gap_plain
    step = next(s for s in fixed.manifest["steps"] if s["step"] == "strategy")
    assert step["strategy"] == "cpm_log1p_combat" and step["params"]["batch_key"] == "batch"
    out = ad.read_h5ad(fixed.processed_path)
    assert "counts" in out.layers  # raw counts survive the correction


def test_combat_strategy_refuses_without_a_batch_column(tmp_path):
    src = tmp_path / "in.h5ad"
    a = _batched(); a.obs = a.obs.drop(columns=["batch"]); a.write_h5ad(src)
    with pytest.raises(PrepError, match="batch"):
        run_prep(src, tmp_path / "run", strategy="cpm_log1p_combat")


def test_combat_is_never_auto_recommended(tmp_path):
    src = tmp_path / "in.h5ad"
    _batched().write_h5ad(src)
    res = run_prep(src, tmp_path / "run", batch_key="batch")
    assert res.strategy != "cpm_log1p_combat"
    assert all(r["strategy"] != "cpm_log1p_combat" for r in res.manifest["recommendations"])


def test_drop_doublets_refuses_an_implausible_call_rate(tmp_path, monkeypatch):
    """On structureless data scrublet called 83% of cells doublets and the run
    went on as if that were a result. Above the ceiling it is a failed
    threshold and the pipeline must say so instead of deleting the dataset."""
    import scanpy as sc

    def half_are_doublets(adata, **k):
        n = adata.n_obs
        adata.obs["doublet_score"] = np.linspace(0, 1, n)
        adata.obs["predicted_doublet"] = np.arange(n) % 2 == 0
    monkeypatch.setattr(sc.pp, "scrublet", half_are_doublets)
    src = tmp_path / "in.h5ad"
    _messy().write_h5ad(src)
    with pytest.raises(PrepError, match="ceiling"):
        run_prep(src, tmp_path / "run", drop_doublets=True)
    # a deliberately raised ceiling lets it through, and the manifest shows it
    from grnbench.qc import Thresholds
    res = run_prep(src, tmp_path / "run2", drop_doublets=True,
                   thresholds=Thresholds(doublet_rate_max=0.6))
    step = next(s for s in res.manifest["steps"] if s["step"] == "drop_doublets")
    assert step["n_removed"] == 300  # every other of 600 cells
    assert res.manifest["thresholds"]["doublet_rate_max"] == 0.6


def test_non_finite_strategy_output_is_refused(tmp_path, monkeypatch):
    import grnbench.pipeline as pl
    real = pl.apply_strategy

    def poisoned(name, adata, **kw):
        out, prov = real(name, adata, **kw)
        out.X = out.X.toarray() if sp.issparse(out.X) else np.asarray(out.X)
        out.X[0, 0] = np.nan
        return out, prov

    monkeypatch.setattr(pl, "apply_strategy", poisoned)
    src = tmp_path / "in.h5ad"
    _batched().write_h5ad(src)
    with pytest.raises(PrepError, match="non-finite"):
        run_prep(src, tmp_path / "run", strategy="cpm_log1p")
    assert not (tmp_path / "run" / "processed.h5ad").exists()
