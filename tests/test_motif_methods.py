"""Motif-prior and cisTarget-pruning methods, plus the RegDiffusion adapter.

A tiny cisTarget resource folder is written per test: a rankings .feather (motifs x
genes, genome-wide ranks) and a motif -> TF annotation .tbl, so no download is needed.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from grnbench.infer import list_methods, run_method
from grnbench.infer import cistarget as ct

pa = pytest.importorskip("pyarrow")
import pyarrow.feather as pf  # noqa: E402

N_BACKGROUND = 400        # genome-wide genes in the fake database
TARGETS = [f"T{i}" for i in range(25)]      # TFA's true targets: motif near them
OTHERS = [f"G{i}" for i in range(25)]       # co-expressed with TFA, no motif


def _write_resources(root, rng):
    genes = ["TFA", "TFB"] + TARGETS + OTHERS + [f"BG{i}" for i in range(N_BACKGROUND)]
    motifs = ["mA"] + [f"m{i}" for i in range(60)]
    ranks = np.array([rng.permutation(len(genes)) for _ in motifs], dtype=np.int16)
    # motif mA ranks TFA's targets at the very top of the genome, OTHERS at the bottom
    order = TARGETS + [g for g in genes if g not in TARGETS and g not in OTHERS] + OTHERS
    ranks[0] = np.array([order.index(g) for g in genes], dtype=np.int16)
    table = pd.DataFrame(ranks, columns=genes)
    table["motifs"] = motifs
    pf.write_feather(table, str(root / "hg38_test.genes_vs_motifs.rankings.feather"))
    ann = pd.DataFrame({"#motif_id": ["mA", "m1"], "gene_name": ["TFA", "TFB"]})
    ann.to_csv(root / "motifs-test.hgnc-m0.001-o0.0.tbl", sep="\t", index=False)


@pytest.fixture()
def resources(tmp_path, monkeypatch):
    _write_resources(tmp_path, np.random.default_rng(0))
    monkeypatch.setenv(ct.ENV_DIR, str(tmp_path))
    ct.load_annotation.cache_clear()
    return tmp_path


@pytest.fixture()
def coexpressed():
    """TFA drives TARGETS and OTHERS equally; only TARGETS carry TFA's motif."""
    rng = np.random.default_rng(1)
    n = 300
    tfa = rng.normal(size=n)
    cols = {"TFA": tfa, "TFB": rng.normal(size=n)}
    for g in TARGETS + OTHERS:
        cols[g] = tfa + rng.normal(scale=0.7, size=n)
    for i in range(30):
        cols[f"BG{i}"] = rng.normal(size=n)
    return pd.DataFrame(cols)


def test_unavailable_without_resources(monkeypatch):
    monkeypatch.delenv(ct.ENV_DIR, raising=False)
    table = {m["name"]: m for m in list_methods()}
    for name in ("motif_prior", "pearson_motif", "pearson_cistarget", "grnboost2_cistarget"):
        assert not table[name]["available"]
        assert ct.ENV_DIR in table[name]["unavailable_reason"] or "arboreto" in table[name]["unavailable_reason"]


def test_guess_species_from_symbol_case():
    assert ct.guess_species(["TP53", "GATA1", "MT-CO1"]) == "human"
    assert ct.guess_species(["Trp53", "Gata1", "mt-Co1"]) == "mouse"


def test_motif_ranks_take_the_best_motif_and_worst_rank_when_missing(resources):
    mr = ct.motif_ranks(["TFA", "TFB", "NOMOTIF"], ["T0", "G0", "ABSENT"], species="human")
    assert mr.loc["TFA", "T0"] < 25                       # planted: targets at the top
    assert (mr.loc["NOMOTIF"] == mr.to_numpy().max()).all()  # no motif -> worst rank
    assert mr.loc["TFA", "ABSENT"] == mr.to_numpy().max()   # gene not in the database


def test_cistarget_keeps_motif_supported_targets_first(resources, coexpressed):
    edges = pd.DataFrame([("TFA", g, 1.0) for g in OTHERS + TARGETS], columns=["source", "target", "weight"])
    edges["weight"] = np.linspace(1.0, 0.5, len(edges))   # OTHERS ranked above TARGETS by the base
    out = ct.cistarget_rerank(edges, coexpressed, species="human")
    assert out.attrs["cistarget_regulons"] == 1
    top = set(out["target"].iloc[:out.attrs["cistarget_supported"]])
    # the leading edge stops where the motif's recovery curve beats random curves the most:
    # (almost) all planted targets, at most a few co-expressed genes that fall there by chance
    assert len(top & set(TARGETS)) >= 20 and len(top & set(OTHERS)) <= 3
    assert len(out) == len(edges)                         # nothing dropped: same universe


def test_motif_methods_run_through_the_registry(resources, coexpressed):
    import anndata as ad
    adata = ad.AnnData(coexpressed.to_numpy(dtype=np.float32) + 5.0,
                       var=pd.DataFrame(index=coexpressed.columns))
    tfs = ["TFA", "TFB"]
    prior, _ = run_method("motif_prior", adata, tf_list=tfs, assume_normalized=True)
    assert set(prior["source"]) == set(tfs)
    assert set(prior[prior.source == "TFA"].head(25)["target"]) == set(TARGETS)
    combined, _ = run_method("pearson_motif", adata, tf_list=tfs, assume_normalized=True)
    tfa = combined[combined.source == "TFA"].head(25)
    # rank product lifts motif-backed co-expressed targets over equally co-expressed ones
    assert len(set(tfa["target"]) & set(TARGETS)) > len(set(tfa["target"]) & set(OTHERS))
    assert combined["weight"].is_monotonic_decreasing


def test_regdiffusion_adapter_returns_canonical_edges():
    pytest.importorskip("regdiffusion")
    import anndata as ad
    rng = np.random.default_rng(0)
    X = rng.poisson(3.0, size=(120, 12)).astype(np.float32)
    adata = ad.AnnData(X, var=pd.DataFrame(index=[f"g{i}" for i in range(12)]))
    edges, prov = run_method("regdiffusion", adata, tf_list=["g0", "g1"], n_steps=20)
    assert list(edges.columns[:3]) == ["source", "target", "weight"]
    assert set(edges["source"]) <= {"g0", "g1"}
    assert (edges["source"] != edges["target"]).all()
    assert prov["normalization"].startswith("raw counts detected")


def _real_bench():
    import importlib.util
    import sys
    from pathlib import Path
    bench = Path(__file__).parents[1] / "benchmarks"
    sys.path.insert(0, str(bench))          # real_bench imports engine_bench as a sibling
    spec = importlib.util.spec_from_file_location("real_bench_models", bench / "real_bench.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_real_bench_combined_names_resolve_to_registered_bases():
    rb = _real_bench()
    assert rb._split_combined("pearson_motif") == ("pearson_signed", "_motif")
    assert rb._split_combined("grnboost2_cistarget") == ("grnboost2", "_cistarget")
    assert rb._split_combined("motif_prior") == ("motif_prior", None)
    registered = {m["name"] for m in list_methods()}
    for a, b in rb.DECLARED_PAIRS:
        assert a in registered and b in registered


def test_real_bench_models_section_pairs_over_datasets_only():
    rb = _real_bench()
    rng = np.random.default_rng(0)
    rows = []
    for i in range(6):                       # 6 datasets = 6 units
        row = {"task": f"cxg__d{i}"}
        for arm in ("auto", "baseline"):
            row[f"m_pearson_signed_{arm}_auprc_ratio"] = 1.2 + rng.normal(0, .01)
            row[f"m_pearson_motif_{arm}_auprc_ratio"] = 1.6 + rng.normal(0, .01)
            row[f"m_pearson_signed_{arm}_epr"] = 2.0
            row[f"m_pearson_motif_{arm}_epr"] = np.nan if i == 0 else 2.5 + rng.normal(0, .01)
            row[f"m_pearson_signed_{arm}_seconds"] = row[f"m_pearson_motif_{arm}_seconds"] = 1.0
        rows.append(row)
    text = "\n".join(rb.models_section(pd.DataFrame(rows)))
    assert "| pearson_motif | pearson_signed | 6 |" in text        # AUPRC: all 6 units
    assert "| pearson_motif | pearson_signed | 5 |" in text        # EPR: the NaN unit dropped
    assert "supported" in text and "| pearson_motif |" in text


def test_regdiffusion_drops_zero_variance_genes(monkeypatch):
    """Issue #13: RegDiffusion refuses constant genes ('1 genes have 0 variance'), which appear
    after the 3,000-cell subsample of real data. The adapter trains on the varying genes only;
    constant genes get no edges. Stubbed regdiffusion: runs without torch."""
    import sys
    import types

    from grnbench.infer.regdiffusion_adapter import regdiffusion_runner

    seen = {}

    class Trainer:
        def __init__(self, X, **kw):
            if (np.asarray(X).std(axis=0) == 0).any():
                raise ValueError("1 genes have 0 variance. Please remove these genes from your data.")
            seen["n_genes"] = X.shape[1]
            self.n = X.shape[1]

        def train(self):
            pass

        def get_adj(self):
            return np.arange(self.n * self.n, dtype=float).reshape(self.n, self.n) + 1.0

    monkeypatch.setitem(sys.modules, "regdiffusion", types.SimpleNamespace(RegDiffusionTrainer=Trainer))
    rng = np.random.default_rng(0)
    expr = pd.DataFrame(rng.normal(size=(50, 5)), columns=["TF1", "TF2", "g1", "g2", "g3"])
    expr["g2"] = 1.5                       # constant after subsampling
    expr["TF2"] = 0.0                      # a constant regulator too
    edges = regdiffusion_runner(expr, tf_names=["TF1", "TF2"], n_steps=1)
    assert seen["n_genes"] == 3
    assert set(edges["source"]) == {"TF1"}
    assert "g2" not in set(edges["target"]) and "TF2" not in set(edges["target"])
    assert (edges["weight"] > 0).all()
