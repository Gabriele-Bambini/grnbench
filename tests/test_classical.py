"""The numpy/sklearn method family.

Each method is tested against the thing it exists to do, not merely that it
returns a DataFrame. A method that runs and produces well-formed garbage passes
a smoke test and fails the user.
"""
from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from grnbench.eval.metrics import evaluate_ranking
from grnbench.infer import list_methods, run_method

NUMPY_ONLY = ["pearson_signed", "mutual_info", "clr"]
SKLEARN = ["ppcor", "lasso", "random_forest"]
ALL_NEW = NUMPY_ONLY + SKLEARN


def _adata(X, genes):
    return ad.AnnData(np.asarray(X, dtype="float32"),
                      obs=pd.DataFrame(index=[f"c{i}" for i in range(len(X))]),
                      var=pd.DataFrame(index=list(genes)))


@pytest.fixture(scope="module")
def planted():
    """25 TFs driving 75 targets through a saturating response, plus noise."""
    rng = np.random.default_rng(0)
    n_cells, n_tf, n_other = 200, 25, 75
    tfs = [f"TF{i:02d}" for i in range(n_tf)]
    others = [f"G{i:03d}" for i in range(n_other)]
    act = rng.lognormal(0.0, 0.7, size=(n_cells, n_tf))
    tgt = np.zeros((n_cells, n_other))
    gold = set()
    for t in range(n_tf):
        for p in rng.choice(n_other, size=3, replace=False):
            tgt[:, p] += 2.0 * np.tanh(act[:, t] - 1.0)
            gold.add((tfs[t], others[p]))
    tgt += rng.normal(0, 0.3, size=tgt.shape)
    X = np.hstack([act, np.maximum(tgt, 0.0)])
    return _adata(X, tfs + others), tfs, gold


# --------------------------------------------------------------------------- #
# contract
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("method", ALL_NEW)
def test_method_is_registered_and_described(method):
    entry = next(m for m in list_methods() if m["name"] == method)
    assert entry["description"], "a registered method must say what it is"
    assert entry["origin"] == "builtin"


@pytest.mark.parametrize("method", ALL_NEW)
def test_edge_list_obeys_the_canonical_contract(method, planted):
    adata, tfs, _ = planted
    edges, prov = run_method(method, adata, tf_list=tfs, assume_normalized=True)

    assert list(edges.columns[:3]) == ["source", "target", "weight"]
    assert not edges.empty
    # ranked descending: every metric in the package assumes it
    w = edges["weight"].to_numpy()
    assert np.all(np.diff(w) <= 1e-12), "edges must arrive sorted by weight"
    assert np.isfinite(w).all(), "no NaN/Inf weights"
    assert (edges["source"] != edges["target"]).all(), "self-loops must be dropped"
    # the regulator whitelist is a restriction, not a suggestion
    assert set(edges["source"]) <= set(tfs)
    assert prov, "run_method must return provenance"


@pytest.mark.parametrize("method", ALL_NEW)
def test_method_is_deterministic(method, planted):
    adata, tfs, _ = planted
    a, _ = run_method(method, adata, tf_list=tfs, assume_normalized=True)
    b, _ = run_method(method, adata, tf_list=tfs, assume_normalized=True)
    pd.testing.assert_frame_equal(a.reset_index(drop=True), b.reset_index(drop=True))


@pytest.mark.parametrize("method", ALL_NEW)
def test_method_beats_chance_on_a_planted_network(method, planted):
    adata, tfs, gold = planted
    edges, _ = run_method(method, adata, tf_list=tfs, assume_normalized=True)
    res = evaluate_ranking(edges, gold, with_curve=False)
    assert res["n_gold_in_list"] > 0, "the gold must be reachable at all"
    assert res["auprc_ratio"] > 1.5, (
        f"{method} scored {res['auprc_ratio']:.2f}x random on a planted network")


# --------------------------------------------------------------------------- #
# each method's reason to exist
# --------------------------------------------------------------------------- #

def test_mutual_info_sees_a_dependence_correlation_is_blind_to():
    """A U-shaped response has ~zero linear correlation but plenty of MI.

    This is the entire justification for carrying an MI method next to the
    correlation baselines; if it fails, `mutual_info` is dead weight.
    """
    rng = np.random.default_rng(1)
    n = 400
    driver = rng.uniform(-2, 2, n)
    response = driver ** 2 + rng.normal(0, 0.15, n)      # symmetric -> r ~ 0
    filler = rng.normal(0, 1, (n, 6))
    X = np.column_stack([driver, response, filler])
    genes = ["TF0", "TARGET"] + [f"N{i}" for i in range(6)]
    adata = _adata(X, genes)

    assert abs(np.corrcoef(driver, response)[0, 1]) < 0.15, "fixture must be non-linear"

    pear, _ = run_method("pearson_signed", adata, tf_list=["TF0"],
                         assume_normalized=True)
    mi, _ = run_method("mutual_info", adata, tf_list=["TF0"], assume_normalized=True)

    def rank_of_target(df):
        return int(np.flatnonzero(df["target"].to_numpy() == "TARGET")[0])

    assert rank_of_target(mi) == 0, "MI must rank the true target first"
    assert rank_of_target(pear) > 0, "correlation is expected to miss it"


def test_clr_demotes_a_hub_that_raw_mi_promotes():
    """CLR exists to normalise away genes that correlate with everything.

    A broadcast gene driving every other gene gets high MI with all of them.
    CLR z-scores each gene against its own MI background, so that hub's edges
    should fall relative to a specific pair's.
    """
    rng = np.random.default_rng(2)
    n, n_bg = 400, 12
    hub = rng.normal(0, 1, n)
    bg = 0.9 * hub[:, None] + rng.normal(0, 0.4, (n, n_bg))   # hub drives all
    a = rng.normal(0, 1, n)
    b = a + rng.normal(0, 0.25, n)                            # one specific pair
    X = np.column_stack([hub, a, b, bg])
    genes = ["HUB", "A", "B"] + [f"BG{i}" for i in range(n_bg)]
    adata = _adata(X, genes)
    tfs = ["HUB", "A"]

    def weight(df, s, t):
        m = (df["source"] == s) & (df["target"] == t)
        return float(df.loc[m, "weight"].iloc[0])

    mi, _ = run_method("mutual_info", adata, tf_list=tfs, assume_normalized=True)
    cl, _ = run_method("clr", adata, tf_list=tfs, assume_normalized=True)

    # ratio of the specific pair to a hub edge, before and after CLR
    mi_ratio = weight(mi, "A", "B") / max(weight(mi, "HUB", "BG0"), 1e-12)
    clr_ratio = weight(cl, "A", "B") / max(weight(cl, "HUB", "BG0"), 1e-12)
    assert clr_ratio > mi_ratio, (
        f"CLR should favour the specific pair over the hub edge "
        f"(MI ratio {mi_ratio:.2f} -> CLR ratio {clr_ratio:.2f})")


def test_ppcor_drops_an_indirect_edge_that_correlation_keeps():
    """A -> B -> C makes A and C correlate; conditioning should break it.

    That is the whole point of a partial-correlation method, and the reason
    BEELINE carries PPCOR alongside plain correlation.
    """
    rng = np.random.default_rng(3)
    n = 500
    A = rng.normal(0, 1, n)
    B = A + rng.normal(0, 0.2, n)
    C = B + rng.normal(0, 0.2, n)          # C depends on A only through B
    noise = rng.normal(0, 1, (n, 5))
    X = np.column_stack([A, B, C, noise])
    genes = ["A", "B", "C"] + [f"N{i}" for i in range(5)]
    adata = _adata(X, genes)

    def weight(df, s, t):
        m = (df["source"] == s) & (df["target"] == t)
        return float(df.loc[m, "weight"].iloc[0])

    pear, _ = run_method("pearson_signed", adata, tf_list=["A", "B"],
                         assume_normalized=True)
    pp, _ = run_method("ppcor", adata, tf_list=["A", "B"], assume_normalized=True)

    # correlation cannot tell the indirect A-C from the direct B-C
    assert weight(pear, "A", "C") > 0.8, "A and C must correlate strongly"
    # partial correlation must weaken A-C far more than it weakens B-C
    assert weight(pp, "A", "C") < weight(pear, "A", "C") / 2
    assert weight(pp, "A", "C") < weight(pp, "B", "C")


def test_random_forest_matches_genie3_where_both_are_available():
    """`random_forest` exists so GENIE3's estimator runs without arboreto.

    They are not bit-identical (different tree implementations), but they must
    agree on which edges matter; if they do not, one of them is wrong.
    """
    avail = {m["name"] for m in list_methods() if m["available"]}
    if "genie3" not in avail:
        pytest.skip("arboreto not installed")

    rng = np.random.default_rng(4)
    n, n_tf, n_tgt = 150, 6, 20
    act = rng.lognormal(0, 0.6, (n, n_tf))
    W = np.zeros((n_tf, n_tgt))
    for t in range(n_tf):
        W[t, rng.choice(n_tgt, 3, replace=False)] = rng.uniform(0.8, 1.5)
    tgt = act @ W + rng.normal(0, 0.2, (n, n_tgt))
    genes = [f"TF{i}" for i in range(n_tf)] + [f"G{i}" for i in range(n_tgt)]
    adata = _adata(np.hstack([act, tgt]), genes)
    tfs = genes[:n_tf]

    rf, _ = run_method("random_forest", adata, tf_list=tfs, assume_normalized=True)
    g3, _ = run_method("genie3", adata, tf_list=tfs, assume_normalized=True)

    k = 20
    top_rf = set(zip(rf.head(k)["source"], rf.head(k)["target"]))
    top_g3 = set(zip(g3.head(k)["source"], g3.head(k)["target"]))
    jaccard = len(top_rf & top_g3) / len(top_rf | top_g3)
    assert jaccard > 0.4, (
        f"sklearn and arboreto random forests disagree on the top {k} edges "
        f"(Jaccard {jaccard:.2f})")


def test_tf_list_none_makes_every_gene_a_regulator(planted):
    adata, _, _ = planted
    edges, _ = run_method("pearson_signed", adata, assume_normalized=True)
    assert set(edges["source"]) == set(adata.var_names)


# --------------------------------------------------------------------------- #
# symmetric statistics: A->B and B->A are bit-identical
# --------------------------------------------------------------------------- #

SYMMETRIC = ["pearson_signed", "spearman_signed", "mutual_info", "clr", "ppcor"]


@pytest.fixture(scope="module")
def dup_data():
    """300 cells x 60 genes (float32), 20 regulators, a few near-duplicate genes."""
    rng = np.random.default_rng(1)
    X = rng.lognormal(0, 1, (300, 60)).astype("float32")
    X[:, 1] = X[:, 0] + rng.normal(0, 0.01, 300).astype("float32")
    X[:, 21] = X[:, 20] * 2 + rng.normal(0, 0.01, 300).astype("float32")
    X[:, 41] = X[:, 3]
    X[rng.random(X.shape) < 0.3] = 0
    genes = [f"G{i}" for i in range(60)]
    return _adata(X, genes), genes[:20]


@pytest.mark.parametrize("method", SYMMETRIC)
def test_symmetric_statistic_gives_identical_weights_both_ways(method, dup_data):
    adata, regs = dup_data
    edges, _ = run_method(method, adata, regs)
    w = edges.set_index(["source", "target"])["weight"]
    pairs = [(a, b) for a in regs for b in regs if a < b and (a, b) in w.index
             and (b, a) in w.index]
    assert len(pairs) > 100
    ab = np.array([w[p] for p in pairs])
    ba = np.array([w[p[::-1]] for p in pairs])
    assert np.array_equal(ab, ba)


def test_pearson_weights_still_equal_abs_corrcoef(dup_data):
    adata, regs = dup_data
    edges, _ = run_method("pearson_signed", adata, regs, assume_normalized=True)
    genes = list(adata.var_names)
    ref = np.abs(np.corrcoef(adata.X.astype(np.float64), rowvar=False))
    ref = np.nan_to_num(ref, nan=0.0)
    got = edges.assign(i=edges["source"].map(genes.index),
                       j=edges["target"].map(genes.index))
    assert np.allclose(got["weight"], ref[got["i"], got["j"]], rtol=0, atol=1e-12)
