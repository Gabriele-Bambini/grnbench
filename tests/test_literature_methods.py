"""Methods added from the GRN literature (docs/models_catalog.md).

Every method is tested against what it exists to do, not only that it returns a frame:
the contract and determinism are parametrised over all of them, and each family has
tests that a wrong implementation would fail (a DPI that prunes nothing, a "robust"
correlation that is not, a community that is not the mean rank, a PCR that is not one).
"""
from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from grnbench.cli_infer import default_method_names
from grnbench.eval.metrics import GoldIndex, evaluate_ranking
from grnbench.infer import MethodError, MethodSpec, list_methods, register_method, run_method
from grnbench.infer import base as infer_base
from grnbench.infer import association as assoc
from grnbench.infer import regression as regr

NEW = ["bicor", "wgcna_tom", "propr_rho", "aracne_dpi", "mrnet", "tigress", "bbsr",
       "pcnet", "extra_trees", "gradient_boosting", "community"]
# the 17 methods registered before these: the set `--methods available` must keep running
OLD = {"spearman_signed", "genie3", "grnboost2", "aracne3", "pearson_signed", "mutual_info", "clr",
       "ppcor", "random_forest", "lasso", "regdiffusion", "motif_prior", "pearson_motif",
       "grnboost2_motif", "regdiffusion_motif", "pearson_cistarget", "grnboost2_cistarget"}
# small settings so the suite stays fast; the defaults are exercised by the smoke run in the PR
PARAMS = {"extra_trees": {"n_estimators": 20}, "gradient_boosting": {"n_estimators": 25},
          "tigress": {"n_resamples": 20}}


def _adata(X, genes):
    return ad.AnnData(np.asarray(X, dtype="float32"),
                      obs=pd.DataFrame(index=[f"c{i}" for i in range(len(X))]),
                      var=pd.DataFrame(index=list(genes)))


def _frame(X, genes):
    return pd.DataFrame(np.asarray(X, dtype=float), columns=list(genes))


def _weight(edges, source, target):
    row = edges[(edges["source"] == source) & (edges["target"] == target)]
    assert len(row) == 1, f"{source}->{target} not reported exactly once"
    return float(row["weight"].iloc[0])


@pytest.fixture(scope="module")
def planted():
    """12 TFs driving 3 targets each through a saturating response, 28 other genes, noise."""
    rng = np.random.default_rng(0)
    n_cells, n_tf, n_other = 160, 12, 28
    tfs = [f"TF{i:02d}" for i in range(n_tf)]
    others = [f"G{i:03d}" for i in range(n_other)]
    act = rng.lognormal(0.0, 0.7, size=(n_cells, n_tf))
    tgt = np.zeros((n_cells, n_other))
    gold = set()
    for t in range(n_tf):
        for p in rng.choice(n_other, size=3, replace=False):
            tgt[:, p] += 2.0 * np.tanh(act[:, t] - 1.0)
            gold.add((tfs[t], others[p]))
    tgt += rng.normal(0, 0.5, size=tgt.shape)
    return _adata(np.hstack([act, np.maximum(tgt, 0.0)]), tfs + others), tfs, gold


_first_run: dict[str, pd.DataFrame] = {}


def _run(method, planted):
    adata, tfs, _ = planted
    return run_method(method, adata, tf_list=tfs, assume_normalized=True, **PARAMS.get(method, {}))[0]


@pytest.fixture(scope="module")
def runs(planted):
    return lambda m: _first_run.setdefault(m, _run(m, planted))


# --------------------------------------------------------------------------- #
# registry: the default set is untouched
# --------------------------------------------------------------------------- #

def test_new_methods_are_registered_extra_and_described():
    table = {m["name"]: m for m in list_methods()}
    for name in NEW:
        entry = table[name]
        assert entry["description"], "a registered method must say what it is"
        assert entry["origin"] == "builtin"
        assert entry["default"] is False
        assert entry["available"], entry["unavailable_reason"]


def test_default_set_is_exactly_the_methods_registered_before():
    """`--methods available` keeps meaning what it meant for every registered experiment."""
    table = [m for m in list_methods() if m["origin"] == "builtin"]
    assert {m["name"] for m in table if m["default"]} == OLD
    chosen = set(default_method_names(table))
    assert chosen <= OLD and not (chosen & set(NEW))
    everything = set(default_method_names(table, everything=True))
    assert set(NEW) <= everything
    assert everything >= chosen


# --------------------------------------------------------------------------- #
# contract
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("method", NEW)
def test_edge_list_obeys_the_canonical_contract(method, runs, planted):
    _, tfs, _ = planted
    edges = runs(method)
    assert list(edges.columns[:3]) == ["source", "target", "weight"]
    assert not edges.empty
    w = edges["weight"].to_numpy()
    assert np.all(np.diff(w) <= 1e-12), "edges must arrive sorted by weight"
    assert np.isfinite(w).all(), "no NaN/Inf weights"
    assert (w >= 0).all(), "weights are evidence strengths, never negative"
    assert (edges["source"] != edges["target"]).all(), "self-loops must be dropped"
    assert set(edges["source"]) <= set(tfs), "the regulator whitelist is a restriction"
    assert not edges.duplicated(["source", "target"]).any()


@pytest.mark.parametrize("method", NEW)
def test_method_is_deterministic(method, runs, planted):
    first = runs(method)
    again = _run(method, planted)
    pd.testing.assert_frame_equal(first.reset_index(drop=True), again.reset_index(drop=True))


@pytest.mark.parametrize("method", [m for m in NEW if m not in ("pcnet",)])
def test_method_recovers_a_planted_network_far_above_chance(method, runs, planted):
    """pcnet is excluded on purpose: three principal components cannot isolate twelve
    independent modules (see its docstring); its correctness is tested against an
    independent PCR below."""
    adata, tfs, gold = planted
    genes = list(adata.var_names)
    universe = pd.DataFrame([(s, g) for s in tfs for g in genes if s != g], columns=["source", "target"])
    result = evaluate_ranking(runs(method), GoldIndex.coerce(gold), with_curve=False, universe=universe)
    chance = len(gold) / len(universe)
    assert result["auprc_ratio"] > 5.0, f"{method}: AUPRC {result['auprc'] if 'auprc' in result else ''} vs chance {chance:.3f}"


# --------------------------------------------------------------------------- #
# bicor
# --------------------------------------------------------------------------- #

def test_bicor_resists_outlier_cells_that_break_pearson():
    rng = np.random.default_rng(1)
    n = 200
    x = rng.normal(size=n)
    y = x + rng.normal(scale=0.3, size=n)
    x[:10], y[:10] = 25.0, -25.0                    # ten cells that anti-correlate violently
    X = np.column_stack([x, y, rng.normal(size=n)])
    expr = _frame(X, ["TF", "G1", "G2"])
    assert np.corrcoef(x, y)[0, 1] < -0.5, "the outliers must flip Pearson's sign in this fixture"
    edges = assoc.bicor(expr, ["TF"])
    row = edges[edges["target"] == "G1"].iloc[0]
    assert row["weight"] > 0.8 and int(row["direction"]) == 1


def test_bicor_falls_back_to_pearson_for_genes_with_zero_mad():
    """More than half the cells at 0 gives MAD 0: no robust scale, so Pearson, not NaN."""
    rng = np.random.default_rng(2)
    n = 200
    tf = rng.normal(size=n)
    sparse = np.where(rng.random(n) < 0.7, 0.0, tf + rng.normal(scale=0.2, size=n))
    expr = _frame(np.column_stack([tf, sparse]), ["TF", "G"])
    got = _weight(assoc.bicor(expr, ["TF"]), "TF", "G")
    assert got == pytest.approx(abs(np.corrcoef(tf, sparse)[0, 1]), abs=1e-9)


def test_bicor_sign_is_reported_as_direction():
    rng = np.random.default_rng(3)
    x = rng.normal(size=100)
    expr = _frame(np.column_stack([x, -x + rng.normal(scale=0.1, size=100)]), ["TF", "G"])
    edges = assoc.bicor(expr, ["TF"])
    assert int(edges[edges["target"] == "G"]["direction"].iloc[0]) == -1


# --------------------------------------------------------------------------- #
# wgcna
# --------------------------------------------------------------------------- #

def _naive_tom(adj: np.ndarray) -> np.ndarray:
    """TOM straight from the definition, a double loop over genes."""
    n = len(adj)
    k = adj.sum(axis=1) - np.diag(adj)
    tom = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            shared = sum(adj[i, u] * adj[u, j] for u in range(n) if u not in (i, j))
            tom[i, j] = (shared + adj[i, j]) / (min(k[i], k[j]) + 1 - adj[i, j])
    return tom


def test_wgcna_tom_matches_the_definition():
    rng = np.random.default_rng(4)
    latent = rng.normal(size=(120, 2))
    X = latent @ rng.normal(size=(2, 9)) + rng.normal(scale=0.7, size=(120, 9))
    genes = [f"g{i}" for i in range(9)]
    expr = _frame(X, genes)
    edges = assoc.wgcna_tom(expr, ["g0", "g1"], power=3)
    adj = np.abs(np.corrcoef(X, rowvar=False)) ** 3
    np.fill_diagonal(adj, 0.0)
    tom = _naive_tom(adj)
    for s, si in (("g0", 0), ("g1", 1)):
        for t in range(2, 9):
            assert _weight(edges, s, genes[t]) == pytest.approx(tom[si, t], abs=1e-10)


def test_scale_free_fit_recovers_an_exact_power_law_and_ignores_flat_connectivity():
    points = np.arange(1, 11)
    k = np.repeat(points, np.round(1e5 / points ** 2).astype(int))      # p(k) proportional to k^-2
    r2, slope = assoc._scale_free_fit(k)
    assert slope == pytest.approx(-2.0, abs=0.05) and r2 > 0.99
    flat = np.repeat(points, 500)                                        # every bin equally full
    assert assoc._scale_free_fit(flat)[0] == pytest.approx(0.0, abs=1e-9)


def test_pick_soft_power_reports_whether_the_cut_was_reached():
    rng = np.random.default_rng(6)
    corr = np.abs(np.corrcoef(rng.normal(size=(80, 60)), rowvar=False))
    power, info = assoc.pick_soft_power(corr, r2_cut=2.0)       # unreachable cut
    assert power in assoc.WGCNA_POWERS and info["cut_reached"] is False
    easy, info2 = assoc.pick_soft_power(corr, r2_cut=-5.0)      # reached by the first power
    assert easy == assoc.WGCNA_POWERS[0] and info2["cut_reached"] is True


def test_wgcna_and_pcnet_refuse_matrices_too_big_for_memory():
    expr = _frame(np.random.default_rng(7).normal(size=(30, 12)), [f"g{i}" for i in range(12)])
    with pytest.raises(MethodError, match="max_genes"):
        assoc.wgcna_tom(expr, ["g0"], max_genes=10)
    with pytest.raises(MethodError, match="max_genes"):
        regr.pcnet(expr, ["g0"], max_genes=10)


def test_default_gene_caps_guard_runtime_and_memory():
    wide = lambda n: _frame(np.zeros((5, n)), [f"g{i}" for i in range(n)])  # noqa: E731
    with pytest.raises(MethodError, match="max_genes=6000"):
        assoc.wgcna_tom(wide(6001), ["g0"])
    with pytest.raises(MethodError, match="max_genes=4000"):
        regr.pcnet(wide(4001), ["g0"])


def test_wgcna_power_is_chosen_without_constant_genes_as_wgcna_does():
    rng = np.random.default_rng(20)
    latent = rng.normal(size=(150, 3))
    X = latent @ rng.normal(size=(3, 40)) + rng.normal(scale=0.8, size=(150, 40))
    genes = [f"g{i}" for i in range(40)]
    plain = assoc.wgcna_tom(_frame(X, genes), ["g0", "g1"])
    with_constant = assoc.wgcna_tom(_frame(np.column_stack([X, np.ones(150)]), genes + ["flat"]), ["g0", "g1"])
    merged = plain.merge(with_constant, on=["source", "target"], suffixes=("_a", "_b"))
    assert len(merged) == len(plain)
    np.testing.assert_allclose(merged["weight_a"], merged["weight_b"], atol=1e-12)


# --------------------------------------------------------------------------- #
# proportionality
# --------------------------------------------------------------------------- #

def test_propr_rho_is_one_minus_scaled_variance_of_the_difference():
    rng = np.random.default_rng(8)
    X = rng.normal(size=(150, 3)) @ np.array([[1.0, 0.4, 0.0], [0.0, 1.0, 0.3], [0.0, 0.0, 1.0]])
    expr = _frame(X, ["TF", "A", "B"])
    edges = assoc.propr_rho(expr, ["TF"])
    for g, j in (("A", 1), ("B", 2)):
        rho = 1.0 - np.var(X[:, 0] - X[:, j], ddof=1) / (np.var(X[:, 0], ddof=1) + np.var(X[:, j], ddof=1))
        assert _weight(edges, "TF", g) == pytest.approx(abs(rho), abs=1e-10)


def test_propr_rho_penalises_a_change_of_scale_that_correlation_ignores():
    rng = np.random.default_rng(9)
    x = rng.normal(size=300)
    expr = _frame(np.column_stack([x, 5.0 * x]), ["TF", "G"])
    rho = _weight(assoc.propr_rho(expr, ["TF"]), "TF", "G")
    pear = _weight(run_method("pearson_signed", _adata(expr.to_numpy(), expr.columns), ["TF"],
                              assume_normalized=True)[0], "TF", "G")
    assert pear == pytest.approx(1.0) and rho < 0.4


# --------------------------------------------------------------------------- #
# ARACNe-style
# --------------------------------------------------------------------------- #

@pytest.fixture(scope="module")
def chain():
    """TF_A -> TF_B -> G (and a lone H): A and G are only associated through B."""
    rng = np.random.default_rng(10)
    n = 800
    a = rng.normal(size=n)
    b = np.tanh(1.5 * a) + rng.normal(scale=0.35, size=n)
    g = np.tanh(1.5 * b) + rng.normal(scale=0.35, size=n)
    h = rng.normal(size=n)
    return _frame(np.column_stack([a, b, g, h]), ["TF_A", "TF_B", "G", "H"])


def test_aracne_dpi_prunes_the_indirect_edge_and_keeps_the_direct_ones(chain):
    edges = assoc.aracne_dpi(chain, ["TF_A", "TF_B"])
    assert _weight(edges, "TF_A", "TF_B") > 0
    assert _weight(edges, "TF_B", "G") > 0
    assert _weight(edges, "TF_A", "G") == 0.0, "A-G is explained by B: the DPI must remove it"
    assert _weight(edges, "TF_B", "H") == 0.0, "an independent gene is not significant"


def test_aracne_dpi_needs_a_third_regulator_to_prune_and_tolerance_loosens_it(chain):
    only_a = assoc.aracne_dpi(chain, ["TF_A"])            # no other TF can be the intermediate
    assert _weight(only_a, "TF_A", "G") > 0
    tolerant = assoc.aracne_dpi(chain, ["TF_A", "TF_B"], tolerance=0.999)
    assert _weight(tolerant, "TF_A", "G") > 0


def test_aracne_dpi_removes_an_edge_that_ties_the_indirect_path():
    """ARACNe prunes on <=: a gene that is an exact copy of TF_B ties A-G with A-B."""
    rng = np.random.default_rng(21)
    a = rng.normal(size=600)
    b = np.tanh(1.5 * a) + rng.normal(scale=0.4, size=600)
    expr = _frame(np.column_stack([a, b, b.copy()]), ["TF_A", "TF_B", "G"])
    edges = assoc.aracne_dpi(expr, ["TF_A", "TF_B"])
    assert _weight(edges, "TF_A", "TF_B") > 0          # no third regulator can explain it
    assert _weight(edges, "TF_A", "G") == 0.0          # equal MI through TF_B: the tie goes


def test_aracne_dpi_threshold_is_a_p_value_not_a_constant(chain):
    loose = assoc.aracne_dpi(chain, ["TF_A", "TF_B"], pvalue=0.5)
    strict = assoc.aracne_dpi(chain, ["TF_A", "TF_B"], pvalue=1e-30)
    assert (loose["weight"] > 0).sum() >= (strict["weight"] > 0).sum()


# --------------------------------------------------------------------------- #
# MRNET
# --------------------------------------------------------------------------- #

def test_mrnet_keeps_one_of_two_redundant_regulators_and_the_complementary_one():
    rng = np.random.default_rng(11)
    n = 1500
    tf1 = rng.normal(size=n)
    tf2 = tf1 + rng.normal(scale=0.05, size=n)            # a copy of TF1
    tf3 = rng.normal(size=n)
    target = np.tanh(tf1) + np.tanh(tf3) + rng.normal(scale=0.3, size=n)
    expr = _frame(np.column_stack([tf1, tf2, tf3, target]), ["TF1", "TF2", "TF3", "T"])
    edges = assoc.mrnet(expr, ["TF1", "TF2", "TF3"])
    w1, w2, w3 = (_weight(edges, f"TF{i}", "T") for i in (1, 2, 3))
    assert w3 > 0
    assert (w1 > 0) != (w2 > 0), "exactly one of the two redundant regulators survives"
    pear = run_method("pearson_signed", _adata(expr.to_numpy(), expr.columns), ["TF1", "TF2", "TF3"],
                      assume_normalized=True)[0]
    assert min(_weight(pear, "TF1", "T"), _weight(pear, "TF2", "T")) > 0.3, \
        "control: plain correlation keeps both"


# --------------------------------------------------------------------------- #
# TIGRESS
# --------------------------------------------------------------------------- #

def test_tigress_scores_true_regulators_high_and_noise_regulators_low():
    rng = np.random.default_rng(12)
    n = 400
    tfs = rng.normal(size=(n, 8))
    target = tfs[:, 0] + 0.8 * tfs[:, 1] + rng.normal(scale=0.5, size=n)
    names = [f"TF{i}" for i in range(8)] + ["T"]
    expr = _frame(np.column_stack([tfs, target]), names)
    edges = regr.tigress(expr, names[:8], n_resamples=40)
    score = {s: _weight(edges, s, "T") for s in names[:8]}
    assert score["TF0"] > 0.7 and score["TF1"] > 0.6
    assert max(score[f"TF{i}"] for i in range(2, 8)) < min(score["TF0"], score["TF1"]) - 0.2
    assert all(0.0 <= v <= 1.0 for v in score.values())


def test_tigress_seed_changes_the_resamples_not_the_conclusion():
    rng = np.random.default_rng(13)
    tfs = rng.normal(size=(300, 6))
    expr = _frame(np.column_stack([tfs, tfs[:, 0] + rng.normal(scale=0.4, size=300)]),
                  [f"TF{i}" for i in range(6)] + ["T"])
    a = regr.tigress(expr, [f"TF{i}" for i in range(6)], n_resamples=20, seed=1)
    b = regr.tigress(expr, [f"TF{i}" for i in range(6)], n_resamples=20, seed=2)
    best = lambda e: e[e["target"] == "T"].sort_values("weight").iloc[-1]["source"]  # noqa: E731
    assert best(a) == best(b) == "TF0"
    assert not np.allclose(a.sort_values(["source", "target"])["weight"],
                           b.sort_values(["source", "target"])["weight"])


# --------------------------------------------------------------------------- #
# BBSR
# --------------------------------------------------------------------------- #

def test_bbsr_selects_the_true_subset_and_nothing_else():
    rng = np.random.default_rng(14)
    n = 600
    tfs = rng.normal(size=(n, 9))
    target = 1.0 * tfs[:, 2] + 0.8 * tfs[:, 5] + rng.normal(scale=0.5, size=n)
    names = [f"TF{i}" for i in range(9)] + ["T"]
    edges = regr.bbsr(_frame(np.column_stack([tfs, target]), names), names[:9])
    chosen = set(edges.loc[(edges["target"] == "T") & (edges["weight"] > 0), "source"])
    assert chosen == {"TF2", "TF5"}


def test_bbsr_weight_of_a_lone_regulator_is_its_explained_variance():
    rng = np.random.default_rng(15)
    n = 500
    x = rng.normal(size=n)
    y = 0.8 * x + rng.normal(scale=0.7, size=n)
    others = rng.normal(size=(n, 3))
    names = ["TF0", "TF1", "TF2", "TF3", "T"]
    edges = regr.bbsr(_frame(np.column_stack([x, others, y]), names), names[:4])
    r2 = float(np.corrcoef(x, y)[0, 1] ** 2)
    assert _weight(edges, "TF0", "T") == pytest.approx(r2, abs=1e-6)


def test_bbsr_refuses_what_it_cannot_do_honestly():
    rng = np.random.default_rng(22)
    names = [f"TF{i}" for i in range(20)] + ["T"]
    wide = _frame(rng.normal(size=(300, 21)), names)
    with pytest.raises(MethodError, match="exceeds 16"):
        regr.bbsr(wide, names[:20], n_preselect=17)
    few_cells = _frame(rng.normal(size=(15, 21)), names)
    with pytest.raises(MethodError, match="at least 2 x n_preselect"):
        regr.bbsr(few_cells, names[:20])


def test_bbsr_reports_nothing_for_a_target_unrelated_to_every_regulator():
    rng = np.random.default_rng(16)
    names = [f"TF{i}" for i in range(5)] + ["T"]
    edges = regr.bbsr(_frame(rng.normal(size=(300, 6)), names), names[:5])
    assert (edges.loc[edges["target"] == "T", "weight"] == 0).all(), "the empty model must win on noise"


# --------------------------------------------------------------------------- #
# pcNet
# --------------------------------------------------------------------------- #

def _svd_pcr(X: np.ndarray, j: int, n_comp: int) -> np.ndarray:
    """Coefficients of gene j regressed on the top PCs of the other (standardised) genes."""
    z, _ = regr._standardize(X)
    others = [i for i in range(z.shape[1]) if i != j]
    _, _, vt = np.linalg.svd(z[:, others], full_matrices=False)
    v = vt[:n_comp].T
    scores = z[:, others] @ v
    return others, v @ ((scores.T @ z[:, j]) / (scores ** 2).sum(axis=0))


@pytest.mark.parametrize("n_genes", [40, 450], ids=["dense_eigh", "lanczos_eigsh"])
def test_pcnet_equals_an_independent_principal_component_regression(n_genes):
    rng = np.random.default_rng(17)
    latent = rng.normal(size=(300, 6))
    X = latent @ rng.normal(size=(6, n_genes)) + rng.normal(scale=0.8, size=(300, n_genes))
    genes = [f"g{i}" for i in range(n_genes)]
    tfs = genes[:5]
    edges = regr.pcnet(_frame(X, genes), tfs, n_comp=3)
    for j in (n_genes - 1, n_genes // 2, 3):
        others, beta = _svd_pcr(X, j, 3)
        for k, g in enumerate(others):
            if genes[g] in tfs:
                expected = abs(beta[k])
                got = float(edges[(edges["source"] == genes[g]) & (edges["target"] == genes[j])]["weight"].iloc[0]) \
                    if genes[g] != genes[j] else None
                assert got == pytest.approx(expected, abs=1e-7)


def test_pcnet_sign_is_the_sign_of_the_coefficient():
    rng = np.random.default_rng(18)
    f = rng.normal(size=(250, 1))
    X = np.hstack([f + rng.normal(scale=0.2, size=(250, 1)), -f + rng.normal(scale=0.2, size=(250, 1)),
                   rng.normal(size=(250, 6))])
    genes = [f"g{i}" for i in range(8)]
    edges = regr.pcnet(_frame(X, genes), ["g0"], n_comp=3)
    assert int(edges[edges["target"] == "g1"]["direction"].iloc[0]) == -1


# --------------------------------------------------------------------------- #
# tree ensembles
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("method", ["extra_trees", "gradient_boosting"])
def test_tree_methods_find_a_non_monotone_dependence_correlation_cannot_see(method):
    rng = np.random.default_rng(19)
    n = 300
    tf = rng.normal(size=(n, 5))
    target = (tf[:, 3] ** 2) + rng.normal(scale=0.2, size=n)      # U-shaped: Pearson ~ 0
    names = [f"TF{i}" for i in range(5)] + ["T"]
    expr = _frame(np.column_stack([tf, target]), names)
    assert abs(np.corrcoef(tf[:, 3], target)[0, 1]) < 0.2
    edges = getattr(regr, method)(expr, names[:5], n_estimators=40)
    top = edges[edges["target"] == "T"].sort_values("weight").iloc[-1]
    assert top["source"] == "TF3"


# --------------------------------------------------------------------------- #
# community
# --------------------------------------------------------------------------- #

def test_community_is_the_mean_rank_of_its_members(planted):
    adata, tfs, _ = planted
    members = ["pearson_signed", "spearman_signed", "mutual_info"]
    got, _ = run_method("community", adata, tf_list=tfs, assume_normalized=True, members=",".join(members))
    ranks = []
    for m in members:
        e, _ = run_method(m, adata, tf_list=tfs, assume_normalized=True)
        ranks.append(e.set_index(["source", "target"])["weight"].rank(ascending=False, method="average"))
    mean_rank = pd.concat(ranks, axis=1).mean(axis=1)
    expected = (1.0 - (mean_rank - 1.0) / len(mean_rank)).rename("expected")
    merged = got.set_index(["source", "target"])["weight"].to_frame().join(expected)
    np.testing.assert_allclose(merged["weight"], merged["expected"], atol=1e-12)
    assert got["weight"].max() <= 1.0 and got["weight"].min() > 0.0


def test_community_refuses_instead_of_silently_dropping_a_member(planted):
    adata, tfs, _ = planted
    with pytest.raises(MethodError, match="not registered"):
        run_method("community", adata, tf_list=tfs, assume_normalized=True, members="pearson_signed,nope")
    with pytest.raises(MethodError, match="itself"):
        run_method("community", adata, tf_list=tfs, assume_normalized=True, members="community,clr")
    with pytest.raises(MethodError, match="at least one"):
        run_method("community", adata, tf_list=tfs, assume_normalized=True, members=" , ")
    ok, why = __import__("grnbench.infer.community", fromlist=["x"]).community_availability("regdiffusion")
    if not ok:                                              # torch is optional: unavailability is data
        assert "regdiffusion" in why


def test_community_forwards_the_seed_to_its_members(planted):
    adata, tfs, _ = planted
    seen = []

    def spy(expr, tf_names=None, seed=None, **_):
        seen.append(seed)
        return run_method("pearson_signed", adata, tf_list=tf_names, assume_normalized=True)[0]

    saved = dict(infer_base.REGISTRY)
    try:
        register_method(MethodSpec("spy_member", "records the seed it receives", spy, default=False))
        run_method("community", adata, tf_list=tfs, assume_normalized=True, members="spy_member", seed=7)
        run_method("community", adata, tf_list=tfs, assume_normalized=True, members="spy_member")
    finally:
        infer_base.REGISTRY.clear()
        infer_base.REGISTRY.update(saved)
    assert seen == [7, None]


def _with_members(planted, outputs, **kw):
    """Run `community` over spy members that each return a fixed edge list."""
    adata, tfs, _ = planted
    saved = dict(infer_base.REGISTRY)
    try:
        names = []
        for i, rows in enumerate(outputs):
            name = f"spy_{i}"
            frame = _frame_edges(rows)
            register_method(MethodSpec(name, "fixed output", lambda expr, tf_names=None, _f=frame, **_: _f.copy(),
                                       default=False))
            names.append(name)
        return run_method("community", adata, tf_list=tfs, assume_normalized=True, members=",".join(names), **kw)
    finally:
        infer_base.REGISTRY.clear()
        infer_base.REGISTRY.update(saved)


def _frame_edges(rows):
    return pd.DataFrame(rows, columns=["source", "target", "weight"])


def test_community_ranks_a_pair_a_sparse_member_leaves_out_last_for_that_member(planted):
    """Review of #69: a member that reports one edge must not give the others rank 1."""
    a = [("TF00", "G000", 1.0)]                                              # sparse member
    b = [("TF00", "G000", 3.0), ("TF00", "G001", 2.0), ("TF00", "G002", 1.0)]
    edges, prov = _with_members(planted, [a, b])
    w = {t: float(edges.loc[edges["target"] == t, "weight"].iloc[0]) for t in ("G000", "G001", "G002")}
    # candidate universe = 3 pairs; ranks: G000 -> (1, 1); G001 -> (3, 2); G002 -> (3, 3)
    assert w["G000"] == pytest.approx(1.0)
    assert w["G001"] == pytest.approx(1.0 - (2.5 - 1.0) / 3)
    assert w["G002"] == pytest.approx(1.0 - (3.0 - 1.0) / 3)
    assert w["G000"] > w["G001"] > w["G002"], "an absent pair must rank below the pairs the member reported"
    assert prov["method_info"]["members"] == ["spy_0", "spy_1"]


def test_community_with_complete_members_is_unchanged_by_the_fill_rule(planted):
    full = [("TF00", f"G00{i}", 3.0 - i) for i in range(3)]
    edges, _ = _with_members(planted, [full, full])
    assert edges["target"].tolist() == ["G000", "G001", "G002"]
    assert edges["weight"].tolist() == pytest.approx([1.0, 1.0 - 1 / 3, 1.0 - 2 / 3])


def test_wgcna_records_the_power_it_used_in_the_provenance(planted):
    adata, tfs, _ = planted
    _, prov = run_method("wgcna_tom", adata, tf_list=tfs, assume_normalized=True)
    info = prov["method_info"]
    assert info["power"] in assoc.WGCNA_POWERS
    assert info["power_selection"].startswith(("scale_free_cut_reached", "fallback_best_signed_r2"))
    assert isinstance(info["signed_scale_free_r2"], float)
    _, given = run_method("wgcna_tom", adata, tf_list=tfs, assume_normalized=True, power=4)
    assert given["method_info"] == {"power": 4, "power_selection": "given"}
    import json
    json.dumps(prov)                                   # the provenance has to be serialisable


def test_community_default_members_are_fixed_and_always_runnable():
    from grnbench.infer.community import DEFAULT_MEMBERS, community_availability
    assert DEFAULT_MEMBERS == ("spearman_signed", "clr", "ppcor", "lasso", "mrnet", "bicor")
    ok, why = community_availability()
    assert ok, why
