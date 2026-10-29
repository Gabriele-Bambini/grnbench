"""Knockdown evaluator on synthetic screens where the right answer is known."""
from itertools import combinations
import json

import numpy as np
import pandas as pd
import pytest

from grnbench.eval.knockdown import (PROP, combine_sets, contrast_test, cross_scores, edge_matrix,
                                     evaluate, main, propensity, responder_matrix)
from grnbench.eval.comparison import holm_adjust
from grnbench.eval.metrics import average_precision, evaluate_ranking, expected_average_precision

N_TF, N_GENES, N_SPECIFIC, N_HUB = 20, 300, 15, 10


def screen(seed=0):
    """20 TFs; each responds on 15 genes of its own plus 10 hub genes shared by every TF."""
    rng = np.random.default_rng(seed)
    genes = [f"G{i}" for i in range(N_GENES)]
    tfs, hubs, others = genes[:N_TF], genes[N_TF:N_TF + N_HUB], genes[N_TF + N_HUB:]
    specific = {t: list(rng.choice(others, N_SPECIFIC, replace=False)) for t in tfs}
    resp = pd.DataFrame([(t, g) for t in tfs for g in specific[t] + hubs], columns=["tf", "gene"])
    return genes, tfs, hubs, specific, resp


def network(genes, tfs, boost, seed):
    """Every TF plus five extra regulators; ``boost(source)`` -> {gene: added weight}."""
    rng = np.random.default_rng(seed)
    sources = tfs + genes[-5:]
    rows = []
    for s in sources:
        w = rng.normal(size=len(genes))
        for g, b in boost(s).items():
            w[genes.index(g)] += b
        rows += [(s, g, w[i]) for i, g in enumerate(genes) if g != s]
    return pd.DataFrame(rows, columns=["source", "target", "weight"])


@pytest.fixture(scope="module")
def result():
    genes, tfs, hubs, specific, resp = screen()
    spec = lambda b, hub=0.0: (lambda s: {**{g: b for g in specific.get(s, [])}, **{g: hub for g in hubs}})
    models = {
        "pearson_signed": network(genes, tfs, spec(1.0), 1),                 # weak TF-specific signal
        "pearson_motif": network(genes, tfs, spec(3.0, hub=5.0), 2),         # motif adds TF-specific signal
        "pearson_motif_genetf": network(genes, tfs, spec(1.0, hub=5.0), 3),  # ... and hubs (its control)
        "grnboost2": network(genes, tfs, spec(1.0), 4),
        "grnboost2_motif": network(genes, tfs, spec(1.0, hub=5.0), 5),       # motif adds hubs only
        "grnboost2_motif_genetf": network(genes, tfs, spec(1.0, hub=5.0), 6),
        "random": network(genes, tfs, spec(0.0), 7),                         # negative control
        "hub_only": network(genes, tfs, lambda s: {g: 5.0 for g in hubs}, 8),
    }
    return evaluate(models, resp, genes, n_perm=2000, seed=0)


def verdicts(res):
    return {r["hypothesis"]: r for r in res["family"]}


def test_positive_control_passes_both_tests(result):
    r = verdicts(result)["H2 pearson_motif > own propensity"]
    assert r["verdict"] == "supported"
    assert r["p_t_greater_holm"] < 0.05 and r["p_perm_greater_holm"] < 0.05


def test_negative_control_fails(result):
    r = verdicts(result)["H2 random > own propensity"]
    assert r["verdict"] == "inconclusive"
    assert r["p_perm_greater"] > 0.05 and r["p_t_greater"] > 0.05


def test_hub_only_network_beats_random_but_not_its_propensity(result):
    # Naively it looks like a strong network: every TF's top targets are true responders (the hubs).
    assert result["arm_means"]["hub_only"] > 3 * result["arm_means"]["random"]
    r = verdicts(result)["H2 hub_only > own propensity"]
    assert r["verdict"] != "supported"
    assert r["p_perm_greater"] > 0.05


def test_motif_lift_from_hubs_is_caught_but_a_real_lift_is_not(result):
    v = verdicts(result)
    assert v["H1a pearson_motif > pearson_signed"]["verdict"] == "supported"
    assert v["H1b pearson_motif > pearson_motif_genetf"]["verdict"] == "supported"
    # the raw lift is fooled by hub genes; the TF-generic motif control is not
    assert v["H1a grnboost2_motif > grnboost2"]["verdict"] == "supported"
    assert v["H1b grnboost2_motif > grnboost2_motif_genetf"]["verdict"] == "inconclusive"
    assert v["H1b grnboost2_motif > grnboost2_motif_genetf"]["p_perm_greater"] > 0.05


def test_family_size_and_holm_over_both_directions(result):
    assert len(result["family"]) == 2 + 2 + 6        # H1a, H1b, H2 (no H2 for the _genetf controls)
    assert result["n_directed_hypotheses"] == 20
    for r in result["family"]:
        assert r["p_t_greater_holm"] >= r["p_t_greater"] and r["p_perm_greater_holm"] >= r["p_perm_greater"]
    json.dumps({k: v for k, v in result.items() if k != "per_unit"})


def test_unit_endpoint_matches_evaluate_ranking():
    genes, tfs, _, _, resp = screen()
    edges = network(genes, tfs, lambda s: {}, 7)
    res = evaluate({"m": edges}, resp, genes, family=[], secondary=[], n_perm=1)
    tf = tfs[3]
    gold = {(tf, g) for g in resp.loc[resp.tf == tf, "gene"]}
    universe = pd.DataFrame({"source": tf, "target": [g for g in genes if g != tf]})
    ref = evaluate_ranking(edges[edges.source == tf], gold, with_curve=False, universe=universe)
    assert res["per_unit"].loc[tf, "m"] == pytest.approx(ref["auprc_ratio"])


def test_missing_edges_share_one_bottom_tie():
    genes = ["A", "B", "C", "D", "E"]
    edges = pd.DataFrame({"source": ["A"], "target": ["C"], "weight": [1.0]})
    mat = edge_matrix(edges, ["A"], genes)
    assert mat[0, 2] == 1.0 and np.isneginf(np.delete(mat[0], 2)).all()
    resp = responder_matrix(pd.DataFrame({"tf": ["A", "A"], "gene": ["C", "A"]}), ["A"], genes)
    assert resp[0].tolist() == [False, False, True, False, False]       # own gene dropped
    c = cross_scores(mat, resp, np.array([0]), "precision_at_k", k=1)
    assert c[0, 0] == 1.0
    c = cross_scores(np.full((1, 5), -np.inf), resp, np.array([0]), "precision_at_k", k=1)
    assert c[0, 0] == pytest.approx(0.25)                               # expected value inside the tie


def test_propensity_is_leave_one_tf_out():
    m = np.array([[3., 2., 1.], [1., 2., 3.], [1., 2., 3.]])
    p = propensity(m, ["X", "Y", "Z"], ["X", "Q"])
    pct = np.array([[1, 2 / 3, 1 / 3], [1 / 3, 2 / 3, 1], [1 / 3, 2 / 3, 1]])
    assert p[0] == pytest.approx(pct[1:].mean(axis=0))      # X's own row left out
    assert p[1] == pytest.approx(pct.mean(axis=0))          # Q is no regulator: all rows


def test_exclusions_are_reported():
    genes, tfs, _, _, resp = screen()
    extra = pd.DataFrame({"tf": ["G5"] * 3 + ["NOT_MEASURED"] * 20,
                          "gene": genes[100:103] + genes[100:120]})
    resp = pd.concat([resp[resp.tf != "G5"], extra])
    res = evaluate({"m": network(genes, tfs, lambda s: {}, 8)}, resp, genes, family=[], secondary=[], n_perm=1)
    assert res["n_units"] == N_TF - 1
    assert "3 responders" in res["excluded"]["G5"] and "NOT_MEASURED" in res["excluded"]


def test_contrast_validation():
    cross = {"a": np.eye(4), "b": np.zeros((4, 4))}
    with pytest.raises(ValueError, match="not scored"):
        contrast_test(cross, {"a": 1, "c": -1})
    with pytest.raises(ValueError, match="n_perm"):
        contrast_test(cross, {"a": 1, "b": -1}, n_perm=0)
    r = contrast_test(cross, {"a": 1, "b": -1}, n_perm=4000)
    assert r["mean_contrast"] == 1 and r["status"] == "zero_variance_unresolved"
    # only the identity permutation (1 in 4! = 24) reaches the observed value
    assert r["p_perm_greater"] == pytest.approx(1 / 24, abs=0.01)


def test_combine_sets():
    d = {"family": [{"hypothesis": h, "verdict": v} for h, v in
                    [("a", "supported"), ("b", "supported"), ("c", "inconclusive"), ("d", "refuted")]]}
    r = {"family": [{"hypothesis": h, "verdict": v} for h, v in
                    [("a", "supported"), ("b", "inconclusive"), ("c", "supported"), ("d", "refuted")]]}
    assert combine_sets(d, r) == {"a": "supported (replicated)", "b": "supported in design, inconclusive in replication",
                                  "c": "inconclusive in design, supported in replication", "d": "refuted (replicated)"}


def test_combine_sets_is_symmetric_and_keeps_hypotheses_of_either_screen():
    fam = lambda **v: {"family": [{"hypothesis": h, "verdict": x} for h, x in v.items()]}
    out = combine_sets(fam(a="inconclusive", b="supported", c="inconclusive"),
                       fam(a="supported", b="inconclusive", c="inconclusive", extra="supported"))
    assert out == {"a": "inconclusive in design, supported in replication",
                   "b": "supported in design, inconclusive in replication",
                   "c": "inconclusive", "extra": "missing in design, supported in replication"}
    assert list(out) == ["a", "b", "c", "extra"]                       # design order first, then replication-only
    assert combine_sets(fam(a="supported"), fam())["a"] == "supported in design, missing in replication"
    assert combine_sets(fam(a="missing"), fam(a="missing"))["a"] == "missing"


def test_expected_average_precision():
    rng = np.random.default_rng(3)
    s = rng.normal(size=300)
    y = (rng.random(300) < 0.1).astype(int)
    order = np.argsort(-s)
    assert expected_average_precision(s[order], y[order]) == pytest.approx(average_precision(y[order]), abs=1e-12)
    assert expected_average_precision(s, y) == pytest.approx(average_precision(y[order]), abs=1e-12)   # unsorted input
    # a long tied tail: closed form equals the mean over random tie-breaks of the plain AP
    n, ranked = 400, 60
    w = np.r_[np.arange(ranked, 0, -1).astype(float), np.full(n - ranked, -np.inf)]
    y = np.zeros(n, dtype=int)
    y[rng.choice(n, 25, replace=False)] = 1
    draws = []
    for _ in range(20000):
        yb = y.copy()
        yb[ranked:] = rng.permutation(y[ranked:])
        draws.append(average_precision(yb))
    assert expected_average_precision(w, y) == pytest.approx(np.mean(draws), abs=0.002)
    assert expected_average_precision(w, y) > average_precision(y, w)      # the block convention is below the expectation
    # all tied: not the prevalence (that is the block convention); with one relevant item, the mean of 1/rank
    one = np.zeros(10, dtype=int)
    one[4] = 1
    assert expected_average_precision(np.zeros(10), one) == pytest.approx(np.mean(1 / np.arange(1, 11)))
    n, r = 6, 2
    labels = np.r_[np.ones(r, dtype=int), np.zeros(n - r, dtype=int)]
    tied = expected_average_precision(np.zeros(n), labels)
    aps = []
    for responders in combinations(range(n), r):
        ordered = np.zeros(n, dtype=int)
        ordered[list(responders)] = 1
        aps.append(np.mean(np.cumsum(ordered)[ordered == 1] / (np.flatnonzero(ordered) + 1)))
    assert tied == pytest.approx(np.mean(aps), abs=1e-12)
    assert tied / (r / n) > 1
    full = np.ones(10, dtype=int)
    assert expected_average_precision(np.zeros(10), full) == pytest.approx(1.0)
    assert np.isnan(expected_average_precision(np.zeros(3), np.zeros(3)))


def test_cli_writes_tables(tmp_path):
    genes, tfs, _, specific, resp = screen()
    edges = network(genes, tfs, lambda s: {g: 3.0 for g in specific.get(s, [])}, 9)
    edges.to_csv(tmp_path / "e.tsv", sep="\t", index=False)
    resp.to_csv(tmp_path / "r.tsv", sep="\t", index=False)
    (tmp_path / "g.txt").write_text("\n".join(genes), encoding="utf-8")
    assert main(["--edges", f"pearson_signed={tmp_path / 'e.tsv'}", "--responders", str(tmp_path / "r.tsv"),
                 "--genes", str(tmp_path / "g.txt"), "--out", str(tmp_path / "o"), "--n-perm", "200"]) == 0
    fam = json.loads((tmp_path / "o" / "family_auprc_ratio.json").read_text(encoding="utf-8"))
    assert fam["family"][0]["verdict"] == "supported"
    per_unit = pd.read_csv(tmp_path / "o" / "per_unit_auprc_ratio.tsv", sep="\t")
    assert {"tf", "n_responders", "pearson_signed", PROP + "pearson_signed", "response_propensity"} <= set(per_unit.columns)


def test_missing_arm_stays_in_the_family_as_p_one_and_a_typo_still_raises():
    genes, tfs, _, _, resp = screen()
    edges = {"m": network(genes, tfs, lambda s: {}, 7)}
    fam = [("H2 m > own propensity", {"m": 1, PROP + "m": -1}), ("H2 gone > own propensity", {"gone": 1, PROP + "gone": -1})]
    res = evaluate(edges, resp, genes, family=fam, secondary=[], missing_arms={"gone"}, n_perm=50)
    assert res["n_directed_hypotheses"] == 4
    gone = res["family"][1]
    assert gone["verdict"] == "missing" and gone["unavailable"] == ["gone"] and np.isnan(gone["mean_contrast"])
    assert gone["p_t_greater_holm"] == gone["p_perm_less_holm"] == 1.0
    m = res["family"][0]
    assert [m["p_t_greater_holm"], m["p_t_less_holm"]] == pytest.approx(
        holm_adjust([m["p_t_greater"], m["p_t_less"], 1.0, 1.0])[:2])     # multiplier of the full family
    with pytest.raises(ValueError, match="not scored"):
        evaluate(edges, resp, genes, family=fam, secondary=[], n_perm=50)
