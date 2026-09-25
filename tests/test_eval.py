"""Evaluation module tests: hand-computed metrics, nulls, gold loader, stability."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from grnbench.eval.gold import GoldStandard, coverage_report, load_gold
from grnbench.eval.metrics import (
    GoldIndex,
    average_precision,
    auc_roc,
    evaluate_ranking,
    identifier_overlap,
    label_edges,
    pr_curve,
    precision_at_k,
    top_k_overlap,
)
from grnbench.eval.nulls import degree_preserving_nulls, null_metrics
from grnbench.eval.stability import subsample_stability


def _ranked_edges():
    """10 edges ranked by weight; gold hits exactly at ranks 1 and 4."""
    rows = [
        ("TF1", "G1", 10.0), ("TF2", "G2", 9.0), ("TF3", "G3", 8.0),
        ("TF1", "G4", 7.0),                              # <- gold hit
        ("TF4", "G5", 6.0), ("TF5", "G6", 5.0), ("TF6", "G7", 4.0),
        ("TF7", "G8", 3.0), ("TF8", "G9", 2.0), ("TF9", "G10", 1.0),
    ]
    return pd.DataFrame(rows, columns=["source", "target", "weight"])


def _gold_pairs():
    return {("TF1", "G1"), ("TF1", "G4")}


def test_label_edges():
    y = label_edges(_ranked_edges(), _gold_pairs())
    assert y.tolist() == [1, 0, 0, 1, 0, 0, 0, 0, 0, 0]


def test_hand_computed_metric_block():
    res = evaluate_ranking(_ranked_edges(), _gold_pairs())
    # AP = (1/1 + 2/4) / 2 = 0.75 ; prevalence = 2/10 ; ratio = 3.75
    assert res["auprc"] == pytest.approx(0.75)
    assert res["auprc_ratio"] == pytest.approx(3.75)
    # k = n_gold = 2 -> precision@2 = 1/2 -> EPR = 0.5 / 0.2 = 2.5
    assert res["precision_at_k"] == pytest.approx(0.5)
    assert res["epr"] == pytest.approx(2.5)
    # AUROC by ranks: pos ranks (asc) 10 and 7 -> U=(17-3)/(2*8)=0.875
    assert res["auroc"] == pytest.approx(0.875)
    assert res["n_edges"] == 10 and res["n_gold_in_list"] == 2


def test_average_precision_and_precision_at_k_edges():
    y = np.array([1, 0, 0, 1, 0])
    assert average_precision(y) == pytest.approx((1.0 + 2 / 4) / 2)
    assert precision_at_k(y, 2) == pytest.approx(0.5)
    assert np.isnan(average_precision(np.zeros(4, dtype=int)))


def test_auc_roc_ties_are_average_ranked():
    w = np.array([5.0, 5.0, 1.0])
    y = np.array([1, 0, 0])
    # ranks asc: both 5.0 share rank 2.5 -> U=(2.5-1)/(1*2)=0.75
    assert auc_roc(w, y) == pytest.approx(0.75)
    assert np.isnan(auc_roc(np.array([1.0]), np.array([1])))


def test_evaluate_empty_raises():
    with pytest.raises(ValueError):
        evaluate_ranking(pd.DataFrame(columns=["source", "target", "weight"]), set())


def test_degree_preserving_null_keeps_out_degrees():
    edges = pd.DataFrame({
        "source": ["A"] * 3 + ["B"] * 2,
        "target": ["x", "y", "z", "w", "v"],
        "weight": [5, 4, 3, 2, 1][:5],
    })
    for null in degree_preserving_nulls(edges, n=3, seed=0):
        assert null["source"].value_counts().to_dict() == {"A": 3, "B": 2}
        assert set(null["target"]) == {"x", "y", "z", "w", "v"}


def test_null_metrics_strong_method_beats_score_null():
    res = null_metrics(_ranked_edges(), _gold_pairs(), n_nulls=30,
                       mode="random", seed=1)
    assert res["observed"] == pytest.approx(0.75)
    assert res["null_mean"] < res["observed"]
    assert 0.0 < res["p_value"] <= 1.0
    # degenerate nulls (all gold hits lost) are excluded and counted, honestly
    assert res["n_nulls_used"] + res["n_nulls_degenerate"] == 30
    assert res["n_nulls_used"] >= 1


def test_load_gold_file_with_alias_headers(tmp_path):
    p = tmp_path / "gold.tsv"
    p.write_text("Regulator\tTarget\tConfidence\nTF1\tG1\tA\nTF1\tG4\tB\n", encoding="utf-8")
    gold = load_gold(str(p))
    assert isinstance(gold, GoldStandard)
    assert gold.pairs == frozenset({("TF1", "G1"), ("TF1", "G4")})


def test_load_gold_reads_the_beeline_gene1_gene2_convention(tmp_path):
    """Regression: every BEELINE ground-truth network was unreadable.

    BEELINE ships its ChIP-seq / STRING / lof-gof networks as `Gene1,Gene2`, and
    that convention propagated to most published GRN ground truths. The alias
    table accepted `gene2` as a target but not `gene1` as a source, so all ten
    BEELINE networks raised "need source/regulator + target columns".
    """
    p = tmp_path / "hESC-ChIP-seq-network.csv"
    p.write_text("Gene1,Gene2\nPOU5F1,NANOG\nSOX2,POU5F1\n", encoding="utf-8")
    gold = load_gold(str(p))
    assert gold.pairs == frozenset({("POU5F1", "NANOG"), ("SOX2", "POU5F1")})

    # some BEELINE networks carry a third Score column, which is not a filter:
    # every listed edge is a positive, exactly as BEELINE treats them
    scored = tmp_path / "mHSC-ChIP-seq-network.csv"
    scored.write_text("Gene1,Gene2,Score\nGATA1,KLF1,880\nGATA1,ZFPM1,12\n",
                      encoding="utf-8")
    assert len(load_gold(str(scored))) == 2


def test_load_gold_does_not_mistake_a_gzipped_csv_for_tsv(tmp_path):
    """`.csv.gz` must be read as CSV, not as one wide tab-separated column."""
    import gzip

    p = tmp_path / "net.csv.gz"
    with gzip.open(p, "wt", encoding="utf-8") as fh:
        fh.write("Gene1,Gene2\nTF1,G1\n")
    assert load_gold(str(p)).pairs == frozenset({("TF1", "G1")})


def test_load_gold_error_names_the_accepted_columns(tmp_path):
    p = tmp_path / "weird.csv"
    p.write_text("colA,colB\nx,y\n", encoding="utf-8")
    with pytest.raises(ValueError, match="regulator"):
        load_gold(str(p))
    try:
        load_gold(str(p))
    except ValueError as e:
        assert "gene1" in str(e), "the error must list what IS accepted"


def test_load_gold_unknown_spec_raises():
    with pytest.raises(FileNotFoundError, match="dorothea"):
        load_gold("mystandards_v9")


def test_coverage_report(tmp_path):
    edges = _ranked_edges()
    cov = coverage_report({"demo": edges}, GoldStandard("g", frozenset(_gold_pairs()), "test"))
    assert cov["demo"]["gold_reachable"] == 2
    assert cov["demo"]["gold_pairs"] == 2


def test_subsample_stability_on_cascade(cascade_adata):
    res = subsample_stability(cascade_adata, "spearman_signed",
                              n_subsamples=3, frac=0.75, top_k=20, seed=5)
    assert 0.0 <= res["jaccard_mean"] <= 1.0
    assert len(res["jaccard_all"]) == 3  # C(3,2) pairs
    assert all(c > 0 for c in res["edge_counts"])


# --- gene-identifier matching -------------------------------------------------
#
# BEELINE ships mESC expression in MGI case ("Pou5f1") and every mouse network
# in upper case ("POU5F1"). Exact matching scored 64,557 reachable gold edges as
# zero hits and reported NaN without a single line of output; these tests pin
# down that this cannot happen silently again.

def _mgi_edges():
    return pd.DataFrame(
        [("Pou5f1", "Nanog", 3.0), ("Sox2", "Pou5f1", 2.0), ("Nanog", "Sox2", 1.0)],
        columns=["source", "target", "weight"])


def _upper_gold():
    return {("POU5F1", "NANOG"), ("SOX2", "POU5F1")}


def test_case_folded_matching_recovers_beeline_mesc_hits():
    assert label_edges(_mgi_edges(), _upper_gold()).tolist() == [1, 1, 0]


def test_case_sensitive_opt_out_still_available():
    y = label_edges(_mgi_edges(), _upper_gold(), case_sensitive=True)
    assert y.tolist() == [0, 0, 0]


def test_evaluate_ranking_scores_across_case():
    res = evaluate_ranking(_mgi_edges(), _upper_gold(), with_curve=False)
    assert res["n_gold_in_list"] == 2
    assert not np.isnan(res["auprc"])


def test_identifier_overlap_flags_case_mismatch():
    d = identifier_overlap(_mgi_edges(), _upper_gold())
    assert d["verdict"] == "case-mismatch"
    assert d["matched_pairs_case_sensitive"] == 0
    assert d["matched_pairs_case_insensitive"] == 2
    assert d["shared_genes_case_insensitive"] == 3
    assert d["example_gold_genes"]  # the message must be able to show both sides


def test_identifier_overlap_flags_namespace_mismatch():
    edges = pd.DataFrame([("ENSG001", "ENSG002", 1.0)],
                         columns=["source", "target", "weight"])
    d = identifier_overlap(edges, _upper_gold())
    assert d["verdict"] == "no-overlap"
    assert d["shared_genes_case_insensitive"] == 0


def test_identifier_overlap_clean_when_cases_agree():
    assert identifier_overlap(_ranked_edges(), _gold_pairs())["verdict"] == "ok"


def test_case_sensitive_index_cannot_hide_the_mismatch():
    """The diagnostic must not inherit the setting it exists to second-guess."""
    gi = GoldIndex(_upper_gold(), case_sensitive=True)
    assert evaluate_ranking(_mgi_edges(), gi, with_curve=False)["n_gold_in_list"] == 0
    assert identifier_overlap(_mgi_edges(), gi)["verdict"] == "case-mismatch"


def test_gold_index_is_reused_not_rebuilt():
    gi = GoldIndex(_upper_gold())
    assert GoldIndex.coerce(gi) is gi


def test_nulls_accept_a_prebuilt_index():
    res = null_metrics(_ranked_edges(), GoldIndex(_gold_pairs()), n_nulls=10,
                       mode="degree", seed=1)
    assert res["n_nulls_used"] + res["n_nulls_degenerate"] == 10


def test_coverage_report_folds_case():
    gold = GoldStandard("g", frozenset(_upper_gold()), "test")
    cov = coverage_report({"demo": _mgi_edges()}, gold)
    assert cov["demo"]["gold_reachable"] == 2
    strict = coverage_report({"demo": _mgi_edges()}, gold, case_sensitive=True)
    assert strict["demo"]["gold_reachable"] == 0


# --------------------------------------------------------------------------- #
# regressions found by the BEELINE backtest / adversarial audit
# --------------------------------------------------------------------------- #

def _misses_every_gold():
    """A ranking whose candidate space contains no gold pair, but whose
    degree-preserving shuffles land on plenty of them."""
    edges = pd.DataFrame({
        "source": [f"TF{i % 5}" for i in range(200)],
        "target": [f"G{i}" for i in range(200)],
        "weight": np.linspace(1.0, 0.0, 200),
    })
    gold = {(f"TF{(i + 1) % 5}", f"G{i}") for i in range(200)}
    return edges, gold


def test_undefined_observation_is_not_a_significant_p_value():
    """A method that reached no gold edge scored nothing; it did not win.

    `valid >= nan` is elementwise False, so the naive estimator returned
    1/(n+1) - the smallest p-value it can produce - and the HTML report
    rendered that as a green 'significant' badge beside empty metric cells.
    """
    edges, gold = _misses_every_gold()
    assert np.isnan(evaluate_ranking(edges, gold, with_curve=False)["auprc"])
    res = null_metrics(edges, gold, n_nulls=45, mode="degree", seed=0)
    assert res["n_nulls_used"] > 0, "the nulls themselves must be scorable"
    assert res["observed_undefined"] is True
    assert np.isnan(res["p_value"])


def test_a_real_observation_still_gets_a_p_value():
    res = null_metrics(_ranked_edges(), _gold_pairs(), n_nulls=20, mode="random",
                       seed=0)
    assert res["observed_undefined"] is False
    assert not np.isnan(res["p_value"])


def test_top_k_overlap_folds_case_like_gold_matching():
    """Two identical networks differing only in symbol case are not disjoint."""
    a = pd.DataFrame({"source": ["Pou5f1", "Sox2", "Nanog"],
                      "target": ["Sox2", "Nanog", "Klf4"],
                      "weight": [3.0, 2.0, 1.0]})
    b = a.assign(source=a["source"].str.upper(), target=a["target"].str.upper())
    assert top_k_overlap({"mine": a, "theirs": b})["jaccard"]["mine"]["theirs"] == 1.0
    with pytest.raises(ValueError, match='identical candidate supports'):
        top_k_overlap({"mine": a, "theirs": b}, case_sensitive=True)


def _tied_ranking(seed: int):
    """400 edges over 4 score levels; only the order WITHIN a level changes."""
    n = 400
    pairs = [(f"TF{i % 20}", f"G{i}") for i in range(n)]
    gold = set(pairs[:40])
    df = pd.DataFrame({"source": [p[0] for p in pairs],
                       "target": [p[1] for p in pairs],
                       "weight": np.repeat([3.0, 2.0, 1.0, 0.0], n // 4)})
    shuffled = (df.sample(frac=1.0, random_state=seed)
                  .sort_values("weight", ascending=False, kind="stable")
                  .reset_index(drop=True))
    return shuffled, gold


def test_metrics_do_not_depend_on_the_order_of_tied_edges():
    """Equal scores assert nothing about which edge outranks which.

    Before this was fixed, six orderings of the SAME scores gave AUPRC 0.319 to
    0.478 and EPR 3.0 to 5.0, while the tie-safe AUROC never moved - proof the
    ranking was identical and only the tie-breaking was being measured.
    """
    results = [evaluate_ranking(*_tied_ranking(s), with_curve=False)
               for s in range(6)]
    for key in ("auprc", "precision_at_k", "epr", "auroc"):
        values = [r[key] for r in results]
        assert max(values) - min(values) < 1e-12, f"{key} still tie-sensitive"


def test_tie_ambiguity_is_reported_not_hidden():
    res = evaluate_ranking(*_tied_ranking(0), with_curve=False)
    lo, hi = res["epr_range"]
    assert res["tie_block_at_k"] > 0, "the cut at k falls inside a tied block"
    assert lo < res["epr"] < hi, "the reported EPR sits inside its own interval"
    # All 40 gold edges sit inside the first tied block of 100, and k is 40, so
    # tie-breaking could put anywhere from 0 to 40 of them above the cut. The
    # reported value is the expectation under a uniform random break -
    # 40 * 40/100 = 16 hits, precision 0.4, EPR 0.4 / 0.1 - which is NOT the
    # midpoint of [0, 10]: the endpoints are extremes, not equally likely.
    assert (lo, hi) == pytest.approx((0.0, 10.0))
    assert res["epr"] == pytest.approx(4.0)


def test_distinct_scores_are_unaffected_by_tie_awareness():
    """No ties -> the tie-aware path must reproduce the textbook definition."""
    edges, gold = _ranked_edges(), _gold_pairs()
    assert edges["weight"].is_unique
    res = evaluate_ranking(edges, gold, with_curve=False)
    y = label_edges(edges, gold)
    assert res["auprc"] == pytest.approx(average_precision(y))
    assert res["precision_at_k"] == pytest.approx(precision_at_k(y, res["epr_k"]))
    assert res["tie_block_at_k"] == 0
    assert res["epr_range"][0] == pytest.approx(res["epr_range"][1])


def test_pr_curve_collapses_ties_without_running_off_the_end():
    edges, gold = _tied_ranking(0)
    c = pr_curve(label_edges(edges, gold), weights=edges["weight"].to_numpy())
    assert len(c["recall"]) == len(c["precision"]) == 4   # one point per level
    assert c["recall"][-1] == pytest.approx(1.0)


def test_label_permutation_null_equals_score_permutation_null_with_ties():
    """The fast null permutes labels over fixed scores; it must equal permuting
    the scores over fixed pairs and re-scoring - including inside tie blocks."""
    from grnbench.eval.metrics import metric_from_labels
    from grnbench.eval.universe import align_ranking
    rng = np.random.default_rng(11)
    n = 60
    edges = pd.DataFrame({"source": [f"TF{i % 7}" for i in range(n)],
                          "target": [f"G{i}" for i in range(n)],
                          "weight": rng.integers(0, 6, size=n).astype(float)})  # heavy ties
    gold = {(f"TF{i % 7}", f"G{i}") for i in rng.choice(n, 12, replace=False)}
    gi = GoldIndex(gold)
    aligned = align_ranking(edges)
    w, y = aligned.weight.to_numpy(float), gi.labels(aligned).astype(np.int64)
    k = 7
    for _ in range(5):
        perm = rng.permutation(n)
        slow = evaluate_ranking(aligned.assign(weight=w[perm]), gi, top_k=k, with_curve=False)
        inv = np.argsort(perm)          # slot j now hosts the pair that received w[j]
        for metric in ("auprc", "auprc_ratio", "epr", "precision_at_k", "auroc"):
            assert metric_from_labels(metric, w, y[inv], k) == pytest.approx(slow[metric])


def test_regulator_null_keeps_label_counts_within_each_source():
    from grnbench.eval.universe import align_ranking
    n = 40
    edges = pd.DataFrame({"source": [f"TF{i % 4}" for i in range(n)],
                          "target": [f"G{i}" for i in range(n)],
                          "weight": np.arange(n, dtype=float)[::-1]})
    gold = {(f"TF{i % 4}", f"G{i}") for i in (0, 1, 2, 5, 9, 13, 17)}   # all in TF0/TF1
    res = null_metrics(edges, gold, n_nulls=25, mode="regulator", seed=3, metric="auroc")
    assert res["inference_available"] and res["n_nulls_used"] == 25
    # every draw keeps the hits inside TF0/TF1, so the null AUROC stays away from
    # the global-permutation value; sanity: observed is a real number
    assert np.isfinite(res["observed"]) and np.isfinite(res["null_mean"])
