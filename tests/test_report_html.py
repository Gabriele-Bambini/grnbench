"""HTML comparison report: content, escaping, and graceful degradation."""
from __future__ import annotations

import math
import re

import pandas as pd
import pytest

from grnbench.eval.metrics import evaluate_ranking, pr_curve, top_k_overlap
from grnbench.html_report import render_html


def _edges(pairs):
    return pd.DataFrame(pairs, columns=["source", "target", "weight"])


def _report():
    a = _edges([("TF1", "G1", 10.0), ("TF2", "G2", 9.0),
                ("TF1", "G4", 8.0), ("TF3", "G3", 7.0)])
    b = _edges([("TF2", "G2", 5.0), ("TF3", "G3", 4.0),
                ("TF1", "G1", 3.0), ("TF1", "G4", 2.0)])
    gold = {("TF1", "G1"), ("TF1", "G4")}
    return {
        "created_utc": "2026-09-05T10:00:00+00:00",
        "gold": "toy-gold", "gold_origin": "test", "gold_pairs": 2,
        "null_mode": "degree", "n_nulls": 10,
        "methods": {
            "alpha": {**evaluate_ranking(a, gold), "file": "alpha.edges.tsv.gz",
                      "null": {"metric": "auprc", "observed": 0.8, "null_mean": 0.3,
                               "null_std": 0.1, "p_value": 0.02,
                               "n_nulls_used": 8, "n_nulls_degenerate": 2}},
            "beta": {**evaluate_ranking(b, gold), "file": "beta.edges.tsv.gz",
                     "null": {"metric": "auprc", "observed": 0.4, "null_mean": 0.35,
                              "null_std": 0.1, "p_value": 0.4,
                              "n_nulls_used": 10, "n_nulls_degenerate": 0},
                     "stability": {"method": "beta", "n_subsamples": 3,
                                   "subsample_frac": 0.8, "top_k": 100,
                                   "jaccard_mean": 0.42, "jaccard_std": 0.05,
                                   "jaccard_all": [0.4, 0.45, 0.41],
                                   "edge_counts": [4, 4, 4]}},
        },
        "coverage": {"alpha": {"gold_pairs": 2, "gold_reachable": 2,
                               "tf_coverage_in_data": 1},
                     "beta": {"gold_pairs": 2, "gold_reachable": 2,
                              "tf_coverage_in_data": 1}},
        "overlap": top_k_overlap({"alpha": a, "beta": b}, k=4),
        "plugins": [{"name": "beta", "origin": "command", "requested_as": "b.json"}],
        "environment": {"grnbench": "0.5.0", "python": "3.11.15"},
    }


def test_report_is_standalone_and_covers_every_section():
    html = render_html(_report())
    assert html.startswith("<!doctype html>")
    for section in ("Verdict", "Precision", "Against the null",
                    "Do the methods agree", "Subsample stability",
                    "Gold reachability", "Preprocessing QC", "Provenance"):
        assert section in html, f"missing section: {section}"
    assert "<svg" in html and "polyline" in html      # PR curves drawn inline
    assert "alpha" in html and "beta" in html
    # nothing may be fetched at render time: fully self-contained
    assert "http://" not in html and "https://" not in html
    assert "<script" not in html


def test_degenerate_nulls_are_surfaced_not_hidden():
    html = render_html(_report())
    assert "Degenerate nulls counted" in html
    assert "2 shuffle(s)" in html


def test_stability_absence_tells_the_user_how_to_get_it():
    rep = _report()
    rep["methods"]["beta"].pop("stability")
    html = render_html(rep)
    assert "--stability-data" in html


def test_prep_manifest_gates_are_rendered_with_status_badges():
    prep = {
        "strategy_chosen": "cpm_log1p",
        "forced_past_gates": True,
        "warnings": ["TF coverage low: 2/100"],
        "qc": {"gates": [
            {"name": "mito_fraction", "status": "fail", "detail": "median 40%"},
            {"name": "sparsity", "status": "warn", "detail": "96%"},
            {"name": "batch_annotation", "status": "pass", "detail": "column 'donor'"},
        ]},
        "steps": [{"step": "ingest", "detail": "format=h5ad"},
                  {"step": "filter_cells", "n_cells_after": 380,
                   "n_genes_after": 600, "detail": "removed 20"}],
        "environment": {"scanpy": "1.11.5"},
    }
    html = render_html(_report(), prep)
    assert 'class="badge fail">FAIL' in html
    assert 'class="badge warn">WARN' in html
    assert "run forced past failing gates" in html
    assert "TF coverage low" in html
    assert "Preprocessing steps" in html


def test_html_escapes_hostile_method_names():
    rep = _report()
    rep["methods"]["<img src=x onerror=alert(1)>"] = rep["methods"].pop("alpha")
    rep["gold"] = "</title><script>alert(2)</script>"
    html = render_html(rep)
    assert "<img src=x" not in html
    assert "<script>alert(2)</script>" not in html
    assert "&lt;img src=x" in html


def test_report_survives_a_method_with_no_gold_hits():
    rep = _report()
    rep["methods"]["ghost"] = {"n_edges": 3, "n_gold_in_list": 0, "prevalence": 0.0,
                               "auprc": float("nan"), "auprc_ratio": float("nan"),
                               "epr_k": 0, "precision_at_k": float("nan"),
                               "epr": float("nan"), "auroc": float("nan"),
                               "pr_curve": pr_curve([0, 0, 0]), "null": None}
    html = render_html(rep)
    assert "ghost" in html and "&ndash;" in html   # NaNs render as an en dash


def test_empty_report_does_not_crash():
    html = render_html({"gold": "none", "gold_pairs": 0, "methods": {}})
    assert "No method produced a scorable edge list" in html


def test_pr_curve_shape_and_baseline():
    y = [1, 0, 0, 1, 0, 0, 0, 0, 0, 0]
    c = pr_curve(y)
    assert c["n_pos"] == 2
    assert c["baseline"] == pytest.approx(0.2)
    assert c["recall"][-1] == pytest.approx(1.0)
    assert c["precision"][0] == pytest.approx(1.0)
    assert len(c["recall"]) == len(c["precision"]) == 10
    assert pr_curve([0, 0, 0])["recall"] == []


def test_pr_curve_is_thinned_for_big_rankings():
    y = [1] + [0] * 4999
    c = pr_curve(y, max_points=50)
    assert len(c["recall"]) <= 50


def test_top_k_overlap_is_symmetric_and_bounded():
    a = _edges([("A", "1", 3.0), ("A", "2", 2.0), ("A", "3", 1.0)])
    b = _edges([("A", "2", 3.0), ("A", "3", 2.0), ("A", "9", 1.0)])
    universe = pd.concat([a[['source', 'target']], b[['source', 'target']]]).drop_duplicates()
    ov = top_k_overlap({"a": a, "b": b}, k=3, universe=universe)
    assert ov["jaccard"]["a"]["b"] == ov["jaccard"]["b"]["a"]
    assert ov["jaccard"]["a"]["a"] == 1.0
    # |{2,3}| / |{1,2,3,9}| = 0.5
    assert ov["jaccard"]["a"]["b"] == pytest.approx(0.5)
    assert ov["top_k_sizes"] == {"a": 3, "b": 3}


def test_heatmap_cells_are_coloured_by_agreement():
    html = render_html(_report())
    assert re.search(r"background:rgba\(0,114,178,[\d.]+\)", html)


# --------------------------------------------------------------------------- #
# regressions: the HTML is the artefact people forward, so what the CLI warns
# about on stderr must also survive into the page
# --------------------------------------------------------------------------- #

def test_no_overlap_verdict_reaches_the_page_before_the_numbers():
    """A vacuous run must announce itself above the verdict table, not below."""
    report = _report()
    report["identifiers"] = {
        "alpha": {"verdict": "no-overlap", "matched_pairs_case_sensitive": 0,
                  "matched_pairs_case_insensitive": 0,
                  "shared_genes_case_insensitive": 0,
                  "example_edge_genes": ["Pou5f1"], "example_gold_genes": ["ENSG1"]},
        "beta": {"verdict": "ok", "matched_pairs_case_sensitive": 2,
                 "matched_pairs_case_insensitive": 2,
                 "shared_genes_case_insensitive": 4},
    }
    html = render_html(report)
    assert "Gene-identifier matching" in html
    assert "NO-OVERLAP" in html
    assert "vacuous" in html
    # the alarm precedes the verdict table it invalidates
    assert html.index("vacuous") < html.index("<h2>Verdict</h2>")


def test_case_mismatch_is_reported_even_though_it_was_recovered():
    report = _report()
    report["identifiers"] = {
        "alpha": {"verdict": "case-mismatch", "matched_pairs_case_sensitive": 0,
                  "matched_pairs_case_insensitive": 61491,
                  "shared_genes_case_insensitive": 1516,
                  "example_edge_genes": ["Pou5f1"], "example_gold_genes": ["POU5F1"]},
    }
    html = render_html(report)
    assert "CASE-MISMATCH" in html
    assert "61,491" in html
    # no page-top alarm: the hits were recovered, so the numbers do mean something
    assert "vacuous" not in html


def test_both_environments_are_rendered_and_labelled():
    """Edge lists get produced on a cluster and scored on a laptop weeks later."""
    report = _report()
    report["environment"] = {"python": "3.11-laptop", "arboreto": "not installed"}
    prep = {"environment": {"python": "3.11-cluster", "arboreto": "0.1.6"}}
    html = render_html(report, prep)
    assert "3.11-laptop" in html and "3.11-cluster" in html
    assert "evaluation" in html and "preprocessing" in html


def test_a_renamed_method_says_so_even_with_no_identifier_block():
    report = _report()
    report["name_collisions"] = [{"file": "runs/b/alpha.edges.tsv.gz",
                                  "derived": "alpha", "renamed_to": "b/alpha"}]
    html = render_html(report)
    assert "b/alpha" in html and "renamed" in html.lower()


def test_undefined_null_is_never_badged_as_significant():
    report = _report()
    report["methods"]["alpha"]["null"] = {
        "metric": "auprc", "observed": float("nan"), "observed_undefined": True,
        "null_mean": 0.48, "null_std": 0.1, "p_value": float("nan"),
        "n_nulls_used": 45, "n_nulls_degenerate": 0}
    html = render_html(report)
    assert '<span class="badge pass">nan' not in html
    assert '<span class="badge fail">n/a</span>' in html


def test_tie_interval_is_shown_next_to_a_coin_flip_epr():
    """A diagnostic that is computed and not rendered is the defect we just fixed."""
    report = _report()
    report["methods"]["alpha"].update(
        {"epr": 4.0, "tie_block_at_k": 100, "epr_range": [0.0, 10.0]})
    html = render_html(report)
    assert "0.00&ndash;10.00" in html
    assert "100 tied" in html


def test_no_tie_interval_when_the_scores_are_distinct():
    report = _report()
    report["methods"]["alpha"].update(
        {"tie_block_at_k": 0, "epr_range": [2.5, 2.5]})
    html = render_html(report)
    assert "tied</span>" not in html


def test_tied_epr_shows_its_interval_not_just_a_point():
    """A coin flip must not be printed as a measurement."""
    report = _report()
    report["methods"]["alpha"].update(epr=4.0, epr_range=[0.0, 10.0],
                                      tie_block_at_k=100)
    html = render_html(report)
    assert "0.00&ndash;10.00" in html and "100 tied" in html


def test_untied_epr_shows_no_interval():
    report = _report()
    report["methods"]["alpha"].update(epr=2.5, epr_range=[2.5, 2.5],
                                      tie_block_at_k=0)
    html = render_html(report)
    assert "tied</span>" not in html


def test_legacy_significance_section_omits_invalid_tiers():
    report = _report()
    report["significance"] = {
        "metric": "epr", "n_boot": 200, "resample": "gold subsample 0.5",
        "methods": ["alpha", "beta"],
        "observed": {"alpha": 3.0, "beta": 2.9},
        "pairs": {"alpha|beta": {"p_a_gt_b": 0.61, "diff_mean": 0.1,
                                 "diff_ci": [-0.4, 0.6], "significant": False}},
        "n_comparisons": 1, "n_significant": 0, "fraction_resolved": 0.0,
        "tiers": [["alpha", "beta"]],
    }
    html = render_html(report)
    assert "Gold-thinning sensitivity" in html
    assert 'Legacy annotation-resampling inference omitted' in html
    assert 'Tiers the data can defend' not in html
    assert '<span class="badge fail">0%</span>' not in html


def test_significance_absent_explains_how_to_get_it():
    html = render_html(_report())
    assert "Gold-thinning sensitivity" in html
    assert "--gold-thinning-sensitivity" in html


def test_thinning_section_uses_scenario_quantiles_without_tiers_or_p_badges():
    from grnbench.eval.significance import gold_thinning_sensitivity
    a = _edges([('TF', 'G1', 3), ('TF', 'G2', 2), ('TF', 'G3', 1)])
    report = _report()
    report['gold_thinning_sensitivity'] = gold_thinning_sensitivity(
        {'a': a, 'b': a.assign(weight=-a.weight)}, {('TF', 'G1')}, n_draws=4)
    html = render_html(report)
    assert 'scenario quantiles' in html and 'fixed k=1' in html
    assert 'not confidence intervals' in html
    assert 'Tiers the data can defend' not in html
    assert 'reporting noise' not in html
    assert 'no better than chance' not in html
