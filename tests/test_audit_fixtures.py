"""Audit 2026-09-13 fixtures are regressions, not just evidence.

``audit/2026-09-13/eval_findings.json`` records six inputs on which the
pre-0.7.0 evaluator misbehaved. Every recorded case must stay mapped to a test
that feeds the same input to the current code and asserts the repaired
behaviour; a new audit case without a regression fails here.
"""
import importlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from grnbench.eval.metrics import average_precision, evaluate_ranking
from grnbench.eval.nulls import degree_preserving_nulls, null_metrics, random_null

AUDIT = Path(__file__).resolve().parents[1] / 'audit' / '2026-09-13'

REGRESSIONS = {
    'default_cutoff_inversion':
        ('test_eval_contract', 'test_thinning_freezes_cutoff_and_makes_no_inferential_claim'),
    'degree_null_duplicates':
        ('test_eval_contract', 'test_dense_degree_rewire_reports_degeneracy'),
    'null_support_and_random_margins':
        ('test_audit_fixtures', 'test_complete_digraph_null_keeps_loop_free_support'),
    'tier_nontransitivity':
        ('test_report_html', 'test_legacy_significance_section_omits_invalid_tiers'),
    'candidate_list_truncation':
        ('test_eval_contract', 'test_truncation_does_not_remove_gold_from_denominator'),
    'finite_random_ap_baseline':
        ('test_audit_fixtures', 'test_prevalence_baseline_is_disclosed_as_not_exact'),
}


def findings():
    data = json.loads((AUDIT / 'eval_findings.json').read_text(encoding='utf-8'))
    return {f['case']: f for f in data['findings']}


def test_every_audit_finding_has_a_named_regression():
    assert set(findings()) == set(REGRESSIONS)
    for module, name in REGRESSIONS.values():
        assert callable(getattr(importlib.import_module(module), name, None)), (module, name)


def test_complete_digraph_null_keeps_loop_free_support():
    f = findings()['null_support_and_random_margins']
    assert f['n_null_self_loops'] > 0      # what the audited code produced
    e = pd.DataFrame(f['original'])[['source', 'target', 'weight']]
    # Every loop-free pair is already present: no simple swap is admissible.
    with pytest.raises(ValueError, match='degenerate'):
        degree_preserving_nulls(e, n=1, seed=0)
    r = null_metrics(e, {('A', 'B')}, mode='degree', n_nulls=5)
    assert r['n_nulls_degenerate'] == 5 and not r['inference_available']
    assert np.isnan(r['p_value'])
    for null in random_null(e, n=5, seed=0):
        assert set(zip(null.source, null.target)) == set(zip(e.source, e.target))


def test_prevalence_baseline_is_disclosed_as_not_exact():
    f = findings()['finite_random_ap_baseline']
    n = f['n_edges']
    exact = np.mean([average_precision(np.eye(n, dtype=int)[i]) for i in range(n)])
    assert exact == pytest.approx(f['exact_expected_random_AP'])
    assert not np.isclose(exact, f['prevalence'])
    e = pd.DataFrame({'source': ['T'] * n, 'target': [f'G{i}' for i in range(n)],
                      'weight': np.arange(n, 0, -1, dtype=float)})
    r = evaluate_ranking(e, {('T', 'G2')}, with_curve=False)
    assert r['prevalence'] == pytest.approx(f['prevalence'])
    assert 'not exact finite-sample random AP' in r['ap_baseline_note']
