"""Annotation thinning is descriptive; no equivalence or superiority claims."""
import numpy as np
import pandas as pd
import pytest
from grnbench.eval.significance import (gold_thinning_sensitivity,
    reference_ordering_differences, render_significance_md, render_contradiction_md)


def world():
    edges = pd.DataFrame({'source': ['TF'] * 100,
                         'target': [f'G{i}' for i in range(100)],
                         'weight': np.arange(100, 0, -1)})
    return edges, {('TF', f'G{i}') for i in range(20)}


def test_ordering_summary_uses_scenario_frequency_not_probability_of_superiority():
    a, gold = world()
    b = a.assign(weight=-a.weight)
    r = gold_thinning_sensitivity({'a': a, 'b': b}, gold, n_draws=40)
    assert r['pairs']['a|b']['fraction_a_gt_b'] == 1
    assert r['pairs']['a|b']['diff_quantiles'][0] > 0
    assert not {'tiers', 'fraction_resolved', 'n_significant'} & r.keys()


def test_identical_methods_have_zero_difference_without_equivalence_claim():
    a, gold = world()
    r = gold_thinning_sensitivity({'a': a, 'b': a}, gold, n_draws=15)
    rec = r['pairs']['a|b']
    assert rec['diff_quantiles'] == [0, 0]
    assert rec['fraction_a_gt_b'] == 0
    assert 'significant' not in rec


@pytest.mark.parametrize('kwargs', [{'frac': 0}, {'frac': 1.01}, {'frac': np.nan},
                                  {'n_draws': 0}, {'n_draws': -1}, {'n_draws': 2.5},
                                  {'top_k': 0}])
def test_invalid_scenario_configuration_rejected(kwargs):
    a, gold = world()
    with pytest.raises(ValueError):
        gold_thinning_sensitivity({'a': a, 'b': a}, gold, **kwargs)


def test_fraction_one_returns_full_observed_scores():
    a, gold = world()
    b = a.assign(weight=-a.weight)
    r = gold_thinning_sensitivity({'a': a, 'b': b}, gold, frac=1, n_draws=5)
    delta = r['observed']['a'] - r['observed']['b']
    assert r['pairs']['a|b']['diff_quantiles'] == [delta, delta]




def test_references_describe_different_tasks_without_confidently_wrong_claims():
    blocks = [{'dataset': 'cells', 'ground_truth': ref,
               'pairs': {'a|b': {'observed_difference': diff}}}
              for ref, diff in [('binding', 1.), ('perturbation', -1.)]]
    r = reference_ordering_differences(blocks)
    assert r['n_reversed'] == 1 and r['reversal_rate'] == 1
    md = render_contradiction_md(r)
    assert '100%' in md and 'different biological relations' in md
    assert 'cannot both be right' not in md


def test_markdown_labels_quantiles_and_fixed_cutoff():
    a, gold = world()
    result = gold_thinning_sensitivity({'a': a, 'b': a}, gold, n_draws=5)
    md = render_significance_md(result)
    assert 'scenario quantiles' in md and 'fixed k=20' in md
    assert 'not confidence intervals' in md and 'Tiers' not in md


def test_explicit_universe_accepts_sparse_predictions_and_keeps_reachable_gold():
    a, gold = world()
    r = gold_thinning_sensitivity({'a': a, 'b': a.head(1)}, gold,
                                  universe=a[['source', 'target']], n_draws=5)
    assert r['n_gold_reachable'] == 20 and r['epr_k'] == 20


def test_case_sensitive_contract_is_respected():
    from grnbench.eval.metrics import GoldIndex
    a, gold = world()
    a = a.assign(source='tf')
    with pytest.raises(ValueError, match='no reachable'):
        gold_thinning_sensitivity({'a': a, 'b': a}, GoldIndex(gold, case_sensitive=True))
