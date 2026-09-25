"""Regression invariants from the dated scientific audit."""
import numpy as np
import pandas as pd
import pytest

from grnbench.eval.metrics import GoldIndex, evaluate_ranking, top_k_overlap
from grnbench.eval.nulls import random_null, degree_preserving_nulls


def frame(n=100):
    return pd.DataFrame({'source': ['TF'] * n,
                         'target': [f'G{i}' for i in range(n)],
                         'weight': np.arange(n, 0, -1, dtype=float)})


def test_truncation_does_not_remove_gold_from_denominator():
    e = frame()
    gold = {('TF', 'G0'), ('TF', 'G99')}
    result = evaluate_ranking(e.head(1), gold, universe=e[['source', 'target']])
    assert result['auprc'] == pytest.approx(.51)
    assert result['n_gold_in_universe'] == 2
    assert result['missing_predictions'] == 99
    assert result['support_origin'] == 'explicit-universe'


@pytest.mark.parametrize('bad', [np.nan, np.inf, -np.inf])
def test_invalid_scores_rejected(bad):
    with pytest.raises(ValueError, match='finite'):
        evaluate_ranking(frame(3).assign(weight=[bad, 2, 1]), {('TF', 'G0')})


def test_duplicates_and_case_collisions_rejected():
    with pytest.raises(ValueError, match='duplicate'):
        evaluate_ranking(pd.concat([frame(3), frame(3)]), {('TF', 'G0')})
    with pytest.raises(ValueError, match='collision'):
        GoldIndex([('TF', 'Gene'), ('TF', 'GENE')])


def test_random_null_keeps_candidate_pairs_exactly():
    e = frame(12).assign(source=[f'TF{i % 3}' for i in range(12)])
    for null in random_null(e, n=5):
        assert set(zip(null.source, null.target)) == set(zip(e.source, e.target))
        assert sorted(null.weight) == sorted(e.weight)


def test_dense_degree_rewire_reports_degeneracy():
    e = pd.DataFrame([('A', 'x', 4), ('A', 'y', 3), ('B', 'x', 2), ('B', 'y', 1)],
                     columns=['source', 'target', 'weight'])
    with pytest.raises(ValueError, match='degenerate'):
        degree_preserving_nulls(e, n=3, seed=0)


def test_overlap_includes_entire_boundary_tie():
    a = frame(8).assign(weight=1.0)
    b = a.iloc[::-1].reset_index(drop=True)
    r = top_k_overlap({'a': a, 'b': b}, k=2)
    assert r['jaccard']['a']['b'] == 1
    assert r['top_k_sizes']['a'] == 8


def test_thinning_freezes_cutoff_and_makes_no_inferential_claim():
    from grnbench.eval.significance import gold_thinning_sensitivity
    # Audit counterexample: full k=100 favours A; recomputing k=50 reverses it.
    e = frame(1000)
    gold = {('TF', f'G{i}') for i in range(100)}
    def ranking(positions):
        positive, negative = iter(range(100)), iter(range(100, 1000))
        order = [next(positive) if i in positions else next(negative) for i in range(1000)]
        return e.iloc[order].assign(weight=np.arange(1000, 0, -1))
    a = ranking(set(range(50, 90)) | set(range(200, 260)))
    b = ranking(set(range(25)) | set(range(300, 375)))
    r = gold_thinning_sensitivity({'a': a, 'b': b}, gold, n_draws=100, seed=1)
    assert r['epr_k'] == 100
    assert r['observed'] == pytest.approx({'a': 4., 'b': 2.5})
    assert r['pairs']['a|b']['diff_quantiles'][0] > 0
    assert 'tiers' not in r and 'n_significant' not in r
    assert 'significant' not in r['pairs']['a|b']


def test_comparison_rejects_different_support_without_universe():
    from grnbench.eval.significance import gold_thinning_sensitivity
    with pytest.raises(ValueError, match='identical.*support'):
        gold_thinning_sensitivity({'a': frame(20), 'b': frame(19)}, {('TF', 'G0')})


def test_repeated_symbols_are_valid_but_diagonal_is_not():
    assert evaluate_ranking(frame(10), {('TF', 'G0')})['n_edges'] == 10
    with pytest.raises(ValueError, match='self-loops'):
        evaluate_ranking(frame(1).assign(target='TF'), {('TF', 'G0')})


def test_empty_output_can_be_scored_only_against_explicit_universe():
    r = evaluate_ranking(frame(0), {('TF', 'G0'), ('TF', 'G9')}, universe=frame(10))
    assert r['auprc'] == pytest.approx(.2)
    assert r['auroc'] == pytest.approx(.5)
    assert r['missing_predictions'] == 10


def test_sparse_degree_swaps_preserve_unique_in_and_out_degrees_and_no_loops():
    e = pd.DataFrame([('A', 'B', 6), ('B', 'C', 5), ('C', 'D', 4),
                      ('D', 'E', 3), ('E', 'F', 2), ('F', 'A', 1)],
                     columns=['source', 'target', 'weight'])
    for null in degree_preserving_nulls(e, n=8):
        assert not null.duplicated(['source', 'target']).any()
        assert not (null.source == null.target).any()
        assert null.groupby('source').target.nunique().to_dict() == e.groupby('source').target.nunique().to_dict()
        assert null.groupby('target').source.nunique().to_dict() == e.groupby('target').source.nunique().to_dict()


def test_regulator_null_is_a_separate_hypothesis_with_fixed_source_score_multisets():
    from grnbench.eval.nulls import regulator_stratified_nulls, null_metrics
    e = frame(12).assign(source=[f'TF{i % 3}' for i in range(12)])
    for null in regulator_stratified_nulls(e, n=4):
        for source in e.source.unique():
            assert sorted(null.loc[null.source == source, 'weight']) == sorted(e.loc[e.source == source, 'weight'])
    r = null_metrics(e, {('TF0', 'G0')}, mode='regulator', n_nulls=4)
    assert 'within each source' in r['null_hypothesis']


def test_degree_diagnostic_never_produces_inferential_p_value():
    from grnbench.eval.nulls import null_metrics
    e = frame(12).assign(source=[f'TF{i % 3}' for i in range(12)])
    r = null_metrics(e, {('TF0', 'G0')}, mode='degree', n_nulls=4)
    assert np.isnan(r['p_value']) and not r['inference_available']


def test_cli_sparse_metadata_and_gold_reachability_use_explicit_universe(tmp_path):
    import json
    from grnbench.cli_eval import main
    u = frame(10)
    u[['source', 'target']].to_csv(tmp_path / 'u.tsv', sep='\t', index=False)
    u.head(3).to_csv(tmp_path / 'a.tsv', sep='\t', index=False)
    (tmp_path / 'gold.tsv').write_text('source\ttarget\nTF\tG9\n')
    assert main([str(tmp_path / 'a.tsv'), '--gold', str(tmp_path / 'gold.tsv'),
                 '--universe', str(tmp_path / 'u.tsv'), '--out', str(tmp_path / 'out'),
                 '--skip-nulls']) == 0
    r = json.loads((tmp_path / 'out/eval_report.json').read_text())
    assert r['methods']['a']['missing_predictions'] == 7
    assert r['methods']['a']['support_origin'] == 'explicit-universe'
    assert r['coverage']['a']['gold_reachable'] == 1


def test_algorithm_seed_is_independent_from_cell_sampling_seed(monkeypatch):
    import anndata as ad
    import grnbench.eval.stability as stability
    calls = []
    def fake_run(method, data, **params):
        calls.append((tuple(data.obs_names), params['seed']))
        return frame(10), {}
    monkeypatch.setattr(stability, 'run_method', fake_run)
    data = ad.AnnData(np.ones((20, 3)))
    stability.subsample_stability(data, 'fake', seed=7, algorithm_seed=99,
                                  n_subsamples=3, frac=.5, top_k=3)
    assert {seed for _, seed in calls} == {99}
    assert len({idx for idx, _ in calls}) > 1
    calls.clear()
    stability.subsample_stability(data, 'fake', seed=7, algorithm_seed=99,
                                  vary_algorithm_seed=True, n_subsamples=3, frac=.5, top_k=3)
    assert [seed for _, seed in calls] == [99, 100, 101]
