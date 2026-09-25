import hashlib
import json
import importlib.util
from pathlib import Path

import pandas as pd
import numpy as np
import anndata as ad
import pytest

from grnbench.report import file_fingerprint


def test_full_digest_detects_same_size_tail_change(tmp_path):
    path = tmp_path / 'large.bin'
    data = b'a' * 70000
    path.write_bytes(data)
    first = file_fingerprint(path, full_below_mb=0)
    assert first.get('sha256') == hashlib.sha256(data).hexdigest()
    path.write_bytes(data[:-1] + b'b')
    assert file_fingerprint(path, full_below_mb=0)['sha256'] != first['sha256']


def test_recursive_directory_inventory(tmp_path):
    (tmp_path / 'nested').mkdir()
    (tmp_path / 'nested' / 'counts.mtx').write_text('123')
    before = file_fingerprint(tmp_path)
    assert before['files'][0]['path'] == 'nested/counts.mtx'
    (tmp_path / 'nested' / 'counts.mtx').write_text('124')
    assert file_fingerprint(tmp_path)['sha256'] != before['sha256']


def backtest_module():
    spec = importlib.util.spec_from_file_location('backtest_reliable', Path(__file__).parents[1] / 'benchmarks/beeline_backtest.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_beeline_subset_filters_then_selects_non_tfs_by_variance(tmp_path):
    genes = ['TF', 'badTF', 'low', 'high', 'bad']
    expr = pd.DataFrame([[1]*5, [2]*5], columns=genes)
    pd.DataFrame({'pValue':[.0001, .1, .0002, .0003, .5],
                  'variance':[1, 900, 2, 20, 1000]}, index=genes).to_csv(tmp_path / 'GeneOrdering.csv')
    selected, tfs, meta = backtest_module().beeline_subset(expr, tmp_path, ['TF','badTF'], 1)
    assert list(selected.columns) == ['TF', 'high']
    assert tfs == ['TF']
    assert meta['p_value_cutoff_effective'] == .01 / 5


def test_infer_retains_partial_failure_and_complete_universe(tmp_path, monkeypatch):
    from grnbench import cli_infer
    path = tmp_path / 'input.h5ad'
    ad.AnnData(np.array([[1.,2.,3.],[2.,4.,2.],[3.,1.,4.]]),
               var=pd.DataFrame(index=['TF','a','b'])).write_h5ad(path)
    tfpath = tmp_path / 'tfs.txt'
    tfpath.write_text('TF\n')
    monkeypatch.setattr(cli_infer, 'list_methods', lambda:[{'name':n,'available':True} for n in ['good','bad']])
    def run(name, *args, **kwargs):
        if name == 'bad':
            raise RuntimeError('deliberate method failure')
        return pd.DataFrame({'source':['TF'],'target':['a'],'weight':[.9]}), {'method':name,'custom_evidence':'kept'}
    monkeypatch.setattr(cli_infer,'run_method',run)
    out = tmp_path / 'out'
    rc = cli_infer.main([str(path),'--out',str(out),'--methods','good,bad','--tf-list',str(tfpath)])
    assert rc != 0
    assert (out / 'inference_manifest.json').exists()
    manifest = json.loads((out / 'inference_manifest.json').read_text())
    assert manifest['runs']['good']['custom_evidence'] == 'kept'
    assert 'deliberate method failure' in manifest['errors']['bad']
    universe = pd.read_csv(out/'candidate_universe.tsv.gz',sep='\t')
    assert set(map(tuple,universe.to_numpy())) == {('TF','a'),('TF','b')}


def test_backtest_unavailable_method_writes_diagnostic_manifest(tmp_path, monkeypatch):
    import grnbench.infer
    monkeypatch.setattr(grnbench.infer,'list_methods',lambda:[])
    out = tmp_path/'out'
    rc = backtest_module().main(['--data',str(tmp_path),'--networks',str(tmp_path),
                               '--out',str(out),'--methods','missing'])
    assert rc != 0
    assert (out/'beeline_backtest.json').exists()
    assert json.loads((out/'beeline_backtest.json').read_text())['status'] == 'failed'


def test_subset_agrees_with_pinned_official_fixture(tmp_path):
    import subprocess
    import sys
    official = Path(__file__).parents[1]/'audit/2026-09-13/implementation/beeline-reference-generateExpInputs.py'
    assert hashlib.sha256(official.read_bytes()).hexdigest() == backtest_module().PROTOCOL['reference_sha256']
    genes = ['TF', 'badTF', 'low', 'high', 'bad']
    expr = pd.DataFrame([[1]*5,[2]*5],columns=genes)
    expr.T.to_csv(tmp_path/'ExpressionData.csv')
    pd.DataFrame({'p':[.0001,.1,.0002,.0003,.5],'variance':[1,900,2,20,1000]},index=genes).to_csv(tmp_path/'GeneOrdering.csv')
    pd.DataFrame({'TF':['TF','badTF']}).to_csv(tmp_path/'tfs.csv',index=False)
    run = subprocess.run([sys.executable,str(official),'-e',str(tmp_path/'ExpressionData.csv'),
                          '-g',str(tmp_path/'GeneOrdering.csv'),'-i',str(tmp_path/'tfs.csv'),
                          '-n','1','-o',str(tmp_path/'reference')],capture_output=True,text=True)
    assert run.returncode == 0, run.stderr
    reference = pd.read_csv(tmp_path/'reference-ExpressionData.csv',index_col=0).T
    actual,_,_ = backtest_module().beeline_subset(expr,tmp_path,['TF','badTF'],1)
    assert list(actual.columns) == list(reference.columns)


def test_output_reuse_refused_without_overwriting(tmp_path):
    from grnbench.report import require_empty_output
    import pytest
    out = require_empty_output(tmp_path/'out')
    artifact = out/'old.edges.tsv'
    artifact.write_text('previous run')
    with pytest.raises(ValueError,match='empty'):
        require_empty_output(out)
    assert artifact.read_text() == 'previous run'


def test_backtest_persists_replayable_edges_and_provenance(tmp_path, monkeypatch):
    import grnbench.infer
    ds = tmp_path/'data'/'hESC'
    ds.mkdir(parents=True)
    expr = pd.DataFrame([[1,2,3],[3,1,5],[2,5,2],[7,3,1]],columns=['TF','A','B'])
    expr.T.to_csv(ds/'ExpressionData.csv')
    pd.DataFrame({'p':[0,0,0],'var':[1,2,3]},index=expr.columns).to_csv(ds/'GeneOrdering.csv')
    net = tmp_path/'networks'
    (net/'Networks'/'human').mkdir(parents=True)
    pd.DataFrame({'TF':['TF']}).to_csv(net/'human-tfs.csv',index=False)
    for name in ['hESC-ChIP-seq-network.csv','Non-specific-ChIP-seq-network.csv','STRING-network.csv']:
        pd.DataFrame({'source':['TF'],'target':['A']}).to_csv(net/'Networks'/'human'/name,index=False)
    out = tmp_path/'out'
    rc = backtest_module().main(['--data',str(tmp_path/'data'),'--networks',str(net),
                                '--out',str(out),'--methods','spearman_signed',
                                '--n-nulls','2','--n-var','2','--data-release','test-fixture-v1'])
    result = json.loads((out/'beeline_backtest.json').read_text())
    assert rc == 0, result['errors']
    assert result['data_release'] == 'test-fixture-v1'
    assert len(result['rows']) == 3
    prov = result['datasets']['hESC']['runs']['spearman_signed']
    assert Path(prov['edges_file']).exists()
    assert prov['edges_fingerprint']['sha256']
    assert result['datasets']['hESC']['effective_genes']['values'] == ['A','B','TF']
    from sklearn.metrics import average_precision_score, roc_auc_score
    persisted = pd.read_csv(prov['edges_file'],sep='\t')
    candidates = pd.read_csv(out/'hESC'/'candidate_universe.tsv.gz',sep='\t')
    aligned = candidates.merge(persisted,on=['source','target'],how='left')
    assert aligned.weight.notna().all()
    labels = (aligned.source.eq('TF') & aligned.target.eq('A')).astype(int)
    for row in result['rows']:
        assert row['auprc'] == pytest.approx(average_precision_score(labels,aligned.weight))
        assert row['auroc'] == pytest.approx(roc_auc_score(labels,aligned.weight))


@pytest.mark.parametrize('fail_stage',['pip','pytest'])
def test_hpc_setup_propagates_install_and_test_failures(tmp_path,fail_stage):
    import os
    import shutil
    import subprocess
    bash = shutil.which('bash') or 'C:/Program Files/Git/bin/bash.exe'
    if not Path(bash).exists():
        pytest.skip('Bash unavailable; HPC failure-path test needs Bash')
    profile = tmp_path/'conda.sh'
    profile.write_text('''conda() { if [ "$1" = "env" ]; then echo test-env; fi; return 0; }
grnbench() { return 0; }
python() { if [ "$2" = "$FAIL_STAGE" ]; then return 23; fi; return 0; }
''')
    env = {**os.environ,'CONDA_SH':profile.as_posix(),'FAIL_STAGE':fail_stage}
    script = Path(__file__).parents[1]/'hpc/setup_env.sh'
    result = subprocess.run([bash,str(script),'test-env',str(Path(__file__).parents[1])],
                            env=env,capture_output=True,text=True)
    assert result.returncode == 23, result.stdout + result.stderr
    assert 'SETUP_DONE' not in result.stdout


def _synthetic_rows(n_datasets=7, methods=('a', 'b', 'c'), seed=0):
    """Datasets x methods x references; method 'a' beats 'b' on EPR by a clear paired margin."""
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n_datasets):
        for gt in ['specific-chip', 'string']:
            base = rng.normal(1.0, .05)
            for m in methods:
                shift = {'a': .5, 'b': 0., 'c': rng.normal(0, .6)}[m]
                rows.append({'dataset': f'd{i}', 'method': m, 'ground_truth': gt,
                             'epr': base + shift + rng.normal(0, .02),
                             'auprc_ratio': base + rng.normal(0, .02)})
    return rows


def test_unit_comparison_uses_datasets_as_units_and_resolves_only_stable_gaps():
    mod = backtest_module()
    res = mod.unit_comparison(_synthetic_rows(), alpha=.05)
    assert res['units'] == 'datasets' and res['unresolved_is_equivalent'] is False
    blocks = {(b['ground_truth'], b['metric']): b for b in res['blocks']}
    assert all(b['n_complete_datasets'] == 7 and b['status'] == 'estimated' for b in blocks.values())
    epr = blocks[('specific-chip', 'epr')]
    assert {(r['winner'], r['loser']) for r in epr['resolved']} >= {('a', 'b')}
    # auprc_ratio carries no real method effect: nothing should resolve there
    assert blocks[('specific-chip', 'auprc_ratio')]['resolved'] == []
    assert blocks[('string', 'auprc_ratio')]['resolved'] == []


def test_unit_comparison_refuses_inference_below_three_units():
    mod = backtest_module()
    res = mod.unit_comparison(_synthetic_rows(n_datasets=2))
    assert res['n_pairs'] == 0 and res['n_resolved'] == 0
    assert {b['status'] for b in res['blocks']} == {'too_few_units'}


def test_unit_comparison_drops_incomplete_datasets_loudly():
    mod = backtest_module()
    rows = [r for r in _synthetic_rows() if not (r['dataset'] == 'd3' and r['method'] == 'c')]
    block = next(b for b in mod.unit_comparison(rows)['blocks']
                 if b['ground_truth'] == 'specific-chip' and b['metric'] == 'epr')
    assert block['n_datasets'] == 7 and block['n_complete_datasets'] == 6
    assert block['dropped_datasets'] == ['d3']


def test_regime_check_only_applies_to_specific_chip_and_flags_far_outliers():
    mod = backtest_module()
    rows = [{'dataset': 'x', 'method': 'm', 'ground_truth': 'specific-chip', 'epr': 1.5, 'auprc_ratio': 4.0},
            {'dataset': 'x', 'method': 'm', 'ground_truth': 'string', 'epr': 9.0, 'auprc_ratio': 9.0},
            {'dataset': 'x', 'method': 'n', 'ground_truth': 'specific-chip', 'epr': float('nan'), 'auprc_ratio': 1.0}]
    rc = mod.regime_check(rows)
    assert rc['applies_to'] == ['specific-chip']
    assert rc['n_checked'] == 3  # string rows ignored, NaN skipped
    flagged = [(c['metric'], c['above_regime']) for c in rc['checks'] if c['method'] == 'm']
    assert dict(flagged) == {'epr': False, 'auprc_ratio': True}
    assert rc['n_above_regime'] == 1


def test_summarize_flag_annotates_existing_json_and_writes_markdown(tmp_path):
    mod = backtest_module()
    src = tmp_path / 'beeline_backtest.json'
    src.write_text(json.dumps({'status': 'complete', 'rows': _synthetic_rows(), 'errors': []}), encoding='utf-8')
    rc = mod.main(['--data', 'x', '--networks', 'x', '--out', 'x', '--summarize', str(src)])
    assert rc == 0
    result = json.loads(src.read_text(encoding='utf-8'))
    assert result['unit_comparison']['n_pairs'] == 4 * 3
    md = (tmp_path / 'summary.md').read_text(encoding='utf-8')
    assert '## Dataset-unit comparison' in md and 'unresolved does not mean equivalent' in md
    assert '### specific-chip / epr - n=7 datasets' in md


def test_regime_check_labels_small_count_epr_excursions():
    """8 hits against 1.8 expected on a 680-pair reference is EPR 4.5 - noise, not a leak."""
    mod = backtest_module()
    rows = [{'dataset': 'mDC', 'method': 'clr', 'ground_truth': 'specific-chip',
             'epr': 4.554, 'auprc_ratio': 1.15, 'epr_k': 680, 'prevalence': 680 / 263220,
             'precision_at_k': 8 / 680}]
    rc = mod.regime_check(rows)
    epr = next(c for c in rc['checks'] if c['metric'] == 'epr')
    assert epr['above_regime'] and epr['small_count']
    assert epr['hits_at_k'] == pytest.approx(8) and epr['expected_random_hits'] == pytest.approx(1.76, abs=.01)
    assert rc['n_above_regime'] == 1 and rc['n_above_regime_small_count'] == 1
