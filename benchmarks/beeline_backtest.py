"""Replayable BEELINE-derived selection and grnbench fixed-universe metrics.

This is not historical leaderboard replication or biological validation.
"""
from __future__ import annotations
import argparse
import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd

PROTOCOL = {
    'name': 'beeline-v1.1-selection-grnbench-fixed-universe-v1',
    'reference_commit': '37464085eb8a95d6cc6a3d3a3c649d36db6052ed',
    'reference_url': 'https://github.com/Murali-group/Beeline/blob/37464085eb8a95d6cc6a3d3a3c649d36db6052ed/utils/generateExpInputs.py',
    'reference_sha256': 'a50fbf3ff9e082f2f07f4f123f4d287ef6864be28ab695392aa54f701fc007c2',
    'selection': 'p<0.01/n tested expression genes; significant TFs plus top N significant non-TFs by variance; alphabetical output',
    'identifier_extension': 'Unambiguous case-insensitive TF matching preserving expression identifiers',
    'evaluation': 'TF x selected genes, no self loops; absent outputs tied last; average precision',
    'scope': 'BEELINE-derived selection and grnbench metrics; not historical leaderboard replication or biological validation',
}
SPECIFIC_NETWORK = {
    'hESC': ('human', 'hESC-ChIP-seq-network.csv'),
    'hHEP': ('human', 'HepG2-ChIP-seq-network.csv'),
    'mDC': ('mouse', 'mDC-ChIP-seq-network.csv'),
    'mESC': ('mouse', 'mESC-ChIP-seq-network.csv'),
    'mHSC-E': ('mouse', 'mHSC-ChIP-seq-network.csv'),
    'mHSC-GM': ('mouse', 'mHSC-ChIP-seq-network.csv'),
    'mHSC-L': ('mouse', 'mHSC-ChIP-seq-network.csv'),
}
GENERIC_NETWORKS = {
    'human': [('nonspecific-chip','Non-specific-ChIP-seq-network.csv'),('string','STRING-network.csv')],
    'mouse': [('nonspecific-chip','Non-Specific-ChIP-seq-network.csv'),('string','STRING-network.csv')],
}
EXTRA_NETWORKS = {'mESC': [('lofgof', 'mouse/mESC-lofgof-network.csv')]}
# Best real-data scores reported by Pratapa et al. 2020 (Nat Methods 17:147),
# against cell-type-specific ChIP-seq ONLY. No comparable published figure
# exists for the other references, so they are not compared to anything.
PUBLISHED_REGIME = {'specific-chip': {'auprc_ratio': 1.10, 'epr': 1.20}}
REGIME_FACTOR = 3.0
SMALL_COUNT_HITS = 5.0   # expected random hits at k below this: EPR is a small count
INFERENCE_METRICS = ('epr', 'auprc_ratio')


def build_parser():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['data','networks','out']:
        p.add_argument('--'+key, required=True)
    p.add_argument('--datasets', default=None)
    p.add_argument('--methods', default='spearman_signed,genie3,grnboost2')
    p.add_argument('--plugin', action='append', default=None)
    p.add_argument('--n-var', type=int, default=500, help='Number of significant non-TF genes; add significant TFs')
    p.add_argument('--n-nulls', type=int, default=25)
    p.add_argument('--n-boot', type=int, default=200, help='Deprecated compatibility argument; no inferential label bootstrap is run')
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--data-release', default='unspecified-local-input', help='Release/DOI declaration; zenodo:19009603:v4 for fetch_beeline.sh')
    p.add_argument('--alpha', type=float, default=.05, help='Family-wise level for the dataset-unit comparison')
    p.add_argument('--summarize', default=None, metavar='JSON',
                   help='Recompute regime check and dataset-unit comparison on an existing beeline_backtest.json and exit')
    return p


def regime_check(rows) -> dict:
    """Descriptive: where do specific-chip scores sit relative to the published best?"""
    checks = []
    for r in rows:
        ceiling = PUBLISHED_REGIME.get(r['ground_truth'])
        if not ceiling:
            continue
        for metric, best in ceiling.items():
            v = r.get(metric)
            if v is None or not np.isfinite(v):
                continue
            check = {'dataset': r['dataset'], 'method': r['method'], 'metric': metric,
                     'value': float(v), 'published_best': best,
                     'ratio_to_published': float(v) / best,
                     'above_regime': bool(v > REGIME_FACTOR * best)}
            k, prev, p_k = r.get('epr_k'), r.get('prevalence'), r.get('precision_at_k')
            if metric == 'epr' and None not in (k, prev, p_k):
                # EPR is a ratio of two counts; on a small reference it is a
                # handful of hits against an expectation of one or two.
                check['hits_at_k'] = float(p_k) * k
                check['expected_random_hits'] = float(prev) * k
                check['small_count'] = bool(prev * k < SMALL_COUNT_HITS)
            checks.append(check)
    above = [c for c in checks if c['above_regime']]
    return {'reference': 'Pratapa et al. 2020, real scRNA-seq vs cell-type-specific ChIP-seq',
            'applies_to': list(PUBLISHED_REGIME), 'factor': REGIME_FACTOR,
            'n_checked': len(checks), 'n_above_regime': len(above),
            'n_above_regime_small_count': sum(bool(c.get('small_count')) for c in above),
            'checks': checks,
            'interpretation': ('Sanity check on the evaluation code, not a leaderboard. A score far above the '
                               'published best on the same reference is either a scoring/matching defect or a '
                               'small-count reference (hits_at_k against expected_random_hits below '
                               f'{SMALL_COUNT_HITS:g}); read the counts before quoting the ratio.')}


def unit_comparison(rows, alpha=.05, metrics=INFERENCE_METRICS) -> dict:
    """Paired t + Holm with datasets as the independent units, one family per reference and metric."""
    from grnbench.eval.comparison import compare_units
    table = pd.DataFrame(rows)
    blocks, resolved, n_pairs = [], 0, 0
    if table.empty:
        return {'units': 'datasets', 'alpha': alpha, 'blocks': blocks, 'n_resolved': 0, 'n_pairs': 0}
    for gt, sub in table.groupby('ground_truth', sort=True):
        for metric in metrics:
            wide = sub.pivot_table(index='dataset', columns='method', values=metric, aggfunc='first')
            complete = wide.dropna(axis=0, how='any')
            block = {'ground_truth': gt, 'metric': metric, 'n_datasets': int(wide.shape[0]),
                     'n_complete_datasets': int(complete.shape[0]),
                     'dropped_datasets': sorted(set(wide.index) - set(complete.index))}
            if complete.shape[0] < 3 or complete.shape[1] < 2:
                block['status'] = 'too_few_units'
            else:
                res = compare_units(complete, alpha=alpha)
                block.update(status='estimated', comparison=res)
                won = [p for p in res['pairs'] if p['winner']]
                resolved += len(won); n_pairs += res['n_pairs']
                block['resolved'] = [{'winner': p['winner'],
                                      'loser': p['b'] if p['winner'] == p['a'] else p['a'],
                                      'difference': p['difference']} for p in won]
            blocks.append(block)
    return {'units': 'datasets', 'alpha': alpha, 'metrics': list(metrics),
            'procedure': 'paired_t_holm per reference and metric; Bonferroni simultaneous intervals',
            'n_resolved': resolved, 'n_pairs': n_pairs, 'blocks': blocks,
            'unresolved_is_equivalent': False,
            'interpretation': 'Datasets are the only independent replicates in this design; gold edges and cells are not units. Unresolved pairs are not shown equivalent.'}


def render_summary_md(result) -> str:
    from grnbench.eval.comparison import render_comparison
    rc, uc = result.get('regime_check', {}), result.get('unit_comparison', {})
    best = PUBLISHED_REGIME['specific-chip']
    lines = ['# BEELINE backtest summary', '',
             f"Status: {result.get('status')}; {len(result.get('rows', []))} scored results; {len(result.get('errors', []))} errors.", '',
             '## Published-regime check (specific-chip only)', '',
             f"{rc.get('n_above_regime', 0)} of {rc.get('n_checked', 0)} scores exceed {rc.get('factor')}x the published best "
             f"(AUPRC ratio {best['auprc_ratio']}, EPR {best['epr']}); "
             f"{rc.get('n_above_regime_small_count', 0)} of those are small-count EPRs.", '',
             rc.get('interpretation', ''), '']
    if rc.get('checks'):
        lines += ['| dataset | method | metric | value | ratio to published best | hits at k / expected at random |',
                  '|---|---|---|---:|---:|---:|']
        for c in rc['checks']:
            counts = (f"{c['hits_at_k']:.0f} / {c['expected_random_hits']:.1f}" + (' (small count)' if c.get('small_count') else '')
                      if 'hits_at_k' in c else '')
            flag = ' **above regime**' if c['above_regime'] else ''
            lines.append(f"| {c['dataset']} | {c['method']} | {c['metric']} | {c['value']:.3f} | {c['ratio_to_published']:.2f}{flag} | {counts} |")
        lines.append('')
    lines += ['## Dataset-unit comparison', '',
              f"Units: datasets. {uc.get('n_resolved', 0)} of {uc.get('n_pairs', 0)} method pairs resolved at alpha={uc.get('alpha')}.", '',
              uc.get('interpretation', ''), '']
    assumptions = None
    for b in uc.get('blocks', []):
        lines.append(f"### {b['ground_truth']} / {b['metric']} - n={b['n_complete_datasets']} datasets")
        if b.get('dropped_datasets'):
            lines.append(f"Dropped (incomplete across methods): {', '.join(b['dropped_datasets'])}")
        if b['status'] != 'estimated':
            lines += ['', f"Too few complete units for inference ({b['n_complete_datasets']} < 3).", '']
            continue
        body = render_comparison(b['comparison']).split('\n', 1)[1]
        body, _, assumptions = body.partition('\n## Assumptions\n')
        lines += ['', body.strip('\n'), '']
    if assumptions:
        lines += ['## Assumptions (all comparisons above)', '', assumptions.strip('\n'), '']
    return '\n'.join(lines) + '\n'


def summarize(result: dict, alpha: float) -> dict:
    result['regime_check'] = regime_check(result['rows'])
    result['unit_comparison'] = unit_comparison(result['rows'], alpha=alpha)
    return result


def find_datasets(root: Path) -> dict[str, Path]:
    canon = {k.casefold(): k for k in SPECIFIC_NETWORK}
    found = {}
    for expr in sorted(root.rglob('ExpressionData.csv')):
        d = expr.parent
        if not (d / 'GeneOrdering.csv').exists():
            continue
        for part in [d.name, *[q.name for q in d.parents[:2]]]:
            key = canon.get(part.casefold())
            if key:
                if key in found and found[key] != d:
                    raise ValueError(f'Multiple inputs for {key}: {found[key]}, {d}; narrow --data')
                found[key] = d
                break
    return found


def load_expression(d):
    expr = pd.read_csv(d / 'ExpressionData.csv', index_col=0).T
    expr.columns = expr.columns.astype(str)
    return expr


def beeline_subset(expr, d, tfs, n_var):
    from grnbench.tfs import resolve_tf_names
    if n_var < 0:
        raise ValueError('n_var must be nonnegative')
    order = pd.read_csv(d / 'GeneOrdering.csv', index_col=0)
    order.index = order.index.astype(str)
    if order.index.has_duplicates or expr.columns.has_duplicates:
        raise ValueError('Duplicate gene identifiers are ambiguous')
    if len(order.columns) < 2:
        raise ValueError('GeneOrdering requires p-value and variance columns')
    order = order.loc[order.index.isin(expr.columns)].copy()
    if order.empty:
        raise ValueError('GeneOrdering has no expression genes')
    values = order.iloc[:, :2].apply(pd.to_numeric, errors='raise')
    if not np.isfinite(values.to_numpy()).all() or not values.iloc[:,0].between(0,1).all():
        raise ValueError('Invalid p-values or variance in GeneOrdering')
    cutoff = .01 / len(order)
    passing = values.sort_values(values.columns[0])
    passing = passing.loc[passing.iloc[:,0] < cutoff]
    all_tfs = resolve_tf_names(expr.columns, tfs)
    present_tfs = [g for g in all_tfs if g in passing.index]
    if not present_tfs:
        raise ValueError('No TF passes the BEELINE p-value filter')
    non_tfs = passing.drop(index=present_tfs)
    top = non_tfs.sort_values(non_tfs.columns[1], ascending=False).iloc[:n_var].index
    keep = sorted(set(top) | set(present_tfs))
    return expr.loc[:,keep], sorted(present_tfs), {
        'n_cells': len(expr), 'n_genes_available': len(expr.columns),
        'n_variable_selected': len(top), 'n_tfs_present': len(present_tfs),
        'n_genes_kept': len(keep), 'n_tested': len(order),
        'p_value_cutoff_effective': cutoff, 'protocol': PROTOCOL,
    }


def main(argv=None):
    args = build_parser().parse_args(argv)
    from grnbench.report import environment_block, file_fingerprint, sequence_identity, require_empty_output, write_json
    if args.summarize:
        import json
        src = Path(args.summarize)
        result = summarize(json.loads(src.read_text(encoding='utf-8')), args.alpha)
        write_json(src, result)
        (src.parent / 'summary.md').write_text(render_summary_md(result), encoding='utf-8')
        rc, uc = result['regime_check'], result['unit_comparison']
        print(f"regime: {rc['n_above_regime']}/{rc['n_checked']} above {REGIME_FACTOR}x published; "
              f"units: {uc['n_resolved']}/{uc['n_pairs']} pairs resolved -> {src.parent / 'summary.md'}")
        return 0
    from grnbench.eval.gold import load_gold
    from grnbench.eval.metrics import GoldIndex, evaluate_ranking
    from grnbench.eval.nulls import null_metrics
    from grnbench.infer import list_methods, run_method, save_edges
    from grnbench.tfs import load_tf_list
    import anndata as ad
    try:
        out = require_empty_output(args.out)
    except ValueError as e:
        print(str(e), file=sys.stderr)
        return 1
    result = {'protocol': PROTOCOL, 'data_release': args.data_release,
              'args': vars(args), 'environment': environment_block(),
              'status':'running', 'datasets':{}, 'rows':[], 'errors':[], 'notes':[]}
    def persist():
        write_json(out / 'beeline_backtest.json', result)
        if result['rows']:
            pd.DataFrame(result['rows']).to_csv(out / 'beeline_backtest.csv', index=False)
    persist()
    try:
        if args.plugin:
            from grnbench.plugins import register_plugins
            result['plugins'] = register_plugins(args.plugin)
        wanted = list(dict.fromkeys(m.strip() for m in args.methods.split(',') if m.strip()))
        available = {m['name'] for m in list_methods() if m['available']}
        methods = [m for m in wanted if m in available]
        for m in wanted:
            if m not in available:
                result['errors'].append({'method':m, 'error':'Requested method unavailable'})
        if not methods:
            raise ValueError('No requested method is available')
        datasets = find_datasets(Path(args.data))
        if args.datasets:
            requested = {d.strip().casefold() for d in args.datasets.split(',')}
            missing = requested - {d.casefold() for d in datasets}
            if missing:
                raise ValueError(f'Requested datasets absent: {sorted(missing)}')
            datasets = {k:v for k,v in datasets.items() if k.casefold() in requested}
        if not datasets:
            raise ValueError('No BEELINE datasets found')
        if args.n_nulls < 1:
            raise ValueError('n_nulls must be positive')
        for name,d in sorted(datasets.items()):
            try:
                dsout = out/name
                dsout.mkdir()
                organism, netfile = SPECIFIC_NETWORK[name]
                tfpath = Path(args.networks)/f'{organism}-tfs.csv'
                expr,tfs,meta = beeline_subset(load_expression(d),d,load_tf_list(tfpath),args.n_var)
                meta.update({'expression_input':file_fingerprint(d/'ExpressionData.csv'),
                             'ordering_input':file_fingerprint(d/'GeneOrdering.csv'),
                             'tf_input':file_fingerprint(tfpath),
                             'effective_genes':sequence_identity(expr.columns),
                             'effective_tfs':sequence_identity(tfs), 'runs':{}, 'golds':{}})
                result['datasets'][name] = meta
                universe = pd.DataFrame(((s,t) for s in tfs for t in expr.columns if s!=t), columns=['source','target'])
                upath = dsout/'candidate_universe.tsv.gz'
                universe.to_csv(upath,sep='\t',index=False,compression='gzip')
                meta['universe'] = file_fingerprint(upath)
                adata = ad.AnnData(expr.to_numpy(dtype=np.float32),
                                  obs=pd.DataFrame(index=expr.index.astype(str)),
                                  var=pd.DataFrame(index=expr.columns))
                apath = dsout/'selected_expression.h5ad'
                adata.write_h5ad(apath, compression='gzip')
                meta['selected_expression'] = file_fingerprint(apath)
                networks = [('specific-chip',f'{organism}/{netfile}')]
                networks += [(k,f'{organism}/{f}') for k,f in GENERIC_NETWORKS[organism]]
                networks += EXTRA_NETWORKS.get(name,[])
                golds = {}
                genes = {g.upper() for g in expr.columns}
                for kind,rel in networks:
                    gp = Path(args.networks)/'Networks'/rel
                    if not gp.exists():
                        result['errors'].append({'dataset':name,'ground_truth':kind,'error':f'Missing network: {gp}'})
                        continue
                    gold = load_gold(str(gp))
                    pairs = sorted((s,t) for s,t in gold.pairs if s.upper() in genes and t.upper() in genes and s.upper()!=t.upper())
                    goldpath = dsout/f'{kind}.gold.tsv.gz'
                    pd.DataFrame(pairs,columns=['source','target']).to_csv(goldpath,sep='\t',index=False,compression='gzip')
                    meta['golds'][kind] = {'input':file_fingerprint(gp),'effective':file_fingerprint(goldpath),'n_pairs':len(pairs)}
                    if not pairs:
                        result['errors'].append({'dataset':name,'ground_truth':kind,'error':'No gold pairs in selected genes'})
                        continue
                    golds[kind] = GoldIndex(frozenset(pairs))
                for method in methods:
                    try:
                        start = time.monotonic()
                        edges,prov = run_method(method,adata,tf_list=tfs,assume_normalized=True,
                                                seed=args.seed,threads=args.threads)
                        edgepath = save_edges(edges,dsout,method)
                        prov.update({'edges_file':str(edgepath),'edges_fingerprint':file_fingerprint(edgepath),
                                     'infer_seconds':time.monotonic()-start})
                        meta['runs'][method] = prov
                        write_json(dsout/f'{method}.provenance.json',prov)
                        persist()
                        for kind,gold in golds.items():
                            try:
                                metrics = evaluate_ranking(edges,gold,universe=universe,with_curve=False)
                                null = null_metrics(edges,gold,n_nulls=args.n_nulls,mode='random',
                                                    seed=args.seed,universe=universe)
                                result['rows'].append({'dataset':name,'method':method,'ground_truth':kind,
                                                       **metrics,'random_ranking_null':null})
                            except Exception as e:
                                result['errors'].append({'dataset':name,'method':method,'ground_truth':kind,'error':f'{type(e).__name__}: {e}'})
                    except Exception as e:
                        result['errors'].append({'dataset':name,'method':method,'error':f'{type(e).__name__}: {e}'})
                    persist()
            except Exception as e:
                result['errors'].append({'dataset':name,'error':f'{type(e).__name__}: {e}'})
            persist()
    except Exception as e:
        result['errors'].append({'stage':'setup','error':f'{type(e).__name__}: {e}'})
    result['status'] = ('partial_failure' if result['rows'] else 'failed') if result['errors'] else 'complete'
    try:
        summarize(result, args.alpha)
        (out / 'summary.md').write_text(render_summary_md(result), encoding='utf-8')
    except Exception as e:
        result['errors'].append({'stage': 'summary', 'error': f'{type(e).__name__}: {e}'})
    persist()
    rc, uc = result.get('regime_check', {}), result.get('unit_comparison', {})
    print(f"{result['status']}: {len(result['rows'])} scored results; {len(result['errors'])} errors -> {out}")
    print(f"regime: {rc.get('n_above_regime', 0)}/{rc.get('n_checked', 0)} above {REGIME_FACTOR}x published; "
          f"units: {uc.get('n_resolved', 0)}/{uc.get('n_pairs', 0)} pairs resolved")
    return 1 if result['errors'] else 0


if __name__ == '__main__':
    raise SystemExit(main())
