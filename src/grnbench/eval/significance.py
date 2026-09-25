"""Descriptive annotation-thinning scenarios, not statistical inference.

Removed annotations become unlabelled. Quantiles describe these scenarios; they
are not confidence intervals. Independent-unit inference is in eval.comparison.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from grnbench.eval.metrics import GoldIndex, _tie_blocks, auc_roc, precision_at_k
from grnbench.eval.universe import align_rankings

_METRICS = ('epr', 'auprc', 'auprc_ratio', 'precision_at_k', 'auroc')


class _MethodState:
    def __init__(self, edges, gold_ids):
        self.weights = edges.weight.to_numpy(dtype=float)
        self.ends = _tie_blocks(self.weights)
        self.gold_id = np.array([gold_ids.get(p, -1) for p in zip(edges.source, edges.target)])

    def score(self, keep, k, metric):
        y = np.zeros(len(self.weights), dtype=np.int8)
        hit = self.gold_id >= 0
        y[hit] = keep[self.gold_id[hit]]
        n_pos = int(y.sum())
        if not n_pos:
            return float('nan')
        prevalence = n_pos / len(y)
        if metric == 'auroc':
            return auc_roc(self.weights, y)
        if metric in ('epr', 'precision_at_k'):
            precision = precision_at_k(y, k, self.weights)
            return precision / prevalence if metric == 'epr' else precision
        hits = np.cumsum(y)[self.ends - 1].astype(float)
        ap = float(np.sum(np.diff(np.r_[0., hits / n_pos]) * hits / self.ends))
        return ap / prevalence if metric == 'auprc_ratio' else ap


def gold_thinning_sensitivity(edge_sets: dict[str, pd.DataFrame], gold,
                              *, metric='epr', n_draws=200, frac=.5,
                              top_k=None, seed=42, universe=None) -> dict:
    """Uniformly retain frac of reachable gold; freeze k before all draws."""
    if metric not in _METRICS:
        raise ValueError(f'metric must be one of {_METRICS}')
    if isinstance(n_draws, bool) or int(n_draws) != n_draws or n_draws < 1:
        raise ValueError('n_draws must be a positive integer')
    if not np.isfinite(frac) or not 0 < frac <= 1:
        raise ValueError('frac must be in (0, 1]')
    if len(edge_sets) < 2:
        raise ValueError('comparing methods needs at least two edge lists')
    if any('|' in str(n) for n in edge_sets):
        raise ValueError('method names cannot contain |')
    gi = GoldIndex.coerce(gold)
    aligned = align_rankings(edge_sets, universe, case_sensitive=gi.case_sensitive)
    first = next(iter(aligned.values()))
    support = set(zip(first.source, first.target))
    gold_pairs = sorted(gi.lookup & support)
    n_gold = len(gold_pairs)
    if n_gold == 0:
        raise ValueError('no reachable gold positives in the candidate universe')
    if top_k is not None and (isinstance(top_k, bool) or int(top_k) != top_k or top_k <= 0):
        raise ValueError('top_k must be a positive integer')
    k = n_gold if top_k is None else min(int(top_k), len(first))
    gold_ids = {p: i for i, p in enumerate(gold_pairs)}
    states = {n: _MethodState(e, gold_ids) for n, e in aligned.items()}
    names = list(states)
    m = max(1, int(round(frac * n_gold)))
    rng = np.random.default_rng(seed)
    draws = {n: np.empty(int(n_draws), dtype=float) for n in names}
    for b in range(int(n_draws)):
        keep = np.zeros(n_gold, dtype=np.int8)
        keep[rng.choice(n_gold, size=m, replace=False)] = 1
        for n in names:
            draws[n][b] = states[n].score(keep, k, metric)
    full = np.ones(n_gold, dtype=np.int8)
    observed = {n: states[n].score(full, k, metric) for n in names}
    pairs = {}
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            d = draws[a] - draws[b]
            d = d[np.isfinite(d)]
            pairs[f'{a}|{b}'] = {
                'fraction_a_gt_b': float(np.mean(d > 0)) if d.size else float('nan'),
                'diff_mean': float(d.mean()) if d.size else float('nan'),
                'diff_quantiles': np.quantile(d, [.025, .975]).tolist() if d.size else [float('nan')] * 2,
                'n_valid_draws': int(d.size),
                'observed_difference': observed[a] - observed[b],
            }
    return {'analysis': 'gold-thinning-sensitivity', 'metric': metric,
            'n_draws': int(n_draws), 'frac': float(frac), 'n_gold_retained': m,
            'n_gold_reachable': n_gold, 'epr_k': k, 'seed': int(seed),
            'support_origin': first.attrs['evaluation_support']['support_origin'],
            'resample': f'uniform gold thinning ({m}/{n_gold} reachable pairs, no replacement)',
            'observed': observed, 'methods': names, 'pairs': pairs,
            'n_comparisons': len(pairs),
            'interpretation': 'Descriptive annotation-removal scenarios; quantiles are not confidence intervals; no superiority or equivalence inference.'}




def reference_ordering_differences(blocks):
    """Describe observed ordering changes between reference-specific tasks."""
    by_pair = {}
    for block in blocks:
        for key, rec in (block.get('pairs') or {}).items():
            a, b = key.split('|', 1)
            diff = rec.get('observed_difference', rec.get('diff_mean', float('nan')))
            if not np.isfinite(diff) or diff == 0:
                continue
            by_pair.setdefault((block.get('dataset', ''), tuple(sorted((a, b)))), {})[
                block.get('ground_truth', '')] = a if diff > 0 else b
    comparable = {k: v for k, v in by_pair.items() if len(v) >= 2}
    reversals = [{'dataset': ds, 'methods': list(pair), 'winner_by_ground_truth': winners}
                 for (ds, pair), winners in comparable.items() if len(set(winners.values())) > 1]
    return {'analysis': 'reference-ordering-differences', 'n_comparable_pairs': len(comparable),
            'n_reversed': len(reversals), 'reversals': reversals,
            'reversal_rate': len(reversals) / len(comparable) if comparable else float('nan'),
            'interpretation': 'Different references may measure different biological relations; ordering reversals do not establish a wrong reference or statistical contradiction.'}



def render_contradiction_md(rep):
    if not rep.get('n_comparable_pairs'):
        return ''
    return ('## Reference-specific ordering differences\n\n'
            f"{rep['n_reversed']} of {rep['n_comparable_pairs']} comparable method pairs "
            f"reverse observed order ({rep['reversal_rate']:.0%}).\n\n"
            + rep['interpretation'] + '\n')


def render_significance_md(res):
    lines = [f"## Gold-thinning sensitivity ({res['metric']})", '',
             res['interpretation'], '',
             f"{res['resample']}; {res['n_draws']} draws; fixed k={res['epr_k']}.", '',
             '| comparison | fraction A > B | mean difference | 2.5% / 97.5% scenario quantiles |',
             '|---|---|---|---|']
    for key, rec in res['pairs'].items():
        q = rec['diff_quantiles']
        lines.append(f"| {key} | {rec['fraction_a_gt_b']:.3f} | {rec['diff_mean']:+.3f} | {q[0]:+.3f} / {q[1]:+.3f} |")
    return '\n'.join(lines) + '\n'
