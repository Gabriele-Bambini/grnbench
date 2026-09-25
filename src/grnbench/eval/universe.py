"""Fixed candidate support and strict identifier/score validation.

Missing predictions receive a single bottom tie; they remain in all denominators.
Each input namespace must map injectively under the configured case convention.
Case differences between separate input namespaces are permitted.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd


def canonical_pairs(frame, *, case_sensitive=False, name='pairs', allow_self_loops=False) -> pd.DataFrame:
    df = (frame[['source', 'target']].copy() if isinstance(frame, pd.DataFrame)
          else pd.DataFrame(list(frame), columns=['source', 'target']))
    if df.isna().any().any():
        raise ValueError(f'{name}: missing gene identifier')
    for col in ('source', 'target'):
        df[col] = df[col].astype(str)
        if df[col].str.strip().eq('').any():
            raise ValueError(f'{name}: empty gene identifier')
    if not case_sensitive:
        originals = pd.concat([df.source, df.target]).unique()
        folded = pd.Index(originals).str.upper()
        if folded.has_duplicates:
            raise ValueError(f'{name}: ambiguous identifier case collision')
        df = df.apply(lambda x: x.str.upper())
    if df.duplicated().any():
        raise ValueError(f'{name}: duplicate canonical pairs')
    if not allow_self_loops and (df.source == df.target).any():
        raise ValueError(f'{name}: self-loops are outside the candidate contract')
    return df.reset_index(drop=True)


def validate_edges(edges, *, case_sensitive=False, allow_empty=False):
    if not {'source', 'target', 'weight'} <= set(edges.columns):
        raise ValueError('edge list requires source, target, weight')
    if edges.empty and not allow_empty:
        raise ValueError('empty edge list')
    pairs = canonical_pairs(edges, case_sensitive=case_sensitive, name='edges')
    if np.iscomplexobj(edges['weight'].to_numpy()):
        raise ValueError('edge weights must be finite real numeric values')
    try:
        weights = edges['weight'].to_numpy(dtype=float)
    except (ValueError, TypeError) as exc:
        raise ValueError('edge weights must be finite numeric values') from exc
    if not np.isfinite(weights).all():
        raise ValueError('edge weights must be finite numeric values')
    pairs['weight'] = weights
    return pairs


def load_universe(path, *, case_sensitive=False):
    """Read a CSV/TSV (optionally gzipped) with source,target headers."""
    path = Path(path)
    sep = ',' if '.csv' in path.suffixes else '\t'
    return canonical_pairs(pd.read_csv(path, sep=sep),
                           case_sensitive=case_sensitive, name='universe')


def align_ranking(edges, universe=None, *, case_sensitive=False):
    clean = validate_edges(edges, case_sensitive=case_sensitive,
                           allow_empty=universe is not None)
    n_input = len(clean)
    if universe is None:
        aligned = clean
        origin = 'conditional-on-output-list'
    else:
        pairs = canonical_pairs(universe, case_sensitive=case_sensitive, name='universe')
        if pairs.empty:
            raise ValueError('empty candidate universe')
        support = set(map(tuple, pairs.to_numpy()))
        if not set(zip(clean.source, clean.target)) <= support:
            raise ValueError('prediction contains pairs outside the explicit universe')
        # Ranks preserve every observed score ordering/tie while leaving finite
        # zero below the weakest prediction, including -maxfloat inputs.
        clean['weight'] = clean.weight.rank(method='dense').to_numpy()
        aligned = pairs.merge(clean, on=['source', 'target'], how='left', validate='one_to_one')
        aligned['weight'] = aligned.weight.fillna(0.0)
        origin = 'explicit-universe'
    aligned = aligned.sort_values('weight', ascending=False, kind='stable').reset_index(drop=True)
    aligned.attrs['evaluation_support'] = {
        'support_origin': origin,
        'n_predictions': n_input,
        'missing_predictions': len(aligned) - n_input,
        'n_candidates': len(aligned),
        'missing_score_policy': 'tied-bottom',
    }
    return aligned


def align_rankings(edge_sets, universe=None, *, case_sensitive=False):
    """Align comparisons on explicit U, or require identical output supports."""
    aligned = {n: align_ranking(e, universe, case_sensitive=case_sensitive)
               for n, e in edge_sets.items()}
    if universe is None and len(aligned) > 1:
        supports = [set(zip(e.source, e.target)) for e in aligned.values()]
        if any(s != supports[0] for s in supports[1:]):
            raise ValueError('comparisons without --universe require identical candidate supports')
        for e in aligned.values():
            e.attrs['evaluation_support']['support_origin'] = 'shared-output-support'
    return aligned


def top_pairs(edges, k, *, case_sensitive=False):
    if isinstance(k, bool) or int(k) != k or k <= 0:
        raise ValueError('top_k must be a positive integer')
    df = validate_edges(edges, case_sensitive=case_sensitive)
    df = df.sort_values('weight', ascending=False, kind='stable')
    threshold = df.weight.iloc[min(int(k), len(df)) - 1]
    head = df[df.weight >= threshold]
    return set(zip(head.source, head.target))
