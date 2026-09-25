"""Explicit conditional randomization diagnostics for ranked edge lists."""
from __future__ import annotations
import numpy as np
import pandas as pd
from grnbench.eval.metrics import GoldIndex, evaluate_ranking, metric_from_labels
from grnbench.eval.universe import align_ranking, validate_edges


def _count(n):
    if isinstance(n, bool) or int(n) != n or n < 1:
        raise ValueError('number of null draws must be a positive integer')
    return int(n)


def random_null(edges: pd.DataFrame, n: int, seed=42, *, case_sensitive=False):
    """Permute scores globally; preserve every candidate pair and score tie."""
    validated = validate_edges(edges, case_sensitive=case_sensitive)
    edges = edges.copy().reset_index(drop=True).assign(weight=validated.weight.to_numpy())
    rng = np.random.default_rng(seed)
    return [edges.assign(weight=rng.permutation(edges.weight.to_numpy()))
            .sort_values('weight', ascending=False, kind='stable').reset_index(drop=True)
            for _ in range(_count(n))]


def regulator_stratified_nulls(edges: pd.DataFrame, n: int, seed=42, *, case_sensitive=False):
    """Permute scores within each source: a distinct conditional hypothesis."""
    validated = validate_edges(edges, case_sensitive=case_sensitive)
    edges = edges.copy().reset_index(drop=True).assign(weight=validated.weight.to_numpy())
    rng = np.random.default_rng(seed)
    out = []
    groups = list(edges.groupby('source', sort=False).indices.values())
    for _ in range(_count(n)):
        weights = edges.weight.to_numpy(copy=True)
        for idx in groups:
            weights[idx] = rng.permutation(weights[idx])
        out.append(edges.assign(weight=weights).sort_values('weight', ascending=False,
                                                          kind='stable').reset_index(drop=True))
    return out


def degree_preserving_nulls(edges: pd.DataFrame, n: int, seed=42, *, case_sensitive=False):
    """Sparse simple directed double-edge swaps; diagnostic only.

    Both unique in/out degrees are preserved. Reject loops/duplicates. This
    bounded chain is not claimed to mix uniformly; immobile graphs raise.
    """
    validated = validate_edges(edges, case_sensitive=case_sensitive)
    edges = edges.copy().reset_index(drop=True).assign(weight=validated.weight.to_numpy())
    if (edges.source == edges.target).any():
        raise ValueError('degree rewiring forbids self-loops')
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(_count(n)):
        shuffled = edges.copy()
        source = shuffled.source.to_numpy()
        target = shuffled.target.to_numpy(copy=True)
        pairs = set(zip(source, target))
        swaps = 0
        desired = max(1, 5 * len(edges))
        for attempt in range(max(100, 50 * len(edges))):
            if len(edges) < 2:
                break
            i, j = rng.choice(len(edges), 2, replace=False)
            s, t, u, v = source[i], target[i], source[j], target[j]
            if s == u or t == v or s == v or u == t or (s, v) in pairs or (u, t) in pairs:
                continue
            pairs.remove((s, t)); pairs.remove((u, v))
            pairs.add((s, v)); pairs.add((u, t))
            target[i], target[j] = v, t
            swaps += 1
            if swaps >= desired:
                break
        if not swaps:
            raise ValueError('degenerate degree rewiring: no admissible simple-edge swap')
        shuffled['target'] = target
        shuffled.attrs['rewiring'] = {'accepted_swaps': swaps,
                                      'mixing_calibrated': False}
        out.append(shuffled)
    return out


def null_metrics(edges, gold, n_nulls=50, mode='random', seed=42,
                 metric='auprc', *, universe=None, top_k=None):
    """Monte Carlo score-permutation test conditional on fixed U and labels.

    Exchangeability is an explicit null hypothesis, not biological validation.
    Degree swaps change support and therefore return diagnostic scores only.

    For the ``random`` and ``regulator`` modes the pairs never change, so a
    permutation of the scores over fixed pairs is realised as the equivalent
    permutation of the labels over the fixed, already-ranked score vector: same
    distribution, same tie blocks, and no identifier re-validation per draw
    (which was ~6 s of the ~6.5 s each draw cost on a 370k-edge list).
    """
    n_nulls = _count(n_nulls)
    gi = GoldIndex.coerce(gold)
    observed = evaluate_ranking(edges, gi, top_k=top_k, with_curve=False,
                                universe=universe)[metric]
    aligned = align_ranking(edges, universe, case_sensitive=gi.case_sensitive)
    generators = {'random': random_null, 'regulator': regulator_stratified_nulls,
                  'degree': degree_preserving_nulls}
    if mode not in generators:
        raise ValueError(f'unknown null mode {mode!r}')
    # Use the observed support exactly even if gold is inaccessible; don't
    # discard selected draws and silently condition the randomization test.
    y = gi.labels(aligned).astype(np.int64)
    fixed_k = top_k if top_k is not None else int(y.sum())
    values = []
    unchanged = 0
    if mode == 'degree':
        original = aligned.set_index(['source', 'target']).weight.sort_index()
    else:
        w = aligned.weight.to_numpy(dtype=float)          # ranking order
        k = min(int(fixed_k), len(y))
        groups = (None if mode == 'random'
                  else list(aligned.groupby('source', sort=False).indices.values()))
    # One draw at a time; materializing n_nulls million-edge frames multiplies
    # peak RAM needlessly. Child RNG streams remain reproducible.
    for child_seed in np.random.SeedSequence(seed).spawn(n_nulls):
        if mode != 'degree':
            rng = np.random.default_rng(child_seed)
            if groups is None:
                y_perm = y[rng.permutation(len(y))]
            else:
                y_perm = y.copy()
                for idx in groups:
                    y_perm[idx] = y[rng.permutation(idx)]
            unchanged += int(np.array_equal(y_perm, y))
            values.append(metric_from_labels(metric, w, y_perm, k))
            continue
        try:
            sample = generators[mode](aligned, 1, seed=child_seed,
                                      case_sensitive=gi.case_sensitive)[0]
        except ValueError as exc:
            if 'degenerate' not in str(exc):
                raise
            continue
        current = sample.set_index(['source', 'target']).weight.sort_index()
        unchanged += int(current.equals(original))
        values.append(evaluate_ranking(sample, gi, top_k=fixed_k or None,
                                       with_curve=False)[metric])
    arr = np.asarray(values, dtype=float)
    valid = arr[np.isfinite(arr)]
    undefined = not np.isfinite(observed)
    inference = mode != 'degree' and not undefined and len(valid) == n_nulls
    p = (int(np.sum(valid >= observed)) + 1) / (n_nulls + 1) if inference else float('nan')
    hypothesis = {'random': 'global score exchangeability on fixed candidate universe',
                  'regulator': 'score exchangeability within each source on fixed candidate universe',
                  'degree': 'uncalibrated simple-degree rewiring diagnostic; changing support'}[mode]
    return {'metric': metric, 'null_mode': mode, 'null_hypothesis': hypothesis,
            'seed': int(seed), 'n_nulls_requested': n_nulls,
            'observed': float(observed), 'observed_undefined': bool(undefined),
            'null_mean': float(valid.mean()) if valid.size else float('nan'),
            'null_std': float(valid.std()) if valid.size else float('nan'),
            'p_value': float(p), 'inference_available': bool(inference),
            'n_nulls_used': int(valid.size), 'n_nulls_degenerate': n_nulls - int(valid.size),
            'n_unchanged_draws': unchanged,
            'interpretation': 'Conditional randomization diagnostic; unadjusted across methods; no biological validation.'}
