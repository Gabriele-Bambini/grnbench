"""Subsampling stability: do top-k edges survive resampling the cells?

For each subsample the method is re-run from scratch; stability is the mean
pairwise Jaccard of the top-k edge sets including boundary ties. Instability
describes sensitivity, not proof of noise or biological incorrectness.
"""
from __future__ import annotations

import numpy as np

from grnbench.infer import run_method
from grnbench.eval.universe import top_pairs


def _top_pairs(edges, k: int) -> set[tuple[str, str]]:
    return top_pairs(edges, k)


def subsample_stability(adata,
                        method: str,
                        tf_list: list[str] | None = None,
                        *,
                        frac: float = 0.8,
                        n_subsamples: int = 5,
                        top_k: int = 100,
                        assume_normalized: bool = False,
                        seed: int = 42,
                        algorithm_seed: int = 42,
                        vary_algorithm_seed: bool = False,
                        **method_params) -> dict:
    if not np.isfinite(frac) or not 0 < frac <= 1:
        raise ValueError('frac must be in (0, 1]')
    if int(n_subsamples) != n_subsamples or n_subsamples < 2:
        raise ValueError('n_subsamples must be an integer >= 2')
    if top_k <= 0 or int(top_k) != top_k:
        raise ValueError('top_k must be a positive integer')
    rng = np.random.default_rng(seed)
    n_cells = adata.n_obs
    if n_cells < 2:
        raise ValueError('stability needs at least two cells')
    n_keep = min(n_cells, max(2, int(round(frac * n_cells))))

    pair_sets: list[set] = []
    sizes: list[int] = []
    for i in range(int(n_subsamples)):
        idx = np.sort(rng.choice(n_cells, size=n_keep, replace=False))
        sub = adata[idx].copy()
        edges, _prov = run_method(method, sub, tf_list=tf_list,
                                  assume_normalized=assume_normalized,
                                  seed=algorithm_seed + i if vary_algorithm_seed else algorithm_seed,
                                  **method_params)
        pair_sets.append(_top_pairs(edges, top_k))
        sizes.append(len(edges))

    jaccards = []
    for i in range(len(pair_sets)):
        for j in range(i + 1, len(pair_sets)):
            union = pair_sets[i] | pair_sets[j]
            inter = pair_sets[i] & pair_sets[j]
            jaccards.append(len(inter) / len(union) if union else float("nan"))
    arr = np.asarray(jaccards, dtype=float)
    return {
        "method": method,
        "cell_resampling_seed": int(seed),
        "algorithm_seed": int(algorithm_seed),
        "vary_algorithm_seed": bool(vary_algorithm_seed),
        "tie_policy": "include-complete-boundary-tie",
        "n_subsamples": int(n_subsamples),
        "subsample_frac": frac,
        "top_k": int(top_k),
        "jaccard_mean": float(arr.mean()) if arr.size else float("nan"),
        "jaccard_std": float(arr.std()) if arr.size else float("nan"),
        "jaccard_all": [round(float(x), 6) for x in arr],
        "edge_counts": sizes,
    }
