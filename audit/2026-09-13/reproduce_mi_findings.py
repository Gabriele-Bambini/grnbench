"""Read-only source audit: MI/CLR must not infer dependence from row indices.

Run from any directory. Loads classical.py directly, without importing scanpy
or the grnbench package. Writes only a JSON evidence file beside this script.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    source = root / "src/grnbench/infer/classical.py"
    spec = importlib.util.spec_from_file_location("audited_classical", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    zero = pd.DataFrame(np.zeros((16, 3)), columns=["A", "B", "C"])
    zero_mi = module.mutual_info(zero, tf_names=["A"], n_bins=8)
    zero_clr = module.clr(zero, tf_names=["A"], n_bins=8)
    # A constant random variable has entropy zero and MI zero with any variable.
    mi_ab = float(zero_mi.loc[zero_mi.target == "B", "weight"].iloc[0])
    clr_ab = float(zero_clr.loc[zero_clr.target == "B", "weight"].iloc[0])

    rng = np.random.default_rng(3)
    sparse = pd.DataFrame(rng.binomial(1, 0.1, (1000, 3)), columns=["A", "B", "C"])
    original = module.mutual_info(sparse, tf_names=["A"], n_bins=8)
    # Apply the SAME row permutation to all genes. The paired observations are
    # unchanged; dependence estimates must be invariant to their ordering.
    permuted = module.mutual_info(sparse.sample(frac=1, random_state=5),
                                  tf_names=["A"], n_bins=8)
    original = original.loc[original.target != "A"].reset_index(drop=True)
    permuted = permuted.loc[permuted.target != "A"].reset_index(drop=True)
    invariant = bool(np.allclose(original.weight, permuted.weight, atol=1e-10))
    original_order = original.sort_values("weight", ascending=False).target.tolist()
    permuted_order = permuted.sort_values("weight", ascending=False).target.tolist()

    evidence = {
        "source": str(source),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "zero_fixture": {
            "shape": [16, 3], "n_bins": 8,
            "expected_mutual_information_A_B": 0.0,
            "actual_mutual_information_A_B": mi_ab,
            "actual_clr_A_B": clr_ab,
            "mi_equals_maximum_bin_entropy": bool(np.isclose(mi_ab, np.log(8))),
            "constant_variable_property_passes": bool(np.isclose(mi_ab, 0)),
        },
        "common_row_permutation_fixture": {
            "shape": [1000, 3], "distribution": "independent Bernoulli(p=0.1)",
            "generation_seed": 3, "permutation_seed": 5,
            "original_edges": original.to_dict(orient="records"),
            "permuted_edges": permuted.to_dict(orient="records"),
            "cell_order_invariance_passes": invariant,
            "original_target_order": original_order,
            "permuted_target_order": permuted_order,
            "ranking_flips": original_order != permuted_order,
        },
        "interpretation": "Stable rank tie breaking shares row-index information across genes; zero ties become artificial dependence. No real-data inference was run.",
    }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(evidence, indent=2), encoding="utf-8")
    print(json.dumps(evidence, indent=2))


if __name__ == "__main__":
    main()
