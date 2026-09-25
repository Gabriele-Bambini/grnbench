"""Independent-unit inference: reference calculations, not gold resampling."""
import json

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from grnbench.eval.comparison import compare_units, holm_adjust, scores_from_long


def test_pair_matches_scipy_and_two_directed_hypotheses():
    frame = pd.DataFrame({"a": [2., 4., 4., 6., 7.], "b": [1., 1., 2., 3., 4.]})
    result = compare_units(frame)
    rec = result["pairs"][0]
    expected = stats.ttest_rel(frame.a, frame.b, alternative="greater")
    assert rec["difference"] == pytest.approx((frame.a - frame.b).mean())
    assert rec["p_a_gt_b"] == pytest.approx(expected.pvalue)
    assert rec["p_a_gt_b_adjusted"] == pytest.approx(min(1., 2 * expected.pvalue))
    expected_ci = stats.ttest_rel(frame.a, frame.b).confidence_interval()
    assert rec["simultaneous_ci"] == pytest.approx([expected_ci.low, expected_ci.high])
    assert result["n_directed_hypotheses"] == 2
    json.dumps(result, allow_nan=False)


def test_holm_matches_known_stepdown_and_original_order():
    assert holm_adjust([.04, .01, .03, .9]) == pytest.approx([.09, .04, .09, .9])
    with pytest.raises(ValueError):
        holm_adjust([np.nan])


def test_zero_variance_cannot_manufacture_certainty():
    scores = pd.DataFrame({"a": [2., 2., 2., 2.], "b": [1., 1., 1., 1.]})
    pair = compare_units(scores)["pairs"][0]
    assert pair["winner"] is None
    assert pair["simultaneous_ci"] == [None, None]
    assert pair["status"] == "zero_variance_unresolved"
    assert pair["p_a_gt_b_adjusted"] == 1.
    assert compare_units(scores.assign(b=scores.a))["pairs"][0]["winner"] is None


def test_practical_margin_is_tested_in_both_directions():
    scores = pd.DataFrame({"a": [1.1, 1.15, 1.2, 1.12, 1.16], "b": np.ones(5)})
    assert compare_units(scores, margin=0)["pairs"][0]["winner"] == "a"
    assert compare_units(scores, margin=.3)["pairs"][0]["winner"] is None
    rev = scores[["b", "a"]]
    assert compare_units(rev, margin=0)["pairs"][0]["winner"] == "a"


@pytest.mark.parametrize("frame", [
    pd.DataFrame({"a": [1, 2], "b": [3, 4]}),
    pd.DataFrame({"a": [1, 2, 3], "b": [3, np.nan, 4]}),
    pd.DataFrame({"a": [1, 2, 3], "b": [3, 5, 4]}, index=["u", "u", "v"]),
])
def test_invalid_experimental_units_rejected(frame):
    with pytest.raises(ValueError):
        compare_units(frame)


def test_missing_method_name_is_not_a_method_called_nan():
    with pytest.raises(ValueError, match="identifier"):
        compare_units(pd.DataFrame([[1, 2], [2, 4], [3, 7]], columns=["a", np.nan]))


def test_long_table_requires_unique_complete_unit_method_cells():
    table = pd.DataFrame({"unit": ["x", "x", "y", "y", "z", "z"],
                          "method": ["a", "b"] * 3, "score": range(6)})
    wide = scores_from_long(table)
    assert wide.shape == (3, 2)
    with pytest.raises(ValueError, match="duplicate"):
        scores_from_long(pd.concat([table, table.head(1)]))
    with pytest.raises(ValueError, match="missing"):
        scores_from_long(table.iloc[:-1])


def test_comparison_cli_requires_unit_declaration(tmp_path):
    from grnbench.cli_compare import main
    path = tmp_path / "scores.csv"
    pd.DataFrame({"unit": ["x", "x", "y", "y", "z", "z"],
                  "method": ["a", "b"] * 3, "score": [1, 2, 2, 4, 3, 5]}).to_csv(path, index=False)
    with pytest.raises(SystemExit):
        main([str(path), "--out", str(tmp_path / "missing")])
    out = tmp_path / "valid"
    assert main([str(path), "--out", str(out), "--unit-description", "independent experiments",
                 "--assume-independent-units"]) == 0
    result = json.loads((out / "comparison.json").read_text())
    assert result["unit_description"] == "independent experiments"
    assert "independent" in result["assumptions"][0]
    assert (out / "scores.csv").exists()
