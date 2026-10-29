"""Cell subsampling is keyed by barcode, not by position (D-024)."""
import numpy as np

from grnbench.screen import stable_cell_order, stable_subsample

K = 3000


def _names(n=5000):
    return np.array([f"AAAC{i:06d}-1" for i in range(n)])


def test_removing_three_cells_changes_at_most_three_sampled_barcodes():
    names = _names()
    cut = np.delete(names, [10, 2500, 4999])
    a = set(names[stable_subsample(names, K, 0)])
    b = set(cut[stable_subsample(cut, K, 0)])
    assert len(a - b) <= 3 and len(b - a) <= 3
    # positional rule of origin/main, for contrast
    old = lambda n: set(n[np.sort(np.random.default_rng(0).choice(len(n), K, replace=False))])
    assert len(old(names) & old(cut)) < K - 500


def test_sample_does_not_depend_on_row_order():
    names = _names()
    shuffled = names[np.random.default_rng(1).permutation(len(names))]
    assert set(names[stable_subsample(names, K, 0)]) == set(shuffled[stable_subsample(shuffled, K, 0)])


def test_deterministic_per_seed_and_differs_across_seeds():
    names = _names()
    assert np.array_equal(stable_subsample(names, K, 0), stable_subsample(names, K, 0))
    assert not np.array_equal(stable_subsample(names, K, 0), stable_subsample(names, K, 1))
    idx = stable_subsample(names, K, 0)
    assert len(idx) == K and np.all(np.diff(idx) > 0)


def test_small_input_returns_all_positions():
    names = _names(200)
    assert np.array_equal(stable_subsample(names, K, 0), np.arange(200))


def test_duplicate_names_do_not_crash():
    names = np.array(["a", "b", "a", "b", "c"])
    order = stable_cell_order(names, 0)
    assert sorted(order) == [0, 1, 2, 3, 4]
    assert len(stable_subsample(names, 3, 0)) == 3


def test_golden_order_is_fixed_across_processes():
    """blake2b64 of "seed|barcode" is pinned: Python's salted hash() or another digest would fail here."""
    names = [f"AAAC{i:06d}-1" for i in range(10)]
    assert stable_cell_order(names, 0).tolist() == [9, 8, 0, 4, 1, 2, 5, 7, 3, 6]
    assert stable_cell_order(names, 1).tolist() == [1, 3, 5, 0, 8, 9, 2, 4, 7, 6]
    assert stable_subsample(names, 4, 0).tolist() == [0, 4, 8, 9]


def test_cell_sampling_rule_name():
    from grnbench.screen import CELL_SAMPLING
    assert CELL_SAMPLING == "blake2b64-seed-barcode-v1"
