"""Power of the paired t over datasets (grnbench.eval.power): checks against the G*Power
manual and against simulation, not against the formula itself."""
import numpy as np
import pytest
from scipy import stats

from grnbench.eval.power import SampleSizeNotReached, min_detectable_dz, n_for_power, paired_t_power


def test_two_sided_power_matches_the_gpower_manual_matched_pairs_example():
    # G*Power 3.1 manual, section 19 (matched pairs): dz 0.421637, N 50, two tails, alpha .05
    assert paired_t_power(0.421637, 50, 0.05, sides=2) == pytest.approx(0.832114, abs=5e-7)


@pytest.mark.parametrize("sides", [1, 2])
def test_power_at_zero_effect_is_the_test_size(sides):
    assert paired_t_power(0.0, 20, 0.05, sides=sides) == pytest.approx(0.05, abs=1e-12)


def test_two_sided_power_does_not_depend_on_the_sign_of_the_effect():
    assert paired_t_power(-0.3, 40, 0.05, sides=2) == pytest.approx(paired_t_power(0.3, 40, 0.05, sides=2))


def test_one_sided_power_agrees_with_simulated_paired_t_tests():
    # the directed test of compare_units: p_a_gt_b = t.sf(mean / se, n - 1)
    rng = np.random.default_rng(42)
    n, dz, alpha, reps = 20, 0.5, 0.05 / 3, 40_000
    diffs = rng.normal(dz, 1.0, size=(reps, n))
    tstat = diffs.mean(axis=1) / (diffs.std(axis=1, ddof=1) / np.sqrt(n))
    simulated = float((stats.t.sf(tstat, n - 1) < alpha).mean())
    expected = paired_t_power(dz, n, alpha, sides=1)
    mc_se = np.sqrt(expected * (1 - expected) / reps)
    assert abs(simulated - expected) < 4 * mc_se


@pytest.mark.parametrize("dz, n, alpha", [(0.5, 200, 0.0125), (2.0, 20, 0.0125), (5.0, 24, 0.0125)])
def test_two_sided_power_stays_finite_for_large_effects(dz, n, alpha):
    # the lower tail at a large noncentrality must not turn the whole power into NaN
    two = paired_t_power(dz, n, alpha, sides=2)
    assert np.isfinite(two) and paired_t_power(dz, n, alpha / 2) <= two <= 1


@pytest.mark.parametrize("dz", [3.0, 5.0, 8.0])
def test_two_sided_n_for_power_is_the_smallest_sufficient_n_for_large_effects(dz):
    n = n_for_power(dz, 0.0125, 0.8, sides=2)
    assert paired_t_power(dz, n, 0.0125, sides=2) >= 0.8
    assert n == 3 or paired_t_power(dz, n - 1, 0.0125, sides=2) < 0.8


def test_two_sided_min_detectable_dz_inverts_power_at_large_n():
    d = min_detectable_dz(160, 0.0125, 0.8, sides=2)
    assert paired_t_power(d, 160, 0.0125, sides=2) == pytest.approx(0.8, abs=1e-9)


def test_n_for_power_returns_the_smallest_sufficient_n():
    for dz, alpha, target in [(0.33889917, 0.05 / 3, 0.8), (0.2817, 0.025, 0.9), (1.2, 0.05, 0.8)]:
        n = n_for_power(dz, alpha, target)
        assert paired_t_power(dz, n, alpha) >= target
        assert n == 3 or paired_t_power(dz, n - 1, alpha) < target


@pytest.mark.parametrize("dz", [8.0, 50.0])
def test_n_for_power_never_plans_fewer_units_than_compare_units_accepts(dz):
    # compare_units refuses fewer than three units, even when two would reach the power
    assert paired_t_power(dz, 2, 0.05) >= 0.8
    assert n_for_power(dz, 0.05, 0.8) == 3


def test_n_for_power_refuses_a_cap_below_three_units():
    with pytest.raises(ValueError, match="n_max"):
        n_for_power(50.0, 0.05, 0.8, n_max=2)


@pytest.mark.parametrize("call", [
    lambda: paired_t_power(0.3, "20", 0.05),
    lambda: min_detectable_dz("20", 0.05, 0.8),
    lambda: n_for_power(0.3, 0.05, 0.8, n_max="500"),
    lambda: paired_t_power(0.3, None, 0.05),
])
def test_a_number_of_units_that_is_not_a_number_is_a_value_error(call):
    with pytest.raises(ValueError, match="integer"):
        call()


def test_min_detectable_dz_inverts_power_and_shrinks_with_n():
    d60 = min_detectable_dz(60, 0.025, 0.8)
    d200 = min_detectable_dz(200, 0.025, 0.8)
    assert paired_t_power(d60, 60, 0.025) == pytest.approx(0.8, abs=1e-9)
    assert paired_t_power(d200, 200, 0.025) == pytest.approx(0.8, abs=1e-9)
    assert 0 < d200 < d60


@pytest.mark.parametrize("call", [
    lambda: paired_t_power(0.3, 1, 0.05),            # one unit: no variance estimate
    lambda: paired_t_power(0.3, 20, 0.0),
    lambda: paired_t_power(0.3, 20, 1.0),
    lambda: paired_t_power(float("nan"), 20, 0.05),
    lambda: paired_t_power(0.3, 20, 0.05, sides=3),
    lambda: paired_t_power(True, 20, 0.05),           # a flag is not an effect size
    lambda: n_for_power(0.0, 0.05, 0.8),              # no positive effect: no finite n
    lambda: n_for_power(-0.2, 0.05, 0.8),
    lambda: n_for_power(0.3, 0.05, 0.05),             # target power must exceed alpha
    lambda: n_for_power(0.3, 0.05, 1.0),
    lambda: min_detectable_dz(20, 0.05, 0.04),
    lambda: min_detectable_dz(1, 0.05, 0.8),
])
def test_invalid_inputs_raise_instead_of_returning_a_number(call):
    with pytest.raises(ValueError):
        call()


def test_n_for_power_refuses_when_the_cap_is_reached():
    with pytest.raises(SampleSizeNotReached, match="n_max"):
        n_for_power(0.01, 0.05, 0.8, n_max=1000)


@pytest.mark.filterwarnings("ignore:.*did not converge:RuntimeWarning")
def test_min_detectable_dz_refuses_a_root_that_does_not_give_the_target_power():
    # at n = 2 and a tiny alpha scipy's noncentral t loses accuracy; a wrong dz must not come back
    with pytest.raises(ValueError, match="does not reach"):
        min_detectable_dz(2, 1e-6, 0.8)
