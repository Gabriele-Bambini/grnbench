"""Power of the paired t over independent units (datasets), for planning experiments.

Model: G*Power's matched-pairs t test (Faul et al. 2007, Behav Res Methods 39:175,
doi:10.3758/BF03193146, Table 3; G*Power 3.1 manual, section 19). With n paired units and
dz = mean paired difference / SD of the differences, the statistic is t(n - 1) under H0
and noncentral t(n - 1, dz * sqrt(n)) under H1. ``sides=1`` is the directed test of
``compare_units`` (``p_a_gt_b``); plan a Holm family by passing the level of the step
the hypothesis is expected to face (alpha / k for the first step of k hypotheses).

These numbers plan a new experiment. Power computed from an effect estimated on the same
data is a function of its p-value and says nothing new about that data; report the
smallest detectable effect at the achieved n instead.
"""
from __future__ import annotations

import math
from numbers import Real

from scipy import optimize, stats


class SampleSizeNotReached(ValueError):
    """The target power needs more units than the cap allows."""


def _finite(x) -> bool:
    return isinstance(x, Real) and not isinstance(x, bool) and math.isfinite(x)


def _check_level(alpha: float, sides: int) -> None:
    if not _finite(alpha) or not 0 < alpha < 1:
        raise ValueError("alpha must be between zero and one")
    if isinstance(sides, bool) or sides not in (1, 2):
        raise ValueError("sides must be 1 (directed) or 2")


MIN_UNITS = 3   # compare_units refuses fewer than three units, so no plan asks for fewer


def _check_n(n: int, minimum: int = 2, name: str = "n") -> int:
    if not _finite(n) or not float(n).is_integer() or n < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum} paired units, got {n!r}")
    return int(n)


def _check_target(power: float, alpha: float) -> None:
    if not _finite(power) or not alpha < power < 1:
        raise ValueError("target power must lie strictly between alpha and one")


def paired_t_power(dz: float, n: int, alpha: float, *, sides: int = 1) -> float:
    """Probability that the paired t rejects H0 at level ``alpha`` when the true effect is ``dz``."""
    _check_level(alpha, sides)
    n = _check_n(n)
    if not _finite(dz):
        raise ValueError("dz must be a finite number")
    df, noncentrality = n - 1, dz * math.sqrt(n)
    if sides == 1:
        power = stats.nct.sf(stats.t.ppf(1 - alpha, df), df, noncentrality)
    else:
        # lower tail P(T < -c | nc) written as P(T > c | -nc): scipy's nct.cdf returns NaN
        # there once the noncentrality is large
        critical = stats.t.ppf(1 - alpha / 2, df)
        power = stats.nct.sf(critical, df, noncentrality) + stats.nct.sf(critical, df, -noncentrality)
    power = float(power)
    if not math.isfinite(power):
        raise ValueError(f"power is not finite for dz={dz}, n={n}, alpha={alpha}, sides={sides}")
    return power


def n_for_power(dz: float, alpha: float, power: float, *, sides: int = 1, n_max: int = 100_000) -> int:
    """Smallest number of paired units, at least ``MIN_UNITS`` (3), whose power for ``dz``
    reaches ``power``: ``compare_units`` refuses fewer than three units, so a larger effect
    still gets 3.

    Raises instead of returning a placeholder when ``dz`` is not positive or ``n_max``
    units are still not enough.
    """
    _check_level(alpha, sides)
    _check_target(power, alpha)
    if not _finite(dz) or dz <= 0:
        raise ValueError("dz must be positive: no finite n detects a null or reversed effect")
    n_max = _check_n(n_max, MIN_UNITS, "n_max")
    reaches = lambda n: paired_t_power(dz, n, alpha, sides=sides) >= power  # noqa: E731
    if reaches(MIN_UNITS):
        return MIN_UNITS
    low, high = MIN_UNITS, min(2 * MIN_UNITS, n_max)
    while not reaches(high):
        if high >= n_max:
            raise SampleSizeNotReached(f"power {power} for dz {dz} needs more than n_max={n_max} units")
        low, high = high, min(2 * high, n_max)
    while high - low > 1:  # reaches(high) and not reaches(low)
        middle = (low + high) // 2
        low, high = (low, middle) if reaches(middle) else (middle, high)
    return high


def min_detectable_dz(n: int, alpha: float, power: float, *, sides: int = 1) -> float:
    """Smallest positive dz that ``n`` paired units detect with probability ``power``."""
    _check_level(alpha, sides)
    _check_target(power, alpha)
    n = _check_n(n)
    gap = lambda dz: paired_t_power(dz, n, alpha, sides=sides) - power  # noqa: E731
    high = 1.0
    while gap(high) < 0:
        high *= 2
        if high > 1e6:
            raise ValueError(f"no dz up to 1e6 reaches power {power} with n={n}")
    dz = float(optimize.brentq(gap, 0.0, high, xtol=1e-12))
    if abs(gap(dz)) > 1e-6:   # scipy's noncentral t loses accuracy at tiny alpha and n
        raise ValueError(f"the root dz={dz} does not reach power {power} with n={n}, alpha={alpha}")
    return dz
