"""P11 for EXP-20261016-LM-1: smallest detectable effect of the paired t over unit TFs, before any score.

Model (G*Power matched pairs, as in EXP-20261020-GB-1's eval/power.py, PR #129): with n units and
dz = mean paired difference / SD of the differences, t ~ noncentral t(n - 1, dz * sqrt(n)) under H1, one-sided.
Holm over the family's directed hypotheses: a true effect faces level alpha / k at worst (first step, k directed
hypotheses) and alpha at best (last step). A result needs the same verdict in both screens, and the screens share
no cell, so the joint power is the product. The permutation p-value (the second condition of "supported") is not
modelled. Unit counts: check P9 outcome (usable TFs under `core`); the scorer re-applies the >= 11 responder rule
inside G, so the analysed n can only be smaller. usage: python p11_power.py
"""
import json
import math

from scipy import optimize, stats

N = {"rpe1": 73, "k562": 109, "both_stratum": 41}
ALPHA, TARGET = 0.05, 0.80
FAMILIES = {"with_regdiffusion": 24, "without_regdiffusion": 16}


def power(dz, n, alpha):
    return float(stats.nct.sf(stats.t.ppf(1 - alpha, n - 1), n - 1, dz * math.sqrt(n)))


def mde(n, alpha, target=TARGET):
    return float(optimize.brentq(lambda d: power(d, n, alpha) - target, 1e-6, 5.0))


def joint_mde(n1, n2, alpha, target=TARGET):
    return float(optimize.brentq(lambda d: power(d, n1, alpha) * power(d, n2, alpha) - target, 1e-6, 5.0))


# instrument checks: G*Power's matched-pairs example (dz 0.5, one-sided alpha 0.05: n 27 gives power 0.8118),
# and no effect gives power alpha
assert abs(power(0.5, 27, 0.05) - 0.8118) < 0.001, power(0.5, 27, 0.05)
assert abs(power(0.0, 73, 0.05) - 0.05) < 1e-9
out = {"target_power": TARGET, "n": N, "mde_dz": {}}
for fam, k in FAMILIES.items():
    for step, a in (("worst_step", ALPHA / k), ("last_step", ALPHA)):
        key = f"{fam}:{step}"
        out["mde_dz"][key] = {s: round(mde(n, a), 3) for s, n in N.items()}
        out["mde_dz"][key]["joint_rpe1_and_k562"] = round(joint_mde(N["rpe1"], N["k562"], a), 3)
print(json.dumps(out, indent=1))
