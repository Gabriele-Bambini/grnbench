"""C4's EPR lift of the motif prior, re-read from the committed rows: why n = 20 did not
resolve it, and how many never-used datasets a new study would need (EXP-20261020-GB-1).

Reads only the committed ``all_rows_real.csv`` and ``REPORT_REAL.md`` of EXP-20260925-LM-2,
whose bytes must match that run's ``SHA256SUMS``. Gate first: the registered code
(``real_bench.models_section``: ``compare_units`` per pair, Holm over the 4 directed hypotheses
of each metric and arm) must print exactly the EPR tables of the committed report, no line more
and no line less, otherwise nothing is computed. Then a descriptive diagnosis of the
per-dataset differences, the registered code rerun with each dataset left out (which single
datasets change C4's verdict) and a prospective power table (``grnbench.eval.power``).

Nothing here tests H1 again or changes C4: its EPR comparison stays not resolved, and not
resolved is not zero. No power is computed from C4's estimate: the achieved n = 20 is
described by the smallest detectable effect, and C4's estimate enters only as the number of
datasets a new study would need if the true effect equalled it.

    python benchmarks/epr_power.py --out benchmarks/results/EXP-20261020-GB-1
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import real_bench  # noqa: E402
from grnbench.eval.power import (SampleSizeNotReached, min_detectable_dz, n_for_power,  # noqa: E402
                                 paired_t_power)

ROOT = Path(__file__).resolve().parents[1]
C4_DIR = ROOT / "benchmarks" / "results" / "EXP-20260925-LM-2"
PAIR = ("pearson_motif", "pearson_signed")
ARMS = ("auto", "baseline")
REQUIRED = ["task", "reference_edges"] + [f"{arm}_n_contexts" for arm in ARMS] + [
    f"m_{model}_{arm}_{metric}" for arm in ARMS for model, metric in
    [("pearson_motif", "epr"), ("pearson_signed", "epr"), ("motif_prior", "epr"),
     ("pearson_motif", "auprc_ratio"), ("pearson_signed", "auprc_ratio")]]
ALPHA = 0.05
# Level the EPR hypothesis faces: first Holm step of C4's family of 4 directed hypotheses;
# C4's family after the other pair (p ~ 1e-4) is rejected first; a family of one pair.
LEVELS = {"alpha/4": ALPHA / 4, "alpha/3": ALPHA / 3, "alpha/2": ALPHA / 2}
GRID_LEVELS = ("alpha/4", "alpha/2")
EFFECTS = (0.05, 0.10, 0.12, 0.15, 0.20)   # EPR units; 0.10 = the plan's proposed smallest effect
GRID_N = (20, 60, 100, 120, 150, 200)
CURVE_N = np.arange(10, 251)   # x axis of the power-by-n figure
BOOTSTRAP_REPS, BOOTSTRAP_SEED = 20_000, 42
WILCOXON_ARE_NORMAL, WILCOXON_ARE_MIN = 3 / math.pi, 0.864   # Hodges & Lehmann 1956


class GateError(ValueError):
    """The committed rows no longer give the registered report: nothing else is computed."""


def load_rows(path: str | Path) -> pd.DataFrame:
    rows = pd.read_csv(path)
    missing = [c for c in REQUIRED if c not in rows]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")
    return rows


def _epr_tables(lines: list[str]) -> list[str]:
    """Headers and rows of the 'Declared comparisons, EPR' tables, in order."""
    out, i = [], 0
    while i < len(lines):
        if lines[i].startswith("#### Declared comparisons, EPR"):
            out.append(lines[i])
            i += 1
            while i < len(lines) and not lines[i].startswith("|"):
                i += 1
            while i < len(lines) and lines[i].startswith("|"):
                out.append(lines[i])
                i += 1
        else:
            i += 1
    return out


def registered_epr_lines(rows: pd.DataFrame) -> list[str]:
    """The EPR declared-comparison tables printed by the registered report code."""
    out = _epr_tables(real_bench.models_section(rows))
    if sum(line.startswith("#### Declared comparisons, EPR") for line in out) != len(ARMS):
        raise GateError(f"expected one EPR comparison table per arm, got: {out}")
    return out


def check_against_report(lines: list[str], report_text: str) -> None:
    """The code's EPR tables and the committed report's must be the same lines, in order."""
    committed = _epr_tables([line.rstrip("\r") for line in report_text.splitlines()])
    if lines != committed:
        raise GateError("EPR tables differ; printed by the code but not in the report: "
                        + " / ".join(x for x in lines if x not in committed)
                        + "; in the report but not printed: " + " / ".join(x for x in committed if x not in lines))


def _pair_lines(epr_lines: list[str]) -> dict[str, str]:
    """The line of PAIR in each arm's EPR table."""
    out, arm = {}, None
    for line in epr_lines:
        if line.startswith("#### Declared comparisons, EPR"):
            arm = next((a for a in ARMS if f", {a} pipeline" in line), None)
        elif arm and line.startswith(f"| {PAIR[0]} | {PAIR[1]} |"):
            out[arm] = line
    return out


def _verdict(line: str) -> str:
    """Last column of a declared-comparison line: the supported model or 'unresolved'."""
    return line.strip().strip("|").split("|")[-1].strip()


def leave_one_out(rows: pd.DataFrame, gate_lines: list[str]) -> dict:
    """C4's registered report code rerun once per dataset on all the others: which single
    datasets change the verdict on PAIR. An influence check of C4's own test, not a new test."""
    full = _pair_lines(gate_lines)
    if set(full) != set(ARMS):
        raise GateError(f"the gate lines lack the {PAIR} line of an arm: {sorted(full)}")
    out = {arm: {"verdict_all": _verdict(full[arm]), "drops": len(rows), "changed": []} for arm in ARMS}
    ids = rows["task"].str.replace("cxg__", "", regex=False)
    for i, dataset in zip(rows.index, ids):
        lines = _pair_lines(registered_epr_lines(rows.drop(index=i)))
        for arm in ARMS:
            if arm not in lines:
                raise GateError(f"without {dataset} the registered code prints no {arm} line for {PAIR}")
            if _verdict(lines[arm]) != out[arm]["verdict_all"]:
                out[arm]["changed"].append({"dropped": dataset, "line": lines[arm]})
    return out


def verify_inputs(rows_path: Path, report_path: Path) -> dict[str, str]:
    """sha256 of the input bytes as committed, checked against the run's SHA256SUMS."""
    sums_path = Path(rows_path).parent / "SHA256SUMS"
    if not sums_path.exists():
        raise GateError(f"{sums_path} not found: the inputs cannot be checked")
    registered = {}
    for line in sums_path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            digest, name = line.split(maxsplit=1)
            registered[name.lstrip("*").removeprefix("./")] = digest
    out = {}
    for path in (Path(rows_path), Path(report_path)):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if registered.get(path.name) != digest:
            raise GateError(f"{path.name}: sha256 {digest} does not match {sums_path} "
                            f"({registered.get(path.name, 'not listed')})")
        out[path.name] = digest
    return out


def _t(x: np.ndarray) -> float:
    return float(x.mean() / (x.std(ddof=1) / math.sqrt(len(x))))


def _dz(x: np.ndarray) -> float:
    return float(x.mean() / x.std(ddof=1))


def _bootstrap_dz(diff: np.ndarray) -> tuple[list[float], int]:
    """10th and 90th percentile of dz over resamples, and how many resamples had zero SD (left out)."""
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    values = []
    for _ in range(BOOTSTRAP_REPS):
        x = rng.choice(diff, len(diff), replace=True)
        sd = x.std(ddof=1)
        if sd > 0:
            values.append(x.mean() / sd)
    if not values:
        raise ValueError("every bootstrap resample has zero SD: dz has no bootstrap band")
    return [float(v) for v in np.percentile(values, [10, 90])], BOOTSTRAP_REPS - len(values)


def diagnose(rows: pd.DataFrame, arm: str) -> dict:
    """Descriptive statistics of pearson_motif - pearson_signed on EPR, one arm."""
    a_name, b_name = PAIR
    mat = real_bench._model_matrix(rows, arm, "epr")[[a_name, b_name]]
    if not mat.index.is_unique:
        raise ValueError(f"{arm}: dataset ids are not unique: {sorted(set(mat.index[mat.index.duplicated()]))}")
    incomplete = mat.index[mat.isna().any(axis=1)].tolist()
    mat = mat.dropna()     # as the registered report: a dataset without both scores leaves this pair
    a, b = mat[a_name].to_numpy(float), mat[b_name].to_numpy(float)
    diff, n = a - b, len(mat)
    if n < 3:
        raise ValueError(f"{arm}: need >= 3 datasets with an EPR for both models")
    sd = float(diff.std(ddof=1))
    se = sd / math.sqrt(n)
    half = float(stats.t.ppf(1 - ALPHA / 2, n - 1)) * se
    loo = [_t(np.delete(diff, i)) for i in range(n)]
    shapiro = stats.shapiro(diff)
    auprc = real_bench._model_matrix(rows, arm, "auprc_ratio")[[a_name, b_name]].loc[mat.index]
    missing_auprc = auprc.index[auprc.isna().any(axis=1)].tolist()
    if missing_auprc:
        raise ValueError(f"{arm}: AUPRC ratio missing for {missing_auprc}, which have both EPR scores")
    auprc_diff = (auprc[a_name] - auprc[b_name]).to_numpy(float)
    positive = bool((a > 0).all() and (b > 0).all())   # the log ratio exists only then
    log_ratio = np.log(a / b) if positive else None
    by_task = rows.set_index(rows["task"].str.replace("cxg__", "", regex=False))
    ref_edges = by_task["reference_edges"]
    chi = stats.chi2.ppf([1 - ALPHA / 2, ALPHA / 2], n - 1)
    contexts = by_task.loc[mat.index, f"{arm}_n_contexts"].value_counts().sort_index()
    boot, zero_sd = _bootstrap_dz(diff)
    return {
        "n": n, "datasets_without_both_scores": incomplete,
        "mean_a": float(a.mean()), "mean_b": float(b.mean()),
        "sd_a": float(a.std(ddof=1)), "sd_b": float(b.std(ddof=1)),
        "corr_a_b": float(np.corrcoef(a, b)[0, 1]),
        "mean_diff": float(diff.mean()), "sd_diff": sd,
        "sd_diff_ci95": [float(math.sqrt((n - 1) * sd ** 2 / c)) for c in chi],
        "se": se, "t": _t(diff), "df": n - 1,
        "p_one_sided": float(stats.t.sf(_t(diff), n - 1)),
        "ci95": [float(diff.mean() - half), float(diff.mean() + half)],
        "dz": _dz(diff), "dz_bootstrap_10_90": boot, "dz_bootstrap_zero_sd_resamples": zero_sd,
        "median_diff": float(np.median(diff)), "min_diff": float(diff.min()), "max_diff": float(diff.max()),
        "n_positive": int((diff > 0).sum()),
        "skew": float(stats.skew(diff)), "excess_kurtosis": float(stats.kurtosis(diff)),
        "shapiro_w": float(shapiro.statistic), "shapiro_p": float(shapiro.pvalue),
        "loo_t_range": [float(min(loo)), float(max(loo))],
        "loo_t_extremes_drop": [mat.index[int(np.argmin(loo))], mat.index[int(np.argmax(loo))]],
        "spearman_diff_reference_edges": float(stats.spearmanr(diff, ref_edges.loc[mat.index]).statistic),
        "n_contexts": {str(k): int(v) for k, v in contexts.items()},
        "auprc_ratio": {"mean_diff": float(auprc_diff.mean()), "sd_diff": float(auprc_diff.std(ddof=1)),
                        "dz": _dz(auprc_diff)},
        "log_ratio": None if log_ratio is None else {
            "mean": float(log_ratio.mean()), "sd": float(log_ratio.std(ddof=1)), "dz": _dz(log_ratio)},
        "per_dataset": [{"dataset": d, "epr_a": float(x), "epr_b": float(y), "diff": float(x - y),
                         "reference_edges": int(ref_edges.loc[d])}
                        for d, x, y in zip(mat.index, a, b)],
    }


def _n_or_reason(dz: float, alpha: float, power: float) -> int | str:
    """n for a data-derived dz, or why there is none (a reason, not a crash of the report)."""
    if dz <= 0:
        return "none (dz <= 0)"
    try:
        return n_for_power(dz, alpha, power)
    except SampleSizeNotReached as exc:
        return f"not reached: {exc}"


def planning(diag: dict) -> dict:
    """Datasets a new study with the same paired t needs; effects in EPR units use this arm's SD."""
    sd, (sd_lo, sd_hi) = diag["sd_diff"], diag["sd_diff_ci95"]
    lo, hi = diag["dz_bootstrap_10_90"]
    effects = {}
    for effect in EFFECTS:
        effects[f"{effect:.2f}"] = {level: {
            "dz": effect / sd, "n80": n_for_power(effect / sd, a, .8), "n90": n_for_power(effect / sd, a, .9),
            "n80_sd_ci95": [n_for_power(effect / sd_lo, a, .8), n_for_power(effect / sd_hi, a, .8)]}
            for level, a in LEVELS.items()}
    n80 = effects["0.10"]["alpha/2"]["n80"]
    return {
        "effects": effects,
        "observed": {level: {"dz": diag["dz"], "n80": _n_or_reason(diag["dz"], a, .8),
                             "n90": _n_or_reason(diag["dz"], a, .9),
                             "n80_dz_bootstrap_10_90": [_n_or_reason(lo, a, .8), _n_or_reason(hi, a, .8)]}
                     for level, a in LEVELS.items()},
        "wilcoxon_n80_at_0.10_alpha/2": {"t": n80, "normal_ARE": math.ceil(n80 / WILCOXON_ARE_NORMAL),
                                         "minimum_ARE": math.ceil(n80 / WILCOXON_ARE_MIN)},
    }


def sensitivity(diag: dict, n: int) -> dict:
    out = {}
    for level, a in LEVELS.items():
        d80, d90 = min_detectable_dz(n, a, .8), min_detectable_dz(n, a, .9)
        out[level] = {"mde80_dz": d80, "mde80_epr": d80 * diag["sd_diff"],
                      "mde90_dz": d90, "mde90_epr": d90 * diag["sd_diff"]}
    return out


def grid(diag: dict) -> dict:
    """Power for +0.10 EPR and the smallest detectable effect by n; never power at C4's estimate."""
    sd = diag["sd_diff"]
    return {level: [{"n": n, "power_at_0.10": paired_t_power(0.10 / sd, n, LEVELS[level]),
                     "mde80_dz": min_detectable_dz(n, LEVELS[level], .8),
                     "mde80_epr": min_detectable_dz(n, LEVELS[level], .8) * sd,
                     "mde90_epr": min_detectable_dz(n, LEVELS[level], .9) * sd} for n in GRID_N]
            for level in GRID_LEVELS}


def analyse(rows: pd.DataFrame, report_text: str) -> dict:
    lines = registered_epr_lines(rows)
    check_against_report(lines, report_text)
    arms = {arm: diagnose(rows, arm) for arm in ARMS}
    for arm, d in arms.items():   # the diagnosed pair is the one the gate reprinted
        prefix = f"| {PAIR[0]} | {PAIR[1]} | {d['n']} | {d['mean_diff']:+.3f} |"
        if not any(line.startswith(prefix) for line in lines):
            raise GateError(f"{arm}: the diagnosed pair {prefix} is not among the gate lines")
    return {"experiment": "EXP-20261020-GB-1", "pair": list(PAIR), "metric": "epr",
            "gate": {"lines": lines, "matches_report": True},
            "arms": arms,
            "leave_one_out": leave_one_out(rows, lines),
            "sensitivity_n20": {arm: sensitivity(d, d["n"]) for arm, d in arms.items()},
            "planning": {arm: planning(d) for arm, d in arms.items()},
            "grid": {arm: grid(d) for arm, d in arms.items()},
            "settings": {"alpha": ALPHA, "levels": LEVELS, "effects_epr": list(EFFECTS), "grid_n": list(GRID_N),
                         "bootstrap": {"reps": BOOTSTRAP_REPS, "seed": BOOTSTRAP_SEED, "percentiles": [10, 90]},
                         "wilcoxon_are": {"normal": WILCOXON_ARE_NORMAL, "minimum": WILCOXON_ARE_MIN}}}


def _f(x, digits: int = 3) -> str:
    return x if isinstance(x, str) else f"{x:.{digits}f}"


def render(result: dict) -> str:
    a_name, b_name = result["pair"]
    L = [f"# {result['experiment']}: EPR of {a_name} vs {b_name} in C4, re-read from the committed rows", "",
         "Generated by `python benchmarks/epr_power.py`. Exploratory: descriptive diagnosis and planning numbers.",
         "C4's registered verdict is unchanged: the EPR comparison is not resolved, and not resolved is not zero.",
         "No number below is a new test of H1. No power is computed from C4's estimate: the achieved n is described",
         "by the smallest detectable effect, and C4's estimate enters only as the datasets a new study would need.", "",
         "## Inputs", ""]
    L += [f"- `{k}`: sha256 `{v}` (bytes as committed; matches that run's `SHA256SUMS`)"
          for k, v in result["inputs"].items()]
    L += ["", "## Gate: the registered code reprints the committed EPR comparisons", "",
          "The lines below are exactly the EPR tables of the committed `REPORT_REAL.md`, none more, none less.", ""]
    for line in result["gate"]["lines"]:
        L += (["", line, ""] if line.startswith("####") else [line])
    L += ["", "## Per-dataset differences (descriptive)", "",
          "| statistic | " + " | ".join(result["arms"]) + " |", "|---|" + "---|" * len(result["arms"])]
    A = result["arms"]
    def row(label, fn):
        L.append(f"| {label} | " + " | ".join(fn(A[arm]) for arm in A) + " |")
    row("n datasets", lambda d: str(d["n"]))
    row("datasets without both EPR scores (left out of this pair)",
        lambda d: ", ".join(d["datasets_without_both_scores"]) or "none")
    row(f"mean EPR {a_name} / {b_name}", lambda d: f"{_f(d['mean_a'])} / {_f(d['mean_b'])}")
    row("correlation of the two models' EPR", lambda d: _f(d["corr_a_b"]))
    row("mean difference [95% interval]", lambda d: f"{d['mean_diff']:+.3f} [{d['ci95'][0]:+.3f}, {d['ci95'][1]:+.3f}]")
    row("SD of differences [95% interval, chi-square]",
        lambda d: f"{_f(d['sd_diff'])} [{_f(d['sd_diff_ci95'][0])}, {_f(d['sd_diff_ci95'][1])}]")
    row("dz = mean / SD [bootstrap 10-90%]",
        lambda d: f"{_f(d['dz'], 2)} [{_f(d['dz_bootstrap_10_90'][0], 2)}, {_f(d['dz_bootstrap_10_90'][1], 2)}]")
    row("median, min, max", lambda d: f"{d['median_diff']:+.3f}, {d['min_diff']:+.3f}, {d['max_diff']:+.3f}")
    row("datasets with a positive difference", lambda d: f"{d['n_positive']}/{d['n']}")
    row("skew, excess kurtosis", lambda d: f"{_f(d['skew'], 2)}, {_f(d['excess_kurtosis'], 2)}")
    row("Shapiro-Wilk W (p)", lambda d: f"{_f(d['shapiro_w'])} ({_f(d['shapiro_p'], 2)})")
    row("t leaving one dataset out: min, max", lambda d: f"{_f(d['loo_t_range'][0], 2)}, {_f(d['loo_t_range'][1], 2)}")
    row("Spearman, difference vs reference_edges (not the EPR's k)",
        lambda d: _f(d["spearman_diff_reference_edges"], 2))
    row("same pair on AUPRC ratio: mean, SD, dz", lambda d: ", ".join(
        _f(d["auprc_ratio"][k], 3 if k != "dz" else 2) for k in ("mean_diff", "sd_diff", "dz")))
    row("log(EPR ratio): mean, SD, dz (scale diagnostic only)", lambda d: "not computed (an EPR is 0)"
        if d["log_ratio"] is None
        else ", ".join(_f(d["log_ratio"][k], 3 if k != "dz" else 2) for k in ("mean", "sd", "dz")))
    row("contexts per dataset (count)", lambda d: ", ".join(f"{k}: {v}" for k, v in d["n_contexts"].items()))
    L += ["", "## Leaving one dataset out with the registered report code", "",
          "`real_bench.models_section` (C4's procedure: paired t, Holm over the 4 directed hypotheses of each",
          "metric and arm) rerun once per dataset on all the other datasets. An influence check of C4's own test,",
          "not a new test.", ""]
    for arm, o in result["leave_one_out"].items():
        L.append(f"- {arm}: verdict with all {o['drops']} datasets `{o['verdict_all']}`; it changes in "
                 f"{len(o['changed'])} of {o['drops']} single-dataset drops" + (":" if o["changed"] else "."))
        L += [f"  - without `{c['dropped']}`: `{c['line']}`" for c in o["changed"]]
    L += ["", f"## Per dataset (EPR {a_name} - {b_name})", "",
          "| dataset | reference_edges | " + " | ".join(f"{arm} diff" for arm in A) + " |",
          "|---|---|" + "---|" * len(A)]
    by_arm = {arm: {r["dataset"]: r for r in A[arm]["per_dataset"]} for arm in A}
    first = next(iter(A))
    for r in sorted(A[first]["per_dataset"], key=lambda r: -r["diff"]):
        L.append(f"| {r['dataset']} | {r['reference_edges']} | " + " | ".join(
            f"{by_arm[arm][r['dataset']]['diff']:+.3f}" if r["dataset"] in by_arm[arm] else "-" for arm in A) + " |")
    L += ["", "## What n = 20 could detect (sensitivity, one-sided paired t)", "",
          "| arm | level | smallest dz at 80% | in EPR units | smallest dz at 90% | in EPR units |", "|---|---|---|---|---|---|"]
    for arm, by_level in result["sensitivity_n20"].items():
        for level, s in by_level.items():
            L.append(f"| {arm} | {level} | {_f(s['mde80_dz'], 2)} | {s['mde80_epr']:+.3f} | "
                     f"{_f(s['mde90_dz'], 2)} | {s['mde90_epr']:+.3f} |")
    L += ["", "## Datasets a new study would need (planning; never-used datasets, same paired t)", "",
          "Effects in EPR units are converted with each arm's SD of differences; the last column uses",
          "the SD's 95% interval instead (n grows with the square of the SD).", "",
          "| arm | EPR effect | level | dz | n for 80% | n for 90% | n for 80% at SD interval |", "|---|---|---|---|---|---|---|"]
    for arm, plan in result["planning"].items():
        for effect, by_level in plan["effects"].items():
            for level, e in by_level.items():
                L.append(f"| {arm} | +{effect} | {level} | {_f(e['dz'], 2)} | {e['n80']} | {e['n90']} | "
                         f"{e['n80_sd_ci95'][0]}-{e['n80_sd_ci95'][1]} |")
    L += ["", "If the true effect equalled C4's estimate (a scenario, with the bootstrap band of dz):", "",
          "| arm | level | dz | n for 80% | n for 90% | n for 80% at dz bootstrap 10% / 90% |", "|---|---|---|---|---|---|"]
    for arm, plan in result["planning"].items():
        for level, o in plan["observed"].items():
            lo, hi = o["n80_dz_bootstrap_10_90"]
            L.append(f"| {arm} | {level} | {_f(o['dz'], 2)} | {o['n80']} | {o['n90']} | {lo} / {hi} |")
    L += ["", "## Power for +0.10 EPR and smallest detectable EPR effect by n", "",
          "| arm | level | n | power at +0.10 | smallest EPR effect, 80% | 90% |",
          "|---|---|---|---|---|---|"]
    for arm, by_level in result["grid"].items():
        for level, rows_ in by_level.items():
            for g in rows_:
                L.append(f"| {arm} | {level} | {g['n']} | {_f(g['power_at_0.10'])} | "
                         f"{g['mde80_epr']:+.3f} | {g['mde90_epr']:+.3f} |")
    w = result["planning"]["auto"]["wilcoxon_n80_at_0.10_alpha/2"]
    L += ["", "## Notes", "",
          f"- Levels: alpha/4 = first Holm step of C4's family of 4 directed hypotheses; alpha/3 = the same family "
          f"once the other pair is rejected first; alpha/2 = a family of one pair (two directions).",
          f"- A signed-rank test instead of the t (not the registered test) would need about {w['normal_ARE']} "
          f"(normal differences) and at most {w['minimum_ARE']} datasets where the t needs {w['t']} "
          "(auto, +0.10, alpha/2; asymptotic relative efficiency 3/pi and its minimum 0.864).",
          "- `reference_edges` counts CollecTRI pairs of the scored TFs over all genes in both halves; it is not "
          "the k of the EPR. k, prevalence, the 1/prevalence ceiling and ties at k need the kept artifacts (step 15-2).",
          "- Figure: `power_epr.png` (left: per-dataset differences with the 95% interval; right: power by n "
          "for +0.10 EPR)."]
    return "\n".join(L) + "\n"


def power_curves(result: dict) -> list[dict]:
    """Right panel of the figure: power by n for +0.10 EPR (auto arm SD); never for C4's estimate."""
    dz = 0.10 / result["arms"]["auto"]["sd_diff"]
    return [{"label": f"+0.10, {level}", "level": level, "dz": dz, "n": CURVE_N,
             "power": [paired_t_power(dz, int(n), LEVELS[level]) for n in CURVE_N]} for level in GRID_LEVELS]


def plot(result: dict, out_png: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from grnbench import viz
    viz._use_style()
    fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4))
    for i, (arm, d) in enumerate(result["arms"].items()):
        diffs = [r["diff"] for r in d["per_dataset"]]
        jitter = np.linspace(-0.15, 0.15, len(diffs))
        left.scatter(np.full(len(diffs), i) + jitter, diffs, s=14, alpha=0.7)
        left.errorbar(i + 0.3, d["mean_diff"], yerr=[[d["mean_diff"] - d["ci95"][0]], [d["ci95"][1] - d["mean_diff"]]],
                      fmt="o", color="black", capsize=3)
    left.axhline(0, color="grey", lw=0.8)
    left.set_xticks(range(len(result["arms"])), list(result["arms"]))
    left.set_ylabel(f"EPR {result['pair'][0]} - {result['pair'][1]}")
    left.set_title(f"per dataset (n = {result['arms']['auto']['n']})")
    for curve in power_curves(result):
        right.plot(curve["n"], curve["power"], label=curve["label"])
    right.axhline(0.8, color="grey", lw=0.8)
    right.set_xlabel("datasets (n)")
    right.set_ylabel("power (auto arm SD)")
    right.legend(fontsize=7)
    fig.tight_layout()
    viz._save(fig, out_png)


def main(argv: list[str] | None = None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--rows", type=Path, default=C4_DIR / "all_rows_real.csv")
    ap.add_argument("--report", type=Path, default=C4_DIR / "REPORT_REAL.md")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    sums = verify_inputs(args.rows, args.report)
    result = analyse(load_rows(args.rows), args.report.read_text(encoding="utf-8"))
    def rel(p: Path) -> str:
        try:
            return p.resolve().relative_to(ROOT).as_posix()
        except ValueError:
            return p.name
    result["inputs"] = {rel(args.rows): sums[args.rows.name], rel(args.report): sums[args.report.name]}
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "power.json").write_text(json.dumps(result, indent=1, allow_nan=False) + "\n",
                                          encoding="utf-8", newline="\n")
    (args.out / "POWER.md").write_text(render(result), encoding="utf-8", newline="\n")
    plot(result, args.out / "power_epr.png")
    for arm, d in result["arms"].items():
        print(f"{arm}: mean diff {d['mean_diff']:+.3f}, SD {d['sd_diff']:.3f}, dz {d['dz']:.2f}; "
              f"n for 80% at +0.10 (alpha/2): {result['planning'][arm]['effects']['0.10']['alpha/2']['n80']}")
    for arm, o in result["leave_one_out"].items():
        print(f"{arm}: leaving one dataset out changes the verdict in {len(o['changed'])} of {o['drops']} drops"
              + "".join(f"; without {c['dropped']}: {c['line']}" for c in o["changed"]))
    return result


if __name__ == "__main__":
    main()
