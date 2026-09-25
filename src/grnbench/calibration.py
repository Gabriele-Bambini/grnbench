"""Score-level calibration pilot for the established paired-t/Holm reference.

Y[u,m] = mu[m] + shared_sd*Z[u] + method_sd*E[u,m]. Independent normal
Z/E imply known mean differences and exact t assumptions. Repeating each unit
several times deliberately violates independence: it must not be called a GRN
simulation or evidence that biological benchmark comparisons are calibrated.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta

from grnbench.eval.comparison import compare_units


def binomial_bounds(events: int, trials: int, confidence: float = .95) -> dict:
    """Separate one-sided exact Clopper-Pearson bounds (not a two-sided 95% CI)."""
    if not isinstance(events, (int, np.integer)) or not isinstance(trials, (int, np.integer)):
        raise ValueError("events and trials must be integers")
    if trials < 1 or not 0 <= events <= trials or not .5 < confidence < 1:
        raise ValueError("require 0 <= events <= trials, trials > 0 and .5 < confidence < 1")
    tail = 1 - confidence
    return {"events": int(events), "trials": int(trials), "rate": events / trials,
            "lower": 0. if events == 0 else float(beta.ppf(tail, events, trials - events + 1)),
            "upper": 1. if events == trials else float(beta.ppf(confidence, events + 1, trials - events)),
            "one_sided_confidence": confidence}


def calibration_gate(false_superiority: int, covered: int, trials: int,
                     alpha: float = .05, tolerance: float = .025,
                     confidence: float = .95) -> dict:
    if not 0 < alpha < 1 or not 0 <= tolerance < 1 - alpha:
        raise ValueError("invalid alpha/tolerance")
    err = binomial_bounds(false_superiority, trials, confidence)
    cov = binomial_bounds(covered, trials, confidence)
    ceiling, floor = alpha + tolerance, 1 - alpha - tolerance
    if err["lower"] > ceiling or cov["upper"] < floor:
        status = "fail"
    elif err["upper"] <= ceiling and cov["lower"] >= floor:
        status = "pass"
    else:
        status = "inconclusive"
    return {"status": status, "familywise_false_superiority": err, "simultaneous_coverage": cov,
            "error_ceiling": ceiling, "coverage_floor": floor,
            "interpretation": "scenario-specific Monte Carlo gate; not a global validation certificate"}


def _integer(value, name, minimum):
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


def _validate_config(config):
    c = json.loads(json.dumps(config, allow_nan=False))
    allowed = {"seed", "replicates", "alpha", "margin", "tolerance", "mc_confidence", "scenarios"}
    if set(c) - allowed:
        raise ValueError(f"unknown config fields: {sorted(set(c) - allowed)}")
    c.setdefault("seed", 42)
    c.setdefault("replicates", 200)
    c.setdefault("alpha", .05)
    c.setdefault("margin", 0.)
    c.setdefault("tolerance", .025)
    c.setdefault("mc_confidence", .95)
    _integer(c["seed"], "seed", 0)
    _integer(c["replicates"], "replicates", 1)
    if not 0 < c["alpha"] < 1 or not 0 <= c["tolerance"] < 1 - c["alpha"]:
        raise ValueError("invalid alpha or tolerance")
    if not .5 < c["mc_confidence"] < 1 or c["margin"] < 0:
        raise ValueError("invalid confidence or margin")
    if not isinstance(c.get("scenarios"), list) or not c["scenarios"]:
        raise ValueError("at least one scenario is required")
    names = set()
    for s in c["scenarios"]:
        if set(s) - {"name", "n_units", "means", "shared_noise", "independent_noise", "repeat_units"}:
            raise ValueError("unknown scenario fields")
        if not isinstance(s.get("name"), str) or not s["name"].strip() or s["name"] in names:
            raise ValueError("scenario names must be nonempty and unique")
        names.add(s["name"])
        n = _integer(s.get("n_units"), "n_units", 3)
        repeat = _integer(s.setdefault("repeat_units", 1), "repeat_units", 1)
        if n % repeat or repeat > n:
            raise ValueError("n_units must be divisible by repeat_units")
        means = np.asarray(s.get("means"), dtype=float)
        if means.ndim != 1 or len(means) < 2 or not np.isfinite(means).all():
            raise ValueError("means must be a finite vector with at least two methods")
        if s.setdefault("shared_noise", .5) < 0 or s.setdefault("independent_noise", 1.) <= 0:
            raise ValueError("shared noise must be nonnegative and independent noise positive")
    return c


def run_calibration(config: dict, out: str | Path) -> dict:
    c = _validate_config(config)
    out = Path(out)
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise FileExistsError("output must be absent or empty; use a fresh run directory")
    try:
        return _run_calibration(c, out)
    except Exception as exc:
        manifest_path = out / "manifest.json"
        if manifest_path.exists():
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest.update(status="failed", error=f"{type(exc).__name__}: {exc}")
            _write_json(manifest_path, manifest)
        raise


def _run_calibration(c: dict, out: Path) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    encoded = json.dumps(c, sort_keys=True, separators=(",", ":"))
    from grnbench.report import environment_block
    manifest = {"config": c, "config_sha256": hashlib.sha256(encoded.encode()).hexdigest(),
                "environment": environment_block(), "status": "running",
                "source_sha256": {str(p.name): hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in [Path(__file__), Path(__file__).parent / "eval/comparison.py"]}}
    _write_json(out / "manifest.json", manifest)
    rows, summaries = [], []
    streams = np.random.SeedSequence(c["seed"]).spawn(len(c["scenarios"]))
    for scenario, stream in zip(c["scenarios"], streams):
        rng = np.random.default_rng(stream)
        mu = np.asarray(scenario["means"], dtype=float)
        names = [f"method_{i}" for i in range(len(mu))]
        n, repeat = scenario["n_units"], scenario["repeat_units"]
        independent_n = n // repeat
        true_deltas = {(names[i], names[j]): float(mu[i] - mu[j])
                       for i in range(len(mu)) for j in range(i + 1, len(mu))}
        for rep in range(c["replicates"]):
            observed = (mu + scenario["shared_noise"] * rng.normal(size=(independent_n, 1))
                        + scenario["independent_noise"] * rng.normal(size=(independent_n, len(mu))))
            observed = np.repeat(observed, repeat, axis=0)
            inference = compare_units(pd.DataFrame(observed, columns=names), c["alpha"], c["margin"])
            false_call, covers, correct_calls, alternatives, calls = False, True, 0, 0, 0
            for p in inference["pairs"]:
                delta = true_deltas[p["a"], p["b"]]
                lo, hi = p["simultaneous_ci"]
                covers &= lo is None or lo <= delta <= hi
                has_alternative = abs(delta) > c["margin"]
                alternatives += has_alternative
                if p["winner"]:
                    calls += 1
                    correct = (delta > c["margin"] if p["winner"] == p["a"]
                               else -delta > c["margin"])
                    correct_calls += correct
                    false_call |= not correct
            rows.append({"scenario": scenario["name"], "replicate": rep,
                         "false_superiority": int(false_call), "covers_all": int(covers),
                         "power": correct_calls / alternatives if alternatives else None,
                         "unresolved_fraction": 1 - calls / len(inference["pairs"])})
        subset = rows[-c["replicates"]:]
        gate = calibration_gate(sum(x["false_superiority"] for x in subset),
                                sum(x["covers_all"] for x in subset), c["replicates"],
                                c["alpha"], c["tolerance"], c["mc_confidence"])
        powers = [x["power"] for x in subset if x["power"] is not None]
        summaries.append({**scenario, "n_independent_units": independent_n,
                          "independence_assumption_holds": repeat == 1,
                          "truth": "known Gaussian score means, not estimated GRN accuracy",
                          "gate": gate, "mean_power": float(np.mean(powers)) if powers else None,
                          "mean_unresolved_fraction": float(np.mean([x["unresolved_fraction"] for x in subset]))})
    result = {"scope": "score_level_gaussian_pilot", "biological_validation": False,
              "novelty_established": False, "procedure": "paired_t_holm",
              "config_sha256": manifest["config_sha256"], "scenarios": summaries,
              "limitations": ["No real expression matrix, GRN simulator, or incomplete-label model is tested.",
                              "Means are known exactly in this score-level model; this does not hold for real algorithm rankings.",
                              "Each gate is scenario-specific; no simultaneous confidence across gates is claimed.",
                              "Repeated-unit scenarios intentionally violate independence.",
                              "Power is averaged across independent experiments; dependent method pairs are not counted as independent trials."]}
    pd.DataFrame(rows).to_csv(out / "experiments.csv", index=False)
    _write_json(out / "calibration.json", result)
    lines = ["# Score-level calibration pilot", "", "This is not biological validation or evidence of novelty.", "",
             "| scenario | independent units / rows | FWER (upper bound) | coverage (lower bound) | power | gate |",
             "|---|---:|---:|---:|---:|---|"]
    for s in summaries:
        e, v = s["gate"]["familywise_false_superiority"], s["gate"]["simultaneous_coverage"]
        power = "n/a (no alternative)" if s["mean_power"] is None else f"{s['mean_power']:.3f}"
        name = s["name"].replace("|", "\\|").replace("\n", " ")
        lines.append(f"| {name} | {s['n_independent_units']} / {s['n_units']} | {e['rate']:.3f} ({e['upper']:.3f}) | {v['rate']:.3f} ({v['lower']:.3f}) | {power} | {s['gate']['status']} |")
    lines.extend(["", f"Bounds: separate one-sided exact {c['mc_confidence']:.0%} Monte Carlo limits.",
                  f"Gate: error upper <= {c['alpha'] + c['tolerance']:.3f} and coverage lower >= {1-c['alpha']-c['tolerance']:.3f}.",
                  "Finite precision may yield an inconclusive gate. A failure in repeated-unit controls is expected.", ""])
    lines.extend("- " + x for x in result["limitations"])
    (out / "calibration.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    manifest["status"] = "completed"
    _write_json(out / "manifest.json", manifest)
    return result


def _write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
