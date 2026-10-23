"""Compare the confirmatory lines of a candidate-study rerun with the committed CONFIRM.json.

    python benchmarks/compare_confirm.py benchmarks/results/EXP-20261003-LM-1/CONFIRM.json \
        RERUN/CONFIRM.json --record experiments/EXP-20261003-LM-1.md [--rtol 1e-9]

``compare_runs.py`` checks the task rows of a ``candidate_study`` run (``--rows-only``) but cannot
read its primary: the registered test lives in ``CONFIRM.json`` (``candidate_study.py confirm``).
This script reads the record's ``primary comparison: A > B | metric | condition`` lines and, for
each, compares the pair ``(a=A, b=B)`` of that condition field by field: ``n``, ``excluded``,
``diff``, both ``ci`` bounds, ``p_a_gt_b``, ``p_b_gt_a`` (Holm-adjusted) and the decision
``supported``. Numbers are compared exactly unless ``--rtol`` is given; text (the decision) is
always compared exactly. The pairs of the other conditions are printed as information only.

Prints one ``primary:`` line per registered primary, then ``verdict: MATCH`` (exit 0: every
primary matches) or ``verdict: DIFFERENT`` (exit 1). A record with no parsable primary, a primary
pair missing from either file and an unreadable file are errors (exit 2), never a pass.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

_PRIMARY_LINE = re.compile(r"^- primary comparison:\s*(.+)$", re.MULTILINE)
_SPEC = re.compile(r"^\s*`?(\w+)\s*>\s*(\w+)\s*\|\s*(\w+)\s*\|\s*(\w+)`?\s*$")


class CompareError(Exception):
    pass


def read_primaries(record: Path) -> list[tuple[str, str, str, str]]:
    """(A, B, metric, condition) of every ``A > B | metric | condition`` in the record's primary line."""
    m = _PRIMARY_LINE.search(record.read_text(encoding="utf-8"))
    if not m:
        raise CompareError(f"{record}: no '- primary comparison:' line")
    out = []
    for spec in m.group(1).split(";"):
        s = _SPEC.match(spec)
        if not s:
            raise CompareError(f"unparsable primary {spec.strip()!r} (expected 'A > B | metric | condition')")
        out.append(s.groups())
    return out


def _load(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as e:
        raise CompareError(f"cannot read {path}: {e}") from e


def _pair(confirm: dict, a: str, b: str, condition: str) -> dict | None:
    return next((p for p in confirm.get(condition, []) if p.get("a") == a and p.get("b") == b), None)


def _same(x, y, rtol: float) -> bool:
    if isinstance(x, bool) or isinstance(y, bool) or not all(isinstance(v, (int, float)) for v in (x, y)):
        return x == y
    if math.isnan(x) and math.isnan(y):
        return True
    return x == y if rtol == 0 else math.isclose(x, y, rel_tol=rtol, abs_tol=0.0)


def compare_pair(c: dict, r: dict, rtol: float) -> list[str]:
    """One line per differing field (``ci`` compared bound by bound)."""
    diffs = []
    for key in sorted(set(c) | set(r)):
        cv, rv = c.get(key), r.get(key)
        if isinstance(cv, list) and isinstance(rv, list) and len(cv) == len(rv):
            diffs += [f"{key}[{i}]: {x!r} -> {y!r}" for i, (x, y) in enumerate(zip(cv, rv))
                      if not _same(x, y, rtol)]
        elif not _same(cv, rv, rtol):
            diffs.append(f"{key}: {cv!r} -> {rv!r}")
    return diffs


def _fmt(p: dict) -> str:
    lo, hi = p["ci"]
    return (f"n={p['n']} excluded={p['excluded']} mean {p['diff']:+.3f} [{lo:+.3f}, {hi:+.3f}] "
            f"p_Holm(A>B) {p['p_a_gt_b']:.3g} p_Holm(B>A) {p['p_b_gt_a']:.3g} supported {p['supported']}")


def run(committed: Path, rerun: Path, record: Path, rtol: float) -> int:
    primaries = read_primaries(record)
    c_all, r_all = _load(committed), _load(rerun)
    tol = "exact" if rtol == 0 else f"tolerance rtol={rtol:g} atol=0"
    print(f"conditions: committed {sorted(c_all)}, rerun {sorted(r_all)}; {tol}")
    ok = True
    for a, b, metric, condition in primaries:
        c, r = _pair(c_all, a, b, condition), _pair(r_all, a, b, condition)
        if c is None or r is None:
            side = "committed" if c is None else "rerun"
            raise CompareError(f"primary {a} > {b} | {condition} missing from the {side} file")
        diffs = compare_pair(c, r, rtol)
        ok &= not diffs
        print(f"primary: {'DIFFERENT' if diffs else 'MATCH'} - {a} > {b} | {metric} | {condition}")
        print(f"  committed: {_fmt(c)}")
        print(f"  rerun:     {_fmt(r)}")
        for d in diffs:
            print(f"    {d}")
    keys = {(a, b, cond) for a, b, _, cond in primaries}
    for condition in sorted(set(c_all) | set(r_all)):
        for p in c_all.get(condition, []):
            if (p.get("a"), p.get("b"), condition) in keys:
                continue
            r = _pair(r_all, p.get("a"), p.get("b"), condition)
            state = "MISSING in rerun" if r is None else ("DIFFERENT" if compare_pair(p, r, rtol) else "MATCH")
            print(f"secondary (information only): {state} - {p.get('a')} > {p.get('b')} | {condition}")
    print(f"verdict: {'MATCH' if ok else 'DIFFERENT'}")
    return 0 if ok else 1


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Compare a candidate-study rerun's CONFIRM.json with the committed one.")
    p.add_argument("committed", type=Path, help="committed CONFIRM.json")
    p.add_argument("rerun", type=Path, help="the rerun's CONFIRM.json")
    p.add_argument("--record", type=Path, required=True,
                   help="experiment record holding 'primary comparison: A > B | metric | condition'")
    p.add_argument("--rtol", type=float, default=0.0)
    args = p.parse_args(argv)
    if not (math.isfinite(args.rtol) and args.rtol >= 0):
        print("error: --rtol must be finite and >= 0", file=sys.stderr)
        return 2
    try:
        return run(args.committed, args.rerun, args.record, args.rtol)
    except CompareError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
