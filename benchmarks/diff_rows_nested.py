"""Where do two rows/ folders differ, at any depth, timings apart? Exploratory diagnosis, not a verdict.

    python benchmarks/diff_rows_nested.py COMMITTED_ROWS RERUN_ROWS [--rtol 0] [--show 20]

``compare_runs.py`` skips timing fields (``seconds``, ``*_seconds``) only at the top level of a row
and compares a list inside a row by plain equality. A ``candidate_study`` row keeps one entry per
candidate in ``candidates``, each with its own ``seconds``, so two runs with identical scores show
every row as DIFFERENT there (first seen in the EXP-20261003-LM-1 rerun, issue #95). This script
walks both rows in parallel through dicts and lists and reports:

- how many rows each side has and which are missing;
- how many leaf values were compared, and how many timing values were skipped (and differ);
- every non-timing difference (numeric or other), with the largest relative numeric difference.

Two NaN values count as equal. A key or list element present on one side only is an "other"
difference. Numbers are compared exactly unless ``--rtol`` is given. Exit 0 when there is no
non-timing difference and no missing row, 1 otherwise, 2 when a folder has no rows or a row is
unreadable. It describes a rerun; the verdict of a reproduction is the rule of its issue.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

MISSING = "<missing>"


def _timing(key: str) -> bool:
    return key == "seconds" or key.endswith("_seconds")


def _number(v) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool)


class _Diff:
    def __init__(self, rtol: float):
        self.rtol = rtol
        self.leaves = self.timing = self.timing_differ = 0
        self.numeric: list[tuple[str, str, object, object, float]] = []
        self.other: list[tuple[str, str, object, object]] = []

    def walk(self, task: str, path: str, a, b, timing: bool = False) -> None:
        if isinstance(a, dict) and isinstance(b, dict):
            for k in list(a) + [k for k in b if k not in a]:
                sub = f"{path}.{k}" if path else k
                if k not in a or k not in b:
                    self.other.append((task, sub, a.get(k, MISSING), b.get(k, MISSING)))
                else:
                    self.walk(task, sub, a[k], b[k], timing or _timing(k))
        elif isinstance(a, list) and isinstance(b, list):
            for i in range(max(len(a), len(b))):
                sub = f"{path}[{i}]"
                if i >= len(a) or i >= len(b):
                    self.other.append((task, sub, a[i] if i < len(a) else MISSING, b[i] if i < len(b) else MISSING))
                else:
                    self.walk(task, sub, a[i], b[i], timing)
        elif timing:
            self.timing += 1
            self.timing_differ += a != b
        else:
            self.leaves += 1
            if _number(a) and _number(b):
                if math.isnan(a) and math.isnan(b):
                    return
                same = a == b if self.rtol == 0 else math.isclose(a, b, rel_tol=self.rtol, abs_tol=0.0)
                if not same:
                    rel = abs(a - b) / max(abs(a), abs(b)) if max(abs(a), abs(b)) > 0 else math.inf
                    self.numeric.append((task, path, a, b, rel))
            elif a != b:
                self.other.append((task, path, a, b))


def _rows(folder: Path) -> dict[str, dict]:
    out = {}
    for f in sorted(folder.glob("*.json")):
        try:
            out[f.stem] = json.loads(f.read_text(encoding="utf-8"))
        except (OSError, ValueError) as e:
            raise ValueError(f"cannot read {f.name}: {e}") from e
    if not out:
        raise ValueError(f"no rows in {folder}")
    return out


def _short(v, n: int = 120) -> str:
    s = v if v == MISSING else repr(v)
    return s if len(s) <= n else s[: n - 3] + "..."


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Diff two rows/ folders at any depth, timing fields apart.")
    p.add_argument("committed", type=Path)
    p.add_argument("rerun", type=Path)
    p.add_argument("--rtol", type=float, default=0.0)
    p.add_argument("--show", type=int, default=20, help="differences listed (default 20)")
    args = p.parse_args(argv)
    try:
        ca, ra = _rows(args.committed), _rows(args.rerun)
    except ValueError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    d = _Diff(args.rtol)
    missing = [t for t in ca if t not in ra]
    extra = [t for t in ra if t not in ca]
    for t in ca:
        if t in ra:
            d.walk(t, "", ca[t], ra[t])
    line = f"rows: {len(ca)} committed, {len(ra)} rerun, {len(missing)} missing in rerun"
    print(line + (f", {len(extra)} extra in rerun" if extra else ""))
    print(f"leaf values compared: {d.leaves}; timing values skipped: {d.timing} ({d.timing_differ} differ); "
          + ("exact" if args.rtol == 0 else f"tolerance rtol={args.rtol:g} atol=0"))
    n = len(d.numeric) + len(d.other)
    maxrel = max((x[4] for x in d.numeric), default=0.0)
    print(f"non-timing differences: {n} (numeric {len(d.numeric)}, other {len(d.other)}); "
          f"max relative numeric difference {maxrel:.3g}")
    shown = sorted(d.numeric, key=lambda x: -x[4])[: args.show] + d.other[: max(0, args.show - len(d.numeric))]
    for task, path, a, b, *_ in shown:
        print(f"  {task}\t{path}\t{_short(a)} -> {_short(b)}")
    for t in missing:
        print(f"  {t}\t(row)\tpresent -> {MISSING}")
    return 0 if n == 0 and not missing else 1


if __name__ == "__main__":
    sys.exit(main())
