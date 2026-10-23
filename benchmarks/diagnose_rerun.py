"""Where does a rerun differ from the committed rows? Exploratory diagnosis, not a verdict.

    python benchmarks/diagnose_rerun.py COMMITTED_ROWS RERUN_ROWS

COMMITTED_ROWS and RERUN_ROWS are two ``rows/`` folders of ``real_bench.py`` (one JSON file per
dataset). The verdict of a reproduction comes from ``benchmarks/compare_runs.py`` and the rule of
its issue; this script only describes the differences that tool reports, for the record's
diagnosis paragraph (first used for EXP-20260925-LM-2 section 4b, issue #24):

- how many datasets have at least one differing field;
- the relative size of every differing numeric field (median, p90, max);
- every differing non-numeric field (planner choices such as ``auto_candidate``), per dataset;
- the primary inputs (``m_<model>_auto_auprc_ratio`` of pearson_motif and pearson_signed): how many
  differ, the smallest and largest relative difference, and the six largest changes.

Timing fields (``seconds``, ``*_seconds``) are skipped, as in ``compare_runs.py``. Two NaN values
count as equal. A numeric field that is NaN, null or missing on one side only is listed with the
non-numeric fields, never folded into the relative differences. Percentiles are nearest-rank. The
relative-difference summary covers every numeric field, counters included (for example
``n_warnings``), not only scores. Exit code 2 when a folder has no rows or the rerun lacks a
committed dataset.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

PRIMARY_MODELS = ("pearson_motif", "pearson_signed")


def _timing(key: str) -> bool:
    return key == "seconds" or key.endswith("_seconds")


def _number(value) -> bool:
    """A finite int or float (not a bool, not NaN)."""
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def _nan(value) -> bool:
    return isinstance(value, float) and math.isnan(value)


def diagnose(committed: Path, rerun: Path) -> list[str]:
    ids = sorted(p.stem for p in committed.glob("*.json"))
    if not ids:
        raise SystemExit(f"no rows in {committed}")
    missing = [i for i in ids if not (rerun / f"{i}.json").exists()]
    if missing:
        print(f"rerun lacks {len(missing)} committed datasets: {missing[:3]}", file=sys.stderr)
        raise SystemExit(2)
    rel, planner, prim, prim_other, changed = [], {}, [], [], set()
    for i in ids:
        a = json.loads((committed / f"{i}.json").read_text(encoding="utf-8"))
        b = json.loads((rerun / f"{i}.json").read_text(encoding="utf-8"))
        for k in sorted(set(a) | set(b)):
            x, y = a.get(k), b.get(k)
            if _timing(k) or x == y or (_nan(x) and _nan(y)):
                continue
            changed.add(i)
            if _number(x) and _number(y):
                rel.append((abs(y - x) / max(abs(x), 1e-12), k, i))
            else:
                planner.setdefault(i, []).append(f"{k}: {x} -> {y}")
        for m in PRIMARY_MODELS:
            k = f"m_{m}_auto_auprc_ratio"
            x, y = a.get(k), b.get(k)
            if _number(x) and _number(y):
                prim.append((abs(y - x) / max(abs(x), 1e-12), m, i, x, y))
            elif not (x == y or (_nan(x) and _nan(y))):
                prim_other.append(f"  not comparable  {m:15s} {i}  {x} -> {y}")
    rel.sort()
    pct = lambda p: rel[min(len(rel) - 1, int(round(p * (len(rel) - 1))))][0]
    out = [f"datasets with at least one differing field: {len(changed)}/{len(ids)}"]
    if rel:
        top, top_key, top_ds = rel[-1]
        out.append(f"numeric fields differing: {len(rel)}; relative difference median {pct(0.5):.2e}, "
                   f"p90 {pct(0.9):.2e}, max {top:.2e} ({top_key} in {top_ds})")
    else:
        out.append("numeric fields differing: 0")
    out.append("non-numeric (planner) fields differing:")
    out += [f"  {i}  {s}" for i, items in planner.items() for s in items]
    prim.sort(reverse=True)
    nonzero = [p[0] for p in prim if p[0] > 0]
    head = (f"primary inputs (auto AUPRC ratio, {len(PRIMARY_MODELS)} models x {len(ids)} datasets): "
            f"{len(nonzero) + len(prim_other)} differ")
    if nonzero:
        head += f"; relative difference min {min(nonzero):.2e}, max {max(nonzero):.2e}"
    if prim_other:
        head += f"; {len(prim_other)} not comparable (NaN, null or missing on one side)"
    out.append(head)
    out += [f"  {p[0]:.2e}  {p[1]:15s} {p[2]}  {p[3]:.6f} -> {p[4]:.6f}" for p in prim[:6] if p[0] > 0]
    out += prim_other
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("committed", type=Path, help="committed rows/ folder")
    ap.add_argument("rerun", type=Path, help="rerun rows/ folder")
    args = ap.parse_args(argv)
    print("\n".join(diagnose(args.committed, args.rerun)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
