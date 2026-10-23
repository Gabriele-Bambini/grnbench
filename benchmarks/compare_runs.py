"""Compare a rerun of a registered real-data experiment with its committed results.

    python benchmarks/compare_runs.py benchmarks/results/EXP-20260925-LM-2 RERUN_DIR \
        --record experiments/EXP-20260925-LM-2.md [--rtol 0 --atol 0] [--json diff.json]

A reproduction is judged from files, by a script fixed before the rerun, not by eye:

1. tasks   - the same task ids (``rows/<id>.json``) on both sides;
2. rows    - every field of every row: same keys, same values. Numbers are compared exactly
             unless ``--rtol``/``--atol`` are given (another OS or BLAS can move the last digits;
             say so in the record if you use them). Timing fields (``seconds``, ``*_seconds``) are
             not compared and are counted apart, so "N numeric fields identical" has one meaning;
3. report  - ``REPORT_REAL.md`` as an aligned line diff (difflib), ignoring line endings (a
             Windows run writes CRLF) and the timing columns (``s``, ``s per dataset ...``).
             Information only, not part of the verdict: the report is generated from the rows,
             so with identical rows a report difference comes from the code that wrote it (a newer
             report format, say), never from the data. To compare reports, regenerate both from
             their rows with the same code;
4. primary - the record's ``primary comparison: A > B | metric | arm``: the per-dataset values
             that enter the registered test must be equal on both sides (same datasets, same
             values within the tolerance). Mean, 95% interval and raw one-sided p(A>B) are
             recomputed with ``compare_units`` on each side for the reader; the primary's row of
             the report's declared-comparisons table (Holm-adjusted p) is shown for information.

A primary comparison is required: pass ``--record`` or ``--primary``, or ``--rows-only`` to compare
tasks and rows without one (the verdict then says so). A record or ``--primary`` that yields no
parsable ``A > B | metric | arm`` is an error (exit 2), never a pass. A primary is ``MATCH`` only if
the same usable units are on both sides and there are at least three of them, the minimum
``compare_units`` accepts; with fewer it is ``NOT EVALUABLE``: the registered test was not reproduced.

A malformed primary (text glued to the arm, an extra ``|`` field, a comparison that does not start
its segment), the same model on both sides, a negative or non-finite tolerance and an unreadable row
file are errors (exit 2). A list inside a row is compared by plain equality and ignores the tolerance
(it can give a false DIFFERENT, never a false MATCH). Read the ``primary:`` line, not only the verdict:
with ``--rows-only`` the verdict can say MATCH with no primary checked.

Prints one line per check, then ``verdict: MATCH`` (exit 0: tasks, rows and primary all match),
``verdict: DIFFERENT`` (exit 1) or ``verdict: NOT EVALUABLE`` (exit 1: nothing differs, but the
primary could not be run on both sides).
``primary: MATCH`` on its own is what a rule such as "the claim becomes reproduced only if the
primary comparison matches" reads. ``--json`` writes every difference, for the record.
"""
from __future__ import annotations

import argparse
import difflib
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from grnbench.eval.comparison import compare_units
from grnbench.texreport import _PRIMARY_RE, record_fields

METRIC_LABELS = {"auprc_ratio": "AUPRC ratio", "epr": "EPR"}  # as printed by real_bench.models_section
SHOWN = 20  # differences printed per check; --json keeps all of them
MIN_UNITS = 3  # compare_units refuses fewer independent units (a test pins this)


def is_timing(key: str) -> bool:
    return key == "seconds" or key.endswith("_seconds")


def _timing_header(cell: str) -> bool:
    return cell == "s" or cell.startswith("s per ") or "seconds" in cell


def load_rows(run: Path) -> dict[str, dict]:
    """rows/<id>.json -> row; an unreadable file or one that is not a JSON object is a ValueError naming it."""
    rows = {}
    for p in sorted((run / "rows").glob("*.json")):
        try:
            row = json.loads(p.read_text(encoding="utf-8"))
        except ValueError as exc:                     # JSONDecodeError and bad UTF-8 are ValueErrors
            raise ValueError(f"{p}: not valid JSON ({exc})") from exc
        if not isinstance(row, dict):
            raise ValueError(f"{p}: a row must be a JSON object, not {type(row).__name__}")
        rows[p.stem] = row
    return rows


def _kind(v) -> str:
    if isinstance(v, bool):
        return "bool"
    if isinstance(v, (int, float)):
        return "numeric"
    return "other"


def same(a, b, rtol: float = 0.0, atol: float = 0.0) -> bool:
    if _kind(a) != _kind(b):
        return False
    if _kind(a) == "numeric":
        if math.isnan(a) and math.isnan(b):
            return True
        return a == b or math.isclose(a, b, rel_tol=rtol, abs_tol=atol)
    return a == b


def compare_rows(old: dict, new: dict, rtol: float, atol: float) -> tuple[dict, list[str]]:
    """Field-by-field over the task ids present on both sides. Counts are by the committed value's type."""
    counts = {"numeric": 0, "bool": 0, "other": 0, "timing": 0}
    diffs = []
    for tid in sorted(set(old) & set(new)):
        a, b = old[tid], new[tid]
        for key in sorted(set(a) | set(b)):
            if is_timing(key):
                counts["timing"] += 1
                continue
            if key not in a or key not in b:
                diffs.append(f"{tid}\t{key}\tonly in {'committed' if key in a else 'rerun'}")
                continue
            counts[_kind(a[key])] += 1
            if not same(a[key], b[key], rtol, atol):
                diffs.append(f"{tid}\t{key}\t{a[key]!r} -> {b[key]!r}")
    return counts, diffs


def _cells(line: str) -> list[str]:
    return [c.strip() for c in line.strip().strip("|").split("|")]


def _is_separator(line: str) -> bool:
    s = line.strip()
    return s.startswith("|") and set(s) <= set("|-: ")


def _mask_timing(lines: list[str]) -> list[str]:
    """Each table row with its timing cells replaced by a placeholder (headers kept)."""
    out, timing = [], set()
    for i, line in enumerate(lines):
        header = line.startswith("|") and i + 1 < len(lines) and _is_separator(lines[i + 1])
        if not line.startswith("|"):
            timing = set()
        elif header:
            timing = {j for j, c in enumerate(_cells(line)) if _timing_header(c)}
        if timing and line.startswith("|") and not header and not _is_separator(line):
            line = "| " + " | ".join("<s>" if j in timing else c for j, c in enumerate(_cells(line))) + " |"
        out.append(line)
    return out


def compare_report(old_text: str, new_text: str) -> list[str]:
    """Aligned line diff (difflib), timing cells ignored: an inserted line shows as one line,
    not as every line after it."""
    a, b = old_text.splitlines(), new_text.splitlines()
    diffs = []
    matcher = difflib.SequenceMatcher(None, _mask_timing(a), _mask_timing(b), autojunk=False)
    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag != "equal":
            diffs += [f"- line {k + 1}: {a[k]!r}" for k in range(i1, i2)]
            diffs += [f"+ line {k + 1}: {b[k]!r}" for k in range(j1, j2)]
    return diffs


def primary_specs(record: Path | None, given: str | None) -> list[tuple[str, str, str, str]]:
    text = given if given is not None else record_fields(record).get("primary comparison", "")
    parts = [s.strip() for s in text.split(";") if s.strip()]
    if given is not None:
        matches = [_PRIMARY_RE.fullmatch(s.strip(" `")) for s in parts]
        if not parts or any(m is None for m in matches):
            return []
    else:
        # Record values can include explanatory continuation prose, including semicolons, after
        # the machine-readable comparison. A comparison must start the segment and end where its arm
        # ends: text glued to the arm ('auto-v2', 'auto.v2') or another field ('| two-sided') makes it
        # malformed. Prose-only segments are ignored; comparison-like segments that do not parse are
        # an error.
        matches = []
        for part in parts:
            match = _PRIMARY_RE.match(part)
            if match is not None and not re.match(r"[-.]\w|`?\s*[|>]", part[match.end():]):
                matches.append(match)
            elif ">" in part or "|" in part:
                return []
        if not matches:
            return []
    return [m.groups() for m in matches]


def primary_pair(rows: dict, a: str, b: str, metric: str, arm: str) -> pd.DataFrame:
    """The units entering the registered test: datasets with a score for both members
    (as real_bench.models_section builds the pair)."""
    ca, cb = f"m_{a}_{arm}_{metric}", f"m_{b}_{arm}_{metric}"
    df = pd.DataFrame({tid: {a: r.get(ca), b: r.get(cb)} for tid, r in rows.items()}, index=[a, b]).T
    return df.apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()


def primary_report_row(text: str, a: str, b: str, metric: str, arm: str) -> str | None:
    want = f"#### Declared comparisons, {METRIC_LABELS.get(metric, metric)}, {arm} pipeline"
    inside = False
    for line in text.splitlines():
        if line.startswith("#"):
            inside = line.startswith(want)
        elif inside and line.startswith("|") and _cells(line)[:2] == [a, b]:
            return line.strip()
    return None


def _summary(pair: pd.DataFrame) -> tuple[str, bool]:
    """One line for the reader, and whether compare_units can run on this pair."""
    if len(pair) < MIN_UNITS:
        return f"n={len(pair)} (too few units for the test)", False
    try:
        res = compare_units(pair)["pairs"][0]
    except ValueError as exc:
        return f"n={len(pair)} (compare_units refused: {exc})", False
    lo, hi = res["simultaneous_ci"]
    ci = "unbounded" if lo is None else f"[{lo:+.3f}, {hi:+.3f}]"
    return f"n={len(pair)} mean {res['difference']:+.3f} {ci} raw one-sided p(A>B) {res['p_a_gt_b']:.3g}", True


def compare_primary(old: dict, new: dict, spec, old_report: str | None, new_report: str | None,
                    rtol: float, atol: float) -> dict:
    a, b, metric, arm = spec
    po, pn = primary_pair(old, *spec), primary_pair(new, *spec)
    diffs = [f"unit {u} only in {'committed' if u in po.index else 'rerun'}"
             for u in sorted(set(po.index) ^ set(pn.index))]
    for u in sorted(set(po.index) & set(pn.index)):
        for m in (a, b):
            if not same(po.at[u, m], pn.at[u, m], rtol, atol):
                diffs.append(f"{u}\t{m}\t{float(po.at[u, m])!r} -> {float(pn.at[u, m])!r}")
    rows = [primary_report_row(t, *spec) if t is not None else None for t in (old_report, new_report)]
    (text_old, ok_old), (text_new, ok_new) = _summary(po), _summary(pn)
    usable = ok_old and ok_new
    status = "DIFFERENT" if diffs else ("MATCH" if usable else "NOT EVALUABLE")
    reason = None if status != "NOT EVALUABLE" else (
        f"the registered test needs at least {MIN_UNITS} usable units on both sides and could not run "
        f"(committed n={len(po)}, rerun n={len(pn)}): it was not reproduced")
    return {"spec": f"{a} > {b} | {metric} | {arm}", "status": status, "match": status == "MATCH",
            "reason": reason, "differences": diffs, "committed": text_old, "rerun": text_new,
            "report_rows": rows, "report_row_same": None if rows[0] is None else rows[0] == rows[1]}


def _tolerance(text: str) -> float:
    value = float(text)
    if not math.isfinite(value) or value < 0:
        raise argparse.ArgumentTypeError(f"{text!r}: a tolerance must be finite and not negative")
    return value


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Compare a rerun with the committed results of a registered experiment.")
    p.add_argument("committed", type=Path, help="committed results folder (rows/, REPORT_REAL.md)")
    p.add_argument("rerun", type=Path, help="the rerun's output folder")
    p.add_argument("--record", type=Path, help="experiment record holding 'primary comparison: A > B | metric | arm'")
    p.add_argument("--primary", help="'A > B | metric | arm' (overrides the record)")
    p.add_argument("--rows-only", action="store_true",
                   help="compare tasks and rows only, with no primary comparison (the verdict says so)")
    p.add_argument("--rtol", type=_tolerance, default=0.0)
    p.add_argument("--atol", type=_tolerance, default=0.0)
    p.add_argument("--json", type=Path, help="write every difference here")
    args = p.parse_args(argv)
    if args.rows_only and (args.record or args.primary is not None):
        p.error("--rows-only cannot be combined with --record or --primary")
    if not args.rows_only and args.primary is None and args.record is not None and not args.record.is_file():
        print(f"experiment record does not exist or is not a file: {args.record}", file=sys.stderr)
        return 2
    for run in (args.committed, args.rerun):
        if not (run / "rows").is_dir():
            print(f"no rows/ folder in {run}", file=sys.stderr)
            return 2
    specs = []
    if not args.rows_only:
        if not args.record and args.primary is None:
            print("no primary comparison given: pass --record or --primary, or --rows-only to compare "
                  "tasks and rows without one", file=sys.stderr)
            return 2
        specs = primary_specs(args.record, args.primary)
        if not specs:
            where = "--primary" if args.primary is not None else f"the record {args.record}"
            print(f"no primary comparison could be parsed from {where}: expected "
                  "'primary comparison: A > B | metric | arm'", file=sys.stderr)
            return 2
        for a, b, _, _ in specs:
            if a == b:
                print(f"the primary comparison names the same model twice ({a} > {b}): nothing to compare",
                      file=sys.stderr)
                return 2

    try:
        old, new = load_rows(args.committed), load_rows(args.rerun)
    except ValueError as exc:
        print(f"unreadable row file: {exc}", file=sys.stderr)
        return 2
    only_old, only_new = sorted(set(old) - set(new)), sorted(set(new) - set(old))
    tasks_ok = bool(old) and not only_old and not only_new
    print(f"tasks: {'same' if tasks_ok else 'DIFFERENT'} ({len(old)} committed, {len(new)} rerun"
          + (f"; only committed: {', '.join(only_old)}" if only_old else "")
          + (f"; only rerun: {', '.join(only_new)}" if only_new else "") + ")")

    counts, row_diffs = compare_rows(old, new, args.rtol, args.atol)
    tol = f"; tolerance rtol={args.rtol:g} atol={args.atol:g}" if args.rtol or args.atol else "; exact"
    print(f"rows: {'identical' if not row_diffs else f'{len(row_diffs)} DIFFERENT'} - compared "
          f"{counts['numeric']} numeric, {counts['bool']} bool and {counts['other']} other fields; "
          f"{counts['timing']} timing fields not compared{tol}")
    for d in row_diffs[:SHOWN]:
        print(f"  {d}")

    reports = [r / "REPORT_REAL.md" for r in (args.committed, args.rerun)]
    texts = [r.read_text(encoding="utf-8") if r.exists() else None for r in reports]
    if None in texts:
        report_diffs = [f"REPORT_REAL.md missing in {r.parent}" for r, t in zip(reports, texts) if t is None]
    else:
        report_diffs = compare_report(*texts)
    print("report: " + ("identical (line endings and timing columns ignored)" if not report_diffs
                        else f"{len(report_diffs)} lines differ (information only, not in the verdict: "
                             "with identical rows a difference comes from the report code)"))
    for d in report_diffs[:SHOWN]:
        print(f"  {d}")

    primaries = [compare_primary(old, new, s, *texts, args.rtol, args.atol) for s in specs]
    if args.rows_only:
        print("primary: not checked (--rows-only)")
    for pr in primaries:
        print(f"primary: {pr['status']} - {pr['spec']}")
        print(f"  committed: {pr['committed']}")
        print(f"  rerun:     {pr['rerun']}")
        if pr["reason"]:
            print(f"  {pr['reason']}")
        if pr["report_row_same"] is False:
            print(f"  report row differs (information): {pr['report_rows'][0]!r} -> {pr['report_rows'][1]!r}")
        for d in pr["differences"][:SHOWN]:
            print(f"  {d}")

    differs = not tasks_ok or bool(row_diffs) or any(pr["status"] == "DIFFERENT" for pr in primaries)
    unevaluable = any(pr["status"] == "NOT EVALUABLE" for pr in primaries)
    verdict = "DIFFERENT" if differs else "NOT EVALUABLE" if unevaluable else "MATCH"
    ok = verdict == "MATCH"
    print(f"verdict: {verdict}" + (" (rows only: no primary comparison checked)" if ok and args.rows_only else ""))
    if args.json:
        args.json.write_text(json.dumps({
            "verdict": verdict, "mode": "rows-only" if args.rows_only else "primary",
            "tolerance": {"rtol": args.rtol, "atol": args.atol},
            "tasks": {"committed": len(old), "rerun": len(new), "only_committed": only_old, "only_rerun": only_new},
            "rows": {"counts": counts, "differences": row_diffs},
            "report": {"differences": report_diffs, "in_verdict": False},
            "primary": primaries}, indent=1), encoding="utf-8")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
