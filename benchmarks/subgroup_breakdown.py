"""Descriptive breakdown of a paired per-dataset difference by organism, assay family and suspension.

Declared in EXP-20261022-GB-1 and EXP-20261022-GB-2 before their runs, so a pooled mean cannot hide a
subgroup that goes the other way. No test and no claim: per subgroup, the number of datasets, the mean
difference and its 95% t interval (none for a single dataset). Datasets are matched to the pool
manifest by the 8-character id suffix of their task id.

    python benchmarks/subgroup_breakdown.py --rows RUN/all_rows_real.csv \\
        --datasets experiments/EXP-20261022-GB-1.datasets.tsv --out RUN/SUBGROUPS.md \\
        [--pair pearson_motif,pearson_signed] [--metric epr] [--arm auto]
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

FULL_LENGTH = ("smart-seq", "patch-seq", "fluidigm")
ORGANISMS = {"Homo sapiens": "human", "Mus musculus": "mouse"}
SUSPENSION = {"cell": "cell", "nucleus": "nucleus", "cell;nucleus": "both"}


def assay_family(label: str) -> str:
    low = label.lower()
    if low.startswith("10x 3'"):
        return "10x 3'"
    if low.startswith("10x 5'"):
        return "10x 5'"
    if "multiome" in low:
        return "10x multiome"
    return "full-length" if any(k in low for k in FULL_LENGTH) else "other"


def load_meta(datasets_tsv) -> dict[str, dict]:
    """8-character id suffix of each folder -> organism, assay family, suspension."""
    with open(datasets_tsv, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE))
    return {r["folder"][-8:]: {"organism": ORGANISMS.get(r["organism"], r["organism"]),
                               "assay family": assay_family(r["assay"]),
                               "suspension": SUSPENSION.get(r["suspension_type"], r["suspension_type"])}
            for r in rows}


def real_differences(rows: pd.DataFrame, pair: tuple[str, str], metric: str, arm: str) -> dict[str, float]:
    """A - B per dataset from a real_bench row table; datasets without both scores are left out."""
    a, b = (rows[f"m_{m}_{arm}_{metric}"] for m in pair)
    keep = a.notna() & b.notna()
    return {str(t)[-8:]: float(x - y) for t, x, y in zip(rows["task"][keep], a[keep], b[keep])}


def breakdown(diffs: dict[str, float], meta: dict[str, dict]) -> pd.DataFrame:
    missing = sorted(k for k in diffs if k not in meta)
    if missing:
        raise ValueError(f"scored datasets not in the manifest: {missing}")
    out = []
    for factor in ("organism", "assay family", "suspension"):
        levels = sorted({meta[k][factor] for k in diffs})
        for level in levels:
            x = np.array([v for k, v in diffs.items() if meta[k][factor] == level])
            lo = hi = math.nan
            if len(x) > 1:
                half = float(stats.t.ppf(0.975, len(x) - 1) * x.std(ddof=1) / math.sqrt(len(x)))
                lo, hi = float(x.mean()) - half, float(x.mean()) + half
            out.append({"factor": factor, "level": level, "n": len(x), "mean": float(x.mean()), "lo": lo, "hi": hi})
    return pd.DataFrame(out)


def render(table: pd.DataFrame, what: str, left_out: int = 0) -> str:
    lines = [f"# Subgroup breakdown: {what}", "",
             "Descriptive, declared before the run; no test and no claim. Mean paired difference per subgroup "
             "with its 95% t interval (none for one dataset).",
             f"{left_out} dataset(s) without both scores left out.", "",
             "| factor | subgroup | n | mean difference | 95% interval |", "|---|---|---|---|---|"]
    for r in table.itertuples():
        ci = "-" if math.isnan(r.lo) else f"[{r.lo:+.3f}, {r.hi:+.3f}]"
        lines.append(f"| {r.factor} | {r.level} | {r.n} | {r.mean:+.3f} | {ci} |")
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--rows", required=True, type=Path, help="all_rows_real.csv of a real_bench run")
    p.add_argument("--datasets", required=True, type=Path, help="the record's datasets.tsv")
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--pair", default="pearson_motif,pearson_signed")
    p.add_argument("--metric", default="epr")
    p.add_argument("--arm", default="auto")
    a = p.parse_args(argv)
    pair = tuple(a.pair.split(","))
    rows = pd.read_csv(a.rows)
    diffs = real_differences(rows, pair, a.metric, a.arm)
    table = breakdown(diffs, load_meta(a.datasets))
    text = render(table, f"{a.metric}, {pair[0]} - {pair[1]}, {a.arm} arm", left_out=len(rows) - len(diffs))
    a.out.write_text(text, encoding="utf-8", newline="\n")
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
