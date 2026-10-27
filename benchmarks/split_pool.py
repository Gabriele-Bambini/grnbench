"""Split a dataset pool into a confirmatory and an exploratory part, fixed by a seed before any look.

Each dataset gets the key sha256("<seed>:<dataset_id>"); the n datasets with the smallest keys
are confirmatory, the others exploratory, each part in key order. The split depends only on the
seed and the dataset ids: not on the row order of the manifest, the data or any result. Both
outputs keep the input columns, so ``fetch_cellxgene.py --from-manifest`` reads them.

    python benchmarks/split_pool.py --manifest POOL.tsv --seed 42 --confirmatory 120 \\
        --out-confirmatory C.datasets.tsv --out-exploratory E.datasets.tsv
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from pathlib import Path


def _key(seed: int, dataset_id: str) -> str:
    return hashlib.sha256(f"{seed}:{dataset_id}".encode()).hexdigest()


def split(rows: list[dict], seed: int, n_confirmatory: int) -> tuple[list[dict], list[dict]]:
    ids = [r["dataset_id"] for r in rows]
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate dataset ids in the pool")
    if not 0 < n_confirmatory <= len(rows):
        raise ValueError(f"cannot take {n_confirmatory} confirmatory datasets from a pool of {len(rows)}")
    ordered = sorted(rows, key=lambda r: _key(seed, r["dataset_id"]))
    return ordered[:n_confirmatory], ordered[n_confirmatory:]


def _write(path: Path, cols: list[str], rows: list[dict]) -> None:
    lines = ["\t".join(cols)] + ["\t".join(str(r[c]) for c in cols) for r in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--seed", required=True, type=int)
    p.add_argument("--confirmatory", required=True, type=int, help="number of confirmatory datasets")
    p.add_argument("--out-confirmatory", required=True, type=Path)
    p.add_argument("--out-exploratory", required=True, type=Path)
    a = p.parse_args(argv)
    paths = [x.resolve() for x in (a.manifest, a.out_confirmatory, a.out_exploratory)]
    if len(set(paths)) != 3:
        p.error("--manifest, --out-confirmatory and --out-exploratory must be three different files")
    with open(a.manifest, encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        rows, cols = list(reader), list(reader.fieldnames or [])
    if "dataset_id" not in cols:
        raise SystemExit(f"{a.manifest}: no dataset_id column")
    short = [i + 2 for i, r in enumerate(rows) if None in r.values() or None in r]
    if short:
        raise SystemExit(f"{a.manifest}: rows with missing or extra columns at lines {short[:5]}")
    conf, expl = split(rows, a.seed, a.confirmatory)
    _write(a.out_confirmatory, cols, conf)
    _write(a.out_exploratory, cols, expl)
    print(f"{len(conf)} confirmatory, {len(expl)} exploratory (seed {a.seed}, pool {len(rows)})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
