"""Download public 10x Genomics raw and filtered feature-barcode matrices (issue #50).

    python benchmarks/fetch_10x_raw.py --out DIR [--only pbmc_1k_v3,heart_1k_v3]

Writes ``DIR/<name>/raw_feature_bc_matrix.h5``, ``DIR/<name>/filtered_feature_bc_matrix.h5`` and
``DIR/<name>/truth.json`` (organism, technology). Licence: 10x dataset terms, CC BY 4.0. Existing
files are kept, so a rerun only fetches what is missing. Record the bytes with
``python benchmarks/data_manifest.py DATA_ROOT`` and check them with ``--check``.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path
from urllib.request import Request, urlopen

BASE = "https://cf.10xgenomics.com/samples/cell-exp"
DATASETS = {  # name -> (Cell Ranger release folder, organism)
    "pbmc_1k_v3": ("3.0.0", "human"),
    "pbmc_10k_v3": ("3.0.0", "human"),
    "5k_pbmc_protein_v3": ("3.0.2", "human"),
    "pbmc_1k_v2": ("3.0.0", "human"),
    "neuron_1k_v3": ("3.0.0", "mouse"),
    "neuron_10k_v3": ("3.0.0", "mouse"),
    "heart_1k_v3": ("3.0.0", "mouse"),
    "heart_10k_v3": ("3.0.0", "mouse"),
}
KINDS = ("raw", "filtered")


def url(name: str, kind: str) -> str:
    release, _ = DATASETS[name]
    return f"{BASE}/{release}/{name}/{name}_{kind}_feature_bc_matrix.h5"


def fetch(out: Path, names: list[str]) -> int:
    for name in names:
        folder = out / name
        folder.mkdir(parents=True, exist_ok=True)
        for kind in KINDS:
            dest = folder / f"{kind}_feature_bc_matrix.h5"
            if dest.exists():
                print(f"have {dest}")
                continue
            tmp = dest.with_suffix(".part")
            # The 10x CDN answers 403 to urllib's default User-Agent.
            req = Request(url(name, kind), headers={"User-Agent": "grnbench-fetch/1.0"})
            with urlopen(req, timeout=120) as resp, open(tmp, "wb") as fh:
                shutil.copyfileobj(resp, fh, 1 << 20)
            tmp.replace(dest)
            print(f"ok {dest} {dest.stat().st_size} bytes")
        truth = {"technology": "umi", "organism": DATASETS[name][1], "source": "10x Genomics public dataset",
                 "url_raw": url(name, "raw"), "url_filtered": url(name, "filtered")}
        (folder / "truth.json").write_text(json.dumps(truth, indent=1), encoding="utf-8")
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", required=True)
    p.add_argument("--only", default=None, help="comma-separated subset of " + ", ".join(DATASETS))
    args = p.parse_args(argv)
    names = [n.strip() for n in args.only.split(",")] if args.only else list(DATASETS)
    bad = [n for n in names if n not in DATASETS]
    if bad:
        p.error(f"unknown datasets {bad}")
    return fetch(Path(args.out), names)


if __name__ == "__main__":
    sys.exit(main())
