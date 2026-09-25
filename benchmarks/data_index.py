"""Write INDEX.md for a grnbench_data folder: one row per dataset.

    python benchmarks/data_index.py DEST
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path


def _size(path: Path) -> float:
    if path.is_file():
        return path.stat().st_size / 1e6
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file()) / 1e6


def main(dest: str) -> int:
    root = Path(dest)
    rows = []
    man = root / "cellxgene" / "manifest.tsv"
    if man.exists():
        with open(man, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                folder = root / "cellxgene" / r["folder"]
                if (folder / "data.h5ad").exists():
                    rows.append((f"cellxgene/{r['folder']}", "CELLxGENE Discover", r["assay"], r["organism"],
                                 r["cell_count"], f"{_size(folder):.0f}", r["title"][:70]))
    for folder, src, proto, org, cells, note in [
            ("10x_pbmc3k/filtered", "10x Genomics", "10x Chromium v1", "Homo sapiens", "~2,700", "Cell Ranger cells"),
            ("10x_pbmc3k/raw", "10x Genomics", "10x Chromium v1", "Homo sapiens", "all barcodes", "unfiltered droplets"),
            ("scanpy_paul15", "scanpy / Paul et al. 2015", "MARS-seq", "Mus musculus", "2,730", "differentiation"),
            ("beeline", "Zenodo 19009603 (BEELINE v1.1)", "various (processed)", "human + mouse", "7 datasets",
             "log-normalized inputs + reference networks")]:
        p = root / folder
        if p.exists():
            rows.append((folder, src, proto, org, cells, f"{_size(p):.0f}", note))
    refs = sorted((root / "references").glob("*.tsv"))
    lines = ["# grnbench_data", "", "Real datasets for the engine mass test "
             "(fetched by hpc/fetch_engine_datasets.sh).", "",
             "| folder | source | protocol | organism | cells | MB | note |", "|---|---|---|---|---|---|---|"]
    lines += ["| " + " | ".join(map(str, r)) + " |" for r in rows]
    lines += ["", "Reference networks: " + (", ".join(f"`references/{r.name}`" for r in refs) or "none"), ""]
    (root / "INDEX.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"INDEX.md: {len(rows)} datasets, {len(refs)} reference networks")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "."))
