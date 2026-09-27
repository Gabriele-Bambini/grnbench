"""Data provenance manifest: one row per data file under a data root, with source and sha256.

    python benchmarks/data_manifest.py DATA_ROOT > experiments/datasets.tsv
    python benchmarks/data_manifest.py DATA_ROOT --check experiments/datasets.tsv

``--check`` verifies every data file under DATA_ROOT against the manifest (same path, same
sha256) and exits 1 if any file differs or is not listed. A text file (.json/.tsv/.txt/.tbl)
whose only difference is LF vs CRLF line endings counts as a match and is reported as such:
fetch scripts write text in the platform's line endings, so a truth.json written on Linux
differs in bytes from the Windows copy that was hashed.

Paths are relative to DATA_ROOT, so the manifest carries no personal path. Commit the TSV;
never commit the data. Experiment records cite datasets by the `path` column.
"""
from __future__ import annotations

import csv
import datetime as dt
import hashlib
import sys
from pathlib import Path

SOURCES = {  # prefix of the relative path -> (source, version note)
    "10x_pbmc3k/": ("https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/", "Cell Ranger 1.1.0"),
    "scanpy_paul15/": ("scanpy.datasets.paul15 (Paul et al. 2015, Cell)", "scanpy version in the lock"),
    "beeline/": ("https://zenodo.org/records/19009603", "BEELINE v1.1 (Zenodo v4)"),
    "references/collectri": ("OmniPath CollecTRI via decoupler", "decoupler version in the lock"),
    "references/dorothea": ("OmniPath DoRothEA ABC via decoupler", "decoupler version in the lock"),
    "cistarget/": ("https://resources.aertslab.org/cistarget/", "cisTarget v10 clustered"),
}
SKIP = {"INDEX.md", "README.md", "manifest.tsv", "desktop.ini"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main(root_arg: str) -> int:
    root = Path(root_arg)
    cxg = {}
    man = root / "cellxgene" / "manifest.tsv"
    if man.exists():
        with open(man, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                cxg[f"cellxgene/{r['folder']}/"] = (r["url"], f"CELLxGENE dataset_id {r['dataset_id']}")
    out = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    out.writerow(["path", "bytes", "sha256", "file_date", "source", "version"])
    for f in sorted(p for p in root.rglob("*") if p.is_file() and p.name not in SKIP):
        rel = f.relative_to(root).as_posix()
        src, ver = next((v for k, v in {**SOURCES, **cxg}.items() if rel.startswith(k)), ("UNKNOWN - add it", ""))
        file_date = dt.date.fromtimestamp(f.stat().st_mtime).isoformat()
        out.writerow([rel, f.stat().st_size, sha256(f), file_date, src, ver])
    return 0


TEXT_SUFFIXES = {".json", ".tsv", ".txt", ".tbl"}


def check(root_arg: str, manifest: str) -> int:
    """Compare the files under DATA_ROOT with the manifest rows; 0 only if all match."""
    root = Path(root_arg)
    with open(manifest, encoding="utf-8") as fh:
        want = {r["path"]: r["sha256"] for r in csv.DictReader(fh, delimiter="\t")}
    counts = {"ok": 0, "ok-eol": 0, "MISMATCH": 0, "UNLISTED": 0}
    for f in sorted(p for p in root.rglob("*") if p.is_file() and p.name not in SKIP):
        rel = f.relative_to(root).as_posix()
        if rel not in want:
            status = "UNLISTED"
        elif sha256(f) == want[rel]:
            status = "ok"
        elif f.suffix in TEXT_SUFFIXES and _eol_sha(f) == want[rel]:
            status = "ok-eol"
        else:
            status = "MISMATCH"
        counts[status] += 1
        if status != "ok":
            print(f"{status}\t{rel}")
    print(" ".join(f"{k}={v}" for k, v in counts.items()))
    return 1 if counts["MISMATCH"] or counts["UNLISTED"] else 0


def _eol_sha(path: Path) -> str:
    """sha256 of the file with its line endings flipped (LF <-> CRLF)."""
    data = path.read_bytes()
    flipped = (data.replace(b"\r\n", b"\n") if b"\r\n" in data
               else data.replace(b"\n", b"\r\n"))
    return hashlib.sha256(flipped).hexdigest()


if __name__ == "__main__":
    if "--check" in sys.argv:
        i = sys.argv.index("--check")
        sys.exit(check(sys.argv[1] if i > 1 else ".", sys.argv[i + 1]))
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "."))
