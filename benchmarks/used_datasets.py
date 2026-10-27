"""Datasets and studies that earlier work used, found from every identifier in records and results.

Scans ``experiments/`` and ``benchmarks/results/`` at the given git refs (default: HEAD and every
remote branch) for CELLxGENE identifiers: 36-character ids, the ``__xxxxxxxx`` suffix of dataset
folders, and lines holding only an 8-character id prefix (as in ``*.excluded_ids.txt``). Each is
resolved against a saved copy of the CELLxGENE listing: a dataset id, the id of a dataset's H5AD
download file, or the 8-character prefix of a dataset id (every dataset with that prefix, to be safe).
Writes the dataset ids with their download-file ids and the unresolved identifiers (for
``fetch_cellxgene.py --exclude-ids``), and the collections of the datasets (``--exclude-collections``).

    python benchmarks/used_datasets.py --listing LISTING.json --skip experiments/EXP-20261022-GB \\
        --out-ids experiments/X.excluded_ids.txt --out-collections experiments/X.excluded_collections.txt
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TOKEN = r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}|__[0-9a-f]{8}"
BARE_PREFIX = r"^[0-9a-f]{8}\s*$"
FILE_ID = re.compile(r"/([0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})\.h5ad", re.I)


def _grep(ref: str, pattern: str, paths, skip) -> list[str]:
    cmd = ["git", "grep", "-o", "-I", "-i", "-E", pattern, ref, "--", *paths, *[f":!{s}" for s in skip]]
    p = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    if p.returncode not in (0, 1):      # 1 = no match
        raise RuntimeError(f"git grep failed on {ref}: {p.stderr.strip()}")
    return [line.rsplit(":", 1)[1].strip().lower().lstrip("_") for line in p.stdout.splitlines()]


def scan(refs, paths=("experiments", "benchmarks/results"), skip=()) -> set[str]:
    tokens: set[str] = set()
    for ref in refs:
        tokens.update(_grep(ref, TOKEN, paths, skip))
        tokens.update(_grep(ref, BARE_PREFIX, paths, skip))
    return tokens


def resolve(tokens, listing) -> tuple[set[str], set[str]]:
    """``(used dataset ids, identifiers that match nothing in the listing)``."""
    ids = {d["dataset_id"].lower() for d in listing}
    by_file = {}
    for d in listing:
        for asset in d.get("assets") or []:
            m = FILE_ID.search(asset.get("url") or "")
            if m:
                by_file[m.group(1).lower()] = d["dataset_id"].lower()
    used, unresolved = set(), set()
    for t in tokens:
        if len(t) == 36:
            hits = {t} if t in ids else {by_file[t]} if t in by_file else set()
        else:
            hits = {i for i in ids if i.startswith(t)}
        if hits:
            used |= hits
        else:
            unresolved.add(t)
    return used, unresolved


def _remote_refs() -> list[str]:
    p = subprocess.run(["git", "for-each-ref", "--format=%(refname)", "refs/remotes/"],
                       capture_output=True, text=True, cwd=ROOT, check=True)
    return [r for r in p.stdout.split() if not r.endswith("/HEAD")]


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--listing", required=True, help="saved JSON of the CELLxGENE curation API listing")
    p.add_argument("--refs", nargs="*", default=None, help="git refs to scan (default: HEAD and every remote branch)")
    p.add_argument("--skip", nargs="*", default=[], help="path prefixes to leave out (the new experiment's own files)")
    p.add_argument("--out-ids", required=True, type=Path)
    p.add_argument("--out-collections", required=True, type=Path)
    a = p.parse_args(argv)
    listing = json.loads(Path(a.listing).read_text(encoding="utf-8"))
    refs = a.refs if a.refs is not None else ["HEAD", *_remote_refs()]
    tokens = scan(refs, skip=[s if s.endswith("*") else s + "*" for s in a.skip])
    used, unresolved = resolve(tokens, listing)
    files = set()
    for d in listing:
        if d["dataset_id"].lower() in used:
            for asset in d.get("assets") or []:
                m = FILE_ID.search(asset.get("url") or "")
                if m:
                    files.add(m.group(1).lower())
    colls = {(d.get("collection_id") or "").lower() for d in listing if d["dataset_id"].lower() in used}
    if "" in colls:
        raise SystemExit("a used dataset has no collection_id in the listing")
    a.out_ids.write_text("\n".join(sorted(used | files | unresolved)) + "\n", encoding="utf-8", newline="\n")
    a.out_collections.write_text("\n".join(sorted(colls)) + "\n", encoding="utf-8", newline="\n")
    print(f"{len(refs)} refs, {len(tokens)} identifiers: {len(used)} datasets in {len(colls)} collections, "
          f"{len(unresolved)} unresolved (kept in the id list)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
