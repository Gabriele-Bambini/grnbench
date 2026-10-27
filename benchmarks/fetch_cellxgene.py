"""Pick and download a diverse set of small CELLxGENE Discover datasets.

Selection (from the public curation API listing): single-cell/nucleus assays
only (no spatial), H5AD <= ``--max-mb``, 500-50,000 cells; round-robin over
assay x organism so protocols and species are all represented, tissues varied.
``manifest.tsv`` records the curated assay / organism / tissue of every
dataset - the ground truth the engine's context detection is checked against.

    python benchmarks/fetch_cellxgene.py --out DIR --max 30 --max-mb 100 [--dry-run]
        [--exclude-ids used_ids.txt] [--exclude-assay targeted]
    python benchmarks/fetch_cellxgene.py --out DIR --from-manifest experiments/EXP-...datasets.tsv

``--from-manifest`` skips the selection and fetches exactly the rows of an earlier
``manifest.tsv`` (same folders, URLs and truth.json), so a rerun gets the same files even if
the public listing has changed since. Check them with ``data_manifest.py --check``. It reads
the listing only for the collection name and citation of each README.md (the engine reads
README.md as a text source), and says so on stderr for each dataset it cannot fill; the
downloads do not depend on it, and ``--dry-run`` does not read it.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from http.client import HTTPException
from pathlib import Path
from urllib.request import Request, urlopen

API = "https://api.cellxgene.cziscience.com/curation/v1/datasets?visibility=PUBLIC"
SPATIAL = ("slide-seq", "visium", "seeker", "merfish", "stereo")
FULL_LENGTH = ("smart-seq", "smart seq", "patch-seq", "fluidigm")


def assay_technology(label: str) -> str:
    low = label.lower()
    return "full_length" if any(k in low for k in FULL_LENGTH) else "umi"


def slug(text: str, n: int = 40) -> str:
    out = "".join(c.lower() if c.isalnum() else "-" for c in text)
    while "--" in out:
        out = out.replace("--", "-")
    return out.strip("-")[:n].strip("-") or "unknown"


def folder_for(d: dict) -> Path:
    """cellxgene/<protocol>/<organism>/<tissue>__<id8>/ - readable and unique."""
    tissue = d["tissue"].split(";")[0] if d["tissue"] else "tissue"
    return Path(slug(d["assay"])) / slug(d["organism"]) / f"{slug(tissue, 30)}__{d['dataset_id'][:8]}"


def _labels(items) -> list[str]:
    return [i.get("label", "") for i in items or []]


def select(datasets: list[dict], n: int, max_mb: float) -> list[dict]:
    pool = []
    for d in datasets:
        h5 = [a for a in d.get("assets", []) if a.get("filetype") == "H5AD"]
        assays, orgs = _labels(d.get("assay")), _labels(d.get("organism"))
        if not h5 or len(assays) != 1 or len(orgs) != 1 or d.get("spatial"):
            continue
        if any(k in assays[0].lower() for k in SPATIAL):
            continue
        size = h5[0].get("filesize") or 0
        cells = d.get("cell_count") or 0
        if not (0 < size <= max_mb * 1e6 and 500 <= cells <= 50_000):
            continue
        pool.append({"dataset_id": d["dataset_id"], "title": d.get("title", ""),
                     "collection": d.get("collection_name", ""), "citation": d.get("citation", ""),
                     "assay": assays[0], "organism": orgs[0],
                     "tissue": ";".join(_labels(d.get("tissue"))[:3]),
                     "suspension_type": ";".join(d.get("suspension_type") or []),
                     "cell_count": cells, "filesize": size, "url": h5[0]["url"]})
    groups: dict[tuple, list] = {}
    for d in sorted(pool, key=lambda d: (d["filesize"], d["dataset_id"])):
        groups.setdefault((d["assay"], d["organism"]), []).append(d)
    picked, seen_tissue = [], set()
    while len(picked) < n and any(groups.values()):
        for key in sorted(groups):
            g = groups[key]
            if not g or len(picked) >= n:
                continue
            fresh = next((d for d in g if d["tissue"] not in seen_tissue), g[0])
            g.remove(fresh)
            picked.append(fresh)
            seen_tissue.add(fresh["tissue"])
    return picked


def _readme_text_from_listing(picked: list[dict]) -> None:
    """Fill ``collection`` and ``citation`` of manifest rows from the listing, by dataset id, so
    README.md is the text the listing route writes. Never fails the download: a dataset the
    listing cannot fill keeps empty values and gets one loud line on stderr."""
    try:
        with urlopen(Request(API, headers={"User-Agent": "grnbench"}), timeout=120) as r:
            listing, why = {d["dataset_id"]: d for d in json.load(r)}, "not in the CELLxGENE listing"
    except (OSError, ValueError, HTTPException, KeyError, TypeError) as e:   # unreachable, not JSON, malformed
        listing, why = {}, f"CELLxGENE listing not fetched ({e})"
    for d in picked:
        hit = listing.get(d["dataset_id"])
        if hit is None:
            print(f"WARNING {d['dataset_id']} ({d['folder']}): {why}; its README.md has no collection "
                  "name or citation, so it differs from the one the listing route writes. The engine reads "
                  "README.md as a text source: planner warnings may differ from the original run.",
                  file=sys.stderr, flush=True)
            continue
        d["collection"], d["citation"] = hit.get("collection_name", ""), hit.get("citation", "")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", required=True)
    p.add_argument("--max", type=int, default=30)
    p.add_argument("--max-mb", type=float, default=100)
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--organisms", default=None,
                   help="comma-separated curated organism labels to keep, e.g. 'Homo sapiens,Mus musculus' "
                        "(only these have CollecTRI references to score against)")
    p.add_argument("--exclude-ids", default=None,
                   help="file with dataset ids (or 8-character id prefixes), one per line, to skip - "
                        "e.g. every dataset used by earlier experiments, for a fresh confirmatory set")
    p.add_argument("--exclude-assay", default=None,
                   help="comma-separated case-insensitive substrings of assay labels to skip, e.g. 'targeted'")
    p.add_argument("--from-manifest", default=None,
                   help="fetch exactly the datasets of this manifest.tsv (no selection; the listing "
                        "only supplies the collection name and citation of README.md)")
    a = p.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    if a.from_manifest:
        with open(a.from_manifest, encoding="utf-8") as fh:
            picked = [dict(r, filesize=int(r["filesize"]), cell_count=int(r["cell_count"]),
                           collection="", citation="") for r in csv.DictReader(fh, delimiter="\t")]
        if not a.dry_run:
            _readme_text_from_listing(picked)
        return _write_and_fetch(out, picked, a.dry_run)
    with urlopen(Request(API, headers={"User-Agent": "grnbench"}), timeout=120) as r:
        datasets = json.load(r)
    if a.organisms:
        keep = {o.strip().lower() for o in a.organisms.split(",")}
        datasets = [d for d in datasets if {x.lower() for x in _labels(d.get("organism"))} <= keep]
    if a.exclude_ids:
        ids = [x.strip().lower() for x in Path(a.exclude_ids).read_text(encoding="utf-8").split() if x.strip()]
        datasets = [d for d in datasets if not any(d["dataset_id"].lower().startswith(i) for i in ids)]
    if a.exclude_assay:
        subs = [x.strip().lower() for x in a.exclude_assay.split(",") if x.strip()]
        datasets = [d for d in datasets if not any(s_ in lab.lower() for lab in _labels(d.get("assay"))
                                                   for s_ in subs)]
    picked = select(datasets, a.max, a.max_mb)
    for d in picked:
        d["folder"] = folder_for(d).as_posix()
    return _write_and_fetch(out, picked, a.dry_run)


def _write_and_fetch(out: Path, picked: list[dict], dry_run: bool) -> int:
    cols = ["folder", "dataset_id", "assay", "organism", "tissue", "suspension_type", "cell_count",
            "filesize", "title", "url"]
    lines = ["\t".join(cols)] + ["\t".join(str(d[c]).replace("\t", " ") for c in cols) for d in picked]
    (out / "manifest.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    total = sum(d["filesize"] for d in picked) / 1e6
    print(f"{len(picked)} datasets, {total:.0f} MB total")
    for d in picked:
        print(f"  {d['assay']:28s} {d['organism']:14s} {d['cell_count']:>7} cells "
              f"{d['filesize'] / 1e6:6.1f} MB  {d['title'][:60]}")
    if dry_run:
        return 0
    for d in picked:
        folder = out / d["folder"]
        folder.mkdir(parents=True, exist_ok=True)
        target = folder / "data.h5ad"
        if not (target.exists() and target.stat().st_size == d["filesize"]):
            with urlopen(Request(d["url"], headers={"User-Agent": "grnbench"}), timeout=600) as r, \
                    open(target, "wb") as fh:
                while chunk := r.read(1 << 20):
                    fh.write(chunk)
            time.sleep(1)
        (folder / "README.md").write_text(
            f"{d['title']}\n\nCollection: {d['collection']}\n{d['citation']}\n", encoding="utf-8")
        (folder / "truth.json").write_text(json.dumps(
            {"assay": d["assay"], "technology": assay_technology(d["assay"]),
             "organism": "human" if "sapiens" in d["organism"].lower() else
             "mouse" if "musculus" in d["organism"].lower() else d["organism"],
             "tissue": d["tissue"]}, indent=1), encoding="utf-8")
        print(f"ok {d['folder']}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
