"""Download CRISPRi Perturb-seq screens from the scPerturb collection (Zenodo 7041849, CC BY 4.0).

    python benchmarks/fetch_perturb.py --out DIR [--only k562_gwps,rpe1] [--parallel 8]

Replogle et al. 2022 knockdowns are the intervention truth for the knockdown experiment; Norman 2019
(CRISPRa) is listed for completeness. Each person's data root is different, so DIR is always an
argument; the files land directly in DIR (nothing is appended), so use ``<data root>/perturb`` for
``data_manifest.py`` to find their source. Files are pinned by name, size and Zenodo md5: a rerun skips files whose md5 already
matches, resumes every unfinished byte range, and refuses a file whose md5 differs. The parts are
hashed in order first; only if that md5 matches are they appended onto part00 one at a time, each deleted
once copied, so joining needs one part of free space, not a second copy of the file. Each file is
fetched as ``--parallel`` byte ranges at once: on the Bocconi login node one connection gets about
0.25 MB/s and eight together about 12 MB/s. Record the bytes afterwards with
``python benchmarks/data_manifest.py DATA_ROOT``. Login node: downloads only.
"""
from __future__ import annotations

import argparse
import hashlib
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import Request, urlopen

RECORD = "https://zenodo.org/api/records/7041849/files"
FILES = {  # key -> (file name, bytes, md5) as published in record 7041849 v1.0
    "k562_gwps": ("ReplogleWeissman2022_K562_gwps.h5ad", 8805466154, "13db594f8f1d2ccb88fec44a13e414dc"),
    "k562_essential": ("ReplogleWeissman2022_K562_essential.h5ad", 1546729675, "d8cba17576d1a8afc0f7d71b79cad0f7"),
    "rpe1": ("ReplogleWeissman2022_rpe1.h5ad", 1236886900, "cc7f1ec50aeb3a3e1b4a6cfa713d80fa"),
    "norman2019": ("NormanWeissman2019_filtered.h5ad", 698680199, "c870e6967d91c017d9da827bab183cd6"),
}
DEFAULT = ("k562_gwps", "rpe1")
CHUNK = 1 << 20
TRIES = 5


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(CHUNK), b""):
            h.update(block)
    return h.hexdigest()


def fetch_range(url: str, part: Path, start: int, end: int) -> None:
    """Bytes ``start..end`` (inclusive) into ``part``, continuing from what ``part`` already holds."""
    for attempt in range(TRIES):
        have = part.stat().st_size if part.exists() else 0
        if start + have > end:
            return
        try:
            req = Request(url, headers={"Range": f"bytes={start + have}-{end}"})
            with urlopen(req, timeout=60) as r, part.open("ab") as f:
                if r.status != 206:
                    raise OSError(f"server ignored the byte range (HTTP {r.status})")
                while block := r.read(CHUNK):
                    f.write(block)
        except OSError as e:
            if attempt == TRIES - 1:
                raise
            print(f"retry {part.name}: {e}", file=sys.stderr, flush=True)


def fetch_one(out: Path, key: str, parallel: int) -> bool:
    name, size, want = FILES[key]
    dest = out / name
    if dest.exists() and md5(dest) == want:
        print(f"have {name}")
        return True
    print(f"=== {name} ({size / 1e9:.2f} GB, {parallel} ranges)", flush=True)
    step = -(-size // parallel)
    ranges = [(i, i * step, min(size, (i + 1) * step) - 1) for i in range(parallel) if i * step < size]
    parts = [dest.with_name(f"{name}.part{i:02d}") for i, _, _ in ranges]
    old = sorted(dest.parent.glob(f"{name}.part*"))
    if old and len(old) != len(parts):
        print(f"found {len(old)} unfinished parts of {name} but this run uses {len(parts)}: resume with the "
              f"same --parallel as the interrupted run ({len(old)}), or delete the .part files to restart",
              file=sys.stderr, flush=True)
    url = f"{RECORD}/{name}/content"
    with ThreadPoolExecutor(len(ranges)) as pool:
        for job in [pool.submit(fetch_range, url, p, s, e) for p, (_, s, e) in zip(parts, ranges)]:
            job.result()
    h = hashlib.md5()
    for p in parts:                      # hash the parts in order without writing a second copy
        with p.open("rb") as src:
            for block in iter(lambda: src.read(CHUNK), b""):
                h.update(block)
    if h.hexdigest() != want:
        print(f"md5 mismatch for {name}: got {h.hexdigest()}, want {want}; kept the .partNN files, "
              "delete them to restart", file=sys.stderr)
        return False
    try:                                 # append onto part00, deleting each part once copied: peak extra space = one part
        with parts[0].open("ab") as f:
            for p in parts[1:]:
                with p.open("rb") as src:
                    shutil.copyfileobj(src, f, CHUNK)
                p.unlink()
        parts[0].replace(dest)
    except OSError as e:
        with parts[0].open("r+b") as f:  # cut part00 back to its own range, so a rerun refetches the deleted parts
            f.truncate(ranges[0][2] - ranges[0][1] + 1)
        print(f"join of {name} failed: {e}; free space and rerun with the same --parallel "
              "(parts already joined are fetched again)", file=sys.stderr, flush=True)
        return False
    for p in dest.parent.glob(f"{name}.part*"):
        p.unlink()
    dest.with_name(name + ".joining").unlink(missing_ok=True)
    print(f"ok {name}", flush=True)
    return True


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--out", type=Path, required=True, help="destination folder, e.g. ~/grnbench_data/perturb")
    p.add_argument("--only", default=",".join(DEFAULT), help=f"comma-separated keys from {sorted(FILES)}")
    p.add_argument("--parallel", type=int, default=8, help="byte ranges fetched at once per file")
    a = p.parse_args(argv)
    keys = [k.strip() for k in a.only.split(",") if k.strip()]
    unknown = sorted(set(keys) - set(FILES))
    if unknown:
        p.error(f"unknown keys {unknown}; choose from {sorted(FILES)}")
    if a.parallel < 1:
        p.error("--parallel must be at least 1")
    out = a.out.expanduser()
    out.mkdir(parents=True, exist_ok=True)
    return 0 if all([fetch_one(out, k, a.parallel) for k in keys]) else 1


if __name__ == "__main__":
    sys.exit(main())
