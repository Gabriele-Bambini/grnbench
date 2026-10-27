"""Data provenance scripts used to put a registered experiment's data on another machine:
fetch_cellxgene --from-manifest (same datasets, no re-selection, same README.md) and
data_manifest --check (same bytes, or the same text with other line endings)."""
from __future__ import annotations

import hashlib
import importlib.util
import io
import json
import subprocess
import sys
from pathlib import Path
from urllib.error import URLError

import pytest

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("fetch_cellxgene_data_scripts",
                                               ROOT / "benchmarks" / "fetch_cellxgene.py")
fx = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fx)

BODY = b"not really an h5ad"
H5AD_URL = "https://example.invalid/4ec57620.h5ad"
LISTED = {"dataset_id": "4ec57620-0000-0000-0000-000000000000", "title": "Mouse M1 nuclei",
          "collection_name": "Evolution of cellular diversity in primary motor cortex",
          "citation": "Publication: https://doi.org/10.0000/example",
          "assay": [{"label": "10x 3' v3"}], "organism": [{"label": "Mus musculus"}],
          "tissue": [{"label": "primary motor cortex"}], "suspension_type": ["nucleus"],
          "cell_count": 1000, "assets": [{"filetype": "H5AD", "filesize": len(BODY), "url": H5AD_URL}]}


@pytest.fixture
def cellxgene(monkeypatch):
    """Fake CELLxGENE: the listing (API) serves ``listing`` (None = unreachable), any other URL BODY."""
    state = {"listing": [LISTED], "calls": []}

    def urlopen(req, timeout):
        state["calls"].append(req.full_url)
        if req.full_url == fx.API:
            if state["listing"] is None:
                raise URLError("no route to host")
            return io.BytesIO(json.dumps(state["listing"]).encode())
        return io.BytesIO(BODY)

    monkeypatch.setattr(fx, "urlopen", urlopen)
    monkeypatch.setattr(fx.time, "sleep", lambda s: None)
    return state


def _sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def _check(root: Path, manifest: Path, *, require_all: bool = False):
    command = [sys.executable, str(ROOT / "benchmarks" / "data_manifest.py"), str(root),
               "--check", str(manifest)]
    if require_all:
        command.append("--require-all")
    return subprocess.run(command, capture_output=True, text=True)


def _manifest(path: Path, rows):
    path.write_text("path\tbytes\tsha256\tfile_date\tsource\tversion\n"
                    + "".join(f"{p}\t0\t{s}\t-\t-\t-\n" for p, s in rows), encoding="utf-8")


def test_check_passes_identical_files_and_crlf_only_text(tmp_path):
    data = tmp_path / "data"
    (data / "cellxgene" / "x").mkdir(parents=True)
    (data / "cellxgene" / "x" / "data.h5ad").write_bytes(b"\x00binary\r\n")
    (data / "cellxgene" / "x" / "truth.json").write_bytes(b'{\n "organism": "human"\n}')
    _manifest(tmp_path / "m.tsv", [("cellxgene/x/data.h5ad", _sha(b"\x00binary\r\n")),
                                   ("cellxgene/x/truth.json", _sha(b'{\r\n "organism": "human"\r\n}'))])
    r = _check(data, tmp_path / "m.tsv")
    assert r.returncode == 0, r.stdout + r.stderr
    assert "ok=1 ok-eol=1 MISMATCH=0 UNLISTED=0" in r.stdout


def test_check_fails_on_changed_binary_and_unlisted_files(tmp_path):
    data = tmp_path / "data"
    (data / "cellxgene" / "x").mkdir(parents=True)
    (data / "cellxgene" / "x" / "data.h5ad").write_bytes(b"\x00binary\n")
    (data / "cellxgene" / "x" / "extra.h5ad").write_bytes(b"new")
    # a binary file never gets the line-ending allowance
    _manifest(tmp_path / "m.tsv", [("cellxgene/x/data.h5ad", _sha(b"\x00binary\r\n"))])
    r = _check(data, tmp_path / "m.tsv")
    assert r.returncode == 1
    assert "MISMATCH\tcellxgene/x/data.h5ad" in r.stdout
    assert "UNLISTED\tcellxgene/x/extra.h5ad" in r.stdout


def test_check_require_all_reports_registered_files_missing_from_root(tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    (data / "present.h5ad").write_bytes(b"present")
    _manifest(tmp_path / "m.tsv", [("present.h5ad", _sha(b"present")),
                                   ("missing.h5ad", _sha(b"missing"))])

    subset = _check(data, tmp_path / "m.tsv")
    assert subset.returncode == 0, subset.stdout + subset.stderr
    assert "MISSING" not in subset.stdout

    exact = _check(data, tmp_path / "m.tsv", require_all=True)
    assert exact.returncode == 1
    assert "MISSING\tmissing.h5ad" in exact.stdout
    assert "ok=1 ok-eol=0 MISMATCH=0 UNLISTED=0 MISSING=1" in exact.stdout


def test_check_cli_flags_are_order_independent(tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    (data / "present.h5ad").write_bytes(b"present")
    manifest = tmp_path / "m.tsv"
    _manifest(manifest, [("present.h5ad", _sha(b"present"))])
    script = str(ROOT / "benchmarks" / "data_manifest.py")

    commands = [
        [sys.executable, script, str(data), "--check", str(manifest), "--require-all"],
        [sys.executable, script, "--require-all", "--check", str(manifest), str(data)],
    ]
    for command in commands:
        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode == 0, result.stdout + result.stderr
        assert "ok=1 ok-eol=0 MISMATCH=0 UNLISTED=0 MISSING=0" in result.stdout


def test_fetch_from_manifest_keeps_the_recorded_datasets(tmp_path):
    cols = ["folder", "dataset_id", "assay", "organism", "tissue", "suspension_type", "cell_count",
            "filesize", "title", "url"]
    row = ["10x-3-v3/homo-sapiens/lung__abcdef12", "abcdef12-0000", "10x 3' v3", "Homo sapiens",
           "lung", "cell", "1000", "123", "A title", "https://example.invalid/a.h5ad"]
    man = tmp_path / "in.tsv"
    man.write_text("\t".join(cols) + "\n" + "\t".join(row) + "\n", encoding="utf-8")
    out = tmp_path / "out"
    r = subprocess.run([sys.executable, str(ROOT / "benchmarks" / "fetch_cellxgene.py"), "--out", str(out),
                        "--from-manifest", str(man), "--dry-run"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "1 datasets" in r.stdout
    written = (out / "manifest.tsv").read_text(encoding="utf-8").splitlines()
    assert written[1].split("\t") == row


def test_fetch_from_manifest_writes_the_readme_of_the_listing_route(tmp_path, cellxgene):
    """The engine reads README.md as a text source (#156: C6's n_warnings 4 -> 3), so a rerun
    fetched with --from-manifest must write the same README as the original listing fetch."""
    first, rerun = tmp_path / "first", tmp_path / "rerun"
    assert fx.main(["--out", str(first), "--max", "1"]) == 0
    assert fx.main(["--out", str(rerun), "--from-manifest", str(first / "manifest.tsv")]) == 0
    folder = "10x-3-v3/mus-musculus/primary-motor-cortex__4ec57620"
    readme = (first / folder / "README.md").read_text(encoding="utf-8")
    assert readme == f"{LISTED['title']}\n\nCollection: {LISTED['collection_name']}\n{LISTED['citation']}\n"
    assert (rerun / folder / "README.md").read_text(encoding="utf-8") == readme
    assert (rerun / folder / "data.h5ad").read_bytes() == BODY
    assert (rerun / folder / "truth.json").read_bytes() == (first / folder / "truth.json").read_bytes()


@pytest.mark.parametrize("listing", [[], None, {"error": "busy"}],
                         ids=["dataset-not-listed", "listing-unreachable", "listing-malformed"])
def test_fetch_from_manifest_without_the_listing_still_downloads_and_says_so(tmp_path, cellxgene,
                                                                            capsys, listing):
    assert fx.main(["--out", str(tmp_path / "first"), "--max", "1"]) == 0
    cellxgene["listing"] = listing
    capsys.readouterr()
    rerun = tmp_path / "rerun"
    assert fx.main(["--out", str(rerun), "--from-manifest", str(tmp_path / "first" / "manifest.tsv")]) == 0
    folder = rerun / "10x-3-v3/mus-musculus/primary-motor-cortex__4ec57620"
    assert (folder / "data.h5ad").read_bytes() == BODY
    assert (folder / "README.md").read_text(encoding="utf-8") == f"{LISTED['title']}\n\nCollection: \n\n"
    err = capsys.readouterr().err.splitlines()
    assert len(err) == 1, err
    assert LISTED["dataset_id"] in err[0] and "README.md" in err[0] and "planner" in err[0]


def test_fetch_from_manifest_dry_run_does_not_read_the_listing(tmp_path, cellxgene):
    assert fx.main(["--out", str(tmp_path / "first"), "--max", "1", "--dry-run"]) == 0
    cellxgene["calls"].clear()
    assert fx.main(["--out", str(tmp_path / "rerun"), "--from-manifest",
                    str(tmp_path / "first" / "manifest.tsv"), "--dry-run"]) == 0
    assert cellxgene["calls"] == []
    assert (tmp_path / "rerun" / "manifest.tsv").read_text(encoding="utf-8") == \
        (tmp_path / "first" / "manifest.tsv").read_text(encoding="utf-8")


def test_lm2_dataset_list_matches_the_provenance_manifest():
    """Every LM-2 dataset has a data.h5ad and a truth.json row with a sha256 in datasets.tsv,
    which is what data_manifest --check verifies on a new machine."""
    import csv
    with open(ROOT / "experiments" / "EXP-20260925-LM-2.datasets.tsv", encoding="utf-8") as fh:
        folders = [r["folder"] for r in csv.DictReader(fh, delimiter="\t")]
    with open(ROOT / "experiments" / "datasets.tsv", encoding="utf-8") as fh:
        paths = {r["path"]: r["sha256"] for r in csv.DictReader(fh, delimiter="\t")}
    assert len(folders) == 20
    for f in folders:
        for name in ("data.h5ad", "truth.json"):
            assert len(paths[f"cellxgene/{f}/{name}"]) == 64
