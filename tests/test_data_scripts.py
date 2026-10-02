"""Data provenance scripts used to put a registered experiment's data on another machine:
fetch_cellxgene --from-manifest (same datasets, no re-selection) and data_manifest --check
(same bytes, or the same text with other line endings)."""
from __future__ import annotations

import hashlib
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


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
