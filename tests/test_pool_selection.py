"""Selection of a large never-used dataset pool (fetch_cellxgene --one-per-collection, --max-cells)
and its seeded split into a confirmatory and an exploratory part (split_pool.py). Fake listings
only; nothing is downloaded."""
from __future__ import annotations

import csv
import io
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "benchmarks"))
import fetch_cellxgene  # noqa: E402
import split_pool  # noqa: E402


def _ds(i, coll, assay="10x 3' v3", org="Homo sapiens", mb=10, cells=2000, tissue="lung"):
    return {"dataset_id": f"{i:08x}-0000-0000-0000-000000000000", "collection_id": coll,
            "collection_name": coll, "title": f"t{i}", "citation": "", "assay": [{"label": assay}],
            "organism": [{"label": org}], "tissue": [{"label": tissue}], "suspension_type": ["cell"],
            "cell_count": cells,
            "assets": [{"filetype": "H5AD", "filesize": int(mb * 1e6), "url": f"https://example.invalid/{i}.h5ad"}]}


LISTING = [_ds(1, "A", mb=30), _ds(2, "A", mb=5), _ds(3, "B", assay="Smart-seq2", cells=800),
           _ds(4, "B", mb=50, tissue="brain"), _ds(5, "C", org="Mus musculus", mb=20, tissue="kidney"),
           _ds(6, "D", mb=40, cells=120_000, tissue="heart"), _ds(7, "E", mb=900, cells=300_000, tissue="liver"),
           _ds(8, "F", assay="10x 5' v1", mb=12, tissue="blood"), _ds(9, "G", mb=8, cells=400)]


def _ids(picked):
    return [d["dataset_id"][:8] for d in picked]


def test_default_selection_is_unchanged():
    # recorded from the selection before --one-per-collection and --max-cells existed
    assert _ids(fetch_cellxgene.select(LISTING, 10, 100)) == [
        "00000002", "00000005", "00000008", "00000003", "00000004", "00000001"]


def test_one_per_collection_keeps_the_smallest_file_of_each_collection():
    picked = fetch_cellxgene.select(LISTING, 10, 100, one_per_collection=True)
    assert sorted(_ids(picked)) == ["00000002", "00000003", "00000005", "00000008"]
    assert len({d["collection_id"] for d in picked}) == len(picked)


def test_max_cells_admits_larger_datasets_and_keeps_the_500_cell_floor():
    picked = fetch_cellxgene.select(LISTING, 10, 1000, max_cells=500_000, one_per_collection=True)
    assert sorted(_ids(picked)) == ["00000002", "00000003", "00000005", "00000006", "00000007", "00000008"]


def test_one_per_collection_refuses_a_dataset_without_a_collection():
    listing = LISTING + [_ds(10, "", mb=3)]
    with pytest.raises(ValueError, match="collection"):
        fetch_cellxgene.select(listing, 10, 100, one_per_collection=True)


def test_cli_dry_run_writes_the_collection_of_each_pick(tmp_path, monkeypatch):
    class _Response(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    monkeypatch.setattr(fetch_cellxgene, "urlopen", lambda *a, **k: _Response(json.dumps(LISTING).encode()))
    assert fetch_cellxgene.main(["--out", str(tmp_path), "--max", "10", "--max-mb", "1000",
                                 "--max-cells", "500000", "--one-per-collection", "--dry-run"]) == 0
    with open(tmp_path / "manifest.tsv", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 6 and {r["collection_id"] for r in rows} == {"A", "B", "C", "D", "E", "F"}


def test_cli_excludes_every_dataset_of_a_listed_collection(tmp_path, monkeypatch):
    class _Response(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    monkeypatch.setattr(fetch_cellxgene, "urlopen", lambda *a, **k: _Response(json.dumps(LISTING).encode()))
    (tmp_path / "colls.txt").write_text("A\nE\n", encoding="utf-8")
    assert fetch_cellxgene.main(["--out", str(tmp_path), "--max", "10", "--max-mb", "1000", "--max-cells", "500000",
                                 "--one-per-collection", "--exclude-collections", str(tmp_path / "colls.txt"),
                                 "--dry-run"]) == 0
    with open(tmp_path / "manifest.tsv", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert {r["collection_id"] for r in rows} == {"B", "C", "D", "F"}


def _fake_api(monkeypatch):
    class _Response(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    monkeypatch.setattr(fetch_cellxgene, "urlopen", lambda *a, **k: _Response(json.dumps(LISTING).encode()))


def test_an_exclusion_file_written_with_a_bom_still_excludes_its_first_entry(tmp_path, monkeypatch, capsys):
    _fake_api(monkeypatch)
    (tmp_path / "colls.txt").write_text("﻿A\r\nE\r\n", encoding="utf-8")   # PowerShell 5 Out-File
    assert fetch_cellxgene.main(["--out", str(tmp_path), "--max", "10", "--max-mb", "1000", "--max-cells", "500000",
                                 "--one-per-collection", "--exclude-collections", str(tmp_path / "colls.txt"),
                                 "--dry-run"]) == 0
    with open(tmp_path / "manifest.tsv", encoding="utf-8") as fh:
        assert {r["collection_id"] for r in csv.DictReader(fh, delimiter="\t")} == {"B", "C", "D", "F"}
    assert "--exclude-collections: 3 datasets dropped" in capsys.readouterr().out


@pytest.mark.parametrize("flag, text", [("--exclude-ids", "0\n"), ("--exclude-ids", "\n"), ("--exclude-collections", " \n")])
def test_an_empty_or_malformed_exclusion_file_is_refused(tmp_path, monkeypatch, flag, text):
    _fake_api(monkeypatch)
    (tmp_path / "x.txt").write_text(text, encoding="utf-8")
    with pytest.raises(SystemExit):
        fetch_cellxgene.main(["--out", str(tmp_path), flag, str(tmp_path / "x.txt"), "--dry-run"])


def test_cli_selects_from_a_saved_listing_without_the_network(tmp_path, monkeypatch):
    def _no_network(*a, **k):
        raise AssertionError("the API must not be called when --listing is given")

    monkeypatch.setattr(fetch_cellxgene, "urlopen", _no_network)
    (tmp_path / "listing.json").write_text(json.dumps(LISTING), encoding="utf-8")
    assert fetch_cellxgene.main(["--out", str(tmp_path), "--listing", str(tmp_path / "listing.json"), "--max", "10",
                                 "--max-mb", "1000", "--max-cells", "500000", "--one-per-collection",
                                 "--dry-run"]) == 0
    with open(tmp_path / "manifest.tsv", encoding="utf-8") as fh:
        assert len(list(csv.DictReader(fh, delimiter="\t"))) == 6


def test_from_manifest_refuses_selection_options(tmp_path):
    (tmp_path / "in.tsv").write_text("folder\tdataset_id\n", encoding="utf-8")
    r = subprocess.run([sys.executable, str(ROOT / "benchmarks" / "fetch_cellxgene.py"), "--out", str(tmp_path / "out"),
                        "--from-manifest", str(tmp_path / "in.tsv"), "--one-per-collection", "--dry-run"],
                       capture_output=True, text=True)
    assert r.returncode == 2 and "--from-manifest" in r.stderr


def test_split_refuses_to_overwrite_its_input(tmp_path):
    (tmp_path / "pool.tsv").write_text("folder\tdataset_id\nx/a\tid-1\nx/b\tid-2\n", encoding="utf-8")
    r = subprocess.run([sys.executable, str(ROOT / "benchmarks" / "split_pool.py"), "--manifest", str(tmp_path / "pool.tsv"),
                        "--seed", "42", "--confirmatory", "1", "--out-confirmatory", str(tmp_path / "pool.tsv"),
                        "--out-exploratory", str(tmp_path / "e.tsv")], capture_output=True, text=True)
    assert r.returncode != 0
    assert (tmp_path / "pool.tsv").read_text(encoding="utf-8").count("\n") == 3


def test_from_manifest_keeps_a_title_that_starts_with_a_quote(tmp_path):
    cols = ["folder", "dataset_id", "assay", "organism", "tissue", "suspension_type", "cell_count",
            "filesize", "title", "url", "collection_id"]
    row = ["10x-3-v3/homo-sapiens/lung__abcdef12", "abcdef12-0000", "10x 3' v3", "Homo sapiens", "lung", "cell",
           "1000", "123", '"Quoted" start', "https://example.invalid/a.h5ad", "c1"]
    (tmp_path / "in.tsv").write_text("\t".join(cols) + "\n" + "\t".join(row) + "\n", encoding="utf-8")
    r = subprocess.run([sys.executable, str(ROOT / "benchmarks" / "fetch_cellxgene.py"), "--out", str(tmp_path / "out"),
                        "--from-manifest", str(tmp_path / "in.tsv"), "--dry-run"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert (tmp_path / "out" / "manifest.tsv").read_text(encoding="utf-8").splitlines()[1].split("\t") == row


def test_split_refuses_a_row_with_missing_columns(tmp_path):
    (tmp_path / "pool.tsv").write_text("folder\tdataset_id\tassay\nx/a\tid-1\t10x\nx/b\tid-2\n", encoding="utf-8")
    r = subprocess.run([sys.executable, str(ROOT / "benchmarks" / "split_pool.py"), "--manifest", str(tmp_path / "pool.tsv"),
                        "--seed", "42", "--confirmatory", "1", "--out-confirmatory", str(tmp_path / "c.tsv"),
                        "--out-exploratory", str(tmp_path / "e.tsv")], capture_output=True, text=True)
    assert r.returncode != 0 and "missing" in (r.stderr + r.stdout)
    assert not (tmp_path / "c.tsv").exists()


def _pool(n):
    return [{"folder": f"x/y/t__{i:08x}", "dataset_id": f"{i:08x}-1111", "collection_id": f"c{i}"} for i in range(n)]


def test_split_is_seeded_and_does_not_depend_on_row_order():
    pool = _pool(30)
    c1, e1 = split_pool.split(pool, seed=42, n_confirmatory=18)
    c2, e2 = split_pool.split(list(reversed(pool)), seed=42, n_confirmatory=18)
    assert [r["dataset_id"] for r in c1] == [r["dataset_id"] for r in c2]
    assert len(c1) == 18 and len(e1) == 12
    assert {r["dataset_id"] for r in c1} | {r["dataset_id"] for r in e1} == {r["dataset_id"] for r in pool}
    assert not {r["dataset_id"] for r in c1} & {r["dataset_id"] for r in e1}
    c3, _ = split_pool.split(pool, seed=7, n_confirmatory=18)
    assert {r["dataset_id"] for r in c3} != {r["dataset_id"] for r in c1}


@pytest.mark.parametrize("pool, n", [(_pool(5), 6), (_pool(5) + _pool(1), 3), (_pool(5), 0)])
def test_split_refuses_impossible_or_ambiguous_requests(pool, n):
    with pytest.raises(ValueError):
        split_pool.split(pool, seed=42, n_confirmatory=n)


def test_split_cli_writes_two_manifests_with_the_input_columns(tmp_path):
    cols = ["folder", "dataset_id", "assay", "collection_id"]
    lines = ["\t".join(cols)] + [f"x/y/t__{i:08x}\t{i:08x}-1111\t10x\tc{i}" for i in range(10)]
    (tmp_path / "pool.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    r = subprocess.run([sys.executable, str(ROOT / "benchmarks" / "split_pool.py"), "--manifest", str(tmp_path / "pool.tsv"),
                        "--seed", "42", "--confirmatory", "6", "--out-confirmatory", str(tmp_path / "c.tsv"),
                        "--out-exploratory", str(tmp_path / "e.tsv")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    c = (tmp_path / "c.tsv").read_text(encoding="utf-8").splitlines()
    e = (tmp_path / "e.tsv").read_text(encoding="utf-8").splitlines()
    assert c[0].split("\t") == cols and e[0].split("\t") == cols
    assert len(c) - 1 == 6 and len(e) - 1 == 4
    assert "6 confirmatory, 4 exploratory" in r.stdout
