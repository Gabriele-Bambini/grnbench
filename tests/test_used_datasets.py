"""used_datasets.py: every CELLxGENE dataset and study earlier records and results used, so a new
pool can exclude them (D-057). The scan runs on this repository at HEAD; the listing is fake."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "benchmarks"))
import used_datasets  # noqa: E402

LISTING = [
    {"dataset_id": "0bcda669-aaaa-4aaa-8aaa-aaaaaaaaaaaa", "collection_id": "coll-ld1",
     "assets": [{"filetype": "H5AD", "url": "https://datasets.example.invalid/11111111-2222-4333-8444-555555555555.h5ad"}]},
    {"dataset_id": "6de332e1-bbbb-4bbb-8bbb-bbbbbbbbbbbb", "collection_id": "coll-ld1b", "assets": []},
    {"dataset_id": "12345678-cccc-4ccc-8ccc-cccccccccccc", "collection_id": "coll-new", "assets": []},
]


def test_resolve_finds_datasets_by_id_by_download_file_and_by_prefix():
    used, unresolved = used_datasets.resolve(
        {"6de332e1-bbbb-4bbb-8bbb-bbbbbbbbbbbb", "11111111-2222-4333-8444-555555555555", "12345678",
         "deadbeef", "99999999-9999-4999-8999-999999999999"}, LISTING)
    assert used == {"0bcda669-aaaa-4aaa-8aaa-aaaaaaaaaaaa", "6de332e1-bbbb-4bbb-8bbb-bbbbbbbbbbbb",
                    "12345678-cccc-4ccc-8ccc-cccccccccccc"}
    assert unresolved == {"deadbeef", "99999999-9999-4999-8999-999999999999"}


def test_scan_reads_folder_suffixes_in_results_and_bare_prefixes_in_exclusion_lists():
    everywhere = used_datasets.scan(["HEAD"])
    assert {"0bcda669", "6de332e1"} <= everywhere                    # EXP-20261001-LD-1 results
    records_only = used_datasets.scan(["HEAD"], paths=("experiments",))
    assert "0bcda669" in records_only                                 # bare line of an excluded_ids file


def test_scan_skips_the_given_paths():
    # the pool's own files are left out too, as the real call does with --skip
    assert "0bcda669" not in used_datasets.scan(
        ["HEAD"], paths=("experiments",),
        skip=("experiments/*.excluded_ids.txt", "experiments/EXP-20261022-GB*"))


def test_scan_refuses_an_unknown_ref():
    with pytest.raises(RuntimeError, match="no-such-ref"):
        used_datasets.scan(["no-such-ref"])


def test_cli_writes_ids_with_their_download_files_and_the_studies(tmp_path):
    (tmp_path / "listing.json").write_text(json.dumps(LISTING), encoding="utf-8")
    assert used_datasets.main(["--listing", str(tmp_path / "listing.json"), "--refs", "HEAD",
                               "--out-ids", str(tmp_path / "ids.txt"),
                               "--out-collections", str(tmp_path / "colls.txt")]) == 0
    ids = set((tmp_path / "ids.txt").read_text(encoding="utf-8").split())
    assert {"0bcda669-aaaa-4aaa-8aaa-aaaaaaaaaaaa", "11111111-2222-4333-8444-555555555555",
            "6de332e1-bbbb-4bbb-8bbb-bbbbbbbbbbbb"} <= ids
    assert "12345678-cccc-4ccc-8ccc-cccccccccccc" not in ids
    assert set((tmp_path / "colls.txt").read_text(encoding="utf-8").split()) == {"coll-ld1", "coll-ld1b"}
