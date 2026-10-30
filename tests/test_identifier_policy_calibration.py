from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd

from benchmarks.identifier_policy_calibration import find_inputs, main


def _write(path, genes):
    ad.AnnData(
        X=np.ones((2, len(genes))),
        var=pd.DataFrame(index=genes),
    ).write_h5ad(path)


def test_calibration_distinguishes_policy_refusal_and_writes_portable_rows(tmp_path, capsys):
    root = tmp_path / "data"
    root.mkdir()
    _write(root / "symbols.h5ad", ["TP53", "MYC", "GAPDH"])
    _write(root / "ensembl.h5ad", [f"ENSG{i:011d}" for i in range(20)])
    table = tmp_path / "result.tsv"

    assert main([str(root), "--tsv", str(table)]) == 2

    output = capsys.readouterr().out
    written = table.read_text(encoding="utf-8")
    assert "symbols.h5ad\taccepted" in written
    assert "ensembl.h5ad\tpolicy_refused" in written
    assert "organism_scope\torganism_source\tcase_fold_action" in written
    assert "policy refused: 1" in output
    assert str(root) not in written


def test_calibration_refuses_an_empty_root(tmp_path):
    assert main([str(tmp_path)]) == 1


def test_input_discovery_ignores_non_10x_h5_cache_files(tmp_path):
    root = tmp_path / "data"
    root.mkdir()
    arbitrary = root / "internal-cache.h5"
    tenx = root / "raw_feature_bc_matrix.h5"
    h5ad = root / "sample.h5ad"
    for path in (arbitrary, tenx, h5ad):
        path.touch()

    assert find_inputs(root) == [tenx, h5ad]
