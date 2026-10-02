"""Scenario simulator: every injected problem is really there and labelled."""
import hashlib

import numpy as np
import pytest

from grnbench.profile import technology_from_counts
from grnbench.simulate import PRESETS, Scenario, random_scenario, simulate, write_scenario


def _small(**kw):
    return Scenario(n_cells=400, n_genes=300, n_tfs=10, targets_per_tf=6, **kw)


def test_counts_are_integer_and_the_network_is_planted():
    adata, tfs, gold = simulate(_small())
    X = adata.X.toarray()
    assert np.allclose(X, np.round(X)) and (X >= 0).all()
    assert set(gold["source"]) <= set(tfs) and set(gold["target"]) <= set(adata.var_names)


def test_problems_are_labelled_per_barcode():
    adata, _, _ = simulate(_small(damaged_frac=0.1, empty_frac=0.1, doublet_frac=0.1, raw_droplets=500))
    counts = adata.obs["simulated_state"].value_counts()
    assert {"ok", "damaged", "empty", "doublet", "raw_empty"} <= set(counts.index)
    assert counts["raw_empty"] == 500


def test_full_length_counts_look_like_reads():
    adata, _, _ = simulate(_small(protocol="full_length"))
    assert technology_from_counts(adata.X)[0] == "full_length"
    umi = Scenario(n_cells=400, n_genes=1000, n_tfs=10, targets_per_tf=6)
    assert technology_from_counts(simulate(umi)[0].X)[0] == "umi"


def test_confounded_batches_follow_cell_type():
    adata, _, _ = simulate(_small(n_batches=3, batch_confounding=1.0, batch_shift_sd=0.3))
    tab = adata.obs.groupby("cell_type")["batch"].nunique()
    assert (tab == 1).all()


@pytest.mark.parametrize("ids,expect", [("mouse", "Tf01"), ("ensembl", "ENSG00000000001")])
def test_identifier_styles(ids, expect):
    adata, tfs, _ = simulate(_small(gene_ids=ids))
    assert tfs[0] == expect and adata.var_names[0] == expect


def test_misleading_readme_contradicts_the_protocol(tmp_path):
    paths = write_scenario(tmp_path, _small(readme="misleading"))
    assert "Smart-seq2" in paths["readme"].read_text()
    assert "readme" not in write_scenario(tmp_path / "b", _small(readme="none"))


def test_every_preset_and_random_scenario_is_valid():
    for scn in list(PRESETS.values()) + [random_scenario(i) for i in range(20)]:
        assert scn.damaged_frac + scn.empty_frac + scn.doublet_frac < 1
        assert isinstance(scn.expected(), dict)


def test_d019_default_simulation_checksum_is_stable():
    """Make future changes to the D-019 default random stream explicit."""
    adata, _, gold = simulate(PRESETS["small"].with_seed(0))
    digest = hashlib.sha256()
    digest.update(np.asarray(adata.X.toarray(), dtype="<f4").tobytes())
    for column in ("cell_type", "simulated_state"):
        digest.update("\0".join(adata.obs[column].astype(str)).encode())
    ordered = gold.sort_values(["source", "target", "sign"]).reset_index(drop=True)
    digest.update(ordered.to_csv(index=False, lineterminator="\n").encode())

    assert digest.hexdigest() == "848b924bbe9726a32270cd810023b6beba198d6794412ad436bac6e51be00870"


def test_context_specific_truth_has_declared_shared_fraction_and_fixed_degree():
    scn = _small(n_cell_types=3, shared_edge_fraction=0.5, type_balance_ratio=1.0)
    adata, tfs, union, truth = simulate(scn, return_context_truth=True)
    by_context = truth["by_context"]
    assert set(by_context) == set(adata.obs["cell_type"].astype(str).unique())
    assert len(truth["shared"]) == len(tfs) * 3       # round(6 * 0.5)
    shared = set(zip(truth["shared"].source, truth["shared"].target))
    context_pairs = []
    for frame in by_context.values():
        assert (frame.groupby("source").size() == scn.targets_per_tf).all()
        pairs = set(zip(frame.source, frame.target))
        assert shared <= pairs
        context_pairs.append(pairs)
    assert len(set.intersection(*context_pairs)) == len(shared)
    assert len(union) > len(next(iter(by_context.values())))


def test_context_truth_files_and_manifest_are_written(tmp_path):
    import json
    paths = write_scenario(tmp_path, _small(n_cell_types=2, shared_edge_fraction=0.5))
    manifest = json.loads(paths["truth_manifest"].read_text())
    assert manifest["split_policy"].startswith("split union edges first")
    assert set(manifest["contexts"]) == set(paths["gold_by_context"])
    assert all(path.exists() for path in paths["gold_by_context"].values())
    assert paths["gold_shared"].exists()


@pytest.mark.parametrize("kwargs", [
    {"shared_edge_fraction": -0.1}, {"shared_edge_fraction": 1.1},
    {"type_balance_ratio": 0.5}, {"n_cell_types": 0},
])
def test_invalid_context_truth_parameters_fail_loudly(kwargs):
    with pytest.raises(ValueError):
        _small(**kwargs)
