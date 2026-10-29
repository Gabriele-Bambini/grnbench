"""benchmarks/knockdown_bench.py on a small synthetic screen with planted responders."""
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
import knockdown_bench as kb  # noqa: E402

N_GENES, N_CTRL, N_PER = 300, 500, 60
GENES = [f"G{i}" for i in range(N_GENES)]
PLANTED = GENES[10:25]            # 15 genes that respond to TFA


def _cells(rng, n, mult=None, kd=None):
    mean = np.full(N_GENES, 5.0)
    for g, m in (mult or {}).items():
        mean[GENES.index(g)] *= m
    for g, f in (kd or {}).items():
        mean[GENES.index(g)] *= f
    return rng.poisson(mean, size=(n, N_GENES)).astype(np.float32)


@pytest.fixture(scope="module")
def screen(tmp_path_factory):
    rng = np.random.default_rng(7)
    up = {g: 2.0 for g in PLANTED}
    blocks, rows = [], []

    def add(n, label, guide, **kw):
        blocks.append(_cells(rng, n, **kw))
        rows.extend([(label, guide)] * n)

    add(N_CTRL, "control", "ctl")
    add(N_PER, "G0", "G0_P1", mult=up, kd={"G0": 0.05})                       # TFA, two constructs, agree
    add(N_PER, "G0", "G0_P2", mult=up, kd={"G0": 0.05})
    add(N_PER, "G1", "G1_P1", mult=up)                                          # TFB: own gene not knocked down
    add(N_PER, "G2", "G2_P1", mult=up, kd={"G2": 0.05})                         # TFD: one construct
    both = {**{g: 2.0 for g in PLANTED}}
    add(N_PER, "G3", "G3_P1", mult=both, kd={"G3": 0.05})                       # TFE: constructs disagree
    add(N_PER, "G3", "G3_P2", mult={g: 1 / 2 for g in PLANTED}, kd={"G3": 0.05})
    obs = pd.DataFrame(rows, columns=["perturbation", "guide_id"], index=[f"c{i:05d}" for i in range(len(rows))])
    obs["nperts"] = 1
    path = tmp_path_factory.mktemp("screen") / "screen.h5ad"
    ad.AnnData(X=np.vstack(blocks), obs=obs, var=pd.DataFrame(index=GENES)).write_h5ad(path)
    tfs = tmp_path_factory.mktemp("ref") / "collectri.tsv"
    pd.DataFrame({"source": ["G0", "G1", "G2", "G3", "NOTMEASURED"], "target": "G50", "weight": 1}).to_csv(
        tfs, sep="\t", index=False)
    out = tmp_path_factory.mktemp("truth")
    assert kb.main(["truth", "--h5ad", str(path), "--tfs", str(tfs), "--out", str(out), "--chunk", "37"]) == 0
    return path, tfs, out


def test_planted_tf_is_usable_and_its_genes_are_found(screen):
    units = pd.read_csv(screen[2] / "units.tsv", sep="\t", keep_default_na=False).set_index("tf")
    resp = pd.read_csv(screen[2] / "responders.tsv", sep="\t")
    for tf in ("G0", "G2"):
        assert units.loc[tf, "usable"] and units.loc[tf, "knockdown"] > 0.9
        found = set(resp.gene[resp.tf == tf])
        assert set(PLANTED) <= found and len(found - set(PLANTED)) <= 10
        planted = resp[(resp.tf == tf) & resp.gene.isin(PLANTED)]
        assert (planted.log2fc > 0.5).all()
        # CP10k is compositional: the extras are mild decreases just past the 0.25 floor
        assert (resp[(resp.tf == tf) & ~resp.gene.isin(PLANTED)].log2fc < 0).all()
    assert units.loc["G0", "n_constructs"] == 2 and units.loc["G0", "sign_test"]
    assert units.loc["G2", "n_constructs"] == 1 and not units.loc["G2", "sign_test"]


def test_no_knockdown_and_unmeasured_tfs_are_excluded_with_the_reason(screen):
    units = pd.read_csv(screen[2] / "units.tsv", sep="\t", keep_default_na=False).set_index("tf")
    assert not units.loc["G1", "usable"] and units.loc["G1", "exclusion_reason"] == "knockdown_below_threshold"
    # a TF that is in the reference but not perturbed is not listed; a perturbed one absent from var is
    assert "NOTMEASURED" not in units.index
    ad_ = ad.read_h5ad(screen[0])
    ad_.obs["perturbation"] = ad_.obs.perturbation.astype(str)
    ad_.obs.loc[ad_.obs.perturbation == "G2", "perturbation"] = "NOTMEASURED"
    ctrl = kb.Controls(ad_.X, np.flatnonzero(ad_.obs.perturbation == "control"), 50)
    cells = np.flatnonzero(ad_.obs.perturbation == "NOTMEASURED")
    unit, resp = kb.tf_unit(ad_.X, ad_.var_names, kb.gene_lookup(ad_.var_names), ctrl, "NOTMEASURED", cells,
                            np.array(["x"] * len(cells)), 50)
    assert unit["exclusion_reason"] == "tf_not_measured" and not unit["usable"] and resp.empty


def test_constructs_that_disagree_in_sign_do_not_make_responders(screen):
    units = pd.read_csv(screen[2] / "units.tsv", sep="\t", keep_default_na=False).set_index("tf")
    assert units.loc["G3", "n_constructs"] == 2
    assert not units.loc["G3", "usable"] and units.loc["G3", "exclusion_reason"] == "too_few_responders"


def test_summary_counts(screen):
    s = kb.summarise(pd.read_csv(screen[2] / "units.tsv", sep="\t", keep_default_na=False))
    assert s["n_perturbed_tfs"] == 4 and s["n_usable"] == 2
    assert s["excluded_by_reason"] == {"knockdown_below_threshold": 1, "too_few_responders": 1}
    assert s["responders_among_usable"]["min"] >= 11
    meta = json.loads((screen[2] / "truth.json").read_text())
    assert "CP10k" in meta["normalisation"] and meta["control_split"]["seed"] == 0


def test_control_halves_are_disjoint_and_deterministic(screen):
    a = ad.read_h5ad(screen[0])
    ctl = np.flatnonzero(a.obs.perturbation == "control")
    i1, t1 = kb.split_controls(a.obs.index, ctl)
    i2, t2 = kb.split_controls(a.obs.index[::-1][::-1], ctl)
    assert not set(i1) & set(t1) and set(i1) | set(t1) == set(ctl) and abs(len(i1) - len(t1)) <= 1
    assert (i1 == i2).all() and (t1 == t2).all()
    assert set(kb.split_controls(a.obs.index, ctl, seed=1)[0]) != set(i1)


def test_controls_command_writes_inference_half_with_all_genes(screen, tmp_path):
    assert kb.main(["controls", "--h5ad", str(screen[0]), "--out", str(tmp_path), "--n", "100", "--chunk", "30"]) == 0
    c = ad.read_h5ad(tmp_path / "controls.h5ad")
    full = ad.read_h5ad(screen[0])
    inf, truth = kb.split_controls(full.obs.index, np.flatnonzero(full.obs.perturbation == "control"))
    assert c.n_obs == 100 and c.n_vars == N_GENES
    assert set(c.obs_names) <= set(full.obs_names[inf]) and not set(c.obs_names) & set(full.obs_names[truth])

    res = kb.baseline_check(tmp_path / "controls.h5ad", tmp_path / "work")
    assert res["error"] is None and res["genes_kept"] == N_GENES and res["cells_kept"] == 100


def test_reader_never_loads_all_rows(screen):
    a = ad.read_h5ad(screen[0])
    asked = []

    class Spy:
        def __getitem__(self, key):
            asked.append(len(key[0]))
            return a.X[key]

    n = a.n_obs
    ctrl = kb.Controls(Spy(), np.flatnonzero(a.obs.perturbation == "control"), 40)
    cells = np.flatnonzero(a.obs.perturbation == "G0")
    kb.tf_unit(Spy(), a.var_names, kb.gene_lookup(a.var_names), ctrl, "G0", cells,
               a.obs.guide_id.to_numpy()[cells], 40)
    assert asked and max(asked) <= 40 < n


def test_subset_reads_each_block_once_and_truth_on_it_is_byte_identical(screen, tmp_path):
    import h5py
    from anndata.io import write_elem
    path, tfs, _ = screen
    full = ad.read_h5ad(path)
    chunked = tmp_path / "chunked.h5ad"                 # dense, gzip, chunked like the K562 screen
    full.write_h5ad(chunked)
    with h5py.File(chunked, "a") as f:
        del f["X"]
        write_elem(f, "X", full.X, dataset_kwargs={"chunks": (37, 50), "compression": "gzip"})
    sub = tmp_path / "sub.h5ad"
    assert kb.main(["subset", "--h5ad", str(chunked), "--tfs", str(tfs), "--out", str(sub)]) == 0
    meta = json.loads((tmp_path / "sub.subset.json").read_text())
    assert meta["block_rows"] == 37 and meta["x_chunks"] == [37, 50] and meta["n_rows_in"] == full.n_obs
    s = ad.read_h5ad(sub)
    assert s.n_obs == meta["n_rows_out"] == full.n_obs and hasattr(s.X, "tocsr")    # every cell here is a control or a TF
    assert (s.obs.index == full.obs.index).all() and (s.X.toarray() == full.X).all()

    # a TF-less reference keeps the controls only, in the original order
    only = tmp_path / "only.tsv"
    pd.DataFrame({"source": ["G0"]}).to_csv(only, sep="\t", index=False)
    assert kb.main(["subset", "--h5ad", str(chunked), "--tfs", str(only), "--out", str(tmp_path / "g0.h5ad"),
                    "--block-rows", "74"]) == 0
    g0 = ad.read_h5ad(tmp_path / "g0.h5ad")
    want = full.obs.index[full.obs.perturbation.isin(["control", "G0"])]
    assert (g0.obs.index == want).all() and g0.n_obs == N_CTRL + 2 * N_PER

    asked = []

    class Spy:
        def __getitem__(self, key):
            asked.append((key[0].start, key[0].stop))
            return full.X[key]

    keep = (full.obs.perturbation != "G1").to_numpy()
    n = sum(c.shape[0] for c in kb.subset_blocks(Spy(), keep, 100))
    assert n == keep.sum() and asked == [(i, i + 100) for i in range(0, full.n_obs, 100)]   # each block once, in order

    for name, h5ad in (("full", chunked), ("sub", sub)):
        assert kb.main(["truth", "--h5ad", str(h5ad), "--tfs", str(tfs), "--out", str(tmp_path / name),
                        "--chunk", "37"]) == 0
    for f in ("responders.tsv", "units.tsv"):
        assert (tmp_path / "full" / f).read_bytes() == (tmp_path / "sub" / f).read_bytes()
    assert len((tmp_path / "sub" / "responders.tsv").read_text().splitlines()) > 1
    for name, h5ad in (("cfull", chunked), ("csub", sub)):
        assert kb.main(["controls", "--h5ad", str(h5ad), "--out", str(tmp_path / name), "--n", "100"]) == 0
    assert (ad.read_h5ad(tmp_path / "cfull" / "controls.h5ad").X == ad.read_h5ad(tmp_path / "csub" / "controls.h5ad").X).all()



# ---------------------------------------------------------------- --norm core
BLOCK = GENES[100:160]            # 60 genes a knockdown raises 6x: a compositional shift under CP10k


@pytest.fixture(scope="module")
def lib_screen(tmp_path_factory):
    """Cells with true library factors; TF G0's knockdown raises BLOCK strongly; G5 raises PLANTED 2x."""
    rng = np.random.default_rng(11)
    blocks, rows, fac = [], [], []

    def add(n, label, guide, mult=None, kd=None):
        f = rng.lognormal(0, 0.3, n)
        mean = np.full(N_GENES, 5.0)
        for g, m in {**(mult or {}), **(kd or {})}.items():
            mean[GENES.index(g)] *= m
        blocks.append(rng.poisson(f[:, None] * mean).astype(np.float32))
        rows.extend([(label, guide)] * n)
        fac.extend(f)

    add(N_CTRL, "control", "ctl")
    add(N_PER, "G0", "G0_P1", {**{g: 6.0 for g in BLOCK}, **{g: 2.0 for g in PLANTED}}, {"G0": 0.05})
    add(N_PER, "G5", "G5_P1", {g: 2.0 for g in PLANTED}, {"G5": 0.05})
    obs = pd.DataFrame(rows, columns=["perturbation", "guide_id"], index=[f"c{i:05d}" for i in range(len(rows))])
    obs["nperts"] = 1
    obs["core_scale_factor"] = fac
    d = tmp_path_factory.mktemp("lib")
    ad.AnnData(X=np.vstack(blocks), obs=obs, var=pd.DataFrame(index=GENES)).write_h5ad(d / "lib.h5ad")
    pd.DataFrame({"source": ["G0", "G5"], "target": "G50", "weight": 1}).to_csv(d / "tfs.tsv", sep="\t", index=False)
    return d


def _truth(d, out, *extra, h5ad="lib.h5ad"):
    assert kb.main(["truth", "--h5ad", str(d / h5ad), "--tfs", str(d / "tfs.tsv"), "--out", str(out),
                    "--chunk", "37", *extra]) == 0
    return pd.read_csv(out / "units.tsv", sep="\t", keep_default_na=False), pd.read_csv(out / "responders.tsv", sep="\t")


def test_default_norm_is_cp10k_and_byte_identical(screen, tmp_path, monkeypatch):
    path, tfs, _ = screen
    args = ["truth", "--h5ad", str(path), "--tfs", str(tfs), "--chunk", "37"]
    assert kb.main([*args, "--out", str(tmp_path / "default")]) == 0
    assert kb.main([*args, "--out", str(tmp_path / "explicit"), "--norm", "cp10k"]) == 0
    seen = []
    monkeypatch.setattr(kb, "normalise", lambda b, i, s: (seen.append(s), kb.cp10k(b))[1])   # the old call, verbatim
    assert kb.main([*args, "--out", str(tmp_path / "old")]) == 0
    assert seen and all(s is None for s in seen)
    for f in ("responders.tsv", "units.tsv"):
        assert (tmp_path / "default" / f).read_bytes() == (tmp_path / "explicit" / f).read_bytes() \
            == (tmp_path / "old" / f).read_bytes() == (screen[2] / f).read_bytes()
    meta = json.loads((tmp_path / "default" / "truth.json").read_text())
    assert meta["normalisation"] == kb.NORMALISATION and "norm" not in meta
    assert meta["rule"]["test"].endswith("on log1p CP10k, BH within TF over tested genes")


def test_core_norm_removes_compositional_decreases(lib_screen, tmp_path):
    u_cp, r_cp = _truth(lib_screen, tmp_path / "cp")
    u_core, r_core = _truth(lib_screen, tmp_path / "core", "--norm", "core")
    assert set(PLANTED) <= set(r_core.gene[r_core.tf == "G5"]) and set(PLANTED) <= set(r_core.gene[r_core.tf == "G0"])
    assert not set(PLANTED) <= set(r_cp.gene[r_cp.tf == "G0"])      # CP10k: the block's 6x swamps the 2x genes
    spurious = lambda r: int(((r.tf == "G0") & ~r.gene.isin(PLANTED + BLOCK) & (r.log2fc < 0)).sum())
    assert spurious(r_core) < spurious(r_cp) and spurious(r_cp) >= 20
    assert set(BLOCK) <= set(r_core.gene[r_core.tf == "G0"])
    meta = json.loads((tmp_path / "core" / "truth.json").read_text())
    assert meta["normalisation"] == kb.NORMALISATION_CORE and meta["norm_obs_column"] == "core_scale_factor"
    assert meta["n_cells_excluded_no_factor"] == 0
    assert "log1p counts / core_scale_factor" in meta["rule"]["test"] and "CP10k" not in meta["rule"]["test"]


def test_missing_or_zero_factor_excludes_the_cell_and_counts_it(lib_screen, tmp_path):
    a = ad.read_h5ad(lib_screen / "lib.h5ad")
    g0 = np.flatnonzero(a.obs.perturbation == "G0")
    ctl = np.flatnonzero(a.obs.perturbation == "control")
    a.obs["core_scale_factor"] = a.obs.core_scale_factor.astype(float)
    a.obs.iloc[g0[:3], a.obs.columns.get_loc("core_scale_factor")] = np.nan
    a.obs.iloc[g0[3:5], a.obs.columns.get_loc("core_scale_factor")] = 0.0
    a.obs.iloc[ctl[:4], a.obs.columns.get_loc("core_scale_factor")] = -1.0
    a.write_h5ad(lib_screen / "holes.h5ad")
    units, _ = _truth(lib_screen, tmp_path / "holes", "--norm", "core", h5ad="holes.h5ad")
    assert int(units.set_index("tf").loc["G0", "n_cells_total"]) == N_PER - 5
    meta = json.loads((tmp_path / "holes" / "truth.json").read_text())
    assert meta["n_cells_excluded_no_factor"] == 9 and meta["control_split"]["n_controls"] == N_CTRL - 4
    with pytest.raises(SystemExit):                        # no factor column at all
        kb.core_scale(pd.DataFrame(index=["a"]), "core")


def test_core_norm_through_subset_keeps_the_factor_column(lib_screen, tmp_path):
    sub = tmp_path / "sub.h5ad"
    assert kb.main(["subset", "--h5ad", str(lib_screen / "lib.h5ad"), "--tfs", str(lib_screen / "tfs.tsv"),
                    "--out", str(sub)]) == 0
    assert "core_scale_factor" in ad.read_h5ad(sub).obs.columns
    _truth(lib_screen, tmp_path / "full", "--norm", "core")
    assert kb.main(["truth", "--h5ad", str(sub), "--tfs", str(lib_screen / "tfs.tsv"), "--out", str(tmp_path / "viasub"),
                    "--chunk", "37", "--norm", "core"]) == 0
    for f in ("responders.tsv", "units.tsv"):
        assert (tmp_path / "full" / f).read_bytes() == (tmp_path / "viasub" / f).read_bytes()


# ---------------------------------------------------------------- total-shift (check P9 sensitivity)
@pytest.fixture(scope="module")
def shift_screen(tmp_path_factory):
    """Controls at 1000 core-adjusted UMIs; G0 at 800 (median of 700/800/900), G1 the same via a factor of 2, G2 at 1.1."""
    d = tmp_path_factory.mktemp("shift")
    cells = [("control", 1000.0, 1.0)] * 20 + [("G0", u, 1.0) for u in (700, 800, 900)] * 5 \
        + [("G1", 2000.0, 2.0)] * 10 + [("G2", 1100.0, 1.0)] * 10 + [("G9", 100.0, 1.0)] * 10
    obs = pd.DataFrame(cells, columns=["perturbation", "UMI_count", "core_scale_factor"],
                       index=[f"c{i:03d}" for i in range(len(cells))])
    obs["nperts"] = 1
    ad.AnnData(X=np.ones((len(cells), 3), dtype=np.float32), obs=obs,
               var=pd.DataFrame(index=["a", "b", "c"])).write_h5ad(d / "s.h5ad")
    pd.DataFrame({"tf": ["G2", "G0", "G1", "G9"], "usable": [True, True, True, False]}).to_csv(
        d / "units.tsv", sep="\t", index=False)
    pd.DataFrame({"source": ["G0", "G1", "G2", "G9"], "target": "a", "weight": 1}).to_csv(d / "tfs.tsv", sep="\t", index=False)
    return d


def test_total_shift_matches_a_hand_computation_and_the_subset_copy(shift_screen, tmp_path):
    d = shift_screen
    out = tmp_path / "ts.tsv"
    assert kb.main(["total-shift", "--h5ad", str(d / "s.h5ad"), "--truth", str(d), "--out", str(out)]) == 0
    t = pd.read_csv(out, sep="\t")
    assert list(t.columns) == ["tf", "n_cells", "rel_core_adjusted_total", "within_0_9_1_1"]
    assert list(t.tf) == ["G0", "G1", "G2"] and list(t.n_cells) == [15, 10, 10]        # sorted, usable only
    assert t.rel_core_adjusted_total.tolist() == pytest.approx([0.8, 1.0, 1.1])
    assert t.within_0_9_1_1.tolist() == [False, True, True]                               # bounds inclusive
    sub = tmp_path / "sub.h5ad"
    assert kb.main(["subset", "--h5ad", str(d / "s.h5ad"), "--tfs", str(d / "tfs.tsv"), "--out", str(sub)]) == 0
    assert "UMI_count" in ad.read_h5ad(sub).obs.columns
    out2 = tmp_path / "ts_sub.tsv"
    assert kb.main(["total-shift", "--h5ad", str(sub), "--truth", str(d), "--out", str(out2)]) == 0
    assert out.read_bytes() == out2.read_bytes()


def test_total_shift_fails_loudly_without_its_columns_or_controls(shift_screen, tmp_path):
    d = shift_screen
    a = ad.read_h5ad(d / "s.h5ad")
    ad.AnnData(X=a.X, obs=a.obs.drop(columns="UMI_count"), var=a.var).write_h5ad(tmp_path / "noumi.h5ad")
    with pytest.raises(SystemExit, match="UMI_count"):
        kb.main(["total-shift", "--h5ad", str(tmp_path / "noumi.h5ad"), "--truth", str(d), "--out", str(tmp_path / "x.tsv")])
    with pytest.raises(SystemExit, match="no cell with"):
        kb.main(["total-shift", "--h5ad", str(d / "s.h5ad"), "--truth", str(d), "--out", str(tmp_path / "x.tsv"),
                 "--control-label", "nonesuch"])
    assert not (tmp_path / "x.tsv").exists()


# ---------------------------------------------------------------- controls and truth cut the same halves
def test_controls_and_truth_halves_are_disjoint_even_with_a_control_lacking_a_factor(lib_screen, tmp_path):
    a = ad.read_h5ad(lib_screen / "lib.h5ad")
    ctl = np.flatnonzero(a.obs.perturbation == "control")
    a.obs["core_scale_factor"] = a.obs.core_scale_factor.astype(float)
    a.obs.iloc[ctl[7], a.obs.columns.get_loc("core_scale_factor")] = np.nan
    a.write_h5ad(lib_screen / "nan.h5ad")
    bad = a.obs_names[ctl[7]]
    _truth(lib_screen, tmp_path / "truth", "--norm", "core", h5ad="nan.h5ad")
    assert kb.main(["controls", "--h5ad", str(lib_screen / "nan.h5ad"), "--out", str(tmp_path / "ctl"), "--norm", "core"]) == 0
    inference = set(ad.read_h5ad(tmp_path / "ctl" / "controls.h5ad").obs_names)
    truth_half = (tmp_path / "truth" / "truth_half.txt").read_text().split()
    assert truth_half == sorted(truth_half) and not inference & set(truth_half)
    assert len(inference) + len(truth_half) == N_CTRL - 1 and bad not in inference | set(truth_half)
    assert ad.read_h5ad(tmp_path / "ctl" / "controls.h5ad").uns["source_h5ad"] == "nan.h5ad"
    # no --n: the whole inference half
    assert kb.main(["controls", "--h5ad", str(lib_screen / "nan.h5ad"), "--out", str(tmp_path / "cp")]) == 0
    assert ad.read_h5ad(tmp_path / "cp" / "controls.h5ad").n_obs == N_CTRL // 2     # no cap by default
