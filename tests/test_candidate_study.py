"""benchmarks/candidate_study.py: degradations, the copy path, task listing and the report."""
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
import candidate_study as cs  # noqa: E402


@pytest.fixture()
def counts():
    rng = np.random.default_rng(1)
    n, g = 600, 80
    names = [f"G{i}" for i in range(g - 4)] + ["MT-CO1", "MT-ND1", "mt-Atp6", "MT-CYB"]
    X = rng.poisson(rng.gamma(2.0, 1.5, size=g), size=(n, g)).astype(np.float32)
    obs = pd.DataFrame({"cell_type": rng.choice(["a", "b", "c"], n)}, index=[f"c{i}" for i in range(n)])
    return ad.AnnData(sp.csr_matrix(X), obs=obs, var=pd.DataFrame(index=names))


def _total(a):
    return float(a.X.sum())


def _is_int(a):
    return np.array_equal(a.X.data, np.rint(a.X.data))


def test_clean_is_identity_and_all_conditions_keep_integers(counts):
    clean = cs.degrade(counts, "clean")
    assert (clean.X != sp.csr_matrix(counts.X)).nnz == 0 and clean.n_obs == counts.n_obs
    for cond in cs.CONDITIONS:
        out = cs.degrade(counts, cond, seed=3)
        assert _is_int(out) and out.X.min() >= 0, cond
        assert out.n_vars == counts.n_vars


def test_unknown_condition_is_refused(counts):
    with pytest.raises(ValueError):
        cs.degrade(counts, "nope")


@pytest.mark.parametrize("cond,keep", [("depth_30", 0.3), ("depth_10", 0.1)])
def test_depth_thinning_keeps_the_stated_fraction(counts, cond, keep):
    out = cs.degrade(counts, cond, seed=0)
    assert out.n_obs == counts.n_obs
    assert _total(out) / _total(counts) == pytest.approx(keep, rel=0.05)
    assert (out.X.toarray() <= counts.X.toarray()).all()


def test_ambient_adds_about_twenty_percent_in_the_global_profile(counts):
    out = cs.degrade(counts, "ambient_20", seed=0)
    assert _total(out) / _total(counts) == pytest.approx(1.2, rel=0.03)
    assert (out.X.toarray() >= counts.X.toarray()).all()
    added = out.X.toarray().sum(0) - counts.X.toarray().sum(0)
    profile = counts.X.toarray().sum(0)
    assert np.corrcoef(added, profile)[0, 1] > 0.95


def test_doublets_append_cells_without_a_marker_column(counts):
    out = cs.degrade(counts, "doublets_20", seed=0)
    k = round(0.2 * counts.n_obs)
    assert out.n_obs == counts.n_obs + k
    assert list(out.obs.columns) == list(counts.obs.columns)
    new = out[counts.n_obs:]
    assert all("-dbl" in n for n in new.obs_names) and out.obs_names.is_unique
    # each doublet is the sum of two distinct original cells, mostly of different types
    X0, Xd = counts.X.toarray(), new.X.toarray()
    diff_type = 0
    for j, name in enumerate(new.obs_names):
        first = counts.obs_names.get_loc(name.rsplit("-dbl", 1)[0])
        rest = Xd[j] - X0[first]
        match = np.flatnonzero((X0 == rest).all(1))
        match = match[match != first]
        assert len(match) >= 1
        diff_type += any(counts.obs["cell_type"].iloc[m] != counts.obs["cell_type"].iloc[first] for m in match)
    assert diff_type / k > 0.9


def test_damaged_thins_only_non_mito_genes_in_a_random_fraction_of_cells(counts):
    out = cs.degrade(counts, "damaged_15", seed=0)
    before, after = counts.X.toarray(), out.X.toarray()
    mito = cs.mito_mask(counts)
    assert mito.sum() == 4
    assert (after[:, mito] == before[:, mito]).all()
    ratio = after[:, ~mito].sum(1) / np.maximum(before[:, ~mito].sum(1), 1)
    hit = ratio < 0.6
    assert hit.mean() == pytest.approx(0.15, abs=0.02)
    assert ratio[hit].mean() == pytest.approx(0.3, abs=0.05)


def test_mito_mask_reads_symbols_when_var_names_are_ensembl():
    var = pd.DataFrame({"feature_name": ["ACTB", "MT-CO1", "GAPDH"]},
                       index=["ENSG00000075624", "ENSG00000198804", "ENSG00000111640"])
    a = ad.AnnData(sp.csr_matrix(np.ones((3, 3), dtype=np.float32)), var=var)
    assert cs.mito_mask(a).tolist() == [False, True, False]


def test_mixed_is_doublets_then_ambient_then_depth(counts):
    out = cs.degrade(counts, "mixed", seed=0)
    assert out.n_obs == counts.n_obs + round(0.2 * counts.n_obs)
    expected = _total(counts) * 1.4 * 1.2 * 0.3   # 20% doublets each carry two cells' reads, ambient +20%, thin 0.3
    assert _total(out) == pytest.approx(expected, rel=0.08)


def test_degradation_is_deterministic_per_seed(counts):
    for cond in cs.CONDITIONS:
        a, b = cs.degrade(counts, cond, seed=5), cs.degrade(counts, cond, seed=5)
        assert (a.X != b.X).nnz == 0 and list(a.obs_names) == list(b.obs_names), cond
    assert (cs.degrade(counts, "depth_30", seed=1).X != cs.degrade(counts, "depth_30", seed=2).X).nnz > 0


def test_clean_copy_gives_the_same_plan_and_profile_as_the_original(tmp_path):
    from grnbench.autoprep import diagnose, plan_preprocessing
    from grnbench.simulate import PRESETS, write_scenario

    scn = PRESETS["clean"].with_seed(1)
    paths = write_scenario(tmp_path / "orig", scn)
    tfs = paths["tfs"].read_text(encoding="utf-8").split()
    copy = cs.write_degraded_copy(paths["counts"], "clean", 0, tmp_path / "copy")
    assert (tmp_path / "copy" / "README.md").exists() == (tmp_path / "orig" / "README.md").exists()
    plans = [plan_preprocessing(diagnose(p), tf_list=tfs) for p in (paths["counts"], copy)]
    assert plans[0].choices == plans[1].choices
    for key in ("technology", "organism"):
        assert plans[0].profile.get(key) == plans[1].profile.get(key)


def test_make_tasks_lists_cellxgene_datasets_times_conditions(tmp_path):
    for sub in ("a/human/x__1", "b/mouse/y__2"):
        d = tmp_path / "cellxgene" / sub
        d.mkdir(parents=True)
        (d / "truth.json").write_text(json.dumps({"organism": "human"}))
    tasks = cs.make_tasks(tmp_path, ["clean", "mixed"])
    assert [t["id"] for t in tasks] == ["cxg__a__human__x__1__clean", "cxg__a__human__x__1__mixed",
                                        "cxg__b__mouse__y__2__clean", "cxg__b__mouse__y__2__mixed"]
    assert all(t["seed"] == 0 and t["path"].endswith("data.h5ad") for t in tasks)
    assert len(cs.make_tasks(tmp_path, ["clean"], max_datasets=1)) == 1


def test_variants_exist_in_the_knowledge_base_and_differ_from_baseline_by_one_change():
    from grnbench.kb import DECISION_INDEX

    for decision, option in cs.VARIANTS:
        DECISION_INDEX[decision].option(option)
        assert cs.BASELINE[decision] != option
    cands, skipped = cs.variant_candidates()
    assert len(cands) + len(skipped) == len(cs.VARIANTS)
    for c in cands:
        assert sum(c["choices"][k] != v for k, v in cs.BASELINE.items()) == 1


def _row(cond, base, auto, noref, variant):
    def cand(name, source, v):
        return {"name": name, "source": source, "eval_mean": v, "choices": {"x": 1}}
    cands = [cand("recommended", "engine", auto), cand("baseline", "baseline", base),
             cand("base+hvg=tf_plus_500", "variant", variant)]
    return {"task": f"t__{cond}", "dataset": "t", "condition": cond, "seed": 0, "error": None,
            "rules": "recommended", "engine_auto": "recommended", "engine_auto_noref": noref,
            "candidates": cands, "n_cells_in": 10, "n_genes_in": 5, "reference_edges": 60}


def test_report_writes_candidates_csv_and_exploratory_summary(tmp_path):
    (tmp_path / "rows").mkdir()
    rows = [_row("clean", 1.0, 1.2, "baseline", 1.5), _row("mixed", 1.0, 0.8, "recommended", 0.5),
            {"task": "u__mixed", "dataset": "u", "condition": "mixed", "error": "Boom: x"}]
    for r in rows:
        (tmp_path / "rows" / f"{r['task']}.json").write_text(json.dumps(r))
    text = cs.report(tmp_path)
    df = pd.read_csv(tmp_path / "candidates.csv")
    assert len(df) == 6 and set(df["source"]) == {"engine", "baseline", "variant"}
    assert "EXPLORATORY" in text and (tmp_path / "SUMMARY.md").read_text(encoding="utf-8") == text
    assert "| clean | 1 | 0 | 0 | 1.000 | 1.200 | 1.000 | 1.200 | 1.200 | 1.500 | 0/1 |" in text
    assert "| mixed | 1 | 0 | 1 | 1.000 | 0.800 | 0.800 | 0.800 | 1.000 | 1.000 | 1/1 |" in text
    assert "base+hvg=tf_plus_500 | +0.500 [1] | -0.500 [1] |" in text


def test_rule_choice_uses_screening_scores_only_and_respects_pool():
    row = {"candidates": [
        {"name": "baseline", "source": "baseline", "eval_mean": 1.0, "sel_mean": 1.0, "alias_of": None},
        {"name": "recommended", "source": "engine", "eval_mean": 0.9, "sel_mean": 1.2, "alias_of": None},
        {"name": "base+x=y", "source": "variant", "eval_mean": 1.5, "sel_mean": 1.3, "alias_of": None},
        {"name": "dup", "source": "variant", "eval_mean": 0.9, "sel_mean": 9.0, "alias_of": "recommended"},
        {"name": "broken", "source": "variant", "eval_mean": float("nan"), "sel_mean": 5.0, "alias_of": None},
    ]}
    sel = lambda c: c["sel_mean"]  # noqa: E731
    assert cs.rule_choice(row, "engine", sel) == "recommended"
    assert cs.rule_choice(row, "all", sel) == "base+x=y"


def test_confirm_family_on_synthetic_rows():
    base = dict(cs.BASELINE)
    rows = []
    for i in range(6):
        cands = [
            {"name": "baseline", "source": "baseline", "choices": base, "eval_mean": 1.30 + 0.01 * i,
             "sel_mean": 1.30, "edge_stability": 0.5, "depth_leakage": 0.1, "batch_leakage": 0.0},
            {"name": "batch=center_per_batch", "source": "engine", "choices": {**base, "batch": "center_per_batch"},
             "eval_mean": 1.10 + 0.02 * (i % 3), "sel_mean": 1.10, "edge_stability": 0.9,
             "depth_leakage": 0.0, "batch_leakage": 0.0},
            {"name": "base+hvg=tf_plus_2000", "source": "variant", "choices": {**base, "hvg": "tf_plus_2000"},
             "eval_mean": 9.0, "sel_mean": 9.0, "edge_stability": 0.1, "depth_leakage": 0.5, "batch_leakage": 0.0},
        ]
        rows.append({"condition": "clean", "candidates": cands})
    recs, table = cs.confirm(rows, "clean")
    by = {(r["a"], r["b"]): r for r in recs}
    assert by[("baseline", "noref_shipped")]["supported"] == "baseline"
    # the gene-set-changing hvg variant is outside the preserving pool, so the rule picks baseline
    assert by[("ref_geneset_preserving", "baseline")]["diff"] == 0
    assert "| baseline | noref_shipped | 6 |" in table


def test_run_task_expands_portable_home_paths(tmp_path, monkeypatch):
    """real_bench.make_tasks writes "~/..." paths since #29; the job must still find the data."""
    from grnbench.simulate import PRESETS, write_scenario

    home = tmp_path / "home"
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("USERPROFILE", str(home))   # Path.home() on Windows
    paths = write_scenario(home / "data", PRESETS["small"].with_seed(0))
    task = {"id": "t", "dataset": "d", "condition": "clean", "seed": 0,
            "path": "~/data/" + paths["counts"].name, "truth": {"organism": "human"}}
    row = cs.run_study_task(task, None, tmp_path / "work")
    assert row["error"] is None, row.get("traceback")
    assert row["skipped"].startswith("reference_edges 0")


# -- shared genes with the baseline ------------------------------------------------------- #

class _Plan:
    def value(self, key):
        return "log_normalized"


def _records(tmp_path, gene_counts):
    """One processed h5ad per candidate over G0..G<n-1> (same cells), first 6 genes are the TFs."""
    rng = np.random.default_rng(0)
    base = rng.normal(size=(60, max(gene_counts.values())))
    base[:, 10:20] += base[:, [0]]                         # G0 drives G10..G19
    recs = {}
    for name, n in gene_counts.items():
        path = tmp_path / f"{name}.h5ad"
        a = ad.AnnData(base[:, :n].astype(np.float32), var=pd.DataFrame(index=[f"G{i}" for i in range(n)]))
        a.obs_names = [f"c{i}" for i in range(60)]
        a.write_h5ad(path)
        recs[name] = [{"processed": str(path)}]
    return recs


TFS = [f"G{i}" for i in range(6)]
PAIRS = {("G0", f"G{i}") for i in range(10, 20)} | {("G1", "G3"), ("G2", "G45"), ("G3", "G12")}


def _score(recs):
    from engine_bench import score_candidates

    genes = {n: cs.candidate_genes(r) for n, r in recs.items()}
    cache, out = {}, {}
    for name in sorted(recs, key=lambda n: n != "baseline"):
        out[name] = score_candidates({name: recs[name]}, {name: _Plan()}, TFS, PAIRS, pair_ref="baseline",
                                     pair_genes=genes, ref_cache=cache)[name]
    return out, genes


def test_pair_fields_baseline_equals_own_and_own_universe_is_unchanged(tmp_path):
    from engine_bench import score_candidates

    recs = _records(tmp_path, {"baseline": 40, "wide": 50, "narrow": 30})
    out, genes = _score(recs)
    for name, row in out.items():
        old = score_candidates({name: recs[name]}, {name: _Plan()}, TFS, PAIRS)[name]
        for k, v in old.items():                                           # exactly equal
            assert row[k] == v or (np.isnan(v) and np.isnan(row[k]))
        assert row["pair_n_genes"] == len(genes[name] & genes["baseline"])
    b = out["baseline"]
    assert b["pair_n_genes"] == 40
    assert b["pair_eval_mean"] == b["eval_mean"] == b["pair_base_eval_mean"]
    assert out["wide"]["pair_n_genes"] == 40               # superset of the baseline: G_C == baseline genes
    assert out["narrow"]["pair_n_genes"] == 30             # subset: the candidate's own genes
    for name in ("wide", "narrow"):
        assert np.isfinite(out[name]["pair_eval_mean"]) and np.isfinite(out[name]["pair_base_eval_mean"])
    # the baseline side is the baseline's own edges on its own genes, whatever the candidate
    assert out["wide"]["pair_base_eval_mean"] == b["eval_mean"]


def test_report_has_shared_genes_table_and_rule(tmp_path):
    (tmp_path / "rows").mkdir()
    row = _row("clean", 1.0, 1.2, "baseline", 1.5)
    for c, v, b in zip(row["candidates"], (1.1, 1.0, 1.3), (1.0, 1.0, 1.0)):
        c.update(pair_eval_mean=v, pair_base_eval_mean=b, sel_mean=v)
    (tmp_path / "rows" / "a.json").write_text(json.dumps(row))
    text = cs.report(tmp_path)
    assert "## Variants vs baseline on shared genes" in text
    assert "base+hvg=tf_plus_500 | +0.300 [1] |" in text.split("on shared genes")[1]
    assert "judged on shared genes | +0.300 (1/0) |" in text
