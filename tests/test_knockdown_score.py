"""benchmarks/knockdown_score.py (infer, score, combine) on a synthetic screen with known answers.

30 TFs; each responds on 12 genes it regulates plus 8 hub genes that respond to every knockdown. Controls are
Poisson counts where each TF drives its 12 targets, so Pearson recovers the regulatory edges for real. The expensive
models are stand-ins (spearman for GRNBoost2, a missing RegDiffusion) and the cisTarget files are replaced by a fake
motif-rank matrix; the motif combination, the gene_tf control, the engine preprocessing, the evaluator and the
permutation null all run their real code.
"""
import contextlib
import shutil
import hashlib
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
import knockdown_score as ks  # noqa: E402

from grnbench.eval import knockdown as kd  # noqa: E402

from grnbench.infer import base as ib  # noqa: E402
from grnbench.infer import cistarget  # noqa: E402
from grnbench.infer.correlation import spearman_signed_abs  # noqa: E402

N_TF, N_SPEC, N_HUB, N_NOISE, N_CELLS = 30, 12, 8, 82, 500
N_GENES = N_TF + N_TF * N_SPEC + N_HUB + N_NOISE
GENES = [f"G{i}" for i in range(N_GENES)]
TFS = GENES[:N_TF]
TARGETS = {t: GENES[N_TF + N_SPEC * i: N_TF + N_SPEC * (i + 1)] for i, t in enumerate(TFS)}
HUBS = GENES[N_TF + N_TF * N_SPEC: N_TF + N_TF * N_SPEC + N_HUB]
WORST = 5000


# ---------------------------------------------------------------- stand-ins
def fake_ranks(tfs, genes, mode):
    """TFs x genes motif ranks. Hubs get good ranks for every TF; in 'informative' mode a TF's own targets too."""
    genes = list(map(str, genes))
    col = {g: i for i, g in enumerate(genes)}
    M = np.empty((len(tfs), len(genes)), dtype=np.int32)
    for r, t in enumerate(map(str, tfs)):
        rng = np.random.default_rng(100 + TFS.index(t))          # a TF's row never depends on the other TFs
        M[r] = rng.integers(1500, WORST, size=len(genes))
        for g in HUBS:
            M[r, col[g]] = rng.integers(1, 800)
        if mode == "informative":
            for g in TARGETS[t]:
                M[r, col[g]] = rng.integers(1, 800)
    return M


def specs(mode):
    def oracle(expr, tf_names=None, **_):                         # the true regulatory edges, nothing else
        rng = np.random.default_rng(1)
        return pd.DataFrame([(t, g, 1.0 + 0.001 * rng.random()) for t in tf_names for g in TARGETS[t]],
                            columns=["source", "target", "weight"])

    def dense(extra):
        def run(expr, tf_names=None, **_):
            rng = np.random.default_rng(3)
            rows = [(t, g, 0.1 * rng.normal() + extra(g)) for t in tf_names for g in expr.columns if g != t]
            return pd.DataFrame(rows, columns=["source", "target", "weight"])
        return run

    def boom(expr, tf_names=None, **_):
        raise RuntimeError("boom: stand-in model crashed")

    mk = lambda name, runner, unavailable=None: ib.MethodSpec(name, "test stand-in", runner, unavailable=unavailable)
    return [mk("oracle", oracle), mk("hubs", dense(lambda g: 10.0 if g in HUBS else 0.0)),
            mk("random_net", dense(lambda g: 0.0)), mk("grnboost2", spearman_signed_abs), mk("boom", boom),
            mk("regdiffusion", None, "regdiffusion not installed (torch missing)")]


@contextlib.contextmanager
def patched(mode):
    with pytest.MonkeyPatch.context() as mp:
        _patch(mp, mode)
        yield


def _patch(mp, mode):
    for s in specs(mode):
        mp.setitem(ib.REGISTRY, s.name, s)
    mp.setattr(cistarget, "motif_ranks", lambda tfs, genes, species=None, root=None: pd.DataFrame(
        fake_ranks(tfs, genes, mode), index=list(map(str, tfs)), columns=list(map(str, genes))))
    mp.setattr(ks.ms, "motif_parts", lambda tfs_e, genes: (
        fake_ranks(tfs_e, genes, mode), np.ones(len(genes), bool), np.ones(len(tfs_e), int), WORST))


# ---------------------------------------------------------------- the screen
def write_truth(path: Path, tfs, perturbed_only=()):
    path.mkdir(parents=True, exist_ok=True)
    units = pd.DataFrame({"tf": list(tfs) + list(perturbed_only), "usable": [True] * len(tfs) + [False] * len(perturbed_only)})
    units["exclusion_reason"] = ""
    resp = pd.DataFrame([(t, g, 0.01, 1.0, 200, 0.9) for t in tfs for g in TARGETS[t] + HUBS],
                        columns=["tf", "gene", "q", "log2fc", "n_cells", "kd"])
    units.to_csv(path / "units.tsv", sep="\t", index=False)
    resp.to_csv(path / "responders.tsv", sep="\t", index=False)
    (path / "truth_half.txt").write_text(" ".join(f"t{i:04d}" for i in range(N_CELLS)), encoding="utf-8")
    sha = lambda f: hashlib.sha256((path / f).read_bytes()).hexdigest()
    (path / "truth.json").write_text(json.dumps({"h5ad": "synthetic.h5ad", "sha256": {
        "units.tsv": sha("units.tsv"), "responders.tsv": sha("responders.tsv")}}))


@pytest.fixture(scope="module")
def screen(tmp_path_factory):
    root = tmp_path_factory.mktemp("kdscore")
    rng = np.random.default_rng(11)
    z = rng.normal(size=(N_CELLS, N_TF))
    logmean = np.full((N_CELLS, N_GENES), 1.0) + rng.normal(0, 0.3, size=(N_CELLS, 1))   # library size
    for i, t in enumerate(TFS):
        logmean[:, i] += 0.5 * z[:, i]
        for g in TARGETS[t]:
            logmean[:, GENES.index(g)] += 0.8 * z[:, i]
    X = rng.poisson(np.exp(logmean)).astype(np.float32)
    obs = pd.DataFrame(index=[f"c{i:04d}" for i in range(N_CELLS)])
    a = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=GENES))
    a.uns["source_h5ad"] = "synthetic.h5ad"
    a.write_h5ad(root / "controls.h5ad")
    pd.DataFrame({"source": TFS + ["NOTHERE1", "NOTHERE2"], "target": "G470", "weight": 1}).to_csv(
        root / "collectri.tsv", sep="\t", index=False)
    pd.DataFrame([(t, g, 1) for t in TFS for g in TARGETS[t]], columns=["source", "target", "weight"]).to_csv(
        root / "reference.tsv", sep="\t", index=False)
    write_truth(root / "truth", TFS)
    write_truth(root / "other", TFS[:15], perturbed_only=TFS[15:18])      # perturbed there: 18; usable there: 15
    return root


def run(screen, tag, models, mode, perms=2000, seed=0, max_cells=3000, extra_score=()):
    root = screen / tag
    shutil.rmtree(root, ignore_errors=True)
    with patched(mode):
        argv = ["infer", "--experiment", ks.EXPERIMENT, "--controls", str(screen / "controls.h5ad"),
                "--tfs", str(screen / "collectri.tsv"), "--models", models, "--out", str(root / "infer"),
                "--max-cells", str(max_cells), "--attempts", str(root / "attempts.json")]
        if ks.main(argv):                       # a failed model: the first inference is incomplete, the job resubmits
            assert json.loads((root / "infer" / "infer.json").read_text())["status"] == "incomplete"
            shutil.rmtree(root / "infer")
            assert ks.main(argv) == 0           # same commit, second failure of the same arms: publishable
    assert ks.main(["score", "--infer", str(root / "infer"), "--truth", str(screen / "truth"), "--out", str(root / "score"),
                    "--perms", str(perms), "--seed", str(seed), *extra_score]) == 0
    return (json.loads((root / "infer" / "infer.json").read_text()), json.loads((root / "score" / "family.json").read_text()))


INFORMATIVE = ("pearson_signed,pearson_motif,pearson_motif_genetf,grnboost2,grnboost2_motif,grnboost2_motif_genetf,"
               "regdiffusion,regdiffusion_motif,regdiffusion_motif_genetf,oracle,random_net,boom")


@pytest.fixture(scope="module")
def informative(screen):
    return run(screen, "informative", INFORMATIVE, "informative",
               extra_score=["--collectri", str(screen / "reference.tsv"), "--essential-from", str(screen / "other"),
                            "--both-with", str(screen / "other")])


@pytest.fixture(scope="module")
def generic(screen):
    return run(screen, "generic", "pearson_signed,pearson_motif,pearson_motif_genetf,hubs,hubs_motif,hubs_motif_genetf",
               "generic")


def fam(result):
    return {r["hypothesis"]: r for r in result[1]["primary"]["family"]}


def supported(r):
    return r["verdict"] == "supported"


# ---------------------------------------------------------------- (a) a real network beats its controls
def test_true_regulatory_network_beats_its_propensity_and_the_permuted_null(informative):
    r = fam(informative)["H2 oracle > own propensity"]
    assert supported(r) and r["mean_contrast"] > 0.5
    assert r["p_t_greater_holm"] < 0.05 and r["p_perm_greater_holm"] < 0.05
    assert r["mean_contrast"] > r["null_mean"] + 5 * r["null_sd"]          # the pairing of TF and knockdown carries it
    means = informative[1]["primary"]["arm_means"]
    assert means["oracle"] > 5 * means["propensity:oracle"]


def test_pearson_run_for_real_recovers_the_edges_and_a_tf_specific_prior_adds_to_it(informative):
    f = fam(informative)
    assert supported(f["H2 pearson_signed > own propensity"])
    assert supported(f["H1a pearson_motif > pearson_signed"]) and supported(f["H1b pearson_motif > pearson_motif_genetf"])
    means = informative[1]["primary"]["arm_means"]
    assert means["pearson_signed"] > 3 and means["pearson_motif"] > means["pearson_motif_genetf"]


# ---------------------------------------------------------------- (b) hub-only networks are caught
def test_hub_only_network_is_caught(generic):
    f = fam(generic)
    h2 = f["H2 hubs > own propensity"]
    assert not supported(h2) and abs(h2["mean_contrast"]) < 0.5
    means = generic[1]["primary"]["arm_means"]
    assert means["hubs"] > 2 and means["hubs"] == pytest.approx(means["propensity:hubs"], rel=0.5)   # it scores well, but only on hubs


def test_generic_prior_lifts_pearson_but_the_lift_needs_no_tf_identity(generic):
    f = fam(generic)
    a = f["H1a pearson_motif > pearson_signed"]                            # hubs respond to every knockdown
    assert a["mean_contrast"] > 1 and a["p_t_greater_holm"] < 0.05         # the paired t sees the lift ...
    assert a["p_perm_greater_holm"] > 0.5 and not supported(a)             # ... the TF-label null reproduces it (more)
    assert a["null_mean"] > a["mean_contrast"]
    assert not supported(f["H1b pearson_motif > pearson_motif_genetf"])    # the gene_tf control gets the same lift
    assert not supported(f["H1b hubs_motif > hubs_motif_genetf"])
    assert supported(f["H2 pearson_signed > own propensity"])


# ---------------------------------------------------------------- (c) random passes nothing
def test_random_network_passes_nothing(informative):
    f = fam(informative)
    r = f["H2 random_net > own propensity"]
    assert r["verdict"] == "inconclusive"
    assert informative[1]["primary"]["arm_means"]["random_net"] < 1.6


# ---------------------------------------------------------------- (d) skipped models are recorded
def test_models_that_cannot_run_are_recorded_with_the_reason(informative):
    infer_meta, family = informative
    sk = infer_meta["skipped"]
    assert set(sk) == {"regdiffusion", "regdiffusion_motif", "regdiffusion_motif_genetf", "boom"}
    assert "torch" in sk["regdiffusion"] and "base model regdiffusion skipped" in sk["regdiffusion_motif"]
    assert "boom: stand-in model crashed" in sk["boom"]
    assert not set(sk) & set(infer_meta["arms"]) and family["skipped_arms"] == sk
    assert family["primary"]["n_directed_hypotheses"] == 2 * len(family["verdicts"])
    assert {h for h, v in family["verdicts"].items() if v == "missing"} == set(family["missing_hypotheses"])


def test_report_lists_skipped_arms_the_family_and_the_decision_rule(screen, informative):
    text = (screen / "informative" / "score" / "REPORT.md").read_text(encoding="utf-8")
    random_ratio = informative[1]["random_auprc_ratio"]
    assert random_ratio > 1
    assert f"random expectation: AUPRC ratio {random_ratio:.3f} (expected AP of a random order, mean over units)" in text
    assert "AUPRC ratio 1," not in text
    assert "## Arms that did not run" in text and "regdiffusion_motif_genetf: base model regdiffusion skipped" in text
    assert "H1b pearson_motif > pearson_motif_genetf" in text and "needs TF-specific motif identity" in text
    assert "If the other screen gives the same verdicts: oracle: its targets predict knockdown responses" in text
    assert "random_net: H2 inconclusive; reported as inconclusive, never as no difference." in text


# ---------------------------------------------------------------- inference record
def test_inference_record_has_cells_regulators_seconds_and_commit(informative):
    m = informative[0]
    assert m["status"] == "complete" and m["code_commit"] and m["n_genes"] == N_GENES
    assert m["n_regulators"] == N_TF and m["tfs_in_reference"] == N_TF + 2 and set(m["regulators"]) == set(TFS)
    assert m["cells_in"] == N_CELLS and m["cells_after_qc"] <= N_CELLS and m["cells_after_cap"] == m["cells_after_qc"]
    assert m["cap_applied"] is False and m["preprocessing"]["choices"]["hvg"] == "all_genes"
    assert all(v["seconds"] >= 0 and v["file"].startswith("edges/") for v in m["arms"].values())
    assert m["arms"]["pearson_motif"]["base"] == "pearson_signed" and m["arms"]["pearson_motif_genetf"]["kind"] == "genetf"


def test_cap_is_applied_after_qc_with_the_real_bench_rule(screen):
    with patched("informative"):
        meta = ks.infer(screen / "controls.h5ad", screen / "collectri.tsv", ["pearson_signed"], screen / "capped", max_cells=400,
                        experiment=ks.EXPERIMENT)
    assert meta["cap_applied"] and meta["cells_after_cap"] == 400 and meta["cells_after_qc"] > 400
    assert "after QC" in meta["cap_rule"] and meta["cap_seed"] == 0
    qc = ad.read_h5ad(json.loads((screen / "capped" / "prep" / "contexts.json").read_text())[0]["processed"])
    keep = np.sort(np.random.default_rng(0).choice(qc.n_obs, 400, replace=False))          # real_bench.score_models
    assert (screen / "capped" / "cells.txt").read_text().split() == list(qc.obs_names[keep])
    assert qc.n_obs == meta["cells_after_qc"] and qc.n_vars == meta["n_genes"] == N_GENES


def test_arm_names_are_validated():
    assert ks.parse_arm("pearson_signed") == ("pearson_signed", "base")
    assert ks.parse_arm("pearson_motif") == ("pearson_signed", "motif")
    assert ks.parse_arm("grnboost2_motif_genetf") == ("grnboost2", "genetf")
    for bad in ("grnboost2_genetf", "pearson_cistarget_genetf", "grnboost2_cistarget"):
        with pytest.raises(ValueError):
            ks.parse_arm(bad)


# ---------------------------------------------------------------- secondary outputs
def test_reference_row_and_strata(informative):
    p = informative[1]
    assert p["primary"]["arm_means"]["collectri"] > 10 * p["primary"]["arm_means"]["propensity:collectri"]
    assert not any("collectri" in h for h in p["verdicts"])                  # descriptive: outside the family
    assert p["strata_members"] == {"essential": 18, "non_essential": 12, "usable_in_both_screens": 15}
    assert p["strata"]["essential"]["n_units"] == 18 and p["strata"]["non_essential"]["n_units"] == 12
    assert set(p["detection"]) >= {"pearson_signed", "oracle"}


def test_scores_tsv_is_long_with_one_row_per_unit_and_arm(screen, informative):
    s = pd.read_csv(screen / "informative" / "score" / "scores.tsv", sep="\t")
    prim = s[(s.endpoint == "auprc_ratio") & (s.stratum == "all")]
    assert set(s.endpoint) == {"auprc_ratio", "top50", "top200"} and {"essential", "non_essential"} <= set(s.stratum)
    assert prim.groupby("arm").tf.nunique().eq(N_TF).all() and "response_propensity" in set(prim.arm)
    assert not prim.duplicated(["tf", "arm"]).any()


def test_truth_tables_must_match_their_hashes(screen, tmp_path):
    bad = tmp_path / "truth"
    write_truth(bad, TFS)
    (bad / "responders.tsv").write_text((bad / "responders.tsv").read_text() + "G0\tG1\t0.01\t1\t200\t0.9\n")
    with pytest.raises(ValueError, match="does not match the sha256"):
        ks.load_truth(bad)


# ---------------------------------------------------------------- (e) determinism
def test_same_seed_same_outputs_and_a_new_seed_moves_only_the_permutation_p(screen):
    models = "pearson_signed,pearson_motif,pearson_motif_genetf"
    a = run(screen, "det_a", models, "informative", perms=300, seed=5, max_cells=400)
    b = run(screen, "det_b", models, "informative", perms=300, seed=5, max_cells=400)
    c = run(screen, "det_c", models, "informative", perms=300, seed=6, max_cells=400)
    for arm in a[0]["arms"]:
        ea = pd.read_csv(screen / "det_a" / "infer" / a[0]["arms"][arm]["file"], sep="\t")
        eb = pd.read_csv(screen / "det_b" / "infer" / b[0]["arms"][arm]["file"], sep="\t")
        pd.testing.assert_frame_equal(ea, eb)
    for f in ("scores.tsv", "family.json"):
        assert (screen / "det_a" / "score" / f).read_bytes() == (screen / "det_b" / "score" / f).read_bytes()
    # REPORT.md also lists each arm's wall-clock seconds (last cell of the arms table, the only 5-cell table),
    # the one field that may differ between runs
    def no_time(tag):
        lines = (screen / tag / "score" / "REPORT.md").read_text(encoding="utf-8").splitlines()
        return [ln.rsplit("|", 2)[0] if ln.startswith("|") and ln.count("|") == 6 else ln for ln in lines]
    ra, rb = no_time("det_a"), no_time("det_b")
    assert ra == rb and any(ln.startswith("| pearson_signed |") for ln in ra)
    assert a[1]["inputs"]["edges_sha256"].keys() == b[1]["inputs"]["edges_sha256"].keys()
    assert (screen / "det_a" / "score" / "scores.tsv").read_bytes() == (screen / "det_c" / "score" / "scores.tsv").read_bytes()
    pa, pc = fam(a)["H2 pearson_signed > own propensity"], fam(c)["H2 pearson_signed > own propensity"]
    assert pa["mean_contrast"] == pc["mean_contrast"] and pa["p_t_greater"] == pc["p_t_greater"]
    assert (pa["null_mean"], pa["p_perm_greater"]) != (pc["null_mean"], pc["p_perm_greater"])


# ---------------------------------------------------------------- decision rule and combine
def test_decision_rule_per_family_and_per_model():
    v = {"H1a pearson_motif > pearson_signed": "supported (replicated)",
         "H1b pearson_motif > pearson_motif_genetf": "inconclusive",
         "H1a grnboost2_motif > grnboost2": "supported (replicated)",
         "H1b grnboost2_motif > grnboost2_motif_genetf": "supported (replicated)",
         "H2 pearson_signed > own propensity": "supported in design, inconclusive in replication",
         "H2 grnboost2 > own propensity": "inconclusive"}
    text = "\n".join(ks.decision_lines(v, combined=True))
    assert "pearson_motif: the motif lift on knockdown truth is consistent with a TF-generic gene propensity" in text
    assert "grnboost2_motif: the motif prior raises recovery" in text and "needs TF-specific motif identity" in text
    assert "pearson_signed: H2 not replicated, no claim" in text and "grnboost2: H2 inconclusive" in text
    none = {"H2 a > own propensity": "inconclusive", "H2 b > own propensity": "inconclusive"}
    assert "no evidence, in either screen" in "\n".join(ks.decision_lines(none, combined=True))
    # one screen inconclusive and the other supported is "not replicated", never the all-inconclusive H2 line
    one = {"H2 a > own propensity": "inconclusive in design, supported in replication",
           "H2 b > own propensity": "supported in design, missing in replication"}
    text = chr(10).join(ks.decision_lines(one, combined=True))
    assert "a: H2 not replicated, no claim" in text and "b: H2 not replicated, no claim" in text
    assert "inconclusive for every model" not in text
    refuted = {"H1a x_motif > x": "refuted (replicated)"}
    assert "lowers recovery of causal responders" in "\n".join(ks.decision_lines(refuted, combined=True))


def copy_score(src: Path, dst: Path, **changes) -> Path:
    """A copy of a score folder whose family.json carries different top-level fields."""
    shutil.copytree(src, dst)
    fj = json.loads((dst / "family.json").read_text(encoding="utf-8"))
    (dst / "family.json").write_text(json.dumps({**fj, **changes}), encoding="utf-8")
    return dst


def test_combine_refuses_two_runs_on_the_same_screen(screen, informative, tmp_path):
    with pytest.raises(SystemExit, match="same screen"):
        ks.main(["combine", "--design", str(screen / "informative" / "score"),
                 "--replication", str(screen / "informative" / "score"), "--out", str(tmp_path / "c")])


def test_score_refuses_overlapping_cells_a_swapped_screen_and_a_missing_truth_json(screen, informative, tmp_path):
    def score(truth, name="s"):
        return ks.main(["score", "--infer", str(screen / "informative" / "infer"), "--truth", str(truth),
                        "--out", str(tmp_path / name), "--perms", "50"])
    bad = tmp_path / "overlap"
    shutil.copytree(screen / "truth", bad)
    cells = (screen / "informative" / "infer" / "cells.txt").read_text().split()
    (bad / "truth_half.txt").write_text(" ".join(["t9999", cells[0]]), encoding="utf-8")
    with pytest.raises(SystemExit, match="in both the inference cells and the truth half"):
        score(bad)
    gone = tmp_path / "nohalf"
    shutil.copytree(screen / "truth", gone)
    (gone / "truth_half.txt").unlink()
    with pytest.raises(SystemExit, match="truth_half.txt"):
        score(gone)
    swapped = tmp_path / "swapped"
    shutil.copytree(screen / "truth", swapped)
    meta = json.loads((swapped / "truth.json").read_text())
    (swapped / "truth.json").write_text(json.dumps({**meta, "h5ad": "other.h5ad"}))
    with pytest.raises(SystemExit, match="other.h5ad"):
        score(swapped)
    (swapped / "truth.json").unlink()
    with pytest.raises(SystemExit, match="truth.json is missing"):
        score(swapped)


def test_registered_hypotheses_that_could_not_be_computed_are_listed_with_the_reason(screen, informative):
    p = informative[1]
    miss = p["missing_hypotheses"]
    assert set(miss) == {h for h, _ in kd.default_family(ks.DEFAULT_MODELS.split(","))[0] if "regdiffusion" in h} \
        | {"H2 boom > own propensity"}
    assert "skipped by infer" in miss["H2 regdiffusion > own propensity"]
    assert p["source_h5ad"] == "synthetic.h5ad" and all(p["verdicts"][h] == "missing" for h in miss)
    text = (screen / "informative" / "score" / "REPORT.md").read_text(encoding="utf-8")
    assert "## Registered hypotheses that could not be computed" in text and "H2 regdiffusion > own propensity" in text
    # a model that was never requested is reported as not requested
    short = run(screen, "short", "pearson_signed", "informative", perms=50)[1]["missing_hypotheses"]
    assert short["H2 grnboost2 > own propensity"] == "grnboost2: not requested"


def test_combine_needs_the_same_verdict_in_both_screens(screen, informative, generic, tmp_path):
    out = tmp_path / "combined"
    other = copy_score(screen / "informative" / "score", tmp_path / "other_score", source_h5ad="other.h5ad")
    assert ks.main(["combine", "--design", str(screen / "informative" / "score"),
                    "--replication", str(other), "--out", str(out)]) == 0
    final = json.loads((out / "combined.json").read_text())
    assert final["H2 oracle > own propensity"] == "supported (replicated)"
    assert final["H2 random_net > own propensity"] == "inconclusive"
    assert "oracle: its targets predict knockdown responses beyond gene propensity" in (out / "COMBINED.md").read_text()
    # design vs replication with different verdicts for the same hypothesis: no claim
    out2 = tmp_path / "mixed"
    assert ks.main(["combine", "--design", str(screen / "informative" / "score"),
                    "--replication", str(copy_score(screen / "generic" / "score", tmp_path / "g", source_h5ad="other.h5ad")),
                    "--out", str(out2)]) == 0
    mixed = json.loads((out2 / "combined.json").read_text())
    assert mixed["H1b pearson_motif > pearson_motif_genetf"].startswith("supported in design, ")


# ---------------------------------------------------------------- --total-shift (check P9 sensitivity)
def test_total_shift_stratum_is_added_and_never_changes_a_primary_verdict(screen, tmp_path):
    models = "pearson_signed,pearson_motif,pearson_motif_genetf"
    plain = run(screen, "shift_a", models, "informative", perms=300, max_cells=400)
    shift = tmp_path / "shift.tsv"
    pd.DataFrame({"tf": TFS, "n_cells": 50, "rel_core_adjusted_total": 1.0,
                  "within_0_9_1_1": [i < 20 for i in range(N_TF)]}).to_csv(shift, sep="\t", index=False)
    withflag = run(screen, "shift_b", models, "informative", perms=300, max_cells=400, extra_score=["--total-shift", str(shift)])
    assert withflag[1]["strata_members"] == {"total_within_0.9_1.1": 20}
    assert withflag[1]["strata"]["total_within_0.9_1.1"]["n_units"] == 20
    assert withflag[1]["verdicts"] == plain[1]["verdicts"] and withflag[1]["primary"] == plain[1]["primary"]
    text = (screen / "shift_b" / "score" / "REPORT.md").read_text(encoding="utf-8")
    assert "sensitivity analysis declared in check P9" in text and "never changes a primary verdict" in text


# ---------------------------------------------------------------- exclusion (c): the retry, the family, the resume guard
def infer_argv(screen, out, attempts, models="pearson_signed,pearson_motif,regdiffusion"):
    return ["infer", "--experiment", ks.EXPERIMENT, "--controls", str(screen / "controls.h5ad"),
            "--tfs", str(screen / "collectri.tsv"), "--models", models, "--out", str(out), "--attempts", str(attempts)]


def test_a_failed_model_fails_the_first_inference_and_publishes_only_after_the_second(screen, tmp_path, monkeypatch):
    attempts = tmp_path / "run" / "rpe1" / "infer_attempts.json"
    with patched("informative"):
        assert ks.main(infer_argv(screen, tmp_path / "a1", attempts)) == 3
        m1 = json.loads((tmp_path / "a1" / "infer.json").read_text())
        assert m1["status"] == "incomplete" and m1["attempt"] == {"regdiffusion": 1}
        with pytest.raises(ValueError, match="not complete"):
            ks.load_arms(tmp_path / "a1")                             # score refuses it
        assert ks.main(infer_argv(screen, tmp_path / "a2", attempts)) == 0
        m2 = json.loads((tmp_path / "a2" / "infer.json").read_text())
        assert m2["status"] == "complete" and m2["attempt"] == {"regdiffusion": 2} and set(m2["skipped"]) == {"regdiffusion"}
        monkeypatch.setattr(ks.ms, "code_commit", lambda: "other")    # another commit counts from one again
        assert ks.main(infer_argv(screen, tmp_path / "a3", attempts)) == 3
    assert json.loads((tmp_path / "a3" / "infer.json").read_text())["attempt"] == {"regdiffusion": 1}


def test_an_inference_without_failures_is_complete_on_the_first_run(screen, tmp_path):
    with patched("informative"):
        assert ks.main(infer_argv(screen, tmp_path / "ok", tmp_path / "att.json", "pearson_signed,pearson_motif")) == 0
    assert json.loads((tmp_path / "ok" / "infer.json").read_text())["status"] == "complete"
    assert json.loads((tmp_path / "att.json").read_text()).get("failures", {}) == {}   # only the log of arms run


def test_a_failed_arm_does_not_shrink_the_holm_family(informative):
    from grnbench.eval.comparison import holm_adjust
    registered = [h for h, _ in kd.default_family(ks.DEFAULT_MODELS.split(","))[0]]
    assert len(registered) == 12
    prim = informative[1]["primary"]                                 # RegDiffusion failed in this run
    fam = prim["family"]
    assert set(registered) <= {r["hypothesis"] for r in fam} and prim["n_directed_hypotheses"] == 2 * len(fam) >= 24
    gone = [r for r in fam if "regdiffusion" in r["hypothesis"]]
    assert len(gone) == 4 and all(r["verdict"] == "missing" and r["p_perm_less"] == 1.0 for r in gone)
    for kind in ("t", "perm"):
        raw = [r[f"p_{kind}_{d}"] for r in fam for d in ("greater", "less")]
        assert [r[f"p_{kind}_{d}_holm"] for r in fam for d in ("greater", "less")] == pytest.approx(holm_adjust(raw))


def identity(**kw):
    return {"experiment": ks.EXPERIMENT, "code_commit": "c1", "snapshot_sha256": "s1", "models": ["a", "b"], "perms": 10000,
            "norm": "core", **kw}


def test_run_identity_writes_once_and_refuses_any_difference(tmp_path):
    run = tmp_path / "run"
    ks.check_run(run, identity())
    assert json.loads((run / "run.json").read_text())["code_commit"] == "c1"
    ks.check_run(run, identity())                                   # same job again: resumes
    for change in ({"code_commit": "c2"}, {"snapshot_sha256": "s2"}, {"models": ["a"]}, {"perms": 200}, {"norm": "cp10k"}):
        with pytest.raises(SystemExit, match=next(iter(change))):
            ks.check_run(run, identity(**change))


def test_a_folder_with_steps_and_no_run_json_is_not_resumed(tmp_path):
    (tmp_path / "run" / "rpe1" / "truth").mkdir(parents=True)
    with pytest.raises(SystemExit, match="no run.json"):
        ks.check_run(tmp_path / "run", identity())
    assert not (tmp_path / "run" / "run.json").exists()


def test_score_folders_of_another_run_are_refused_before_combine(screen, informative, tmp_path):
    run = tmp_path / "run"
    for s in ks.SCREENS:
        shutil.copytree(screen / "truth", run / s / "truth")
        shutil.copytree(screen / "informative" / "infer", run / s / "infer")
        shutil.copytree(screen / "informative" / "score", run / s / "score")
    ident = identity(code_commit=informative[0]["code_commit"], models=informative[0]["models_requested"], perms=2000)
    with pytest.raises(SystemExit, match="no run.json"):
        ks.check_run(run, ident)                                    # steps nobody vouches for
    (run / "run.json").write_text(json.dumps(ident))
    ks.check_run(run, ident)                                        # both score folders verified
    fam = run / "k562" / "score" / "family.json"
    d = json.loads(fam.read_text())
    fam.write_text(json.dumps({**d, "perms": 50}))
    with pytest.raises(SystemExit, match="k562/score"):
        ks.check_run(run, ident)
    fam.write_text(json.dumps({**d, "inputs": {**d["inputs"], "code_commit": "old"}}))
    with pytest.raises(SystemExit, match="k562/score"):
        ks.check_run(run, ident)
    fam.write_text(json.dumps(d))
    (run / "rpe1" / "truth" / "truth.json").write_text(json.dumps({"h5ad": "x", "sha256": {"units.tsv": "0"}}))
    with pytest.raises(SystemExit, match="rpe1/score"):
        ks.check_run(run, ident)


def test_run_identity_command_hashes_the_snapshot(tmp_path):
    snap = tmp_path / "snap.tsv"
    snap.write_text("a\n")
    args = ["run-identity", "--run", str(tmp_path / "run"), "--snapshot", str(snap), "--perms", "10", "--norm", "core"]
    assert ks.main(args) == 0 and ks.main(args) == 0
    snap.write_text("b\n")
    with pytest.raises(SystemExit, match="snapshot_sha256"):
        ks.main(args)


# ---------------------------------------------------------------- exclusion (c): a job killed mid-arm counts too
def die_in_regdiffusion(monkeypatch):
    """A process killed inside an arm (Slurm timeout): nothing after run_method runs, as with a BaseException."""
    real = ks.run_method

    def run_method(base, *a, **kw):
        if base == "regdiffusion":
            raise KeyboardInterrupt("simulated Slurm kill")
        return real(base, *a, **kw)
    monkeypatch.setattr(ks, "run_method", run_method)


def killed_run(screen, out, attempts, monkeypatch):
    with monkeypatch.context() as mp:
        die_in_regdiffusion(mp)
        with pytest.raises(KeyboardInterrupt):
            ks.main(infer_argv(screen, out, attempts))


def test_an_arm_killed_mid_run_counts_as_one_failed_attempt_on_the_next_run(screen, tmp_path, monkeypatch):
    attempts = tmp_path / "run" / "rpe1" / "infer_attempts.json"
    with patched("informative"):
        killed_run(screen, tmp_path / "k1", attempts, monkeypatch)
        left = json.loads(attempts.read_text())
        assert left["arms"]["regdiffusion"]["status"] == "started" and left["arms"]["regdiffusion"]["started_at"]
        assert left["arms"]["pearson_signed"]["status"] == "done"
        assert ks.main(infer_argv(screen, tmp_path / "k2", attempts)) == 0          # kill (1) + caught failure (2): published
    m = json.loads((tmp_path / "k2" / "infer.json").read_text())
    assert m["status"] == "complete" and m["interrupted"] == ["regdiffusion"]
    assert m["attempt"] == {"regdiffusion": 2}


def test_two_killed_runs_on_one_commit_reach_the_second_failure_disposition(screen, tmp_path, monkeypatch):
    attempts = tmp_path / "run" / "rpe1" / "infer_attempts.json"
    with patched("informative"):
        killed_run(screen, tmp_path / "k1", attempts, monkeypatch)
        killed_run(screen, tmp_path / "k2", attempts, monkeypatch)                  # counts the first kill, dies again
        assert ks.main(infer_argv(screen, tmp_path / "k3", attempts)) == 0          # counts the second: published
    m = json.loads((tmp_path / "k3" / "infer.json").read_text())
    assert m["status"] == "complete" and m["attempt"] == {"regdiffusion": 2} and m["interrupted"] == ["regdiffusion"]
    assert "did not finish" in m["skipped"]["regdiffusion"] and "regdiffusion" not in m["arms"]
    assert json.loads(attempts.read_text())["failures"] == {"regdiffusion": 2}      # not rerun a third time


def test_a_killed_run_is_not_counted_against_another_commit(screen, tmp_path, monkeypatch):
    attempts = tmp_path / "run" / "rpe1" / "infer_attempts.json"
    with patched("informative"):
        killed_run(screen, tmp_path / "k1", attempts, monkeypatch)
        monkeypatch.setattr(ks.ms, "code_commit", lambda: "other")
        assert ks.main(infer_argv(screen, tmp_path / "k2", attempts)) == 3
    m = json.loads((tmp_path / "k2" / "infer.json").read_text())
    assert m["attempt"] == {"regdiffusion": 1} and "interrupted" not in m
