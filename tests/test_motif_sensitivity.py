"""benchmarks/motif_sensitivity.py (EXP-20261003-GB-1, stage 1): synthetic checks only.

The mock run folders are written by ``real_bench.score_models`` itself (adapted from
``audit/2026-09-27/proposal/mock_run.py``), so they have the exact ``--keep-artifacts`` layout
and the "committed" rows the calibration gate must reproduce were computed by the registered
code, not by the script under test. A fake cisTarget folder (two hg38 databases + an hgnc table)
and a fake ``collectri_human.tsv`` stand in for the registered files.

- planted mock: half of every TF's targets carry its motif at the top of the motif rankings
  (target-level signal);
- null mock: motif rankings are random (no target information), the number of motifs per TF
  tracks its reference degree (TF-level effect), and Ensembl-ID genes plus symbol genes absent
  from the databases are tied at the worst rank (membership effect). Target sets are drawn per
  dataset, so a single random database cannot line up with them by chance across datasets.
- gene mock (from the agent review of the draft record, 2026-10-05: an agent with a fresh context
  started from the session that wrote the script, not an independent review in the sense of
  AGENTS.md section 5): every motif of every TF shares one ranking ordered by a per-gene
  propensity, and the same propensity decides how many TFs regulate a gene in the reference (hub
  targets). The prior knows which genes are regulated, not by which TF: no TF-specific target
  information. Before T3 it was called target-level.
- multi mock: the planted mock with two contexts per arm (the committed rows have 1-4), so the
  dataset mean over contexts and over each permutation seed is exercised; a one-dataset mouse mock
  (mm10 databases, an mgi table, collectri_mouse.tsv) runs through run_unit; tiny contexts give the
  degenerate shapes (no motif, one TF with a motif, too few positives, no gene in the databases).
"""
from __future__ import annotations

import json
import os
import shutil
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

pf = pytest.importorskip("pyarrow.feather")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "benchmarks"))
import motif_sensitivity as ms  # noqa: E402
import real_bench  # noqa: E402
from grnbench.infer import cistarget as ct  # noqa: E402

N_PERM = 8             # the registered run uses 100; the mocks only need the mechanics
CAP = 3000             # real_bench.ARTIFACT_TOP_EDGES in the mocks, so kept files are truncated
EXP = "MOCK-1"
N_TASKS = 5
N_TF, N_SYM, N_ENS, N_BG, N_OUT_DB, N_NO_MOTIF = 20, 220, 30, 1000, 25, 4


class _Plan:
    def __init__(self, strategy):
        self.strategy = strategy

    def value(self, key, default=None):
        return self.strategy if key == "strategy" else default


def build_mock(root: Path, *, plant: bool, n_tasks: int = N_TASKS, seed: int = 7, eff: float = 0.12,
               gene: bool = False, models=("pearson_signed", "pearson_motif", "motif_prior"), n_ctx: int = 1,
               organism: str = "human") -> dict:
    """Run folder (rows + artifacts), committed folder (rows + REPORT_REAL.md), cisTarget, refs.

    n_ctx contexts per arm (the committed rows have 1-4). organism "mouse" uses title-case symbols,
    mm10 databases, an mgi motif table and collectri_mouse.tsv (12 of the 40 registered datasets are
    mouse); the defaults give the same mock as before, name for name and number for number."""
    rng = np.random.default_rng(seed)
    mouse = organism == "mouse"
    n_tf, n_sym, n_ens = N_TF, N_SYM, N_ENS
    tfs = [f"Tf{i}" if mouse else f"TF{i}" for i in range(n_tf)]
    sym = [f"Gn{i}" if mouse else f"G{i}" for i in range(n_sym)]
    ens = [f"ENSMUSG{i:011d}" if mouse else f"ENSG{i:011d}" for i in range(n_ens)]
    genes = tfs + sym + ens
    K = np.r_[np.zeros(N_NO_MOTIF, int), rng.integers(1, 9, n_tf - N_NO_MOTIF)]   # motif-less TFs
    not_in_db = set(rng.choice(sym, N_OUT_DB, replace=False))           # symbols absent from the DBs
    db_genes = [g for g in tfs + sym if g not in not_in_db] + [f"Bg{i}" if mouse else f"BG{i}" for i in range(N_BG)]
    if gene:                       # one shared ranking by gene propensity; hub targets in the reference
        prop = dict(zip(sym, rng.normal(size=n_sym)))
        w = np.exp(1.5 * np.array([prop[g] for g in sym]))
        truths = [{t: list(rng.choice(sym, 15, replace=False, p=w / w.sum())) for t in tfs}
                  for _ in range(n_tasks)]
        score = np.array([(prop[g] if g in prop else rng.normal()) + 0.5 * rng.normal() for g in db_genes])
        shared = np.empty(len(db_genes), np.int32)
        shared[np.argsort(-score, kind="stable")] = np.arange(len(db_genes), dtype=np.int32)
    elif plant:
        fixed = {t: list(rng.choice(sym, 20, replace=False)) for t in tfs}
        truths = [fixed] * n_tasks
    else:                                                               # degree grows with K
        truths = [{t: list(rng.choice(sym, int(6 + 3 * k), replace=False)) for t, k in zip(tfs, K)}
                  for _ in range(n_tasks)]

    cis = root / "cistarget"
    cis.mkdir(parents=True)
    pos = {g: k for k, g in enumerate(db_genes)}
    ann = []
    for db in ("mm10_a_fake", "mm10_b_fake") if mouse else ("hg38_a_fake", "hg38_b_fake"):
        ids, rows = [], []
        for i, t in enumerate(tfs):
            for j in range(K[i]):
                r = shared.copy() if gene else rng.permutation(len(db_genes)).astype(np.int32)
                if plant:                                               # half the targets carry the motif
                    for g in truths[0][t][:10]:
                        if g in pos:
                            a = np.flatnonzero(r == rng.integers(0, 150))[0]
                            r[a], r[pos[g]] = r[pos[g]], r[a]
                ids.append(f"{db}_m{i}_{j}")
                rows.append(r)
                ann.append((f"{db}_m{i}_{j}", t))
        for j in range(30):
            ids.append(f"{db}_bg{j}")
            rows.append(rng.permutation(len(db_genes)).astype(np.int32))
        tab = pd.DataFrame(np.vstack(rows), columns=db_genes)
        tab.insert(0, "motifs", ids)
        pf.write_feather(tab, str(cis / f"{db}.genes_vs_motifs.rankings.feather"))
    pd.DataFrame(ann, columns=["#motif_id", "gene_name"]).to_csv(
        cis / f"motifs-fake.{'mgi' if mouse else 'hgnc'}-m0.001-o0.0.tbl", sep="\t", index=False)

    refs = root / "references"
    refs.mkdir()
    ref_rows = sorted({(t, g) for tr in truths for t in tfs for g in tr[t]}) + \
        [(t, f"OTHER{k}") for t in tfs for k in range(3)]
    pd.DataFrame(ref_rows, columns=["source", "target"]).to_csv(refs / f"collectri_{organism}.tsv", sep="\t",
                                                                index=False)

    def context(true, n_cells, cseed, counts):
        r = np.random.default_rng(cseed)
        act = r.normal(size=(n_cells, n_tf))
        X = {t: act[:, i] for i, t in enumerate(tfs)}
        for g in sym + ens:
            X[g] = r.normal(size=n_cells) * 0.8
        for i, t in enumerate(tfs):
            for g in true[t]:
                X[g] = X[g] + eff * act[:, i]
        for g in ens:
            X[g] = X[g] + eff * act[:, r.integers(n_tf)]
        lam = np.exp(0.5 * pd.DataFrame(X)[genes].to_numpy() + 1.0)
        mat = r.poisson(lam).astype(np.float32)
        mat = mat if counts else np.log1p(mat).astype(np.float32)
        return ad.AnnData(mat, var=pd.DataFrame(index=genes),
                          obs=pd.DataFrame(index=[f"c{k}" for k in range(n_cells)]))

    run, work = root / "run", root / "work"
    (run / "rows").mkdir(parents=True)
    work.mkdir()
    committed = root / "results" / EXP
    (committed / "rows").mkdir(parents=True)
    old_cap, old_env = real_bench.ARTIFACT_TOP_EDGES, os.environ.get(ct.ENV_DIR)
    real_bench.ARTIFACT_TOP_EDGES = CAP
    os.environ[ct.ENV_DIR] = str(cis)
    try:
        for ti in range(n_tasks):
            task = f"cxg__mock__{'mus-musculus' if mouse else 'homo-sapiens'}__t{ti}"
            pairs = {(t, g) for t in tfs for g in truths[ti][t]}
            sel, ev = real_bench_split(pairs)
            recs = {"baseline": [], "candX": [{"processed": None, "context": "skipped"}]}  # auto: from ctx1
            for cand, counts in (("baseline", False), ("candX", True)):
                for k in range(n_ctx):
                    p = work / (f"{task}_{cand}" + (f"_{k}" if k else "") + ".h5ad")
                    context(truths[ti], 220, 100 * ti + (cand == "candX") + 7 * k, counts).write_h5ad(p)
                    recs[cand].append({"processed": str(p), "context": "all" + (str(k) if k else "")})
            ta = run / "artifacts" / task
            ta.mkdir(parents=True)
            for half, pr in (("selection", sel), ("evaluation", ev)):
                pd.DataFrame(sorted(pr), columns=["source", "target"]).to_csv(
                    ta / f"reference_{half}.tsv", sep="\t", index=False)
            (ta / "decisions.json").write_text(json.dumps({
                "candidates": {"baseline": "baseline", "auto": "candX"},
                "chosen_plans": {"baseline": {"candidate": "baseline", "choices": {"normalization": "cpm_log1p"}},
                                 "auto": {"candidate": "candX", "choices": {"normalization": "raw_counts"}}}}))
            row = {"task": task, "truth": json.dumps({"organism": organism}), "reference_edges": len(pairs)}
            row.update(real_bench.score_models(
                recs, {"baseline": _Plan("cpm_log1p"), "candX": _Plan("raw_counts")},
                {"baseline": "baseline", "auto": "candX"}, list(models), sorted({s for s, _ in pairs}), ev,
                artifacts_dir=ta))
            text = json.dumps(row, default=float)
            (run / "rows" / f"{task}.json").write_text(text, encoding="utf-8")
            (committed / "rows" / f"{task}.json").write_text(text, encoding="utf-8")
    finally:
        real_bench.ARTIFACT_TOP_EDGES = old_cap
        if old_env is None:
            os.environ.pop(ct.ENV_DIR, None)
        else:
            os.environ[ct.ENV_DIR] = old_env
    df = pd.DataFrame([json.loads(p.read_text()) for p in sorted((committed / "rows").glob("*.json"))])
    (committed / "REPORT_REAL.md").write_text("\n".join(real_bench.models_section(df)) + "\n", encoding="utf-8")
    return {"root": root, "run": run, "cis": cis, "refs": refs, "results": root / "results"}


def real_bench_split(pairs):
    from grnbench.screen import split_reference
    return split_reference(pairs, 0)


def run_pipeline(mock: dict, out: Path, n_perm: int = N_PERM) -> dict:
    """plan -> run -> report through main(), as on the cluster."""
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(real_bench, "ARTIFACT_TOP_EDGES", CAP)
        mp.setenv(ct.ENV_DIR, str(mock["cis"]))
        tasks = out / "sens_tasks.json"
        assert ms.main(["plan", "--run", f"{EXP}={mock['run']}", "--out", str(tasks),
                        "--results-root", str(mock["results"])]) == 0
        assert ms.main(["run", "--tasks", str(tasks), "--out", str(out), "--refs", str(mock["refs"]),
                        "--n-perm", str(n_perm), "--workroot", str(mock["root"] / "work")]) == 0
        assert ms.main(["report", "--out", str(out)]) == 0
    return json.loads((out / "SENSITIVITY.json").read_text())


@pytest.fixture(scope="module")
def planted(tmp_path_factory):
    root = tmp_path_factory.mktemp("planted")
    mock = build_mock(root, plant=True)
    mock["summary"] = run_pipeline(mock, root / "sens", n_perm=20)  # rank p can reach 1/21 < 0.05
    mock["out"] = root / "sens"
    return mock


@pytest.fixture(scope="module")
def gene(tmp_path_factory):
    root = tmp_path_factory.mktemp("gene")
    mock = build_mock(root, plant=False, gene=True)
    mock["summary"] = run_pipeline(mock, root / "sens", n_perm=20)
    mock["out"] = root / "sens"
    return mock


@pytest.fixture(scope="module")
def null(tmp_path_factory):
    root = tmp_path_factory.mktemp("null")
    mock = build_mock(root, plant=False, seed=11)
    mock["summary"] = run_pipeline(mock, root / "sens")
    mock["out"] = root / "sens"
    return mock


@pytest.fixture(scope="module")
def multi(tmp_path_factory):
    """Two contexts per arm (the committed rows have 1-4): the mean over contexts is exercised."""
    root = tmp_path_factory.mktemp("multi")
    mock = build_mock(root, plant=True, n_tasks=3, n_ctx=2)
    mock["summary"] = run_pipeline(mock, root / "sens", n_perm=4)
    mock["out"] = root / "sens"
    return mock


def _arm(mock, arm="auto"):
    return mock["summary"]["sets"][EXP]["arms"][arm]


def _test(arm_summary, name):
    return next(t for t in arm_summary["tests"] if t["test"] == name)


# (a) planted target-level signal is called target-level ------------------------------------- #

def test_planted_motif_mock_is_called_target_level(planted):
    s = _arm(planted)
    t1, t2, t3 = _test(s, "T1"), _test(s, "T2"), _test(s, "T3")
    assert t1["holm_family"] == 6 and t1["n"] == N_TASKS
    assert t1["difference"] > 0 and t1["p_a_gt_b_holm"] < 0.05          # T1 supported
    assert t2["difference"] > 0 and t2["p_a_gt_b_holm"] < 0.05          # T2 supported
    assert t3["difference"] > 0 and t3["p_a_gt_b_holm"] < 0.05          # T3 supported: TF-specific
    assert s["L"] > 0 and s["delta_target"] >= 0.5 * s["L"] and s["L_gene"] < 0.5 * s["L"]
    assert s["verdict"] == "target-level" and s["shortfalls"] == []
    assert s["p_perm_T1"] < 0.05 and s["t1_robust_to_database"] is True   # robustness condition met
    assert s["downgraded_by_t1_robustness"] is False and s["secondary_disagreements"] == []
    assert s["secondary_not_compared"] == []                            # every pair compared on a healthy run
    text = (planted["out"] / "SENSITIVITY.md").read_text()
    assert "verdict (auto): supports target-level evidence" in text
    assert "Secondary checks (sign agreement; they never change the verdict): all agree." in text


# (a') a TF-generic gene propensity is not --------------------------------------------------- #

def test_gene_level_mock_is_not_called_target_level(gene):
    s = _arm(gene)
    t1, t2, t3 = _test(s, "T1"), _test(s, "T2"), _test(s, "T3")
    assert s["L"] > 0.5                                       # a large pooled lift ...
    assert t1["p_a_gt_b_holm"] < 0.05 and t2["p_a_gt_b_holm"] < 0.05   # ... that passes T1 and T2
    assert t3["difference"] == 0 and not t3["p_a_gt_b_holm"] < 0.05   # gene_tf IS the prior here
    assert s["L_gene"] == s["L"]                              # the gene propensity reaches all of it
    assert s["verdict"] == "mixed"
    assert [x[:2] for x in s["shortfalls"]][:1] == ["T3"] and any(x.startswith("L_gene") for x in s["shortfalls"])
    text = (gene["out"] / "SENSITIVITY.md").read_text()
    assert "verdict (auto): mixed" in text and "consistent with a TF-generic gene propensity" in text


# (b) TF-level + membership structure alone is not ------------------------------------------- #

def test_no_motif_mock_is_not_called_target_level(null):
    s = _arm(null)
    t1, t2 = _test(s, "T1"), _test(s, "T2")
    assert s["L"] > 0.1                                # a pooled lift from a prior with no target information
    assert s["components"]["tf_mean"] > 0.5 * s["L"]   # ... reached by TF-level reweighting alone
    assert not t1["p_a_gt_b_holm"] < 0.05              # T1 not supported
    assert abs(s["delta_target"]) < 0.5 * s["L"]        # target identity: below the rule's half
    assert not t2["p_a_gt_b_holm"] < 0.05              # T2 not supported
    assert s["macro_R_lift"] < 0.05                    # no within-TF lift on U'' (random ranks add noise)
    assert s["verdict"] == "design effects only"
    assert "verdict (auto): consistent with design effects only" in (null["out"] / "SENSITIVITY.md").read_text()


# (c) calibration gate ----------------------------------------------------------------------- #

def test_gate_passes_on_a_run_folder_written_by_real_bench(planted):
    st = planted["summary"]["sets"][EXP]
    assert st["gate_pass"] == st["units"] == N_TASKS and st["excluded"] == []
    assert st["gate2"]["pass"] and st["gate2"]["substituted"] == N_TASKS
    heads = {(c["a"], c["b"]) for c in st["gate2"]["checks"]}
    assert {("pearson_motif", "pearson_signed"), ("pearson_motif", "motif_prior")} <= heads
    for p in sorted((planted["out"] / "rows").glob("*.json")):
        row = json.loads(p.read_text())
        assert row["flags"] == [] and max(row["gate"].values()) == 0.0
        assert all(e["match"] and e["rows"] == CAP for e in row["edges"].values())   # truncated, matched
        assert sorted({c["ctx"] for c in row["ctx"] if c["arm"] == "auto"}) == ["ctx1"]


def _altered_unit(planted, tmp_path, alter) -> dict:
    run = tmp_path / "run"
    shutil.copytree(planted["run"], run)
    alter(run / "artifacts" / "cxg__mock__homo-sapiens__t0")
    units = ms.plan_units([f"{EXP}={run}"], planted["results"])
    return units[0]


def _run_one(planted, unit, monkeypatch):
    monkeypatch.setattr(real_bench, "ARTIFACT_TOP_EDGES", CAP)
    monkeypatch.setenv(ct.ENV_DIR, str(planted["cis"]))
    row, _ = ms.run_unit(unit, str(planted["refs"]), n_perm=2, workroot=str(planted["root"] / "work"))
    return row


def test_gate_fails_when_a_kept_processed_file_is_altered(planted, tmp_path, monkeypatch):
    def alter(tdir):
        path = tdir / "auto" / "ctx1" / "processed.h5ad"
        a = ad.read_h5ad(path)
        X = np.asarray(a.X).copy()
        X[:, :5] = np.random.default_rng(0).permutation(X[:, :5])       # scramble 5 TF columns
        a.X = X
        a.write_h5ad(path, compression="gzip")
    row = _run_one(planted, _altered_unit(planted, tmp_path, alter), monkeypatch)
    assert any(f.startswith("GATE score pearson_signed/auto") for f in row["flags"])
    assert not any("baseline" in f for f in row["flags"])               # the untouched arm still passes


def test_gate_fails_when_a_kept_edge_file_is_altered(planted, tmp_path, monkeypatch):
    def alter(tdir):
        path = tdir / "baseline" / "ctx0" / "pearson_signed.tsv.gz"
        e = ms.read_kept(path)
        e.loc[7, "weight"] += 1e-9
        e.to_csv(path, sep="\t", index=False, compression="gzip")
    row = _run_one(planted, _altered_unit(planted, tmp_path, alter), monkeypatch)
    assert [f for f in row["flags"] if f.startswith("GATE")] and \
        all(f.startswith("GATE edges baseline/ctx0/pearson_signed") for f in row["flags"])


def test_plan_notes_a_run_row_that_is_not_the_committed_row(planted, tmp_path):
    run = tmp_path / "run"
    shutil.copytree(planted["run"], run, ignore=shutil.ignore_patterns("artifacts"))
    p = run / "rows" / "cxg__mock__homo-sapiens__t1.json"
    p.write_text(p.read_text().replace('"reference_edges"', '"reference_edges_x"'))
    (run / "rows" / "cxg__mock__homo-sapiens__t2.json").unlink()        # gate 1 reads the committed row
    units = {u["task"]: u for u in ms.plan_units([f"{EXP}={run}"], planted["results"])}
    assert [n[:3] for n in units["cxg__mock__homo-sapiens__t1"]["plan_notes"]] == ["ROW"]
    assert [n[:6] for n in units["cxg__mock__homo-sapiens__t2"]["plan_notes"]] == ["no row"]
    assert not any(f.startswith("MISSING row") for f in units["cxg__mock__homo-sapiens__t2"]["plan_flags"])
    assert units["cxg__mock__homo-sapiens__t0"]["plan_notes"] == []
    assert all(any(f.startswith("MISSING artifacts") for f in u["plan_flags"]) for u in units.values())


def test_excluded_datasets_are_reported_not_dropped(planted, tmp_path):
    out = tmp_path / "sens"
    shutil.copytree(planted["out"], out)
    p = out / "rows" / f"{EXP}__cxg__mock__homo-sapiens__t2.json"
    row = json.loads(p.read_text())
    row["flags"] = ["GATE score pearson_signed/auto: |recomputed - committed| = 0.1"]
    p.write_text(json.dumps(row))
    text = ms.report(out)
    st = json.loads((out / "SENSITIVITY.json").read_text())["sets"][EXP]
    assert st["gate_pass"] == N_TASKS - 1 and [e["id"] for e in st["excluded"]] == [row["id"]]
    assert _test(st["arms"]["auto"], "T1")["n"] == N_TASKS - 1
    assert ("\n  - cxg__mock__homo-sapiens__t2 (completed): 1 reason(s)\n"
            "    - GATE score pearson_signed/auto: |recomputed - committed| = 0.1\n") in text
    # gate 2 cannot pass when a dataset fails gate 1 (record), but the tests run on the others
    assert st["gate2"]["pass"] is False and "Gate 2 cannot pass with 1 dataset(s)" in text
    assert st["not_interpretable"] is False and st["arms"]["auto"]["verdict"] != "not interpretable"


def test_every_gate_1_failure_is_listed_with_its_size_and_a_cut_is_marked(planted, tmp_path):
    """Record, platform note: the report prints the size of every gate-1 failure. Six score and six
    edge-file failures (the formats run_unit writes) are each printed in full; only free text longer
    than 400 characters is cut, and the cut says so."""
    out = tmp_path / "sens"
    shutil.copytree(planted["out"], out)
    rng = np.random.default_rng(0)
    gate = [f"GATE score {m}/{arm}: |recomputed - committed| = {float(rng.uniform(1e-12, 9e-12))}"
            for m in ms.RECOMPUTABLE for arm in ms.ARMS]
    gate += [f"GATE edges {arm}/ctx{i}/{m}: " + str({"match": False, "rows": 3000, "same_pairs_in_order": True,
                                                       "max_abs_weight_diff": float(rng.uniform(1e-12, 9e-12))})
             for i, arm in enumerate(ms.ARMS) for m in ms.RECOMPUTABLE]
    error = "ERROR MemoryError: " + "x" * 481                                  # 500 characters, as _run_and_save cuts
    for task, status, flags in (("t2", "completed", gate), ("t3", "error", [error])):
        p = out / "rows" / f"{EXP}__cxg__mock__homo-sapiens__{task}.json"
        row = json.loads(p.read_text())
        row["status"], row["flags"] = status, flags
        p.write_text(json.dumps(row))
    text = ms.report(out)
    assert len(gate) == 12 and all(len(f) < 400 for f in gate) and len("; ".join(gate)) > 400
    lines = text.splitlines()
    i = lines.index("  - cxg__mock__homo-sapiens__t2 (completed): 12 reason(s)")
    assert lines[i + 1:i + 13] == [f"    - {f}" for f in gate]                 # every size, in full, in order
    j = lines.index("  - cxg__mock__homo-sapiens__t3 (error): 1 reason(s)")
    assert lines[j + 1] == f"    - {error[:400]} [... 100 characters cut; full text in SENSITIVITY.json]"
    reasons = {e["id"].rsplit("__", 1)[1]: e["reasons"] for e in
               json.loads((out / "SENSITIVITY.json").read_text())["sets"][EXP]["excluded"]}
    assert reasons == {"t2": gate, "t3": [error]}                             # the JSON keeps every reason whole


def test_gate2_failure_with_every_dataset_through_gate_1_is_not_interpretable(planted, tmp_path, monkeypatch):
    out = tmp_path / "sens"
    shutil.copytree(planted["out"], out)
    bad = {"heading": "#### Declared comparisons, AUPRC ratio, auto pipeline", "a": "pearson_motif",
           "b": "pearson_signed", "committed": "| x |", "recomputed": "| y |", "equal": False}
    monkeypatch.setattr(ms, "gate2", lambda rows, committed_dir: ([bad], N_TASKS, N_TASKS))
    text = ms.report(out)
    st = json.loads((out / "SENSITIVITY.json").read_text())["sets"][EXP]
    assert st["gate_pass"] == N_TASKS and st["not_interpretable"] is True
    assert {a["verdict"] for a in st["arms"].values()} == {"not interpretable"}
    assert "verdict (auto): not interpretable" in text and "supports target-level" not in text


def test_a_planned_unit_without_a_row_is_reported_as_missing(planted, tmp_path):
    """A killed worker or a job cut by its time limit writes no row: exclusion rule (c) reports and
    excludes the unit; it is never dropped and never blamed on the report code path (gate 2)."""
    assert planted["summary"]["missing_units"] == []
    assert "- units without a row: 0\n" in (planted["out"] / "SENSITIVITY.md").read_text()
    out = tmp_path / "sens"
    shutil.copytree(planted["out"], out)
    uid = f"{EXP}__cxg__mock__homo-sapiens__t2"
    (out / "rows" / f"{uid}.json").unlink()
    text = ms.report(out)
    summary = json.loads((out / "SENSITIVITY.json").read_text())
    st = summary["sets"][EXP]
    assert summary["missing_units"] == [uid]
    assert st["units"] == N_TASKS and st["gate_pass"] == N_TASKS - 1
    assert [(e["id"], e["status"]) for e in st["excluded"]] == [(uid, "missing")]
    assert st["excluded"][0]["reasons"][0].startswith(f"MISSING row: no rows/{uid}.json was written")
    assert st["not_interpretable"] is False and st["arms"]["auto"]["verdict"] != "not interpretable"
    assert _test(st["arms"]["auto"], "T1")["n"] == N_TASKS - 1
    assert text.startswith(f"# {ms.EXPERIMENT}") and "\n**INCOMPLETE: 1 planned unit(s) have no row" in text
    assert f"- units planned: {N_TASKS}; rows: {N_TASKS - 1};" in text
    assert f"- units without a row: 1 ({uid})\n" in text
    assert "  - cxg__mock__homo-sapiens__t2 (missing): 1 reason(s)\n    - MISSING row: no rows/" in text
    assert "excluded (reported, never silently): 1" in text and "Gate 2 cannot pass with 1 dataset(s)" in text
    assert "does not reproduce" not in text
    # a set whose every row is missing is still reported, with every unit excluded
    for p in (out / "rows").glob("*.json"):
        p.unlink()
    text = ms.report(out)
    st = json.loads((out / "SENSITIVITY.json").read_text())["sets"][EXP]
    assert st["units"] == N_TASKS and st["gate_pass"] == 0 and st["not_interpretable"] is False
    assert {e["status"] for e in st["excluded"]} == {"missing"} and len(st["excluded"]) == N_TASKS
    assert "Not tested: no dataset of this set passed gate 1." in text
    # without sens_tasks.json the number of planned units is unknown, and nothing is invented
    (out / "sens_tasks.json").unlink()
    text = ms.report(out)
    assert "- units planned: unknown; rows: 0;" in text and "INCOMPLETE" not in text
    assert json.loads((out / "SENSITIVITY.json").read_text())["missing_units"] is None


def test_a_failed_code_check_makes_the_set_not_interpretable(planted, tmp_path):
    out = tmp_path / "sens"
    shutil.copytree(planted["out"], out)
    p = out / "rows" / f"{EXP}__cxg__mock__homo-sapiens__t3.json"
    row = json.loads(p.read_text())
    p.write_text(json.dumps({k: row[k] for k in ("id", "exp", "task")} |
                            {"status": "check_failed", "flags": ["CHECK CheckFailed: built-in checks: x"]}))
    text = ms.report(out)
    st = json.loads((out / "SENSITIVITY.json").read_text())["sets"][EXP]
    assert st["not_interpretable"] is True and st["arms"]["auto"]["verdict"] == "not interpretable"
    assert "A code check failed on 1 dataset(s)" in text
    assert ("check_failed units are not resubmitted; after the fix every unit is rerun at the fixed commit into "
            "a new results folder and the fix is a dated deviation naming both commits") in text
    assert "no resubmission with the same commit" not in text


def test_rows_carry_no_home_path_and_check_failures_are_not_rerun(tmp_path, monkeypatch):
    home = tmp_path / "home" / "someone"
    home.mkdir(parents=True)
    monkeypatch.setattr(ms.Path, "home", lambda: home)
    unit = {"id": "u1", "exp": EXP, "task": "t"}

    def boom(exc):
        def run_unit(*a, **k):
            raise exc
        return run_unit
    monkeypatch.setattr(ms, "run_unit", boom(FileNotFoundError(f"{home}/grnbench_artifacts/x/processed.h5ad")))
    job = (unit, None, str(tmp_path / "o"), 2, None, "c0ffee0")
    assert ms._run_and_save(job)[1].startswith("error")
    text = (tmp_path / "o" / "rows" / "u1.json").read_text()
    assert str(home) not in text and "~/grnbench_artifacts/x/processed.h5ad" in text
    assert json.loads(text)["code_commit"] == "c0ffee0"
    monkeypatch.setattr(ms, "run_unit", boom(ms.CheckFailed("self-consistency failed")))
    assert ms._run_and_save(job)[1].startswith("check_failed")
    assert json.loads((tmp_path / "o" / "rows" / "u1.json").read_text())["code_commit"] == "c0ffee0"
    assert ms._run_and_save(job)[1] == "skipped (check_failed)"


def _fake_run(planted, tmp_path, monkeypatch, kinds):
    """Units u_<kind> run through main() with a stand-in run_unit (no computation): ok completes, err
    raises a runtime error while state['err_fails'], bad fails a code check, excl is excluded by the
    plan (the real run_unit). state['commit'] is the code commit the job runs."""
    base = {"exp": EXP, "pairs": [list(p) for p in ms.STAGE1_PAIRS], "recompute": [], "plan_notes": [],
            "committed": str(tmp_path / "no_committed_rows")}
    units = [{**base, "id": f"u_{k}", "task": f"t_{k}", "plan_flags": ["MISSING artifacts: x"] if k == "excl" else []}
             for k in kinds]
    state, real_run_unit = {"commit": "aaaaaaa", "err_fails": True}, ms.run_unit

    def run_unit(unit, refs, n_perm=ms.N_PERM, workroot=None):
        if unit["id"] == "u_excl":
            return real_run_unit(unit, refs, n_perm, workroot)
        if unit["id"] == "u_bad":
            raise ms.CheckFailed("fast scorer != evaluate_ranking")
        if unit["id"] == "u_err" and state["err_fails"]:
            raise OSError("disk quota exceeded")
        return {"id": unit["id"], "exp": EXP, "task": unit["task"], "status": "completed", "flags": [],
                "n_perm": n_perm}, {}
    monkeypatch.setattr(ms, "run_unit", run_unit)
    monkeypatch.setattr(ms, "code_commit", lambda: state["commit"])
    monkeypatch.setenv(ct.ENV_DIR, str(planted["cis"]))
    tasks = tmp_path / f"tasks_{'_'.join(kinds)}.json"
    tasks.write_text(json.dumps(units))

    def main(out, refs=planted["refs"], n_perm=2):
        return ms.main(["run", "--tasks", str(tasks), "--out", str(out), "--refs", str(refs), "--n-perm", str(n_perm)])
    return state, main


def _rows_of(out):
    return {p.stem: json.loads(p.read_text()) for p in sorted((out / "rows").glob("*.json"))}


def _bytes(out):
    return {p.relative_to(out).as_posix(): p.read_bytes() for p in sorted(out.rglob("*")) if p.is_file()}


def test_every_row_records_its_commit_and_run_refuses_to_mix_commits_in_one_out(planted, tmp_path, monkeypatch,
                                                                                  capsys):
    """Rule (c) resubmits with the same commit: a resubmission at that commit reruns the error row; a
    job at another commit (a re-sync, or a code fix) is refused before any row or RUN_INPUTS.json is
    touched, so one results folder never holds rows of two commits."""
    state, run = _fake_run(planted, tmp_path, monkeypatch, ("ok", "err", "excl"))
    out = tmp_path / "out"
    assert run(out) == 2                                              # u_err: runtime error
    rows = _rows_of(out)
    assert {k: (r["status"], r["code_commit"]) for k, r in rows.items()} == {
        "u_ok": ("completed", "aaaaaaa"), "u_err": ("error", "aaaaaaa"), "u_excl": ("excluded", "aaaaaaa")}
    inputs = json.loads((out / "RUN_INPUTS.json").read_text())
    assert [i["code_commit"] for i in inputs["invocations"]] == ["aaaaaaa"]
    assert inputs["invocations"][0]["platform"] == inputs["platform"]
    capsys.readouterr()
    state.update(commit="bbbbbbb", err_fails=False)                   # code re-synced: refused, nothing touched
    before = _bytes(out)
    assert run(out) == 2
    err = capsys.readouterr().err
    assert "holds 3 row(s) from commit aaaaaaa" in err and "this job runs bbbbbbb" in err
    assert "rerun every unit at the fixed commit into a new, empty --out" in err
    assert "has already copied this job's GIT_COMMIT to OUT/CODE_COMMIT" in err    # the sbatch ran first
    assert "echo aaaaaaa > OUT/CODE_COMMIT" in err
    assert _bytes(out) == before
    state["commit"] = "aaaaaaa"                                       # the same commit: the error row is rerun
    assert run(out) == 0
    rows = _rows_of(out)
    assert {k: (r["status"], r["code_commit"]) for k, r in rows.items()} == {
        "u_ok": ("completed", "aaaaaaa"), "u_err": ("completed", "aaaaaaa"), "u_excl": ("excluded", "aaaaaaa")}
    again = json.loads((out / "RUN_INPUTS.json").read_text())
    assert [i["code_commit"] for i in again["invocations"]] == ["aaaaaaa", "aaaaaaa"]
    assert {k: again[k] for k in ("references", "cistarget", "platform")} == \
        {k: inputs[k] for k in ("references", "cistarget", "platform")}  # kept, not overwritten
    p = out / "rows" / "u_ok.json"                                    # a row written before rows had a commit
    p.write_text(json.dumps({k: v for k, v in json.loads(p.read_text()).items() if k != "code_commit"}))
    before = _bytes(out)
    capsys.readouterr()
    err = ""
    assert run(out) == 2
    err = capsys.readouterr().err
    assert "holds 1 row(s) from commit none recorded" in err
    assert "OUT/CODE_COMMIT" not in err          # no commit to restore when the old row recorded none
    assert _bytes(out) == before


def test_a_skipped_check_failed_unit_exits_2_with_the_rerun_rule(planted, tmp_path, monkeypatch, capsys):
    state, run = _fake_run(planted, tmp_path, monkeypatch, ("ok", "bad"))
    out = tmp_path / "out"
    assert run(out) == 2
    assert _rows_of(out)["u_bad"]["status"] == "check_failed" and _rows_of(out)["u_bad"]["code_commit"] == "aaaaaaa"
    capsys.readouterr()
    assert run(out) == 2                                              # never 0 as if the fix had run
    cap = capsys.readouterr()
    assert "u_bad: skipped (check_failed)" in cap.out
    assert "a check_failed unit is never resubmitted" in cap.err and "naming both commits" in cap.err


def test_run_inputs_are_kept_on_resubmission_and_refused_when_they_change(planted, tmp_path, monkeypatch, capsys):
    state, run = _fake_run(planted, tmp_path, monkeypatch, ("ok",))
    out = tmp_path / "out"
    assert run(out) == 0
    first = (out / "RUN_INPUTS.json").read_bytes()
    refs = tmp_path / "refs_changed"
    shutil.copytree(planted["refs"], refs)
    with open(refs / "collectri_human.tsv", "a") as fh:
        fh.write("TFX\tGX\n")
    capsys.readouterr()
    assert run(out, refs=refs) == 2
    err = capsys.readouterr().err
    assert "inputs changed between submissions: references/collectri_human.tsv differ" in err
    assert (out / "RUN_INPUTS.json").read_bytes() == first
    assert run(out) == 0                                              # the same inputs: kept, one more invocation
    inputs = json.loads((out / "RUN_INPUTS.json").read_text())
    assert len(inputs["invocations"]) == 2 and inputs["references"] == json.loads(first)["references"]


def test_report_names_the_code_commits_and_the_registered_inputs_it_lacks(planted, tmp_path):
    text = (planted["out"] / "SENSITIVITY.md").read_text()
    st = planted["summary"]["sets"][EXP]
    commit = {json.loads(p.read_text())["code_commit"] for p in (planted["out"] / "rows").glob("*.json")}
    assert st["code_commits"] == sorted(commit) and len(commit) == 1 and st["mixed_code_commits"] is False
    assert f"- code commit of the rows: {commit.pop()}\n" in text and "MIXED CODE COMMITS" not in text
    assert "- CODE_COMMIT: not in this folder\n" in text and planted["summary"]["code_commit_file"] is None
    snapshot = ms._snapshot_inputs()          # the mock reads a collectri_human.tsv (another hash) and nothing else
    assert len(snapshot) == 8 and planted["summary"]["missing_inputs"] == \
        [p for p in snapshot if p != "references/collectri_human.tsv"]
    for p in planted["summary"]["missing_inputs"]:
        assert f"- input {p}: MISSING from RUN_INPUTS.json" in text
    assert "- input references/collectri_human.tsv: sha256 " in text and "invocations" not in text
    out = tmp_path / "sens"
    shutil.copytree(planted["out"], out)
    inputs = json.loads((out / "RUN_INPUTS.json").read_text())
    inputs["references"]["collectri_mouse.tsv"] = "0" * 64           # recorded (with another hash): not MISSING
    (out / "RUN_INPUTS.json").write_text(json.dumps(inputs))
    p = out / "rows" / f"{EXP}__cxg__mock__homo-sapiens__t2.json"
    row = json.loads(p.read_text())
    p.write_text(json.dumps({**row, "code_commit": "0000000"}))
    text = ms.report(out)
    summary = json.loads((out / "SENSITIVITY.json").read_text())
    assert summary["sets"][EXP]["mixed_code_commits"] is True and len(summary["sets"][EXP]["code_commits"]) == 2
    assert "- **MIXED CODE COMMITS**: the rows of this set come from 0000000, " in text
    assert "references/collectri_mouse.tsv" not in summary["missing_inputs"] and len(summary["missing_inputs"]) == 6
    assert "- input references/collectri_mouse.tsv: sha256 0000000000000000... (NOT the registered hash)" in text
    p.write_text(json.dumps(row))                                     # one commit, but CODE_COMMIT names another
    (out / "CODE_COMMIT").write_text("1111111\n")
    text = ms.report(out)
    assert "- CODE_COMMIT (copied by hpc/engine_bench.sbatch at the last submission): 1111111" in text
    assert "MIXED CODE COMMITS" in text and " and CODE_COMMIT is 1111111." in text


# (d) the built-in checks of PLAN.md section 2 ------------------------------------------------ #

def _contexts(mock):
    for p in sorted((mock["out"] / "rows").glob("*.json")):
        yield from json.loads(p.read_text())["ctx"]


@pytest.mark.parametrize("which", ["planted", "null", "gene"])
def test_constant_row_and_tie_controls_cannot_move_their_metrics(which, request):
    mock = request.getfixturevalue(which)
    n = 0
    for c in _contexts(mock):
        n += 1
        for v in ("tf_mean", "tf_kexp", "tf_presence"):                 # constant rows: no macro lift
            for m in ("macro", "macro_T", "macro_R"):
                assert c[f"{v}_{m}"] == c[f"base_{m}"], (v, m)
        for v in ("tf_presence", "db_only", "worst_tie"):               # ties on U'': no lift at all
            for m in ("pooled_R", "macro_R"):
                assert c[f"{v}_{m}"] == c[f"base_{m}"], (v, m)
        for v in ("tf_mean", "tfmean_db"):
            assert c[f"{v}_macro_R"] == c["base_macro_R"], v
        assert c["n_tfs_no_motif"] == N_NO_MOTIF and c["genes_ensembl"] == N_ENS and c["gold_in_worst_tie"] >= 0
        assert c["builtin_checks"] == "held" and c["n_tfs_macro_R"] <= c["n_tfs_macro"] <= c["n_tfs"]
    assert n == 2 * N_TASKS                                            # datasets x 2 arms


def test_builtin_check_violations_are_code_errors():
    rec = {f"{v}_{m}": 1.0 for v in ("base",) + ms.FIXED for m in ms.METRICS}
    assert ms.builtin_check_violations(rec) == []
    assert ms.builtin_check_violations({**rec, "db_only_pooled_R": 1.5}) == ["db_only_pooled_R 1.5 != base_pooled_R 1.0"]
    assert ms.builtin_check_violations({**rec, "worst_tie_macro_R": 2.0}, worst_tie_applies=False) == []
    assert ms.builtin_check_violations({**rec, "base_macro_R": np.nan, "tf_mean_macro_R": np.nan,
                                        "tf_kexp_macro_R": np.nan, "tf_presence_macro_R": np.nan,
                                        "db_only_macro_R": np.nan, "worst_tie_macro_R": np.nan,
                                        "tfmean_db_macro_R": np.nan}) == []          # NaN equals NaN


def test_keyed_permutation_is_one_random_database_shared_by_every_dataset():
    rng = np.random.default_rng(0)
    tfs = ["TF1", "TF2", "TF3"]
    pool = [f"G{i}" for i in range(40)] + [f"H{i}" for i in range(15)]
    value = {(t, g): float(v) for t in tfs for g, v in zip(pool, rng.permutation(len(pool)) + 1)}
    genes_a = pool[:40]
    genes_b = list(rng.permutation(pool[10:]))                         # 30 shared genes, another order
    K, worst = np.array([2, 1, 0]), 1000

    def perm(genes, seed, name="perm_db"):
        M = np.array([[value[(t, g)] for g in genes] for t in tfs])
        in_db = np.array([not g.endswith("9") for g in genes])
        P = ms.variant_matrix(name, M, in_db, K, worst, seed, (ms.name_keys(tfs), ms.name_keys(genes)))
        return M, in_db, P
    for seed in range(3):
        Ma, ia, Pa = perm(genes_a, seed)
        Mb, ib, Pb = perm(genes_b, seed)
        for M, in_db, P in ((Ma, ia, Pa), (Mb, ib, Pb)):
            assert np.array_equal(P[:, ~in_db], M[:, ~in_db])                         # mask kept
            assert np.array_equal(np.sort(P[:, in_db], 1), np.sort(M[:, in_db], 1))   # TF values kept
            assert not np.array_equal(P, M)
        shared = [g for g in genes_a if g in genes_b and not g.endswith("9")]
        for i in range(len(tfs)):                                      # the same order in both datasets
            ra = [Pa[i, genes_a.index(g)] for g in shared]
            rb = [Pb[i, genes_b.index(g)] for g in shared]
            assert np.array_equal(np.argsort(ra), np.argsort(rb))
    assert not np.array_equal(perm(genes_a, 0)[2], perm(genes_a, 1)[2])            # seeds differ
    M, _, P = perm(genes_a, 0, "perm_all")
    assert np.array_equal(np.sort(P, 1), np.sort(M, 1))
    with pytest.raises(ValueError, match="seed"):
        ms.variant_matrix("perm_db", M, np.ones(40, bool), K, worst, 0)


def test_gene_controls_keep_the_tf_values_and_follow_the_other_tfs():
    M = np.array([[1.0, 5.0, 3.0, 9.0], [2.0, 6.0, 4.0, 9.0], [7.0, 8.0, 1.0, 9.0], [9.0, 9.0, 9.0, 9.0]])
    in_db, K = np.array([True, True, True, False]), np.array([1, 2, 1, 0])
    G = ms.variant_matrix("gene_tf", M, in_db, K, 9)
    assert np.array_equal(np.sort(G[:, :3], 1)[:3], np.sort(M[:3, :3], 1))       # each TF keeps its values
    assert (G[:, 3] == 9).all() and (G[3] == 9).all()                          # mask and motif-less TF
    # TF 0's leave-one-out gene means over TFs 1, 2: [4.5, 7.0, 2.5] -> its best value goes to gene 2
    assert list(G[0, :3]) == [3.0, 5.0, 1.0]
    assert list(ms.variant_matrix("gene_mean", M, in_db, K, 9)[0, :3]) == [4.5, 7.0, 2.5]


# script mechanics --------------------------------------------------------------------------- #

def test_local_rank_product_equals_the_registered_one(planted, monkeypatch):
    monkeypatch.setenv(ct.ENV_DIR, str(planted["cis"]))
    from grnbench.infer import prepare_expression, run_method
    tdir = planted["run"] / "artifacts" / "cxg__mock__homo-sapiens__t3"
    adata = ad.read_h5ad(tdir / "baseline" / "ctx0" / "processed.h5ad")
    tfs = sorted(set(pd.read_csv(tdir / "reference_evaluation.tsv", sep="\t")["source"]))
    E, _ = run_method("pearson_signed", adata, tf_list=tfs, assume_normalized=True)
    expr, _ = prepare_expression(adata, assume_normalized=True)
    tfs_e = list(dict.fromkeys(E["source"]))
    M, in_db, K, worst = ms.motif_parts(tfs_e, tuple(expr.columns))
    m = M[E["source"].map({t: i for i, t in enumerate(tfs_e)}).to_numpy(),
          E["target"].map({g: i for i, g in enumerate(expr.columns)}).to_numpy()].astype(float)
    mine, ref = ms.rank_product(E, m), ct.motif_rank_product(E, expr)
    assert mine[["source", "target"]].equals(ref[["source", "target"]])
    assert np.array_equal(mine["weight"].to_numpy(), ref["weight"].to_numpy())
    assert worst == N_BG + N_SYM + N_TF - N_OUT_DB and int((~in_db).sum()) == N_OUT_DB + N_ENS
    assert (K == 0).sum() == N_NO_MOTIF


def test_stage2_models_are_refused_with_the_reason(planted, tmp_path, capsys):
    with pytest.raises(ms.Stage2Refused, match="GRNBoost2 calibration gate"):
        ms.plan_units([f"{EXP}={planted['run']}"], planted["results"], [("grnboost2_motif", "grnboost2")])
    assert ms.main(["plan", "--run", f"{EXP}={planted['run']}", "--out", str(tmp_path / "t.json"),
                    "--results-root", str(planted["results"]), "--pairs", "regdiffusion_motif:regdiffusion"]) == 2
    assert "stage 2" in capsys.readouterr().err and not (tmp_path / "t.json").exists()
    unit = {**ms.plan_units([f"{EXP}={planted['run']}"], planted["results"])[0],
            "pairs": [["grnboost2_motif", "grnboost2"]]}
    with pytest.raises(ms.Stage2Refused):
        ms.run_unit(unit, str(planted["refs"]))
    tasks = tmp_path / "tasks.json"
    tasks.write_text(json.dumps([unit]))
    assert ms.main(["run", "--tasks", str(tasks), "--out", str(tmp_path / "o"), "--refs", str(planted["refs"])]) == 2
    assert ms.main(["run", "--tasks", str(tasks), "--out", str(tmp_path / "o"), "--refs", str(planted["refs"]),
                    "--models", "regdiffusion+regdiffusion_motif"]) == 2
    assert not (tmp_path / "o" / "rows").exists()


def test_run_refuses_before_any_row_when_the_cistarget_folder_or_a_reference_is_missing(planted, tmp_path,
                                                                                        monkeypatch, capsys):
    """A typo in the --export of the submit line (GRNBENCH_CISTARGET_DIR or REFS) stops `run` before
    any unit runs: no row, no RUN_INPUTS.json, and an existing RUN_INPUTS.json keeps its bytes."""
    tasks = str(planted["out"] / "sens_tasks.json")
    monkeypatch.setattr(real_bench, "ARTIFACT_TOP_EDGES", CAP)

    def run(out, refs=planted["refs"]):
        return ms.main(["run", "--tasks", tasks, "--out", str(out), "--refs", str(refs), "--n-perm", "2",
                        "--workroot", str(planted["root"] / "work")])
    monkeypatch.delenv(ct.ENV_DIR, raising=False)
    assert run(tmp_path / "fresh") == 2
    err = capsys.readouterr().err
    assert "preflight failed (no unit was run, no row written)" in err and f"{ct.ENV_DIR} is unset" in err
    assert not (tmp_path / "fresh").exists()
    out = tmp_path / "resumed"                                     # a finished run, resubmitted with the typo
    shutil.copytree(planted["out"], out)
    before = {p.name: p.read_bytes() for p in [out / "RUN_INPUTS.json", *sorted((out / "rows").glob("*.json"))]}
    monkeypatch.setenv(ct.ENV_DIR, str(tmp_path / "no_such_folder"))
    assert run(out) == 2 and f"{ct.ENV_DIR} is unset or not a folder" in capsys.readouterr().err
    assert before == {p.name: p.read_bytes() for p in [out / "RUN_INPUTS.json", *sorted((out / "rows").glob("*.json"))]}
    empty = tmp_path / "empty"
    empty.mkdir()
    monkeypatch.setenv(ct.ENV_DIR, str(empty))                     # a folder without the human files
    assert run(tmp_path / "o2") == 2
    assert "no cisTarget rankings database and motif table for human" in capsys.readouterr().err
    monkeypatch.setenv(ct.ENV_DIR, str(planted["cis"]))
    assert run(tmp_path / "o3", refs=empty) == 2                   # --refs without collectri_human.tsv
    assert "no collectri_human.tsv" in capsys.readouterr().err
    assert not any((tmp_path / d).exists() for d in ("o2", "o3"))
    assert ms.preflight_problems(json.loads(Path(tasks).read_text()), str(planted["refs"])) == []


def test_verdict_rule_precedence():
    def tests(t1, t1_rev, t2, t3=True, t2_rev=False):
        mk = lambda name, a, b: {"test": name, "res": {}, "p_a_gt_b_holm": a, "p_b_gt_a_holm": b}  # noqa: E731
        return [mk("T1", 0.01 if t1 else 0.5, 0.01 if t1_rev else 1.0),
                mk("T2", 0.01 if t2 else 0.5, 0.01 if t2_rev else 1.0), mk("T3", 0.01 if t3 else 0.5, 1.0)]
    def v(tests_, L, d, art, gene=0.02, robust=True):
        return ms.stage1_verdict(tests_, L, d, art, gene, t1_robust=robust)
    assert v(tests(True, False, True), 0.2, 0.15, 0.05) == "target-level"
    assert v(tests(True, False, True), 0.2, 0.15, 0.12) == "mixed"            # half reached without targets
    assert v(tests(True, False, True), 0.2, 0.05, 0.05) == "mixed"            # target identity < half
    assert v(tests(False, False, True), 0.2, 0.05, 0.12) == "mixed"
    assert v(tests(False, False, True), 0.2, 0.05, 0.05) == "inconclusive"
    assert v(tests(False, True, False), 0.2, -0.1, 0.0) == "design effects only"
    assert v(tests(False, False, False), 0.2, 0.0, 0.15) == "design effects only"
    assert v(tests(False, False, False), 0.2, 0.0, 0.05) == "inconclusive"
    assert v(tests(True, False, True), -0.01, 0.15, 0.0) == "inconclusive"    # no lift to decompose
    # robustness condition (record, GB 2026-10-04): target-level needs T1's permutation rank p < 0.05
    assert v(tests(True, False, True), 0.2, 0.15, 0.05, robust=False) == "mixed"
    assert v(tests(False, False, False), 0.2, 0.0, 0.15, robust=False) == "design effects only"  # unaffected
    # gene-level condition (record, GB 2026-10-05): target-level also needs T3 and L_gene < 0.5 L
    assert v(tests(True, False, True, t3=False), 0.2, 0.15, 0.05) == "mixed"
    assert v(tests(True, False, True), 0.2, 0.15, 0.05, gene=0.12) == "mixed"
    assert v(tests(True, False, True), 0.2, 0.15, 0.05, gene=float("nan")) == "mixed"
    assert v(tests(False, False, True, t3=False), 0.2, 0.05, 0.05) == "inconclusive"   # T3 only in rule 3
    # T2's reverse joins T1's: the prior lowers within-TF ranking where no design effect can act
    assert v(tests(True, False, False, t2_rev=True), 0.2, 0.15, 0.0) == "design effects only"   # T1 and T3 too
    # ... reported when T1 or T3 is supported, with the verdict the rules give without T2's reverse
    def via_t2(tests_, art):
        return ms.t2_reverse_only(tests_, 0.2, 0.15, art, 0.02, t1_robust=True)
    assert via_t2(tests(True, False, False, t2_rev=True), 0.0) == (True, "inconclusive")    # rule 4, L_art < 0.5 L
    assert via_t2(tests(True, False, False, t2_rev=True), 0.15) == (True, "mixed")          # rule 4, L_art >= 0.5 L
    assert via_t2(tests(False, False, False, t2_rev=True), 0.15) == (True, "design effects only")  # T3 alone: rule 5
    assert via_t2(tests(False, False, False, t3=False, t2_rev=True), 0.15) == (False, None)  # nothing supported
    assert via_t2(tests(False, True, False, t2_rev=True), 0.0) == (False, None)              # T1 reversed as well
    assert via_t2(tests(True, False, True), 0.0) == (False, None)                           # rule 2 did not fire
    short =ms.target_level_shortfalls(tests(True, False, True, t3=False), 0.2, 0.15, 0.05, 0.12, t1_robust=False)
    assert [x.split(" ")[0] for x in short] == ["T3", "L_gene", "T1"]


def test_gate2_reproduces_the_committed_declared_lines_of_c4_and_c6():
    """On the real committed rows: substituting each dataset's own committed values must print the
    committed REPORT_REAL.md lines (checks the gate's parsing and the registered report code)."""
    for exp, pairs in (("EXP-20261002-LM-1", {("pearson_motif", "pearson_signed")}),
                       ("EXP-20260925-LM-2", {("pearson_motif", "pearson_signed"), ("pearson_motif", "motif_prior")})):
        committed = ROOT / "benchmarks" / "results" / exp
        recompute = ["pearson_signed", "pearson_motif"] + (["motif_prior"] if exp.endswith("LM-2") else [])
        rows = []
        for p in sorted((committed / "rows").glob("*.json")):
            r = json.loads(p.read_text(encoding="utf-8"))
            rows.append({"task": r["task"], "recompute": recompute,
                         "registered": {k: v for k, v in r.items()
                                        if k.startswith(tuple(f"m_{m}_" for m in recompute))}})
        checks, n_sub, n_set = ms.gate2(rows, committed)
        assert n_sub == n_set == 20 and checks and all(c["equal"] for c in checks), exp
        assert pairs <= {(c["a"], c["b"]) for c in checks}
        auto = next(c for c in checks if c["a"] == "pearson_motif" and c["b"] == "pearson_signed"
                    and "auto pipeline" in c["heading"])
        assert "| 20 |" in auto["committed"]
        rows[0]["registered"]["m_pearson_motif_auto_auprc_ratio"] += 0.01
        assert not all(c["equal"] for c in ms.gate2(rows, committed)[0])


# (e) the declared tests and decision quantities are the record's definitions ---------------- #

def _rows(mock):
    return [json.loads(p.read_text()) for p in sorted((mock["out"] / "rows").glob("*.json"))]


def _mean_diff(rows, arm, a, b, combo="pearson_motif"):
    d = [r["means"][f"{arm}|{combo}|{a}"] - r["means"][f"{arm}|{combo}|{b}"] for r in rows]
    return float(np.mean(d))


def _holm_by_hand(ps):
    order = sorted(range(len(ps)), key=lambda i: ps[i])
    out, run = [0.0] * len(ps), 0.0
    for rank, i in enumerate(order):
        run = max(run, (len(ps) - rank) * ps[i])
        out[i] = min(1.0, run)
    return out


def test_holm_runs_over_the_six_directed_hypotheses_and_shrinks_by_two_per_untested_test():
    from grnbench.eval.comparison import compare_units
    rng = np.random.default_rng(4)
    n = 7
    cols = {"real_pooled": 1.5 + rng.normal(0, 0.1, n), "perm_db_pooled": 1.4 + rng.normal(0, 0.1, n),
            "real_macro_R": 1.30 + rng.normal(0, 0.1, n), "base_macro_R": 1.20 + rng.normal(0, 0.1, n),
            "gene_tf_macro_R": 1.25 + rng.normal(0, 0.1, n)}

    def mk(drop=()):                       # rows; the metrics in `drop` have a value in 2 datasets only
        return [{"task": f"t{i}", "means": {f"auto|pm|{k}": (float("nan") if k in drop and i >= 2 else float(v[i]))
                                            for k, v in cols.items()}} for i in range(n)]

    def raw(a, b):
        r = compare_units(pd.DataFrame({"a": cols[a], "b": cols[b]}))["pairs"][0]
        return [r["p_a_gt_b"], r["p_b_gt_a"]]

    def flat(got):
        return [x for t in got if t["res"] is not None for x in (t["p_a_gt_b_holm"], t["p_b_gt_a_holm"])]
    ps = raw("real_pooled", "perm_db_pooled") + raw("real_macro_R", "base_macro_R") + raw("real_macro_R", "gene_tf_macro_R")
    want = _holm_by_hand(ps)
    per_test = sum((_holm_by_hand(ps[i:i + 2]) for i in range(0, 6, 2)), [])
    assert not np.allclose(want, per_test)                 # the data can tell a family of 6 from three of 2
    got = ms.declared_tests(mk(), "auto", "pm")
    assert flat(got) == pytest.approx(want, abs=1e-12) and {t["holm_family"] for t in got} == {6}
    # T3 with only 2 datasets with a value: not run, the family is the 4 hypotheses of T1 and T2
    got4 = ms.declared_tests(mk(drop=("gene_tf_macro_R",)), "auto", "pm")
    assert [t["res"] is None for t in got4] == [False, False, True] and {t["holm_family"] for t in got4} == {4}
    assert flat(got4) == pytest.approx(_holm_by_hand(ps[:4]), abs=1e-12)
    # T2 and T3 untested (macro_R of the real prior in 2 datasets only): T1's 2 hypotheses alone
    got2 = ms.declared_tests(mk(drop=("real_macro_R",)), "auto", "pm")
    assert [t["res"] is None for t in got2] == [False, True, True] and {t["holm_family"] for t in got2} == {2}
    assert flat(got2) == pytest.approx(_holm_by_hand(ps[:2]), abs=1e-12)


@pytest.mark.parametrize("which", ["planted", "gene", "null"])
def test_declared_tests_and_decision_quantities_are_the_registered_definitions(which, request):
    """T1-T3, L, D, L_art, L_gene, the components and the ratios, recomputed from the per-dataset
    means of the rows with the definitions of section 1 of the record."""
    mock = request.getfixturevalue(which)
    rows = _rows(mock)
    for arm in ms.ARMS:
        s = _arm(mock, arm)
        for name, a, b in (("T1", "real_pooled", "perm_db_pooled"), ("T2", "real_macro_R", "base_macro_R"),
                           ("T3", "real_macro_R", "gene_tf_macro_R")):
            assert _test(s, name)["difference"] == pytest.approx(_mean_diff(rows, arm, a, b), abs=1e-12), name
        L = _mean_diff(rows, arm, "real_pooled", "base_pooled")
        assert s["L"] == pytest.approx(L, abs=1e-12)
        assert s["delta_target"] == pytest.approx(_mean_diff(rows, arm, "real_pooled", "perm_db_pooled"), abs=1e-12)
        assert s["L_art"] == pytest.approx(_mean_diff(rows, arm, "tfmean_db_pooled", "base_pooled"), abs=1e-12)
        assert s["L_gene"] == pytest.approx(_mean_diff(rows, arm, "gene_tf_pooled", "base_pooled"), abs=1e-12)
        for v in ms.COMPONENTS:
            assert s["components"][v] == pytest.approx(_mean_diff(rows, arm, f"{v}_pooled", "base_pooled"), abs=1e-12), v
        assert s["ratios"]["D/L"][0] == pytest.approx(s["delta_target"] / L, abs=1e-9)
        assert s["ratios"]["L_art/L"][0] == pytest.approx(s["L_art"] / L, abs=1e-9)
        assert s["ratios"]["L_gene/L"][0] == pytest.approx(s["L_gene"] / L, abs=1e-9)


def test_permutation_rank_p_adds_one_and_counts_ties_against_the_real_prior():
    def rows(real, perm, n=2):
        return [{"means": {"auto|c|real_pooled": real}, "seed_means": {"auto|c|perm_db_pooled": perm}}] * n
    assert ms.permutation_rank_p(rows(1.0, [0.5] * 20), "auto", "c") == pytest.approx((1 / 21, 2, 20))
    assert ms.permutation_rank_p(rows(1.0, [0.5] * 19 + [1.0]), "auto", "c")[0] == pytest.approx(2 / 21)
    assert ms.permutation_rank_p(rows(1.0, [2.0] * 20), "auto", "c")[0] == pytest.approx(1.0)
    # a dataset with a non-finite seed mean (or real value) leaves the rank and the returned count
    bad_seed = {"means": {"auto|c|real_pooled": 9.0}, "seed_means": {"auto|c|perm_db_pooled": [0.5] * 19 + [np.nan]}}
    bad_real = {"means": {"auto|c|real_pooled": np.nan}, "seed_means": {"auto|c|perm_db_pooled": [9.0] * 20}}
    assert ms.permutation_rank_p(rows(1.0, [0.5] * 19 + [1.0], 3) + [bad_seed, bad_real], "auto", "c") == \
        pytest.approx((2 / 21, 3, 20))
    assert np.isnan(ms.permutation_rank_p([bad_seed, bad_real], "auto", "c")[0])


@pytest.mark.parametrize("case", ["p_perm_1", "19_seeds"])
def test_t1_permutation_rank_failure_downgrades_the_verdict_in_the_report(planted, tmp_path, case):
    """The planted mock is target-level; only T1's robustness condition is made to fail: by a
    permuted prior above the real one in every seed, or by 19 seeds (1/20 cannot reach 0.05)."""
    out = tmp_path / "sens"
    shutil.copytree(planted["out"], out)
    key = "auto|pearson_motif|perm_db_pooled"
    for p in (out / "rows").glob("*.json"):
        row = json.loads(p.read_text())
        row["seed_means"][key] = [99.0] * len(row["seed_means"][key]) if case == "p_perm_1" else [-99.0] * 19
        p.write_text(json.dumps(row))
    text = ms.report(out)
    s = json.loads((out / "SENSITIVITY.json").read_text())["sets"][EXP]["arms"]["auto"]
    assert s["verdict"] == "mixed" and s["downgraded_by_t1_robustness"] is True
    assert [x[:2] for x in s["shortfalls"]] == ["T1"] and s["shortfalls"][0].startswith("T1 not robust")
    assert "downgraded from 'supports target-level evidence' to 'mixed'" in text
    if case == "p_perm_1":
        assert s["p_perm_T1"] == pytest.approx(1.0) and s["t1_robust_to_database"] is False
        assert "20 seeds): p = 1\n" in text and "(not assessable: too few seeds)" not in text
    else:                                                  # p_perm = 1/20: not assessable, not merely >= 0.05
        assert s["p_perm_T1"] == pytest.approx(1 / 20) and s["t1_robust_to_database"] is None
        assert "19 seeds): p = 0.05 (not assessable: too few seeds)" in text


def test_a_resubmission_with_another_n_perm_is_refused_and_rows_are_kept(planted, tmp_path, monkeypatch, capsys):
    """The planted folder was made with 20 seeds: a resume with --n-perm 8 would skip its completed rows
    and leave a set that mixes seed counts, so it is refused before anything is written."""
    out = tmp_path / "sens"
    shutil.copytree(planted["out"], out)
    (commit,) = {json.loads(p.read_text())["code_commit"] for p in (out / "rows").glob("*.json")}
    monkeypatch.setattr(ms, "code_commit", lambda: commit)
    monkeypatch.setenv(ct.ENV_DIR, str(planted["cis"]))
    argv = ["run", "--tasks", str(out / "sens_tasks.json"), "--out", str(out), "--refs", str(planted["refs"]),
            "--workroot", str(tmp_path / "work")]
    before = _bytes(out)
    capsys.readouterr()
    assert ms.main(argv + ["--n-perm", "8"]) == 2
    err = capsys.readouterr().err
    assert f"and {N_TASKS - 1} more row(s): existing row has n_perm=20, this run asks 8" in err
    assert "use a fresh --out or rerun them with one --n-perm (registered: 100)" in err and "from commit" not in err
    assert _bytes(out) == before
    assert ms.main(argv + ["--n-perm", "20"]) == 0                    # the same seed count: every unit skipped
    assert {k: v for k, v in _bytes(out).items() if k != "RUN_INPUTS.json"} == \
        {k: v for k, v in before.items() if k != "RUN_INPUTS.json"}


def test_rows_with_mixed_seed_counts_stop_the_report_with_a_clear_message(planted, tmp_path, capsys):
    out = tmp_path / "sens"
    shutil.copytree(planted["out"], out)
    p = out / "rows" / f"{EXP}__cxg__mock__homo-sapiens__t1.json"
    row = json.loads(p.read_text())
    key = "auto|pearson_motif|perm_db_pooled"
    row["seed_means"][key] = row["seed_means"][key][:8]               # as a resume with --n-perm 8 would leave it
    p.write_text(json.dumps(row))
    capsys.readouterr()
    assert ms.main(["report", "--out", str(out)]) == 2
    err = capsys.readouterr().err
    assert f"auto|pearson_motif: rows mix permutation seed counts [8, 20] (8 seeds: 1 row(s), e.g. {row['id']}; " \
           f"20 seeds: {N_TASKS - 1} row(s)" in err
    assert "rerun them with one --n-perm" in err and "inhomogeneous" not in err
    with pytest.raises(ValueError, match=r"mix permutation seed counts \[19, 20\]"):
        ms.permutation_rank_p([{"means": {"auto|c|real_pooled": 1.0}, "seed_means": {"auto|c|perm_db_pooled": s}}
                               for s in ([0.5] * 20, [0.5] * 19)], "auto", "c")


def test_a_seed_count_other_than_the_registered_100_is_reported_as_a_deviation(planted, tmp_path, monkeypatch):
    line = "DEVIATION: permutation seeds 20, registered 100 (p_perm resolution 1/21)"
    assert f"\n- {line}\n" in (planted["out"] / "SENSITIVITY.md").read_text()
    assert planted["summary"]["deviations"] == [line]
    out = tmp_path / "sens"
    shutil.copytree(planted["out"], out)
    monkeypatch.setattr(ms, "N_PERM", 20)                              # as if 20 were the registered count
    assert "DEVIATION" not in ms.report(out)
    assert json.loads((out / "SENSITIVITY.json").read_text())["deviations"] == []


def test_ratio_interval():
    def rows(real, base, num):
        return [{"task": f"t{i}", "means": {"auto|c|real_pooled": r, "auto|c|base_pooled": b, "auto|c|tfmean_db_pooled": x}}
                for i, (r, b, x) in enumerate(zip(real, base, num))]
    real, base, num = [1.5, 1.4, 1.7, 1.2, 1.6], [1.0, 1.1, 1.2, 1.0, 1.1], [1.2, 1.1, 1.5, 1.1, 1.3]
    point, lo, hi = ms.ratio_interval(rows(real, base, num), "auto", "c", ("tfmean_db", "base"))
    want = np.mean(np.subtract(num, base)) / np.mean(np.subtract(real, base))   # mean(num) / mean(L), not mean(num / L)
    assert point == pytest.approx(want, abs=1e-12) and lo <= point <= hi and lo < hi
    assert point != pytest.approx(np.mean(np.subtract(num, base) / np.subtract(real, base)), abs=1e-6)
    # fewer than 2 datasets with a value (NaN rows are dropped), or mean L == 0: no ratio
    one = rows(real[:2], base[:2], [num[0], np.nan])
    assert all(np.isnan(v) for v in ms.ratio_interval(one, "auto", "c", ("tfmean_db", "base")))
    flat = rows([1.5, 0.5], [1.0, 1.0], [1.2, 1.1])                    # L = +0.5 and -0.5
    assert all(np.isnan(v) for v in ms.ratio_interval(flat, "auto", "c", ("tfmean_db", "base")))


def test_verdict_cuts_at_exactly_half_of_the_lift():
    mk = lambda name, a, b: {"test": name, "res": {}, "p_a_gt_b_holm": a, "p_b_gt_a_holm": b}  # noqa: E731
    sup = [mk("T1", 0.01, 1.0), mk("T2", 0.01, 1.0), mk("T3", 0.01, 1.0)]
    none = [mk("T1", 0.5, 1.0), mk("T2", 0.5, 1.0), mk("T3", 0.5, 1.0)]
    v = lambda t, L, d, art, gene=0.02: ms.stage1_verdict(t, L, d, art, gene, t1_robust=True)  # noqa: E731
    assert v(sup, 0.2, 0.1, 0.05) == "target-level"          # D == 0.5 L passes ("D >= 0.5 L")
    assert v(sup, 0.2, 0.15, 0.1) == "mixed"                 # L_art == 0.5 L fails ("L_art < 0.5 L")
    assert v(sup, 0.2, 0.15, 0.05, gene=0.1) == "mixed"      # L_gene == 0.5 L fails
    assert v(none, 0.2, 0.0, 0.1) == "design effects only"   # rule 5: L_art >= 0.5 L


# (e') the verdict line says why: rule 4, tests not run, the one shared database --------------- #

NOISE = (0.010, -0.012, 0.004, -0.006, 0.004)            # mean 0: a difference made of it is unresolved


def _fabricated(planted, tmp_path, t1=0.0, t2=0.0, t3=0.0, art=0.75, robust=True) -> Path:
    """Copy of the planted report folder with fabricated decision means for the stage-1 pair, both
    arms alike: L = 0.2, L_art = art * L and L_gene = 0.1 L in every dataset; the T1, T2 and T3
    differences are the given mean plus NOISE (0: unresolved, +0.1: supported, -0.1: reversed).
    Gates 1 and 2 stay the planted mock's (they read other fields). robust=False puts every
    permuted seed above the real prior (p_perm = 1)."""
    out = tmp_path / "sens"
    shutil.copytree(planted["out"], out)
    for i, p in enumerate(sorted((out / "rows").glob("*.json"))):
        row = json.loads(p.read_text())
        e, base, base_r = NOISE[i], 1.3 + 0.05 * i, 1.2 + 0.03 * i
        vals = {"base_pooled": base, "real_pooled": base + 0.2, "perm_db_pooled": base + 0.2 - (t1 + e),
                "tfmean_db_pooled": base + art * 0.2, "gene_tf_pooled": base + 0.02,
                "base_macro_R": base_r, "real_macro_R": base_r + t2 + e, "gene_tf_macro_R": base_r + t2 - t3}
        for arm in ms.ARMS:
            row["means"].update({f"{arm}|pearson_motif|{k}": v for k, v in vals.items()})
            key = f"{arm}|pearson_motif|perm_db_pooled"
            row["seed_means"][key] = [vals["perm_db_pooled"] if robust else 99.0] * len(row["seed_means"][key])
        p.write_text(json.dumps(row))
    return out


def _verdict_line(text, arm="auto"):
    return next(x for x in text.splitlines() if x.startswith(f"- **verdict ({arm}): "))


@pytest.mark.parametrize("supported, missing, robust", [("T1", "T2", True), ("T2", "T1", True), ("T1", "T2", False)])
def test_a_rule_4_mixed_names_its_reasons_on_the_verdict_line(planted, tmp_path, supported, missing, robust):
    """Exactly one of T1 and T2 supported and L_art = 0.75 L: record rule 4 says mixed. The line names
    the test that is not supported and the L_art cut; the rule-3 list and its TF-generic sentence stay
    out, and T1's robustness is reported only."""
    out = _fabricated(planted, tmp_path, robust=robust, **{supported.lower(): 0.1})
    text = ms.report(out)
    s = json.loads((out / "SENSITIVITY.json").read_text())["sets"][EXP]["arms"]["auto"]
    assert [_test(s, n)["p_a_gt_b_holm"] < 0.05 for n in ("T1", "T2")] == [supported == "T1", supported == "T2"]
    assert s["verdict"] == "mixed" and s["shortfalls"] == [] and s["downgraded_by_t1_robustness"] is False
    assert s["rule4_reasons"] == [f"{missing} not supported",
                                  "L_art >= 0.5 L: the design-effect control reaches at least half of the lift"] + \
        ([] if robust else [ms.T1_NOT_ROBUST])
    assert s["untested"] == [] and s["n_analysed"] == N_TASKS
    line = _verdict_line(text)
    assert line.startswith("- **verdict (auto): mixed**. T1 and/or T2 is supported by the paired t test")
    assert (f" Record rule 4: only one of T1 and T2 is supported ({missing} is not) and the design-effect "
            "control reaches at least half of the lift (L_art/L 0.75 >= 0.5).") in line
    assert ("T1 not robust to the realization of the one shared database (permutation rank p >= 0.05 or not "
            "assessable): reported only, it does not change this verdict (record).") in line or robust
    assert "not robust" not in line or not robust
    assert "Conditions of the strongest verdict" not in line and "TF-generic gene propensity" not in line
    assert "Evidence that target identity matters exists" not in text and "named below" not in text
    assert "Record rule 4" in _verdict_line(text, "baseline")


def test_tests_that_were_not_run_are_called_not_run_and_rule_5_is_unchanged(null, tmp_path):
    """Two datasets through gate 1: no declared test runs (fewer than 3). The verdict follows the record
    (rule 5 on the point estimates), and the line says the tests were not run, not that they failed."""
    out = tmp_path / "sens"
    shutil.copytree(null["out"], out)
    for p in sorted((out / "rows").glob("*.json"))[2:]:
        row = json.loads(p.read_text())
        row["flags"] = ["GATE score pearson_signed/auto: |recomputed - committed| = 0.1"]
        p.write_text(json.dumps(row))
    text = ms.report(out)
    st = json.loads((out / "SENSITIVITY.json").read_text())["sets"][EXP]
    s = st["arms"]["auto"]
    assert st["gate_pass"] == 2 and st["not_interpretable"] is False
    assert s["verdict"] == "design effects only"                    # rule 5 on the point estimates, unchanged
    assert s["untested"] == ["T1", "T2", "T3"] and s["n_analysed"] == 2 and s["tests"][0]["holm_family"] == 0
    line = _verdict_line(text)
    assert ("T1, T2, T3 not tested (fewer than 3 analysed datasets with a value; the family shrinks, record): "
            "'not supported' here means not tested, and the verdict rests on the point estimates L and L_art/L "
            "only.") in line
    assert "(c) neither T1 nor T2 was run and" in line and "neither T1 nor T2 is supported" not in line
    # a test run in neither arm has no sign: not compared, never a disagreement, never "all agree"
    assert "DISAGREE" not in line and s["secondary_disagreements"] == []
    assert [x.split(" mean")[0] for x in s["secondary_not_compared"]] == [f"baseline arm T{i}" for i in (1, 2, 3)]
    assert ("Secondary checks (sign agreement; they never change the verdict): no compared pair disagrees. Not "
            "compared (no finite value in one or both arms): baseline arm T1 mean difference +nan vs +nan;") in line
    full = (null["out"] / "SENSITIVITY.md").read_text()            # every dataset analysed: nothing added
    assert "not tested (fewer than 3" not in full and "(c) neither T1 nor T2 is supported and" in full
    assert _arm(null)["untested"] == [] and _arm(null)["n_analysed"] == N_TASKS


@pytest.mark.parametrize("art, without", [(0.1, "inconclusive"), (0.75, "mixed")])
def test_rule_2_through_t2_reverse_alone_is_named_and_the_verdict_stands(planted, tmp_path, art, without):
    """T1 and T3 supported, T2 reversed: registered rule 2 says 'design effects only' and that stands.
    The line adds that it fired through T2's reverse alone and what the rules give without it."""
    out = _fabricated(planted, tmp_path, t1=0.1, t2=-0.1, t3=0.1, art=art)
    text = ms.report(out)
    s = json.loads((out / "SENSITIVITY.json").read_text())["sets"][EXP]["arms"]["auto"]
    assert [(_test(s, n)["p_a_gt_b_holm"] < 0.05, _test(s, n)["p_b_gt_a_holm"] < 0.05) for n in ("T1", "T2", "T3")] \
        == [(True, False), (False, True), (True, False)]
    assert s["verdict"] == "design effects only"                    # rule 2, as registered
    assert s["rule2_via_t2_reverse_only"] is True and s["verdict_without_t2_reverse"] == without
    line = _verdict_line(text)
    assert line.startswith("- **verdict (auto): consistent with design effects only**. One of three routes (record "
                           "rules 2 and 5): (a) T1 reversed: ")
    assert line.count("C4/C6 may say that the motif prior raises the pooled AUPRC ratio, not that sequence evidence "
                      "at the target does; a confirmatory test with macro_R as primary metric needs a D-xxx.") == 1
    assert (" Rule 2 fired through T2's reverse only: the rank-product combination lowers macro_R below the base "
            "within TF. T1 and T3 are supported. Without T2's reverse the registered rules would give "
            f"'{ms.VERDICTS[without][0]}'. The registered verdict and its consequence stand; quote this together "
            "with D/L, L_art/L and L_gene/L.") in line
    assert "carries target identity" not in line         # facts only: no reading that softens rule 2
    assert ("The design-effect control reaches less than half of the lift." in line) == (art < 0.5)
    assert _arm(planted)["rule2_via_t2_reverse_only"] is False and _arm(planted)["verdict_without_t2_reverse"] is None


def test_secondary_checks_name_every_sign_disagreement_and_never_call_a_missing_value_one():
    def arm(L, epr, *diffs):
        return {"L": L, "epr_lift": epr, "tests": [{"test": n, "res": None if d is None else {"difference": d}}
                                                   for n, d in zip(("T1", "T2", "T3"), diffs)]}
    prim = arm(0.2, 0.1, 0.1, 0.1, 0.1)
    assert ms.secondary_checks(prim, None) == ([], [])
    assert ms.secondary_checks(prim, arm(0.05, 0.0, 0.3, 0.2, 0.01)) == ([], [])     # same signs, any size
    assert ms.secondary_text([], []) == "all agree."
    # a real sign flip is a disagreement; a test run in one arm only is not compared
    dis, nc = ms.secondary_checks(prim, arm(-0.1, 0.0, 0.1, -0.1, None))
    assert dis == ["baseline arm L -0.1000 vs +0.2000", "baseline arm T2 mean difference -0.1000 vs +0.1000"]
    assert nc == ["baseline arm T3 mean difference +nan vs +0.1000"]
    assert ms.secondary_text(dis, nc) == (
        "DISAGREE: baseline arm L -0.1000 vs +0.2000; baseline arm T2 mean difference -0.1000 vs +0.1000. "
        "Not compared (no finite value in one or both arms): baseline arm T3 mean difference +nan vs +0.1000.")
    # EPR against L: a sign flip disagrees, a missing EPR is not compared
    assert ms.secondary_checks(arm(0.2, -0.05, 0.1, 0.1, 0.1), None) == \
        (["EPR(real) - EPR(base) -0.0500 vs L +0.2000"], [])
    assert ms.secondary_checks(arm(0.2, float("nan"), 0.1, 0.1, 0.1), None) == \
        ([], ["EPR(real) - EPR(base) +nan vs L +0.2000"])
    # neither arm ran the tests (fewer than 3 datasets): three pairs not compared, and the text says so
    dis, nc = ms.secondary_checks(arm(0.2, 0.1, None, None, None), arm(0.15, 0.0, None, None, None))
    assert dis == [] and len(nc) == 3 and all(x.endswith("+nan vs +nan") for x in nc)
    assert ms.secondary_text(dis, nc).startswith("no compared pair disagrees. Not compared (no finite value in one "
                                                 "or both arms): baseline arm T1 mean difference +nan vs +nan; ")
    # no dataset at all: EPR and L are not compared either
    empty = arm(float("nan"), float("nan"), None, None, None)
    assert ms.secondary_checks(empty, empty)[0] == [] and len(ms.secondary_checks(empty, empty)[1]) == 5


def test_every_table_of_declared_tests_says_its_p_values_are_conditional_on_the_one_database(planted):
    lines = (planted["out"] / "SENSITIVITY.md").read_text().splitlines()
    heads = [i for i, x in enumerate(lines) if x.startswith("Paired t over datasets")]
    assert len(heads) == len(ms.STAGE1_PAIRS) * len(ms.ARMS)
    sentence = ("All datasets of this set share one cisTarget database and one reference: these p-values measure "
                "consistency over datasets given that database, not inference over database realizations. Only "
                "the T1 permutation rank below accounts for it; T2 and T3 have no such check in this record.")
    assert all(lines[i + 1] == sentence and lines[i + 2] == "" and lines[i + 3].startswith("| test |")
               for i in heads)
    mixed = ms.VERDICTS["mixed"][1]
    assert "not evidence that target identity matters" in mixed and "named on this line" in mixed
    assert "Evidence that target identity matters exists" not in mixed and "named below" not in mixed


# (f) every registered code check and gate-1 condition fires ----------------------------------- #

def _first_unit(planted):
    return ms.plan_units([f"{EXP}={planted['run']}"], planted["results"])[0]


@pytest.mark.parametrize("what, match", [("builtin", "built-in checks"), ("selfcheck", "self-consistency"),
                                         ("fast", "fast scorer"), ("internal", "internal consistency")])
def test_each_code_check_stops_the_unit_with_check_failed(planted, monkeypatch, what, match):
    """One fault per registered code check: run_unit must raise CheckFailed (status check_failed)."""
    if what == "builtin":
        monkeypatch.setattr(ms, "builtin_check_violations", lambda rec, applies=True: ["tf_mean_macro 1 != base_macro 2"])
    elif what == "selfcheck":
        orig = ms.rank_product
        monkeypatch.setattr(ms, "rank_product", lambda E, m: orig(E, -np.asarray(m, float)))
    elif what == "fast":
        orig_score = ms.Scorer.score

        def off(self, w):
            r = orig_score(self, w)
            r["pooled"] += 1e-6
            return r
        monkeypatch.setattr(ms.Scorer, "score", off)
    else:
        orig_mean = ms._finite_mean
        monkeypatch.setattr(ms, "_finite_mean", lambda values: orig_mean(values) + 1e-6)
    with pytest.raises(ms.CheckFailed, match=match):
        _run_one(planted, _first_unit(planted), monkeypatch)


def test_an_exception_inside_score_models_is_a_runtime_error_and_is_rerun(planted, tmp_path, monkeypatch):
    """score_models records a model's exception (m_<model>_<arm>_error) instead of raising. A
    transient fault there (I/O, memory) is a runtime error, exclusion rule (c): status error with the
    original message, rerun by the resubmission. Never check_failed (rule 0, never rerun)."""
    import grnbench.infer as infer
    orig, budget = infer.run_method, {"faults": 1}

    def flaky(*a, **k):
        if budget["faults"] and sys._getframe(1).f_code.co_name == "score_models":
            budget["faults"] -= 1
            raise MemoryError("simulated transient out-of-memory")
        return orig(*a, **k)
    monkeypatch.setattr(infer, "run_method", flaky)
    monkeypatch.setattr(real_bench, "ARTIFACT_TOP_EDGES", CAP)
    monkeypatch.setenv(ct.ENV_DIR, str(planted["cis"]))
    unit, tasks, out = _first_unit(planted), tmp_path / "tasks.json", tmp_path / "sens"
    tasks.write_text(json.dumps([unit]))
    argv = ["run", "--tasks", str(tasks), "--out", str(out), "--refs", str(planted["refs"]), "--n-perm", "2",
            "--workroot", str(planted["root"] / "work")]
    path = out / "rows" / f"{unit['id']}.json"
    assert ms.main(argv) == 2 and budget["faults"] == 0
    row = json.loads(path.read_text())
    assert row["status"] == "error" and len(row["flags"]) == 1
    assert row["flags"][0].startswith("ERROR ScoreModelsError: score_models pearson_signed_auto: "
                                      "MemoryError: simulated transient out-of-memory")
    assert ms.main(argv) == 0                                           # the resubmission reruns it
    row = json.loads(path.read_text())
    assert row["status"] == "completed" and row["flags"] == []


def test_a_finite_disagreement_with_score_models_is_still_a_failed_code_check(planted, monkeypatch):
    """Only an exception recorded by score_models is a runtime error: a finite pooled value that
    differs from the script's per-context mean by 1e-6 is the internal-consistency code check."""
    orig = ms.registered_rescore

    def off(*a, **k):
        reg, edges = orig(*a, **k)
        assert not any(key.endswith("_error") for key in reg)
        reg["m_pearson_motif_auto_auprc_ratio"] += 1e-6
        return reg, edges
    monkeypatch.setattr(ms, "registered_rescore", off)
    with pytest.raises(ms.CheckFailed, match="internal consistency: auto/real pooled"):
        _run_one(planted, _first_unit(planted), monkeypatch)


def _no_gold_copy(planted, src: Path) -> Path:
    """The kept auto context of t0 with every non-TF symbol renamed (G* -> Z*): no evaluation pair
    (all TF -> G*) lies in its universe."""
    a = ad.read_h5ad(planted["run"] / "artifacts" / "cxg__mock__homo-sapiens__t0" / "auto" / "ctx1" / "processed.h5ad")
    a.var_names = [("Z" + g[1:]) if g.startswith("G") else g for g in map(str, a.var_names)]
    a.write_h5ad(src / "nogold.h5ad")
    return src / "nogold.h5ad"


def _unit_with_auto_contexts(planted, root: Path, auto: list[Path], monkeypatch) -> tuple[dict, dict]:
    """t0 of the planted mock re-kept with the given auto contexts (ctx1, ctx2, ...): the registered
    score_models writes the kept files and the committed row, as in build_mock."""
    task, tdir0 = "cxg__mock__homo-sapiens__t0", planted["run"] / "artifacts" / "cxg__mock__homo-sapiens__t0"
    run, results, src = root / "run", root / "results", root / "src"
    ta = run / "artifacts" / task
    for p in (ta, run / "rows", results / EXP / "rows", src):
        p.mkdir(parents=True)
    for name in ("reference_selection.tsv", "reference_evaluation.tsv", "decisions.json"):
        shutil.copy(tdir0 / name, ta / name)
    shutil.copy(planted["results"] / EXP / "REPORT_REAL.md", results / EXP / "REPORT_REAL.md")
    shutil.copy(tdir0 / "baseline" / "ctx0" / "processed.h5ad", src / "baseline.h5ad")
    sel, ev = ms.read_pairs(ta / "reference_selection.tsv"), ms.read_pairs(ta / "reference_evaluation.tsv")
    ev_pairs = set(zip(ev["source"], ev["target"]))
    recs = {"baseline": [{"processed": str(src / "baseline.h5ad"), "context": "all"}],
            "candX": [{"processed": None, "context": "skipped"}] + [{"processed": str(p), "context": p.stem} for p in auto]}
    monkeypatch.setattr(real_bench, "ARTIFACT_TOP_EDGES", CAP)
    monkeypatch.setenv(ct.ENV_DIR, str(planted["cis"]))
    row = {"task": task, "truth": json.dumps({"organism": "human"}), "reference_edges": len(sel) + len(ev)}
    row.update(real_bench.score_models(recs, {"baseline": _Plan("cpm_log1p"), "candX": _Plan("raw_counts")},
                                       {"baseline": "baseline", "auto": "candX"},
                                       ["pearson_signed", "pearson_motif", "motif_prior"],
                                       sorted(set(sel["source"]) | set(ev["source"])), ev_pairs, artifacts_dir=ta))
    assert not any(k.endswith("_error") for k in row) and {t[:1] for _, t in ev_pairs} == {"G"}
    for d in (run / "rows", results / EXP / "rows"):
        (d / f"{task}.json").write_text(json.dumps(row, default=float), encoding="utf-8")
    (unit,) = ms.plan_units([f"{EXP}={run}"], results)
    assert unit["plan_flags"] == []
    return unit, row


def test_a_context_without_evaluation_gold_is_nan_and_skipped(planted, tmp_path, monkeypatch):
    """A context whose universe holds no evaluation pair has an undefined (NaN) AUPRC ratio, in the
    fast scorer and in evaluate_ranking alike: agreement, not a failed code check. The dataset value
    is the mean over the contexts with a finite value, as score_models (record section 1); the
    context's macro metrics are NaN and leave the tests that need them (exclusion rule d)."""
    from grnbench.eval.metrics import GoldIndex
    from grnbench.infer import prepare_expression, run_method
    from grnbench.tfs import resolve_tf_names
    monkeypatch.setenv(ct.ENV_DIR, str(planted["cis"]))
    # decompose alone: gold made of pairs outside the data
    tdir = planted["run"] / "artifacts" / "cxg__mock__homo-sapiens__t0"
    adata = ad.read_h5ad(tdir / "baseline" / "ctx0" / "processed.h5ad")
    genes = pd.Index(adata.var_names).astype(str)
    tfs = [f"TF{i}" for i in range(N_TF)]
    E, _ = run_method("pearson_signed", adata, tf_list=tfs, assume_normalized=True)
    expr, _ = prepare_expression(adata, assume_normalized=True)
    universe = pd.DataFrame([(s, g) for s in resolve_tf_names(genes, tfs) for g in genes if s != g],
                            columns=["source", "target"])
    labelable = set(pd.read_csv(planted["refs"] / "collectri_human.tsv", sep="\t")["target"].str.upper())
    gi = GoldIndex.coerce({("TF0", "NOWHERE0"), ("TF5", "NOWHERE1")})
    rec, _ = ms.decompose(E, expr, universe, gi, labelable, n_perm=2)
    assert rec["gold_in_universe"] == 0 and rec["n_tfs_macro"] == rec["n_tfs_macro_R"] == 0
    assert all(np.isnan(rec[f"{v}_{m}"]) for v in ("base",) + ms.VARIANTS for m in ms.METRICS)
    assert rec["builtin_checks"].startswith("held")
    # run_unit: auto arm = ctx1 (gold) + ctx2 (no gold), kept and committed by score_models itself
    src = tmp_path / "nogold"
    src.mkdir()
    ctx1 = tdir / "auto" / "ctx1" / "processed.h5ad"
    unit, committed = _unit_with_auto_contexts(planted, tmp_path / "m", [ctx1, _no_gold_copy(planted, src)], monkeypatch)
    row = _run_one(planted, unit, monkeypatch)
    assert row["status"] == "completed" and row["flags"] == [] and max(row["gate"].values()) == 0.0
    ctx = {(c["arm"], c["ctx"]): c for c in row["ctx"]}
    assert sorted(ctx) == [("auto", "ctx1"), ("auto", "ctx2"), ("baseline", "ctx0")]
    c1, c2 = ctx[("auto", "ctx1")], ctx[("auto", "ctx2")]
    assert c2["gold_in_universe"] == 0 < c1["gold_in_universe"] and c2["n_tfs_macro"] == 0
    assert all(np.isnan(c2[f"{v}_{m}"]) for v in ("base", "real") for m in ("pooled", "macro", "macro_R"))
    means = row["means"]
    assert means["auto|pearson_motif|real_pooled"] == c1["real_pooled"] == committed["m_pearson_motif_auto_auprc_ratio"]
    assert means["auto|pearson_motif|base_pooled"] == c1["base_pooled"] == committed["m_pearson_signed_auto_auprc_ratio"]


def test_an_arm_without_evaluation_gold_fails_gate_1_not_a_code_check(planted, tmp_path, monkeypatch):
    """Every auto context without gold: the dataset value is NaN in the script and in score_models
    (agreement, no CheckFailed), and gate 1 stays strict: a NaN never equals the committed value, so
    the dataset is flagged and excluded, never silently passed."""
    src = tmp_path / "nogold"
    src.mkdir()
    unit, committed = _unit_with_auto_contexts(planted, tmp_path / "m", [_no_gold_copy(planted, src)], monkeypatch)
    assert np.isnan(committed["m_pearson_motif_auto_auprc_ratio"])
    row = _run_one(planted, unit, monkeypatch)
    assert row["status"] == "completed" and np.isnan(row["means"]["auto|pearson_motif|real_pooled"])
    assert sorted(f.split(":")[0] for f in row["flags"]) == sorted(f"GATE score {m}/auto" for m in ms.RECOMPUTABLE)


def test_gate_fails_loudly_when_the_committed_value_is_missing(planted, tmp_path, monkeypatch):
    results = tmp_path / "results"
    shutil.copytree(planted["results"], results)
    p = results / EXP / "rows" / "cxg__mock__homo-sapiens__t0.json"
    row = json.loads(p.read_text())
    del row["m_pearson_motif_baseline_auprc_ratio"]          # NaN difference: must flag, never pass
    p.write_text(json.dumps(row))
    unit = next(u for u in ms.plan_units([f"{EXP}={planted['run']}"], results)
                if u["task"] == "cxg__mock__homo-sapiens__t0")
    flags = _run_one(planted, unit, monkeypatch)["flags"]
    assert any(f.startswith("GATE score pearson_motif/baseline") for f in flags)


def test_gate_fails_when_kept_edge_targets_are_swapped(planted, tmp_path, monkeypatch):
    def alter(tdir):
        path = tdir / "baseline" / "ctx0" / "pearson_signed.tsv.gz"
        e = ms.read_kept(path)
        i = next(k for k in range(10) if e.loc[k, "target"] != e.loc[k + 1, "target"])
        e.loc[i, "target"], e.loc[i + 1, "target"] = e.loc[i + 1, "target"], e.loc[i, "target"]   # same weights
        e.to_csv(path, sep="\t", index=False, compression="gzip")
    row = _run_one(planted, _altered_unit(planted, tmp_path, alter), monkeypatch)
    assert any(f.startswith("GATE edges baseline/ctx0/pearson_signed") for f in row["flags"])


def test_worst_tie_check_is_skipped_when_an_in_database_gene_has_the_worst_rank(planted, monkeypatch):
    orig = ms.motif_parts

    def parts(tfs_e, genes):                       # lower the worst rank to the median: many in-database
        M, in_db, K, worst = orig(tfs_e, genes)    # genes of a TF with a motif now sit at or past W
        return M, in_db, K, float(np.median(M[np.ix_(K > 0, in_db)]))
    monkeypatch.setattr(ms, "motif_parts", parts)
    row = _run_one(planted, _first_unit(planted), monkeypatch)
    assert row["ctx"] and all("worst_tie on U'' not applicable" in c["builtin_checks"] for c in row["ctx"])


def test_report_counts_the_contexts_where_worst_tie_was_not_applicable(planted, tmp_path):
    out = tmp_path / "sens"
    shutil.copytree(planted["out"], out)
    p = out / "rows" / f"{EXP}__cxg__mock__homo-sapiens__t1.json"
    row = json.loads(p.read_text())
    row["ctx"][0]["builtin_checks"] = ("held (worst_tie on U'' not applicable: an in-database gene has the "
                                       "worst rank for a TF with a motif)")
    p.write_text(json.dumps(row))
    text = ms.report(out)
    assert f"held on {2 * N_TASKS}/{2 * N_TASKS} contexts; worst_tie on U'' not applicable on 1\n" in text
    assert "not applicable" not in (planted["out"] / "SENSITIVITY.md").read_text()


# (g) the keyed permutation, the gene propensity, U', U'' and the macro threshold, by definition --- #

def test_keyed_permutation_gives_every_tf_its_own_order_and_ignores_the_case_of_names():
    genes = [f"G{i}" for i in range(80)]
    tfs = ["TFA", "TFB", "TFC"]
    M = np.tile(np.arange(1.0, 81.0), (3, 1))                            # three TFs, the same value profile
    in_db, K = np.ones(80, bool), np.array([1, 1, 1])

    def perm(tf_names, gene_names, rows=slice(None)):
        keys = (ms.name_keys(tf_names), ms.name_keys(gene_names))
        return ms.variant_matrix("perm_db", M[rows], in_db, K[rows], 1000, 0, keys)
    P = perm(tfs, genes)
    assert not np.array_equal(P[0], P[1]) and not np.array_equal(P[1], P[2]) and not np.array_equal(P[0], P[2])
    assert np.array_equal(perm(["TFB"], genes, slice(1, 2))[0], P[1])    # a TF's row does not depend on the others
    assert np.array_equal(perm([t.lower() for t in tfs], [g.lower() for g in genes]), P)   # case-insensitive keys


def test_gene_tf_is_the_leave_one_out_propensity_of_the_other_tfs_with_motifs():
    rng = np.random.default_rng(5)
    n_tf, n_g, W = 9, 30, 500.0
    in_db = rng.random(n_g) > 0.25
    K = np.array([0, 3, 1, 2, 0, 5, 1, 4, 2])
    M = np.full((n_tf, n_g), W)
    for i in np.flatnonzero(K):
        M[i, in_db] = rng.permutation(np.arange(1.0, in_db.sum() + 1))
    G = ms.variant_matrix("gene_tf", M, in_db, K, W)
    motif_tfs, cols = np.flatnonzero(K > 0), np.flatnonzero(in_db)
    assert (G[K == 0] == W).all() and (G[:, ~in_db] == W).all()
    differs = False
    for i in motif_tfs:
        prop = M[np.ix_([j for j in motif_tfs if j != i], cols)].mean(0)
        expect = np.full(n_g, W)
        expect[cols[np.argsort(prop, kind="stable")]] = np.sort(M[i, cols])
        assert np.array_equal(G[i], expect), i
        full = M[np.ix_(motif_tfs, cols)].mean(0)
        differs |= not np.array_equal(np.argsort(full, kind="stable"), np.argsort(prop, kind="stable"))
    assert differs                                  # the matrix can tell leave-one-out from the full mean


def test_macro_needs_at_least_five_positives_in_the_scored_row():
    from grnbench.eval.metrics import GoldIndex
    targets = [f"G{i}" for i in range(30)]
    n_pos = {"A": 5, "B": 4, "C": 8}
    uni = pd.DataFrame([(t, g) for t in n_pos for g in targets], columns=["source", "target"])
    gi = GoldIndex.coerce({(t, g) for t, k in n_pos.items() for g in targets[:k]})
    edges = uni.assign(weight=np.random.default_rng(1).random(len(uni)))
    sc = ms.Scorer(uni, gi, {}, edges)
    w = sc.universe_weights(edges["weight"].to_numpy())
    ap = {t: ms.Scorer._ap_ratio(sc.y[b], w[b]) for t, b in zip(sc.tf_names, sc.blocks)}
    assert ms.MIN_POS_TF == 5
    assert sc.macro(w) == pytest.approx((ap["A"] + ap["C"]) / 2)            # B (4) out, A (5) in
    mask = (~(((uni["source"] == "A") & (uni["target"] == "G0"))             # A: 5 -> 4 positives
              | ((uni["source"] == "C") & uni["target"].isin(["G0", "G1"])))).to_numpy()  # C: 8 -> 6
    bc = sc.blocks[sc.tf_names.index("C")]
    bc = bc[mask[bc]]
    assert sc.macro(w, mask) == pytest.approx(ms.Scorer._ap_ratio(sc.y[bc], w[bc]))   # counted in the masked row
    assert [len(b) for b in sc.macro_blocks()] == [30, 30] and [len(b) for b in sc.macro_blocks(mask)] == [28]


def test_scorer_matches_evaluate_ranking():
    """The fast scorer against the registered evaluate_ranking, per TF row (macro) and per mask (T, R),
    on random universes with heavy weight ties and partial edge lists. TF0 has exactly MIN_POS_TF
    evaluation positives in its row, TF1 one fewer."""
    from grnbench.eval.metrics import GoldIndex, evaluate_ranking

    def ref_macro(E, gold, U):
        vals = []
        for tf, sub in U.groupby("source", sort=False):
            if sum((tf, t) in gold for t in sub["target"]) < 5:
                continue
            e = E.merge(sub[["source", "target"]], on=["source", "target"])
            vals.append(evaluate_ranking(e, gold, with_curve=False, universe=sub)["auprc_ratio"])
        return float(np.mean(vals)) if vals else np.nan
    n = 0
    for seed in range(5):
        rng = np.random.default_rng(seed)
        tfs, genes = [f"TF{i}" for i in range(6)], [f"G{j}" for j in range(40)]
        U = pd.DataFrame([(s, t) for s in tfs for t in genes], columns=["source", "target"])
        n_pos = [5, 4] + [int(k) for k in rng.integers(6, 14, len(tfs) - 2)]
        gold = {(s, str(t)) for s, k in zip(tfs, n_pos) for t in rng.choice(genes, k, replace=False)}
        gi = GoldIndex.coerce(gold)
        E = U.sample(frac=0.7, random_state=seed).reset_index(drop=True)
        E = E.assign(weight=rng.integers(0, 5, len(E)).astype(float))       # heavy ties
        tmask = U["target"].isin(rng.choice(genes, 30, replace=False)).to_numpy()
        rmask = tmask & U["source"].isin(tfs[:-1]).to_numpy()
        fast = ms.Scorer(U, gi, {"T": tmask, "R": rmask}, E).score(E["weight"].to_numpy())
        ref = {"pooled": evaluate_ranking(E, gi, with_curve=False, universe=U)["auprc_ratio"],
               "macro": ref_macro(E, gold, U)}
        for tag, m in (("T", tmask), ("R", rmask)):
            Um = U[m].reset_index(drop=True)
            Em = E.merge(Um, on=["source", "target"])
            ref[f"pooled_{tag}"] = evaluate_ranking(Em, gi, with_curve=False, universe=Um)["auprc_ratio"]
            ref[f"macro_{tag}"] = ref_macro(Em, gold, Um)
        assert set(fast) == set(ref) == set(ms.METRICS)
        for k in ms.METRICS:
            assert np.isfinite(ref[k]) and abs(fast[k] - ref[k]) <= ms.TOL, (seed, k, fast[k], ref[k])
            n += 1
    assert n == 5 * 6


@pytest.mark.parametrize("which", ["planted", "null"])
def test_restricted_universes_match_an_independent_count(which, request):
    """U' = target in a database AND labelable; U'' = U' x TFs with a motif; macro needs 5 positives
    in the row: recounted from the mock's feather, reference and processed files."""
    mock = request.getfixturevalue(which)
    db = set(pf.read_table(str(mock["cis"] / "hg38_a_fake.genes_vs_motifs.rankings.feather")).column_names) - {"motifs"}
    labelable = set(pd.read_csv(mock["refs"] / "collectri_human.tsv", sep="\t")["target"].str.upper())
    with_motif = {f"TF{i}" for i in range(N_NO_MOTIF, N_TF)}
    n = separated = 0
    for row in _rows(mock):
        tdir = mock["run"] / "artifacts" / row["task"]
        ev = pd.read_csv(tdir / "reference_evaluation.tsv", sep="\t")
        gold = set(zip(ev["source"], ev["target"]))
        for c in row["ctx"]:
            n += 1
            genes = list(ad.read_h5ad(tdir / c["arm"] / c["ctx"] / "processed.h5ad").var_names)
            pairs = [(s, g) for s in genes if s.startswith("TF") for g in genes if g != s]
            ok_t = [g in db and g.upper() in labelable for _, g in pairs]       # U'
            ok_r = [t and s in with_motif for t, (s, _) in zip(ok_t, pairs)]    # U''
            assert c["n_universe"] == len(pairs)
            assert c["frac_universe_T"] == pytest.approx(np.mean(ok_t))
            assert c["frac_universe_R"] == pytest.approx(np.mean(ok_r))
            assert c["gold_in_T"] == sum(p in gold for p, k in zip(pairs, ok_t) if k)
            assert c["gold_in_R"] == sum(p in gold for p, k in zip(pairs, ok_r) if k)
            assert 0 < c["frac_universe_R"] < c["frac_universe_T"] < 1
            pos = pd.Series([s for (s, g) in pairs if (s, g) in gold]).value_counts()
            pos_r = pd.Series([s for (s, g), k in zip(pairs, ok_r) if k and (s, g) in gold]).value_counts()
            assert c["n_tfs_macro"] == int((pos >= 5).sum()) and c["n_tfs_macro_R"] == int((pos_r >= 5).sum())
            separated += c["n_tfs_macro_R"] < c["n_tfs_macro"]
    assert n == 2 * N_TASKS and separated                # the counts can tell U'' rows from full rows


# (h) shapes the mocks above lack: several contexts, run-level exclusion and skip, the report
#     fields that decide nothing, degenerate contexts and a mouse dataset ------------------------ #

def test_dataset_values_average_every_context_and_every_permutation_seed(multi):
    """Record section 1: the dataset value is the mean over the dataset's contexts with a finite
    value (as score_models). Each permutation seed's dataset mean averages that seed over the
    contexts: those means feed T1's permutation rank, which no runtime check compares."""
    st = multi["summary"]["sets"][EXP]
    assert st["gate_pass"] == st["units"] == 3 and st["gate2"]["pass"]       # gates 1 and 2 with 2 contexts
    rows = _rows(multi)
    assert len(rows) == 3 and all(row["flags"] == [] and row["status"] == "completed" for row in rows)
    for row in rows:
        seeds = json.loads((multi["out"] / "perm_seeds" / f"{row['id']}.json").read_text())
        for arm in ms.ARMS:
            ctx = [c for c in row["ctx"] if c["arm"] == arm]
            assert len(ctx) == 2 and ctx[0]["real_pooled"] != ctx[1]["real_pooled"]
            for col in ("real_pooled", "base_pooled", "gene_tf_macro_R", "perm_db_pooled", "tfmean_db_pooled"):
                assert row["means"][f"{arm}|pearson_motif|{col}"] == \
                    pytest.approx(np.mean([c[col] for c in ctx]), abs=1e-12), col
            for k in ms.SEED_KEYS:
                per = [seeds[f"{arm}|{c['ctx']}|pearson_motif"][f"perm_db_{k}"] for c in ctx]
                assert len(per[0]) == 4 and not np.allclose(per[0], per[1])   # the contexts can be told apart
                assert row["seed_means"][f"{arm}|pearson_motif|perm_db_{k}"] == \
                    pytest.approx(np.mean(per, axis=0), abs=1e-12), k


def test_finite_mean_ignores_non_finite_contexts_as_score_models_does():
    assert ms._finite_mean([1.0, float("nan"), 3.0, None]) == 2.0
    assert ms._finite_mean([1.0, np.inf, 3.0]) == 2.0                   # score_models keeps np.isfinite values
    assert np.isnan(ms._finite_mean([float("nan"), None])) and np.isnan(ms._finite_mean([]))


def test_units_without_a_committed_score_or_artifacts_are_excluded_not_run(planted, tmp_path):
    """Exclusion rule (a): missing kept artifacts, or a committed row without a stage-1 score, is
    reported and excluded by the plan. run_unit runs nothing for such a unit (never an error row,
    never a failed code check), and the run loop does not count it as a failure."""
    results = tmp_path / "results"
    shutil.copytree(planted["results"], results)
    p = results / EXP / "rows" / "cxg__mock__homo-sapiens__t0.json"
    row = json.loads(p.read_text())
    del row["m_pearson_motif_auto_auprc_ratio"]
    p.write_text(json.dumps(row))
    bare = tmp_path / "bare"
    shutil.copytree(planted["run"], bare, ignore=shutil.ignore_patterns("artifacts"))
    no_score = {u["task"]: u for u in ms.plan_units([f"{EXP}={planted['run']}"], results)}["cxg__mock__homo-sapiens__t0"]
    no_files = {u["task"]: u for u in ms.plan_units([f"{EXP}={bare}"], planted["results"])}["cxg__mock__homo-sapiens__t1"]
    assert no_score["plan_flags"] == ["MISSING committed scores for ['pearson_motif']"]
    assert no_files["plan_flags"] and no_files["plan_flags"][0].startswith("MISSING artifacts")
    for unit in (no_score, no_files):
        out, seeds = ms.run_unit(unit, str(planted["refs"]), n_perm=2)
        assert out["status"] == "excluded" and seeds == {} and "ctx" not in out and out["flags"] == unit["plan_flags"]
        status = ms._run_and_save((unit, str(planted["refs"]), str(tmp_path / "o"), 2, None, "c0ffee0"))[1]
        assert status.startswith("excluded: MISSING") and not status.startswith(("error", "check_failed"))


def test_resubmission_skips_finished_rows_and_retries_runtime_errors(tmp_path, monkeypatch):
    """A resubmission with the same commit (rule c) skips completed and excluded rows and reruns a
    row with a runtime error (check_failed: test_rows_carry_no_home_path_and_check_failures_are_not_rerun)."""
    calls, status = [], ["completed"]

    def fake(unit, refs, n_perm, workroot):
        calls.append(unit["id"])
        return {"id": unit["id"], "exp": unit["exp"], "task": unit["task"], "status": status[0], "flags": []}, {}
    monkeypatch.setattr(ms, "run_unit", fake)
    args = ({"id": "u1", "exp": EXP, "task": "t"}, None, str(tmp_path / "o"), 2, None, "c0ffee0")
    for first in ("completed", "excluded"):
        (tmp_path / "o" / "rows" / "u1.json").unlink(missing_ok=True)
        calls.clear()
        status[0] = first
        assert ms._run_and_save(args)[1] == first and ms._run_and_save(args)[1] == f"skipped ({first})"
        assert calls == ["u1"]                                         # the second call did not run the unit
    (tmp_path / "o" / "rows" / "u1.json").write_text(json.dumps({"id": "u1", "status": "error", "flags": ["ERROR x"]}))
    calls.clear()
    status[0] = "completed"
    assert ms._run_and_save(args)[1] == "completed" and calls == ["u1"]  # a runtime error is retried


def test_auto_is_baseline_follows_the_committed_candidates(planted, tmp_path, monkeypatch):
    """Record: the report prints on how many datasets the auto pipeline is the baseline pipeline
    (13/20 C6, 12/20 C4 in the committed rows); each row takes it from its committed candidates."""
    results = tmp_path / "results"
    shutil.copytree(planted["results"], results)
    files = sorted((results / EXP / "rows").glob("*.json"))[:2]
    for p, (auto, base) in zip(files, (("baseline", "baseline"), ("context_scope=selected_type", "baseline"))):
        row = json.loads(p.read_text())
        row["auto_candidate"], row["baseline_candidate"] = auto, base
        p.write_text(json.dumps(row))
    units = {u["task"]: u for u in ms.plan_units([f"{EXP}={planted['run']}"], results)}
    assert [_run_one(planted, units[p.stem], monkeypatch)["auto_is_baseline"] for p in files] == [True, False]
    assert not any("auto_is_baseline" in r for r in _rows(planted))    # the mock rows carry no candidates


def test_report_counts_the_datasets_where_auto_is_the_baseline(planted, tmp_path):
    out = tmp_path / "sens"
    shutil.copytree(planted["out"], out)
    for p, flag in zip(sorted((out / "rows").glob("*.json")), (True, True, False, None, None)):
        row = json.loads(p.read_text())
        row.pop("auto_is_baseline", None)
        if flag is not None:
            row["auto_is_baseline"] = flag
        p.write_text(json.dumps(row))
    text = ms.report(out)
    assert json.loads((out / "SENSITIVITY.json").read_text())["sets"][EXP]["auto_is_baseline"] == [2, 3]
    assert "\n- datasets where the auto pipeline is the baseline pipeline (committed rows): 2/3; " in text


def test_run_inputs_and_report_carry_the_platform_without_host_or_home(planted):
    """Record, platform note: RUN_INPUTS.json holds the OS, CPU model, Python, numpy and BLAS of the
    run (no host name, no home path); the report prints them and SENSITIVITY.json keeps them."""
    import platform as _platform
    inputs = json.loads((planted["out"] / "RUN_INPUTS.json").read_text())
    assert {"references", "cistarget", "platform", "invocations"} <= set(inputs)
    plat = inputs["platform"]
    assert {"os", "cpu", "python", "numpy", "blas"} <= set(plat) and plat["os"]
    assert plat["numpy"] == np.__version__ and plat["python"] == _platform.python_version()
    assert [i["platform"] for i in inputs["invocations"]] == [plat]
    node, home, blob = _platform.node(), str(Path.home()), json.dumps(inputs)
    assert "node" not in plat and node not in plat.values()            # no host name (AGENTS.md section 11)
    assert len(node) < 6 or node not in blob                            # a short name can occur by chance
    assert home in ("", "/") or home not in blob
    assert json.loads((planted["out"] / "SENSITIVITY.json").read_text())["platform"] == plat
    assert f"\n- platform: {plat['os']}; CPU " in (planted["out"] / "SENSITIVITY.md").read_text()


def test_descriptive_share_of_the_lift_matches_the_decision_ratios(planted):
    """The descriptive 'share of the real lift' of tfmean_db and gene_tf is the ratio the decision
    rule reads (L_art/L, L_gene/L: mean lift of the control over mean lift of the real prior)."""
    import re
    text = (planted["out"] / "SENSITIVITY.md").read_text()
    for arm in ms.ARMS:
        block = text.split(f"### Descriptive, {arm} pipeline")[1].split("\n### ")[0]
        for variant, key in (("tfmean_db", "L_art/L"), ("gene_tf", "L_gene/L")):
            m = re.search(rf"\| pooled \| {variant} \| \d+ \| [^|]+\| [^|]+\| [^|]+\| [^|]+\| (-?[\d.]+) \[", block)
            assert m and float(m.group(1)) == pytest.approx(_arm(planted, arm)["ratios"][key][0], abs=0.006), \
                (arm, variant)


def _tiny_context(root: Path, shape: str):
    """One tiny context for decompose: 6 TFs x 86 genes, every TF x gene pair an edge, and a fake
    cisTarget folder (two databases and a motif table: hg38/hgnc, or mm10/mgi for title-case mouse
    symbols). Shapes: no_motif (no TF has an annotated motif), one_motif (K >= 1 for one TF only),
    few_gold (4 evaluation positives per TF, MIN_POS_TF is 5), no_gene_in_db (no expression gene in
    the databases), mouse (K >= 1 for every TF), shared (TF0/TF1 and TF2/TF3 share a motif, TF4 and TF5
    end with the same set)."""
    from grnbench.eval.metrics import GoldIndex
    rng = np.random.default_rng(3)
    mouse = shape == "mouse"
    tfs = [f"Tf{i}" if mouse else f"TF{i}" for i in range(6)]
    targets = [f"Gn{j}" if mouse else f"G{j}" for j in range(80)]
    genes = tfs + targets
    K = {"no_motif": [0] * 6, "one_motif": [2, 0, 0, 0, 0, 0]}.get(shape, [1, 2, 3, 1, 2, 3])
    db_genes = ([] if shape == "no_gene_in_db" else genes) + [f"Bg{i}" if mouse else f"BG{i}" for i in range(40)]
    cis = root / "cistarget"
    cis.mkdir(parents=True)
    ann = []
    for db in ("mm10_a_fake", "mm10_b_fake") if mouse else ("hg38_a_fake", "hg38_b_fake"):
        motifs = [(f"{db}_m{i}_{j}", t) for i, t in enumerate(tfs) for j in range(K[i])] + [(f"{db}_bg", "NOTINDATA")]
        extra = ([(f"{db}_m0_0", tfs[1]), (f"{db}_m2_0", tfs[3]),                   # TF0/TF1 and TF2/TF3 share one motif
                  (f"{db}_m4_0", tfs[5]), (f"{db}_m4_1", tfs[5])]                   # TF4 and TF5 end with the same set
                 + [(f"{db}_m5_{j}", tfs[4]) for j in range(3)]) if shape == "shared" else []
        ann += motifs + extra
        tab = pd.DataFrame(np.vstack([rng.permutation(len(db_genes)) for _ in motifs]).astype(np.int32),
                           columns=db_genes)
        tab.insert(0, "motifs", [m for m, _ in motifs])
        pf.write_feather(tab, str(cis / f"{db}.genes_vs_motifs.rankings.feather"))
    pd.DataFrame(ann, columns=["#motif_id", "gene_name"]).to_csv(
        cis / f"motifs-fake.{'mgi' if mouse else 'hgnc'}-m0.001-o0.0.tbl", sep="\t", index=False)
    n_gold = ms.MIN_POS_TF - 1 if shape == "few_gold" else 10
    gold = {(t, str(g)) for t in tfs for g in rng.choice(targets, n_gold, replace=False)}
    E = pd.DataFrame([(t, g) for t in tfs for g in genes if g != t], columns=["source", "target"])
    E = E.assign(weight=rng.random(len(E)))
    expr = pd.DataFrame(rng.normal(size=(30, len(genes))), columns=genes)
    return cis, E, expr, GoldIndex.coerce(gold), {g.upper() for g in targets}


@pytest.mark.parametrize("shape", ["no_motif", "one_motif", "few_gold", "no_gene_in_db", "mouse"])
def test_degenerate_contexts_give_nan_or_equal_values_not_errors(shape, tmp_path, monkeypatch):
    """Shapes a real context can take: decompose completes, the built-in checks hold (NaN equals NaN),
    and each shape gives what section 1 of the record says (exclusion rule d, the gene_tf bullet)."""
    from grnbench.infer.cistarget import guess_species
    cis, E, expr, gi, labelable = _tiny_context(tmp_path, shape)
    monkeypatch.setenv(ct.ENV_DIR, str(cis))
    ct.load_annotation.cache_clear()
    rec, seeds = ms.decompose(E, expr, E[["source", "target"]], gi, labelable, n_perm=2)
    ct.load_annotation.cache_clear()
    assert rec["builtin_checks"] == "held" and rec["n_tfs"] == 6 and rec["n_genes"] == 86
    assert len(seeds["perm_db_pooled"]) == 2

    def same(a, b):
        return rec[a] == pytest.approx(rec[b], abs=1e-12)
    if shape == "no_motif":              # the prior is W everywhere: no lift, no U''
        assert rec["n_tfs_no_motif"] == 6 and rec["n_tfs_macro_R"] == 0 and rec["gold_in_R"] == 0
        assert all(same(f"{v}_pooled", "base_pooled") for v in ("real", "perm_db", "gene_tf", "tfmean_db"))
        assert np.isnan(rec["base_macro_R"]) and np.isnan(rec["real_macro_R"])
        assert np.isfinite(rec["real_macro"]) and rec["n_tfs_macro"] == 6
    elif shape == "one_motif":           # record: with one TF with K >= 1, gene_tf equals the real prior
        assert rec["n_tfs_no_motif"] == 5 and rec["n_tfs_macro_R"] == 1
        assert rec["gene_tf_pooled"] == rec["real_pooled"] and rec["gene_tf_macro_R"] == rec["real_macro_R"]
        assert np.isfinite(rec["real_macro_R"])
    elif shape == "few_gold":            # rule (d): no TF with 5 positives, every macro metric undefined
        assert rec["n_tfs_macro"] == rec["n_tfs_macro_R"] == 0 and rec["gold_in_universe"] == 6 * (ms.MIN_POS_TF - 1)
        assert all(np.isnan(rec[f"{v}_{m}"]) for v in ("base", "real", "gene_tf") for m in ("macro", "macro_T", "macro_R"))
        assert np.isfinite(rec["real_pooled"])
    elif shape == "no_gene_in_db":       # every gene tied at W: the real prior adds nothing, U' is empty
        assert rec["genes_not_in_db"] == 86 and rec["gold_in_T"] == 0
        assert same("real_pooled", "base_pooled") and np.isnan(rec["real_macro_R"]) and np.isnan(rec["real_pooled_T"])
    else:                                # mouse: guess_species reads the mm10 databases and the mgi table
        assert guess_species(expr.columns) == "mouse" and rec["n_tfs_no_motif"] == 0
        assert rec["n_tfs_macro_R"] == 6 and np.isfinite(rec["real_macro_R"])
        assert rec["real_pooled"] != pytest.approx(rec["base_pooled"], abs=1e-6)


def test_a_mouse_dataset_runs_through_the_mouse_reference_and_databases(tmp_path, monkeypatch):
    """12 of the 40 registered datasets are mouse. run_unit takes the organism from the committed row,
    reads collectri_mouse.tsv for the labelable targets and the mm10/mgi files (guess_species, the
    registered behaviour), and gate 1 holds on a mock written by score_models with them."""
    mock = build_mock(tmp_path, plant=True, n_tasks=1, organism="mouse")
    (unit,) = ms.plan_units([f"{EXP}={mock['run']}"], mock["results"])
    assert unit["plan_flags"] == [] and unit["task"] == "cxg__mock__mus-musculus__t0"
    row = _run_one(mock, unit, monkeypatch)
    assert row["status"] == "completed" and row["flags"] == [] and max(row["gate"].values()) == 0.0
    assert list(row["reference_sha256"]) == ["collectri_mouse.tsv"]
    assert len(row["ctx"]) == 2 and all(c["builtin_checks"] == "held" for c in row["ctx"])
    assert all(c["n_tfs_no_motif"] == N_NO_MOTIF and c["genes_ensembl"] == N_ENS for c in row["ctx"])
    assert all(c["gold_in_R"] > 0 and np.isfinite(c["real_macro_R"]) for c in row["ctx"])
    assert ms.preflight_problems([unit], str(mock["refs"])) == []
    human_only = tmp_path / "refs_human"
    human_only.mkdir()
    shutil.copy(mock["refs"] / "collectri_mouse.tsv", human_only / "collectri_human.tsv")
    assert ms.preflight_problems([unit], str(human_only)) == [f"--refs {ms._portable(human_only)}: "
                                                              "no collectri_mouse.tsv"]


# motif-sharing diagnostics and gene_tf_np (record, dated note 2026-10-06) ------------------------- #

def test_shared_motifs_are_read_from_the_annotation_and_the_matrix_is_symmetric(tmp_path, monkeypatch):
    cis, E, expr, gi, labelable = _tiny_context(tmp_path, "shared")
    monkeypatch.setenv(ct.ENV_DIR, str(cis))
    ct.load_annotation.cache_clear()
    tfs_e = [f"TF{i}" for i in range(6)]
    sets = ms.motif_sets(tfs_e, tuple(expr.columns))
    M, in_db, K, worst = ms.motif_parts(tfs_e, tuple(expr.columns))
    ct.load_annotation.cache_clear()
    assert [bool(st) for st in sets] == [bool(k) for k in K] and all(K > 0)         # K > 0 exactly when not empty
    assert sets[0] <= sets[1] and sets[2] & sets[3] == {"hg38_a_fake_m2_0", "hg38_b_fake_m2_0"}
    assert sets[4] == sets[5] and len(sets[4]) == 10
    sh = ms.shared_matrix(sets)
    want = np.zeros((6, 6), bool)
    for a, b in ((0, 1), (2, 3), (4, 5)):
        want[a, b] = want[b, a] = True
    assert np.array_equal(sh, want) and not sh.diagonal().any()
    assert ms.shared_matrix([frozenset(), frozenset("a"), frozenset("ab"), frozenset()]).tolist() == \
        [[False] * 4, [False, False, True, False], [False, True, False, False], [False] * 4]


def test_gene_tf_np_leaves_the_tf_and_its_motif_partners_out_of_the_propensity():
    rng = np.random.default_rng(5)
    M = rng.permutation(8 * 12).reshape(8, 12).astype(float) + 1
    in_db = np.array([True] * 10 + [False] * 2)
    K = np.array([1, 2, 1, 3, 1, 2, 0, 1])                                    # TF 6 has no motif
    sets = [frozenset({"a"}), frozenset({"a", "b"}), frozenset({"c"}), frozenset({"b", "d"}),
            frozenset({"e"}), frozenset({"e", "f"}), frozenset(), frozenset({"d"})]
    shared = ms.shared_matrix(sets)
    W = 99.0
    G = ms.variant_matrix("gene_tf_np", M, in_db, K, W, shared=shared)
    cols, rows = np.flatnonzero(in_db), np.flatnonzero(K > 0)
    for i in rows:                                                            # independent brute force
        others = [j for j in rows if j != i and not shared[i, j]]
        prop = M[np.ix_(others, cols)].mean(0)
        want = np.full(M.shape[1], W)
        want[cols[np.argsort(prop, kind="stable")]] = np.sort(M[i, cols], kind="stable")
        assert np.array_equal(G[i], want), i
    assert (G[6] == W).all() and (G[:, ~in_db] == W).all()                    # motif-less TF and the mask
    # TF 0's partner (TF 1) is out of its propensity; gene_tf would still average it in
    assert not np.array_equal(G[0], ms.variant_matrix("gene_tf", M, in_db, K, W)[0])


def test_gene_tf_np_is_gene_tf_without_partners_and_the_real_row_without_anyone_else():
    M = np.array([[1.0, 5.0, 3.0, 9.0], [2.0, 6.0, 4.0, 9.0], [7.0, 8.0, 1.0, 9.0], [9.0, 9.0, 9.0, 9.0]])
    in_db, K = np.array([True, True, True, False]), np.array([1, 2, 1, 0])
    none = np.zeros((4, 4), bool)
    assert np.array_equal(ms.variant_matrix("gene_tf_np", M, in_db, K, 9, shared=none),
                          ms.variant_matrix("gene_tf", M, in_db, K, 9))     # no partner: identical
    everyone = np.ones((4, 4), bool) & ~np.eye(4, dtype=bool)
    G = ms.variant_matrix("gene_tf_np", M, in_db, K, 9, shared=everyone)
    assert np.array_equal(G[:3, :3], M[:3, :3])                               # nobody left: the real row stays
    with pytest.raises(ValueError, match="shared-motif matrix"):
        ms.variant_matrix("gene_tf_np", M, in_db, K, 9)


def test_decompose_reports_the_share_of_tfs_with_a_motif_partner(tmp_path, monkeypatch):
    out = {}
    for shape in ("shared", "mouse", "no_motif"):
        cis, E, expr, gi, labelable = _tiny_context(tmp_path / shape, shape)
        monkeypatch.setenv(ct.ENV_DIR, str(cis))
        ct.load_annotation.cache_clear()
        out[shape], _ = ms.decompose(E, expr, E[["source", "target"]], gi, labelable, n_perm=2)
        ct.load_annotation.cache_clear()
    rec = out["shared"]
    assert rec["builtin_checks"] == "held" and rec["n_tfs_motif"] == 6
    assert rec["frac_tfs_shared_motifs"] == 1.0 and rec["frac_tfs_same_motifs"] == pytest.approx(2 / 6)
    assert np.isfinite(rec["gene_tf_np_pooled"]) and np.isfinite(rec["gene_tf_np_macro_R"])
    rec = out["mouse"]                                                        # each TF has its own motifs
    assert rec["frac_tfs_shared_motifs"] == 0 and rec["frac_tfs_same_motifs"] == 0
    assert all(rec[f"gene_tf_np_{m}"] == pytest.approx(rec[f"gene_tf_{m}"], abs=1e-12) for m in ms.METRICS)
    rec = out["no_motif"]
    assert rec["n_tfs_motif"] == 0 and np.isnan(rec["frac_tfs_shared_motifs"]) and np.isnan(rec["frac_tfs_same_motifs"])


def test_report_prints_the_motif_sharing_diagnostics_gene_tf_np_and_the_n_of_each_quantity(planted):
    text = ms.report(planted["out"])
    for key in ("n_tfs_motif", "frac_tfs_same_motifs", "frac_tfs_shared_motifs"):
        assert f"| {key} |" in text
    assert "| pooled | gene_tf_np |" in text and "gene_tf_np " in text.split("components (pooled lift of each control)")[1]
    for quantity in ("registered lift L", "target-identity gap (T1) D", "design-effect control L_art",
                     "gene-propensity control L_gene"):
        line = next(l for l in text.splitlines() if l.startswith(f"- {quantity}"))
        assert line.endswith(f"(n = {N_TASKS})"), line
    arm = _arm(planted, "auto")
    assert arm["n_quantities"] == {"L": N_TASKS, "D": N_TASKS, "L_art": N_TASKS, "L_gene": N_TASKS}
