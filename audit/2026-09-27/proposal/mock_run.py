"""Build a synthetic run folder with the exact layout real_bench.py --keep-artifacts writes
(artifacts/<task>/{reference_*.tsv, decisions.json, <arm>/ctx<i>/{processed.h5ad, context.json,
<model>.tsv.gz}} + rows/<task>.json), by calling real_bench.score_models itself, so the
sensitivity script can be checked against 'registered' numbers it did not compute.
Fake cisTarget dir (two hg38_* DBs + hgnc .tbl), fake collectri_human.tsv."""
import json, os, sys
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, pyarrow.feather as pf

ROOT = Path(sys.argv[1]); ROOT.mkdir(parents=True, exist_ok=True)
N_TASKS = int(sys.argv[2]) if len(sys.argv) > 2 else 4
MODELS = (sys.argv[3] if len(sys.argv) > 3 else "pearson_signed,pearson_motif").split(",")
REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "benchmarks"))
import real_bench  # noqa: E402
real_bench.ARTIFACT_TOP_EDGES = 3000          # force truncation, as the real 20,000 cap does
from grnbench.screen import split_reference  # noqa: E402

rng = np.random.default_rng(7)
n_tf, n_sym, n_ens = 25, 400, 40
tfs = [f"TF{i}" for i in range(n_tf)]
sym = [f"G{i}" for i in range(n_sym)]
ens = [f"ENSG{i:011d}" for i in range(n_ens)]
genes = tfs + sym + ens
K = np.r_[np.zeros(5, int), rng.integers(1, 9, n_tf - 5)]            # 5 motif-less TFs
not_in_db = set(rng.choice(sym, 30, replace=False))                   # symbol genes absent from DB
db_genes = [g for g in tfs + sym if g not in not_in_db] + [f"BG{i}" for i in range(2000)]
PLANT = int(os.environ.get("PLANT", 1)); EFF = float(os.environ.get("EFF", 0.35))
deg = [25] * n_tf if PLANT else [int(8 + 4 * k) for k in K]          # null mock: degree grows with K
true = {t: list(rng.choice(sym, d, replace=False)) for t, d in zip(tfs, deg)}

cis = ROOT / "cistarget"; cis.mkdir(exist_ok=True)
ann = []
for db in ("hg38_a_fake", "hg38_b_fake"):
    ids, rows = [], []
    for i, t in enumerate(tfs):
        for j in range(K[i]):
            r = rng.permutation(len(db_genes)).astype(np.int32)
            pos = {g: k for k, g in enumerate(db_genes)}
            hits = [pos[g] for g in true[t][:12] if g in pos] if PLANT else []   # half the targets carry the motif
            for h in hits:                                                # swap planted targets to top ranks
                a = np.where(r == rng.integers(0, 200))[0][0]; r[a], r[h] = r[h], r[a]
            ids.append(f"{db}_m{i}_{j}"); rows.append(r); ann.append((f"{db}_m{i}_{j}", t))
    for j in range(40):
        ids.append(f"{db}_bg{j}"); rows.append(rng.permutation(len(db_genes)).astype(np.int32))
    tab = pd.DataFrame(np.vstack(rows), columns=db_genes); tab.insert(0, "motifs", ids)
    pf.write_feather(tab, str(cis / f"{db}.genes_vs_motifs.rankings.feather"))
pd.DataFrame(ann, columns=["#motif_id", "gene_name"]).to_csv(cis / "motifs-fake.hgnc-m0.001-o0.0.tbl", sep="\t", index=False)
os.environ["GRNBENCH_CISTARGET_DIR"] = str(cis)

refs = ROOT / "references"; refs.mkdir(exist_ok=True)
ref_rows = [(t, g) for t in tfs for g in true[t]] + [(t, f"OTHER{k}") for t in tfs for k in range(3)]
pd.DataFrame(ref_rows, columns=["source", "target"]).to_csv(refs / "collectri_human.tsv", sep="\t", index=False)


class FakePlan:
    def __init__(self, strategy): self.s = strategy
    def value(self, key, default=None): return self.s if key == "strategy" else default


def make_context(n_cells, seed, counts):
    r = np.random.default_rng(seed)
    act = r.normal(size=(n_cells, n_tf))
    X = {t: act[:, i] for i, t in enumerate(tfs)}
    for g in sym + ens:
        X[g] = r.normal(size=n_cells) * 0.8
    for i, t in enumerate(tfs):
        for g in true[t]:
            X[g] = X[g] + EFF * act[:, i]
    for g in ens:
        X[g] = X[g] + EFF * act[:, r.integers(n_tf)]
    lin = pd.DataFrame(X)[genes].to_numpy()
    if counts:
        mat = r.poisson(np.exp(0.5 * lin + 1.0)).astype(np.float32)
    else:
        mat = np.log1p(r.poisson(np.exp(0.5 * lin + 1.0))).astype(np.float32)
    return ad.AnnData(mat, var=pd.DataFrame(index=genes), obs=pd.DataFrame(index=[f"c{k}" for k in range(n_cells)]))


out = ROOT / "run"; (out / "rows").mkdir(parents=True, exist_ok=True)
work = ROOT / "work"; work.mkdir(exist_ok=True)
tasks = []
for ti in range(N_TASKS):
    task = f"cxg__mock__homo-sapiens__t{ti}"
    pairs = {(t, g) for t in tfs for g in true[t]}
    sel, ev = split_reference(pairs, 0)
    recs = {"baseline": [], "candX": []}
    for cand, n_ctx, counts in (("baseline", 1, False), ("candX", 2, True)):
        for c in range(n_ctx):
            p = work / f"{task}_{cand}_{c}.h5ad"
            make_context(300 + 50 * c, 100 * ti + 10 * c + (cand == "candX"), counts).write_h5ad(p)
            recs[cand].append({"processed": str(p), "context": f"ctx{c}"})
    plans = {"baseline": FakePlan("cpm_log1p"), "candX": FakePlan("raw_counts")}
    ta = out / "artifacts" / task; ta.mkdir(parents=True, exist_ok=True)
    for half, pr in (("selection", sel), ("evaluation", ev)):
        pd.DataFrame(sorted(pr), columns=["source", "target"]).to_csv(ta / f"reference_{half}.tsv", sep="\t", index=False)
    (ta / "decisions.json").write_text(json.dumps({"candidates": {"baseline": "baseline", "auto": "candX"},
        "chosen_plans": {"baseline": {"candidate": "baseline", "choices": {"normalization": "cpm_log1p"}},
                         "auto": {"candidate": "candX", "choices": {"normalization": "raw_counts"}}}}))
    row = {"task": task, "reference_edges": len(pairs)}
    row.update(real_bench.score_models(recs, plans, {"baseline": "baseline", "auto": "candX"}, MODELS,
                                       sorted({s for s, _ in pairs}), ev, artifacts_dir=ta))
    (out / "rows" / f"{task}.json").write_text(json.dumps(row, default=float))
    tasks.append({"id": task})
    print(task, {k: round(v, 4) for k, v in row.items() if k.endswith("auto_auprc_ratio")}, flush=True)
(out / "real_tasks.json").write_text(json.dumps(tasks))
