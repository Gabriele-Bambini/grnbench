"""Does the min-over-K motif rank + pooled rank product create a TF-level reweighting that
moves pooled AUPRC with zero edge-level motif information? Uses the real grnbench code."""
import os, sys, tempfile
from pathlib import Path
import numpy as np, pandas as pd
import pyarrow.feather as pf
from grnbench.infer import cistarget as ct
from grnbench.eval.metrics import GoldIndex, evaluate_ranking, average_precision

N_TF, N_TGT, N_DB_EXTRA, N_CELLS = 40, 400, 4000, 400
COEF = float(os.environ.get('COEF', '0.03'))

def build(rng, K, edge_info=False, gold=None):
    root = Path(tempfile.mkdtemp())
    tfs = [f"TF{i}" for i in range(N_TF)]
    tg = [f"G{j}" for j in range(N_TGT)]
    genes = tfs + tg + [f"BG{j}" for j in range(N_DB_EXTRA)]
    n = len(genes)
    motifs, ann, rows = [], [], []
    for i, tf in enumerate(tfs):
        for k in range(K[i]):
            m = f"m{i}_{k}"
            r = rng.permutation(n)
            if edge_info and k == 0:
                # put this TF's gold targets in the top 20% of the genome (one real motif)
                gset = [genes.index(t) for (s, t) in gold if s == tf]
                top = rng.choice(int(0.2 * n), len(gset), replace=False)
                rest = np.setdiff1d(np.arange(n), top); rng.shuffle(rest)
                r = np.empty(n, int); others = np.setdiff1d(np.arange(n), gset)
                r[gset] = top; r[others] = rest
            motifs.append(m); ann.append((m, tf)); rows.append(r)
    t = pd.DataFrame(np.array(rows, dtype=np.int32), columns=genes); t["motifs"] = motifs
    pf.write_feather(t, str(root / "hg38_fake.genes_vs_motifs.rankings.feather"))
    pd.DataFrame(ann, columns=["#motif_id", "gene_name"]).to_csv(
        root / "motifs-fake.hgnc-m0.001-o0.0.tbl", sep="\t", index=False)
    ct.load_annotation.cache_clear()
    return root, tfs, tg

def one(seed, K, deg, edge_info=False):
    rng = np.random.default_rng(seed)
    tfs = [f"TF{i}" for i in range(N_TF)]; tg = [f"G{j}" for j in range(N_TGT)]
    gold = set()
    for i, tf in enumerate(tfs):
        for j in rng.choice(N_TGT, deg[i], replace=False):
            gold.add((tf, tg[j]))
    root, _, _ = build(rng, K, edge_info, gold)
    # expression: each gold edge gets a weak Pearson signal
    X = {tf: rng.normal(size=N_CELLS) for tf in tfs}
    Y = {g: rng.normal(size=N_CELLS) for g in tg}
    for s, t in gold:
        Y[t] = Y[t] + COEF * X[s]
    expr = pd.DataFrame({**X, **Y})
    Z = (expr - expr.mean()) / expr.std(ddof=0)
    C = (Z[tfs].T @ Z[tg]) / N_CELLS
    edges = C.stack().rename_axis(["source", "target"]).reset_index(name="w")
    edges["weight"] = edges["w"].abs(); edges = edges[["source", "target", "weight"]]
    universe = edges[["source", "target"]]
    gi = GoldIndex.coerce(gold)
    comb = ct.motif_rank_product(edges, expr, species="human", root=root)
    prior = ct.motif_prior(expr[tfs + tg], tf_names=tfs, species="human", root=root)
    prior = prior[prior.target.isin(tg)]
    res = {}
    for name, e in (("pearson", edges), ("pearson_motif", comb), ("motif_prior", prior)):
        r = evaluate_ranking(e, gi, with_curve=False, universe=universe)
        # per-TF (macro) AUPRC ratio: within-TF ranking only
        ee = e.copy(); ee["y"] = [(s, t) in gold for s, t in zip(ee.source, ee.target)]
        macro = np.mean([average_precision(g.y.to_numpy(), g.weight.to_numpy()) / g.y.mean()
                         for _, g in ((k, v.sort_values('weight', ascending=False, kind='stable')) for k, v in ee.groupby('source')) if g.y.any()])
        res[name] = (r["auprc_ratio"], macro)
    return res

def scenario(label, Kfun, degfun, edge_info=False, reps=10):
    pooled = {"pearson": [], "pearson_motif": [], "motif_prior": []}
    macro = {k: [] for k in pooled}
    for s in range(reps):
        rng = np.random.default_rng(1000 + s)
        deg = degfun(rng); K = Kfun(rng, deg)
        r = one(s, K, deg, edge_info)
        for k, (p, m) in r.items():
            pooled[k].append(p); macro[k].append(m)
    lift_p = np.array(pooled["pearson_motif"]) - np.array(pooled["pearson"])
    lift_m = np.array(macro["pearson_motif"]) - np.array(macro["pearson"])
    print(f"{label:55s} pooled: pearson {np.mean(pooled['pearson']):.3f} motif {np.mean(pooled['pearson_motif']):.3f} "
          f"prior {np.mean(pooled['motif_prior']):.3f} lift {lift_p.mean():+.3f} (pos {int((lift_p>0).sum())}/{reps}) | "
          f"per-TF lift {lift_m.mean():+.3f} (pos {int((lift_m>0).sum())}/{reps})", flush=True)

deg_skew = lambda rng: np.clip(np.round(rng.lognormal(2.5, 0.9, N_TF)), 3, 200).astype(int)
# K rises with degree (rank corr ~ 1), K independent, K falls with degree
K_corr = lambda rng, d: np.clip(1 + (np.argsort(np.argsort(d)) * 30) // N_TF, 1, 40)
K_indep = lambda rng, d: rng.permutation(np.clip(1 + (np.argsort(np.argsort(d)) * 30) // N_TF, 1, 40))
K_anti = lambda rng, d: np.clip(1 + ((N_TF - 1 - np.argsort(np.argsort(d))) * 30) // N_TF, 1, 40)
scenario("no edge info, K increases with TF degree", K_corr, deg_skew)
scenario("no edge info, K independent of TF degree", K_indep, deg_skew)
scenario("no edge info, K decreases with TF degree", K_anti, deg_skew)
scenario("edge info (1 true motif/TF), K independent", K_indep, deg_skew, edge_info=True)
