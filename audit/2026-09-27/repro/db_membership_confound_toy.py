"""Motif prior with NO motif information (random ranks) still 'lifts' Pearson when the
universe contains targets that have no symbol (absent from the motif DB and from the
symbol-keyed reference)."""
import os, tempfile
from pathlib import Path
import numpy as np, pandas as pd, pyarrow.feather as pf
from grnbench.infer import cistarget as ct
from grnbench.infer.classical import pearson_signed
from grnbench.eval.metrics import evaluate_ranking

EFF = float(os.environ.get("EFF", "0.6"))
def run(seed, frac_ensg=0.15):
    rng = np.random.default_rng(seed)
    n, n_tf, n_sym = 600, 15, 400
    n_ens = int(round(frac_ensg * (n_sym) / (1 - frac_ensg)))
    tfs = [f"TF{i}" for i in range(n_tf)]
    sym = [f"GENE{i}" for i in range(n_sym)]
    ens = [f"ENSG{i:011d}" for i in range(n_ens)]
    Z = rng.normal(size=(n, n_tf))
    X = {t: Z[:, i] for i, t in enumerate(tfs)}
    gold = set()
    # each TF drives 12 symbol targets; every no-symbol gene is also driven by one TF
    for g in sym:
        X[g] = rng.normal(size=n)
    for i, t in enumerate(tfs):
        for g in rng.choice(sym, 12, replace=False):
            X[g] = X[g] + EFF * Z[:, i]; gold.add((t, g))
    for g in ens:
        X[g] = rng.normal(size=n) + EFF * Z[:, rng.integers(n_tf)]
    expr = pd.DataFrame(X)
    root = Path(tempfile.mkdtemp())
    db_genes = tfs + sym + [f"BG{i}" for i in range(2000)]     # symbol genes only; no ENSG ids
    motifs = [f"m{i}" for i in range(40)]
    ranks = np.array([rng.permutation(len(db_genes)) for _ in motifs], dtype=np.int32)  # NO information
    tab = pd.DataFrame(ranks, columns=db_genes); tab["motifs"] = motifs
    pf.write_feather(tab, str(root / "hg38_rand.genes_vs_motifs.rankings.feather"))
    pd.DataFrame({"#motif_id": motifs[:n_tf], "gene_name": tfs}).to_csv(
        root / "motifs-rand.hgnc-m0.001-o0.0.tbl", sep="\t", index=False)
    ct.load_annotation.cache_clear()
    genes = list(expr.columns)
    uni = pd.DataFrame([(s, g) for s in tfs for g in genes if s != g], columns=["source", "target"])
    base = pearson_signed(expr, tf_names=tfs)
    base = base[base.source != base.target].rename(columns=dict(zip(base.columns[:3], ["source","target","weight"])))
    comb = ct.motif_rank_product(base[["source","target","weight"]], expr, species="human", root=root)
    prior = ct.motif_prior(expr, tf_names=tfs, species="human", root=root)
    r = lambda e: evaluate_ranking(e[["source","target","weight"]], gold, with_curve=False, universe=uni)["auprc_ratio"]
    return r(base), r(comb), r(prior)

for frac in (0.0, 0.15):
    res = np.array([run(s, frac) for s in range(10)])
    d = res[:, 1] - res[:, 0]
    print(f"no-symbol share {frac:.2f}: pearson {res[:,0].mean():.3f}  pearson_motif(random motifs) {res[:,1].mean():.3f}  "
          f"motif_prior(random) {res[:,2].mean():.3f}  lift {d.mean():+.3f} (min {d.min():+.3f}, max {d.max():+.3f})")
