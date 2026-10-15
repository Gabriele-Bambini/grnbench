"""score_models: a model that fails in ONE context of a dataset keeps a finite dataset score
(mean over the other contexts), so the paired comparison uses different context sets per arm
instead of dropping the dataset from that pair (registered exclusion rule (c))."""
import sys, tempfile
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad
REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "benchmarks"))
import real_bench as rb
import grnbench.infer as gi

rng = np.random.default_rng(0)
tmp = Path(tempfile.mkdtemp())
genes = [f"TF{i}" for i in range(5)] + [f"G{i}" for i in range(60)]
recs = []
for c in range(2):
    X = rng.poisson(3, size=(200, len(genes))).astype(np.float32)
    X[:, 5:25] += X[:, [0]] * (1 + c)          # context-specific signal
    p = tmp / f"ctx{c}.h5ad"
    ad.AnnData(X, var=pd.DataFrame(index=genes)).write_h5ad(p)
    recs.append({"context": f"type{c}", "processed": str(p)})
class Plan:
    def value(self, k): return "raw_counts"
gold = {(f"TF0", f"G{i}") for i in range(20)} | {("TF1", "G40"), ("TF2", "G41")}
real_run = gi.run_method
calls = {"n": 0}
def flaky(name, adata, **kw):
    if name == "spearman_signed":
        calls["n"] += 1
        if calls["n"] == 2:                        # fails in the second context only
            raise RuntimeError("context-specific failure")
    return real_run(name, adata, **kw)
gi.run_method = flaky
out = rb.score_models({"cand": recs}, {"cand": Plan()}, {"auto": "cand"},
                      ["pearson_signed", "spearman_signed"], [f"TF{i}" for i in range(5)], gold)
for k in sorted(out):
    print(k, out[k])
gi.run_method = real_run
for c in range(2):
    o = rb.score_models({"cand": [recs[c]]}, {"cand": Plan()}, {"auto": "cand"},
                        ["pearson_signed", "spearman_signed"], [f"TF{i}" for i in range(5)], gold)
    print(f"ctx{c} alone:", {k: round(v, 3) for k, v in o.items() if k.endswith("auprc_ratio")})
