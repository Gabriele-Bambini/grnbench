"""P9 diagnosis: is the K562 shift to decreases under `core` a global shrink of transcription?
For each usable TF (core truth): median UMI_count / core_scale_factor of its cells relative to the
controls' median, and the share of its responder rows that are decreases. obs only, no X read."""
import csv
import json
import sys

import h5py
import numpy as np

h5, truth = sys.argv[1], sys.argv[2]


def col(g, name):
    o = g[name]
    if isinstance(o, h5py.Group):
        cats = o["categories"][()].astype(str)
        return cats[o["codes"][()]]
    return o[()]


with h5py.File(h5, "r") as f:
    obs = f["obs"]
    pert = col(obs, "perturbation")
    umi = col(obs, "UMI_count").astype(float)
    sf = col(obs, "core_scale_factor").astype(float)
adj = umi / sf
ctrl = np.median(adj[pert == "control"])
units = [u for u in csv.DictReader(open(f"{truth}/units.tsv"), delimiter="\t") if u["usable"] == "True"]
neg, tot = {}, {}
for r in csv.DictReader(open(f"{truth}/responders.tsv"), delimiter="\t"):
    tot[r["tf"]] = tot.get(r["tf"], 0) + 1
    neg[r["tf"]] = neg.get(r["tf"], 0) + (float(r["log2fc"]) < 0)
ratio, share = [], []
for u in units:
    m = pert == u["tf"]
    ratio.append(np.median(adj[m]) / ctrl)
    share.append(neg.get(u["tf"], 0) / max(tot.get(u["tf"], 1), 1))
ratio, share = np.array(ratio), np.array(share)
from scipy.stats import spearmanr
print(json.dumps({
    "usable_tfs": len(units),
    "median_rel_core_adjusted_umi": round(float(np.median(ratio)), 3),
    "tfs_with_rel_below_0.9": int((ratio < 0.9).sum()),
    "median_decreased_share_per_tf": round(float(np.median(share)), 3),
    "spearman_rel_vs_decreased_share": [round(float(x), 4) for x in spearmanr(ratio, share)],
}))
