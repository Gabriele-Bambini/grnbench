"""Read-only metadata checks for EXP-20261016-LM-1 (P1, P5, P7a). h5py only, no X read."""
import collections
import json
import os
import sys

import h5py
import numpy as np

ROOT = os.path.expanduser("~/grnbench_data_perturb")
FILES = {"rpe1": "ReplogleWeissman2022_rpe1.h5ad", "k562": "ReplogleWeissman2022_K562_gwps.h5ad"}
COLLECTRI = os.path.expanduser("~/grnbench_data_fresh1/references/collectri_human.tsv")


def col(group, name):
    obj = group[name]
    if isinstance(obj, h5py.Group) and "categories" in obj:
        cats = obj["categories"][()].astype(str)
        codes = obj["codes"][()]
        return np.where(codes >= 0, cats[np.clip(codes, 0, None)], "NA")
    v = obj[()]
    return v.astype(str) if v.dtype.kind in "SO" else v


def index_of(group):
    key = group.attrs.get("_index", "_index")
    return col(group, key if isinstance(key, str) else key.decode())


sources = set()
with open(COLLECTRI) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    si = header.index("source") if "source" in header else 0
    for line in fh:
        sources.add(line.split("\t")[si].upper())

out = {"collectri_sources": len(sources)}
perturbed = {}
for key, fname in FILES.items():
    with h5py.File(os.path.join(ROOT, fname), "r") as f:
        obs, var = f["obs"], f["var"]
        X = f["X"]
        info = {
            "n_obs": int(f["obs"][obs.attrs.get("_index", "_index")].shape[0]),
            "obs_columns": list(obs.keys()),
            "var_columns": list(var.keys()),
            "X": dict(X.attrs) if isinstance(X, h5py.Group) else str(X.dtype),
            "layers": list(f["layers"].keys()) if "layers" in f else [],
        }
        info["X"] = {k: str(v) for k, v in info["X"].items()} if isinstance(info["X"], dict) else info["X"]
        if isinstance(X, h5py.Group) and "data" in X:
            d = X["data"][: 200000]
            info["X_data_dtype"] = str(d.dtype)
            info["X_first200k_all_integer"] = bool(np.all(np.mod(d, 1) == 0))
        elif isinstance(X, h5py.Dataset):
            d = X[:200, :]
            info["X_dense_dtype"] = str(X.dtype)
            info["X_dense_shape"] = list(X.shape)
            info["X_first200rows_all_integer"] = bool(np.all(np.mod(d, 1) == 0))
        genes = index_of(var)
        info["n_vars"] = int(len(genes))
        sym = None
        for c in ("gene_name", "gene_symbol", "symbol", "gene"):
            if c in var:
                sym = col(var, c)
                info["var_symbol_column"] = c
                break
        symbols = set(np.char.upper((sym if sym is not None else genes).astype(str)))
        info["var_index_examples"] = genes[:3].tolist()
        info["collectri_sources_in_genes"] = len(sources & symbols)
        for c in obs.keys():
            if c in ("__categories",):
                continue
            try:
                v = col(obs, c)
            except Exception as e:  # noqa: BLE001
                info.setdefault("obs_unreadable", []).append(f"{c}: {e}")
                continue
            if v.dtype.kind in "U":
                cnt = collections.Counter(v.tolist())
                info.setdefault("obs_text", {})[c] = {
                    "n_unique": len(cnt),
                    "top": cnt.most_common(4),
                }
            else:
                info.setdefault("obs_numeric", {})[c] = {
                    "min": float(np.nanmin(v)), "max": float(np.nanmax(v)), "dtype": str(v.dtype)}
        out[key] = info
        for c in ("perturbation", "gene", "target", "perturbation_gene"):
            if c in obs:
                v = col(obs, c)
                perturbed[key] = set(np.char.upper(v.astype(str)).tolist())
                out[key]["perturbation_column"] = c
                break
        if key in perturbed:
            pt = perturbed[key]
            out[key]["n_perturbed_labels"] = len(pt)
            out[key]["perturbed_collectri_tfs"] = len(pt & sources)
            out[key]["perturbed_tfs_in_genes"] = len(pt & sources & symbols)
if len(perturbed) == 2:
    a, b = perturbed["rpe1"], perturbed["k562"]
    out["overlap"] = {
        "perturbed_both": len(a & b),
        "rpe1_not_in_k562": len(a - b),
        "collectri_tfs_perturbed_both": len(a & b & sources),
    }
json.dump(out, sys.stdout, indent=1, default=str)
