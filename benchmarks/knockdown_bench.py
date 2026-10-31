"""Knockdown truth for EXP-20261016-LM-1, phase 1 (section 1 of experiments/EXP-20261016-LM-1.md).

    python benchmarks/knockdown_bench.py truth --h5ad F --tfs collectri.tsv --out DIR [--seed 0] [--workers 1]
    python benchmarks/knockdown_bench.py units --truth DIR
    python benchmarks/knockdown_bench.py controls --h5ad F --out DIR [--n 3000]
    python benchmarks/knockdown_bench.py subset --h5ad F --tfs collectri.tsv --out OUT.h5ad [--block-rows N]
    python benchmarks/knockdown_bench.py baseline-check --h5ad DIR/controls.h5ad
    python benchmarks/knockdown_bench.py total-shift --h5ad F --truth DIR --out OUT.tsv

The h5ad is opened in backed mode; X (dense, integer counts) is read in sorted row chunks, never whole.
Peak memory is the truth controls (3,000 x genes float32) plus one chunk plus one TF's 200 test cells.
No network is built or scored here.

`subset` is for screens whose dense X is chunked across many rows (K562: 3,886 x 33), where truth's scattered row
reads decompress every chunk many times. It reads X once, in blocks of whole chunk rows, and writes the control cells
plus every cell perturbing a CollecTRI TF as an uncompressed CSR h5ad (original order, obs and var kept). truth and
controls read it unchanged in backed mode (row reads of an uncompressed CSR) and give identical output.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import mannwhitneyu

# Section 1 constants.
SPLIT_SEED, CONTROL_SEED, CELL_CAP_SEED = 0, 1, 2     # split / truth controls / TF cell cap
N_TRUTH_CONTROLS, MAX_CELLS_PER_TF, MAX_INFERENCE_CELLS = 3000, 200, 3000
MIN_KD, MIN_CELLS_PER_CONSTRUCT, MIN_DETECTION = 0.70, 20, 0.05
FDR, MIN_ABS_LOG2FC, PSEUDOCOUNT, MIN_RESPONDERS = 0.05, 0.25, 0.1, 11
NORMALISATION = "counts / total counts per cell x 10,000 (CP10k); log1p for the rank test"
NORMALISATION_CORE = ("counts / core_scale_factor per cell (obs column; Replogle 2022 size factor on core genes); "
                      "log1p for the rank test")
CORE_COL = "core_scale_factor"
CHUNK = 2000


# ---------------------------------------------------------------- reading
def read_rows(X, indices, chunk: int = CHUNK):
    """Yield (row indices, dense float32 block) for sorted ``indices``, ``chunk`` rows per read."""
    ordered = np.sort(np.asarray(indices, dtype=np.int64))
    for start in range(0, len(ordered), chunk):
        idx = ordered[start:start + chunk]
        block = X[idx, :]
        block = block.toarray() if sp.issparse(block) else np.asarray(block)
        yield idx, block.astype(np.float32, copy=False)


def cp10k(block: np.ndarray) -> np.ndarray:
    total = block.sum(axis=1, keepdims=True, dtype=np.float64)
    return (block * np.divide(1e4, total, out=np.zeros_like(total), where=total > 0)).astype(np.float32)


def normalise(block: np.ndarray, idx: np.ndarray, scale) -> np.ndarray:
    """CP10k when ``scale`` is None, else counts / scale[idx] (``scale``: per-row factors of the whole screen)."""
    if scale is None:
        return cp10k(block)
    return (block / scale[idx][:, None]).astype(np.float32)


def core_scale(obs: pd.DataFrame, norm: str):
    """(per-row factors or None, mask of rows with a usable factor). cp10k needs no factor: every row passes."""
    if norm == "cp10k":
        return None, np.ones(len(obs), dtype=bool)
    if CORE_COL not in obs.columns:
        raise SystemExit(f"--norm core needs obs column {CORE_COL!r}")
    f = pd.to_numeric(obs[CORE_COL], errors="coerce").to_numpy(dtype=np.float64)
    good = np.isfinite(f) & (f > 0)
    return np.where(good, f, 1.0), good


# ---------------------------------------------------------------- control split
def control_cells(obs: pd.DataFrame, perturbation_col: str, control_label: str, norm: str) -> np.ndarray:
    """Boolean mask of the control cells the split is cut from. ``truth`` and ``controls`` both call this, so
    their halves cannot diverge: under ``core`` a control without a positive finite factor is dropped first."""
    _, good = core_scale(obs, norm)
    return usable_cells(obs, perturbation_col) & good & (obs[perturbation_col].astype(str) == control_label).to_numpy()


def usable_cells(obs: pd.DataFrame, perturbation_col: str) -> np.ndarray:
    """Boolean mask of cells that carry at most one perturbation (nperts > 1 dropped when the column exists)."""
    ok = obs[perturbation_col].notna().to_numpy()
    if "nperts" in obs.columns:
        ok = ok & (pd.to_numeric(obs["nperts"], errors="coerce").fillna(1) <= 1).to_numpy()
    return ok


def split_controls(obs_names, control_idx: np.ndarray, seed: int = SPLIT_SEED):
    """Two disjoint halves of the controls: barcodes sorted, permuted with default_rng(seed).

    Returns (inference_half, truth_half) as sorted row indices.
    """
    names = np.asarray(obs_names).astype(str)
    ordered = control_idx[np.argsort(names[control_idx], kind="stable")]
    permuted = ordered[np.random.default_rng(seed).permutation(len(ordered))]
    half = len(permuted) // 2
    return np.sort(permuted[:half]), np.sort(permuted[half:])


def cap(indices: np.ndarray, n: int, seed: int) -> np.ndarray:
    """At most ``n`` of ``indices``: ``np.sort(default_rng(seed).choice(...))`` (real_bench rule)."""
    if len(indices) <= n:
        return np.sort(indices)
    return indices[np.sort(np.random.default_rng(seed).choice(len(indices), n, replace=False))]


def bh(p: np.ndarray) -> np.ndarray:
    q = np.full(len(p), np.nan)
    ok = np.flatnonzero(np.isfinite(p))
    if len(ok):
        order = ok[np.argsort(p[ok])]
        ranked = p[order] * len(ok) / np.arange(1, len(ok) + 1)
        q[order] = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1.0)
    return q


# ---------------------------------------------------------------- truth
class Controls:
    """The truth-half controls, read once: log values for the rank test and per-gene summaries."""

    def __init__(self, X, truth_half: np.ndarray, chunk: int, scale=None):
        self.idx = cap(truth_half, N_TRUTH_CONTROLS, CONTROL_SEED)
        rows = [normalise(b, i, scale) for i, b in read_rows(X, self.idx, chunk)]
        self.cp = np.vstack(rows)
        self.mean = self.cp.mean(axis=0, dtype=np.float64)
        self.detected = (self.cp > 0).mean(axis=0)
        self.log = np.log1p(self.cp)


def gene_lookup(var_names) -> tuple[dict, dict]:
    exact, folded = {}, {}
    for i, g in enumerate(map(str, var_names)):
        exact.setdefault(g, i)
        folded.setdefault(g.upper(), []).append(i)
    return exact, {k: v[0] for k, v in folded.items() if len(v) == 1}


def find_gene(lookup, tf: str):
    exact, folded = lookup
    return exact.get(tf, folded.get(tf.upper()))


def tf_unit(X, var_names, lookup, ctrl: Controls, tf: str, cells: np.ndarray, constructs: np.ndarray,
            chunk: int, scale=None) -> tuple[dict, pd.DataFrame]:
    """One TF: construct filter, cell cap, per-cell Wilcoxon, BH, fold-change floor, sign agreement.

    Pass 1 streams all of the TF's cells and keeps only per-construct sums; pass 2 reads the capped cells.
    """
    n_genes = len(var_names)
    unit = {"tf": tf, "n_cells": 0, "n_constructs": 0, "knockdown": np.nan, "n_responders": 0,
            "usable": False, "exclusion_reason": "", "n_cells_total": len(cells),
            "n_constructs_total": int(len(set(constructs))), "sign_test": False}
    empty = pd.DataFrame(columns=["tf", "gene", "log2fc", "q", "n_cells", "kd"])
    t = find_gene(lookup, tf)
    if t is None:
        unit["exclusion_reason"] = "tf_not_measured"
        return unit, empty
    if ctrl.mean[t] <= 0:
        unit["exclusion_reason"] = "tf_not_detected_in_controls"
        return unit, empty

    pos = {int(c): i for i, c in enumerate(cells)}
    sums, counts = {}, {}
    for idx, block in read_rows(X, cells, chunk):
        cp = normalise(block, idx, scale)
        con = constructs[[pos[int(i)] for i in idx]]
        for c in np.unique(con):
            m = con == c
            sums[c] = sums.get(c, 0.0) + cp[m].sum(axis=0, dtype=np.float64)
            counts[c] = counts.get(c, 0) + int(m.sum())
    kd = {c: 1.0 - sums[c][t] / counts[c] / ctrl.mean[t] for c in sums}
    big = [c for c in sums if counts[c] >= MIN_CELLS_PER_CONSTRUCT]
    kept = sorted(c for c in big if kd[c] >= MIN_KD)
    if kept:
        unit["knockdown"] = float(np.mean([kd[c] for c in kept]))
    elif big:
        unit["knockdown"] = float(max(kd[c] for c in big))
    if not kept:
        unit["exclusion_reason"] = "knockdown_below_threshold" if big else "too_few_cells_per_construct"
        return unit, empty

    kept_cells = cells[np.isin(constructs, kept)]
    use = cap(kept_cells, MAX_CELLS_PER_TF, CELL_CAP_SEED)
    unit.update(n_cells=len(use), n_constructs=len(kept), sign_test=len(kept) > 1)

    cp_rows, det = [], np.zeros(n_genes)
    for idx, block in read_rows(X, use, chunk):
        cp_rows.append(normalise(block, idx, scale))
        det += (block > 0).sum(axis=0)
    cp = np.vstack(cp_rows)
    det /= len(use)
    mean_t = cp.mean(axis=0, dtype=np.float64)
    tested = ~((det < MIN_DETECTION) & (ctrl.detected < MIN_DETECTION))
    tested[t] = False
    p = np.full(n_genes, np.nan)
    cols = np.flatnonzero(tested)
    if len(cols):
        p[cols] = mannwhitneyu(np.log1p(cp[:, cols]), ctrl.log[:, cols], axis=0, method="asymptotic").pvalue
    q = bh(p)
    log2fc = np.log2((mean_t + PSEUDOCOUNT) / (ctrl.mean + PSEUDOCOUNT))
    hit = tested & (q < FDR) & (np.abs(log2fc) >= MIN_ABS_LOG2FC)
    if len(kept) > 1:
        signs = np.array([np.sign(np.log2((sums[c] / counts[c] + PSEUDOCOUNT) / (ctrl.mean + PSEUDOCOUNT)))
                          for c in kept])
        hit &= (signs == signs[0]).all(axis=0) & (signs[0] != 0)
    g = np.flatnonzero(hit)
    unit["n_responders"] = len(g)
    out = pd.DataFrame({"tf": tf, "gene": np.asarray(var_names).astype(str)[g], "log2fc": log2fc[g], "q": q[g],
                        "n_cells": len(use), "kd": unit["knockdown"]})
    if len(g) < MIN_RESPONDERS:
        unit["exclusion_reason"] = "too_few_responders"
        return unit, empty
    unit["usable"] = True
    return unit, out


_STATE = {}


def _init_worker(h5ad, split_seed, chunk, perturbation_col, control_label, norm="cp10k"):
    import anndata as ad
    a = ad.read_h5ad(h5ad, backed="r")
    obs = a.obs
    scale, good = core_scale(obs, norm)
    control = control_cells(obs, perturbation_col, control_label, norm)
    _, truth_half = split_controls(obs.index, np.flatnonzero(control), split_seed)
    _STATE.update(a=a, X=a.X, var_names=a.var_names, lookup=gene_lookup(a.var_names), chunk=chunk, scale=scale,
                  ctrl=Controls(a.X, truth_half, chunk, scale))


def _work(job):
    tf, cells, constructs = job
    s = _STATE
    return tf_unit(s["X"], s["var_names"], s["lookup"], s["ctrl"], tf, cells, constructs, s["chunk"], s["scale"])


def cmd_truth(args) -> int:
    import anndata as ad
    t0 = time.time()
    a = ad.read_h5ad(args.h5ad, backed="r")
    try:
        obs, col = a.obs, args.perturbation_col
        if col not in obs.columns:
            raise SystemExit(f"obs has no column {col!r}; pass --perturbation-col")
        label = obs[col].astype(str)
        scale, good = core_scale(obs, args.norm)
        ok = usable_cells(obs, col)
        lab = label.to_numpy()
        tfs = sorted(pd.read_csv(args.tfs, sep="\t", dtype=str)["source"].dropna().unique())
        used = ok & ((lab == args.control_label) | np.isin(lab, tfs))
        n_no_factor = int((used & ~good).sum())
        ok &= good
        control = control_cells(obs, col, args.control_label, args.norm)
        if control.sum() < 2 * N_TRUTH_CONTROLS:
            print(f"warning: {int(control.sum())} controls, fewer than two halves of {N_TRUTH_CONTROLS}",
                  file=sys.stderr)
        inf_half, truth_half = split_controls(obs.index, np.flatnonzero(control), args.seed)
        construct = obs[args.construct_col].astype(str).to_numpy() if args.construct_col in obs.columns else label.to_numpy()
        wanted = set(tfs)
        by_tf = {}
        for i in np.flatnonzero(ok & np.isin(lab, list(wanted))):
            by_tf.setdefault(lab[i], []).append(i)
        jobs = [(tf, np.array(by_tf[tf]), construct[by_tf[tf]]) for tf in tfs if tf in by_tf]
        n_var, n_obs = a.n_vars, a.n_obs
        initargs = (args.h5ad, args.seed, args.chunk, col, args.control_label, args.norm)
        truth_names = sorted(map(str, obs.index[truth_half]))
    finally:
        a.file.close()

    if args.workers > 1:
        with Pool(args.workers, _init_worker, initargs) as pool:
            results = pool.map(_work, jobs, chunksize=1)
    else:
        _init_worker(*initargs)
        results = []
        for k, job in enumerate(jobs, 1):
            results.append(_work(job))
            print(f"{k}/{len(jobs)} {job[0]} usable={results[-1][0]['usable']} {time.time() - t0:.0f}s", flush=True)
        _STATE["a"].file.close()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "truth_half.txt").write_text("\n".join(truth_names) + "\n", encoding="utf-8")
    cols = ["tf", "n_cells", "n_constructs", "knockdown", "n_responders", "usable", "exclusion_reason",
            "n_cells_total", "n_constructs_total", "sign_test"]
    units = pd.DataFrame([u for u, _ in results], columns=cols)
    resp = pd.concat([r for _, r in results if len(r)], ignore_index=True) if any(len(r) for _, r in results) else results[0][1] if results else pd.DataFrame()
    resp = resp[resp.tf.isin(units.tf[units.usable])] if len(resp) else resp
    units.to_csv(out / "units.tsv", sep="\t", index=False, float_format="%.6g")
    resp.to_csv(out / "responders.tsv", sep="\t", index=False, float_format="%.6g")
    digest = lambda p: hashlib.sha256((out / p).read_bytes()).hexdigest()
    meta = {
        "h5ad": Path(args.h5ad).name, "tfs_file": Path(args.tfs).name, "n_cells": n_obs, "n_genes": n_var,
        "perturbation_col": col, "control_label": args.control_label, "construct_col": args.construct_col,
        "constructs_present": args.construct_col in obs.columns,
        "normalisation": NORMALISATION if args.norm == "cp10k" else NORMALISATION_CORE,
        "control_split": {"method": "barcodes sorted, permuted with np.random.default_rng(seed), first half "
                                    "= inference, second half = truth", "seed": args.seed,
                          "n_controls": int(control.sum()), "n_inference": len(inf_half), "n_truth": len(truth_half)},
        "seeds": {"split": args.seed, "truth_controls": CONTROL_SEED, "tf_cell_cap": CELL_CAP_SEED},
        "rule": {"n_truth_controls": N_TRUTH_CONTROLS, "max_cells_per_tf": MAX_CELLS_PER_TF, "min_kd": MIN_KD,
                 "min_cells_per_construct": MIN_CELLS_PER_CONSTRUCT, "min_detection": MIN_DETECTION,
                 "fdr": FDR, "min_abs_log2fc": MIN_ABS_LOG2FC, "pseudocount": PSEUDOCOUNT,
                 "min_responders": MIN_RESPONDERS,
                 "test": "scipy.stats.mannwhitneyu asymptotic (continuity, tie-corrected) per gene on log1p "
                         + ("CP10k" if args.norm == "cp10k" else f"counts / {CORE_COL}") + ", BH within TF over tested genes",
                 "sign_agreement": "per-construct pseudobulk log2FC (all cells of kept constructs), same sign in "
                                   "every kept construct; only tested when >= 2 kept constructs (units.tsv sign_test)"},
        "n_tfs_in_reference": len(tfs), "n_tfs_perturbed": len(jobs),
        **({"norm": "core", "norm_obs_column": CORE_COL,
            "n_cells_excluded_no_factor": n_no_factor} if args.norm == "core" else {}),
        "sha256": {"units.tsv": digest("units.tsv"), "responders.tsv": digest("responders.tsv")},
        "seconds": round(time.time() - t0, 1),
    }
    (out / "truth.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print(f"{len(jobs)} perturbed TFs, {int(units.usable.sum())} usable, {len(resp)} responder rows -> {out}")
    return 0


# ---------------------------------------------------------------- subset
def subset_blocks(X, keep: np.ndarray, block: int):
    """Yield the kept rows as CSR float32, one full-width read of ``block`` rows of X at a time, in file order."""
    for start in range(0, len(keep), block):
        rows = np.flatnonzero(keep[start:start + block])
        if len(rows):
            b = X[start:start + block, :]
            yield sp.csr_matrix((b if sp.issparse(b) else np.asarray(b))[rows], dtype=np.float32)


def cmd_subset(args) -> int:
    import anndata as ad
    import h5py
    try:
        from anndata.io import sparse_dataset, write_elem
    except ImportError:             # anndata 0.10
        from anndata.experimental import sparse_dataset, write_elem
    t0 = time.time()
    a = ad.read_h5ad(args.h5ad, backed="r")
    try:
        obs, col = a.obs, args.perturbation_col
        if col not in obs.columns:
            raise SystemExit(f"obs has no column {col!r}; pass --perturbation-col")
        lab = obs[col].astype(str).to_numpy()
        ok = usable_cells(obs, col)
        tfs = sorted(pd.read_csv(args.tfs, sep="\t", dtype=str)["source"].dropna().unique())
        keep = ok & ((lab == args.control_label) | np.isin(lab, tfs))
        chunks = getattr(a.X, "chunks", None)
        block = args.block_rows or (chunks[0] if chunks else CHUNK)
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        sub = ad.AnnData(obs=obs[keep].copy(), var=a.var.copy())
        sub.write_h5ad(out)
        n_kept, nnz = int(keep.sum()), 0
        with h5py.File(out, "a") as f:
            for csr in subset_blocks(a.X, keep, block):
                if "X" in f:
                    sparse_dataset(f["X"]).append(csr)
                else:
                    write_elem(f, "X", csr)
                nnz += csr.nnz
        n_in, n_var = a.n_obs, a.n_vars
    finally:
        a.file.close()
    size = out.stat().st_size
    meta = {"input": str(args.h5ad), "input_bytes": Path(args.h5ad).stat().st_size, "n_rows_in": n_in,
            "n_rows_out": n_kept, "n_genes": n_var, "nnz": nnz, "block_rows": block, "x_chunks": list(chunks or []),
            "tfs_file": Path(args.tfs).name, "perturbation_col": col, "control_label": args.control_label,
            "output_bytes": size, "seconds": round(time.time() - t0, 1),
            "selection": "rows with usable_cells (nperts <= 1) whose perturbation is the control label or a "
                         "CollecTRI source; original row order; all obs columns and the obs index kept"}
    out.with_suffix(".subset.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print(f"{n_kept}/{n_in} rows kept, nnz {nnz}, in-memory CSR about {(8 * nnz + 8 * n_kept) / 1e9:.2f} GB, "
          f"{meta['seconds']}s, {size / 1e9:.2f} GB -> {out}")
    return 0


# ---------------------------------------------------------------- total shift (check P9 sensitivity)
def obs_column(g, name: str) -> np.ndarray:
    """One obs column straight from an h5ad's obs group (categorical or plain), without loading X."""
    o = g[name]
    v = o["categories"][()][o["codes"][()]] if hasattr(o, "keys") else o[()]
    return np.array([x.decode() if isinstance(x, bytes) else x for x in v]) if v.dtype.kind in "OS" else v


def cmd_total_shift(args) -> int:
    """Per usable TF: median(UMI_count / core_scale_factor) of its cells over the same median of the controls."""
    import h5py
    with h5py.File(args.h5ad, "r") as f:
        obs = f["obs"]
        missing = [c for c in (args.perturbation_col, "UMI_count", CORE_COL) if c not in obs]
        if missing:
            raise SystemExit(f"obs has no column {missing}; total-shift needs {args.perturbation_col!r}, 'UMI_count' and {CORE_COL!r}")
        pert = obs_column(obs, args.perturbation_col).astype(str)
        adj = obs_column(obs, "UMI_count").astype(float) / obs_column(obs, CORE_COL).astype(float)
    if not (pert == args.control_label).any():
        raise SystemExit(f"no cell with {args.perturbation_col} == {args.control_label!r}")
    ctrl = np.median(adj[pert == args.control_label])
    units = pd.read_csv(Path(args.truth) / "units.tsv", sep="\t", keep_default_na=False)
    rows = []
    for tf in sorted(units.tf[units.usable.astype(str) == "True"].astype(str)):
        m = pert == tf
        rel = float(np.median(adj[m]) / ctrl) if m.any() else float("nan")
        rows.append((tf, int(m.sum()), rel, bool(0.9 <= rel <= 1.1)))
    pd.DataFrame(rows, columns=["tf", "n_cells", "rel_core_adjusted_total", "within_0_9_1_1"]).to_csv(
        args.out, sep="\t", index=False, float_format="%.6g")
    print(f"{len(rows)} usable TFs, {sum(r[3] for r in rows)} within [0.9, 1.1] -> {args.out}")
    return 0


# ---------------------------------------------------------------- units
def summarise(units: pd.DataFrame) -> dict:
    usable = units[units.usable.astype(str) == "True"]
    q = np.percentile(usable.n_responders, [0, 25, 50, 75, 100]) if len(usable) else [None] * 5
    reasons = units.loc[units.usable.astype(str) != "True", "exclusion_reason"].value_counts()
    return {"n_perturbed_tfs": int(len(units)), "n_usable": int(len(usable)),
            "excluded_by_reason": {k: int(v) for k, v in reasons.items()},
            "responders_among_usable": dict(zip(["min", "p25", "median", "p75", "max"],
                                                [None if v is None else float(v) for v in q]))}


def cmd_units(args) -> int:
    print(json.dumps(summarise(pd.read_csv(Path(args.truth) / "units.tsv", sep="\t", keep_default_na=False)), indent=2))
    return 0


# ---------------------------------------------------------------- controls and baseline
def cmd_controls(args) -> int:
    import anndata as ad
    a = ad.read_h5ad(args.h5ad, backed="r")
    try:
        obs, col = a.obs, args.perturbation_col
        control = control_cells(obs, col, args.control_label, args.norm)
        inf_half, _ = split_controls(obs.index, np.flatnonzero(control), args.seed)
        pick = cap(inf_half, args.n, SPLIT_SEED) if args.n else inf_half
        X = np.vstack([b for _, b in read_rows(a.X, pick, args.chunk)])
        out = ad.AnnData(X=X, obs=pd.DataFrame(index=obs.index[pick].astype(str)), var=a.var.copy())
        out.uns["split"] = {"seed": args.seed, "n_inference_half": int(len(inf_half)), "cap_seed": SPLIT_SEED}
        out.uns["source_h5ad"] = Path(args.h5ad).name
    finally:
        a.file.close()
    Path(args.out).mkdir(parents=True, exist_ok=True)
    path = Path(args.out) / "controls.h5ad"
    out.write_h5ad(path)
    print(f"{out.n_obs} inference-half control cells x {out.n_vars} genes -> {path}")
    return 0


def baseline_check(h5ad, workdir=None) -> dict:
    """Run the fixed baseline (D-021) with hvg = all_genes; report genes/cells kept, seconds, error."""
    import tempfile

    import anndata as ad
    from grnbench.autoprep import diagnose, execute_plan, plan_preprocessing
    from grnbench.profile import Profile
    from grnbench.screen import BASELINE_CHOICES
    from grnbench.selector import select

    t0 = time.time()
    res = {"genes_in": None, "genes_kept": None, "cells_kept": None, "seconds": None, "error": None}
    try:
        diag = diagnose(h5ad)
        res["genes_in"] = int(diag.adata.n_vars)
        selection = select(Profile(**plan_preprocessing(diag).profile))
        plan = plan_preprocessing(diag, selection=selection,
                                  choices={**selection.chosen(), **BASELINE_CHOICES, "hvg": "all_genes"})
        work = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="kd_baseline_"))
        records = execute_plan(plan, h5ad, work)
        for r in records:
            if r["error"]:
                res["error"] = f"{r['context']}: {r['error']}"
            elif r["processed"]:
                p = ad.read_h5ad(r["processed"], backed="r")
                res.update(genes_kept=int(p.n_vars), cells_kept=int(p.n_obs))
                p.file.close()
        if res["genes_kept"] is None and not res["error"]:
            res["error"] = "no processed output"
    except Exception as exc:        # the check reports the failure, it does not hide it
        res["error"] = f"{type(exc).__name__}: {exc}"
    res["seconds"] = round(time.time() - t0, 1)
    return res


def cmd_baseline(args) -> int:
    res = baseline_check(args.h5ad, args.workdir)
    print(json.dumps(res, indent=2))
    return 1 if res["error"] else 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    def screen_args(s):
        s.add_argument("--h5ad", required=True)
        s.add_argument("--out", required=True)
        s.add_argument("--seed", type=int, default=SPLIT_SEED, help="control split seed (record: 0)")
        s.add_argument("--chunk", type=int, default=CHUNK, help="rows per read of X")
        s.add_argument("--perturbation-col", default="perturbation")
        s.add_argument("--control-label", default="control")

    t = sub.add_parser("truth")
    screen_args(t)
    t.add_argument("--tfs", required=True, help="CollecTRI tsv with a 'source' column")
    t.add_argument("--construct-col", default="guide_id")
    t.add_argument("--workers", type=int, default=1)
    norm_help = f"cp10k (default) or counts / obs[{CORE_COL!r}] (a subset file keeps that column)"
    t.add_argument("--norm", choices=["cp10k", "core"], default="cp10k", help=norm_help)
    t.set_defaults(fn=cmd_truth)
    u = sub.add_parser("units")
    u.add_argument("--truth", required=True)
    u.set_defaults(fn=cmd_units)
    c = sub.add_parser("controls")
    screen_args(c)
    c.add_argument("--norm", choices=["cp10k", "core"], default="cp10k", help="same control mask as truth --norm")
    c.add_argument("--n", type=int, help="cap on the inference-half cells (default: the whole half; the record caps after QC)")
    c.set_defaults(fn=cmd_controls)
    sb = sub.add_parser("subset")
    sb.add_argument("--h5ad", required=True)
    sb.add_argument("--tfs", required=True, help="CollecTRI tsv with a 'source' column")
    sb.add_argument("--out", required=True, help="output .h5ad (X as CSR)")
    sb.add_argument("--block-rows", type=int, help="rows per read of X (default: the dataset's chunk rows)")
    sb.add_argument("--perturbation-col", default="perturbation")
    sb.add_argument("--control-label", default="control")
    sb.set_defaults(fn=cmd_subset)
    b = sub.add_parser("baseline-check")
    b.add_argument("--h5ad", required=True)
    b.add_argument("--workdir")
    b.set_defaults(fn=cmd_baseline)
    ts = sub.add_parser("total-shift")
    ts.add_argument("--h5ad", required=True)
    ts.add_argument("--truth", required=True)
    ts.add_argument("--out", required=True, help="output .tsv")
    ts.add_argument("--perturbation-col", default="perturbation")
    ts.add_argument("--control-label", default="control")
    ts.set_defaults(fn=cmd_total_shift)
    args = p.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
