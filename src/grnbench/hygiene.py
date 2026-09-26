"""Data-hygiene utilities for the messy realities of scRNA-seq metadata.

Everything here is diagnostic-or-cosmetic and LOUD: label normalization reports
every change, ID-type/organism detection produces warnings instead of silent
TF-matching failures, and structural fixes (sparsify, downcast) are logged as
pipeline steps.
"""
from __future__ import annotations

import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

ENSEMBL_RE = re.compile(r"^ENS[A-Z]{0,4}\d+(\.\d+)?$")
HUMAN_SYMBOL_RE = re.compile(r"^[A-Z][A-Z0-9@.\-]{1,14}$")
MOUSE_SYMBOL_RE = re.compile(r"^[A-Z][a-z]{1,10}\d+[A-Za-z]?$")
BATCH_COL_RE = re.compile(r"(batch|donor|sample_id|run|lane|library|lib_|experiment|patient)",
                          re.IGNORECASE)

_LABEL_MISSING_TOKENS = {"", "nan", "none", "null", "na", "unknown", "unk", "?"}


# --------------------------------------------------------------------------- #
# labels
# --------------------------------------------------------------------------- #

def missing_label_mask(series: pd.Series) -> np.ndarray:
    """True where a label is absent: NaN or one of the usual missing tokens.

    One definition, shared by label hygiene, the obs audit and stratified
    subsampling, so a cell counted as unlabelled in the report is the same cell
    that is excluded from a stratum.
    """
    raw = series.astype("object")
    mask = raw.isna() | raw.astype(str).str.strip().str.lower().isin(_LABEL_MISSING_TOKENS)
    return mask.to_numpy(dtype=bool)


def normalize_labels(series: pd.Series,
                     mapping: dict[str, str] | None = None,
                     case: str | None = None) -> tuple[pd.Series, dict]:
    """Clean a categorical annotation column; report every transformation.

    Steps: missing-token detection -> whitespace strip/collapse -> optional
    case unification -> synonym mapping. NaN stays NaN. Returns (clean series,
    report dict with change counts and old->new samples).
    """
    raw = series.astype("object")
    is_missing = pd.Series(missing_label_mask(raw), index=raw.index)
    cleaned = raw.copy()
    cleaned[is_missing] = np.nan  # missing tokens become real NaN for downstream drops
    cleaned[~is_missing] = (
        raw[~is_missing].astype(str).str.strip().str.replace(r"\s+", " ", regex=True)
    )
    n_ws = int((cleaned[~is_missing].astype(str) != raw[~is_missing].astype(str)).sum())

    if case in {"upper", "lower", "title"}:
        before = cleaned[~is_missing].copy()
        cleaned[~is_missing] = getattr(cleaned[~is_missing].astype(str).str, case)()
        n_case = int((before != cleaned[~is_missing]).sum())
    else:
        n_case = 0

    n_mapped = 0
    if mapping:
        mapped = cleaned[~is_missing].astype(str).map(
            lambda v: mapping.get(v, mapping.get(v.lower(), v)))
        n_mapped = int((mapped != cleaned[~is_missing]).sum())
        cleaned[~is_missing] = mapped

    report = {
        "n_total": int(len(raw)),
        "n_missing": int(is_missing.sum()),
        "missing_tokens_seen": sorted(set(raw[is_missing].astype(str).str.lower()) - {""})[:8],
        "n_whitespace_fixed": n_ws,
        "n_case_changed": n_case,
        "n_synonym_mapped": n_mapped,
        "n_unique_after": int(cleaned[~is_missing].nunique()),
        "unique_values_after": sorted(map(str, cleaned[~is_missing].unique()))[:25],
    }
    return cleaned, report


def load_label_mapping(path: str | Path) -> dict[str, str]:
    """Two-column file (old<TAB|comma>new or header with from/to aliases)."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    sep = "\t" if p.suffix == ".tsv" else ","
    df = pd.read_csv(p, sep=sep, dtype=str)
    cols = {c.strip().lower(): c for c in df.columns}
    src = next((cols[k] for k in ("from", "old", "original", "raw") if k in cols),
               df.columns[0])
    dst = next((cols[k] for k in ("to", "new", "clean", "standard") if k in cols),
               df.columns[-1])
    pairs = df[[src, dst]].dropna().astype(str)
    return {r.strip(): s.strip() for r, s in zip(pairs[src], pairs[dst])}


def audit_obs(adata: AnnData) -> list[dict]:
    """Per-column stats for every object/category obs column."""
    out = []
    for col in adata.obs.columns:
        s = adata.obs[col]
        if not (pd.api.types.is_object_dtype(s)
                or isinstance(s.dtype, pd.CategoricalDtype)
                or pd.api.types.is_string_dtype(s)):
            continue
        non_missing = ~missing_label_mask(s)
        counts = s[non_missing].astype(str).value_counts()
        if len(counts):
            ratio = float(counts.max() / counts.min()) if counts.min() > 0 else float("inf")
        else:
            ratio = float("nan")
        out.append({
            "column": col,
            "n_unique_nonmissing": int(len(counts)),
            "n_missing": int((~non_missing).sum()),
            "largest_class_frac": round(float(counts.max()) / max(int(non_missing.sum()), 1), 4),
            "imbalance_ratio": ratio,
            "smallest_class_cells": int(counts.min()) if len(counts) else 0,
            "batch_like_name": bool(BATCH_COL_RE.search(str(col))),
        })
    return out


def detect_batch_column(adata: AnnData, explicit: str | None = None) -> tuple[str | None, str]:
    """Explicit key wins; else first batch-like-named column with 2..200 levels."""
    if explicit:
        if explicit not in adata.obs.columns:
            raise KeyError(f"--batch-key '{explicit}' not in obs")
        return explicit, "explicit"
    for col in adata.obs.columns:
        if BATCH_COL_RE.search(str(col)):
            n = adata.obs[col].nunique()
            if 2 <= n <= 200:
                return col, "auto-detected"
    return None, "none found"


# --------------------------------------------------------------------------- #
# gene identifiers / organism
# --------------------------------------------------------------------------- #

def detect_id_type(var_names) -> dict:
    names = pd.Index(var_names).astype(str)
    ens = names.str.match(ENSEMBL_RE).mean()
    human = names.str.match(HUMAN_SYMBOL_RE).mean()
    if ens > 0.9:
        verdict = "ensembl_only"
    elif ens > 0.05:
        verdict = "mixed"
    else:
        verdict = "symbols"
    return {"ensembl_frac": round(float(ens), 4),
            "symbol_like_frac": round(float(human), 4),
            "verdict": verdict}


def organism_hint(var_names) -> str:
    names = pd.Index(var_names).astype(str)
    mouse = names.str.match(MOUSE_SYMBOL_RE).mean()
    lower_mito = names.str.startswith(("mt-", "Mt-")).any()
    human_upper = names.str.match(HUMAN_SYMBOL_RE).mean()
    upper_mito = names.str.upper().str.startswith("MT-").any()
    if mouse > 0.3 or (lower_mito and human_upper < 0.5):
        return "mouse-like"
    if human_upper > 0.7 and upper_mito:
        return "human-like"
    return "ambiguous"


def knee_cell_threshold(totals, min_counts: int = 10) -> dict | None:
    """Counts threshold separating cells from empty droplets on the barcode-rank curve.

    Log-rank vs log-total, interpolated and smoothed: the *knee* is the point
    farthest above the chord between the first and last barcode (end of the cell
    plateau); the threshold is the *inflection* after it (steepest descent),
    halfway down the cliff to the ambient plateau - the same two landmarks as
    DropletUtils' ``barcodeRanks``. Returns ``None`` when the curve has no clear
    cliff (fewer than 50 or more than 90% of barcodes would be called).
    """
    t = np.sort(np.asarray(totals, dtype=float))[::-1]
    t = t[t >= min_counts]
    if t.size < 200:
        return None
    x = np.log10(np.arange(1, t.size + 1))
    y = np.log10(t)
    xs = np.linspace(x[0], x[-1], 400)
    ys = np.interp(xs, x, y)
    k = 9
    ys = np.convolve(np.pad(ys, k // 2, mode="edge"), np.ones(k) / k, mode="valid")
    chord = ys[0] + (ys[-1] - ys[0]) * (xs - xs[0]) / (xs[-1] - xs[0])
    i_knee = int(np.argmax(ys - chord))
    d = np.gradient(ys, xs)
    i_inf = i_knee + int(np.argmin(d[i_knee:]))
    threshold = float(10 ** ys[i_inf])
    n_called = int(np.sum(np.asarray(totals) >= threshold))
    if n_called < 50 or n_called > 0.9 * np.asarray(totals).size:
        return None
    return {"method": "knee", "threshold": threshold, "knee_counts": float(10 ** ys[i_knee]),
            "n_called": n_called, "n_barcodes": int(np.asarray(totals).size)}


def duplicate_barcodes_report(adata: AnnData) -> int:
    return int(pd.Index(adata.obs_names).duplicated().sum())


# --------------------------------------------------------------------------- #
# doublets (heuristic), memory hygiene
# --------------------------------------------------------------------------- #

# Scrublet's defaults need more cells than a PCA of 30 components (hard floor
# 31) and enough simulated doublets to place its automatic threshold. Measured
# on the four EXP-20260925-LM-2 datasets that fell back (per-batch scrublet on
# 101 batches, issue #6): batches of <= 30 cells failed; 31-50 cells called a
# median 74% doublets; 51-100 cells 14% (max 33%); 101-200 cells 3.9% (max
# 10%); larger batches 1-2%. Real rates are ~1-10%, so 100 is where the
# threshold becomes usable. Pooling the batches under 100 cells gave pools of
# 206-556 cells calling 0-3%, in line with the large batches of the same data.
SCRUBLET_MIN_CELLS = 100
_POOL = "__pooled_small_batches__"


def scrublet_groups(labels, min_cells: int = SCRUBLET_MIN_CELLS) -> tuple[np.ndarray, dict]:
    """Group cells for scrublet: returns ``(group per cell or None, summary)``.

    Batches with at least ``min_cells`` cells are scored on their own (scrublet
    is meant to run per sample). Smaller batches, including missing labels as
    one level, are pooled and scored together: a heterotypic doublet still
    looks like a mix of two programs across samples, and a 13-cell batch
    cannot fit scrublet's PCA at all. A pool still under ``min_cells`` is not
    scored (None): its cells are kept and counted in ``unscored_cells``.
    """
    lab = pd.Series(labels).astype(object).where(pd.Series(labels).notna(), "<missing>")
    lab = lab.astype(str).to_numpy()
    levels, counts = np.unique(lab, return_counts=True)
    large = set(levels[counts >= min_cells])
    small = counts < min_cells
    groups = np.array([x if x in large else _POOL for x in lab], dtype=object)
    pooled_cells = int(counts[small].sum())
    unscored = 0 < pooled_cells < min_cells
    if unscored:
        groups[groups == _POOL] = None
    return groups, {"min_cells": int(min_cells), "per_batch": len(large),
                    "pooled_batches": int(small.sum()), "pooled_cells": pooled_cells,
                    "unscored_cells": pooled_cells if unscored else 0}


def detect_doublets(adata: AnnData, batch_key: str | None = None,
                    seed: int = 0) -> tuple[str, np.ndarray, str | None, dict]:
    """Score and call doublets; returns ``(method, is_doublet, warning, groups)``.

    Uses scanpy's built-in scrublet on the raw counts in ``adata.X``, per batch
    when ``batch_key`` is given, with batches under ``SCRUBLET_MIN_CELLS``
    pooled (``scrublet_groups``; ``groups`` is its summary). Writes
    ``obs['doublet_score']`` (NaN for unscored cells) and
    ``obs['predicted_doublet']``. When no group reaches the minimum, or
    scrublet fails for any other reason, the joint-extreme heuristic below is
    used instead and the caller gets a warning to record - the method that
    actually ran is never left implicit.
    """
    import scanpy as sc

    labels = adata.obs[batch_key].to_numpy() if batch_key else np.zeros(adata.n_obs)
    groups, info = scrublet_groups(labels)
    if not batch_key:   # one group of all cells: nothing was pooled
        info.update(per_batch=0, pooled_batches=0, pooled_cells=0)
    info["batch_key"] = batch_key
    scored = np.array([g is not None for g in groups])
    notes = []
    try:
        if not scored.any():
            raise ValueError(f"{adata.n_obs} cells, fewer than the {SCRUBLET_MIN_CELLS} "
                             "scrublet needs (SCRUBLET_MIN_CELLS)")
        sub = adata if scored.all() else adata[scored].copy()   # no extra copy when all scored
        sub.obs["_scrublet_group"] = groups[scored].astype(str)
        n_groups = sub.obs["_scrublet_group"].nunique()
        sc.pp.scrublet(sub, batch_key="_scrublet_group" if n_groups > 1 else None,
                       random_state=int(seed), verbose=False)
        del sub.obs["_scrublet_group"]
        if sub is not adata:
            adata.uns["scrublet"] = sub.uns["scrublet"]   # thresholds stay auditable
        score = np.full(adata.n_obs, np.nan)
        score[scored] = sub.obs["doublet_score"].to_numpy(float)
        pred = sub.obs["predicted_doublet"]
        if pred.isna().any():
            notes.append(f"scrublet could not set a doublet threshold for {int(pred.isna().sum())} "
                         "cells (automatic threshold detection failed); they are kept")
        is_dbl = np.zeros(adata.n_obs, dtype=bool)
        is_dbl[scored] = pred.fillna(False).to_numpy(dtype=bool)
        adata.obs["doublet_score"] = score
        adata.obs["predicted_doublet"] = is_dbl
        if batch_key and info["pooled_batches"]:
            what = (f"{info['pooled_batches']} batches under {SCRUBLET_MIN_CELLS} cells "
                    f"({info['pooled_cells']} cells) were pooled")
            notes.append(what + (f" and left unscored (pool under {SCRUBLET_MIN_CELLS} cells; "
                                 f"{info['unscored_cells']} cells kept without a doublet call)"
                                 if info["unscored_cells"] else " and scored together by scrublet"))
        return "scrublet", is_dbl, "; ".join(notes) or None, info
    except Exception as exc:  # noqa: BLE001 - any scrublet failure -> heuristic, loudly
        _frac, flagged = doublet_heuristic_fraction(adata)
        adata.obs["doublet_score"] = np.nan
        adata.obs["predicted_doublet"] = flagged
        hint = (" Install scikit-image (it is in requirements/lock-py311.txt): scanpy's scrublet "
                "needs it to pick the doublet threshold." if "skimage" in str(exc) else "")
        msg = (f"DOUBLET FALLBACK: scrublet failed ({type(exc).__name__}: {exc}); doublets were "
               "called with the joint top-1% genes+counts heuristic instead - a much cruder call."
               + hint)
        warnings.warn(msg, RuntimeWarning, stacklevel=2)
        return "heuristic", flagged, msg, {**info, "unscored_cells": int(adata.n_obs)}


def doublet_heuristic_fraction(adata: AnnData,
                               percentile: float = 99.0) -> tuple[float, np.ndarray]:
    """Fraction of cells jointly extreme in detected genes AND total counts.

    A cheap multiplet proxy - not scrublet; flagged cells are reported, never
    dropped without an explicit user decision.
    """
    X = adata.X
    total = np.asarray(X.sum(axis=1)).ravel()
    n_genes = X.getnnz(axis=1) if sp.issparse(X) else (np.asarray(X) > 0).sum(axis=1)
    t_cut = np.percentile(total, percentile)
    g_cut = np.percentile(n_genes, percentile)
    flagged = (total >= t_cut) & (n_genes >= g_cut)
    return float(flagged.mean()), flagged


def sparsify_if_dense(adata: AnnData, zero_threshold: float = 0.8) -> bool:
    if sp.issparse(adata.X):
        return False
    zero_frac = 1.0 - np.count_nonzero(adata.X) / adata.X.size
    if zero_frac >= zero_threshold:
        adata.X = sp.csr_matrix(adata.X)
        for k, v in list(adata.layers.items()):
            if not sp.issparse(v) and v.ndim == 2:
                adata.layers[k] = sp.csr_matrix(v)
        return True
    return False


def downcast_floats(adata: AnnData) -> bool:
    """float64 -> float32 on X and every layer; returns whether anything moved."""
    changed = False

    def _cast(m):
        nonlocal changed
        if sp.issparse(m):
            if m.dtype == np.float64:
                changed = True
                return m.astype(np.float32)
            return m
        if isinstance(m, np.ndarray) and m.dtype == np.float64:
            changed = True
            return m.astype(np.float32)
        return m

    adata.X = _cast(adata.X)
    for k in list(adata.layers.keys()):
        adata.layers[k] = _cast(adata.layers[k])
    return changed
