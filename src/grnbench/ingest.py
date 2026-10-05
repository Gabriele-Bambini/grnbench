"""Dataset ingestion: common scRNA-seq formats -> counts-backed AnnData.

Contract: after :func:`load_any`, ``adata.X`` holds **raw counts** (original
counts are preserved in ``adata.layers["counts"]``). Anything that does not
look like counts is rejected loudly — silent normalization guessing is how
benchmarks rot — with advice on what the data look like and where counts may be.
Integer input that still looks non-raw (near-identical cell totals, 0/1 only) is
warned about and recorded in the report (:func:`count_plausibility`, issue #42).
Gene-name hygiene (version suffixes, duplicates aggregated by sum) is applied and
logged in ``adata.uns["grnbench_ingest"]``.
"""
from __future__ import annotations

import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

ENSEMBL_VERSION_RE = re.compile(r"^(ENSG\d+)\.\d+$")
ENSEMBL_ANY_RE = re.compile(r"^ENS[A-Z]{0,4}G\d+(\.\d+)?$")
# var columns that carry gene symbols next to Ensembl IDs (CELLxGENE, 10x, scanpy)
SYMBOL_COLUMNS = ("feature_name", "gene_symbols", "gene_symbol", "gene_name", "gene_names",
                  "symbol", "Symbol", "SYMBOL", "features")
GENE_SYMBOL_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]{0,19}$")
# 10x-style cell barcodes: an optional sample prefix, 8+ ACGT bases, optional
# lane suffix. Deliberately narrow - this is used to OVERRIDE gene-likeness,
# so a false positive here would flip a matrix that was already correct.
BARCODE_RE = re.compile(r"^(?:[A-Za-z0-9]+[:_])?[ACGTNacgtn]{8,}(?:[-_]\d+)?$")
MITO_PREFIXES = ("MT-", "MT.", "Mt-", "mt-")
RIBO_PREFIXES = ("RPS", "RPL", "rps", "rpl")

# Count plausibility (issue #42). An integer, non-negative matrix is not always raw
# counts: rounded normalize_total/CPM output and binarized 0/1 matrices pass
# looks_like_counts. PROVISIONAL thresholds, not yet calibrated on the approved real
# datasets (issue #42: benchmarks/count_plausibility.py over the data root must show
# zero false alarms first). Until then the check only warns and records; it never refuses.
PLAUSIBILITY_MIN_CELLS = 20        # fewer cells with counts: totals signal not evaluated
PLAUSIBILITY_MIN_NONZERO = 1000    # fewer non-zero values: value signals not evaluated
TOTALS_CV_NEAR_IDENTICAL = 0.05    # CV of per-cell totals below this: near-identical depth (warns alone)
TOTALS_CV_TIGHT = 0.15             # CV below this: tight totals (warns only together with "no 1s")
NO_ONES_MAX_SHARE = 0.001          # share of 1s among non-zero values at or below this: "no 1s"
# Refusal advice only (they describe, never decide): non-integer values whose per-cell
# totals agree this closely are normalize_total/CPM output; non-integer values never above
# LOG_RANGE_MAX are in the log range (log1p of 1e4-normalized data peaks near 9, log2 CPM near 20).
NORMALIZED_TOTALS_CV = 0.01
LOG_RANGE_MAX = 20.0


class IngestError(ValueError):
    """Raised when input data cannot be safely interpreted as raw counts."""


@dataclass
class IngestReport:
    path: str
    fmt: str
    n_duplicates_merged: int = 0
    version_suffixes_stripped: int = 0
    transposed: bool = False
    counts_source: str = "X"  # "X" or "layers['counts']"
    notes: list[str] = field(default_factory=list)
    n_case_duplicates_merged: int = 0
    # signals from count_plausibility() on the counts actually used (issue #42)
    count_plausibility: dict = field(default_factory=dict)
    # loud-but-not-fatal findings; run_prep copies them into manifest["warnings"]
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return self.__dict__.copy()


def load_any(path: str | Path,
             orientation: str = "auto") -> tuple[ad.AnnData, IngestReport]:
    """Load ``path`` into a counts-backed AnnData.

    Supported: ``.h5ad``, 10x ``.h5``, 10x MTX directory, ``.loom``,
    ``.csv``/``.tsv``/``.txt`` expression tables.

    ``orientation`` applies to the table formats only: ``"auto"`` detects
    whether rows are cells or genes and REFUSES to guess when the evidence is
    ambiguous; pass ``"genes_x_cells"`` (rows are genes — the usual GEO and
    BEELINE layout) or ``"cells_x_genes"`` to state it outright.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    if p.is_dir():
        adata, rep = _load_10x_mtx_dir(p)
    else:
        sfx = p.suffix.lower()
        if sfx == ".h5ad":
            adata, rep = _load_h5ad(p)
        elif sfx == ".h5":
            # NO var_names_make_unique() here: CellRanger references do contain
            # duplicated symbols (e.g. TBCE), and renaming them to SYMBOL-1
            # would make _clean_gene_names' aggregation branch dead - the same
            # data would then yield a different matrix than the MTX path.
            adata = sc.read_10x_h5(p)
            rep = IngestReport(path=str(p), fmt="10x_h5")
        elif sfx == ".loom":
            adata = sc.read_loom(p)
            rep = IngestReport(path=str(p), fmt="loom")
        elif sfx in {".csv", ".tsv", ".txt"}:
            adata, rep = _load_table(p, orientation)
        else:
            raise IngestError(
                f"Unsupported format '{sfx}'. Supported: .h5ad, 10x .h5, "
                "10x MTX directory, .loom, .csv/.tsv/.txt"
            )
    _ensure_counts_backed(adata, rep)
    _assert_finite(adata.X)
    _check_count_plausibility(adata.X, rep)
    adata = _clean_gene_names(adata, rep)
    return adata, rep


# --------------------------------------------------------------------------- #
# loaders
# --------------------------------------------------------------------------- #

def _load_h5ad(p: Path) -> tuple[ad.AnnData, IngestReport]:
    adata = sc.read_h5ad(p)
    rep = IngestReport(path=str(p), fmt="h5ad")
    if "counts" in adata.layers and adata.layers["counts"] is not None:
        layer = adata.layers["counts"]
        _validate_counts_layer(adata.X, layer, adata=adata)
        rep.counts_source = "layers['counts']"
        rep.notes.append(
            "X was replaced by validated layers['counts']; the previous X is "
            "preserved in layers['x_as_loaded']."
        )
        adata.layers["x_as_loaded"] = adata.X.copy()
        adata.X = layer.copy()
    elif adata.raw is not None and not looks_like_counts(adata.X)[0]:
        # CELLxGENE and many scanpy exports keep normalized values in X and the
        # raw counts in .raw - use them rather than refusing the file.
        ok, why = looks_like_counts(adata.raw.X)
        if ok:
            raw = adata.raw
            adata = ad.AnnData(X=raw.X.copy(), obs=adata.obs.copy(), var=raw.var.copy(),
                               uns=dict(adata.uns), obsm=dict(adata.obsm))
            rep.counts_source = "raw.X"
            rep.notes.append(f"X is not raw counts; the counts in .raw ({adata.n_vars} genes) were used.")
    return adata, rep


def _validate_counts_layer(x_as_loaded, layer, adata: ad.AnnData | None = None) -> None:
    """Fail loudly when the trusted counts layer is itself broken.

    ``adata`` (optional) lets the refusal list the file's other count-like
    matrices as suggestions; without it only ``x_as_loaded`` is considered.
    """
    if layer.shape != x_as_loaded.shape:
        if layer.shape == x_as_loaded.shape[::-1]:
            raise IngestError(
                "layers['counts'] appears TRANSPOSED relative to X "
                f"({layer.shape} vs {x_as_loaded.shape}); fix the writer - "
                "refusing to guess."
            )
        raise IngestError(
            f"layers['counts'] shape {layer.shape} != X shape {x_as_loaded.shape}"
        )
    ok, why = looks_like_counts(layer)
    if not ok:
        if adata is not None:
            others = _count_like_sources(x=adata.X, layers=adata.layers, raw=adata.raw,
                                         skip_layer="counts")
        else:
            others = _count_like_sources(x=x_as_loaded)
        raise IngestError(
            f"layers['counts'] does not look like raw counts: {why}. "
            + _refusal_advice(layer, *others,
                              fix="Fix or drop layers['counts'] and re-export the dataset "
                                  "with raw counts")
        )
    _assert_finite(layer)


def _assert_finite(x) -> None:
    data = x.data if sp.issparse(x) else np.asarray(x).ravel()
    if data.size and not np.all(np.isfinite(data)):
        n_bad = int(np.sum(~np.isfinite(data)))
        raise IngestError(
            f"matrix contains {n_bad} NaN/Inf entries; clean or impute-free "
            "filter them upstream - GRN preprocessing refuses to guess."
        )


def _load_10x_mtx_dir(p: Path) -> tuple[ad.AnnData, IngestReport]:
    adata = sc.read_10x_mtx(p, var_names="gene_symbols", make_unique=False)
    rep = IngestReport(path=str(p), fmt="10x_mtx")
    return adata, rep


def _load_table(p: Path, orientation: str = "auto") -> tuple[ad.AnnData, IngestReport]:
    sep = "\t" if p.suffix.lower() in {".tsv", ".txt"} else ","
    df = pd.read_csv(p, sep=sep, index_col=0)
    try:
        df = df.apply(pd.to_numeric, errors="raise")
    except (ValueError, TypeError) as e:
        raise IngestError(f"Table {p.name} is not fully numeric: {e}") from e

    transposed, why = _decide_orientation(df.index.astype(str),
                                          df.columns.astype(str),
                                          orientation, p.name)
    if transposed:
        df = df.T  # rows were genes -> transpose so obs=cells, var=genes
    rep = IngestReport(path=str(p), fmt="table", transposed=transposed)
    rep.notes.append(f"orientation: {why}")
    return ad.AnnData(X=df.to_numpy(dtype=np.float32), obs=pd.DataFrame(index=df.index.astype(str)),
                      var=pd.DataFrame(index=df.columns.astype(str))), rep


def _decide_orientation(rows, cols, orientation: str,
                        name: str) -> tuple[bool, str]:
    """Are the rows genes (transpose needed) or cells?

    Gene-likeness alone cannot answer this. A 10x barcode such as
    ``TGGCCAAAATGTGGTG-1`` matches ``GENE_SYMBOL_RE`` perfectly, so on the
    standard GEO layout (genes as rows, barcodes as columns) both axes scored
    ~1.0, the 0.1 margin was never cleared, and a 4000x600 matrix was ingested
    as 4000 cells x 600 genes — reporting ``transposed: False`` as though that
    were a determination. Every QC metric was then computed on the wrong axis:
    median genes/cell 331 instead of 2204, and ``n_mito_genes`` 0 instead of 2,
    which silently disables the mitochondrial gate.

    So decide on the signal that actually separates the two: barcodes are long
    runs of ACGT and gene symbols are not. Only fall back to gene-likeness when
    neither axis looks like barcodes, and refuse to guess when nothing decides.
    """
    if orientation not in {"auto", "genes_x_cells", "cells_x_genes"}:
        raise IngestError(
            f"orientation must be 'auto', 'genes_x_cells' or 'cells_x_genes'; "
            f"got {orientation!r}")
    if orientation == "genes_x_cells":
        return True, "declared genes x cells by the caller"
    if orientation == "cells_x_genes":
        return False, "declared cells x genes by the caller"

    row_bc, col_bc = _barcode_likeness(rows), _barcode_likeness(cols)
    row_gene, col_gene = _gene_likeness(rows), _gene_likeness(cols)
    scores = (f"rows[barcode {row_bc:.2f}, gene {row_gene:.2f}] "
              f"cols[barcode {col_bc:.2f}, gene {col_gene:.2f}]")

    # 1. barcodes are unmistakable; whichever axis holds them is the cell axis
    if col_bc > 0.6 and row_bc < 0.2:
        return True, f"columns are cell barcodes; {scores}"
    if row_bc > 0.6 and col_bc < 0.2:
        return False, f"rows are cell barcodes; {scores}"

    # 2. a real transcriptome axis carries MT-/RPS/RPL genes; a barcode axis
    #    never does, and this holds even when both axes look "symbol-like"
    row_mk, col_mk = _marker_fraction(rows), _marker_fraction(cols)
    if row_mk >= 0.001 and col_mk == 0.0:
        return True, f"rows carry mito/ribo gene markers; {scores}"
    if col_mk >= 0.001 and row_mk == 0.0:
        return False, f"columns carry mito/ribo gene markers; {scores}"

    # 3. fall back to the old margin, but only when it is genuinely decisive
    if row_gene > col_gene + 0.1:
        return True, f"rows score higher on gene-likeness; {scores}"
    if col_gene > row_gene + 0.1:
        return False, f"columns score higher on gene-likeness; {scores}"

    raise IngestError(
        f"{name}: cannot tell whether rows are cells or genes ({scores}). "
        "Guessing here silently transposes every QC metric, so it is refused. "
        "Pass the layout explicitly: orientation='genes_x_cells' (rows are "
        "genes, the usual GEO/BEELINE layout) or 'cells_x_genes'."
    )


def _gene_likeness(names) -> float:
    """Score how much a name collection looks like gene symbols (0..1)."""
    names = list(names)
    if not names:
        return 0.0
    sym = sum(bool(GENE_SYMBOL_RE.match(n)) for n in names) / len(names)
    numeric = sum(n.replace(".", "", 1).isdigit() for n in names) / len(names)
    return max(0.0, sym - numeric) + 0.5 * min(_marker_fraction(names), 0.1)


def _marker_fraction(names) -> float:
    """Fraction of names that are mitochondrial or ribosomal genes."""
    names = list(names)
    if not names:
        return 0.0
    return sum(n.upper().startswith(MITO_PREFIXES + RIBO_PREFIXES)
               for n in names) / len(names)


def _barcode_likeness(names) -> float:
    """Fraction of names shaped like cell barcodes rather than gene symbols.

    Matches the common forms: bare ``ACGT``-runs of 8+ bases, an optional
    ``-1``/``_1`` lane suffix, and an optional ``sample:``/``sample_`` prefix.
    Gene symbols essentially never look like this — the handful that are pure
    ACGT letters (none reach 8 characters) cannot reach the 0.6 threshold.
    """
    names = list(names)
    if not names:
        return 0.0
    return sum(bool(BARCODE_RE.match(str(n))) for n in names) / len(names)


# --------------------------------------------------------------------------- #
# counts integrity + hygiene
# --------------------------------------------------------------------------- #

def looks_like_counts(X) -> tuple[bool, str]:
    """Heuristic: integer-valued, non-negative matrix (dense or sparse)."""
    data = X.data if sp.issparse(X) else np.asarray(X).ravel()
    if data.size == 0:
        return False, "empty matrix"
    if np.any(data < 0):
        return False, "negative values (scaled/centered data?)"
    frac = np.mean(np.abs(data - np.round(data)) > 1e-6)
    if frac > 1e-4:
        return False, f"{frac:.1%} non-integer values (normalized/log data?)"
    return True, "integer, non-negative"


def _ensure_counts_backed(adata: ad.AnnData, rep: IngestReport) -> None:
    ok, why = looks_like_counts(adata.X)
    if ok:
        adata.layers["counts"] = adata.X.copy()
        return
    others = _count_like_sources(layers=adata.layers, raw=adata.raw)
    raise IngestError(
        "Input does not look like raw counts: " + why + ". "
        + _refusal_advice(adata.X, *others,
                          fix="GRN preprocessing must start from counts: re-export the "
                              "dataset with raw counts (e.g. layers['counts'] in the h5ad)")
    )


def _row_totals(X) -> np.ndarray:
    """Per-cell totals in float64 (float32 sums drift on deep cells); one pass, no copy."""
    if sp.issparse(X):
        return np.asarray(X.sum(axis=1, dtype=np.float64)).ravel()
    return np.asarray(X).sum(axis=1, dtype=np.float64)


def _cv(values: np.ndarray) -> float:
    return float(values.std() / values.mean())


def describe_non_counts(X) -> str:
    """Best guess, for a refusal message, of what a matrix that is not counts holds."""
    data = X.data if sp.issparse(X) else np.asarray(X).ravel()
    if data.size == 0:
        return "an empty matrix"
    if np.any(data < 0):
        return ("scaled or centered values (negative entries, as sc.pp.scale, regress_out or "
                "residuals produce); counts cannot be recovered from them")
    totals = _row_totals(X)
    totals = totals[totals > 0]
    vmax = float(data.max())
    if totals.size >= 2 and _cv(totals) < NORMALIZED_TOTALS_CV:
        return (f"per-cell normalized values without log: every cell totals about "
                f"{float(np.median(totals)):.4g} (CV {_cv(totals):.2g}), as normalize_total or "
                "CPM produce")
    if vmax <= LOG_RANGE_MAX:
        return (f"log-normalized values (non-integer, max {vmax:.3g}), such as "
                "log1p(normalize_total)")
    return (f"non-integer values up to {vmax:.4g}: normalized (TPM/FPKM), expected counts "
            "(RSEM, salmon) or imputed values")


def _count_like_sources(x=None, layers=None, raw=None,
                        skip_layer: str | None = None) -> tuple[list[str], str | None]:
    """Other matrices in the file that pass :func:`looks_like_counts`, plus a note on ``.raw``.

    Suggestions for the refusal message only: never used automatically, because an
    integer layer need not hold total counts (velocyto's 'spliced' does not).
    """
    found: list[str] = []
    if x is not None and looks_like_counts(x)[0]:
        found.append("X")
    for key in (layers.keys() if layers is not None else ()):
        if key is None or key == skip_layer:  # anndata>=0.13 aliases X as layers[None]
            continue
        if looks_like_counts(layers[key])[0]:
            found.append(f"layers['{key}']")
    raw_note = None
    if raw is not None:
        ok, why = looks_like_counts(raw.X)
        if ok:
            found.append(".raw.X")
        else:
            raw_note = f".raw is present but is not counts either ({why})."
    return found, raw_note


def _refusal_advice(X, found: list[str], raw_note: str | None, *, fix: str) -> str:
    """The readable part of a not-counts refusal: what X looks like, candidates, where to look."""
    parts = [f"It looks like {describe_non_counts(X)}."]
    if found:
        parts.append(
            "Integer, non-negative matrices in this file: " + ", ".join(found) + ". None is "
            "used automatically: check that one holds total raw counts (velocyto "
            "'spliced'/'unspliced' do not), store it as layers['counts'] and load again.")
    if raw_note:
        parts.append(raw_note)
    parts.append(
        fix + ". Raw counts are usually in layers['counts'], .raw.X, the 10x "
        "filtered_feature_bc_matrix (.h5 or MTX folder) or a counts file in the GEO "
        "supplementary data.")
    return " ".join(parts)


def count_plausibility(X) -> dict:
    """Signals that an integer, non-negative matrix is still not raw counts (issue #42).

    :func:`looks_like_counts` proves integrality only, and rounded normalize_total/CPM
    output or a binarized 0/1 matrix passes it. Three signals:

    - near-identical per-cell totals, the strongest: sequencing depth always varies
      between real cells, and per-cell normalization removes exactly that variation;
    - only 0/1 values: a detected/not-detected matrix;
    - no 1s among the non-zero values: rounding scaled counts skips 1. Weak on its own
      (full-length protocols such as Smart-seq have few 1s), so it only warns together
      with tight totals (normalized, a gene subset taken, then rounded).

    Cost: one row sum and one pass over the stored values; sparse input is never
    densified or copied. Returns the signals, the ``flags`` that fired, whether they
    ``warn`` and the provisional thresholds used (recorded in the manifest, D-005).
    """
    data = X.data if sp.issparse(X) else np.asarray(X).ravel()
    totals = _row_totals(X)
    totals = totals[totals > 0]
    n_nonzero = int(np.count_nonzero(data))
    sig: dict = {"n_cells_with_counts": int(totals.size), "n_nonzero": n_nonzero}
    if totals.size >= 2:
        sig["cv_cell_totals"] = round(_cv(totals), 6)
    if n_nonzero:
        sig["share_ones"] = round(int(np.count_nonzero(data == 1)) / n_nonzero, 6)
        sig["max_value"] = float(data.max())
    sig["checked_totals"] = bool(totals.size >= PLAUSIBILITY_MIN_CELLS)
    sig["checked_values"] = bool(n_nonzero >= PLAUSIBILITY_MIN_NONZERO)
    flags: list[str] = []
    if sig["checked_totals"]:
        if sig["cv_cell_totals"] < TOTALS_CV_NEAR_IDENTICAL:
            flags.append("near_identical_totals")
        elif sig["cv_cell_totals"] < TOTALS_CV_TIGHT:
            flags.append("tight_totals")
    if sig["checked_values"]:
        if sig["max_value"] <= 1:
            flags.append("binary")
        elif sig["share_ones"] <= NO_ONES_MAX_SHARE:
            flags.append("no_ones")
    sig["flags"] = flags
    sig["warn"] = bool({"near_identical_totals", "binary"} & set(flags)
                       or {"no_ones", "tight_totals"} <= set(flags))
    sig["thresholds"] = {
        "status": "provisional, not calibrated on real data (issue #42)",
        "min_cells": PLAUSIBILITY_MIN_CELLS, "min_nonzero": PLAUSIBILITY_MIN_NONZERO,
        "totals_cv_near_identical": TOTALS_CV_NEAR_IDENTICAL,
        "totals_cv_tight": TOTALS_CV_TIGHT, "no_ones_max_share": NO_ONES_MAX_SHARE,
    }
    return sig


def plausibility_warning(sig: dict) -> str | None:
    """The loud warning text for :func:`count_plausibility` signals, or None."""
    if not sig.get("warn"):
        return None
    flags = set(sig["flags"])
    cv, n_cells, nnz = sig.get("cv_cell_totals"), sig["n_cells_with_counts"], sig["n_nonzero"]
    findings = []
    if "near_identical_totals" in flags:
        findings.append(
            f"per-cell totals are near-identical (CV {cv:.2g} over {n_cells} cells)"
            + (" and no non-zero value is 1" if "no_ones" in flags else "")
            + "; real sequencing depth always varies between cells, so this looks like rounded "
              "normalize_total/CPM output or counts downsampled to one fixed depth per cell")
    elif {"no_ones", "tight_totals"} <= flags:
        findings.append(
            f"no non-zero value is 1 and per-cell totals are unusually tight (CV {cv:.2g} over "
            f"{n_cells} cells); this looks like normalized values rounded after some genes "
            "were removed")
    if "binary" in flags:
        findings.append(f"every non-zero value is 1 ({nnz} values): a binarized "
                        "detected/not-detected matrix, not counts")
    return ("COUNTS LOOK NON-RAW: " + "; ".join(findings) + ". The matrix passed the integer "
            "check, so QC, doublet calling and every model will treat it as raw counts. If it "
            "is not raw, re-export the raw counts (layers['counts'], .raw.X or the 10x matrix). "
            "Warning only: the thresholds are provisional until calibrated on real data "
            "(issue #42).")


def _check_count_plausibility(X, rep: IngestReport) -> None:
    """Record the plausibility signals in ``rep``; warn loudly, never refuse (issue #42)."""
    rep.count_plausibility = count_plausibility(X)
    msg = plausibility_warning(rep.count_plausibility)
    if msg:
        rep.warnings.append(msg)
        warnings.warn(msg, RuntimeWarning, stacklevel=3)


def _clean_gene_names(adata: ad.AnnData, rep: IngestReport) -> ad.AnnData:
    """Strip Ensembl versions and aggregate duplicate genes by sum.

    Returns the (possibly rebuilt) AnnData: duplicates change the matrix
    shape, so the object is reconstructed rather than partially mutated. Gene
    matching elsewhere in grnbench is case-insensitive, therefore identifiers
    that differ only by case are duplicates under the same contract and are
    aggregated here too. The first loaded spelling is retained and the
    case-fold merge is reported separately. Every layer uses the same mapping.
    """
    names = pd.Index(adata.var_names.astype(str))
    if names.str.match(ENSEMBL_ANY_RE).mean() > 0.9:
        for col in SYMBOL_COLUMNS:
            if col not in adata.var.columns:
                continue
            sym = adata.var[col].astype(str)
            usable = (sym.str.len() > 0) & ~sym.str.lower().isin({"nan", "none", ""}) \
                & ~sym.str.match(ENSEMBL_ANY_RE)
            if usable.mean() > 0.5:
                adata.var["gene_id_as_loaded"] = names.to_numpy()
                names = pd.Index(np.where(usable.to_numpy(), sym.to_numpy(), names.to_numpy()))
                rep.notes.append(f"Ensembl IDs mapped to symbols from var['{col}'] "
                                 f"({int(usable.sum())}/{len(names)} genes; IDs kept in "
                                 "var['gene_id_as_loaded'])")
                break
    stripped = names.str.match(ENSEMBL_VERSION_RE)
    if stripped.any():
        names = pd.Index([
            ENSEMBL_VERSION_RE.match(n).group(1) if m else n
            for n, m in zip(names, stripped)
        ])
        rep.version_suffixes_stripped = int(stripped.sum())

    exact_duplicates = int(names.duplicated(keep="first").sum())
    folded = pd.Index([name.casefold() for name in names])
    folded_dup_mask = folded.duplicated(keep="first")
    if not folded_dup_mask.any():
        adata.var.index = names
        if rep.version_suffixes_stripped:
            rep.notes.append(
                f"{rep.version_suffixes_stripped} Ensembl version suffixes stripped."
            )
        adata.uns["grnbench_ingest"] = rep.to_dict()
        return adata

    unique_folded = pd.Index(pd.unique(folded))
    codes = pd.Categorical(folded, categories=unique_folded).codes
    ind = sp.csr_matrix(
        (
            np.ones(len(names), dtype=np.float32),
            (codes, np.arange(len(names))),
        ),
        shape=(len(unique_folded), len(names)),
    )

    def _agg(m):
        return ((ind @ sp.csr_matrix(m).T).T).tocsr()

    first = ~folded_dup_mask
    new_var = adata.var.loc[first].copy()
    representative_names = names[first]
    new_var.index = pd.Index(representative_names, name=new_var.index.name)

    rebuilt = ad.AnnData(
        X=_agg(adata.X),
        obs=adata.obs.copy(),
        var=new_var,
        # anndata>=0.13 aliases X as layers[None]; skip it when rebuilding
        layers={k: _agg(v) for k, v in adata.layers.items() if k is not None},
    )
    rebuilt.uns = dict(adata.uns)
    rep.n_duplicates_merged = int(folded_dup_mask.sum())
    rep.n_case_duplicates_merged = rep.n_duplicates_merged - exact_duplicates
    if exact_duplicates:
        rep.notes.append(f"{exact_duplicates} duplicate gene names aggregated by sum.")
    if rep.n_case_duplicates_merged:
        rep.notes.append(
            f"{rep.n_case_duplicates_merged} case-insensitive duplicate gene names "
            "aggregated by sum; first loaded spelling retained."
        )
    if rep.version_suffixes_stripped:
        rep.notes.append(
            f"{rep.version_suffixes_stripped} Ensembl version suffixes stripped."
        )
    rebuilt.uns["grnbench_ingest"] = rep.to_dict()
    return rebuilt
