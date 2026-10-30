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

import hashlib
import re
import warnings
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

ENSEMBL_VERSION_RE = re.compile(r"^((?:ENSG|ENSMUSG)\d+)\.\d+$")
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
CASE_COLLISION_MAX_FRACTION = 0.01
CROSS_CONVENTION_MAX_FRACTION = 0.001
CROSS_CONVENTION_MIN_GROUPS = 2
_ORF_TOKEN_RE = re.compile(r"(?<=[0-9XY])orf(?=[0-9])")
_MAPPING_ID_HEADERS = {"id", "gene_id", "ensembl_id", "feature_id"}
_MAPPING_SYMBOL_HEADERS = {"symbol", "gene_symbol", "gene_name", "feature_name"}
_ORGANISM_ID_COLUMNS = (
    "gene_ids", "gene_id", "feature_id", "ensembl_id", "gene_id_as_loaded",
)

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


class IdentifierPolicyError(IngestError):
    """Raised when D-028 refuses an unsafe identifier set."""


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
    # D-027: portable, path-free source-to-representative audit.
    identifier_provenance: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return self.__dict__.copy()


def load_any(path: str | Path,
             orientation: str = "auto",
             gene_mapping_path: str | Path | None = None) -> tuple[ad.AnnData, IngestReport]:
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
    try:
        adata, rep = _read_any(p, orientation)
    except IngestError:
        raise
    except (OSError, ValueError, KeyError, TypeError) as exc:    # the reader's own text is jargon: say what to do
        raise _unreadable(p, exc) from exc
    gene_mapping = mapping_source = None
    if gene_mapping_path is not None:
        gene_mapping, mapping_source = load_gene_mapping(gene_mapping_path)
    _ensure_counts_backed(adata, rep)
    _assert_finite(adata.X)
    _check_count_plausibility(adata.X, rep)
    adata = _clean_gene_names(
        adata, rep, gene_mapping=gene_mapping, mapping_source=mapping_source)
    return adata, rep


def _read_any(p: Path, orientation: str) -> tuple[ad.AnnData, IngestReport]:
    if p.is_dir():
        return _load_10x_mtx_dir(p)
    sfx = p.suffix.lower()
    if sfx == ".h5ad":
        return _load_h5ad(p)
    if sfx == ".h5":
        # NO var_names_make_unique() here: CellRanger references do contain
        # duplicated symbols (e.g. TBCE), and renaming them to SYMBOL-1
        # would make _clean_gene_names' aggregation branch dead - the same
        # data would then yield a different matrix than the MTX path.
        return sc.read_10x_h5(p), IngestReport(path=str(p), fmt="10x_h5")
    if sfx == ".loom":
        try:
            adata = sc.read_loom(p)       # anndata imports loompy lazily, inside read_loom
        except ImportError as e:
            raise IngestError("Reading .loom files needs the optional 'loompy' package: "
                              "pip install loompy (from a clone of the repository: "
                              "pip install -e '.[loom]')") from e
        return adata, IngestReport(path=str(p), fmt="loom")
    if sfx in {".csv", ".tsv", ".txt"}:
        return _load_table(p, orientation)
    raise IngestError(
        f"Unsupported format '{sfx}'. Supported: .h5ad, 10x .h5, "
        "10x MTX directory, .loom, .csv/.tsv/.txt"
    )


def _unreadable(p: Path, exc: Exception) -> IngestError:
    """An IngestError that names the file and what to do, for a reader that failed with library jargon."""
    if p.is_dir():
        msg = (f"{p.name} is not a 10x MTX folder: it needs matrix.mtx(.gz), barcodes.tsv(.gz) and "
               "features.tsv(.gz) (or genes.tsv)")
    elif p.stat().st_size == 0:
        msg = f"{p.name} is empty (0 bytes): re-download or re-export it"
    elif isinstance(exc, UnicodeDecodeError):
        msg = f"{p.name} looks like a binary file, not a text table"
    else:
        msg = f"{p.name} is not a valid {p.suffix.lower() or 'data'} file (damaged, cut short or another format)"
    return IngestError(msg)


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
    # pandas renames a repeated header field (GENE, GENE -> GENE, GENE.1); put the file's own header back so
    # _clean_gene_names aggregates the duplicates like every other format (blank fields keep pandas' names).
    header = pd.read_csv(p, sep=sep, header=None, nrows=1, dtype=str, keep_default_na=False)
    if header.shape[1] == df.shape[1] + 1:
        df.columns = [raw or kept for raw, kept in zip(header.iloc[0, 1:], df.columns)]
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


def load_gene_mapping(path: str | Path) -> tuple[dict[str, str], dict]:
    """Read a deterministic two-column ``id -> symbol`` table.

    The table may be comma- or tab-separated and may have an ``id/symbol``
    header.  Conflicting targets for one identifier are refused.  Only the
    content digest and row counts are returned as provenance; the local path is
    deliberately never copied into a manifest or portable diagnostic.
    """
    source = Path(path)
    if not source.is_file():
        raise IngestError(f"gene mapping table does not exist: {source}")
    raw = source.read_bytes()
    try:
        frame = pd.read_csv(source, sep=None, engine="python", header=None, dtype=str,
                            keep_default_na=False, comment="#")
    except Exception as exc:  # noqa: BLE001 - turn parser internals into user advice
        raise IngestError(f"could not read gene mapping table: {exc}") from exc
    if frame.shape[1] != 2:
        raise IngestError(
            f"gene mapping table must have exactly two columns (id, symbol); found {frame.shape[1]}")
    frame = frame.apply(lambda column: column.astype(str).str.strip())
    if len(frame):
        first = (frame.iat[0, 0].casefold(), frame.iat[0, 1].casefold())
        if first[0] in _MAPPING_ID_HEADERS and first[1] in _MAPPING_SYMBOL_HEADERS:
            frame = frame.iloc[1:].reset_index(drop=True)
    if frame.empty or (frame == "").any(axis=None):
        raise IngestError("gene mapping table contains no mappings or has an empty id/symbol")

    mapping: dict[str, str] = {}
    for source_id, symbol in frame.itertuples(index=False, name=None):
        previous = mapping.get(source_id)
        if previous is not None and previous != symbol:
            raise IngestError(
                "gene mapping table maps one identifier to multiple symbols: "
                f"{source_id!r} -> {previous!r}/{symbol!r}")
        mapping[source_id] = symbol
    metadata = {
        "kind": "user_table",
        "sha256": hashlib.sha256(raw).hexdigest(),
        "rows": int(len(frame)),
        "unique_ids": int(len(mapping)),
    }
    return mapping, metadata


def _strip_supported_ensembl_version(name: str) -> str:
    match = ENSEMBL_VERSION_RE.match(name)
    return match.group(1) if match else name


def _symbol_is_usable(value: str) -> bool:
    return bool(value) and value.casefold() not in {"nan", "none"} \
        and not bool(ENSEMBL_ANY_RE.match(value))


def _name_convention(name: str) -> str | None:
    """Measured D-027 spelling class after removing the ``C1orf112`` token."""
    without_orf = _ORF_TOKEN_RE.sub("", name)
    letters = "".join(char for char in without_orf if char.isalpha())
    if not letters:
        return None
    if letters.isupper():
        return "upper"
    if letters.islower():
        return "lower"
    return "mixed"


def _organism_label(value: object) -> str:
    """Portable organism label for the identifier policy, including ontology IDs."""
    label = " ".join(str(value).strip().lower().split())
    if label in {"human", "homo sapiens", "ncbitaxon:9606", "9606"}:
        return "human"
    if label in {"mouse", "mus musculus", "ncbitaxon:10090", "10090"}:
        return "mouse"
    return label


def _identifier_organism_policy(adata: ad.AnnData, names: pd.Index) -> dict:
    """Determine whether human/mouse case-fold merging is safe under D-028.

    Curated organism metadata wins.  Human/mouse Ensembl identifiers are the
    next strongest evidence, followed by the existing symbol-casing hint.  A
    named third organism or ambiguous input is never coerced into the
    human/mouse namespace: distinct case-only spellings are preserved.
    """
    if "organism" in adata.obs.columns and len(adata.obs):
        values = {
            _organism_label(value)
            for value in pd.unique(adata.obs["organism"].dropna().astype(str))
            if _organism_label(value) not in {"", "nan", "none"}
        }
        if values:
            if values == {"human"}:
                scope, label = "human", "human"
            elif values == {"mouse"}:
                scope, label = "mouse", "mouse"
            elif values <= {"human", "mouse"}:
                scope, label = "human_mouse", "human+mouse"
            else:
                scope, label = "other", ", ".join(sorted(values))[:120]
            return {
                "scope": scope,
                "label": label,
                "source": "metadata:obs['organism']",
                "case_fold_action": "merge_and_guard" if scope != "other" else "preserve_and_warn",
            }

    id_collections = [names]
    for column in _ORGANISM_ID_COLUMNS:
        if column in adata.var.columns:
            id_collections.append(pd.Index(adata.var[column].astype(str)))
    for identifiers in id_collections:
        human = float(identifiers.str.match(r"^ENSG\d+(?:\.\d+)?$").mean())
        mouse = float(identifiers.str.match(r"^ENSMUSG\d+(?:\.\d+)?$").mean())
        other_ensembl = float(
            (identifiers.str.match(ENSEMBL_ANY_RE)
             & ~identifiers.str.match(r"^(?:ENSG|ENSMUSG)\d+(?:\.\d+)?$")).mean())
        if human > 0.5:
            return {"scope": "human", "label": "human", "source": "gene_identifiers",
                    "case_fold_action": "merge_and_guard"}
        if mouse > 0.5:
            return {"scope": "mouse", "label": "mouse", "source": "gene_identifiers",
                    "case_fold_action": "merge_and_guard"}
        if other_ensembl > 0.5:
            return {"scope": "other", "label": "non-human/mouse Ensembl",
                    "source": "gene_identifiers", "case_fold_action": "preserve_and_warn"}

    from grnbench.hygiene import organism_hint
    hint = organism_hint(names)
    if hint in {"human-like", "mouse-like"}:
        scope = hint.removesuffix("-like")
        return {"scope": scope, "label": scope, "source": "symbol_casing",
                "case_fold_action": "merge_and_guard"}
    return {"scope": "unknown", "label": None, "source": "insufficient_evidence",
            "case_fold_action": "preserve_and_warn"}


def _portable_ingest_summary(rep: IngestReport) -> dict:
    """Small AnnData ``uns`` copy; complete mappings live in run artifacts."""
    value = rep.to_dict()
    provenance = value.pop("identifier_provenance", {}) or {}
    if provenance:
        value["identifier_summary"] = {
            "policy": provenance.get("policy"),
            "input_gene_count": provenance.get("input_gene_count"),
            "output_gene_count": provenance.get("output_gene_count"),
            "organism_policy": provenance.get("organism_policy"),
            "mapping_source": provenance.get("mapping_source"),
            "mapping_outcome": provenance.get("mapping_outcome"),
            "collision_summary": provenance.get("collision_summary"),
        }
    return value


def _clean_gene_names(
    adata: ad.AnnData,
    rep: IngestReport,
    *,
    gene_mapping: dict[str, str] | None = None,
    mapping_source: dict | None = None,
) -> ad.AnnData:
    """Map, audit and safely aggregate gene identifiers under D-028.

    Input-carried symbols or an explicit, fingerprinted two-column mapping are
    the only mapping sources.  Exact duplicates are always summed.  Distinct
    spellings are summed only for confidently identified human/mouse inputs and
    while the measured case-fold and cross-convention guardrails remain below
    their refusal thresholds.  For another or unknown organism, distinct
    case-only spellings are preserved and reported.  Exact duplicates are
    always summed, and every layer uses the same mapping.
    """
    loaded_names = pd.Index(adata.var_names.astype(str))
    names = loaded_names.copy()
    input_ensembl_fraction = float(loaded_names.str.match(ENSEMBL_ANY_RE).mean())
    applied_mapping = np.zeros(len(names), dtype=bool)
    source_info: dict = {"kind": "none"}
    mapping_transformation = None

    if gene_mapping is not None:
        mapped_names = []
        for name in loaded_names:
            target = gene_mapping.get(name)
            if target is None:
                target = gene_mapping.get(_strip_supported_ensembl_version(name))
            if target is not None and _symbol_is_usable(target):
                mapped_names.append(target)
                applied_mapping[len(mapped_names) - 1] = target != name
            else:
                mapped_names.append(name)
        names = pd.Index(mapped_names)
        source_info = dict(mapping_source or {"kind": "user_table"})
        mapping_transformation = "user_mapping"
    elif input_ensembl_fraction > 0.9:
        # Several producers carry more than one symbol-like column.  A sparse
        # early column must not hide a complete later one merely because it has
        # one usable value: select the column with the greatest usable share.
        candidates = []
        for column in SYMBOL_COLUMNS:
            if column in adata.var.columns:
                symbols = adata.var[column].astype(str).str.strip()
                usable = symbols.map(_symbol_is_usable).to_numpy(dtype=bool)
                candidates.append((int(usable.sum()), column, symbols, usable))
        if candidates:
            _, column, symbols, usable = max(candidates, key=lambda item: item[0])
            if usable.any():
                names = pd.Index(np.where(usable, symbols.to_numpy(), loaded_names.to_numpy()))
                applied_mapping = usable & (names.to_numpy() != loaded_names.to_numpy())
                source_info = {"kind": "input_column", "column": str(column)}
                mapping_transformation = "input_symbol_mapping"
    padded = names.str.strip() != names
    if padded.any():            # ' TP53' and 'TP53\xa0' match no TF list and no gene set: say so and fix them
        rep.notes.append(f"{int(padded.sum())} gene name(s) had leading or trailing whitespace and were "
                         f"stripped (e.g. {names[padded][0]!r})")
        names = pd.Index(names.str.strip())

    if applied_mapping.any():
        adata.var["gene_id_as_loaded"] = loaded_names.to_numpy()
        source_label = ("a supplied table" if source_info["kind"] == "user_table"
                        else f"var['{source_info['column']}']")
        rep.notes.append(
            f"gene identifiers mapped to symbols from {source_label} "
            f"({int(applied_mapping.sum())}/{len(names)} features; IDs kept in "
            "var['gene_id_as_loaded'])")

    remaining_ensembl_fraction = float(names.str.match(ENSEMBL_ANY_RE).mean())
    if input_ensembl_fraction > 0.9 and remaining_ensembl_fraction > 0.9:
        raise IdentifierPolicyError(
            "more than 90% of input gene names are Ensembl IDs and no usable symbol source "
            "mapped enough of them. Supply a two-column id-to-symbol table with "
            "--gene-mapping, or export the input with a feature_name/gene_symbol column; "
            "keeping Ensembl IDs would produce no overlap with symbol-based references.")

    before_strip = names.copy()
    names = pd.Index([_strip_supported_ensembl_version(name) for name in names])
    stripped_mask = names.to_numpy() != before_strip.to_numpy()
    rep.version_suffixes_stripped = int(stripped_mask.sum())

    organism_policy = _identifier_organism_policy(adata, names)
    merge_casefold = organism_policy["case_fold_action"] == "merge_and_guard"
    spelling_counts = Counter(map(str, names))
    folded_spellings: dict[str, list[str]] = {}
    for name in map(str, names):
        spellings = folded_spellings.setdefault(name.casefold(), [])
        if name not in spellings:
            spellings.append(name)
    case_groups_raw = [spellings for spellings in folded_spellings.values() if len(spellings) > 1]
    case_groups = []
    n_cross = 0
    for spellings in case_groups_raw:
        conventions = {_name_convention(name) for name in spellings} - {None}
        cross = len(conventions) >= 2
        n_cross += int(cross)
        case_groups.append({
            "representative": spellings[0] if merge_casefold else None,
            "spellings": spellings,
            "method": "sum" if merge_casefold else "preserve",
            "cross_convention": cross,
        })

    n_features = len(names)
    cross_convention_limit = max(
        CROSS_CONVENTION_MIN_GROUPS,
        CROSS_CONVENTION_MAX_FRACTION * n_features,
    )
    if merge_casefold and n_cross > cross_convention_limit:
        raise IdentifierPolicyError(
            f"identifier hygiene found {n_cross} cross-convention case-fold groups among "
            f"{n_features} genes, above the max(2 groups, 0.1%) safety limit. This can indicate mixed "
            "human/mouse annotation; provide one organism's corrected identifiers.")
    if merge_casefold and len(case_groups) > CASE_COLLISION_MAX_FRACTION * n_features:
        raise IdentifierPolicyError(
            f"identifier hygiene found {len(case_groups)} distinct-spelling case-fold groups "
            f"among {n_features} genes, above the 1% safety limit. This looks like an "
            "annotation collision set; provide a corrected id-to-symbol table instead of "
            "silently merging it.")

    aggregation_keys = pd.Index(
        [name.casefold() for name in names] if merge_casefold else list(map(str, names)))
    duplicate_mask = aggregation_keys.duplicated(keep="first")
    unique_keys = pd.Index(pd.unique(aggregation_keys))
    representative_by_key: dict[str, str] = {}
    for name, key in zip(names, aggregation_keys):
        representative_by_key.setdefault(str(key), str(name))

    mapping_rows: dict[tuple[str, str, tuple[str, ...]], int] = {}
    for position, (loaded, mapped, stripped, fold) in enumerate(
            zip(loaded_names, before_strip, names, aggregation_keys)):
        transformations = []
        if applied_mapping[position] and mapping_transformation:
            transformations.append(mapping_transformation)
        if mapped != stripped:
            transformations.append("strip_ensembl_version")
        representative = representative_by_key[str(fold)]
        if merge_casefold and str(stripped) != representative:
            transformations.append("case_fold_merge")
        key = (str(loaded), representative, tuple(transformations))
        mapping_rows[key] = mapping_rows.get(key, 0) + 1

    exact_groups = [
        {"representative": name, "occurrences": int(count), "method": "sum"}
        for name, count in spelling_counts.items() if count > 1
    ]
    exact_duplicates = sum(group["occurrences"] - 1 for group in exact_groups)
    rep.n_duplicates_merged = int(duplicate_mask.sum())
    rep.n_case_duplicates_merged = (
        rep.n_duplicates_merged - int(exact_duplicates) if merge_casefold else 0)
    rep.identifier_provenance = {
        "policy": "D-028",
        "input_gene_count": int(len(loaded_names)),
        "output_gene_count": int(len(unique_keys)),
        "aggregation_method": "sum",
        "organism_policy": organism_policy,
        "mapping_source": source_info,
        "mapping_outcome": {
            "mapped_features": int(applied_mapping.sum()),
            "unchanged_features": int(len(names) - applied_mapping.sum()),
            "input_ensembl_fraction": round(input_ensembl_fraction, 6),
            "remaining_ensembl_fraction": round(remaining_ensembl_fraction, 6),
        },
        "version_suffixes_stripped": rep.version_suffixes_stripped,
        "source_to_representative": [
            {
                "as_loaded": loaded,
                "representative": representative,
                "occurrences": occurrences,
                "transformations": list(transformations),
            }
            for (loaded, representative, transformations), occurrences in mapping_rows.items()
        ],
        "exact_duplicate_groups": exact_groups,
        "case_fold_groups": case_groups,
        "collision_summary": {
            "case_fold_groups": int(len(case_groups)),
            "case_fold_fraction": round(len(case_groups) / n_features, 6),
            "case_fold_limit": CASE_COLLISION_MAX_FRACTION,
            "cross_convention_groups": int(n_cross),
            "cross_convention_fraction": round(n_cross / n_features, 6),
            "cross_convention_limit": CROSS_CONVENTION_MAX_FRACTION,
            "cross_convention_min_groups": CROSS_CONVENTION_MIN_GROUPS,
            "case_fold_action": organism_policy["case_fold_action"],
        },
    }

    if exact_duplicates:
        rep.notes.append(f"{exact_duplicates} exact duplicate gene names aggregated by sum.")
    if rep.n_case_duplicates_merged:
        rep.notes.append(
            f"{rep.n_case_duplicates_merged} case-insensitive duplicate gene names "
            "aggregated by sum; first loaded spelling retained.")
    elif case_groups and not merge_casefold:
        message = (
            f"{len(case_groups)} distinct-spelling case-fold group(s) preserved because "
            f"the organism scope is {organism_policy['scope']!r}, not confidently human/mouse.")
        rep.notes.append(message)
        rep.warnings.append(message)
    if rep.version_suffixes_stripped:
        rep.notes.append(
            f"{rep.version_suffixes_stripped} human/mouse Ensembl version suffixes stripped.")

    if not duplicate_mask.any():
        adata.var.index = names
        adata.uns["grnbench_ingest"] = _portable_ingest_summary(rep)
        return adata

    codes = pd.Categorical(aggregation_keys, categories=unique_keys).codes
    indicator = sp.csr_matrix(
        (np.ones(len(names), dtype=np.float32), (codes, np.arange(len(names)))),
        shape=(len(unique_keys), len(names)),
    )

    def _aggregate(matrix):
        return ((indicator @ sp.csr_matrix(matrix).T).T).tocsr()

    first = ~duplicate_mask
    new_var = adata.var.loc[first].copy()
    new_var.index = pd.Index(names[first], name=new_var.index.name)
    rebuilt = ad.AnnData(
        X=_aggregate(adata.X),
        obs=adata.obs.copy(),
        var=new_var,
        # anndata>=0.13 aliases X as layers[None]; skip it when rebuilding
        layers={key: _aggregate(value) for key, value in adata.layers.items() if key is not None},
    )
    rebuilt.uns = dict(adata.uns)
    rebuilt.uns["grnbench_ingest"] = _portable_ingest_summary(rep)
    return rebuilt
