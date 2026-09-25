"""Dataset ingestion: common scRNA-seq formats -> counts-backed AnnData.

Contract: after :func:`load_any`, ``adata.X`` holds **raw counts** (original
counts are preserved in ``adata.layers["counts"]``). Anything that does not
look like counts is rejected loudly — silent normalization guessing is how
benchmarks rot. Gene-name hygiene (version suffixes, duplicates aggregated by
sum) is applied and logged in ``adata.uns["grnbench_ingest"]``.
"""
from __future__ import annotations

import re
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
        _validate_counts_layer(adata.X, layer)
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


def _validate_counts_layer(x_as_loaded, layer) -> None:
    """Fail loudly when the trusted counts layer is itself broken."""
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
        raise IngestError(f"layers['counts'] does not look like raw counts: {why}")
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
    msg = (
        "Input does not look like raw counts: " + why + ". GRN preprocessing "
        "must start from counts; re-export the dataset with raw counts (e.g. "
        "layers['counts'] in the h5ad) or pass assume_counts=True if you are "
        "certain (not recommended)."
    )
    raise IngestError(msg)


def _clean_gene_names(adata: ad.AnnData, rep: IngestReport) -> ad.AnnData:
    """Strip Ensembl version suffixes and aggregate duplicate genes by sum.

    Returns the (possibly rebuilt) AnnData: duplicates change the matrix
    shape, so the object is reconstructed rather than partially mutated.
    Every layer is aggregated with the same mapping.
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

    dup_mask = names.duplicated(keep="first")
    if not dup_mask.any():
        adata.var.index = names
        if rep.version_suffixes_stripped:
            rep.notes.append(
                f"{rep.version_suffixes_stripped} Ensembl version suffixes stripped."
            )
        adata.uns["grnbench_ingest"] = rep.to_dict()
        return adata

    uniq = pd.Index(pd.unique(names))
    codes = pd.Categorical(names, categories=uniq).codes
    ind = sp.csr_matrix(
        (
            np.ones(len(names), dtype=np.float32),
            (codes, np.arange(len(names))),
        ),
        shape=(len(uniq), len(names)),
    )

    def _agg(m):
        return ((ind @ sp.csr_matrix(m).T).T).tocsr()

    first = ~dup_mask
    new_var = adata.var.loc[first].copy()
    new_var.index = pd.Index([str(u) for u in uniq], name=new_var.index.name)

    rebuilt = ad.AnnData(
        X=_agg(adata.X),
        obs=adata.obs.copy(),
        var=new_var,
        # anndata>=0.13 aliases X as layers[None]; skip it when rebuilding
        layers={k: _agg(v) for k, v in adata.layers.items() if k is not None},
    )
    rebuilt.uns = dict(adata.uns)
    rep.n_duplicates_merged = int(dup_mask.sum())
    rep.notes.append(f"{rep.n_duplicates_merged} duplicate gene names aggregated by sum.")
    if rep.version_suffixes_stripped:
        rep.notes.append(
            f"{rep.version_suffixes_stripped} Ensembl version suffixes stripped."
        )
    rebuilt.uns["grnbench_ingest"] = rep.to_dict()
    return rebuilt
