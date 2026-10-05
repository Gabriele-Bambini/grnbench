"""Count plausibility calibration: would ingest's non-raw warning fire on real raw data?

    python benchmarks/count_plausibility.py DATA_ROOT
    python benchmarks/count_plausibility.py DATA_ROOT --tsv plausibility.tsv

Walks DATA_ROOT and loads every ``*.h5ad``, every 10x ``*.h5`` and every 10x MTX folder
(a folder holding ``matrix.mtx[.gz]``) with grnbench's own ingest (``load_any``), so the
counts are taken exactly where a run would take them (X, ``layers['counts']`` or ``.raw``).
Per file it prints the signals of ``grnbench.ingest.count_plausibility``: CV of per-cell
totals, share of 1s among non-zero values, maximum, binary or not, the flags and whether
the warning would fire. Inputs ingest refuses (not integer counts) are listed as ``refused``
with the start of the message; they never count as a warning. Text tables are not walked:
their layout needs an orientation decision per file.

Issue #42: the thresholds in ``grnbench.ingest`` are provisional. Run this on the HPC over
the approved data root; the warning may turn into a refusal only after a record shows
``would warn: 0`` on real raw data. Paths are printed relative to DATA_ROOT and messages
have DATA_ROOT replaced by ``<root>``, so the output carries no personal path. Exit code 1
when nothing loadable was found (a wrong root must not read as "zero false alarms").
"""
from __future__ import annotations

import argparse
import gc
import sys
import warnings
from pathlib import Path

COLUMNS = ("path", "status", "counts_source", "n_cells", "n_genes", "cv_cell_totals",
           "share_ones", "max_value", "binary", "flags", "message")


def find_inputs(root: Path) -> list[Path]:
    """Every ingest-readable input under ``root``, sorted: h5ad, 10x h5, 10x MTX folders."""
    found = {p for pattern in ("*.h5ad", "*.h5") for p in root.rglob(pattern) if p.is_file()}
    found |= {p.parent for name in ("matrix.mtx", "matrix.mtx.gz") for p in root.rglob(name)}
    return sorted(found, key=lambda p: p.relative_to(root).as_posix())


def _fmt(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def inspect_file(path: Path, root: Path) -> dict:
    """One output row for ``path``: the plausibility signals, or why it could not be read."""
    from grnbench.ingest import IngestError, load_any

    row = {"path": path.relative_to(root).as_posix()}
    hide = sorted({str(root.resolve()), str(root.absolute())}, key=len, reverse=True)

    def scrub(text: str) -> str:
        for h in hide:
            text = text.replace(h, "<root>")
        text = text.replace(str(Path.home()), "~")
        return " ".join(text.split())[:300]

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # the warning text is read from the report below
            adata, rep = load_any(path)
    except IngestError as exc:
        return {**row, "status": "refused", "message": scrub(str(exc))}
    except Exception as exc:  # noqa: BLE001 - one unreadable file must not stop the walk
        return {**row, "status": "error", "message": scrub(f"{type(exc).__name__}: {exc}")}
    sig = rep.count_plausibility
    row.update(
        status="warn" if sig.get("warn") else "ok",
        counts_source=rep.counts_source, n_cells=int(adata.n_obs), n_genes=int(adata.n_vars),
        cv_cell_totals=sig.get("cv_cell_totals"), share_ones=sig.get("share_ones"),
        max_value=sig.get("max_value"),
        binary="binary" in sig.get("flags", []),
        flags=",".join(sig.get("flags", [])) or "-",
        message=scrub(" ".join(rep.warnings)),
    )
    del adata
    gc.collect()
    return row


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("data_root")
    parser.add_argument("--tsv", metavar="OUT", help="also write the table to this TSV file")
    args = parser.parse_args(argv)

    from grnbench import ingest

    root = Path(args.data_root)
    if not root.is_dir():
        print(f"not a folder: {args.data_root}", file=sys.stderr)
        return 1
    print(f"# thresholds (provisional, issue #42): totals CV near-identical < "
          f"{ingest.TOTALS_CV_NEAR_IDENTICAL}, tight < {ingest.TOTALS_CV_TIGHT}, no 1s <= "
          f"{ingest.NO_ONES_MAX_SHARE}; min cells {ingest.PLAUSIBILITY_MIN_CELLS}, "
          f"min non-zero values {ingest.PLAUSIBILITY_MIN_NONZERO}")
    lines = ["\t".join(COLUMNS)]
    print(lines[0], flush=True)
    rows = []
    for path in find_inputs(root):
        row = inspect_file(path, root)
        rows.append(row)
        lines.append("\t".join(_fmt(row.get(c)) for c in COLUMNS))
        print(lines[-1], flush=True)
    if args.tsv:
        Path(args.tsv).write_text("\n".join(lines) + "\n", encoding="utf-8")
    count = {s: sum(r["status"] == s for r in rows) for s in ("ok", "warn", "refused", "error")}
    print(f"files: {len(rows)}, would warn: {count['warn']}, ok: {count['ok']}, "
          f"refused: {count['refused']}, errors: {count['error']}")
    return 0 if rows else 1


if __name__ == "__main__":
    sys.exit(main())
