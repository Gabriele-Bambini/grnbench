"""Calibrate D-027 identifier guards on the approved real-data root.

    python benchmarks/identifier_policy_calibration.py DATA_ROOT
    python benchmarks/identifier_policy_calibration.py DATA_ROOT --tsv identifier-policy.tsv

The script walks only formats already approved by ``hpc/fetch_engine_datasets.sh``:
``*.h5ad``, named 10x feature-barcode ``*.h5`` files and 10x MTX folders.  It runs the production ingest
path, reports collision counts and distinguishes a D-027 policy refusal from an
unrelated ingest refusal.  Paths are relative to DATA_ROOT and messages are
scrubbed, so the output can be attached to an issue or PR.

Exit 0 means at least one input was inspected and D-027 refused none.  Exit 2
means one or more policy refusals need review; exit 1 means no input was found.
No expression values, barcodes or local absolute paths are written.
"""
from __future__ import annotations

import argparse
import gc
import sys
import warnings
from pathlib import Path

COLUMNS = (
    "path", "status", "input_genes", "output_genes", "mapping_source",
    "organism_scope", "organism_source", "case_fold_action",
    "exact_duplicate_groups", "case_fold_groups", "cross_convention_groups", "message",
)

TENX_H5_NAMES = {"raw_feature_bc_matrix.h5", "filtered_feature_bc_matrix.h5"}


def find_inputs(root: Path) -> list[Path]:
    """Approved ingest-readable inputs, in stable relative-path order."""
    found = {path for path in root.rglob("*.h5ad") if path.is_file()}
    found |= {path for path in root.rglob("*.h5")
              if path.is_file() and path.name in TENX_H5_NAMES}
    found |= {path.parent for name in ("matrix.mtx", "matrix.mtx.gz")
              for path in root.rglob(name)}
    return sorted(found, key=lambda path: path.relative_to(root).as_posix())


def _scrubber(root: Path):
    hidden = sorted({str(root.resolve()), str(root.absolute()), str(Path.home())},
                    key=len, reverse=True)

    def scrub(message: object) -> str:
        text = str(message)
        for value in hidden:
            text = text.replace(value, "<root>" if value != str(Path.home()) else "~")
        return " ".join(text.split())[:300]

    return scrub


def inspect_file(path: Path, root: Path) -> dict:
    """Inspect one file through production ingest and return a portable row."""
    from grnbench.ingest import IdentifierPolicyError, IngestError, load_any

    row = {"path": path.relative_to(root).as_posix()}
    scrub = _scrubber(root)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            adata, report = load_any(path)
    except IdentifierPolicyError as exc:
        return {**row, "status": "policy_refused", "message": scrub(exc)}
    except IngestError as exc:
        return {**row, "status": "ingest_refused", "message": scrub(exc)}
    except Exception as exc:  # noqa: BLE001 - one broken file must not hide the others
        return {**row, "status": "error", "message": scrub(f"{type(exc).__name__}: {exc}")}

    provenance = report.identifier_provenance
    collisions = provenance.get("collision_summary") or {}
    organism_policy = provenance.get("organism_policy") or {}
    result = {
        **row,
        "status": "accepted",
        "input_genes": provenance.get("input_gene_count"),
        "output_genes": provenance.get("output_gene_count"),
        "mapping_source": (provenance.get("mapping_source") or {}).get("kind"),
        "organism_scope": organism_policy.get("scope"),
        "organism_source": organism_policy.get("source"),
        "case_fold_action": organism_policy.get("case_fold_action"),
        "exact_duplicate_groups": len(provenance.get("exact_duplicate_groups") or []),
        "case_fold_groups": collisions.get("case_fold_groups"),
        "cross_convention_groups": collisions.get("cross_convention_groups"),
        "message": "",
    }
    del adata
    gc.collect()
    return result


def _format(value) -> str:
    return "" if value is None else str(value)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("data_root")
    parser.add_argument("--tsv", metavar="OUT", help="also write the portable table to this TSV")
    args = parser.parse_args(argv)

    root = Path(args.data_root)
    if not root.is_dir():
        print(f"not a folder: {args.data_root}", file=sys.stderr)
        return 1
    inputs = find_inputs(root)
    if not inputs:
        print("no approved ingest-readable inputs found", file=sys.stderr)
        return 1

    lines = ["\t".join(COLUMNS)]
    print(lines[0], flush=True)
    rows = []
    for path in inputs:
        row = inspect_file(path, root)
        rows.append(row)
        line = "\t".join(_format(row.get(column)) for column in COLUMNS)
        lines.append(line)
        print(line, flush=True)
    if args.tsv:
        Path(args.tsv).write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")

    counts = {status: sum(row["status"] == status for row in rows) for status in (
        "accepted", "policy_refused", "ingest_refused", "error")}
    print(
        f"files: {len(rows)}, accepted: {counts['accepted']}, "
        f"policy refused: {counts['policy_refused']}, "
        f"other ingest refused: {counts['ingest_refused']}, errors: {counts['error']}")
    return 2 if counts["policy_refused"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
