"""Compare scores across explicitly declared independent experimental units."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import pandas as pd

from grnbench.eval.comparison import compare_units, render_comparison, scores_from_long


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="grnbench compare", description=__doc__)
    p.add_argument("scores", help="CSV/TSV with unit,method,score columns")
    p.add_argument("--out", required=True)
    p.add_argument("--unit-description", required=True, help="what makes one row an independent experimental unit")
    p.add_argument("--assume-independent-units", action="store_true", required=True,
                   help="explicitly assert the independent-unit assumption; does not verify it")
    p.add_argument("--alpha", type=float, default=.05)
    p.add_argument("--margin", type=float, default=0., help="prespecified meaningful difference in score units")
    args = p.parse_args(argv)
    try:
        if not args.unit_description.strip():
            raise ValueError("unit description cannot be empty")
        table = pd.read_csv(args.scores, sep=None, engine="python",
                            dtype={"unit": str, "method": str})
        result = compare_units(scores_from_long(table), args.alpha, args.margin)
        from grnbench.report import environment_block, file_fingerprint
        result.update(unit_description=args.unit_description,
                      independence_asserted_by_user=True,
                      input=file_fingerprint(args.scores), environment=environment_block())
        out = Path(args.out)
        if out.exists() and (not out.is_dir() or any(out.iterdir())):
            raise FileExistsError("output must be absent or empty; use a new run directory")
        out.mkdir(parents=True, exist_ok=True)
        table.to_csv(out / "scores.csv", index=False)
        (out / "comparison.json").write_text(json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")
        (out / "comparison.md").write_text(render_comparison(result), encoding="utf-8")
        print(out / "comparison.md")
        return 0
    except (ValueError, OSError) as exc:
        print(f"grnbench compare: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
