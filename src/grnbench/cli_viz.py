"""``grnbench viz`` - build figures + an index.md from a real_bench.py run directory."""
from __future__ import annotations

import argparse
import sys


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="grnbench viz", description=__doc__)
    p.add_argument("run_dir", help="a benchmarks/real_bench.py output directory (has rows/, "
                                   "optionally artifacts/ from --keep-artifacts)")
    p.add_argument("--out", default=None, help="output directory (default: RUN_DIR/viz)")
    p.add_argument("--tasks", default=None, help="tasks.json mapping task id -> raw data path, "
                                                  "for the QC/embeddings figures")
    p.add_argument("--max-datasets", type=int, default=None,
                   help="cap the number of per-dataset sections built")
    args = p.parse_args(argv)

    from grnbench.viz import build_run

    try:
        out_dir = build_run(args.run_dir, out_dir=args.out, tasks_json=args.tasks,
                            max_datasets=args.max_datasets)
    except (OSError, ValueError) as exc:
        print(f"grnbench viz: {exc}", file=sys.stderr)
        return 2
    print(out_dir / "index.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
