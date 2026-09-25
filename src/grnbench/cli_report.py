"""``grnbench report`` - LaTeX report generated from result files (see grnbench.texreport)."""
from __future__ import annotations

import argparse
import sys


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="grnbench report", description=__doc__)
    p.add_argument("target", nargs="?", help="a real_bench.py run directory or an experiment id (EXP-...)")
    p.add_argument("--exp", default=None, help="experiment id, when TARGET is a folder not named after it")
    p.add_argument("--out", default=None, help="output folder (default: RUN_DIR/latex, or paper/build)")
    p.add_argument("--paper", action="store_true",
                   help="assemble the confirmed/reproduced experiments of paper/CLAIMS.md into one draft")
    p.add_argument("--no-compile", action="store_true", help="write .tex only")
    a = p.parse_args(argv)
    from grnbench.texreport import build_paper, build_report

    try:
        if a.paper:
            main_tex = build_paper(out=a.out, compile_pdf=not a.no_compile)
        elif a.target:
            main_tex = build_report(a.target, exp=a.exp, out=a.out, compile_pdf=not a.no_compile)
        else:
            p.error("give a run directory / EXP id, or --paper")
    except (OSError, ValueError) as exc:
        print(f"grnbench report: {exc}", file=sys.stderr)
        return 2
    print(main_tex)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
