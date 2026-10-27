"""``grnbench`` — umbrella command dispatching to the stage CLIs.

    grnbench run    --data d.h5ad --out runs/x --gold dorotheaABC
    grnbench prep   d.h5ad --out runs/x/prep
    grnbench infer  runs/x/prep/processed.h5ad --out runs/x/infer
    grnbench eval   runs/x/infer/*.edges.tsv.gz --gold g.tsv --out runs/x/eval
    grnbench methods | strategies | version

Each subcommand forwards argv untouched to the same ``main()`` the dedicated
``grnbench-*`` entry point uses, so the two spellings can never drift apart.
"""
from __future__ import annotations

import sys

_USAGE = """grnbench <command> [args]

Commands:
  run          preprocess + infer + evaluate in one go (grnbench-run)
  prep         preprocessing only                      (grnbench-prep)
  infer        inference only                          (grnbench-infer)
  eval         evaluation only                         (grnbench-eval)
  compare      paired scores across declared independent experimental units
  calibrate    run a configured score-level calibration pilot
  viz          build figures + an index.md from a real_bench.py run directory
  report       LaTeX report of one experiment (run dir or EXP id), or --paper for the draft
  board        what is running, open PRs, where things left off (run it first)
  methods      list inference methods and availability
  strategies   list preprocessing strategies
  kb           preprocessing knowledge base: options, availability, docs (--write-docs PATH)
  env          optional tools: R and its packages, Python and command-line extras (--json)
  demo         write the synthetic demo dataset (--out DIR)
  ui           open the graphical interface in the browser
  version      print the installed version

Run `grnbench <command> --help` for a command's own options."""


def _methods() -> int:
    from grnbench.infer import list_methods
    rows = list_methods()
    width = max((len(m["name"]) for m in rows), default=8)
    for m in rows:
        status = "available" if m["available"] else f"UNAVAILABLE ({m['unavailable_reason']})"
        if not m["default"]:
            status += "  [extra: runs with --methods all or by name]"
        print(f"{m['name']:{width}s}  {m['origin']:12s}  {status}")
        print(f"{'':{width}s}  {'':12s}  {m['description']}")
    return 0


def _strategies() -> int:
    from grnbench.strategies import list_strategies
    rows = list_strategies()
    width = max((len(s["name"]) for s in rows), default=8)
    for s in rows:
        status = "available" if s["available"] else f"UNAVAILABLE ({s['unavailable_reason']})"
        print(f"{s['name']:{width}s}  {status:52s}  {s['description']}")
    return 0


def _kb(rest: list[str]) -> int:
    import argparse
    from grnbench import kb
    ap = argparse.ArgumentParser(prog="grnbench kb", description="Print the preprocessing knowledge base.")
    ap.add_argument("--write-docs", nargs="?", const="docs/knowledge/preprocessing_kb.md", metavar="PATH",
                    help="write it as Markdown to PATH (default docs/knowledge/preprocessing_kb.md)")
    try:
        args = ap.parse_args(rest)
    except SystemExit as exc:      # argparse prints its own usage; hand back the code instead of exiting
        return int(exc.code or 0)
    if args.write_docs is not None:
        from pathlib import Path
        target = Path(args.write_docs)
        try:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(kb.render_markdown(), encoding="utf-8")
        except OSError as exc:     # '--write-docs ""' or a folder: say so instead of a traceback
            print(f"grnbench kb: cannot write {args.write_docs!r}: {exc}", file=sys.stderr)
            return 1
        print(f"wrote {target}")
        return 0
    for row in kb.option_table():
        status = "available" if row["available"] else f"UNAVAILABLE ({row['unavailable_reason']})"
        mark = "*" if row["default"] else " "
        print(f"{row['decision']:20s}{mark}{row['option']:26s} {row['evidence']:9s} {status}")
    return 0


_BUILTIN_HELP = {
    "methods": "grnbench methods\n\nList every inference method with its origin and whether it can run here.",
    "strategies": "grnbench strategies\n\nList the preprocessing strategies and whether each can run here.",
    "kb": "grnbench kb [--write-docs [PATH]]\n\nPrint the preprocessing knowledge base (decision, option, evidence level,\n"
          "availability), or write it as Markdown to PATH (default docs/knowledge/preprocessing_kb.md).",
    "env": "grnbench env [--rscript PATH] [--json]\n\nList the optional tools of the cleaning engine: R (found through --rscript,\n"
           "GRNBENCH_RSCRIPT, the conda env grnbench-r, R_HOME or PATH), each R package, Python module and\n"
           "command-line tool the knowledge base names, what needs it, and the command to install what is missing.\n"
           "Nothing is installed. GRNBENCH_RSCRIPT=none turns R off. Exits 1 when looking for R is an error (an R\n"
           "you named does not work, GRNBENCH_R_ENV names no env, an existing env has no working R), else 0.\n"
           "Details: docs/R_BRIDGE.md.",
    "demo": "grnbench demo [--out DIR]\n\nWrite a small synthetic dataset with a known network (counts.h5ad, tfs.txt,\n"
            "gold.tsv, README.md) into DIR (default grnbench_demo). DIR must be new or empty: nothing is overwritten.",
    "ui": "grnbench ui [streamlit options]\n\nOpen the graphical interface in the browser. It listens on this computer only;\n"
          "give --server.address 0.0.0.0 yourself to share it (anyone who can reach it can run code and read files).",
    "version": "grnbench version\n\nPrint the installed version.",
}


def _demo(rest: list[str]) -> int:
    import argparse
    from grnbench.demo import write_demo
    from grnbench.report import require_empty_output
    ap = argparse.ArgumentParser(prog="grnbench demo", description="Write the synthetic demo dataset.")
    ap.add_argument("--out", default="grnbench_demo", help="a new or empty folder (default grnbench_demo)")
    try:
        args = ap.parse_args(rest)
    except SystemExit as exc:      # argparse prints its own usage; hand back the code instead of exiting
        return int(exc.code or 0)
    try:
        # refuse before the 3 s build: never replace someone's files, nor keep a failed run's diagnostic
        require_empty_output(args.out, allow_failed_diagnostic=False)
    except (ValueError, OSError) as exc:
        print(f"grnbench demo: {exc}", file=sys.stderr)
        return 1
    for name, path in write_demo(args.out).items():
        print(f"{name:7s} {path}")
    return 0


def _env(rest: list[str]) -> int:
    import argparse
    import json
    from grnbench import rbridge
    ap = argparse.ArgumentParser(prog="grnbench env", description="List the optional tools and how to install them.")
    ap.add_argument("--rscript", metavar="PATH", help="use this Rscript instead of looking for one")
    ap.add_argument("--json", action="store_true", help="print the rows as JSON")
    try:
        args = ap.parse_args(rest)
    except SystemExit as exc:      # argparse prints its own usage; hand back the code instead of exiting
        return int(exc.code or 0)
    # A diagnostic: ask R, not the cache. R being absent (or turned off) is an answer, not a failure; looking
    # for R being an error (an R you named does not work, GRNBENCH_R_ENV names no env, an existing env has no
    # working R) is a failure, as it is for run_r, so scripts can check it.
    rows, failed = rbridge.environment_status(args.rscript, fresh=True)
    code = 1 if failed else 0
    if args.json:
        print(json.dumps(rows, indent=2))
        return code
    from grnbench import kb
    state = "on" if kb.r_bridge_on() else "off (GRNBENCH_R_BRIDGE=1 lets the planner see R packages)"
    print(f"R in planning: {state}")
    width = max(len(r["tool"]) for r in rows)
    print(f"{'tool':{width}s}  found  {'version':10s}  source")
    for r in rows:
        print(f"{r['tool']:{width}s}  {'yes' if r['found'] else 'no':5s}  {r['version'] or '-':10s}  {r['source'] or '-'}")
        if r["used_for"]:
            print(f"{'':{width}s}  used for: {', '.join(r['used_for'])}")
        for line in r["install"].splitlines() if r["install"] else ():
            print(f"{'':{width}s}  {'install: ' if not r['found'] and not line.startswith('#') else ''}{line}")
    return code


def _ui(rest: list[str]) -> int:
    import subprocess
    from pathlib import Path
    app = Path(__file__).with_name("ui") / "app.py"
    try:
        import streamlit  # noqa: F401
    except ImportError:
        print("the interface needs streamlit and plotly: from a clone of the repository run "
              "pip install -e '.[ui]'", file=sys.stderr)
        return 1
    from grnbench.ui import theme
    # The interface can run uploaded Python and read any folder of this computer: listen here only, unless asked.
    local_only = [] if any(a.startswith("--server.address") for a in rest) else ["--server.address", "localhost"]
    # The product look (accent, radii, chart colours) rides on the command line, so it needs no config file
    # next to the user's data; an option the user gives themselves wins.
    return subprocess.call([sys.executable, "-m", "streamlit", "run", str(app), *local_only,
                            *theme.streamlit_flags(rest), *rest])


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in {"-h", "--help", "help"}:
        print(_USAGE)
        return 0

    cmd, rest = argv[0], argv[1:]
    if cmd in _BUILTIN_HELP and rest[:1] in (["-h"], ["--help"]):
        print(_BUILTIN_HELP[cmd])
        return 0
    if cmd in {"methods", "strategies", "version"} and rest:
        print(f"grnbench {cmd}: takes no arguments (got {' '.join(rest)}); see `grnbench {cmd} --help`", file=sys.stderr)
        return 2
    if cmd in {"-V", "--version", "version"}:
        import grnbench
        print(f"grnbench {grnbench.__version__}")
        return 0
    if cmd == "methods":
        return _methods()
    if cmd == "strategies":
        return _strategies()
    if cmd == "kb":
        return _kb(rest)
    if cmd == "env":
        return _env(rest)
    if cmd == "demo":
        return _demo(rest)
    if cmd == "ui":
        return _ui(rest)

    dispatch = {
        "run": "grnbench.cli_run",
        "prep": "grnbench.cli",
        "infer": "grnbench.cli_infer",
        "eval": "grnbench.cli_eval",
        "compare": "grnbench.cli_compare",
        "calibrate": "grnbench.cli_calibrate",
        "viz": "grnbench.cli_viz",
        "report": "grnbench.cli_report",
        "board": "grnbench.board",
    }
    if cmd not in dispatch:
        print(f"grnbench: unknown command '{cmd}'\n", file=sys.stderr)
        print(_USAGE, file=sys.stderr)
        return 1
    import importlib
    return importlib.import_module(dispatch[cmd]).main(rest)


if __name__ == "__main__":
    raise SystemExit(main())
