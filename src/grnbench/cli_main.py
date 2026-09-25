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
    from grnbench import kb
    if rest[:1] == ["--write-docs"]:
        from pathlib import Path
        target = Path(rest[1] if len(rest) > 1 else "docs/knowledge/preprocessing_kb.md")
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(kb.render_markdown(), encoding="utf-8")
        print(f"wrote {target}")
        return 0
    for row in kb.option_table():
        status = "available" if row["available"] else f"UNAVAILABLE ({row['unavailable_reason']})"
        mark = "*" if row["default"] else " "
        print(f"{row['decision']:20s}{mark}{row['option']:26s} {row['evidence']:9s} {status}")
    return 0


def _demo(rest: list[str]) -> int:
    from grnbench.demo import write_demo
    out = rest[rest.index("--out") + 1] if "--out" in rest else "grnbench_demo"
    for name, path in write_demo(out).items():
        print(f"{name:7s} {path}")
    return 0


def _ui(rest: list[str]) -> int:
    import subprocess
    from pathlib import Path
    app = Path(__file__).with_name("ui") / "app.py"
    try:
        import streamlit  # noqa: F401
    except ImportError:
        print("the interface needs streamlit: pip install 'grnbench[ui]'", file=sys.stderr)
        return 1
    return subprocess.call([sys.executable, "-m", "streamlit", "run", str(app), *rest])


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in {"-h", "--help", "help"}:
        print(_USAGE)
        return 0

    cmd, rest = argv[0], argv[1:]
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
