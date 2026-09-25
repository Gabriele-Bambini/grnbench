"""Run a configured score-level calibration pilot."""
import argparse
import json
from pathlib import Path
import sys

from grnbench.calibration import run_calibration


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="grnbench calibrate", description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)
    try:
        config = json.loads(Path(args.config).read_text(encoding="utf-8"))
        result = run_calibration(config, args.out)
    except (ValueError, OSError, TypeError) as exc:
        print(f"grnbench calibrate: {exc}", file=sys.stderr)
        return 2
    for s in result["scenarios"]:
        print(f"{s['name']}: {s['gate']['status']}")
    print(Path(args.out) / "calibration.md")
    # A negative scientific result is not an execution failure.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
