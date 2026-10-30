"""Export pySCENIC's default modules in a small, portable JSON file.

This helper is intentionally executed by the isolated pySCENIC environment.  It does not
reimplement or patch pySCENIC: it only loads the same adjacency/expression files passed to
``pyscenic ctx``, calls the published ``modules_from_adjacencies`` defaults, and serialises the
module membership needed to diagnose any downstream disagreement.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from pyscenic.cli.utils import load_adjacencies, load_exp_matrix
from pyscenic.utils import modules_from_adjacencies


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--adjacencies", required=True)
    parser.add_argument("--expression", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)

    adjacencies = load_adjacencies(args.adjacencies)
    expression = load_exp_matrix(args.expression, False, False, "CellID", "Gene")
    modules = modules_from_adjacencies(adjacencies, expression)
    payload = []
    for module in modules:
        tf = str(module.transcription_factor)
        genes = sorted(str(g) for g in module.gene2weight)
        payload.append({"name": str(module.name), "tf": tf, "genes": genes})
    Path(args.out).write_text(json.dumps(payload, indent=1, sort_keys=True), encoding="utf-8")
    print(f"modules={len(payload)} out={args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
