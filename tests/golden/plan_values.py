"""Golden values for the option-availability A/B (issue #128).

    PYTHONPATH=<checkout>/src python tests/golden/plan_values.py OUT.json

With ``kb._have`` fixed to the environment of the registered runs (``REGISTERED_ENV``: no R, no
scvi, no magic, no CellBender), so the values do not depend on the machine, it writes:

* ``option_table``: ``kb.option_table()`` (one row per option);
* ``plans``: for every simulator preset and the demo data, with and without a TF list, the whole
  ``plan.to_dict()`` (warnings, reasons, sources, selection, profile; ``option_notes`` left out)
  and the screening candidates ``select(profile).candidates(max_candidates=4, always=("context_scope",))``;
* ``variant_candidates``: ``benchmarks/candidate_study.variant_candidates()``.

``option_availability_main.json`` was written this way from a read-only copy of main 6917c88;
``tests/test_option_availability.py`` compares the current code with it. ``benchmarks/`` is found
next to ``src/`` (``PYTHONPATH=src`` or an editable install, not a site-packages install).
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

REGISTERED_ENV = frozenset({"py:scanpy", "py:decoupler"})


def _plain(obj):
    return json.loads(json.dumps(obj, default=str))


def datasets(workdir) -> dict:
    """``{name: (diag, tfs)}`` for every simulator preset and the demo data."""
    from grnbench.autoprep import diagnose
    from grnbench.demo import write_demo
    from grnbench.simulate import PRESETS, write_scenario

    work = Path(workdir)
    paths = {name: write_scenario(work / name, scn) for name, scn in PRESETS.items()}
    paths["demo"] = write_demo(work / "demo")
    return {name: (diagnose(p["counts"]), p["tfs"].read_text(encoding="utf-8").split())
            for name, p in paths.items()}


def plan_values(data: dict) -> dict:
    """``{"<dataset>/<tfs|no_tfs>": {"plan": plan.to_dict(), "candidates": [...]}}``."""
    from grnbench.autoprep import plan_preprocessing
    from grnbench.profile import build_profile
    from grnbench.selector import select

    out = {}
    for name, (diag, tfs) in data.items():
        for label, tf_list in (("tfs", tfs), ("no_tfs", None)):
            plan = plan_preprocessing(diag, tf_list=tf_list, methods=[])
            doc = plan.to_dict()
            doc.pop("option_notes", None)
            profile = build_profile(diag.adata, diag.qc.metrics, diag.context, id_verdict=diag.id_verdict,
                                    tf_list=tf_list, raw_droplets=diag.raw_droplets)
            cands = select(profile).candidates(max_candidates=4, always=("context_scope",))
            out[f"{name}/{label}"] = _plain({"plan": doc, "candidates": cands})
    return out


def option_table() -> list:
    from grnbench import kb

    return _plain(kb.option_table())


def variant_candidates() -> list:
    import grnbench

    benchmarks = str(Path(grnbench.__file__).resolve().parents[2] / "benchmarks")
    if benchmarks not in sys.path:
        sys.path.insert(0, benchmarks)
    import candidate_study

    return _plain(candidate_study.variant_candidates())


def main(out: str) -> None:
    from grnbench import kb

    kb._have = lambda req: req in REGISTERED_ENV
    with tempfile.TemporaryDirectory() as tmp:
        doc = {"environment": sorted(REGISTERED_ENV), "option_table": option_table(),
               "plans": plan_values(datasets(tmp)), "variant_candidates": variant_candidates()}
    Path(out).write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main(sys.argv[1])
