"""GRN-Bench: explicit preprocessing and evaluation contracts for GRN inference.

Four stages, one edge schema, full provenance:

    grnbench.pipeline    ingest -> QC gates -> strategy -> features -> processed.h5ad
    grnbench.infer       adapters (built-in + plugins) -> [source, target, weight]
    grnbench.eval        gold standards, metrics, nulls, subsample stability
    grnbench.html_report self-contained comparison report

Bring your own model with :mod:`grnbench.plugins` — a Python object exposing
``infer(adata, tf_list, **params)`` (declared as a ``grnbench.methods`` entry
point, or passed as ``--plugin module:Class``) or an external command / Docker
image reading ``expr.parquet`` and writing ``edges.parquet``.
"""
from __future__ import annotations

from grnbench.pipeline import PrepResult, run_prep
from grnbench.qc import Thresholds, run_qc
from grnbench.recommend import recommend_strategy
from grnbench.strategies import list_strategies

__version__ = "0.7.0"

__all__ = [
    "PrepResult",
    "Thresholds",
    "run_prep",
    "run_qc",
    "recommend_strategy",
    "list_strategies",
    "__version__",
]
