"""GRN inference adapters (Parte 2). Public API mirrors grnbench.strategies."""
from __future__ import annotations

from grnbench.infer.base import (
    EDGE_COLUMNS,
    MethodError,
    MethodSpec,
    get_method,
    list_methods,
    load_edges,
    prepare_adata,
    prepare_expression,
    register_method,
    run_method,
    save_edges,
    standardize_edges,
)

__all__ = [
    "EDGE_COLUMNS",
    "MethodError",
    "MethodSpec",
    "get_method",
    "list_methods",
    "load_edges",
    "prepare_adata",
    "prepare_expression",
    "register_method",
    "run_method",
    "save_edges",
    "standardize_edges",
]
