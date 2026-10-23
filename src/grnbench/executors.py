"""Which (decision, option) pairs of the knowledge base grnbench can run.

``PLAN_STEPS`` is the source of truth: the options :func:`grnbench.autoprep.plan_preprocessing`
maps onto :func:`grnbench.pipeline.run_prep`. :func:`executor_for` reads it and so does the plan,
so availability (:meth:`grnbench.kb.DecisionPoint.availability_full`) and the plan cannot drift
apart. An option grnbench cannot run is never chosen by the selector, and forcing one fails loudly
(issue #128). Step 2b adds R executors with :func:`register_executor`, which extends
``PLAN_STEPS``; each one also needs its mapping in the plan.

``RECORDED_ONLY`` are the options that run nothing: the plan records them and no step reads them.
The regulators options are among them: the TF list that runs is always the caller's (issue #131).

Rules every executor obeys:

* D-002: smoothing and imputation never become executable (``FORBIDDEN``);
* D-001: a matrix an executor returns is raw integer counts (:func:`check_integer_counts`).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
import scipy.sparse as sp

KINDS = ("python", "r", "cli")
# D-002: expression smoothed or imputed between cells (MAGIC, scVI's decoded expression).
FORBIDDEN = frozenset({("sparsity", "imputation"), ("batch", "scvi_denoised")})
FORBIDDEN_REASON = "not offered: D-002, no imputation or smoothing"
BUILT_IN = "built into grnbench: run by pipeline.run_prep as autoprep.plan_preprocessing maps it"

PLAN_STEPS: dict[str, tuple[str, ...]] = {
    "identifier_mapping": ("keep",), "cell_calling": ("as_provided", "knee"), "ambient_rna": ("none",),
    "cell_qc_thresholds": ("fixed", "mad", "mad_per_batch"), "doublets": ("none", "scrublet"),
    "gene_filter": ("min_cells", "min_cells_drop_qc_genes"),
    "context_scope": ("pooled", "per_cell_type", "selected_type"),
    "normalization": ("cpm_log1p", "pearson_residuals", "raw_counts"),
    "batch": ("none", "center_per_batch", "combat"), "confounders": ("none",),
    "sparsity": ("single_cells", "metacells", "pseudobulk"),
    "hvg": ("tf_plus_1000", "tf_plus_500", "tf_plus_2000", "all_genes"),
    "regulators": ("user_file", "dorothea", "all_genes"),
    "subsample": ("none", "stratified", "random"), "pseudotime": ("none", "provided"),
}
RECORDED_ONLY = frozenset({
    ("identifier_mapping", "keep"), ("cell_calling", "as_provided"), ("ambient_rna", "none"),
    ("confounders", "none"), ("pseudotime", "none"),
    ("pseudotime", "provided"),                         # no built-in model reads pseudotime
    ("regulators", "user_file"), ("regulators", "dorothea"), ("regulators", "all_genes"),   # #131
})


def _run_prep(*args, **kwargs):
    """The executor of every built-in step: ``pipeline.run_prep`` (imported only when it runs)."""
    from grnbench.pipeline import run_prep

    return run_prep(*args, **kwargs)


@dataclass(frozen=True)
class Executor:
    decision: str
    option: str
    kind: str                 # python | r | cli
    fn: Callable
    detail: str = ""

    def __post_init__(self):
        if self.kind not in KINDS:
            raise ValueError(f"executor kind {self.kind!r}; choose from {KINDS}")
        if not callable(self.fn):
            raise TypeError(f"executor for {self.decision}={self.option} must be callable, got {self.fn!r}")


_REGISTRY: dict[tuple[str, str], Executor] = {}


def _kb_option(decision: str, option: str):
    from grnbench.kb import DECISION_INDEX

    if decision not in DECISION_INDEX:
        raise KeyError(f"unknown decision {decision!r}")
    return DECISION_INDEX[decision].option(option)


def executor_for(decision: str, option: str) -> Executor | None:
    """The executor of ``decision=option``, or None when grnbench cannot run it.

    Raises KeyError for a pair the knowledge base does not know."""
    _kb_option(decision, option)
    if (decision, option) in FORBIDDEN:
        return None
    if (decision, option) in _REGISTRY:
        return _REGISTRY[(decision, option)]
    if option in PLAN_STEPS.get(decision, ()):
        return Executor(decision, option, "python", _run_prep, BUILT_IN)
    return None


def has_executor(decision: str, option: str) -> bool:
    return executor_for(decision, option) is not None


def register_executor(decision: str, option: str, fn: Callable, *, kind: str,
                      detail: str = "") -> Executor:
    """Let grnbench run ``decision=option`` with ``fn`` and add it to ``PLAN_STEPS``.

    Refuses smoothing/imputation (D-002), a pair the knowledge base does not know and a pair that
    already has an executor."""
    if (decision, option) in FORBIDDEN:
        raise ValueError(f"{decision}={option} smooths or imputes expression: grnbench never "
                         "executes it (D-002)")
    ex = Executor(decision, option, kind, fn, detail)
    if has_executor(decision, option):
        raise ValueError(f"{decision}={option} already has an executor")
    _REGISTRY[(decision, option)] = ex
    PLAN_STEPS[decision] = PLAN_STEPS.get(decision, ()) + (option,)
    return ex


def check_integer_counts(matrix, *, what: str = "executor output") -> None:
    """Refuse a matrix that is not raw integer counts (D-001): non-finite, negative or fractional.

    Stricter than ingest's ``looks_like_counts`` (which tolerates a few float-noise values in
    loaded files): an executor writes its own output, so every value must be an integer.
    """
    data = np.asarray(matrix.data if sp.issparse(matrix) else matrix, dtype=float).ravel()
    if not np.isfinite(data).all():
        raise ValueError(f"{what}: non-finite values; grnbench keeps raw integer counts only (D-001)")
    if (data < 0).any():
        raise ValueError(f"{what}: negative values; grnbench keeps raw integer counts only (D-001)")
    frac = np.abs(data - np.round(data)) > 1e-6
    if frac.any():
        raise ValueError(f"{what}: {int(frac.sum())} non-integer values; grnbench keeps raw integer "
                         "counts only (D-001)")
