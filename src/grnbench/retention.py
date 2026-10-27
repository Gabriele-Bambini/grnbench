"""Retention ledger: what each cleaning step removed, with the regulators named.

Accounting only: ``run_prep`` keeps exactly the cells and genes it kept before. For every step
that can remove cells or genes the ledger records the counts before and after and, when a TF list
is given, how many TFs of that list were in the matrix before and after the step and which ones it
removed (names capped at :data:`TF_NAMES_CAP`, counts always complete, names in case-insensitive
alphabetical order).

D-003 protects TFs at feature selection only. A gene filter that runs earlier (minimum cells per
gene, QC-gene exclusion, the per-context gene filter) can drop a TF before the union is taken. A TF
it drops while the TF was *detected* (non-zero in at least one cell of the context, before the
step) is a lost regulator, and the step raises one warning naming it. Whether those filters should
keep detected TFs is a default change, not made here.

TF matching is case-insensitive (``str.casefold``), like :func:`grnbench.tfs.resolve_tf_names`,
but it never raises: zero matches and case collisions are reported by the code paths that already
handle them, so the ledger adds no new failure.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

TF_NAMES_CAP = 20          # names kept per ledger row; the counts are always complete
WARN_NAMES = 5             # names quoted in a warning
WARNING_PREFIX = "TF protection:"

#: The fields of one ledger row, in order (the ``tfs_*`` ones only when a TF list is given).
CELL_GENE_FIELDS = ("cells_before", "cells_after", "cells_removed",
                    "genes_before", "genes_after", "genes_removed", "aggregated")
TF_FIELDS = ("tfs_before", "tfs_after", "tfs_removed", "tfs_removed_detected",
             "tfs_removed_names", "tfs_removed_detected_names")


def sorted_names(names) -> list[str]:
    """Stable order for gene names: case-insensitive, ties broken by the exact spelling."""
    return sorted((str(n) for n in names), key=lambda n: (n.casefold(), n))


def _detected_cells(X, rows: np.ndarray | None, cols: np.ndarray) -> np.ndarray:
    """Per column of ``cols``: in how many of ``rows`` (all rows when None) it is non-zero."""
    sub = X[:, cols]                      # the few removed TF columns first, then the context rows
    if rows is not None:
        sub = sub[rows]
    if sp.issparse(sub):
        return np.asarray((sub != 0).sum(axis=0)).ravel().astype(int)
    return (np.asarray(sub) != 0).sum(axis=0).astype(int)


def _with_cells(name: str, n_cells: int) -> str:
    return f"{name} ({n_cells} cell{'' if n_cells == 1 else 's'})"


class RetentionLedger:
    """Collects one row per cell- or gene-removing step of ``run_prep``.

    ``subset_key``/``subset_value`` name the context; before the context subset a TF counts as
    detected only when it is non-zero in a cell of that context, after it every cell is the context.
    """

    def __init__(self, tf_list=None, *, subset_key=None, subset_value=None):
        self._tfs = (None if tf_list is None else
                     frozenset(str(t).strip().casefold() for t in tf_list if str(t).strip()))
        self._subset_key = subset_key
        self._subset_value = subset_value
        self.rows: list[dict] = []

    def _tf_names(self, var_names) -> list[str]:
        return [str(g) for g in var_names if str(g).casefold() in self._tfs]

    def _context_rows(self, adata) -> np.ndarray | None:
        key = self._subset_key
        if key is None or key not in adata.obs.columns:
            return None
        return adata.obs[key].astype(str).to_numpy() == str(self._subset_value)

    def record(self, step: str, before, after, *, aggregated: bool = False) -> list[str]:
        """Append the row of ``step`` (``before`` -> ``after``); return its warnings (0 or 1)."""
        cells_b, cells_a = int(before.n_obs), int(after.n_obs)
        genes_b, genes_a = int(before.n_vars), int(after.n_vars)
        row: dict = {
            "step": str(step),
            "cells_before": cells_b, "cells_after": cells_a,
            # aggregating cells into pseudobulk samples or metacells is not removing them
            "cells_removed": None if aggregated else cells_b - cells_a,
            "genes_before": genes_b, "genes_after": genes_a, "genes_removed": genes_b - genes_a,
            "aggregated": bool(aggregated),
        }
        warnings: list[str] = []
        if self._tfs is not None:
            tf_before = self._tf_names(before.var_names)
            tf_after = self._tf_names(after.var_names)
            kept = set(tf_after)
            removed = sorted_names(n for n in tf_before if n not in kept)
            detected: list[str] = []
            n_cells: dict[str, int] = {}          # detected TF -> cells of the context it is non-zero in
            if removed:
                cols = np.flatnonzero(np.asarray(before.var_names.isin(removed)))
                hits = _detected_cells(before.X, self._context_rows(before), cols)
                n_cells = {str(before.var_names[c]): int(h) for c, h in zip(cols, hits) if h > 0}
                detected = sorted_names(n_cells)
            row.update({
                "tfs_before": len(tf_before), "tfs_after": len(tf_after),
                "tfs_removed": len(removed), "tfs_removed_detected": len(detected),
                "tfs_removed_names": removed[:TF_NAMES_CAP],
                "tfs_removed_detected_names": detected[:TF_NAMES_CAP],
            })
            if detected:
                more = (f" and {len(detected) - WARN_NAMES} more (manifest 'retention')"
                        if len(detected) > WARN_NAMES else "")
                why = ("feature selection must keep every TF in the data (D-003)"
                       if step == "feature_selection" else
                       "the TF union of feature selection (D-003) comes later and cannot bring them back")
                warnings.append(
                    f"{WARNING_PREFIX} step '{step}' removed {len(detected)} TF(s) of the TF list that "
                    f"were detected in this context (non-zero in at least one cell before the step): "
                    f"{', '.join(_with_cells(n, n_cells[n]) for n in detected[:WARN_NAMES])}{more}. "
                    "They are missing from the processed "
                    f"matrix; {why}.")
        self.rows.append(row)
        return warnings
