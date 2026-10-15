"""Rank aggregation of several methods: the "community" network.

Marbach et al., *Nat Methods* 9:796 (2012), "Wisdom of crowds for robust gene network
inference": averaging the edge ranks of many methods was more robust across data sets
than any single method. This adapter runs the member methods on the same normalised
data and orders every TF-gene pair by its mean rank (Borda count).

Members are explicit: the default set is fixed and a member that cannot run is an
error, not a silent omission, because a community built from whichever methods happened
to be installed would mean something different on every machine.
"""
from __future__ import annotations

import pandas as pd

DEFAULT_MEMBERS = ("spearman_signed", "clr", "ppcor", "lasso", "mrnet", "bicor")


def _members(value) -> list[str]:
    if value is None:
        return list(DEFAULT_MEMBERS)
    names = [m.strip() for m in value.split(",")] if isinstance(value, str) else [str(m).strip() for m in value]
    return [m for m in names if m]


def community_availability(members=None) -> tuple[bool, str]:
    """(available, reason): every member must be registered and runnable here."""
    from grnbench.infer.base import REGISTRY, unavailability_reason

    for name in _members(members):
        spec = REGISTRY.get(name)
        if spec is None:
            return False, f"member '{name}' is not registered"
        reason = unavailability_reason(spec)
        if reason is not None:
            return False, f"member '{name}' unavailable: {reason}"
    return True, ""


def community(adata, tf_list, members=None, seed=None, **_) -> pd.DataFrame:
    """Mean-rank aggregate of ``members`` (default :data:`DEFAULT_MEMBERS`).

    ``weight`` is ``1 - (mean_rank - 1) / n_pairs``, so 1 is the best possible
    aggregate rank. A pair a member does not report ranks last for that member: its rank is
    the size of the candidate universe (the union of every member's pairs), below every pair
    the member did report, so a sparse member cannot lift the pairs it leaves out.
    ``seed``, when given, is passed to every member (the default members are deterministic).
    """
    from grnbench.infer.base import MethodError, run_method

    names = _members(members)
    if not names:
        raise MethodError("community needs at least one member method")
    if "community" in names:
        raise MethodError("community cannot be a member of itself")
    ok, why = community_availability(names)
    if not ok:
        raise MethodError(f"community unavailable: {why}")

    extra = {} if seed is None else {"seed": seed}
    ranks = []
    for name in names:
        edges, _ = run_method(name, adata, tf_list=tf_list, assume_normalized=True, **extra)
        key = pd.MultiIndex.from_frame(edges[["source", "target"]])
        ranks.append(pd.Series(edges["weight"].to_numpy(), index=key)
                     .rank(ascending=False, method="average"))
    table = pd.concat(ranks, axis=1, keys=names)
    table = table.fillna(float(len(table)))      # a pair a member does not report ranks last for that member
    mean_rank = table.mean(axis=1)
    out = mean_rank.rename("mean_rank").reset_index()
    out.columns = ["source", "target", "mean_rank"]
    out["weight"] = 1.0 - (out["mean_rank"] - 1.0) / len(out)
    out.attrs["members"] = names
    return out[["source", "target", "weight"]]
