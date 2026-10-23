"""GENIE3 / GRNBoost2 adapters over the arboreto package.

We drive ``arboreto.core.create_graph`` directly instead of the ``algo``
wrappers: this arboreto release builds a meta dask-dataframe from an empty
delayed list whenever ``include_meta=False`` (its default), which raises
"Must supply at least one delayed object" on modern dask-expr. Requesting the
meta output too avoids the dead branch; we discard it after compute.

A dashboard-less local dask client is injected via ``client_or_address``:
newer distributed/bokeh combos crash at teardown ("Cannot synchronously wait
on a running event loop") when a dashboard server was started.
"""
from __future__ import annotations

import importlib.util
from typing import Callable

import pandas as pd

ARBORETO_UNAVAILABLE = (
    None if importlib.util.find_spec("arboreto") is not None
    else "arboreto not installed - pip install 'grnbench[infer]'"
)


def _load():
    if ARBORETO_UNAVAILABLE is not None:
        from grnbench.infer.base import MethodError
        raise MethodError(f"unavailable: {ARBORETO_UNAVAILABLE}")
    from arboreto import core
    from arboreto.algo import RF_KWARGS, SGBM_KWARGS
    return core, dict(RF_KWARGS), dict(SGBM_KWARGS)


def _quiet_client(threads: int | None = None):
    # an unsized client takes every core on the machine, not the --threads budget
    sized = {} if threads is None else {"n_workers": 1, "threads_per_worker": int(threads)}
    # no quiet fallback: _run_arboreto needs a client, and returning None turned
    # a startup failure into an unrelated AttributeError on None.compute
    try:
        from dask.distributed import Client
        return Client(processes=False, dashboard_address=None, **sized)
    except Exception as e:
        raise RuntimeError(f"could not start a dask Client for arboreto (threads={threads})") from e


def _run_arboreto(expr: pd.DataFrame,
                  tf_names: list[str] | None,
                  seed: int,
                  *,
                  regressor_type: str,
                  early_stop_window_length: int,
                  n_estimators: int | None = None,
                  limit: int | None = None,
                  threads: int | None = None) -> pd.DataFrame:
    core, rf_kwargs, sgbm_kwargs = _load()
    kwargs = dict(rf_kwargs if regressor_type.upper() == "RF" else sgbm_kwargs)
    if n_estimators is not None:
        kwargs["n_estimators"] = int(n_estimators)

    matrix = expr.to_numpy(dtype="float32")
    gene_names = [str(g) for g in expr.columns]
    tf_list = [str(t) for t in tf_names] if tf_names else gene_names

    client = _quiet_client(threads)
    try:
        links_graph, _meta_graph = core.create_graph(
            expression_matrix=matrix,
            gene_names=gene_names,
            tf_names=tf_list,
            regressor_type=regressor_type,
            regressor_kwargs=kwargs,
            client=client,
            include_meta=True,  # True dodges arboreto's empty-meta from_delayed bug
            early_stop_window_length=int(early_stop_window_length),
            limit=limit,
            seed=int(seed) if seed is not None else core.DEMON_SEED,
        )
        links = client.compute(links_graph, sync=True)
    finally:
        client.close()
    return links.sort_values("importance", ascending=False).reset_index(drop=True)


def genie3_runner(expr: pd.DataFrame, tf_names: list[str] | None = None,
                  seed: int = 42, **params) -> pd.DataFrame:
    return _run_arboreto(
        expr, tf_names, seed,
        regressor_type="RF",
        early_stop_window_length=0,   # unused by RF profile (no OOB heuristic)
        n_estimators=params.pop("n_estimators", None),
        limit=params.pop("limit", None),
        threads=params.pop("threads", None),
    )


def grnboost2_runner(expr: pd.DataFrame, tf_names: list[str] | None = None,
                     seed: int = 42, **params) -> pd.DataFrame:
    return _run_arboreto(
        expr, tf_names, seed,
        regressor_type="GBM",
        early_stop_window_length=25,
        n_estimators=params.pop("n_estimators", None),
        limit=params.pop("limit", None),
        threads=params.pop("threads", None),
    )
