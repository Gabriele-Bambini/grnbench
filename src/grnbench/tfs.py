"""Transcription-factor list utilities.

The platform never lets HVG selection silently drop TFs; these helpers load a
TF list and report its coverage in the dataset. DoRothEA fetching is optional
(via ``decoupler``).
"""
from __future__ import annotations

from collections.abc import Iterable, Sequence
from pathlib import Path

import pandas as pd

_TF_COL_HINTS = ("tf", "symbol", "gene", "name", "source")


class TFListError(ValueError):
    pass


def resolve_tf_names(
    gene_names: Sequence[object] | Iterable[object],
    tf_names: Iterable[object] | None,
) -> list[str]:
    """Resolve requested TFs to the original expression identifiers.

    Matching uses Unicode ``casefold`` so, for example, ``POU5F1`` resolves to
    the expression column ``Pou5f1``.  Returned names retain expression order
    and spelling.  Case-insensitive collisions and explicit TF selections that
    match nothing are errors because either would make the regulator universe
    ambiguous or silently expand it to all genes in legacy adapters.

    ``tf_names=None`` means every expression gene is a candidate regulator.
    """
    genes = [str(name) for name in gene_names]
    by_fold: dict[str, list[str]] = {}
    for gene in genes:
        by_fold.setdefault(gene.casefold(), []).append(gene)

    collisions = {fold: names for fold, names in by_fold.items() if len(names) > 1}
    if collisions:
        examples = ", ".join("/".join(names) for names in list(collisions.values())[:3])
        raise TFListError(
            "ambiguous expression identifiers under case-insensitive matching: "
            f"{examples}"
        )

    if tf_names is None:
        return genes
    requested = [str(name).strip() for name in tf_names]
    requested = [name for name in requested if name]
    if not requested:
        raise TFListError("empty TF list: provide at least one regulator identifier")

    wanted = {name.casefold() for name in requested}
    resolved = [gene for gene in genes if gene.casefold() in wanted]
    if not resolved:
        raise TFListError("none of the requested TF identifiers match expression genes")
    return resolved


def load_tf_list(path: str | Path) -> list[str]:
    """Load TF symbols from ``.txt`` (one per line) or ``.csv``/``.tsv``.

    For tabular input, the first column matching a known alias (tf/symbol/
    gene/name/source, case-insensitive) is used, else the first column.
    Returns uppercase, deduplicated symbols preserving order.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    sfx = p.suffix.lower()
    if sfx == ".txt":
        names = p.read_text(encoding="utf-8-sig").splitlines()
    elif sfx in {".csv", ".tsv"}:
        sep = "\t" if sfx == ".tsv" else ","
        df = pd.read_csv(p, sep=sep, dtype=str)
        col = next((c for c in df.columns if str(c).strip().lower() in _TF_COL_HINTS),
                   df.columns[0])
        names = df[col].tolist()
    else:
        raise TFListError(f"Unsupported TF list format: {sfx}")
    out: list[str] = []
    for n in names:
        s = str(n).strip().upper()
        if s and not s.startswith("#") and s not in out:
            out.append(s)
    if not out:
        raise TFListError(f"No TF symbols parsed from {p}")
    return out


def match_tfs_in_data(tf_list: list[str], var_names) -> tuple[list[str], int]:
    """Return (matched TFs present in data, total requested).

    Matching is case-insensitive; human symbols conventionally arrive already
    uppercased by :func:`load_tf_list`.
    """
    resolved = resolve_tf_names(var_names, tf_list)
    present = {name.casefold() for name in resolved}
    # Keep the requested spelling for the preprocessing feature-mask contract;
    # inference uses ``resolve_tf_names`` directly and gets expression spelling.
    matched = [str(tf).upper() for tf in tf_list if str(tf).casefold() in present]
    return matched, len(tf_list)


def decoupler_network(kind: str = "dorothea", organism: str = "human", levels: str = "ABC"):
    """DoRothEA or CollecTRI from decoupler, across its 1.x (``dc.get_*``) and
    2.x (``dc.op.*``) APIs. Returns the decoupler frame (source, target, ...)."""
    try:
        import decoupler as dc
    except ImportError as e:
        raise RuntimeError(
            "decoupler is required for DoRothEA/CollecTRI: "
            "pip install 'grnbench[dorothea]' or pip install decoupler"
        ) from e
    op = getattr(dc, "op", None)
    if kind == "dorothea":
        if op is not None and hasattr(op, "dorothea"):
            return op.dorothea(organism=organism, levels=list(levels or "ABC"))
        return dc.get_dorothea(organism=organism, levels=levels or "ABC")
    if kind == "collectri":
        if op is not None and hasattr(op, "collectri"):
            return op.collectri(organism=organism)
        return dc.get_collectri(organism=organism)
    raise ValueError(f"unknown decoupler network {kind!r}")


def fetch_dorothea(organism: str = "human", levels: str = "ABC"):
    """Fetch DoRothEA regulons via decoupler (optional dependency).

    Returns a DataFrame with columns ``[source(set), target, confidence]``.
    """
    reg = decoupler_network("dorothea", organism=organism, levels=levels)
    cols = {c.lower(): c for c in reg.columns}
    src = cols.get("source") or reg.columns[0]
    tgt = cols.get("target") or reg.columns[1]
    conf = cols.get("confidence") or cols.get("confidence_level")
    out = reg[[src, tgt] + ([conf] if conf else [])].copy()
    out.columns = ["source", "target"] + (["confidence"] if conf else [])
    return out
