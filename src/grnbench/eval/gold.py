"""Gold-standard loading for GRN evaluation.

Two sources:
- a local file (``.tsv``/``.csv``/``.txt``) with regulator/target columns
  (header aliases matched case-insensitively);
- DoRothEA via decoupler (optional dependency), spec format ``dorothea<levels>``
  e.g. ``dorotheaABC`` (confidence A, B and C).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

# "Gene1,Gene2" is the BEELINE convention and, through it, the de-facto format of
# most published GRN ground truths (ChIP-seq, STRING and lof/gof networks alike).
# Accepting "gene2" as a target while omitting "gene1" as a source made every one
# of them unreadable, which is how this asymmetry was found.
_COL_ALIASES = {
    "source": {"source", "regulator", "tf", "tf1", "gene1", "g1", "from",
               "transcription_factor", "regulator_symbol"},
    "target": {"target", "target_gene", "gene2", "g2", "to", "target_symbol"},
    "confidence": {"confidence", "confidence_level", "grade"},
}


@dataclass(frozen=True)
class GoldStandard:
    name: str
    pairs: frozenset[tuple[str, str]]
    origin: str

    def __len__(self) -> int:
        return len(self.pairs)


def _canon(col: str) -> str | None:
    c = str(col).strip().lower()
    for canon, aliases in _COL_ALIASES.items():
        if c in aliases:
            return canon
    return None


def load_gold(spec: str) -> GoldStandard:
    """Load a gold standard from a path or a ``dorothea<levels>`` spec."""
    candidate = Path(spec)
    if candidate.exists():
        return _load_file(candidate)
    low = spec.lower()
    if low.startswith("dorothea"):
        return _load_dorothea(levels=spec[len("dorothea"):])
    raise FileNotFoundError(
        f"gold standard '{spec}': not a file and not a dorothea* spec"
    )


def _load_file(p: Path) -> GoldStandard:
    # decide the separator on the name with any .gz stripped, so "net.csv.gz" is
    # read as CSV rather than silently parsed as one wide tab-separated column
    stem = p.name[:-3] if p.name.lower().endswith(".gz") else p.name
    sep = "\t" if stem.lower().endswith((".tsv", ".txt")) or "tsv" in stem.lower() else ","
    df = pd.read_csv(p, sep=sep)
    mapping = {canon: col for col in df.columns if (canon := _canon(col))}
    if "source" not in mapping or "target" not in mapping:
        raise ValueError(
            f"{p.name}: need regulator + target columns "
            f"(accepted names: {sorted(_COL_ALIASES['source'])} / "
            f"{sorted(_COL_ALIASES['target'])}); found {list(df.columns)}"
        )
    pairs = frozenset(zip(df[mapping["source"]].astype(str),
                          df[mapping["target"]].astype(str)))
    return GoldStandard(name=p.stem, pairs=pairs, origin=str(p))


def _load_dorothea(levels: str = "ABC") -> GoldStandard:
    from grnbench.tfs import decoupler_network
    reg = decoupler_network("dorothea", organism="human", levels=levels or "ABC")
    cols = {c.lower(): c for c in reg.columns}
    src = cols.get("source", reg.columns[0])
    tgt = cols.get("target", reg.columns[1])
    pairs = frozenset(zip(reg[src].astype(str), reg[tgt].astype(str)))
    return GoldStandard(name=f"dorothea{levels}", pairs=pairs,
                        origin=f"decoupler dorothea(levels={levels!r})")


def coverage_report(edges_files: dict[str, pd.DataFrame], gold: GoldStandard,
                    case_sensitive: bool = False) -> dict:
    """How much of the gold can each method's candidate space even see?

    Case-folded by default, for the same reason scoring is: an expression matrix
    in MGI case against an upper-case ChIP-seq network reports "0 / 64557 gold
    reachable", which reads as a genuinely tiny candidate space rather than the
    identifier mismatch it actually is.
    """
    def norm(s: str) -> str:
        return s if case_sensitive else str(s).upper()

    report = {}
    gold_pairs = {(norm(s), norm(t)) for s, t in gold.pairs}
    gold_sources = {s for s, _ in gold_pairs}
    for name, edges in edges_files.items():
        candidates = {(norm(s), norm(t))
                      for s, t in zip(edges["source"], edges["target"])}
        report[name] = {
            "gold_pairs": len(gold),
            "gold_reachable": len(gold_pairs & candidates),
            "tf_coverage_in_data": len({s for s, _ in candidates} & gold_sources),
        }
    return report
