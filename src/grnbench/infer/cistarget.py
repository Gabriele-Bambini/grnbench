"""cisTarget motif pruning for any co-expression method (pySCENIC step 2).

Co-expression alone cannot tell a direct target from a gene that merely moves
with the regulator. SCENIC (Aibar et al., *Nat Methods* 14:1083, 2017; Van de
Sande et al., *Nat Protoc* 15:2247, 2020) keeps a TF -> target edge only when
the TF's binding motif is enriched near the co-expressed targets. This module
implements that step with numpy on the public cisTarget resources, so any
co-expression adapter (``pearson_signed``, ``grnboost2``, ``regdiffusion``)
can be motif-pruned:

1. **Modules** (``modules_from_edges``, pySCENIC defaults): per TF, activating
   targets only (Pearson rho > 0.03); top 50 targets, targets above the 75th
   and 90th weight percentile, targets for which the TF is among their top
   5/10/50 regulators; modules with fewer than 20 genes are dropped.
2. **Enrichment** (``prune``): for every module and every motif of the ranking
   database, the area under the recovery curve of the module genes in the
   motif's genome-wide ranking up to 5% of the genome (``auc_threshold``);
   NES = (AUC - mean) / sd over all motifs. Motifs annotated to the module's TF
   with NES >= 3 are enriched.
3. **Leading edge**: module genes ranked before the point where the motif's
   recovery curve exceeds the mean + 2 sd curve of all motifs the most
   (``rank_threshold`` = 5000, computed on a 5-rank grid).

Output (``cistarget_rerank``): the base method's full ranking, re-ordered so
that motif-supported edges come first (each block keeps the base order). The
evaluation universe therefore stays the same; nothing is dropped.

A second, parameter-free use of the same resources (``motif_prior``,
``motif_rank_product``): score every TF x gene pair by the best genome-wide
rank of the gene among the TF's motifs, and combine it with a co-expression
ranking by rank product (geometric mean of the two ranks). ``motif_prior``
alone is the sequence-only control.

Prior knowledge used: genome sequence (motif rankings) and motif -> TF
annotations (motif collections, motif similarity, orthology). No curated
TF -> target collection (CollecTRI, DoRothEA, TRRUST) enters here, so scoring
against CollecTRI stays non-circular.

Resources: ``GRNBENCH_CISTARGET_DIR`` must hold at least one
``*.genes_vs_motifs.rankings.feather`` (hg38_* for human, mm10_* for mouse)
and the matching motif annotation ``*.hgnc-*.tbl`` / ``*.mgi-*.tbl`` from
https://resources.aertslab.org/cistarget/ (v10 clustered motifs used here).
"""
from __future__ import annotations

import os
import warnings
from functools import lru_cache
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

ENV_DIR = "GRNBENCH_CISTARGET_DIR"
SPECIES_PREFIX = {"human": ("hg38", "hg19", "hgnc"), "mouse": ("mm10", "mm9", "mgi")}


# --------------------------------------------------------------------------- #
# resources
# --------------------------------------------------------------------------- #

def cistarget_dir() -> Path | None:
    raw = os.environ.get(ENV_DIR, "")
    return Path(raw) if raw and Path(raw).is_dir() else None


def resource_files(species: str, root: Path | None = None) -> tuple[list[Path], Path | None]:
    """(ranking databases, motif annotation table) for one species."""
    root = root or cistarget_dir()
    if root is None:
        return [], None
    *db_prefixes, ann_tag = SPECIES_PREFIX[species]
    dbs = sorted(p for p in root.glob("*.genes_vs_motifs.rankings.feather")
                 if p.name.startswith(tuple(db_prefixes)))
    ann = next(iter(sorted(root.glob(f"*.{ann_tag}-*.tbl"))), None)
    return dbs, ann


def cistarget_availability() -> tuple[bool, str]:
    root = cistarget_dir()
    if root is None:
        return False, (f"{ENV_DIR} not set to a folder with cisTarget rankings (*.rankings.feather) "
                       "and motif annotations (*.tbl) from resources.aertslab.org/cistarget")
    try:
        import pyarrow  # noqa: F401
    except ImportError:
        return False, "pyarrow not installed - pip install pyarrow"
    usable = [s for s in SPECIES_PREFIX if all(resource_files(s, root))]
    if not usable:
        return False, f"{root}: no species with both a rankings .feather and a motif .tbl"
    return True, ""


@lru_cache(maxsize=4)
def load_annotation(path: str) -> dict[str, frozenset]:
    """Motif annotation table -> {TF symbol (casefolded): motif ids}."""
    ann = pd.read_csv(path, sep="\t", usecols=["#motif_id", "gene_name"])
    out: dict[str, set] = {}
    for motif, tf in zip(ann["#motif_id"].astype(str), ann["gene_name"].astype(str)):
        out.setdefault(tf.casefold(), set()).add(motif)
    return {k: frozenset(v) for k, v in out.items()}


def _schema_names(path: str) -> list[str]:
    import pyarrow.ipc as ipc
    with open(path, "rb") as fh:
        return list(ipc.open_file(fh).schema.names)


def load_rankings(path: str, genes: tuple[str, ...]) -> tuple[np.ndarray, list[str], np.ndarray, int]:
    """(motif ids, matched gene names in ``genes`` spelling, ranks motifs x genes, total genes in DB)."""
    import pyarrow.feather as pf
    cols = [c for c in _schema_names(path) if c != "motifs"]
    by_fold = {c.casefold(): c for c in cols}
    matched = [(g, by_fold[g.casefold()]) for g in genes if g.casefold() in by_fold]
    table = pf.read_table(path, columns=["motifs"] + [db for _, db in matched], memory_map=True)
    motifs = np.asarray(table.column("motifs").to_pylist(), dtype=object)
    ranks = np.column_stack([table.column(db).to_numpy() for _, db in matched]).astype(np.int32) \
        if matched else np.zeros((len(motifs), 0), dtype=np.int32)
    return motifs, [g for g, _ in matched], ranks, len(cols)


def coverage_report(regulators, genes, found, motif_ids, tf_motifs, dbs) -> dict:
    """Share of ``genes`` present in the rankings databases and of ``regulators`` with an annotated
    motif present in them. Warns (scores untouched) when identifiers probably do not match."""
    regs = list(dict.fromkeys(map(str, regulators)))
    with_motif = sum(1 for t in regs if tf_motifs.get(t.casefold(), frozenset()) & motif_ids)
    n_genes = len(set(map(str, genes)))
    cov = {"genes_total": n_genes, "genes_in_db": len(found),
           "gene_coverage": round(len(found) / n_genes, 4) if n_genes else 0.0,
           "regulators_total": len(regs), "regulators_with_motif": with_motif,
           "regulator_coverage": round(with_motif / len(regs), 4) if regs else 0.0,
           "databases": [Path(d).name for d in dbs]}
    if cov["gene_coverage"] < 0.5 or with_motif == 0:
        warnings.warn(
            f"cisTarget coverage is low: {len(found)}/{n_genes} genes ({cov['gene_coverage']:.0%}) are in "
            f"the rankings and {with_motif}/{len(regs)} regulators have an annotated motif present; the motif "
            "prior is constant or neutral for the rest. Likely cause: gene identifiers do not match the "
            "cisTarget symbols (Ensembl IDs?)", RuntimeWarning, stacklevel=3)
    return cov


def guess_species(genes) -> str:
    """Human symbols are upper case (TP53), mouse symbols title case (Trp53)."""
    alpha = [g for g in map(str, genes) if any(c.isalpha() for c in g)]
    if not alpha:
        return "human"
    upper = sum(g == g.upper() for g in alpha) / len(alpha)
    return "human" if upper >= 0.5 else "mouse"


# --------------------------------------------------------------------------- #
# algorithm
# --------------------------------------------------------------------------- #

def _pearson_for_edges(expr: pd.DataFrame, edges: pd.DataFrame) -> np.ndarray:
    X = expr.to_numpy(dtype=np.float64)
    X = X - X.mean(0)
    sd = np.sqrt((X ** 2).sum(0))
    sd[sd == 0] = np.inf
    Z = X / sd
    col = {g: i for i, g in enumerate(map(str, expr.columns))}
    s = edges["source"].map(col).to_numpy()
    t = edges["target"].map(col).to_numpy()
    ok = ~(pd.isna(s) | pd.isna(t))
    rho = np.full(len(edges), np.nan)
    si, ti = s[ok].astype(int), t[ok].astype(int)
    regs, local = np.unique(si, return_inverse=True)
    corr = Z[:, regs].T @ Z                     # regulators x genes, never edges x cells
    rho[ok] = corr[local, ti]
    return rho


def modules_from_edges(edges: pd.DataFrame, expr: pd.DataFrame, *,
                       top_n_targets=(50,), percentiles=(0.75, 0.90),
                       top_n_regulators=(5, 10, 50), min_genes: int = 20,
                       rho_threshold: float = 0.03) -> list[tuple[str, frozenset]]:
    """pySCENIC ``modules_from_adjacencies`` defaults, activating modules only."""
    e = edges[["source", "target", "weight"]].copy()
    e = e[e["source"] != e["target"]]
    e["rho"] = _pearson_for_edges(expr, e)
    e = e[e["rho"] > rho_threshold]
    if e.empty:
        return []
    modules = []
    e = e.sort_values("weight", ascending=False, kind="stable")
    for tf, g in e.groupby("source", sort=True):
        for n in top_n_targets:
            modules.append((tf, frozenset(g["target"].iloc[:n])))
        for q in percentiles:
            modules.append((tf, frozenset(g.loc[g["weight"] >= g["weight"].quantile(q), "target"])))
    rank_in_target = e.groupby("target")["weight"].rank(method="first", ascending=False)
    for n in top_n_regulators:
        sub = e[rank_in_target <= n]
        for tf, g in sub.groupby("source", sort=True):
            modules.append((tf, frozenset(g["target"])))
    seen, out = set(), []
    for tf, genes in modules:
        if len(genes) >= min_genes and (tf, genes) not in seen:
            seen.add((tf, genes))
            out.append((tf, genes))
    return out


def prune(modules, motifs: np.ndarray, genes: list[str], ranks: np.ndarray, n_total: int,
          tf_motifs: dict[str, frozenset], *, rank_threshold: int = 5000,
          auc_threshold: float = 0.05, nes_threshold: float = 3.0,
          grid_step: int = 5, min_present: int = 10) -> dict[tuple[str, str], float]:
    """Motif-supported (TF, target) pairs with their best NES."""
    col = {g: i for i, g in enumerate(genes)}
    motif_row = {m: i for i, m in enumerate(motifs)}
    R = max(1, int(round(auc_threshold * n_total)))
    n_bins = max(1, rank_threshold // grid_step)
    support: dict[tuple[str, str], float] = {}
    for tf, module in modules:
        rows = [motif_row[m] for m in tf_motifs.get(str(tf).casefold(), ()) if m in motif_row]
        idx = np.array([col[g] for g in module if g in col], dtype=int)
        if not rows or len(idx) < min_present:
            continue
        sub = ranks[:, idx]
        auc = np.clip(R - sub, 0, None).sum(1, dtype=np.float64) / (R * len(idx))
        sd = auc.std()
        if sd == 0:
            continue
        nes = (auc - auc.mean()) / sd
        enriched = [r for r in rows if nes[r] >= nes_threshold]
        if not enriched:
            continue
        # recovery curves of all motifs on a coarse grid: rec[m, b] = #genes ranked < (b+1)*step
        b = np.minimum(sub // grid_step, n_bins)
        flat = (np.arange(sub.shape[0])[:, None] * (n_bins + 1) + b).ravel()
        hist = np.bincount(flat, minlength=sub.shape[0] * (n_bins + 1)).reshape(sub.shape[0], n_bins + 1)
        rec = np.cumsum(hist[:, :n_bins], axis=1, dtype=np.float64)
        thr = rec.mean(0) + 2 * rec.std(0)
        module_genes = np.array([genes[i] for i in idx], dtype=object)
        for r in enriched:
            crit = (int(np.argmax(rec[r] - thr)) + 1) * grid_step
            for g in module_genes[sub[r] < crit]:
                key = (str(tf), str(g))
                support[key] = max(support.get(key, -np.inf), float(nes[r]))
    return support


def cistarget_rerank(edges: pd.DataFrame, expr: pd.DataFrame, species: str | None = None,
                     root: Path | None = None, **params) -> pd.DataFrame:
    """Re-rank ``edges``: motif-supported edges first, base order kept inside each block."""
    species = species or guess_species(expr.columns)
    dbs, ann = resource_files(species, root)
    if not dbs or ann is None:
        from grnbench.infer.base import MethodError
        raise MethodError(f"no cisTarget resources for {species} in {root or cistarget_dir()}")
    tf_motifs = load_annotation(str(ann))
    modules = modules_from_edges(edges, expr)
    support: dict[tuple[str, str], float] = {}
    genes = tuple(map(str, expr.columns))
    found, motif_ids = set(), set()
    for db in dbs:
        motifs, matched, ranks, n_total = load_rankings(str(db), genes)
        found.update(matched)
        motif_ids.update(motifs)
        for k, v in prune(modules, motifs, matched, ranks, n_total, tf_motifs, **params).items():
            support[k] = max(support.get(k, -np.inf), v)
    out = edges.sort_values("weight", ascending=False, kind="stable").reset_index(drop=True)
    base_score = 1.0 - np.arange(len(out)) / max(len(out), 1)          # (0, 1], base order
    supported = np.fromiter(((s, t) in support for s, t in zip(out["source"].astype(str),
                                                                out["target"].astype(str))),
                            dtype=bool, count=len(out))
    out["weight"] = base_score + supported.astype(float)
    out = out.sort_values("weight", ascending=False, kind="stable").reset_index(drop=True)
    out.attrs.update(cistarget_species=species, cistarget_supported=int(supported.sum()),
                     cistarget_modules=len(modules),
                     cistarget_regulons=len({s for s, _ in support}),
                     motif_coverage=coverage_report(out["source"], genes, found, motif_ids, tf_motifs, dbs))
    return out


def with_cistarget(base: Callable[..., pd.DataFrame]) -> Callable[..., pd.DataFrame]:
    """Wrap an ``expr`` runner: base co-expression edges, then motif re-ranking."""
    def runner(expr: pd.DataFrame, tf_names: list[str] | None = None, species: str | None = None,
               **params) -> pd.DataFrame:
        edges = base(expr, tf_names=tf_names, **params)
        edges = edges.rename(columns=dict(zip(edges.columns[:3], ["source", "target", "weight"])))
        return cistarget_rerank(edges, expr, species=species)
    runner.__name__ = f"{getattr(base, '__name__', 'base')}_cistarget"
    return runner


# --------------------------------------------------------------------------- #
# motif prior: sequence evidence as a ranking, and its rank-product with a base
# --------------------------------------------------------------------------- #

def motif_ranks(tfs, genes, species: str | None = None, root: Path | None = None) -> pd.DataFrame:
    """TFs x genes: best (lowest) genome-wide rank of ``gene`` over the TF's motifs, over all
    ranking databases of the species. TFs without an annotated motif and genes missing from
    the databases get the worst rank (the database size)."""
    genes = tuple(map(str, genes))
    species = species or guess_species(genes)
    dbs, ann = resource_files(species, root)
    if not dbs or ann is None:
        from grnbench.infer.base import MethodError
        raise MethodError(f"no cisTarget resources for {species} in {root or cistarget_dir()}")
    tf_motifs = load_annotation(str(ann))
    out = None
    found, motif_ids = set(), set()
    for db in dbs:
        motifs, matched, ranks, n_total = load_rankings(str(db), genes)
        found.update(matched)
        motif_ids.update(motifs)
        row = {m: i for i, m in enumerate(motifs)}
        block = np.full((len(tfs), len(genes)), n_total, dtype=np.int32)
        pos = {g: i for i, g in enumerate(genes)}
        cols = [pos[g] for g in matched]
        for i, tf in enumerate(tfs):
            rr = [row[m] for m in tf_motifs.get(str(tf).casefold(), ()) if m in row]
            if rr and cols:
                block[i, cols] = ranks[rr].min(0)
        out = block if out is None else np.minimum(out, block)
    res = pd.DataFrame(out, index=list(map(str, tfs)), columns=list(genes))
    res.attrs["motif_coverage"] = coverage_report(tfs, genes, found, motif_ids, tf_motifs, dbs)
    return res


def _motif_rank_of_edges(edges: pd.DataFrame, expr: pd.DataFrame, species, root) -> np.ndarray:
    tfs = list(dict.fromkeys(edges["source"].astype(str)))
    mr = motif_ranks(tfs, expr.columns, species, root)
    if "motif_coverage" in mr.attrs:            # ranks a caller built itself carry no coverage report
        edges.attrs["motif_coverage"] = mr.attrs["motif_coverage"]
    r = mr.to_numpy()
    ri = edges["source"].astype(str).map({t: i for i, t in enumerate(mr.index)}).to_numpy()
    ci = edges["target"].astype(str).map({g: i for i, g in enumerate(mr.columns)}).to_numpy()
    return r[ri.astype(int), ci.astype(int)].astype(float)


def motif_prior(expr: pd.DataFrame, tf_names: list[str] | None = None, species: str | None = None,
                root: Path | None = None, **_) -> pd.DataFrame:
    """Sequence-only control: every TF x gene pair scored by its motif rank (no expression used
    beyond the gene list). Tells how much of a combined score the data themselves add."""
    from grnbench.tfs import resolve_tf_names
    tfs = resolve_tf_names(expr.columns, tf_names)
    mr = motif_ranks(tfs, expr.columns, species, root)
    long = mr.stack().rename_axis(["source", "target"]).reset_index(name="rank")
    long = long[long["source"] != long["target"]]
    out = pd.DataFrame({"source": long["source"], "target": long["target"],
                        "weight": -long["rank"].astype(float)})
    if "motif_coverage" in mr.attrs:
        out.attrs["motif_coverage"] = mr.attrs["motif_coverage"]
    return out


def motif_rank_product(edges: pd.DataFrame, expr: pd.DataFrame, species: str | None = None,
                       root: Path | None = None) -> pd.DataFrame:
    """Parameter-free combination: geometric mean of the base rank and the motif rank
    (rank product, Breitling et al. 2004). Ties share their average rank."""
    out = edges.copy()
    base_rank = out["weight"].rank(ascending=False, method="average").to_numpy()
    mot_rank = pd.Series(_motif_rank_of_edges(out, expr, species, root)).rank(
        ascending=True, method="average").to_numpy()
    out["weight"] = -(np.log(base_rank) + np.log(mot_rank)) / 2
    cov = out.attrs.get("motif_coverage")
    out = out.sort_values("weight", ascending=False, kind="stable").reset_index(drop=True)
    if cov is not None:
        out.attrs["motif_coverage"] = cov
    return out


def with_motif_prior(base: Callable[..., pd.DataFrame]) -> Callable[..., pd.DataFrame]:
    """Wrap an ``expr`` runner: base co-expression edges, then rank product with the motif prior."""
    def runner(expr: pd.DataFrame, tf_names: list[str] | None = None, species: str | None = None,
               **params) -> pd.DataFrame:
        edges = base(expr, tf_names=tf_names, **params)
        edges = edges.rename(columns=dict(zip(edges.columns[:3], ["source", "target", "weight"])))
        edges = edges[edges["source"] != edges["target"]]
        return motif_rank_product(edges, expr, species=species)
    runner.__name__ = f"{getattr(base, '__name__', 'base')}_motif"
    return runner
