"""Knockdown truth: score each knocked-down TF's predicted targets against its measured responders.

The unit is the TF (D-022, PROPOSED). Three pieces, all computed from (a) one edge list per model
(``source, target, weight``) and (b) a responders table (``tf, gene``) on a fixed gene universe:

* endpoint: for unit TF ``i``, rank every universe gene except ``i`` by the weight of edge
  ``i -> gene`` (missing edges share one bottom tie) and score that ranking against the responders
  of ``i`` (AUPRC ratio: expected AP over a random order inside each tied block; or precision at k),
  with the tie-safe metrics of ``grnbench.eval.metrics``;
* gene-propensity control: a ranking that carries no TF-specific information. For a network it is the
  mean within-regulator percentile of each gene over every regulator except the unit's own TF
  (leave-one-TF-out in-strength), so a network that puts the same hub genes on top for every TF
  scores the same as its own propensity. ``response_propensity`` does the same from the responders
  of the other units (how often a gene responds to any knockdown);
* dependence-keeping null: one TF-label permutation applied jointly to every unit and every network
  arm (TF ``pi(i)``'s ranking scored against TF ``i``'s responders); control arms stay on their own
  unit. It keeps each network whole (hub genes, correlated TF rows) and each responder set whole
  (shared responders), and breaks only the pairing of a TF with its own knockdown.

A declared hypothesis is a contrast over arms (``{arm: coefficient}``), e.g. ``{"m": 1,
"propensity:m": -1}``. Each is tested in both directions with the paired t of ``compare_units`` and
with the permutation null; Holm runs over all directed hypotheses of the family, separately for each
test. ``supported`` needs both adjusted p-values below alpha (``refuted``: both, other direction).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from grnbench.eval.comparison import compare_units, holm_adjust
from grnbench.eval.metrics import expected_average_precision, metric_from_labels
from grnbench.eval.universe import validate_edges

METRICS = ("auprc_ratio", "precision_at_k")
PROP = "propensity:"          # arm name prefix of a network's own propensity control
RESPONSE_PROP = "response_propensity"
GENETF = "_genetf"   # X_motif_genetf: X rank-producted with the TF-generic motif propensity (GB-1 gene_tf)


def _genes(genes) -> pd.Index:
    idx = pd.Index([str(g).upper() for g in genes])
    if idx.has_duplicates or len(idx) < 3:
        raise ValueError("gene universe must hold at least three distinct genes (case-insensitive)")
    return idx


def edge_matrix(edges: pd.DataFrame, sources, genes) -> np.ndarray:
    """Weights as a (sources x genes) array; a missing edge or a self pair is -inf (bottom tie)."""
    genes, sources = _genes(genes), pd.Index([str(s).upper() for s in sources])
    e = validate_edges(edges)
    mat = np.full((len(sources), len(genes)), -np.inf)
    r, c = sources.get_indexer(e.source), genes.get_indexer(e.target)
    keep = (r >= 0) & (c >= 0)
    mat[r[keep], c[keep]] = e.weight.to_numpy()[keep]
    return mat


def responder_matrix(responders: pd.DataFrame, units, genes) -> np.ndarray:
    """Boolean (units x genes); responders outside the universe are dropped, the TF's own gene too."""
    if not {"tf", "gene"} <= set(responders.columns):
        raise ValueError("responders table needs tf and gene columns")
    genes, units = _genes(genes), pd.Index([str(u).upper() for u in units])
    tf, gene = responders.tf.astype(str).str.upper(), responders.gene.astype(str).str.upper()
    r, c = units.get_indexer(tf), genes.get_indexer(gene)
    keep = (r >= 0) & (c >= 0) & (tf.to_numpy() != gene.to_numpy())
    mat = np.zeros((len(units), len(genes)), dtype=bool)
    mat[r[keep], c[keep]] = True
    return mat


def propensity(matrix: np.ndarray, sources, units) -> np.ndarray:
    """Leave-one-TF-out mean within-row percentile, one row per unit (the TF-generic ranking).

    Each row of ``matrix`` (a regulator's weights, or a unit's responder indicator) becomes average
    ranks scaled to (0, 1]; ties, including the -inf bottom tie, share their rank. Unit ``u`` gets the
    mean over every row except its own (when ``u`` is one of ``sources``).
    """
    sources = [str(s).upper() for s in sources]
    if len(sources) != matrix.shape[0] or len(sources) < 2:
        raise ValueError("propensity needs one row per source and at least two sources")
    pct = rankdata(np.where(np.isfinite(matrix), matrix, -np.finfo(float).max), axis=1) / matrix.shape[1]
    total, pos = pct.sum(axis=0), {s: i for i, s in enumerate(sources)}
    out = np.empty((len(units), matrix.shape[1]))
    for k, u in enumerate(str(u).upper() for u in units):
        j = pos.get(u)
        out[k] = total / len(sources) if j is None else (total - pct[j]) / (len(sources) - 1)
    return out


def response_propensity(resp: np.ndarray, units) -> np.ndarray:
    """How often each gene responds to the other units' knockdowns (leave-one-TF-out)."""
    return propensity(resp.astype(float), units, units)


def cross_scores(rows: np.ndarray, resp: np.ndarray, unit_gene: np.ndarray, metric: str = "auprc_ratio",
                 k: int = 50) -> np.ndarray:
    """C[j, i] = score of ranking row ``j`` on unit ``i``'s universe (all genes but ``i``'s own) and
    responders. The diagonal is the observed per-unit endpoint; the rest feeds the permutation null."""
    if metric not in METRICS:
        raise ValueError(f"metric must be one of {METRICS}")
    n = rows.shape[0]
    if resp.shape != rows.shape or len(unit_gene) != n:
        raise ValueError("rows, responders and unit genes must align")
    out = np.empty((n, n))
    for j in range(n):
        order = np.argsort(-rows[j], kind="stable")
        w = rows[j][order]
        for i in range(n):
            keep = order != unit_gene[i] if unit_gene[i] >= 0 else slice(None)
            y = resp[i][order][keep].astype(np.int8)
            if metric == "auprc_ratio":       # tied blocks (missing edges) are scored as a random order, in expectation
                out[j, i] = expected_average_precision(w[keep], y) / (y.sum() / len(y)) if y.any() else np.nan
            else:
                out[j, i] = metric_from_labels(metric, w[keep], y, k)
    return out


def is_control(arm: str) -> bool:
    """TF-generic control arms (a network's propensity, the response propensity): their row for a unit
    already leaves that unit out, so the permutation never moves them (a moved row would carry the
    unit's own TF back in)."""
    return arm.startswith(PROP) or arm == RESPONSE_PROP


def contrast_test(cross: dict[str, np.ndarray], coefs: dict[str, float], n_perm: int = 10000,
                  seed: int = 0) -> dict:
    """Paired t (``compare_units``) and joint TF-label permutation for one declared contrast.

    One permutation ``pi`` of the units is drawn per null sample and applied to every network arm
    at once (TF ``pi(i)``'s ranking against TF ``i``'s responders); control arms stay on their own
    unit. The statistic is the mean contrast over units."""
    if isinstance(n_perm, bool) or int(n_perm) != n_perm or n_perm < 1:
        raise ValueError("n_perm must be a positive integer")
    missing = set(coefs) - set(cross)
    if missing:
        raise ValueError(f"contrast uses arms that were not scored: {sorted(missing)}")
    n = next(iter(cross.values())).shape[0]
    diag = {a: np.diag(cross[a]) for a in coefs}
    pos = sum((c * diag[a] for a, c in coefs.items() if c > 0), np.zeros(n))
    neg = sum((-c * diag[a] for a, c in coefs.items() if c < 0), np.zeros(n))
    pair = compare_units(pd.DataFrame({"pos": pos, "neg": neg}))["pairs"][0]
    moved = sum((c * cross[a] for a, c in coefs.items() if not is_control(a)), np.zeros((n, n)))
    fixed = float(np.mean(sum((c * diag[a] for a, c in coefs.items() if is_control(a)), np.zeros(n))))
    observed = float(np.mean(pos - neg))
    rng = np.random.default_rng(seed)
    cols = np.arange(n)
    null = np.array([moved[rng.permutation(n), cols].mean() for _ in range(int(n_perm))]) + fixed
    tol = 1e-12 * max(1.0, abs(observed))
    return {"coefs": dict(coefs), "n_units": n, "mean_contrast": observed,
            "ci95": pair["simultaneous_ci"], "status": pair["status"],
            "p_t_greater": pair["p_a_gt_b"], "p_t_less": pair["p_b_gt_a"],
            "p_perm_greater": float((1 + np.sum(null >= observed - tol)) / (n_perm + 1)),
            "p_perm_less": float((1 + np.sum(null <= observed + tol)) / (n_perm + 1)),
            "null_mean": float(null.mean()), "null_sd": float(null.std()), "n_perm": int(n_perm), "seed": int(seed)}


def default_family(models) -> tuple[list, list]:
    """(family, secondary) as declared in the knockdown record.

    Family, for each ``X`` with an ``X_motif`` arm: H1a ``X_motif > X`` (the motif lift on knockdown
    truth); H1b ``X_motif > X_motif_genetf`` (the lift needs TF-specific motif identity, not only a
    gene propensity shared by all TFs; the ``gene_tf`` control of EXP-20261003-GB-1). H2, for every
    model except the ``_genetf`` controls: the model beats its own gene propensity. Secondary
    (reported, outside Holm): every model against the response propensity.
    """
    models = [m for m in models]
    base_of = lambda m: "pearson_signed" if m == "pearson_motif" else m.removesuffix("_motif")
    motif = [m for m in models if m.endswith("_motif") and base_of(m) in models]
    family = [(f"H1a {m} > {base_of(m)}", {m: 1, base_of(m): -1}) for m in motif]
    family += [(f"H1b {m} > {m}{GENETF}", {m: 1, m + GENETF: -1}) for m in motif if m + GENETF in models]
    real = [m for m in models if not m.endswith(GENETF)]
    family += [(f"H2 {m} > own propensity", {m: 1, PROP + m: -1}) for m in real]
    secondary = [(f"S {m} > response propensity", {m: 1, RESPONSE_PROP: -1}) for m in real]
    return family, secondary


def evaluate(edges_by_model: dict[str, pd.DataFrame], responders: pd.DataFrame, genes, *,
             metric: str = "auprc_ratio", k: int = 50, min_responders: int = 11, family=None, secondary=None,
             missing_arms=(), n_perm: int = 10000, seed: int = 0, alpha: float = 0.05) -> dict:
    """Score every model and control arm per unit TF and test the declared family.

    ``missing_arms``: models of the registered family that could not be built. A hypothesis that needs one stays
    in the family (Holm's multiplier does not shrink), with p = 1 in both directions and verdict ``missing``."""
    genes = _genes(genes)
    tfs = sorted(set(responders.tf.astype(str).str.upper()))
    units = [t for t in tfs if t in genes]
    resp = responder_matrix(responders, units, genes)
    n_resp = resp.sum(axis=1)
    excluded = {t: "TF not in the gene universe" for t in tfs if t not in genes}
    excluded.update({u: f"{int(c)} responders in the universe (< {min_responders})"
                     for u, c in zip(units, n_resp) if c < min_responders})
    keep = n_resp >= min_responders
    units, resp = [u for u, ok in zip(units, keep) if ok], resp[keep]
    if len(units) < 3:
        raise ValueError(f"only {len(units)} unit TFs left after exclusions; need at least three")
    unit_gene = genes.get_indexer(units)
    rows = {RESPONSE_PROP: response_propensity(resp, units)}
    for name, edges in edges_by_model.items():
        sources = sorted(set(validate_edges(edges).source))
        mat = edge_matrix(edges, sources, genes)
        own = pd.Index(sources).get_indexer(units)
        rows[name] = np.where(own[:, None] >= 0, mat[np.maximum(own, 0)], -np.inf)
        rows[PROP + name] = propensity(mat, sources, units)
    cross = {a: cross_scores(r, resp, unit_gene, metric, k) for a, r in rows.items()}
    if family is None and secondary is None:
        family, secondary = default_family(edges_by_model)
    gone = lambda c: {a.removeprefix(PROP) for a in set(c) - set(cross)}
    tests = [(h, c, contrast_test(cross, c, n_perm, seed) if not gone(c) or not gone(c) <= set(missing_arms) else
              {"coefs": dict(c), "n_units": len(units), "mean_contrast": float("nan"), "ci95": None, "status": "missing",
               "unavailable": sorted(gone(c)), "p_t_greater": 1.0, "p_t_less": 1.0, "p_perm_greater": 1.0,
               "p_perm_less": 1.0}) for h, c in (family or [])]
    t_adj = holm_adjust([p for _, _, r in tests for p in (r["p_t_greater"], r["p_t_less"])]) if tests else []
    perm_adj = holm_adjust([p for _, _, r in tests for p in (r["p_perm_greater"], r["p_perm_less"])]) if tests else []
    fam = []
    for i, (h, _, r) in enumerate(tests):
        r = {**r, "hypothesis": h, "p_t_greater_holm": t_adj[2 * i], "p_t_less_holm": t_adj[2 * i + 1],
             "p_perm_greater_holm": perm_adj[2 * i], "p_perm_less_holm": perm_adj[2 * i + 1]}
        r["verdict"] = ("missing" if r["status"] == "missing" else
                        "supported" if max(r["p_t_greater_holm"], r["p_perm_greater_holm"]) < alpha else
                        "refuted" if max(r["p_t_less_holm"], r["p_perm_less_holm"]) < alpha else "inconclusive")
        fam.append(r)
    sec = [{**contrast_test(cross, c, n_perm, seed), "hypothesis": h} for h, c in (secondary or [])]
    per_unit = pd.DataFrame({a: np.diag(c) for a, c in cross.items()}, index=pd.Index(units, name="tf"))
    per_unit.insert(0, "n_responders", n_resp[keep])
    return {"metric": metric, "k": k if metric == "precision_at_k" else None, "n_units": len(units),
            "n_genes": len(genes), "alpha": alpha, "n_directed_hypotheses": 2 * len(fam),
            "excluded": excluded, "per_unit": per_unit, "family": fam, "secondary": sec,
            "arm_means": per_unit.drop(columns="n_responders").mean().to_dict()}


def combine_sets(design: dict, replication: dict) -> dict[str, str]:
    """Final verdict per hypothesis: a claim needs the same verdict in both screens."""
    des = {r["hypothesis"]: r["verdict"] for r in design["family"]}
    rep = {r["hypothesis"]: r["verdict"] for r in replication["family"]}
    out = {}
    for h in [*des, *(h for h in rep if h not in des)]:
        a, b = des.get(h, "missing"), rep.get(h, "missing")
        out[h] = (f"{a} (replicated)" if a == b and a not in ("inconclusive", "missing") else
                  a if a == b else f"{a} in design, {b} in replication")
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--edges", nargs="+", required=True, help="model=path.tsv (source, target, weight)")
    ap.add_argument("--responders", required=True, type=Path, help="TSV with tf, gene columns")
    ap.add_argument("--genes", required=True, type=Path, help="gene universe, one per line")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--metric", default="auprc_ratio", choices=METRICS)
    ap.add_argument("--k", type=int, default=50)
    ap.add_argument("--min-responders", type=int, default=11)
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)
    edges = {}
    for spec in a.edges:
        name, _, path = spec.partition("=")
        if not path:
            ap.error(f"--edges takes model=path, got {spec!r}")
        edges[name] = pd.read_csv(path, sep="\t")
    genes = [g for g in Path(a.genes).read_text(encoding="utf-8").split() if g]
    res = evaluate(edges, pd.read_csv(a.responders, sep="\t"), genes, metric=a.metric, k=a.k,
                   min_responders=a.min_responders, n_perm=a.n_perm, seed=a.seed)
    a.out.mkdir(parents=True, exist_ok=True)
    res.pop("per_unit").to_csv(a.out / f"per_unit_{a.metric}.tsv", sep="\t")
    (a.out / f"family_{a.metric}.json").write_text(json.dumps(res, indent=1, allow_nan=True), encoding="utf-8")
    for r in res["family"]:
        print(f"{r['hypothesis']}: {r['mean_contrast']:+.4f}, p_t_holm {r['p_t_greater_holm']:.3g}, "
              f"p_perm_holm {r['p_perm_greater_holm']:.3g} -> {r['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
