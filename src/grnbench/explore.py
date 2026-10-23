"""EXPLORATORY descriptive exploration + per-TF explainability over saved run artifacts.

Reads an artifacts folder written by ``benchmarks/real_bench.py --keep-artifacts`` (one folder per dataset:
``reference_evaluation.tsv`` plus ``<arm>/ctx0/{processed.h5ad, <model>.tsv.gz}``) and writes, all labelled
EXPLORATORY (descriptive, no claim, no test):

* ``tf_atlas.tsv``: per (dataset, TF) detection rate, mean expression, reference targets, per-model per-TF AUPRC ratio;
* ``fig_tf_detectability``: does the motif prior help most for TFs whose mRNA is barely detected? (ranked-in-both pairs,
  coverage, mean lift with dataset-cluster bootstrap);
* ``fig_promoted_edges``: which top-1000 edges the prior keeps, promotes and demotes, and how many are reference edges;
* ``card_<TF>.md`` / ``.png``: one TF's evidence across datasets;
* ``README.md``: the question, Spearman rho per family, medians.

    python -m grnbench.explore --artifacts DIR --out OUT [--arm auto] [--tf TP53]

Scoring follows ``grnbench.eval``: edges are ranked by descending ``weight`` (signed, as stored) and scored with
``evaluate_ranking`` on an explicit universe (the TF against every other gene of the dataset); edges missing from
the stored list (lists are truncated by the run) share one bottom tie. A TF absent from a model's list gets NaN.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

from grnbench.eval.metrics import GoldIndex, evaluate_ranking  # noqa: E402

LABEL = "EXPLORATORY (descriptive, no claim)"
FAMILIES = {"pearson": ("pearson_signed", "pearson_motif"),
            "grnboost2": ("grnboost2", "grnboost2_motif"),
            "regdiffusion": ("regdiffusion", "regdiffusion_motif")}
MODELS = [m for pair in FAMILIES.values() for m in pair]
FAMILY_LABEL = {"pearson": "Pearson", "grnboost2": "GRNBoost2", "regdiffusion": "RegDiffusion"}
OKABE_ITO = {"pearson": "#0072B2", "grnboost2": "#E69F00", "regdiffusion": "#009E73"}
VERMILLION, GREY = "#D55E00", "#999999"
MIN_TARGETS = 5
TOP_FRAC = 0.01
SEED = 0
PROMO_K = 1000
MIN_PRESENT = 3
N_BOOT = 1000
N_BINS = 5
DET_FLOOR = 1e-3  # detection 0 is drawn at this value on the log axis only
QUESTION = "Does the motif prior help most for TFs whose mRNA is barely detected?"
MM = 1 / 25.4
FIG_STYLE = {"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"], "font.size": 7,
             "axes.labelsize": 7, "axes.titlesize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
             "legend.fontsize": 7, "axes.linewidth": 0.6, "lines.linewidth": 1.0, "axes.grid": False,
             "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42,
             "savefig.dpi": 300, "savefig.bbox": "standard", "figure.dpi": 100}


def _short(name: str) -> str:
    p = name.split("__")
    return f"{p[3]}/{p[4][:4]}" if len(p) >= 5 else name


# ------------------------------------------------------------------------------------------ loading
def load_dataset(ddir: Path, arm: str = "auto") -> dict:
    import anndata as ad

    ctx = Path(ddir) / arm / "ctx0"
    adata = ad.read_h5ad(ctx / "processed.h5ad")
    X = adata.X
    # processed is normalised/log-transformed (cpm_log1p for the registered runs); zeros stay zero under
    # log1p and CPM, so nonzero still means "mRNA seen in this cell". Scaling/centering would break this.
    det = np.asarray((X != 0).mean(axis=0)).ravel()
    mean = np.asarray(X.mean(axis=0)).ravel()
    genes = pd.Index(adata.var_names.astype(str).str.upper())
    keep = ~genes.duplicated()
    ref = pd.read_csv(Path(ddir) / "reference_evaluation.tsv", sep="\t").astype(str).apply(lambda c: c.str.upper())
    edges = {}
    for m in MODELS:
        f = ctx / f"{m}.tsv.gz"
        if f.exists():
            e = pd.read_csv(f, sep="\t")
            e[["source", "target"]] = e[["source", "target"]].astype(str).apply(lambda c: c.str.upper())
            edges[m] = e.sort_values("weight", ascending=False, kind="stable").drop_duplicates(["source", "target"])
    return {"genes": genes[keep], "det": pd.Series(det[keep], genes[keep]), "mean": pd.Series(mean[keep], genes[keep]),
            "ref": ref, "edges": edges, "by_src": {m: dict(tuple(e.groupby("source"))) for m, e in edges.items()}}


def load_all(artifacts: Path, arm: str = "auto") -> dict:
    return {d.name: load_dataset(d, arm) for d in sorted(Path(artifacts).iterdir())
            if (d / arm / "ctx0" / "processed.h5ad").exists() and (d / "reference_evaluation.tsv").exists()}


# ------------------------------------------------------------------------------------------ atlas
def _tf_row(name: str, d: dict, tf: str, min_targets: int = MIN_TARGETS) -> dict | None:
    genes = d["genes"]
    if tf not in genes:
        return None
    targets = set(d["ref"].loc[d["ref"].source == tf, "target"]) & set(genes) - {tf}
    if len(targets) < min_targets:
        return None
    gold = GoldIndex({(tf, t) for t in targets})
    universe = pd.DataFrame({"source": tf, "target": [g for g in genes if g != tf]})
    row = {"dataset": name, "tf": tf, "detection": float(d["det"][tf]), "mean_expr": float(d["mean"][tf]),
           "n_ref_targets": len(targets), "status": "EXPLORATORY"}
    for m in MODELS:
        e = d["by_src"].get(m, {}).get(tf)
        row[f"n_stored_{m}"] = 0 if e is None else len(e)
        row[f"auprc_ratio_{m}"] = (float(evaluate_ranking(e, gold, with_curve=False, universe=universe)["auprc_ratio"])
                                   if e is not None and targets else np.nan)
    return row


def build_atlas(data: dict) -> pd.DataFrame:
    rows = [r for name, d in data.items() for tf in sorted(set(d["ref"].source))
            if (r := _tf_row(name, d, tf)) is not None]
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------ figures
def _save(fig, stem: Path) -> list[Path]:
    fig.text(0.995, 0.995, "EXPLORATORY", ha="right", va="top", fontsize=6, color=VERMILLION)
    fig.savefig(stem.with_suffix(".pdf"), metadata={"CreationDate": None})
    fig.savefig(stem.with_suffix(".png"), metadata={"Software": None})
    plt.close(fig)
    return [stem.with_suffix(".pdf"), stem.with_suffix(".png")]


def _bin_ids(x: np.ndarray) -> tuple[np.ndarray, int]:
    """Equal-count bin index of each x."""
    edges = np.unique(np.quantile(x, np.linspace(0, 1, N_BINS + 1)))
    return np.clip(np.digitize(x, edges[1:-1]), 0, len(edges) - 2), len(edges) - 1


def _binned_median(x, y, idx, nb, rng):
    """Median of y per bin with a point bootstrap 95% CI (within bin)."""
    xs, med, ci = [], [], []
    for b in range(nb):
        v = y[idx == b]
        if len(v):
            boot = np.median(rng.choice(v, (N_BOOT, len(v))), axis=1)
            xs.append(np.median(x[idx == b])); med.append(np.median(v)); ci.append(np.percentile(boot, [2.5, 97.5]))
    return np.array(xs), np.array(med), np.array(ci).reshape(-1, 2)


def _cluster_mean(vals, ds, idx, nb, rng):
    """Per-bin mean of vals; 95% CI from resampling whole datasets (the unit), not points."""
    codes, uniq = pd.factorize(ds)
    S, C = np.zeros((len(uniq), nb)), np.zeros((len(uniq), nb))
    np.add.at(S, (codes, idx), np.asarray(vals, float)); np.add.at(C, (codes, idx), 1)
    draws = rng.integers(0, len(uniq), (N_BOOT, len(uniq)))
    with np.errstate(invalid="ignore", divide="ignore"):
        point = S.sum(0) / C.sum(0)
        boot = S[draws].sum(1) / C[draws].sum(1)
    return point, np.nanpercentile(boot, [2.5, 97.5], axis=0).T


def _family_view(atlas: pd.DataFrame, fam: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    a, b = FAMILIES[fam]
    both = atlas[(atlas[f"n_stored_{a}"] >= 1) & (atlas[f"n_stored_{b}"] >= 1)
                 & atlas[f"auprc_ratio_{a}"].notna() & atlas[f"auprc_ratio_{b}"].notna()].copy()
    both["lift"] = both[f"auprc_ratio_{b}"] - both[f"auprc_ratio_{a}"]
    return atlas, both


def detectability_stats(atlas: pd.DataFrame) -> dict:
    """Per family, on (dataset, TF) pairs ranked in both lists of the family (>= 1 stored edge each)."""
    rng = np.random.default_rng(SEED)
    out = {}
    for fam, (a, b) in FAMILIES.items():
        allp, both = _family_view(atlas, fam)
        ds_rho = [spearmanr(g.detection, g.lift)[0] for _, g in both.groupby("dataset") if len(g) >= 10]
        mean_lift = ci = [np.nan, np.nan]
        if len(both):
            m, c = _cluster_mean(both.lift, both.dataset, np.zeros(len(both), int), 1, rng)
            mean_lift, ci = float(m[0]), [float(v) for v in c[0]]
        out[fam] = {"n_rows": len(allp), "n_ranked_both": len(both), "excluded": len(allp) - len(both),
                    "rho": float(spearmanr(both.detection, both.lift)[0]) if len(both) > 2 else np.nan,
                    "rho_per_dataset_median": float(np.median(ds_rho)) if ds_rho else np.nan, "n_datasets_rho": len(ds_rho),
                    "mean_lift": mean_lift, "mean_lift_ci95_dataset_bootstrap": ci,
                    "frac_lift_positive": float((both.lift > 0).mean()) if len(both) else np.nan,
                    "median_without": float(both[f"auprc_ratio_{a}"].median()),
                    "median_with": float(both[f"auprc_ratio_{b}"].median()),
                    "frac_ratio_exactly_1_without": float((both[f"auprc_ratio_{a}"] == 1.0).mean()) if len(both) else np.nan,
                    "frac_ratio_exactly_1_with": float((both[f"auprc_ratio_{b}"] == 1.0).mean()) if len(both) else np.nan,
                    "coverage_without": float((allp[f"n_stored_{a}"] >= 1).mean()),
                    "coverage_with": float((allp[f"n_stored_{b}"] >= 1).mean())}
    return out


def fig_tf_detectability(atlas: pd.DataFrame, out: Path) -> tuple[list[Path], dict]:
    """Rows: ratio (ranked in both), coverage, mean lift, fraction of TFs with positive lift. Returns paths and the
    per-bin numbers behind the figure."""
    rng = np.random.default_rng(SEED)
    stats = detectability_stats(atlas)
    bins = {}
    with plt.rc_context(FIG_STYLE):
        fig, axes = plt.subplots(4, 3, figsize=(183 * MM, 170 * MM), sharex=True)
        for j, (fam, (a, b)) in enumerate(FAMILIES.items()):
            allp, both = _family_view(atlas, fam)
            col = OKABE_ITO[fam]
            xb = both.detection.to_numpy(); ib, nbb = _bin_ids(xb)
            xa = allp.detection.to_numpy(); ia, nba = _bin_ids(xa)
            clip = lambda v: np.clip(v, DET_FLOOR, None)  # noqa: E731 (log axis only)
            ya, yb = both[f"auprc_ratio_{a}"].to_numpy(), both[f"auprc_ratio_{b}"].to_numpy()
            r0, r1, r2, r3 = axes[:, j]
            for y, c, lab in ((ya, GREY, "without prior"), (yb, col, "with motif prior")):
                r0.scatter(clip(xb), y, s=2, color=c, alpha=0.12, linewidths=0, rasterized=True)
                bx, med, ci = _binned_median(xb, y, ib, nbb, rng)
                r0.errorbar(clip(bx), med, yerr=[med - ci[:, 0], ci[:, 1] - med], color=c, marker="o", ms=3, capsize=2, label=lab)
            xm = np.array([np.median(xa[ia == k]) for k in range(nba)])
            cov = {}
            for key, c, lab in ((a, GREY, "without prior"), (b, col, "with motif prior")):
                p, ci = _cluster_mean((allp[f"n_stored_{key}"] >= 1).to_numpy(), allp.dataset, ia, nba, rng)
                cov[key] = (p, ci)
                r1.errorbar(clip(xm), p, yerr=[p - ci[:, 0], ci[:, 1] - p], color=c, marker="o", ms=3, capsize=2, label=lab)
            xl = np.array([np.median(xb[ib == k]) for k in range(nbb)])
            lp, lci = _cluster_mean(both.lift, both.dataset, ib, nbb, rng)
            pp, pci = _cluster_mean((both.lift > 0).to_numpy(), both.dataset, ib, nbb, rng)
            r2.errorbar(clip(xl), lp, yerr=[lp - lci[:, 0], lci[:, 1] - lp], color="black", marker="o", ms=3, capsize=2)
            r3.errorbar(clip(xl), pp, yerr=[pp - pci[:, 0], pci[:, 1] - pp], color="black", marker="o", ms=3, capsize=2)
            r2.axhline(0, color=GREY, lw=0.6); r3.axhline(0.5, color=GREY, lw=0.6)
            bins[fam] = {"coverage": [{"detection_median": float(xm[k]), "without": float(cov[a][0][k]), "with": float(cov[b][0][k]),
                                       "n": int((ia == k).sum())} for k in range(nba)],
                         "ranked_both": [{"detection_median": float(xl[k]), "mean_lift": float(lp[k]),
                                          "ci": [float(v) for v in lci[k]], "frac_positive": float(pp[k]), "n": int((ib == k).sum())}
                                         for k in range(nbb)]}
            s = stats[fam]
            r0.set_title(f"{FAMILY_LABEL[fam]} (n = {s['n_ranked_both']}/{s['n_rows']})")
            r0.set_xscale("log"); r0.set_yscale("log")
            r3.set_xlabel("TF detection rate (fraction of cells > 0)")
            if j == 0:
                r0.set_ylabel("per-TF AUPRC ratio"); r1.set_ylabel("TFs with >= 1 stored edge")
                r2.set_ylabel("mean lift (with - without)"); r3.set_ylabel("fraction of TFs, lift > 0")
                r0.legend(frameon=False, loc="upper left"); r1.legend(frameon=False, loc="lower right")
            r1.set_ylim(0, 1.02); r3.set_ylim(0, 1)
        fig.suptitle("TF detectability and the motif prior", fontsize=8, x=0.01, ha="left", y=0.995)
        fig.tight_layout(rect=(0, 0, 1, 0.975))
        return _save(fig, out / "fig_tf_detectability"), bins


def _pct(df: pd.DataFrame) -> pd.Series:
    """Rank percentile of each stored edge (rank / list length; 1 = last stored edge)."""
    return df.weight.rank(method="min", ascending=False) / len(df)


def promotion_table(data: dict, k: int = PROMO_K) -> pd.DataFrame:
    """Per dataset and family: share of CollecTRI evaluation edges among kept / promoted / demoted top-k edges."""
    rows = []
    for name, d in data.items():
        ref = set(zip(d["ref"].source, d["ref"].target))
        for fam, (a, b) in FAMILIES.items():
            if a not in d["edges"] or b not in d["edges"]:
                continue
            ea, eb = d["edges"][a], d["edges"][b]
            ta, tb = set(zip(ea.source.iloc[:k], ea.target.iloc[:k])), set(zip(eb.source.iloc[:k], eb.target.iloc[:k]))
            base = np.mean([(s, t) in ref for s, t in zip(ea.source, ea.target)])
            for cat, s in (("kept", ta & tb), ("promoted", tb - ta), ("demoted", ta - tb)):
                n_ref = sum(e in ref for e in s)
                rows.append({"dataset": name, "family": fam, "category": cat, "n_edges": len(s), "n_ref": n_ref,
                             "frac_ref": n_ref / len(s) if s else np.nan, "frac_ref_all_stored_without": float(base)})
    return pd.DataFrame(rows)


def fig_promotion(data: dict, out: Path) -> tuple[list[Path], pd.DataFrame]:
    tab = promotion_table(data)
    rng = np.random.default_rng(SEED)
    cats = ["kept", "promoted", "demoted"]
    with plt.rc_context(FIG_STYLE):
        fig, axes = plt.subplots(1, 3, figsize=(183 * MM, 70 * MM), sharey=True)
        for ax, (fam, col) in zip(axes, OKABE_ITO.items()):
            t = tab[tab.family == fam]
            for i, c in enumerate(cats):
                v = t.loc[t.category == c, "frac_ref"].dropna().to_numpy()
                if not len(v):
                    continue
                ax.scatter(i + rng.uniform(-0.18, 0.18, len(v)), v, s=6, color=col, alpha=0.5, linewidths=0)
                boot = np.median(rng.choice(v, (N_BOOT, len(v))), axis=1)  # resamples datasets: one point each
                lo, hi = np.percentile(boot, [2.5, 97.5])
                ax.errorbar([i], [np.median(v)], yerr=[[np.median(v) - lo], [hi - np.median(v)]], color="black", marker="_", ms=10, capsize=3)
            ax.axhline(t.frac_ref_all_stored_without.median(), color=GREY, lw=0.6, ls="--")
            ax.set_xticks(range(3), cats); ax.set_title(FAMILY_LABEL[fam])
        axes[0].set_ylabel(f"share in CollecTRI eval. half (top-{PROMO_K})")
        fig.suptitle("Which edges does the prior promote?", fontsize=8, x=0.01, ha="left", y=0.995)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        return _save(fig, out / "fig_promoted_edges"), tab


# ------------------------------------------------------------------------------------------ card
def _cut(d: dict, model: str) -> float:
    """Rank percentile of the last edge inside the top 1% of the stored list (at least one edge)."""
    n = len(d["edges"][model])
    return np.ceil(TOP_FRAC * n) / n


def default_tf(atlas: pd.DataFrame) -> str:
    n = atlas.groupby("tf").dataset.nunique()
    return sorted(n[n == n.max()].index)[0]


def evidence_card(data: dict, tf: str, out: Path | None = None, model: str = "pearson_motif") -> str:
    """Markdown card for ``tf``; with ``out`` also writes card_<TF>.md and the neighbourhood png."""
    tf = tf.upper()
    rows = [r for n, d in data.items() if (r := _tf_row(n, d, tf, min_targets=0)) is not None]
    if not rows:
        return f"# Evidence card: {tf}\n\n{LABEL}\n\nTF not found in any dataset.\n"
    names = [r["dataset"] for r in rows]
    base = next(a for a, b in FAMILIES.values() if b == model)

    def pcts(n, m):  # target -> rank percentile of the TF's edges in the stored list of model m
        e = data[n]["edges"].get(m)
        return pd.Series(dtype=float) if e is None else _pct(e)[e.source == tf].set_axis(e.target[e.source == tf])

    pa, pb = {n: pcts(n, base) for n in names}, {n: pcts(n, model) for n in names}
    ref_t = {t for d in data.values() for s, t in zip(d["ref"].source, d["ref"].target) if s == tf}

    def stat(t):
        # a target counts in a dataset only if both genes are in its universe (TFs are everywhere: no pooled-mean bias)
        pres = [n for n in names if t in data[n]["genes"] and t != tf]
        fill = lambda p, n: p[n].get(t, 1.0)  # noqa: E731 (stored list lacks the edge: bottom)
        hb = [n for n in pres if fill(pb, n) <= _cut(data[n], model)]
        ha = [n for n in pres if fill(pa, n) <= _cut(data[n], base)]
        return {"target": t, "present": len(pres), "hits": len(hb), "hits_without": len(ha),
                "hit_rate": len(hb) / len(pres) if pres else np.nan, "hits_names": [_short(n) for n in hb],
                "med_with": float(np.median([fill(pb, n) for n in pres])) if pres else np.nan,
                "med_without": float(np.median([fill(pa, n) for n in pres])) if pres else np.nan, "in_ref": t in ref_t}

    cands = [stat(t) for t in sorted({t for s in pb.values() for t in s.index})]
    top = sorted([c for c in cands if c["present"] >= MIN_PRESENT], key=lambda c: (-c["hit_rate"], c["med_with"], c["target"]))[:10]
    known = sorted([c for c in map(stat, sorted(ref_t)) if c["present"] >= MIN_PRESENT], key=lambda c: (c["med_with"], c["target"]))
    md = [f"# Evidence card: {tf}", "", f"**{LABEL}.** Auto arm; {len(names)} datasets. Ratios are per-TF AUPRC ratios "
          f"(TF's predicted targets vs its CollecTRI evaluation targets, universe = dataset genes); NaN = TF absent from "
          f"that model's stored list or no reference target. n stored = the TF's edges in the stored (truncated) list.", "",
          "## Per dataset", "", "| dataset | detection | n ref targets | " + " | ".join(MODELS) + " | n stored (" + base + " / " + model + ") |",
          "|" + "---|" * (4 + len(MODELS))]
    for r in rows:
        md.append(f"| {_short(r['dataset'])} | {r['detection']:.3f} | {r['n_ref_targets']} | "
                  + " | ".join(f"{r[f'auprc_ratio_{m}']:.2f}" for m in MODELS) + f" | {r[f'n_stored_{base}']} / {r[f'n_stored_{model}']} |")
    md += ["", f"## Top 10 targets ({model})", "",
           f"Per target: present = datasets where both genes are in the universe (>= {MIN_PRESENT} required); hits = datasets with "
           f"the edge in the first {TOP_FRAC:.0%} of the stored list; hit rate = hits / present; ties broken by the median rank "
           "percentile over present datasets (edge missing from the stored list counts 1.0). In reference = in the CollecTRI "
           "evaluation half of any dataset.", "",
           "| target | present | hits | hit rate | datasets in top 1% | median pct without prior | median pct with prior | in reference |",
           "|---|---|---|---|---|---|---|---|"]
    for c in top:
        md.append(f"| {c['target']} | {c['present']} | {c['hits']} | {c['hit_rate']:.2f} | {', '.join(c['hits_names']) or '-'} | "
                  f"{c['med_without']:.4f} | {c['med_with']:.4f} | {'y' if c['in_ref'] else 'n'} |")
    md += ["", "## Known targets (CollecTRI evaluation half)", "",
           f"{len(known)} reference targets present in >= {MIN_PRESENT} datasets, best median percentile with prior first. "
           "Hits = datasets with the edge in the top 1% (without / with prior).", "",
           "| target | present | hits without | hits with | median pct without prior | median pct with prior |",
           "|---|---|---|---|---|---|"]
    md += [f"| {c['target']} | {c['present']} | {c['hits_without']} | {c['hits']} | {c['med_without']:.4f} | {c['med_with']:.4f} |"
           for c in known]
    text = "\n".join(md) + "\n"
    if out is not None:
        out = Path(out); out.mkdir(parents=True, exist_ok=True)
        (out / f"card_{tf}.md").write_text(text, encoding="utf-8", newline="\n")
        _card_figure(tf, top, out / f"card_{tf}")
    return text


def _card_figure(tf: str, top: list[dict], stem: Path) -> None:
    ang = np.linspace(0, 2 * np.pi, len(top), endpoint=False)
    all0 = all(c["hit_rate"] == 0 for c in top)
    with plt.rc_context(FIG_STYLE):
        fig, ax = plt.subplots(figsize=(89 * MM, 89 * MM))
        for c, a in zip(top, ang):
            x, y = np.cos(a), np.sin(a)
            ax.plot([0, x], [0, y], color=VERMILLION if c["in_ref"] else GREY, lw=0.5 + 5 * c["hit_rate"], zorder=1)
            ax.scatter([x], [y], s=60, color="white", edgecolor="black", linewidths=0.6, zorder=2)
            ax.text(x * 1.25, y * 1.25, c["target"], ha="center", va="center")
        ax.scatter([0], [0], s=700, color=OKABE_ITO["pearson"], zorder=2)
        ax.text(0, 0, tf, ha="center", va="center", color="white", zorder=3)
        ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.5, 1.5); ax.axis("off")
        ax.set_title(f"{tf}: top-{len(top)} targets by top-1% hit rate"
                     + ("\nall hit rates are 0 (no edge in any top 1%)" if all0 else "") + "\nwidth = hit rate, vermilion = in reference",
                     fontsize=6)
        fig.tight_layout()
        fig.text(0.995, 0.995, "EXPLORATORY", ha="right", va="top", fontsize=6, color=VERMILLION)
        fig.savefig(stem.with_suffix(".png"), metadata={"Software": None})
        plt.close(fig)


# ------------------------------------------------------------------------------------------ driver
def write_readme(atlas: pd.DataFrame, stats: dict, bins: dict, promo: pd.DataFrame, out: Path, arm: str, cards: list[str]) -> None:
    f3 = lambda v: f"{v:+.3f}"  # noqa: E731
    lines = [f"# grnbench explore: {LABEL}", "", f"Arm: {arm}. {atlas.dataset.nunique()} datasets, {len(atlas)} (dataset, TF) rows "
             f"with >= {MIN_TARGETS} reference targets in the dataset. Everything here is descriptive: no test, no claim.", "",
             "## fig_tf_detectability", "", f"Question: {QUESTION}", "",
             "Stored edge lists are truncated (top 20,000), so a TF with no stored edge is unranked, and a TF whose stored edges "
             "hit no reference target scores exactly 1.0. Ratios and lifts use only (dataset, TF) pairs ranked in both lists of the "
             "family (>= 1 stored edge each); excluded pairs are counted. Rows: per-TF AUPRC ratio (medians of 5 equal-count "
             "detection bins, point bootstrap 95% CI); coverage = fraction of all TFs with >= 1 stored edge (the prior may let "
             "low-detection TFs into the ranking); mean lift and fraction of TFs with positive lift (CIs resample datasets, not points; "
             f"seed {SEED}). Lift = ratio with motif prior minus without.", "",
             "| family | ranked in both / rows | excluded | Spearman rho (pairs) | median of per-dataset rho (n datasets) | mean lift [95% CI] | "
             "frac lift > 0 | median ratio without / with | share of ratios exactly 1.0 without / with | coverage without / with |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for f, s in stats.items():
        lines.append(f"| {FAMILY_LABEL[f]} | {s['n_ranked_both']} / {s['n_rows']} | {s['excluded']} | {f3(s['rho'])} | "
                     f"{f3(s['rho_per_dataset_median'])} ({s['n_datasets_rho']}) | {f3(s['mean_lift'])} "
                     f"[{s['mean_lift_ci95_dataset_bootstrap'][0]:+.3f}, {s['mean_lift_ci95_dataset_bootstrap'][1]:+.3f}] | "
                     f"{s['frac_lift_positive']:.3f} | {s['median_without']:.3f} / {s['median_with']:.3f} | "
                     f"{s['frac_ratio_exactly_1_without']:.3f} / {s['frac_ratio_exactly_1_with']:.3f} | "
                     f"{s['coverage_without']:.3f} / {s['coverage_with']:.3f} |")
    lines += ["", "Coverage per detection bin (without -> with prior):", ""]
    for f, b in bins.items():
        lines.append(f"- {FAMILY_LABEL[f]}: " + "; ".join(f"det~{c['detection_median']:.3f} (n={c['n']}): {c['without']:.2f} -> {c['with']:.2f}"
                                                         for c in b["coverage"]))
    lines += ["", f"## fig_promoted_edges", "", f"Which edges does the prior promote? Per dataset and family, top-{PROMO_K} of the stored lists: "
              "kept = in both, promoted = with-prior only, demoted = without-prior only. One point per dataset = share that are "
              "CollecTRI evaluation edges; black = median with 95% CI (resampling datasets); dashed = share among all stored "
              "without-prior edges.", "",
              "| family | category | median share | datasets |", "|---|---|---|---|"]
    for (f, c), g in promo.groupby(["family", "category"], sort=False):
        lines.append(f"| {FAMILY_LABEL[f]} | {c} | {g.frac_ref.median():.4f} | {g.frac_ref.notna().sum()} |")
    lines += ["", "## Cards", "", "Top-10 targets ordered by top-1% hit rate across datasets where both genes are present (>= "
              f"{MIN_PRESENT}); second table lists the TF's reference targets, i.e. what the networks miss.", "",
              "Files: tf_atlas.tsv, promotion.tsv, stats.json, fig_tf_detectability.{pdf,png}, fig_promoted_edges.{pdf,png}, "
              + ", ".join(f"card_{t}.{{md,png}}" for t in cards) + "."]
    (out / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def run(artifacts: Path, out: Path, arm: str = "auto", tf: str | None = None) -> dict:
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    data = load_all(artifacts, arm)
    atlas = build_atlas(data)
    atlas.to_csv(out / "tf_atlas.tsv", sep="\t", index=False, float_format="%.6g", lineterminator="\n")
    stats = detectability_stats(atlas)
    _, bins = fig_tf_detectability(atlas, out)
    _, promo = fig_promotion(data, out)
    promo.to_csv(out / "promotion.tsv", sep="\t", index=False, float_format="%.6g", lineterminator="\n")
    chosen = [tf.upper()] if tf else [default_tf(atlas)] + (["TP53"] if "TP53" in set(atlas.tf) else [])
    chosen = list(dict.fromkeys(chosen))
    for t in chosen:
        evidence_card(data, t, out)
    write_readme(atlas, stats, bins, promo, out, arm, chosen)
    summary = {"rows": len(atlas), "stats": stats, "bins": bins, "cards": chosen,
               "promotion_median": {f"{f}/{c}": float(g.frac_ref.median()) for (f, c), g in promo.groupby(["family", "category"])}}
    (out / "stats.json").write_text(json.dumps(summary, indent=1), encoding="utf-8", newline="\n")
    return summary


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="grnbench.explore", description=__doc__.split("\n")[0])
    p.add_argument("--artifacts", required=True, help="folder with one sub-folder per dataset")
    p.add_argument("--out", required=True)
    p.add_argument("--arm", default="auto")
    p.add_argument("--tf", default=None, help="card TF (default: the TF in most datasets, plus TP53 if present)")
    a = p.parse_args(argv)
    res = run(Path(a.artifacts), Path(a.out), a.arm, a.tf)
    print(json.dumps(res["stats"], indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
