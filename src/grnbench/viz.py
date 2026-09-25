"""Visualization / explainability scaffold for grnbench runs.

Reads what ``benchmarks/real_bench.py`` wrote (``rows/*.json`` always, plus
``artifacts/<task>/...`` when it ran with ``--keep-artifacts``) and turns it
into figures a person can look at to answer a specific question, never a
generic dashboard. Matplotlib is used with the ``Agg`` backend only (no
display needed) and the shared style in ``viz.mplstyle``.

Every figure function returns ``(path, None)`` on success or ``(None, reason)``
when it could not be made (missing inputs, too few points, an optional
dependency not installed). :func:`build_run` walks a whole run directory,
calls what it can, and writes ``index.md`` listing every figure made and
every one skipped, with the reason.

Learned-model attribution (:func:`gradient_x_input`, :func:`integrated_gradients`)
needs ``torch`` and a ``torch.nn.Module``; grnbench's own inference methods are
not torch modules today, so the CLI never calls these - they are here for the
next learned model to use directly.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

_STYLE = Path(__file__).with_name("viz.mplstyle")
DPI = 120
SEED = 0
MAX_EMBED_CELLS = 3000

_BATCH_KEYS = ("batch", "donor_id", "sample", "assay")
_TYPE_KEYS = ("cell_type", "cell_type_ontology_term_id", "celltype", "louvain", "leiden")


def _use_style():
    if str(_STYLE) not in plt.style.available and _STYLE.exists():
        plt.style.use(str(_STYLE))
    else:
        plt.style.use(str(_STYLE))


def _n_label(n: int, unit: str) -> str:
    return f"n = {n:,} {unit}"


def _save(fig, out_png, dpi: int = DPI) -> Path:
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    _use_style()
    fig.savefig(out_png, dpi=dpi)
    plt.close(fig)
    return out_png


def _missing(reason: str):
    return None, reason


def _read_edges(path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="gzip" if str(path).endswith(".gz") else None)


# --------------------------------------------------------------------------- #
# QC / embeddings on raw vs. processed data
# --------------------------------------------------------------------------- #

def _counts(adata):
    """(counts matrix as CSR, gene symbols). Uses ``raw.X`` when ``X`` is normalized (CELLxGENE
    keeps counts there) and ``var['feature_name']`` when var_names are Ensembl ids."""
    import scipy.sparse as sp

    X, var = adata.X, adata.var
    if adata.raw is not None and adata.raw.shape == adata.shape:
        X, var = adata.raw.X, adata.raw.var
    X = sp.csr_matrix(X)
    names = var["feature_name"].astype(str) if "feature_name" in var else pd.Index(adata.var_names).astype(str)
    return X, pd.Index(names)


def _qc_frame(adata) -> pd.DataFrame:
    X, genes = _counts(adata)
    total = np.asarray(X.sum(axis=1)).ravel()
    n_genes = np.diff(X.indptr)
    mito = np.asarray(genes.str.upper().str.startswith("MT-"), dtype=bool)
    mito_sum = np.asarray(X[:, mito].sum(axis=1)).ravel() if mito.any() else np.zeros(X.shape[0])
    with np.errstate(divide="ignore", invalid="ignore"):
        pct = np.where(total > 0, mito_sum / total * 100.0, 0.0)
    return pd.DataFrame({"total_counts": total, "n_genes": n_genes, "pct_mito": pct})


def qc_before_after(raw_adata, processed_adata, out_png):
    """Per-cell total counts, genes detected and % mito, before vs. after cleaning."""
    if raw_adata is None:
        return _missing("no raw AnnData given")
    before = _qc_frame(raw_adata)
    after = _qc_frame(processed_adata) if processed_adata is not None else None
    metrics = [("total_counts", "total counts / cell"), ("n_genes", "genes detected / cell"),
               ("pct_mito", "% mitochondrial counts")]
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6))
    for ax, (col, label) in zip(axes, metrics):
        data = [before[col]] + ([after[col]] if after is not None else [])
        labels = ["before"] + (["after"] if after is not None else [])
        ax.boxplot(data, tick_labels=labels, showfliers=False)
        ax.set_ylabel(label)
    n_after = f", {_n_label(processed_adata.n_obs, 'cells after')}" if processed_adata is not None else \
        ", after: not kept by this run"
    fig.suptitle(f"Does cleaning change per-cell QC? ({_n_label(raw_adata.n_obs, 'cells before')}{n_after})")
    fig.tight_layout()
    return _save(fig, out_png), None


def _first_present(cols, candidates):
    for c in candidates:
        if c in cols:
            return c
    return None


def _pca_2d(adata, seed=SEED) -> np.ndarray:
    """PCA of log1p CPM on the 2,000 most variable genes, from counts."""
    from sklearn.decomposition import PCA

    X, _ = _counts(adata)
    total = np.asarray(X.sum(axis=1)).ravel()
    total[total == 0] = 1
    import scipy.sparse as sp
    L = sp.diags(1e6 / total) @ X
    L = L.log1p() if hasattr(L, "log1p") else np.log1p(L)
    mean = np.asarray(L.mean(axis=0)).ravel()
    var = np.asarray(L.multiply(L).mean(axis=0)).ravel() - mean ** 2
    top = np.argsort(-np.nan_to_num(var))[: min(2000, L.shape[1])]
    D = np.asarray(L[:, top].todense()) - mean[top]
    n = max(2, min(20, D.shape[0] - 1, D.shape[1] - 1))
    return PCA(n_components=n, random_state=seed).fit_transform(D)


def embeddings(adata_raw, adata_processed, color_keys, out_png):
    """PCA (+ UMAP when umap-learn is installed) before vs. after cleaning, coloured by cell type
    and by batch (one row each, when the column exists)."""
    if adata_raw is None:
        return _missing("no raw AnnData given")
    rng = np.random.default_rng(SEED)
    panels = [("before", adata_raw)] + ([("after", adata_processed)] if adata_processed is not None else [])
    keys = [k for k in (_first_present(adata_raw.obs.columns, _TYPE_KEYS),
                        _first_present(adata_raw.obs.columns, _BATCH_KEYS)) if k]
    keys = keys or [None]
    try:
        import umap
    except ImportError:
        umap = None
    rows = [("PCA", k) for k in keys] + ([("UMAP", keys[0])] if umap is not None else [])
    fig, axes = plt.subplots(len(rows), len(panels), figsize=(4.4 * len(panels), 3.8 * len(rows)),
                             squeeze=False)
    total_cells = 0
    for j, (label, ad_) in enumerate(panels):
        sub = ad_
        if sub.n_obs > MAX_EMBED_CELLS:
            sub = sub[np.sort(rng.choice(sub.n_obs, MAX_EMBED_CELLS, replace=False))].copy()
        total_cells += sub.n_obs
        pcs = _pca_2d(sub)
        emb = umap.UMAP(random_state=SEED).fit_transform(pcs) if umap is not None else None
        for i, (kind, key) in enumerate(rows):
            colors = sub.obs[key].astype(str).to_numpy() if key and key in sub.obs else None
            xy = pcs[:, :2] if kind == "PCA" else emb
            _scatter(axes[i, j], xy, colors, f"{kind} {label}" + (f", colour = {key}" if key else ""))
    fig.suptitle(f"Do cell types and batches separate before vs. after cleaning? "
                 f"({_n_label(total_cells, 'cells shown')})")
    fig.tight_layout()
    return _save(fig, out_png, dpi=80), None   # scatter-heavy: lower dpi keeps results small


def _scatter(ax, xy, colors, title):
    if colors is None:
        ax.scatter(xy[:, 0], xy[:, 1], s=6, alpha=0.7)
    else:
        # 8 colour-blind-safe colours: the 7 largest groups get one each, the rest are grey
        counts = pd.Series(colors).value_counts()
        top = list(counts.index[:7]) if len(counts) > 8 else list(counts.index)
        other = ~pd.Series(colors).isin(top).to_numpy()
        if other.any():
            ax.scatter(xy[other, 0], xy[other, 1], s=5, alpha=0.5, color="#BBBBBB",
                       label=f"other ({len(counts) - len(top)} groups)")
        for c in top:
            m = colors == c
            ax.scatter(xy[m, 0], xy[m, 1], s=5, alpha=0.7, label=str(c)[:24])
        ax.legend(markerscale=2, fontsize=7, loc="best", frameon=False)
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])


# --------------------------------------------------------------------------- #
# Decisions
# --------------------------------------------------------------------------- #

def decisions_markdown(row_or_decisions) -> str:
    """Markdown table of the pipeline decisions the engine took, with reasons."""
    d = row_or_decisions
    if isinstance(d, (str, Path)):
        d = json.loads(Path(d).read_text(encoding="utf-8"))
    lines = ["| decision | value | reason |", "|---|---|---|"]
    choices = {}
    if isinstance(d, dict) and "rules_choices" in d:
        rc = d["rules_choices"]
        choices = json.loads(rc) if isinstance(rc, str) else (rc or {})
    for key, val in choices.items():
        lines.append(f"| {key} | {val} | rules recommendation |")
    if isinstance(d, dict) and "chosen_plans" in d:
        for arm, plan in d["chosen_plans"].items():
            lines.append(f"| candidate ({arm}) | {plan.get('candidate')} | chosen for the {arm} arm |")
            for w in plan.get("warnings") or []:
                lines.append(f"| warning ({arm}) | - | {w} |")
    elif isinstance(d, dict) and "candidates" in d:
        for arm, cand in d["candidates"].items():
            lines.append(f"| candidate ({arm}) | {cand} | chosen for the {arm} arm |")
    if len(lines) == 2:
        lines.append("| (none found) | - | - |")
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- #
# Scores across datasets
# --------------------------------------------------------------------------- #

def _model_cols(rows_df: pd.DataFrame, arm: str, metric: str) -> dict[str, str]:
    suffix = f"_{arm}_{metric}"
    return {c[2:-len(suffix)]: c for c in rows_df.columns
            if c.startswith("m_") and c.endswith(suffix)}


def score_distributions(rows_df: pd.DataFrame, metric: str, arm: str, out_png):
    """Per-model distribution of a metric across datasets (strip + box), 1 = random."""
    cols = _model_cols(rows_df, arm, metric)
    if not cols:
        return _missing(f"no m_*_{arm}_{metric} columns in rows_df")
    data = {m: rows_df[c].dropna() for m, c in cols.items()}
    data = {m: v for m, v in data.items() if len(v)}
    if not data:
        return _missing("no finite values for any model")
    order = sorted(data, key=lambda m: -data[m].mean())
    fig, ax = plt.subplots(figsize=(max(5, 1.1 * len(order)), 4.2))
    ax.boxplot([data[m] for m in order], tick_labels=order, showfliers=False)
    rng = np.random.default_rng(SEED)
    for i, m in enumerate(order, start=1):
        x = rng.normal(i, 0.05, size=len(data[m]))
        ax.scatter(x, data[m], s=14, alpha=0.6, zorder=3)
    ax.axhline(1.0, color="grey", linestyle="--", linewidth=1, label="1 = random")
    ax.set_ylabel(metric)
    ax.tick_params(axis="x", rotation=30)
    n_ds = max(len(v) for v in data.values())
    fig.suptitle(f"Which model wins on {metric} ({arm} pipeline)? ({_n_label(n_ds, 'datasets')})")
    ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    return _save(fig, out_png), None


def paired_lift(rows_df: pd.DataFrame, a: str, b: str, metric: str, arm: str, out_png):
    """Per-dataset paired lines a vs. b for one metric; mean difference in the title."""
    cols_a, cols_b = _model_cols(rows_df, arm, metric), _model_cols(rows_df, arm, metric)
    ca, cb = cols_a.get(a), cols_b.get(b)
    if ca is None or cb is None:
        return _missing(f"model column missing for {a!r} or {b!r} ({arm}/{metric})")
    pair = rows_df[[ca, cb]].dropna()
    if len(pair) < 2:
        return _missing(f"fewer than 2 datasets with both {a} and {b} scored")
    diff = (pair[ca] - pair[cb])
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    for _, r in pair.iterrows():
        ax.plot([0, 1], [r[ca], r[cb]], color="0.6", alpha=0.6, marker="o", markersize=4)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([a, b])
    ax.set_ylabel(metric)
    fig.suptitle(f"Does {a} beat {b} on {metric} ({arm})? mean diff {diff.mean():+.3f} "
                f"({_n_label(len(pair), 'datasets')})")
    fig.tight_layout()
    return _save(fig, out_png), None


# --------------------------------------------------------------------------- #
# Per-dataset edge-level figures
# --------------------------------------------------------------------------- #

def _universe_prevalence(edges: pd.DataFrame, ref_pairs) -> tuple[np.ndarray, float]:
    ref = {(str(s).upper(), str(t).upper()) for s, t in ref_pairs}
    labels = np.fromiter(((str(s).upper(), str(t).upper()) in ref
                          for s, t in zip(edges["source"], edges["target"])),
                         dtype=bool, count=len(edges))
    prevalence = labels.mean() if len(labels) else 0.0
    return labels, float(prevalence)


def pr_curves(edges_by_model: dict, ref_pairs, out_png):
    """PR curve per model on its own edge universe; dashed line = prevalence."""
    edges_by_model = {m: e for m, e in (edges_by_model or {}).items() if e is not None and len(e)}
    ref_pairs = list(ref_pairs or [])
    if not edges_by_model or not ref_pairs:
        return _missing("no edges or no reference pairs given")
    fig, ax = plt.subplots(figsize=(5, 4.6))
    n_edges = 0
    prevalences = []
    for m, edges in edges_by_model.items():
        ordered = edges.sort_values("weight", ascending=False, kind="stable")
        labels, prevalence = _universe_prevalence(ordered, ref_pairs)
        prevalences.append(prevalence)
        n_edges = max(n_edges, len(ordered))
        if labels.sum() == 0:
            continue
        tp = np.cumsum(labels)
        rank = np.arange(1, len(labels) + 1)
        precision = tp / rank
        recall = tp / labels.sum()
        ax.plot(recall, precision, label=m)
    if prevalences:
        ax.axhline(float(np.mean(prevalences)), color="grey", linestyle="--", linewidth=1,
                   label="prevalence (random)")
    ax.set_xlabel("recall")
    ax.set_ylabel("precision")
    ax.legend(frameon=False, fontsize=8)
    fig.suptitle(f"Which model ranks reference edges highest? ({_n_label(n_edges, 'edges scored')})")
    fig.tight_layout()
    return _save(fig, out_png), None


def topk_network(edges: pd.DataFrame, ref_pairs, k: int, out_png):
    """networkx drawing of the top-k predicted edges, TP vs. FP, TFs distinct."""
    try:
        import networkx as nx
    except ImportError:
        return _missing("networkx is not installed")
    if edges is None or not len(edges):
        return _missing("no edges given")
    ref = {(str(s).upper(), str(t).upper()) for s, t in (ref_pairs or [])}
    top = edges.sort_values("weight", ascending=False, kind="stable").head(k)
    if not len(top):
        return _missing("no edges after ranking")
    g = nx.DiGraph()
    tfs = set(top["source"].astype(str))
    for _, r in top.iterrows():
        g.add_edge(str(r["source"]), str(r["target"]),
                   tp=(str(r["source"]).upper(), str(r["target"]).upper()) in ref)
    pos = nx.spring_layout(g, seed=SEED)
    fig, ax = plt.subplots(figsize=(6.5, 6.5))
    node_colors = ["#E69F00" if n in tfs else "#0072B2" for n in g.nodes]
    nx.draw_networkx_nodes(g, pos, ax=ax, node_size=60, node_color=node_colors)
    tp_edges = [(u, v) for u, v, d in g.edges(data=True) if d["tp"]]
    fp_edges = [(u, v) for u, v, d in g.edges(data=True) if not d["tp"]]
    nx.draw_networkx_edges(g, pos, edgelist=fp_edges, ax=ax, edge_color="#D55E00", alpha=0.5, width=0.8)
    nx.draw_networkx_edges(g, pos, edgelist=tp_edges, ax=ax, edge_color="#009E73", alpha=0.9, width=1.4)
    n_labels = min(len(g.nodes), 30)
    labels = {n: n for n in list(g.nodes)[:n_labels]}
    nx.draw_networkx_labels(g, pos, labels=labels, ax=ax, font_size=6)
    ax.axis("off")
    n_tp = len(tp_edges)
    fig.suptitle(f"Which top-{len(top)} predicted edges are known regulation? "
                f"({_n_label(len(top), 'edges')}, {n_tp} match the reference, TF nodes in orange)")
    fig.tight_layout()
    return _save(fig, out_png), None


def tf_recovery(edges_by_model: dict, ref_pairs, k_per_tf: int, out_png, out_tsv):
    """Per-TF hits in its top-k targets per model: heatmap + TSV with false positives."""
    edges_by_model = {m: e for m, e in (edges_by_model or {}).items() if e is not None and len(e)}
    ref_pairs = list(ref_pairs or [])
    if not edges_by_model or not ref_pairs:
        return _missing("no edges or no reference pairs given")
    ref_df = pd.DataFrame(ref_pairs, columns=["source", "target"]).astype(str)
    ref_df["source"], ref_df["target"] = ref_df["source"].str.upper(), ref_df["target"].str.upper()
    ref_by_tf = ref_df.groupby("source")["target"].apply(set)
    top_tfs = ref_by_tf.apply(len).sort_values(ascending=False).head(30).index
    if not len(top_tfs):
        return _missing("no TF has reference targets")
    hit_rows, fp_rows = [], []
    heat = pd.DataFrame(index=top_tfs, columns=list(edges_by_model), dtype=float)
    for m, edges in edges_by_model.items():
        e = edges.copy()
        e["source"], e["target"] = e["source"].astype(str).str.upper(), e["target"].astype(str).str.upper()
        for tf in top_tfs:
            targets = ref_by_tf.get(tf, set())
            sub = e[e["source"] == tf].sort_values("weight", ascending=False, kind="stable").head(k_per_tf)
            hits = int(sub["target"].isin(targets).sum())
            heat.loc[tf, m] = hits / max(len(targets), 1)
            hit_rows.append({"tf": tf, "model": m, "hits": hits, "n_ref": len(targets)})
            fp = sub[~sub["target"].isin(targets)].head(3)
            for _, r in fp.iterrows():
                fp_rows.append({"tf": tf, "model": m, "false_positive_target": r["target"],
                                "weight": r["weight"]})
    fig, ax = plt.subplots(figsize=(1.1 * len(edges_by_model) + 2, 0.28 * len(top_tfs) + 1.5))
    im = ax.imshow(heat.to_numpy(dtype=float), aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(len(heat.columns)))
    ax.set_xticklabels(heat.columns, rotation=30, ha="right")
    ax.set_yticks(range(len(heat.index)))
    ax.set_yticklabels(heat.index, fontsize=6)
    fig.colorbar(im, ax=ax, label="fraction of reference targets recovered")
    fig.suptitle(f"Which TFs does each model recover in its top-{k_per_tf} targets? "
                f"({_n_label(len(top_tfs), 'TFs')})")
    fig.tight_layout()
    path = _save(fig, out_png)
    out_tsv = Path(out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(hit_rows).to_csv(out_tsv, sep="\t", index=False)
    if fp_rows:
        fp_path = out_tsv.with_name(out_tsv.stem + "_false_positives.tsv")
        pd.DataFrame(fp_rows).to_csv(fp_path, sep="\t", index=False)
    return path, None


def motif_lift_by_tf(base_edges: pd.DataFrame, motif_edges: pd.DataFrame, ref_pairs, out_png, out_tsv,
                     k_per_tf: int = 20):
    """Per-TF change in hits@k between a base model and its motif version."""
    if base_edges is None or motif_edges is None or not len(base_edges) or not len(motif_edges):
        return _missing("base or motif edges missing")
    ref_pairs = list(ref_pairs or [])
    if not ref_pairs:
        return _missing("no reference pairs given")
    ref_df = pd.DataFrame(ref_pairs, columns=["source", "target"]).astype(str)
    ref_df["source"], ref_df["target"] = ref_df["source"].str.upper(), ref_df["target"].str.upper()
    ref_by_tf = ref_df.groupby("source")["target"].apply(set)

    def hits_per_tf(edges):
        e = edges.copy()
        e["source"], e["target"] = e["source"].astype(str).str.upper(), e["target"].astype(str).str.upper()
        out = {}
        for tf, targets in ref_by_tf.items():
            sub = e[e["source"] == tf].sort_values("weight", ascending=False, kind="stable").head(k_per_tf)
            out[tf] = int(sub["target"].isin(targets).sum())
        return pd.Series(out)

    base_hits, motif_hits = hits_per_tf(base_edges), hits_per_tf(motif_edges)
    lift = (motif_hits - base_hits).dropna().sort_values()
    if not len(lift):
        return _missing("no shared TFs between base and motif edges")
    top = pd.concat([lift.head(10), lift.tail(10)]).drop_duplicates()
    top = top.sort_values()
    fig, ax = plt.subplots(figsize=(5.5, max(3, 0.3 * len(top))))
    colors = ["#D55E00" if v < 0 else "#009E73" for v in top]
    ax.barh(top.index, top.to_numpy(), color=colors)
    ax.axvline(0, color="grey", linewidth=1)
    ax.set_xlabel(f"hits@{k_per_tf} (motif - base)")
    fig.suptitle(f"Which TFs gain or lose from the motif prior? ({_n_label(len(lift), 'TFs compared')})")
    fig.tight_layout()
    path = _save(fig, out_png)
    out_tsv = Path(out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"tf": lift.index, "hits_base": base_hits.reindex(lift.index),
                 "hits_motif": motif_hits.reindex(lift.index), "lift": lift.to_numpy()}
                ).to_csv(out_tsv, sep="\t", index=False)
    return path, None


def calibration(edges: pd.DataFrame, ref_pairs, out_png, n_deciles: int = 10):
    """Fraction of reference edges per score-rank decile, with the prevalence line."""
    if edges is None or not len(edges):
        return _missing("no edges given")
    ref_pairs = list(ref_pairs or [])
    if not ref_pairs:
        return _missing("no reference pairs given")
    ordered = edges.sort_values("weight", ascending=False, kind="stable").reset_index(drop=True)
    labels, prevalence = _universe_prevalence(ordered, ref_pairs)
    decile = np.minimum((np.arange(len(labels)) * n_deciles) // max(len(labels), 1), n_deciles - 1)
    frac = pd.Series(labels).groupby(decile).mean()
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.bar(frac.index + 1, frac.to_numpy(), color="#0072B2")
    ax.axhline(prevalence, color="grey", linestyle="--", linewidth=1, label="prevalence (random)")
    ax.set_xlabel("score decile (1 = highest-scored edges)")
    ax.set_ylabel("fraction that are reference edges")
    ax.legend(frameon=False)
    fig.suptitle(f"Are top-scored edges really more likely correct? ({_n_label(len(ordered), 'edges')})")
    fig.tight_layout()
    return _save(fig, out_png), None


def stability(edge_lists: dict, k: int, out_png):
    """Pairwise Jaccard of top-k edges across runs/contexts/seeds."""
    edge_lists = {name: e for name, e in (edge_lists or {}).items() if e is not None and len(e)}
    if len(edge_lists) < 2:
        return _missing("fewer than 2 edge lists given")
    names = list(edge_lists)
    tops = {n: set(zip(edge_lists[n].sort_values("weight", ascending=False, kind="stable")
                       .head(k)["source"].astype(str).str.upper(),
                       edge_lists[n].sort_values("weight", ascending=False, kind="stable")
                       .head(k)["target"].astype(str).str.upper()))
           for n in names}
    mat = pd.DataFrame(index=names, columns=names, dtype=float)
    for a in names:
        for b in names:
            union = tops[a] | tops[b]
            mat.loc[a, b] = len(tops[a] & tops[b]) / len(union) if union else np.nan
    fig, ax = plt.subplots(figsize=(1.0 * len(names) + 2, 1.0 * len(names) + 1.5))
    im = ax.imshow(mat.to_numpy(dtype=float), vmin=0, vmax=1, cmap="viridis")
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=30, ha="right")
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names)
    fig.colorbar(im, ax=ax, label="Jaccard of top-k edges")
    fig.suptitle(f"How stable is the top-{k} edge list across runs? "
                f"({_n_label(len(names), 'edge lists compared')})")
    fig.tight_layout()
    return _save(fig, out_png), None


# --------------------------------------------------------------------------- #
# Learned-model attribution (opt-in, needs torch + a torch.nn.Module)
# --------------------------------------------------------------------------- #

def gradient_x_input(module, x):
    """(x * d out.sum()/dx) as numpy. Needs torch and a torch.nn.Module; not
    called by the CLI (no learned model artifact ships in the repo today)."""
    try:
        import torch
    except ImportError as exc:
        raise ImportError("gradient_x_input needs torch: pip install torch") from exc
    if not isinstance(module, torch.nn.Module):
        raise ImportError("gradient_x_input needs a torch.nn.Module")
    x = x.clone().detach().requires_grad_(True) if torch.is_tensor(x) else torch.tensor(
        np.asarray(x, dtype=np.float32), requires_grad=True)
    out = module(x)
    out.sum().backward()
    return (x * x.grad).detach().cpu().numpy()


def integrated_gradients(module, x, baseline=None, steps: int = 32):
    """Integrated gradients (Sundararajan et al. 2017). Same requirements as
    :func:`gradient_x_input`."""
    try:
        import torch
    except ImportError as exc:
        raise ImportError("integrated_gradients needs torch: pip install torch") from exc
    if not isinstance(module, torch.nn.Module):
        raise ImportError("integrated_gradients needs a torch.nn.Module")
    x_t = x if torch.is_tensor(x) else torch.tensor(np.asarray(x, dtype=np.float32))
    base = torch.zeros_like(x_t) if baseline is None else (
        baseline if torch.is_tensor(baseline) else torch.tensor(np.asarray(baseline, dtype=np.float32)))
    total = torch.zeros_like(x_t)
    for i in range(1, steps + 1):
        alpha = i / steps
        pt = (base + alpha * (x_t - base)).clone().detach().requires_grad_(True)
        out = module(pt)
        out.sum().backward()
        total = total + pt.grad
    avg_grad = total / steps
    return ((x_t - base) * avg_grad).detach().cpu().numpy()


# --------------------------------------------------------------------------- #
# build_run: walk a real_bench.py output directory, write every figure it can
# --------------------------------------------------------------------------- #

def _load_rows(run_dir: Path) -> pd.DataFrame:
    rows_dir = Path(run_dir) / "rows"
    if not rows_dir.is_dir():
        return pd.DataFrame()
    rows = [json.loads(p.read_text(encoding="utf-8")) for p in sorted(rows_dir.glob("*.json"))]
    return pd.DataFrame(rows)


def _task_artifacts(run_dir: Path, task: str) -> Path | None:
    d = Path(run_dir) / "artifacts" / task
    return d if d.is_dir() else None


def _load_edges_by_model(ctx_dir: Path) -> dict[str, pd.DataFrame]:
    out = {}
    for f in sorted(ctx_dir.glob("*.tsv.gz")):
        try:
            out[f.name[: -len(".tsv.gz")]] = _read_edges(f)
        except Exception:  # noqa: BLE001 - a corrupt artifact is skipped, not fatal
            continue
    return out


_BASE_ALIAS = {"pearson": "pearson_signed", "spearman": "spearman_signed"}


def _motif_pairs(model_names) -> list[tuple[str, str]]:
    """(base, combined) pairs: 'pearson_motif' -> ('pearson_signed', 'pearson_motif')."""
    names = set(model_names)
    pairs = []
    for m in sorted(names):
        for suffix in ("_motif", "_cistarget"):
            if m.endswith(suffix) and m != "motif_prior":
                base = m[: -len(suffix)]
                base = _BASE_ALIAS.get(base, base)
                if base in names:
                    pairs.append((base, m))
    return pairs


def _raw_qc(task, prefix, tasks_map, ctx_dirs, out_dir, note, skipped) -> None:
    """QC and embeddings from the raw data (tasks.json path), with the processed data when kept."""
    raw_path = tasks_map.get(task)
    proc_h5 = ctx_dirs[0] / "processed.h5ad" if ctx_dirs else None
    if raw_path and Path(raw_path).exists():
        import anndata as ad
        raw = ad.read_h5ad(raw_path) if str(raw_path).endswith(".h5ad") else None
        processed = ad.read_h5ad(proc_h5) if proc_h5 and proc_h5.exists() else None
        if raw is not None:
            p, r = qc_before_after(raw, processed, out_dir / f"qc_{task}.png")
            note(prefix, p, r, "tasks.json (raw data)", n=f"n = {raw.n_obs} cells")
            p, r = embeddings(raw, processed, list(_TYPE_KEYS + _BATCH_KEYS),
                              out_dir / f"embed_{task}.png")
            note(prefix, p, r, "tasks.json (raw data)", n=f"n = {raw.n_obs} cells")
        else:
            skipped.append(f"{prefix}: raw data path is not an .h5ad, QC/embeddings skipped")
    else:
        skipped.append(f"{prefix}: no QC/embeddings (no tasks.json path for the raw data)")


def build_run(run_dir, out_dir=None, tasks_json=None, max_datasets=None) -> Path:
    """Build every figure that can be made from ``run_dir`` and write ``index.md``.

    ``run_dir`` is a ``benchmarks/real_bench.py`` output: ``rows/*.json``
    always, ``artifacts/<task>/...`` only with ``--keep-artifacts``.
    ``tasks_json`` optionally maps task id -> raw data path, for the QC and
    embeddings figures (there is otherwise no way to know where the raw data
    lives). Returns the ``out_dir`` written to.
    """
    run_dir = Path(run_dir)
    out_dir = Path(out_dir) if out_dir else run_dir / "viz"
    out_dir.mkdir(parents=True, exist_ok=True)
    made: list[dict] = []
    skipped: list[str] = []

    def note(section, path, reason, src, n=None):
        if path is not None:
            entry = {"section": section, "path": Path(path).relative_to(out_dir).as_posix(),
                     "source": src}
            if n is not None:
                entry["n"] = n
            made.append(entry)
        else:
            skipped.append(f"{section}: {reason} (from {src})")

    rows_df = _load_rows(run_dir)
    tasks_map = {}
    if tasks_json is not None and Path(tasks_json).exists():
        tasks_map = json.loads(Path(tasks_json).read_text(encoding="utf-8"))
        if isinstance(tasks_map, list):
            tasks_map = {t["id"]: t.get("path") for t in tasks_map if "id" in t}
        # paths in committed task files are written as ~/... (never a personal absolute path)
        tasks_map = {k: str(Path(v).expanduser()) if v else v for k, v in tasks_map.items()}

    if len(rows_df):
        for metric in ("auprc_ratio", "epr"):
            for arm in ("auto", "baseline"):
                path, reason = score_distributions(rows_df, metric, arm,
                                                    out_dir / f"scores_{metric}_{arm}.png")
                n_ds = int(rows_df[list(_model_cols(rows_df, arm, metric).values())]
                          .notna().any(axis=1).sum())
                note("Scores across datasets", path, reason, f"rows/*.json ({arm}/{metric})",
                    n=f"n = {n_ds} datasets")
        m_cols = [c for c in rows_df.columns if c.startswith("m_") and c.endswith("_auto_auprc_ratio")]
        models_present = sorted({c[2:-len("_auto_auprc_ratio")] for c in m_cols})
        for base, motif in _motif_pairs(models_present):
            path, reason = paired_lift(rows_df, motif, base, "auprc_ratio", "auto",
                                       out_dir / f"lift_{motif}_vs_{base}.png")
            cols = _model_cols(rows_df, "auto", "auprc_ratio")
            n_pair = int(rows_df[[cols[motif], cols[base]]].dropna().shape[0]) if motif in cols and base in cols else 0
            note("Scores across datasets", path, reason, f"rows/*.json ({motif} vs {base})",
                n=f"n = {n_pair} datasets")
    else:
        skipped.append("Scores across datasets: no rows/*.json in this run directory")

    tasks = sorted({r["task"] for r in rows_df.to_dict("records")}) if len(rows_df) else []
    if max_datasets:
        tasks = tasks[: int(max_datasets)]
    if not tasks:
        skipped.append("Per dataset: no datasets found in rows/*.json")
    for task in tasks:
        row = rows_df[rows_df["task"] == task].iloc[0].to_dict() if len(rows_df) else {}
        art = _task_artifacts(run_dir, task)
        prefix = f"Per dataset: {task}"

        md_path = out_dir / f"decisions_{task}.md"
        md_path.write_text(decisions_markdown(
            json.loads((art / "decisions.json").read_text(encoding="utf-8")) if art and
            (art / "decisions.json").exists() else row), encoding="utf-8")
        made.append({"section": prefix, "path": md_path.relative_to(out_dir).as_posix(),
                    "source": "decisions.json" if art else "rows/*.json (rules_choices)"})

        if art is None:
            skipped.append(f"{prefix}: no edge artifacts (run was made without --keep-artifacts): "
                           "no PR curves, network, TF recovery, motif lift, calibration, stability, "
                           "and QC/embeddings show the raw data only (no processed.h5ad kept)")
            _raw_qc(task, prefix, tasks_map, [], out_dir, note, skipped)
            continue

        ref_eval_path = art / "reference_evaluation.tsv"
        ref_pairs = []
        if ref_eval_path.exists():
            ref_pairs = list(pd.read_csv(ref_eval_path, sep="\t").itertuples(index=False, name=None))
        else:
            skipped.append(f"{prefix}: no reference_evaluation.tsv")

        ctx_dirs = sorted((art / "auto").glob("ctx*")) if (art / "auto").is_dir() else []
        if not ctx_dirs:
            ctx_dirs = sorted((art / "baseline").glob("ctx*")) if (art / "baseline").is_dir() else []
        edges_by_model = _load_edges_by_model(ctx_dirs[0]) if ctx_dirs else {}
        arm_used = ctx_dirs[0].parent.name if ctx_dirs else None

        if ref_pairs and edges_by_model:
            n_edges_total = max((len(e) for e in edges_by_model.values()), default=0)
            p, r = pr_curves(edges_by_model, ref_pairs, out_dir / f"pr_{task}.png")
            note(prefix, p, r, f"artifacts/{task}/{arm_used}/{ctx_dirs[0].name if ctx_dirs else ''}",
                n=f"n = {n_edges_total} edges")
            top_model = max(edges_by_model, key=lambda m: len(edges_by_model[m]))
            p, r = topk_network(edges_by_model[top_model], ref_pairs, 60, out_dir / f"network_{task}.png")
            note(prefix, p, r, f"{top_model}.tsv.gz", n="n = 60 edges")
            p, r = tf_recovery(edges_by_model, ref_pairs, 20, out_dir / f"tf_recovery_{task}.png",
                               out_dir / f"tf_recovery_{task}.tsv")
            note(prefix, p, r, f"artifacts/{task}/{arm_used}", n="n = up to 30 TFs")
            p, r = calibration(edges_by_model[top_model], ref_pairs, out_dir / f"calibration_{task}.png")
            note(prefix, p, r, f"{top_model}.tsv.gz", n=f"n = {len(edges_by_model[top_model])} edges")
            for base, motif in _motif_pairs(edges_by_model):
                p, r = motif_lift_by_tf(edges_by_model[base], edges_by_model[motif], ref_pairs,
                                        out_dir / f"motif_lift_{task}_{base}.png",
                                        out_dir / f"motif_lift_{task}_{base}.tsv")
                note(prefix, p, r, f"{base}.tsv.gz vs {motif}.tsv.gz", n="n = shared TFs")
        elif not ref_pairs:
            pass  # already noted above
        else:
            skipped.append(f"{prefix}: no edge files under artifacts/{task}/auto or baseline")

        if len(ctx_dirs) >= 2 and edges_by_model:
            for m in edges_by_model:
                lists = {}
                for cd in ctx_dirs:
                    f = cd / f"{m}.tsv.gz"
                    if f.exists():
                        lists[cd.name] = _read_edges(f)
                p, r = stability(lists, 500, out_dir / f"stability_{task}_{m}.png")
                note(prefix, p, r, f"artifacts/{task}/{arm_used}/*/{m}.tsv.gz",
                    n=f"n = {len(lists)} edge lists")
        else:
            skipped.append(f"{prefix}: no stability figure (fewer than 2 contexts with edges)")

        _raw_qc(task, prefix, tasks_map, ctx_dirs, out_dir, note, skipped)

    _write_index(out_dir, made, skipped)
    return out_dir


def _write_index(out_dir: Path, made: list[dict], skipped: list[str]) -> None:
    sections: dict[str, list[dict]] = {}
    for item in made:
        sections.setdefault(item["section"], []).append(item)
    lines = ["# grnbench run - visualization index", ""]
    for section, items in sections.items():
        lines.append(f"## {section}")
        lines.append("")
        for it in items:
            n_part = f" ({it['n']})" if it.get("n") else ""
            lines.append(f"- [{it['path']}]({it['path']}){n_part} - built from `{it['source']}`")
        lines.append("")
    lines.append("## Not produced")
    lines.append("")
    if skipped:
        for s in skipped:
            lines.append(f"- {s}")
    else:
        lines.append("- (nothing skipped)")
    lines.append("")
    (out_dir / "index.md").write_text("\n".join(lines), encoding="utf-8")
