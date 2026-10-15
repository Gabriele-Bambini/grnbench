"""Exploratory: how much of the motif-prior lift (C4, C5) is TF-level or DB-membership?

Reads the --keep-artifacts tree of a finished real_bench run; never touches its folder.
Per dataset x arm x context x declared pair (X_motif, X):

  full base edges   pearson_signed recomputed (closed form); regdiffusion / grnboost2 refit with the
                    adapter defaults (seeds 0 / 42) unless the kept file is complete (< 20,000 rows).
                    Gate: the refit's top 20,000 must equal the kept <X>.tsv.gz exactly.
  motif matrix M    cistarget.motif_ranks (the registered one), plus DB membership and K per TF.
  variants          real M and controls (below), each rank-producted with the SAME code as
                    motif_rank_product (asserted equal for real M).
  metrics           registered pooled AUPRC ratio (evaluate_ranking, universe = TF x gene) and EPR;
                    macro within-TF AUPRC ratio; both again on U'' = motif-TF x (DB & labelable) targets.
  gate              mean over contexts of the base and real-motif pooled AUPRC ratio must equal the
                    committed row (|diff| <= 1e-12), else the dataset is flagged and left out.

    python motif_sensitivity.py plan --run EXP=RUN_DIR [--run ...] --out OUT/sens_tasks.json
    python motif_sensitivity.py run  --tasks OUT/sens_tasks.json --out OUT --refs REFS --jobs 8
    python motif_sensitivity.py report --out OUT
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

CAP = int(os.environ.get("SENS_CAP", 20000))    # real_bench.ARTIFACT_TOP_EDGES (mock test: 3000)
PAIRS = {"EXP-20260925-LM-2": [("pearson_motif", "pearson_signed")],
         "EXP-20260926-LM-1": [("regdiffusion_motif", "regdiffusion"), ("grnboost2_motif", "grnboost2")]}
if os.environ.get("SENS_MOCK"):
    PAIRS = {"MOCK": [("pearson_motif", "pearson_signed"), ("grnboost2_motif", "grnboost2")]}
EXTRA = {"EXP-20260925-LM-2": ["motif_prior"], "MOCK": [m for m in os.environ.get("SENS_EXTRA", "").split(",") if m]}    # other members of the registered family
ARMS = ("auto", "baseline")
N_SEEDS = int(os.environ.get("SENS_SEEDS", 100))   # permutation seeds 0..99
MIN_POS_TF = 5
TOL = 1e-12
FIXED = ("real", "tf_mean", "tf_kexp", "tf_presence", "db_only", "worst_tie", "tfmean_db")
PERMS = ("perm_all", "perm_db")
METRICS = ("pooled", "macro", "pooled_T", "macro_T", "pooled_R", "macro_R")


# ----------------------------------------------------------------------------- helpers
def ctx_dirs(task_dir: Path, arm: str) -> list[Path]:
    d = task_dir / arm
    return sorted((p for p in d.glob("ctx*") if (p / "processed.h5ad").exists()),
                  key=lambda p: int(p.name[3:]))           # score_models order: record index i


def read_kept(path: Path) -> pd.DataFrame | None:
    return pd.read_csv(path, sep="\t", float_precision="round_trip") if path.exists() else None   # exact floats


def top_like_artifact(edges: pd.DataFrame) -> pd.DataFrame:
    top = edges.sort_values("weight", ascending=False, kind="stable").head(CAP)
    return top[["source", "target", "weight"]].reset_index(drop=True)


def rank_product(E: pd.DataFrame, m: np.ndarray) -> pd.DataFrame:
    """cistarget.motif_rank_product with the motif ranks passed in (same operations)."""
    out = E.copy()
    base_rank = out["weight"].rank(ascending=False, method="average").to_numpy()
    mot_rank = pd.Series(m).rank(ascending=True, method="average").to_numpy()
    out["weight"] = -(np.log(base_rank) + np.log(mot_rank)) / 2
    return out.sort_values("weight", ascending=False, kind="stable").reset_index(drop=True)


def motif_parts(tfs_e: list[str], genes: tuple[str, ...]):
    """M (TFs x genes, registered), gene-in-DB mask, K motifs per TF, worst rank value."""
    from grnbench.infer.cistarget import (guess_species, load_annotation, load_rankings,
                                          motif_ranks, resource_files)
    M = motif_ranks(tfs_e, genes).to_numpy()
    dbs, ann = resource_files(guess_species(genes))
    tf_motifs = load_annotation(str(ann))
    in_db, K, worst = np.zeros(len(genes), bool), np.zeros(len(tfs_e), int), None
    for db in dbs:
        motifs, matched, _, n_total = load_rankings(str(db), genes)
        in_db |= np.isin(np.asarray(genes, dtype=object), np.asarray(matched, dtype=object))
        present = set(motifs.tolist())
        K += np.array([sum(m in present for m in tf_motifs.get(t.casefold(), ())) for t in tfs_e])
        worst = n_total if worst is None else min(worst, n_total)
    return M, in_db, K, worst


def variant_matrix(name: str, M, in_db, K, worst, seed: int | None = None) -> np.ndarray:
    W = float(worst)
    has = K > 0
    if name == "real":
        return M
    if name in PERMS:
        rng = np.random.default_rng(seed)
        out = M.copy()
        cols = np.flatnonzero(in_db) if name == "perm_db" else np.arange(M.shape[1])
        for i in range(M.shape[0]):
            out[i, cols] = M[i, cols[rng.permutation(len(cols))]]
        return out
    tf_mean = np.array([M[i, in_db].mean() if has[i] and in_db.any() else W for i in range(len(K))])
    if name == "tf_mean":
        return np.repeat(tf_mean[:, None], M.shape[1], 1)
    if name == "tf_kexp":
        v = np.where(has, (W + 1) / (K + 1), W)
        return np.repeat(v[:, None], M.shape[1], 1)
    if name == "tf_presence":
        return np.repeat(np.where(has, 1.0, W)[:, None], M.shape[1], 1)
    if name == "db_only":
        return np.repeat(np.where(in_db, 1.0, W)[None, :], M.shape[0], 0)
    if name == "worst_tie":
        return np.where(M < worst, 1.0, W)
    if name == "tfmean_db":
        return np.where(has[:, None] & in_db[None, :], tf_mean[:, None], W)
    raise ValueError(name)


class Scorer:
    """Registered pooled metric via evaluate_ranking; fast exact equivalents for the rest."""

    def __init__(self, universe: pd.DataFrame, gi, masks: dict[str, np.ndarray]):
        from grnbench.eval.universe import canonical_pairs
        self.U, self.gi = universe, gi
        cu = canonical_pairs(universe, name="universe")
        self.key = pd.Series(np.arange(len(cu)), index=pd.MultiIndex.from_frame(cu))
        self.y = gi.labels(cu).astype(int)
        src = cu["source"].to_numpy()
        self.blocks = [np.flatnonzero(src == s) for s in pd.unique(src)]
        self.masks = masks          # {'T': DB & labelable targets (U'), 'R': same x motif TFs (U'')}

    def weights(self, edges: pd.DataFrame) -> np.ndarray:
        e = edges[["source", "target"]].astype(str).apply(lambda c: c.str.upper())
        pos = self.key.reindex(pd.MultiIndex.from_frame(e)).to_numpy()
        w = np.zeros(len(self.y))
        w[pos] = edges["weight"].rank(method="dense").to_numpy()     # as align_ranking
        return w

    @staticmethod
    def _ap_ratio(y, w):
        from grnbench.eval.metrics import average_precision
        if y.sum() == 0:
            return np.nan
        o = np.argsort(-w, kind="stable")
        return average_precision(y[o], w[o]) / (y.sum() / len(y))

    def macro(self, w, mask=None):
        vals = []
        for b in self.blocks:
            if mask is not None:
                b = b[mask[b]]
            if len(b) and self.y[b].sum() >= MIN_POS_TF:
                vals.append(self._ap_ratio(self.y[b], w[b]))
        return float(np.mean(vals)) if vals else np.nan

    def all(self, edges: pd.DataFrame, registered: bool = False) -> dict:
        from grnbench.eval.metrics import evaluate_ranking
        w = self.weights(edges)
        out = {"pooled": self._ap_ratio(self.y, w), "macro": self.macro(w)}
        for tag, m in self.masks.items():
            out[f"pooled_{tag}"] = self._ap_ratio(self.y[m], w[m])
            out[f"macro_{tag}"] = self.macro(w, m)
        if registered:                                           # the registered code path, as a check
            r = evaluate_ranking(edges, self.gi, with_curve=False, universe=self.U)
            assert abs(r["auprc_ratio"] - out["pooled"]) <= TOL, (r["auprc_ratio"], out["pooled"])
            out["pooled"], out["epr"] = float(r["auprc_ratio"]), float(r["epr"])
        return out


# ----------------------------------------------------------------------------- one dataset
def run_unit(unit: dict, refs: str | None, refit: bool = True) -> dict:
    import anndata as ad
    from grnbench.autoprep import resolve_strategy
    from grnbench.eval.metrics import GoldIndex
    from grnbench.infer import prepare_expression, run_method
    from grnbench.infer.cistarget import motif_rank_product
    from grnbench.tfs import resolve_tf_names

    tdir = Path(unit["run"]) / "artifacts" / unit["task"]
    row = json.loads((Path(unit["run"]) / "rows" / f"{unit['task']}.json").read_text(encoding="utf-8"))
    rd = lambda h: pd.read_csv(tdir / f"reference_{h}.tsv", sep="\t")  # noqa: E731
    sel, ev = rd("selection"), rd("evaluation")
    tfs = sorted(set(sel["source"].astype(str)) | set(ev["source"].astype(str)))   # = run_real_task tfs
    gi = GoldIndex.coerce(set(zip(ev["source"].astype(str), ev["target"].astype(str))))
    labelable = None
    if refs:
        org = json.loads(row["truth"]).get("organism") if "truth" in row else "human"   # run_real_task truth
        ref = pd.read_csv(Path(refs) / f"collectri_{org}.tsv", sep="\t")
        labelable = set(ref["target"].astype(str).str.upper())
    dec = json.loads((tdir / "decisions.json").read_text(encoding="utf-8"))
    out = {"id": unit["id"], "exp": unit["exp"], "task": unit["task"], "flags": [], "ctx": []}
    for arm in ARMS:
        choices = dec["chosen_plans"][arm]["choices"]
        assume = resolve_strategy(choices)[0] != "raw_counts"
        for cdir in ctx_dirs(tdir, arm):
            adata = ad.read_h5ad(cdir / "processed.h5ad")
            genes = pd.Index(adata.var_names).astype(str)
            sources = resolve_tf_names(genes, tfs)
            U = pd.DataFrame([(s, g) for s in sources for g in genes if s != g], columns=["source", "target"])
            expr, _ = prepare_expression(adata, assume_normalized=assume)
            for m in EXTRA.get(unit["exp"], []):                        # e.g. motif_prior (C4 family)
                from grnbench.eval.metrics import evaluate_ranking
                from grnbench.infer.cistarget import motif_prior
                r = evaluate_ranking(motif_prior(expr, tf_names=tfs), gi, with_curve=False, universe=U)
                out["ctx"].append({"arm": arm, "ctx": cdir.name, "pair": m, "extra_pooled": float(r["auprc_ratio"]),
                                   "extra_epr": float(r["epr"])})
            for combo, base in unit["pairs"]:
                kept, kept_m = read_kept(cdir / f"{base}.tsv.gz"), read_kept(cdir / f"{combo}.tsv.gz")
                rec = {"arm": arm, "ctx": cdir.name, "pair": combo, "n_universe": len(U)}
                if kept is None or kept_m is None:                   # model errored (rule c)
                    rec["missing"] = True
                    out["ctx"].append(rec)
                    continue
                if base == "pearson_signed" or (refit and len(kept) >= CAP):
                    E, _ = run_method(base, adata, tf_list=tfs, assume_normalized=assume)
                    same = top_like_artifact(E).equals(kept[["source", "target", "weight"]])
                    rec["edges"] = "recomputed" if same else "RECOMPUTED_MISMATCH"
                else:
                    E = kept
                    rec["edges"] = "kept_complete" if len(kept) < CAP else "kept_TRUNCATED"
                genes_t = tuple(map(str, expr.columns))
                tfs_e = list(dict.fromkeys(E["source"].astype(str)))
                M, in_db, K, worst = motif_parts(tfs_e, genes_t)
                ti = E["source"].astype(str).map({t: i for i, t in enumerate(tfs_e)}).to_numpy()
                gix = E["target"].astype(str).map({g: i for i, g in enumerate(genes_t)}).to_numpy()
                ref_comb = motif_rank_product(E, expr)
                mine = rank_product(E, M[ti, gix].astype(float))
                assert mine[["source", "target"]].equals(ref_comb[["source", "target"]]) and \
                    np.array_equal(mine["weight"].to_numpy(), ref_comb["weight"].to_numpy())
                rec["combo_top_matches_kept"] = bool(top_like_artifact(mine).equals(kept_m[["source", "target", "weight"]]))
                if rec["edges"] != "recomputed" and rec["edges"] != "kept_complete" or not rec["combo_top_matches_kept"]:
                    out["flags"].append(f"EDGES {arm}/{cdir.name}/{combo}: {rec['edges']}, "
                                        f"combo top matches kept: {rec['combo_top_matches_kept']}")
                # U'' = TFs with >= 1 motif x targets in the DB and labelable in the full reference
                up = {g.upper(): i for i, g in enumerate(genes_t)}
                tf_has = {t.upper(): k > 0 for t, k in zip(tfs_e, K)}
                cu_src, cu_tgt = U["source"].str.upper().to_numpy(), U["target"].str.upper().to_numpy()
                tgt_ok = np.array([in_db[up[t]] and (labelable is None or t in labelable) for t in cu_tgt])
                restrict = tgt_ok & np.array([tf_has.get(s, False) for s in cu_src])
                sc = Scorer(U, gi, {"T": tgt_ok, "R": restrict})
                rec.update({f"base_{k}": v for k, v in sc.all(E, registered=True).items()})
                for v in FIXED:
                    Mv = variant_matrix(v, M, in_db, K, worst)
                    res = sc.all(rank_product(E, Mv[ti, gix].astype(float)), registered=(v == "real"))
                    rec.update({f"{v}_{k}": x for k, x in res.items()})
                for v in PERMS:
                    per = [sc.all(rank_product(E, variant_matrix(v, M, in_db, K, worst, s)[ti, gix].astype(float)))
                           for s in range(N_SEEDS)]
                    for k in per[0]:
                        rec[f"{v}_{k}"] = float(np.nanmean([p[k] for p in per]))
                        rec[f"{v}_{k}_seeds"] = [p[k] for p in per]
                y_u = sc.y
                worst_u = np.array([not tf_has.get(s, False) or not in_db[up[t]] for s, t in zip(cu_src, cu_tgt)])
                rec.update(n_tfs=len(tfs_e), n_tfs_no_motif=int((K == 0).sum()), genes_not_in_db=int((~in_db).sum()),
                           genes_ensembl=int(pd.Index(genes_t).str.match(r"^ENS[A-Z]{0,4}G\d+").sum()),
                           frac_universe_worst_tie=float(worst_u.mean()),
                           gold_in_worst_tie=int(y_u[worst_u].sum()), gold_in_universe=int(y_u.sum()),
                           frac_universe_R=float(restrict.mean()), gold_in_R=int(y_u[restrict].sum()),
                           frac_universe_T=float(tgt_ok.mean()), gold_in_T=int(y_u[tgt_ok].sum()),
                           spearman_K_gold=float(pd.Series(K).corr(pd.Series(
                               [int(gi.labels(U[U.source == t]).sum()) for t in tfs_e]), method="spearman")))
                out["ctx"].append(rec)
    # dataset level: mean over contexts (finite only), as score_models; then the gate vs committed rows
    df = pd.DataFrame([r for r in out["ctx"] if not r.get("missing") and "extra_pooled" not in r])
    for (arm, combo), g in (df.groupby(["arm", "pair"]) if len(df) else []):
        base = dict((c, b) for c, b in unit["pairs"])[combo]
        for col in [c for c in g if c.endswith(METRICS + ("epr",))]:
            x = [v for v in g[col].astype(float) if np.isfinite(v)]          # score_models: np.mean of the finite
            out[f"{arm}|{combo}|{col}"] = float(np.mean(x)) if x else np.nan
        for mine, reg in ((f"{arm}|{combo}|base_pooled", f"m_{base}_{arm}_auprc_ratio"),
                          (f"{arm}|{combo}|real_pooled", f"m_{combo}_{arm}_auprc_ratio")):
            d = abs(out[mine] - float(row.get(reg, np.nan)))
            out[f"gate|{mine}"] = d
            if not d <= TOL:
                out["flags"].append(f"GATE FAIL {mine} vs {reg}: |diff| = {d}")
    for arm in ARMS:
        for m in EXTRA.get(unit["exp"], []):
            x = [r["extra_pooled"] for r in out["ctx"] if r["arm"] == arm and r["pair"] == m]
            x = [v for v in x if np.isfinite(v)]
            out[f"{arm}|{m}|extra_pooled"] = float(np.mean(x)) if x else np.nan
            y = [r["extra_epr"] for r in out["ctx"] if r["arm"] == arm and r["pair"] == m]
            y = [v for v in y if np.isfinite(v)]
            out[f"{arm}|{m}|extra_epr"] = float(np.mean(y)) if y else np.nan
            d = abs(out[f"{arm}|{m}|extra_pooled"] - float(row.get(f"m_{m}_{arm}_auprc_ratio", np.nan)))
            out[f"gate|{arm}|{m}"] = d
            if not d <= TOL:
                out["flags"].append(f"GATE FAIL {arm}|{m}: |diff| = {d}")
    return out


def _run_and_save(args):
    unit, refs, outdir, refit = args
    path = Path(outdir) / "rows" / f"{unit['id']}.json"
    if path.exists():
        return unit["id"], "skipped"
    res = run_unit(unit, refs, refit)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(res, default=float), encoding="utf-8")
    return unit["id"], "; ".join(res["flags"]) or "ok"


# ----------------------------------------------------------------------------- report
def report(outdir: Path) -> str:
    from grnbench.eval.comparison import compare_units, holm_adjust
    rows = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((outdir / "rows").glob("*.json"))]
    bad = [r["id"] for r in rows if r["flags"]]
    lines = ["# Motif-prior lift decomposition (EXPLORATORY; unit = dataset)", "",
             f"- datasets: {len(rows)}; gate failures (left out): {len(bad)} {bad}", ""]
    rows = [r for r in rows if not r["flags"]]
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from real_bench import models_section                     # the registered report code, unchanged
    for exp, pairs in PAIRS.items():
        recs = []
        for r in (r for r in rows if r["exp"] == exp):
            d = {"task": r["task"]}
            for arm in ARMS:
                for combo, base in pairs:
                    for mdl, v in ((base, "base"), (combo, "real")):
                        d[f"m_{mdl}_{arm}_auprc_ratio"] = r.get(f"{arm}|{combo}|{v}_pooled", np.nan)
                        d[f"m_{mdl}_{arm}_epr"] = r.get(f"{arm}|{combo}|{v}_epr", np.nan)
                        d[f"m_{mdl}_{arm}_seconds"] = np.nan
                for m in EXTRA.get(exp, []):
                    d[f"m_{m}_{arm}_auprc_ratio"] = r.get(f"{arm}|{m}|extra_pooled", np.nan)
                    d[f"m_{m}_{arm}_epr"] = r.get(f"{arm}|{m}|extra_epr", np.nan)
                    d[f"m_{m}_{arm}_seconds"] = np.nan
            recs.append(d)
        if len(recs) >= 3:
            sec = models_section(pd.DataFrame(recs))
            keep = [l for l in sec if l.startswith(("####", "| A ", "|---|---|---|---|---|---|---|---|"))
                    or (l.startswith("| ") and l.count("|") == 9)]
            lines += [f"## Check (i): {exp} declared comparisons recomputed (must equal the committed REPORT_REAL.md lines)", ""] + keep + [""]
    boot = np.random.default_rng(0)
    primary = []
    for exp, pairs in PAIRS.items():
        R = [r for r in rows if r["exp"] == exp]
        for arm in ARMS:
            def col(combo, v, metric):
                return pd.Series({r["task"]: r.get(f"{arm}|{combo}|{v}_{metric}", np.nan) for r in R}, dtype=float)
            lines += [f"## {exp}, {arm} arm: lift = variant - base (paired t over datasets; Holm over the "
                      f"{2 * len(pairs)} directed hypotheses of each variant x metric, as models_section)", "",
                      "| metric | variant | pair | n | mean lift | 95% interval | p(>0) Holm | p(<0) Holm | "
                      "share of real lift [bootstrap 95%] |", "|---|---|---|---|---|---|---|---|---|"]
            for metric in METRICS:
                for v in FIXED + PERMS:
                    recs, raw = [], []
                    for combo, base in pairs:
                        pair = pd.DataFrame({"v": col(combo, v, metric), "b": col(combo, "base", metric),
                                             "r": col(combo, "real", metric)}).dropna()
                        if len(pair) < 3:
                            continue
                        res = compare_units(pair[["v", "b"]])["pairs"][0]
                        recs.append((combo, pair, res)); raw += [res["p_a_gt_b"], res["p_b_gt_a"]]
                    adj = holm_adjust(raw) if raw else []
                    for i, (combo, pair, res) in enumerate(recs):
                        lo, hi = res["simultaneous_ci"]
                        ci = "unbounded" if lo is None else f"[{lo:+.3f}, {hi:+.3f}]"
                        lv, lr = (pair.v - pair.b).to_numpy(), (pair.r - pair.b).to_numpy()
                        idx = boot.integers(0, len(lv), (2000, len(lv)))
                        bs = lv[idx].mean(1) / lr[idx].mean(1)
                        share = f"{lv.mean() / lr.mean():.2f} [{np.nanpercentile(bs, 2.5):.2f}, {np.nanpercentile(bs, 97.5):.2f}]" \
                            if lr.mean() != 0 else "-"
                        lines.append(f"| {metric} | {v} | {combo} | {len(pair)} | {res['difference']:+.3f} | {ci} | "
                                     f"{adj[2*i]:.3g} | {adj[2*i+1]:.3g} | {share} |")
                        if arm == "auto" and metric == "pooled" and v == "perm_db":
                            primary.append((exp, combo, "T1 real > perm_db (pooled)", pair[["r", "v"]]))
                        if arm == "auto" and metric == "macro_R" and v == "real":
                            primary.append((exp, combo, "T2 real > base (macro, U'')", pair[["r", "b"]]))
            lines.append("")
    lines += ["## Declared primary tests of this analysis (auto arm; paired t, Holm over all rows x 2 directions)", "",
              "| experiment | pair | test | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm |",
              "|---|---|---|---|---|---|---|---|"]
    recs, raw = [], []
    for exp, combo, name, pair in primary:
        res = compare_units(pair)["pairs"][0]
        recs.append((exp, combo, name, len(pair), res)); raw += [res["p_a_gt_b"], res["p_b_gt_a"]]
    adj = holm_adjust(raw) if raw else []
    for i, (exp, combo, name, n, res) in enumerate(recs):
        lo, hi = res["simultaneous_ci"]
        ci = "unbounded" if lo is None else f"[{lo:+.4f}, {hi:+.4f}]"
        lines.append(f"| {exp} | {combo} | {name} | {n} | {res['difference']:+.4f} | {ci} | {adj[2*i]:.3g} | {adj[2*i+1]:.3g} |")
    text = "\n".join(lines) + "\n"
    (outdir / "SENSITIVITY.md").write_text(text, encoding="utf-8")
    return text


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("plan"); a.add_argument("--run", action="append", required=True); a.add_argument("--out", required=True)
    b = sub.add_parser("run")
    for k in ("--tasks", "--out"):
        b.add_argument(k, required=True)
    b.add_argument("--refs"); b.add_argument("--jobs", type=int, default=1)
    b.add_argument("--index", type=int); b.add_argument("--chunk", type=int, default=1)
    b.add_argument("--workroot"); b.add_argument("--models"); b.add_argument("--keep-artifacts", action="store_true")
    b.add_argument("--no-refit", action="store_true", help="Tier A: kept top-20,000 edges only (truncated)")
    c = sub.add_parser("report"); c.add_argument("--out", required=True)
    args = p.parse_args(argv)
    if args.cmd == "plan":
        import hashlib
        units = []
        for spec in args.run:
            exp, run = spec.split("=", 1)
            repo_root = next(  # walk up from this file to find the repo root (has pyproject.toml)
                d for d in Path(__file__).resolve().parents if (d / "pyproject.toml").exists()
            )
            sums = repo_root / "benchmarks" / "results" / exp / "SHA256SUMS"
            committed = dict(l.split(" *./")[::-1] for l in sums.read_text().splitlines()) if sums.exists() else {}
            for r in sorted((Path(run).expanduser() / "rows").glob("*.json")):
                if committed and committed.get(f"rows/{r.name}") != hashlib.sha256(r.read_bytes()).hexdigest():
                    raise SystemExit(f"{r}: not the committed row of {exp} (SHA256SUMS)")   # gate reference = committed
                units.append({"id": f"{exp}__{r.stem}", "exp": exp, "task": r.stem,
                              "run": str(Path(run).expanduser()), "pairs": PAIRS[exp]})
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(units, indent=0), encoding="utf-8")
        print(f"{len(units)} units -> {args.out}")
        return 0
    if args.cmd == "report":
        print(report(Path(args.out)))
        return 0
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    units = json.loads(Path(args.tasks).read_text(encoding="utf-8"))
    if args.index is not None:
        units = units[args.index * args.chunk:(args.index + 1) * args.chunk]
    jobs = [(u, args.refs, args.out, not args.no_refit) for u in units]
    if args.jobs > 1:
        from multiprocessing import get_context
        with get_context("spawn").Pool(args.jobs, maxtasksperchild=1) as pool:
            for uid, status in pool.imap_unordered(_run_and_save, jobs):
                print(f"{uid}: {status}", flush=True)
    else:
        for j in jobs:
            print("%s: %s" % _run_and_save(j), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
