"""Knockdown validation for EXP-20261016-LM-1, phase 2: network inference and scoring.

    python benchmarks/knockdown_score.py infer --controls controls.h5ad --tfs collectri.tsv \\
        --experiment EXP-ID --models LIST --out DIR [--seed 0] [--max-cells 3000] [--attempts FILE]
    python benchmarks/knockdown_score.py score --infer DIR --truth DIR --out DIR [--perms 10000] [--seed 0] \\
        [--collectri collectri.tsv] [--essential-from TRUTH_DIR] [--both-with TRUTH_DIR] [--total-shift TSV]
    python benchmarks/knockdown_score.py combine --design DIR --replication DIR --out DIR
    python benchmarks/knockdown_score.py run-identity --run DIR --snapshot TSV --perms N --norm core [--models LIST]

`infer` (record section 1, "Inference"): the engine's fixed baseline (`screen.BASELINE_CHOICES`) with
`hvg = all_genes` on the control cells in `controls.h5ad` (the inference half; no perturbed cell is read), then
the cell cap AFTER the pipeline's QC with the `real_bench` rule `np.sort(default_rng(seed).choice(n, cap,
replace=False))` (cap = `engine_bench.MAX_SCORE_CELLS` = 3,000, seed 0; a no-op when QC leaves <= 3,000 cells).
Regulators = every CollecTRI source present in the processed genes. Each model of `--models` runs through the code
`real_bench.score_models` uses (`run_method`, `_split_combined`, `cistarget.motif_rank_product`); `X_motif_genetf`
is `X`'s edges rank-producted with `motif_sensitivity.variant_matrix("gene_tf")` by the same rank-product rule
(`motif_sensitivity.rank_product`, checked equal to the registered `motif_rank_product` on the unchanged prior).
Writes `edges/<arm>.edges.tsv.gz` (source, target, weight), `genes.txt`, `regulators.tsv` and `infer.json`
(genes, cells before and after QC and cap, regulators, per-arm seconds, commit, skipped arms and why). A model
that cannot run is recorded in `infer.json["skipped"]`; its `_motif` and `_genetf` arms are skipped with it.
Record exclusion (c), one resubmission: an inference with a skipped arm is `status = "incomplete"` and exits 3, so
the job publishes nothing. `--attempts FILE` keeps the failures per arm and commit across resubmissions; the
inference is `complete` with arms still missing only once every skipped arm has failed twice on one commit. A failure
is a caught error or an arm left `started` in the file by a run that was killed (Slurm timeout): the file is written
(atomic replace) before each arm starts, and the next run of the commit counts the unfinished arm as one failure.

`score` runs `grnbench.eval.knockdown` on every arm against the usable units' responders (auprc_ratio is the
primary endpoint; top-50 and top-200 hit rates are secondary, each with its own Holm family) and writes
`scores.tsv`, `family.json` and `REPORT.md`. The family is always the registered 12 hypotheses / 24 directions
(plus any extra model requested): one that needs a missing arm has p = 1 and verdict `missing`.
`combine` joins the design and replication `family.json` with
`knockdown.combine_sets` and maps the result to the record's decision rule. `run-identity` writes `run.json` in a
run folder (commit, snapshot hash, models, permutations, normalisation) and refuses a folder whose `run.json` or
finished steps do not match this job.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import platform
import sys
import time
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import motif_sensitivity as ms  # noqa: E402  (imports real_bench and engine_bench)
import real_bench  # noqa: E402
from engine_bench import MAX_SCORE_CELLS  # noqa: E402

from grnbench.eval import knockdown as kd  # noqa: E402
from grnbench.eval.metrics import expected_average_precision  # noqa: E402
from grnbench.infer import EDGE_COLUMNS, cistarget, prepare_expression, run_method, save_edges  # noqa: E402
from grnbench.tfs import resolve_tf_names  # noqa: E402

EXPERIMENT = "EXP-20261016-LM-1"
CAP_SEED = 0
DEFAULT_MODELS = ("pearson_signed,pearson_motif,pearson_motif_genetf,grnboost2,grnboost2_motif,"
                  "grnboost2_motif_genetf,regdiffusion,regdiffusion_motif,regdiffusion_motif_genetf")
ENDPOINTS = (("auprc_ratio", "auprc_ratio", 50), ("top50", "precision_at_k", 50), ("top200", "precision_at_k", 200))
REFERENCE_ARM = "collectri"          # descriptive reference row, outside the family


# ---------------------------------------------------------------- arms
def parse_arm(name: str) -> tuple[str, str]:
    """``'pearson_motif_genetf'`` -> ``('pearson_signed', 'genetf')``; kind is base, motif or genetf."""
    core = name.removesuffix(kd.GENETF)
    base, combiner = real_bench._split_combined(core)
    if combiner not in (None, "_motif") or (name != core and combiner != "_motif") or not base:
        raise ValueError(f"{name!r}: arms are X, X_motif and X_motif_genetf")
    return base, ("genetf" if name != core else "motif" if combiner else "base")


# ---------------------------------------------------------------- preprocessing
def preprocess(controls, workdir, sources=None) -> dict:
    """Fixed baseline with ``hvg = all_genes``; the plan gets the reference TFs present in the genes as
    ``tf_list``, as ``real_bench`` does for C4-C7 (``knockdown_bench.baseline_check`` plans without it)."""
    from grnbench.autoprep import diagnose, execute_plan, plan_preprocessing
    from grnbench.profile import Profile
    from grnbench.screen import BASELINE_CHOICES
    from grnbench.selector import select

    t0 = time.time()
    diag = diagnose(controls)
    cells_in, genes_in = int(diag.adata.n_obs), int(diag.adata.n_vars)
    tf_list = resolve_tf_names([str(g) for g in diag.adata.var_names], sources) or None if sources else None
    selection = select(Profile(**plan_preprocessing(diag, tf_list=tf_list).profile))
    plan = plan_preprocessing(diag, tf_list=tf_list, selection=selection,
                              choices={**selection.chosen(), **BASELINE_CHOICES, "hvg": "all_genes"})
    records = execute_plan(plan, controls, workdir)
    failed = [r for r in records if r["error"] or not r["processed"]]
    if failed or len(records) != 1:
        raise RuntimeError("baseline preprocessing did not give one processed context: "
                           + "; ".join(f"{r['context']}: {r['error']}" for r in records))
    adata = ad.read_h5ad(records[0]["processed"])
    return {"adata": adata, "plan": plan, "cells_in": cells_in, "genes_in": genes_in,
            "seconds": round(time.time() - t0, 1)}


def cap_cells(adata, max_cells: int, seed: int):
    """The ``real_bench.score_models`` cap: ``np.sort(default_rng(seed).choice(n, max_cells, replace=False))``."""
    if adata.n_obs <= max_cells:
        return adata, False
    keep = np.sort(np.random.default_rng(seed).choice(adata.n_obs, max_cells, replace=False))
    return adata[keep].copy(), True


# ---------------------------------------------------------------- derived arms
def genetf_arm(E: pd.DataFrame, expr: pd.DataFrame, motif_edges: pd.DataFrame) -> pd.DataFrame:
    """``E`` rank-producted with the TF-generic ``gene_tf`` propensity of the registered motif matrix.

    Same operations as ``motif_sensitivity.decompose``: M, in-database mask, K and the worst rank from
    ``motif_parts`` over the regulators present in ``E``; the rank product with the unchanged M must equal the
    registered ``motif_rank_product`` output (``motif_edges``) exactly, else the arm is not written."""
    genes = tuple(map(str, expr.columns))
    tfs_e = list(dict.fromkeys(E["source"].astype(str)))
    M, in_db, K, worst = ms.motif_parts(tfs_e, genes)
    ti = E["source"].astype(str).map({t: i for i, t in enumerate(tfs_e)}).to_numpy().astype(int)
    gix = E["target"].astype(str).map({g: i for i, g in enumerate(genes)}).to_numpy().astype(int)
    mine = ms.rank_product(E, M[ti, gix].astype(float))
    if not (mine[["source", "target"]].equals(motif_edges[["source", "target"]])
            and np.array_equal(mine["weight"].to_numpy(), motif_edges["weight"].to_numpy())):
        raise ms.CheckFailed("self-consistency failed: rank_product != cistarget.motif_rank_product")
    del mine
    return ms.rank_product(E, ms.variant_matrix("gene_tf", M, in_db, K, worst)[ti, gix].astype(float))


def _write_meta(out: Path, meta: dict) -> None:
    (out / "infer.json").write_text(json.dumps(meta, indent=1, default=str) + "\n", encoding="utf-8")


def _attempt_state(attempts, commit: str) -> dict:
    """The attempts file for this commit: ``failures`` per arm and ``arms``, the log of every arm's last attempt
    (``started_at``, ``status`` started | done | failed). A file of another commit counts for nothing."""
    path = Path(attempts) if attempts else None
    prev = json.loads(path.read_text(encoding="utf-8")) if path and path.exists() else {}
    if prev.get("code_commit") != commit:
        prev = {}
    return {"code_commit": commit, "failures": dict(prev.get("failures", {})), "arms": dict(prev.get("arms", {}))}


def _store_attempts(attempts, state: dict) -> None:
    """Atomic (temp file + os.replace): a kill during the write leaves the previous file whole."""
    if not attempts:
        return
    path = Path(attempts)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(state, indent=1) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def _mark_arms(attempts, state: dict, arms, status: str) -> None:
    now = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    for arm in arms:
        old = state["arms"].get(arm, {})
        state["arms"][arm] = {"code_commit": state["code_commit"], "status": status,
                              "started_at": now if status == "started" else old.get("started_at")}
    _store_attempts(attempts, state)


def record_failures(meta: dict, attempts) -> str:
    """``complete`` when no arm was skipped or every skipped arm has now failed twice on this commit (record
    exclusion (c)); else ``incomplete``, with the failures per arm in ``meta["attempt"]``. A failure is a caught
    error here or an arm left ``started`` by an earlier run that was killed (counted by ``infer`` at its start);
    an arm skipped unrun because it already failed twice is not counted again. ``attempts`` is a json file that
    outlives the job; without it nothing can be counted, so any skipped arm is incomplete."""
    if not meta["skipped"]:
        return "complete"
    state = _attempt_state(attempts, meta["code_commit"])
    meta["attempt"] = {a: state["failures"].get(a, 0) + (a not in meta.get("not_run", ())) for a in meta["skipped"]}
    state["failures"] = {**state["failures"], **meta["attempt"]}
    if attempts:
        _store_attempts(attempts, state)
    return "complete" if min(meta["attempt"].values()) >= 2 else "incomplete"


def infer(controls, tfs_path, models, out, seed: int = CAP_SEED, max_cells: int = MAX_SCORE_CELLS,
          device: str = "cpu", attempts=None, *, experiment: str) -> dict:
    """Run the preprocessing and every model; returns the ``infer.json`` dict."""
    t_start = time.time()
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    wanted = list(dict.fromkeys(m.strip() for m in models if m.strip()))
    parsed = {m: parse_arm(m) for m in wanted}            # validates before any work
    sources = pd.read_csv(tfs_path, sep="\t", dtype=str)["source"].dropna().unique().tolist()

    src = ad.read_h5ad(controls, backed="r")
    source_h5ad = str(src.uns["source_h5ad"]) if "source_h5ad" in src.uns else None
    src.file.close()
    prep = preprocess(controls, out / "prep", sources)
    adata, plan = prep["adata"], prep["plan"]
    genes = [str(g) for g in adata.var_names]
    cells_qc = int(adata.n_obs)
    adata, capped = cap_cells(adata, max_cells, seed)
    regs = resolve_tf_names(genes, sources)
    assume = plan.value("strategy") != "raw_counts"
    expr, norm_note = prepare_expression(adata, assume_normalized=assume)
    (out / "genes.txt").write_text("\n".join(genes) + "\n", encoding="utf-8")
    (out / "cells.txt").write_text("\n".join(map(str, adata.obs_names)) + "\n", encoding="utf-8")
    X = adata.X
    det = np.asarray((X[:, [adata.var_names.get_loc(r) for r in regs]] > 0).mean(axis=0)).ravel()
    pd.DataFrame({"tf": regs, "detect_frac": det}).to_csv(out / "regulators.tsv", sep="\t", index=False)

    meta = {
        "experiment": experiment, "status": "running", "code_commit": ms.code_commit(),
        "controls": Path(controls).name, "source_h5ad": source_h5ad, "tfs_file": Path(tfs_path).name,
        "tfs_sha256": hashlib.sha256(Path(tfs_path).read_bytes()).hexdigest(),
        "preprocessing": {"choices": dict(plan.choices), "strategy": plan.value("strategy"),
                          "normalization_note": norm_note, "seconds": prep["seconds"]},
        "cells_in": prep["cells_in"], "cells_after_qc": cells_qc, "max_cells": int(max_cells),
        "cap_applied": bool(capped), "cells_after_cap": int(adata.n_obs), "cap_seed": int(seed),
        "cap_rule": "after QC: np.sort(np.random.default_rng(seed).choice(n, max_cells, replace=False)) "
                    "(real_bench.score_models); no-op when n <= max_cells",
        "genes_in": prep["genes_in"], "n_genes": len(genes), "n_regulators": len(regs),
        "tfs_in_reference": len(sources), "regulators": regs,
        "models_requested": wanted, "regdiffusion_device": device, "arms": {}, "skipped": {},
        "versions": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
                     "anndata": version("anndata")},
    }
    _write_meta(out, meta)
    state = _attempt_state(attempts, meta["code_commit"])
    dead = sorted(a for a, e in state["arms"].items() if e.get("status") == "started")
    if dead:                                              # an earlier run died inside these arms (Slurm timeout)
        for a in dead:
            state["failures"][a] = state["failures"].get(a, 0) + 1
            state["arms"][a]["status"] = "failed"
        meta["interrupted"] = dead
        _store_attempts(attempts, state)

    def save(arm, E, seconds, base, kind, base_seconds):
        path = save_edges(E[EDGE_COLUMNS], out / "edges", arm)
        meta["arms"][arm] = {"file": path.relative_to(out).as_posix(), "base": base, "kind": kind,
                             "n_edges": int(len(E)), "n_sources": int(E["source"].nunique()),
                             "regulators_without_edges": int(len(set(regs) - set(E["source"]))),
                             "seconds": round(seconds, 1), "seconds_including_base": round(seconds + base_seconds, 1)}
        _write_meta(out, meta)
        print(f"{arm}: {len(E)} edges, {seconds:.1f}s  (elapsed {time.time() - t_start:.0f}s)", flush=True)

    bases = list(dict.fromkeys(b for b, _ in parsed.values()))
    for bi, base in enumerate(bases, 1):
        arms = {kind: arm for arm, (b, kind) in parsed.items() if b == base}
        if all(state["failures"].get(arm, 0) >= 2 for arm in arms.values()):    # nothing left to try
            why = "did not finish within the job (2 attempts on this commit)"
            meta["skipped"].update({arm: why for arm in arms.values()})
            meta.setdefault("not_run", []).extend(arms.values())
            _write_meta(out, meta)
            print(f"[{bi}/{len(bases)}] {base} not run: {why}", flush=True)
            continue
        _mark_arms(attempts, state, arms.values(), "started")
        t0 = time.time()
        try:
            E, _ = run_method(base, adata, tf_list=regs, assume_normalized=assume,
                              **({"device": device} if base == "regdiffusion" else {}))
        except Exception as exc:  # noqa: BLE001 - a model that cannot run is recorded, never dropped silently
            why = f"{type(exc).__name__}: {exc}"[:500]
            meta["skipped"].update({arm: why if kind == "base" else f"base model {base} skipped: {why}"
                                    for kind, arm in arms.items()})
            _write_meta(out, meta)
            _mark_arms(attempts, state, arms.values(), "failed")
            print(f"[{bi}/{len(bases)}] {base} skipped: {why}", flush=True)
            continue
        base_s = time.time() - t0
        if "base" in arms:
            save(arms["base"], E, base_s, base, "base", 0.0)
        motif = None
        if "motif" in arms or "genetf" in arms:
            t0 = time.time()
            try:
                motif = cistarget.motif_rank_product(E, expr)
            except Exception as exc:  # noqa: BLE001
                why = f"{type(exc).__name__}: {exc}"[:500]
                meta["skipped"].update({arms[k]: f"motif prior unavailable: {why}" for k in ("motif", "genetf") if k in arms})
            else:
                if "motif" in arms:
                    save(arms["motif"], motif, time.time() - t0, base, "motif", base_s)
        if motif is not None and "genetf" in arms:
            t0 = time.time()
            try:
                save(arms["genetf"], genetf_arm(E, expr, motif), time.time() - t0, base, "genetf", base_s)
            except Exception as exc:  # noqa: BLE001
                meta["skipped"][arms["genetf"]] = f"{type(exc).__name__}: {exc}"[:500]
        _write_meta(out, meta)
        for arm in arms.values():
            _mark_arms(attempts, state, [arm], "failed" if arm in meta["skipped"] else "done")
        del E, motif
        print(f"[{bi}/{len(bases)}] {base} done", flush=True)
    meta["status"] = record_failures(meta, attempts)
    meta["seconds_total"] = round(time.time() - t_start, 1)
    _write_meta(out, meta)
    return meta


# ---------------------------------------------------------------- scoring
def _sha(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _sha_content(path: Path) -> str:
    """sha256 of the decompressed text: a gzip header carries the write time, the content does not."""
    h = hashlib.sha256()
    with gzip.open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_truth(truth: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Usable units and their responders; the tables must match the hashes written by ``truth``."""
    truth = Path(truth)
    if not (truth / "truth.json").exists():
        raise SystemExit(f"{truth / 'truth.json'} is missing: the truth tables cannot be checked against their hashes")
    meta = json.loads((truth / "truth.json").read_text(encoding="utf-8"))
    for name, digest in (meta.get("sha256") or {}).items():
        if _sha(truth / name) != digest:
            raise ValueError(f"{truth / name} does not match the sha256 in truth.json: the screen is not analysed")
    units = pd.read_csv(truth / "units.tsv", sep="\t", keep_default_na=False)
    usable = units[units["usable"].astype(str) == "True"]
    resp = pd.read_csv(truth / "responders.tsv", sep="\t")
    return usable, resp[resp.tf.isin(usable.tf)], meta


def perturbed_tfs(truth: Path) -> set[str]:
    return set(pd.read_csv(Path(truth) / "units.tsv", sep="\t", keep_default_na=False)["tf"].astype(str).str.upper())


def load_arms(infer_dir: Path) -> tuple[dict, dict[str, pd.DataFrame], list[str]]:
    meta = json.loads((infer_dir / "infer.json").read_text(encoding="utf-8"))
    if meta.get("status") != "complete":
        raise ValueError(f"{infer_dir / 'infer.json'}: inference is not complete (status {meta.get('status')!r})")
    edges = {a: pd.read_csv(infer_dir / v["file"], sep="\t") for a, v in meta["arms"].items()}
    genes = [g for g in (infer_dir / "genes.txt").read_text(encoding="utf-8").split() if g]
    return meta, edges, genes


def reference_edges(path) -> pd.DataFrame:
    ref = pd.read_csv(path, sep="\t", dtype=str)[["source", "target"]].dropna().drop_duplicates()
    ref = ref[ref.source.str.upper() != ref.target.str.upper()]
    ref = ref[~ref.assign(s=ref.source.str.upper(), t=ref.target.str.upper()).duplicated(["s", "t"])]
    return ref.assign(weight=1.0).reset_index(drop=True)


def run_endpoints(edges, resp, genes, perms, seed, family, secondary, strata, missing_arms=()):
    """Evaluate every endpoint on all units and (primary endpoint only) on each stratum."""
    out = {"endpoints": {}, "strata": {}}
    for key, metric, k in ENDPOINTS:
        out["endpoints"][key] = kd.evaluate(edges, resp, genes, metric=metric, k=k, n_perm=perms, seed=seed,
                                            family=family, secondary=secondary, missing_arms=missing_arms)
        print(f"{key}: {out['endpoints'][key]['n_units']} units", flush=True)
    for name, members in strata.items():
        sub = resp[resp.tf.str.upper().isin(members)]
        try:
            out["strata"][name] = kd.evaluate(edges, sub, genes, metric="auprc_ratio", n_perm=perms, seed=seed,
                                              family=family, secondary=secondary, missing_arms=missing_arms)
        except ValueError as exc:            # fewer than three units: reported, not hidden
            out["strata"][name] = {"skipped": str(exc), "per_unit": None}
    return out


def detection(per_unit: pd.DataFrame, regs: pd.DataFrame) -> dict:
    """Spearman between a unit TF's detection in the inference cells and its per-TF score, per arm."""
    frac = regs.assign(tf=regs.tf.str.upper()).set_index("tf").detect_frac
    d = frac.reindex(per_unit.index)
    ok = d.notna()
    out = {}
    for arm in per_unit.columns.drop("n_responders"):
        v = per_unit.loc[ok, arm]
        out[arm] = {"n": int(ok.sum()), "spearman": None if v.nunique() < 2 or d[ok].nunique() < 2
                    else float(spearmanr(d[ok], v)[0])}
    return out


def long_scores(results: dict) -> pd.DataFrame:
    frames = []
    for key, res in results["endpoints"].items():
        frames.append(("all", key, res["per_unit"]))
    for name, res in results["strata"].items():
        if res.get("per_unit") is not None:
            frames.append((name, "auprc_ratio", res["per_unit"]))
    rows = []
    for stratum, endpoint, pu in frames:
        long = pu.reset_index().melt(id_vars=["tf", "n_responders"], var_name="arm", value_name="value")
        rows.append(long.assign(endpoint=endpoint, stratum=stratum))
    out = pd.concat(rows, ignore_index=True)[["endpoint", "stratum", "tf", "n_responders", "arm", "value"]]
    return out.sort_values(["endpoint", "stratum", "arm", "tf"], kind="stable", ignore_index=True)


# ---------------------------------------------------------------- decision rule
RULE = {
    "H1": "the motif prior raises recovery of knockdown responders, and the gain needs TF-specific motif identity",
    "H1a_only": "the motif lift on knockdown truth is consistent with a TF-generic gene propensity "
                "(no TF-specific motif claim)",
    "H1a_refuted": "the prior lowers recovery of causal responders; C6 is limited to agreement with the literature",
    "H2": "its targets predict knockdown responses beyond gene propensity",
    "H2_none": "no evidence, in either screen, that networks inferred from control cells predict knockdown "
               "responses beyond gene propensity (reported as inconclusive, never as no difference)",
}


def _state(v: str) -> str:
    """A verdict, from one screen ('supported') or from ``combine_sets``, as one of: both_supported, both_refuted,
    one_screen, inconclusive, missing."""
    if v == "supported (replicated)":
        return "both_supported"
    if v == "refuted (replicated)":
        return "both_refuted"
    if v == "missing":
        return "missing"
    if " in design, " in v:       # a verdict from one screen only, including a hypothesis absent from the other
        return "one_screen"
    return "inconclusive"


def decision_lines(verdicts: dict[str, str], combined: bool) -> list[str]:
    """The record's decision rule applied per model family (H1a/H1b) and per model (H2).

    ``combined`` False: verdicts are from one screen, so a supported or refuted verdict only becomes a result when
    the other screen gives the same verdict; the line says what the rule would read if it does."""
    st = {h: (_state(v) if combined else {"supported": "both_supported", "refuted": "both_refuted",
                                          "missing": "missing"}.get(v, "inconclusive"))
          for h, v in verdicts.items()}
    cond = "" if combined else "If the other screen gives the same verdicts: "
    lines = []
    for h in verdicts:
        if not h.startswith("H1a "):
            continue
        m = h[4:].split(" > ")[0]
        h1b = next((x for x in verdicts if x.startswith(f"H1b {m} > ")), None)
        a, b = st[h], st.get(h1b, "missing")
        if a == "both_supported" and b == "both_supported":
            text = f"{m}: {RULE['H1']}."
        elif a == "both_supported":
            text = f"{m}: {RULE['H1a_only']}. H1b is {verdicts.get(h1b, 'missing')}."
        elif a == "both_refuted":
            text = f"{m}: {RULE['H1a_refuted']}."
        elif a == "one_screen":
            text = f"{m}: H1a not replicated, no claim ({verdicts[h]})."
        elif a == "missing":
            text = f"{m}: H1a missing, no claim."
        else:
            text = f"{m}: H1a inconclusive; reported as inconclusive, never as no difference."
        lines.append(("- " + cond + text) if a in ("both_supported", "both_refuted") else "- " + text)
    h2 = {h.split()[1]: st[h] for h in verdicts if h.startswith("H2 ")}
    for m, a in h2.items():
        full = next(x for x in verdicts if x.startswith(f"H2 {m} "))
        if a == "both_supported":
            lines.append(f"- {cond}{m}: {RULE['H2']}.")
        elif a == "both_refuted":
            lines.append(f"- {cond}{m}: its targets rank responders worse than its own propensity (not covered by the "
                         "rule; reported).")
        elif a == "one_screen":
            lines.append(f"- {m}: H2 not replicated, no claim ({verdicts[full]}).")
        elif a == "missing":
            lines.append(f"- {m}: H2 missing, no claim.")
        else:
            lines.append(f"- {m}: H2 inconclusive; reported as inconclusive, never as no difference.")
    if h2 and all(a == "inconclusive" for a in h2.values()):
        lines.append(f"- H2 is inconclusive for every model: {RULE['H2_none']}.")
    return lines


# ---------------------------------------------------------------- report
def _p(x) -> str:
    return "nan" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.3g}"


def _table(header, rows) -> list[str]:
    return ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)] + ["| " + " | ".join(map(str, r)) + " |" for r in rows]


def _family_rows(res):
    rows = []
    for r in res["family"]:
        ci = r["ci95"]
        ci = "unbounded" if not ci or ci[0] is None else f"[{ci[0]:+.3f}, {ci[1]:+.3f}]"
        rows.append((r["hypothesis"], f"{r['mean_contrast']:+.4f}", ci, _p(r["p_t_greater_holm"]),
                     _p(r["p_perm_greater_holm"]), _p(r["p_t_less_holm"]), _p(r["p_perm_less_holm"]), r["verdict"]))
    return rows


FAMILY_HEADER = ("hypothesis", "mean contrast", "95% interval", "Holm p_t", "Holm p_perm", "Holm p_t (reverse)",
                 "Holm p_perm (reverse)", "verdict")


def random_auprc_ratio(prim) -> float:
    n = prim["n_genes"] - 1
    ratios = []
    for r in prim["per_unit"].n_responders:
        labels = np.r_[np.ones(int(r), dtype=int), np.zeros(n - int(r), dtype=int)]
        ratios.append(expected_average_precision(np.zeros(n), labels) / (r / n))
    return float(np.mean(ratios))


def build_report(meta, truth_meta, results, det, perms, seed, strata_note, random_ratio, missing=None) -> str:
    prim = results["endpoints"]["auprc_ratio"]
    n_genes = prim["n_genes"]
    pu = prim["per_unit"]
    prevalence = float((pu.n_responders / (n_genes - 1)).mean())
    L = [f"# {EXPERIMENT} knockdown scores", "",
         f"- screen: {truth_meta.get('h5ad', 'unknown')}; truth tables verified against truth.json: "
         f"{'yes' if truth_meta.get('sha256') else 'NO truth.json hashes found'}",
         f"- inference: code {meta['code_commit']}; {meta['cells_in']} control cells, {meta['cells_after_qc']} after QC, "
         f"{meta['cells_after_cap']} after the cap ({meta['max_cells']}, seed {meta['cap_seed']}, "
         f"{'applied' if meta['cap_applied'] else 'not needed'}); {meta['n_genes']} genes; "
         f"{meta['n_regulators']} regulators (CollecTRI sources in the genes)",
         f"- units analysed: {prim['n_units']}; excluded: {len(prim['excluded'])}"
         + (" (" + "; ".join(f"{t}: {r}" for t, r in list(prim['excluded'].items())[:5]) + (" ..." if len(prim['excluded']) > 5 else "") + ")"
            if prim["excluded"] else ""),
         f"- permutations: {perms} joint TF-label draws, seed {seed}; Holm over {prim['n_directed_hypotheses']} directed "
         "hypotheses, separately for the t and for the permutation p-values, alpha 0.05", ""]
    if meta["skipped"]:
        L += ["## Arms that did not run", ""] + [f"- {a}: {w}" for a, w in meta["skipped"].items()] + [""]
    if missing:
        L += ["## Registered hypotheses that could not be computed", ""] + [f"- {h}: {w}" for h, w in missing.items()] + [""]
    L += [f"## Arms (mean over units; random expectation: AUPRC ratio {random_ratio:.3f} "
          f"(expected AP of a random order, mean over units), hit rate = prevalence {prevalence:.4f})", ""]
    arms = list(prim["arm_means"])
    sec = {a: meta["arms"][a]["seconds_including_base"] for a in meta["arms"]}
    L += _table(("arm", "AUPRC ratio", "top-50 hit rate", "top-200 hit rate", "seconds"),
                [(a, f"{prim['arm_means'][a]:.3f}", f"{results['endpoints']['top50']['arm_means'][a]:.4f}",
                  f"{results['endpoints']['top200']['arm_means'][a]:.4f}", sec.get(a, "")) for a in arms]) + [""]
    L += ["## Primary family (AUPRC ratio)", ""] + _table(FAMILY_HEADER, _family_rows(prim)) + [""]
    verdicts = {r["hypothesis"]: r["verdict"] for r in prim["family"]}
    L += ["## Decision rule, this screen", "",
          "A result needs the same verdict in the design and the replication screen (`combine`).", ""]
    L += decision_lines(verdicts, combined=False) + [""]
    L += ["## Secondary (outside Holm; unadjusted p)", ""]
    L += _table(("contrast", "mean contrast", "p_t", "p_perm"),
                [(r["hypothesis"], f"{r['mean_contrast']:+.4f}", _p(r["p_t_greater"]), _p(r["p_perm_greater"]))
                 for r in prim["secondary"]]) + [""]
    for key, label in (("top50", "top-50 hit rate"), ("top200", "top-200 hit rate")):
        L += [f"### {label} family", ""] + _table(FAMILY_HEADER, _family_rows(results["endpoints"][key])) + [""]
    for name, res in results["strata"].items():
        L += [f"### Stratum {name} (AUPRC ratio)", ""]
        if res.get("per_unit") is None:
            L += [f"skipped: {res['skipped']}", ""]
        else:
            L += [f"{res['n_units']} units.", ""] + _table(FAMILY_HEADER, _family_rows(res)) + [""]
    if strata_note:
        L += [strata_note, ""]
    rhos = [(a, d["spearman"], d["n"]) for a, d in det.items()]
    L += ["### TF detection in the inference cells vs the unit's score (Spearman, descriptive)", ""]
    L += _table(("arm", "n", "rho"), [(a, n, "nan" if r is None else f"{r:+.3f}") for a, r, n in rhos]) + [""]
    L += ["Not measured here: effect of the preprocessing choices, of the seeds, and of the responder rule on the verdicts.", ""]
    return "\n".join(L)


def missing_hypotheses(meta: dict, family: list) -> dict[str, str]:
    """Hypotheses of the scored family (``evaluate()["family"]``) with verdict ``missing``, each with the reason."""
    return {r["hypothesis"]: "; ".join(f"{a}: skipped by infer ({meta['skipped'][a]})" if a in meta["skipped"]
                                       else f"{a}: not requested" for a in r["unavailable"])
            for r in family if r["verdict"] == "missing"}


def check_same_screen(infer_dir: Path, meta: dict, truth: Path, truth_meta: dict) -> None:
    """Refuse a swapped screen, and any control cell used for both inference and truth."""
    if meta.get("source_h5ad") != truth_meta.get("h5ad"):
        raise SystemExit(f"inference was run on {meta.get('source_h5ad')!r} but the truth is from {truth_meta.get('h5ad')!r}")
    files = [infer_dir / "cells.txt", truth / "truth_half.txt"]
    gone = [str(f) for f in files if not f.exists()]
    if gone:
        raise SystemExit(f"cannot check that no control cell is used twice, missing: {', '.join(gone)}")
    shared = set(files[0].read_text(encoding="utf-8").split()) & set(files[1].read_text(encoding="utf-8").split())
    if shared:
        raise SystemExit(f"{len(shared)} control cells are in both the inference cells and the truth half, e.g. {sorted(shared)[:3]}")


def cmd_score(args) -> int:
    infer_dir, out = Path(args.infer), Path(args.out)
    meta, edges, genes = load_arms(infer_dir)
    usable, resp, truth_meta = load_truth(args.truth)
    check_same_screen(infer_dir, meta, Path(args.truth), truth_meta)
    regs = pd.read_csv(infer_dir / "regulators.tsv", sep="\t")
    model_arms = list(edges)
    # Holm always runs over the registered family (plus any extra model requested), never over what finished
    registered = list(dict.fromkeys(DEFAULT_MODELS.split(",") + meta["models_requested"]))
    family, secondary = kd.default_family(registered)[0], kd.default_family(model_arms)[1]
    missing_arms = set(registered) - set(model_arms)
    if args.collectri:
        edges = {**edges, REFERENCE_ARM: reference_edges(args.collectri)}
        secondary = secondary + [(f"S {REFERENCE_ARM} > own propensity", {REFERENCE_ARM: 1, kd.PROP + REFERENCE_ARM: -1})]
    strata, note = {}, ""
    usable_tfs = set(usable.tf.astype(str).str.upper())
    if args.essential_from:
        essential = perturbed_tfs(args.essential_from)
        strata["essential"] = usable_tfs & essential
        strata["non_essential"] = usable_tfs - essential
        note = ("Essential = the TF is perturbed in the design screen (check P5); the strata run the same family on "
                "their own units, with the propensity controls built from all regulators of each network.")
    if args.both_with:
        strata["usable_in_both_screens"] = usable_tfs & set(load_truth(args.both_with)[0].tf.astype(str).str.upper())
    if args.total_shift:
        shift = pd.read_csv(args.total_shift, sep="\t", keep_default_na=False)
        strata["total_within_0.9_1.1"] = usable_tfs & set(shift.tf[shift.within_0_9_1_1.astype(str) == "True"].astype(str).str.upper())
        note = (note + "\n\n" if note else "") + ("Stratum total_within_0.9_1.1 is the sensitivity analysis declared in check P9 "
                "(the family without the TFs whose relative core-adjusted total is outside [0.9, 1.1]); descriptive, "
                "it never changes a primary verdict.")
    results = run_endpoints(edges, resp, genes, args.perms, args.seed, family, secondary, strata, missing_arms)
    prim = results["endpoints"]["auprc_ratio"]
    missing = missing_hypotheses(meta, prim["family"])
    random_ratio = random_auprc_ratio(prim)
    det = detection(prim["per_unit"], regs)

    out.mkdir(parents=True, exist_ok=True)
    long_scores(results).to_csv(out / "scores.tsv", sep="\t", index=False, float_format="%.10g")
    dump = {"experiment": EXPERIMENT, "perms": args.perms, "seed": args.seed,
            "random_auprc_ratio": random_ratio,
            "inputs": {"truth": truth_meta.get("sha256"), "tfs_sha256": meta["tfs_sha256"],
                       "edges_sha256": {a: _sha_content(infer_dir / v["file"]) for a, v in meta["arms"].items()},
                       "code_commit": meta["code_commit"]},
            "source_h5ad": meta.get("source_h5ad"), "missing_hypotheses": missing,
            "skipped_arms": meta["skipped"], "detection": det, "strata_members": {k: len(v) for k, v in strata.items()}}
    for key, res in results["endpoints"].items():
        dump["primary" if key == "auprc_ratio" else key] = {k: v for k, v in res.items() if k != "per_unit"}
    dump["strata"] = {k: {x: y for x, y in v.items() if x != "per_unit"} for k, v in results["strata"].items()}
    dump["verdicts"] = {r["hypothesis"]: r["verdict"] for r in dump["primary"]["family"]}
    (out / "family.json").write_text(json.dumps(dump, indent=1, allow_nan=True) + "\n", encoding="utf-8")
    (out / "REPORT.md").write_text(build_report(meta, truth_meta, results, det, args.perms, args.seed, note,
                                                    random_ratio, missing),
                                   encoding="utf-8")
    for r in dump["primary"]["family"]:
        print(f"{r['hypothesis']}: {r['mean_contrast']:+.4f}, p_t_holm {r['p_t_greater_holm']:.3g}, "
              f"p_perm_holm {r['p_perm_greater_holm']:.3g} -> {r['verdict']}")
    return 0


def cmd_combine(args) -> int:
    design = json.loads((Path(args.design) / "family.json").read_text(encoding="utf-8"))
    repl = json.loads((Path(args.replication) / "family.json").read_text(encoding="utf-8"))
    if design.get("source_h5ad") == repl.get("source_h5ad"):
        raise SystemExit(f"design and replication come from the same screen ({design.get('source_h5ad')!r}); "
                         "a claim needs two different screens")
    final = kd.combine_sets(design["primary"], repl["primary"])
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "combined.json").write_text(json.dumps(final, indent=1) + "\n", encoding="utf-8")
    lines = [f"# {EXPERIMENT} combined verdicts", "",
             "A hypothesis is a result only if it has the same verdict in the design and the replication screen.", ""]
    lines += _table(("hypothesis", "design", "replication", "final"),
                    [(h, design["verdicts"].get(h, "missing"), repl["verdicts"].get(h, "missing"), v)
                     for h, v in final.items()]) + ["", "## What the paper may say", ""]
    lines += decision_lines(final, combined=True) + [""]
    (out / "COMBINED.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(f"{h}: {v}" for h, v in final.items()))
    return 0


# ---------------------------------------------------------------- run identity (resume guard)
SCREENS = ("rpe1", "k562")


def run_identity(snapshot, models, perms: int, norm: str) -> dict:
    """What makes two jobs the same run: code, data (the D-020 snapshot pins it), models, permutations, normalisation."""
    return {"experiment": EXPERIMENT, "code_commit": ms.code_commit(), "snapshot_sha256": _sha(Path(snapshot)),
            "models": list(dict.fromkeys(m.strip() for m in models if m.strip())), "perms": int(perms), "norm": norm}


def check_run(run, ident: dict) -> None:
    """First call writes ``run.json`` in the run folder; later calls refuse (SystemExit) a different identity, a folder
    with finished steps and no ``run.json``, and any finished infer or score folder not built by this run."""
    run = Path(run)
    run.mkdir(parents=True, exist_ok=True)
    sentinel = run / "run.json"
    steps = [str(d.relative_to(run)) for d in [*run.glob("*/*"), run / "combined"]
             if d.is_dir() and d.name in ("truth", "infer", "score", "combined")]
    if not sentinel.exists():
        if steps:
            raise SystemExit(f"{run} holds finished steps ({', '.join(steps)}) but no run.json: not resumed")
        sentinel.write_text(json.dumps(ident, indent=1) + "\n", encoding="utf-8")
        return
    have = json.loads(sentinel.read_text(encoding="utf-8"))
    bad = [f"{k}: run has {have.get(k)!r}, this job has {v!r}" for k, v in ident.items() if have.get(k) != v]
    for s in SCREENS:
        d = run / s
        arms = {}
        if (d / "infer").exists():
            m = json.loads((d / "infer" / "infer.json").read_text(encoding="utf-8"))
            arms = {a: _sha_content(d / "infer" / v["file"]) for a, v in m["arms"].items()}
            if (m["code_commit"], m["models_requested"]) != (ident["code_commit"], ident["models"]):
                bad.append(f"{s}/infer was built by commit {m['code_commit']} with models {m['models_requested']}")
        if (d / "score").exists():
            f = json.loads((d / "score" / "family.json").read_text(encoding="utf-8"))
            truth = json.loads((d / "truth" / "truth.json").read_text(encoding="utf-8")).get("sha256")
            want = {"perms": ident["perms"], "inputs": {**f["inputs"], "code_commit": ident["code_commit"],
                                                        "truth": truth, "edges_sha256": arms}}
            if {"perms": f["perms"], "inputs": f["inputs"]} != want:
                bad.append(f"{s}/score was not built by this commit, permutation count, truth and infer folder")
    if bad:
        raise SystemExit(f"{run} is not this run (run.json), refusing to resume:\n  " + "\n  ".join(bad))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    i = sub.add_parser("infer")
    i.add_argument("--controls", required=True, help="inference-half control cells (knockdown_bench controls)")
    i.add_argument("--tfs", required=True, help="CollecTRI tsv with a 'source' column")
    i.add_argument("--experiment", required=True, help="experiment ID written into infer.json, e.g. EXP-20261016-LM-1")
    i.add_argument("--models", default=DEFAULT_MODELS, help="comma-separated arms (default: the record's nine)")
    i.add_argument("--out", required=True)
    i.add_argument("--attempts", help="json file that keeps failures per arm and commit across resubmissions (exclusion (c))")
    i.add_argument("--seed", type=int, default=CAP_SEED, help="cell cap seed (record: 0)")
    i.add_argument("--max-cells", type=int, default=MAX_SCORE_CELLS, help="cell cap after QC (record: 3000)")
    i.add_argument("--device", default="cpu", help="torch device for RegDiffusion (cpu or cuda); recorded in infer.json")
    s = sub.add_parser("score")
    s.add_argument("--infer", required=True)
    s.add_argument("--truth", required=True, help="knockdown_bench truth directory")
    s.add_argument("--out", required=True)
    s.add_argument("--perms", type=int, default=10000)
    s.add_argument("--seed", type=int, default=0, help="permutation seed (record: 0)")
    s.add_argument("--collectri", help="CollecTRI tsv: adds the descriptive reference row (outside the family)")
    s.add_argument("--essential-from", help="truth dir of the design screen (RPE1): TFs perturbed there are essential (P5)")
    s.add_argument("--both-with", help="truth dir of the other screen: adds the stratum of TFs usable in both")
    s.add_argument("--total-shift", help="knockdown_bench total-shift tsv: adds the P9 sensitivity stratum (descriptive)")
    c = sub.add_parser("combine")
    c.add_argument("--design", required=True)
    c.add_argument("--replication", required=True)
    c.add_argument("--out", required=True)
    r = sub.add_parser("run-identity", help="write or check the run's run.json (resume guard)")
    r.add_argument("--run", required=True)
    r.add_argument("--snapshot", required=True)
    r.add_argument("--models", default=DEFAULT_MODELS)
    r.add_argument("--perms", type=int, required=True)
    r.add_argument("--norm", required=True)
    a = p.parse_args(argv)
    if a.cmd == "run-identity":
        check_run(a.run, run_identity(a.snapshot, a.models.split(","), a.perms, a.norm))
        return 0
    if a.cmd == "infer":
        meta = infer(a.controls, a.tfs, a.models.split(","), a.out, a.seed, a.max_cells, a.device, a.attempts,
                     experiment=a.experiment)
        print(f"{len(meta['arms'])} arms, {len(meta['skipped'])} skipped -> {a.out}")
        if meta["status"] != "complete":
            print(f"inference incomplete: arms failed (attempt counts {meta['attempt']}); resubmit the same commit once "
                  "(record exclusion (c))", file=sys.stderr)
            return 3
        return 0
    return cmd_score(a) if a.cmd == "score" else cmd_combine(a)


if __name__ == "__main__":
    sys.exit(main())
