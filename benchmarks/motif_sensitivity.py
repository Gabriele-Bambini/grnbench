"""Motif-prior lift decomposition, stage 1: Pearson recomputed from kept artifacts (exploratory).

EXP-20261003-GB-1 (issue #52). Adapted from ``audit/2026-09-27/proposal/motif_sensitivity.py``,
the draft of the 2026-09-27 audit (``audit/2026-09-27/AUDIT.md`` sections 3 and 6, design in
``audit/2026-09-27/proposal/PLAN.md``; branch ``lm/audit-2026-09-27``, PR #45), to the APIs on
``main`` and to stage 1 only. Changes from the draft: calibration calls
``real_bench.score_models`` itself on the kept inputs instead of a copy of its loop; the refit
paths (GRNBoost2, RegDiffusion) refuse instead of running; the Holm family is per set
(1 pair x 3 tests x 2 directions = 6); the verdict rule of PLAN.md section 5 is computed here.

Question: how much of the pooled AUPRC-ratio lift of ``pearson_motif`` over ``pearson_signed``
(claims C4, C6) needs the motif prior's target identity, and how much comes from the two design
effects the audit found (TF-level reweighting; cisTarget database membership)?

Reads the ``--keep-artifacts`` tree of a finished ``real_bench.py`` run and never writes into it.
Per dataset x arm (auto, baseline) x context:

  calibration  ``real_bench.score_models`` re-scores every kept ``processed.h5ad`` with the
               closed-form models the run scored (pearson_signed, pearson_motif, motif_prior).
               Nothing is refit (tier B of PLAN.md). Gate 1: each model's per-dataset AUPRC ratio
               equals the committed row (|diff| <= 1e-12) and its re-written top edges equal the
               kept ``<model>.tsv.gz`` (same pairs in the same order, weights within 1e-12).
               Gate 2 (report): ``real_bench.models_section`` on the committed rows, with the
               recomputed values substituted, prints the committed declared-comparison lines.
  prior        the registered motif matrix M (``cistarget.motif_ranks``), the genes present in
               the databases, K (annotated motifs per TF) and W (the worst rank, the tie value).
  variants     M and its controls, each rank-producted with the base edges by a line-for-line
               copy of ``cistarget.motif_rank_product`` (checked identical for the unchanged M):
               perm_db (primary null: each TF's in-database values re-assigned to its
               in-database genes in a keyed random order, seeds 0-99; the key depends on seed, TF
               and gene only, so one seed is one random database shared by every dataset),
               perm_all, tf_mean, tf_kexp, tf_presence, db_only, worst_tie, tfmean_db, gene_mean
               and gene_tf (a TF-generic gene propensity: the leave-one-out mean of M over the
               other TFs with a motif). Built-in checks (constant-row and tie controls cannot move
               their metrics) are asserted on every context.
  metrics      pooled (the registered AUPRC ratio), macro (mean within-TF AUPRC ratio over TFs
               with >= 5 evaluation positives), pooled_T / macro_T (targets in the databases and
               labelable in the full CollecTRI), pooled_R / macro_R (same, TFs with >= 1 motif).
  tests        per set, auto arm: T1 pooled(real) > pooled(perm_db), T2 macro_R(real) >
               macro_R(base), T3 macro_R(real) > macro_R(gene_tf); paired t over datasets
               (compare_units), Holm over those 6 directed hypotheses (holm_adjust). Baseline arm:
               the same, secondary.

    python benchmarks/motif_sensitivity.py plan   --run EXP=RUN_DIR [--run ...] --out OUT/sens_tasks.json
    python benchmarks/motif_sensitivity.py run    --tasks OUT/sens_tasks.json --refs REFS --out OUT --jobs 8
    python benchmarks/motif_sensitivity.py report --out OUT

``run`` accepts what ``hpc/engine_bench.sbatch`` passes (``--refs --index --chunk --jobs
--workroot --max-attempts``) and needs ``GRNBENCH_CISTARGET_DIR``. It checks first (before any
unit runs, any row or ``RUN_INPUTS.json`` is written) that ``--refs`` holds
``collectri_<organism>.tsv`` and that the cisTarget folder holds a rankings database and a motif
table for every organism of the runnable units (``preflight_problems``), and exits 2 otherwise.
Every row records the code commit that wrote it (``code_commit``). Finished datasets are skipped on
a resubmission with the same commit and the same inputs; an OUT holding rows of another commit, or
a ``RUN_INPUTS.json`` with other input hashes, is refused (``resume_problems``,
``run_inputs_payload``), and a skipped ``check_failed`` unit exits 2. Stage 2 (refitting GRNBoost2
or RegDiffusion) is refused: see STAGE2_REFUSAL.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shutil
import sys
import tempfile
import traceback
from concurrent.futures.process import BrokenProcessPool
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import real_bench  # noqa: E402
from engine_bench import _atomic_json, pool_results  # noqa: E402

EXPERIMENT = "EXP-20261003-GB-1"
REPO_ROOT = Path(__file__).resolve().parents[1]
REGISTERED_SETS = {"EXP-20261002-LM-1": "primary set (C6, 20 never-used datasets)",
                   "EXP-20260925-LM-2": "secondary set (C4)"}
STAGE1_PAIRS = (("pearson_motif", "pearson_signed"),)
RECOMPUTABLE = ("pearson_signed", "pearson_motif", "motif_prior")   # closed form: no fit, no seed
REFIT_BASES = ("grnboost2", "regdiffusion")
STAGE2_REFUSAL = (
    "{model}: stage 2 (refitting GRNBoost2 or RegDiffusion) is out of scope of EXP-20261003-GB-1. "
    "The audit's GRNBoost2 calibration gate failed on a mock (audit/2026-09-27/AUDIT.md section 6, "
    "proposal/PLAN.md section 0) and has not been fixed and re-validated. Fix it, validate it "
    "against the kept C5/C6 artifacts, and register stage 2 as a dated amendment or a new record "
    "before running it. Stage 1 runs pearson_motif vs pearson_signed only.")
ARMS = ("auto", "baseline")
PRIMARY_ARM = "auto"
N_PERM = 100                    # permutation seeds 0..99
MIN_POS_TF = 5                  # evaluation positives a TF needs to enter a macro metric
TOL = 1e-12                     # calibration tolerance (absolute)
FIXED = ("real", "tf_mean", "tf_kexp", "tf_presence", "db_only", "worst_tie", "tfmean_db",
         "gene_mean", "gene_tf", "gene_tf_np")     # gene_tf_np: descriptive (record 2026-10-06 note)
PERMS = ("perm_all", "perm_db")
VARIANTS = FIXED + PERMS
METRICS = ("pooled", "macro", "pooled_T", "macro_T", "pooled_R", "macro_R")
SEED_KEYS = ("pooled", "macro_R")     # per-seed dataset means kept in the row (T1 robustness)
N_BOOT = 2000
ALPHA = 0.05
VERDICTS = {
    "target-level": (
        "supports target-level evidence",
        "T1, T2 and T3 are supported; the gap to the random assignment (D) is at least half of the "
        "registered lift L; the design-effect control (tfmean_db) and the TF-generic gene "
        "propensity (gene_tf) each reach less than half of L; T1's permutation rank p < 0.05. "
        "The caveat on C4/C6 can be narrowed for Pearson on this set (exploratory; no claim "
        "status changes)."),
    "mixed": (
        "mixed",
        "T1 and/or T2 is supported by the paired t test over datasets, but the strongest verdict is "
        "not reached (conditions named on this line). These t tests are conditional on the one "
        "cisTarget database and reference shared by every dataset; on their own they are not "
        "evidence that target identity matters (only T1's permutation rank accounts for the "
        "database realization). The caveat on C4/C6 stays and names the components with their "
        "sizes."),
    "design effects only": (
        "consistent with design effects only",
        "One of three routes (record rules 2 and 5): (a) T1 reversed: a random target assignment "
        "beats the real prior; (b) T2 reversed: the real prior lowers within-TF ranking where "
        "neither design effect can act; (c) neither T1 nor T2 is supported and the design-effect "
        "control reaches at least half of the lift. "
        # the consequence sentence, unchanged
        "C4/C6 may say that the motif prior raises the pooled AUPRC ratio, not that sequence "
        "evidence at the target does; a confirmatory test with macro_R as primary metric needs a "
        "D-xxx."),
    "inconclusive": (
        "inconclusive",
        "Inconclusive is not equivalence: the caveat on C4/C6 stays as worded."),
    "not interpretable": (
        "not interpretable",
        "A code check failed, or gate 2 failed although every dataset passed gate 1. No result of "
        "this set is read until the cause is found, fixed and recorded as a dated deviation."),
}
T1_NOT_ROBUST = ("T1 not robust to the realization of the one shared database (permutation rank "
                 "p >= 0.05 or not assessable)")
# printed under every table of declared tests (record section 1, robustness condition on T1)
SHARED_DATABASE = ("All datasets of this set share one cisTarget database and one reference: these "
                   "p-values measure consistency over datasets given that database, not inference over "
                   "database realizations. Only the T1 permutation rank below accounts for it; T2 and T3 "
                   "have no such check in this record.")


class Stage2Refused(RuntimeError):
    """A stage-2 model (refit of GRNBoost2 / RegDiffusion) was requested in stage 1."""


class CheckFailed(RuntimeError):
    """A code check failed (self-consistency, fast scorer, internal consistency, built-in checks).
    Not a transient error: the unit is ``check_failed``, it is not resubmitted, and the set it
    belongs to is not interpretable until the code is fixed."""


class ScoreModelsError(RuntimeError):
    """``score_models`` recorded an exception for a model (``m_<model>_<arm>_error``) instead of
    raising: a runtime error (exclusion rule c; I/O, memory), resubmitted with the same commit and
    excluded if it persists. Never a failed code check: those are for finite values that disagree."""


def check_stage1_pairs(pairs) -> None:
    for combo, base in pairs:
        for model in (combo, base):
            if any(model.startswith(b) for b in REFIT_BASES):
                raise Stage2Refused(STAGE2_REFUSAL.format(model=model))
        if (combo, base) not in STAGE1_PAIRS:
            raise ValueError(f"{combo} vs {base}: not a stage-1 pair ({STAGE1_PAIRS})")


# ----------------------------------------------------------------------------- inputs
def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _portable(path: Path) -> str:
    """Repo-relative inside the repo, ``~/...`` below HOME: never a personal absolute path."""
    p = Path(path).expanduser().resolve()
    try:
        return p.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return real_bench._portable_user_path(p)


def _resolve(path: str) -> Path:
    p = Path(path).expanduser()
    return p if p.is_absolute() else REPO_ROOT / p


def ctx_dirs(task_dir: Path, arm: str) -> list[Path]:
    """Kept contexts of one arm, in score_models order (``ctx<i>``, i = record index)."""
    d = Path(task_dir) / arm
    return sorted((p for p in d.glob("ctx*") if (p / "processed.h5ad").exists()),
                  key=lambda p: int(p.name[3:]))


def read_pairs(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def read_kept(path: Path) -> pd.DataFrame | None:
    """A kept edge file with exact floats (``round_trip``; the default parser is off by ~1e-16)."""
    if not Path(path).exists():
        return None
    return pd.read_csv(path, sep="\t", dtype={"source": str, "target": str},
                       keep_default_na=False, float_precision="round_trip")


def compare_edge_files(kept: Path, recomputed: Path) -> dict:
    a, b = read_kept(kept), read_kept(recomputed)
    if a is None or b is None:
        return {"match": False, "why": "kept file missing" if a is None else "not recomputed"}
    same = (len(a) == len(b) and np.array_equal(a["source"].to_numpy(), b["source"].to_numpy())
            and np.array_equal(a["target"].to_numpy(), b["target"].to_numpy()))
    diff = float(np.max(np.abs(a["weight"].to_numpy() - b["weight"].to_numpy()))) if same and len(a) else np.nan
    return {"match": bool(same and diff <= TOL), "rows": int(len(a)), "same_pairs_in_order": bool(same),
            "max_abs_weight_diff": diff}


def input_fingerprints(refs: str | None) -> dict:
    """sha256 of the reference networks and of the cisTarget files this run reads."""
    from grnbench.infer.cistarget import cistarget_dir
    root = cistarget_dir()
    cis = {}
    if root is not None:
        for p in sorted(root.glob("*.genes_vs_motifs.rankings.feather")) + sorted(root.glob("*.tbl")):
            cis[p.name] = _sha256(p)
    return {"references": real_bench.reference_fingerprints(refs), "cistarget": cis}


def preflight_problems(units: list[dict], refs: str) -> list[str]:
    """Why ``run`` cannot start (empty: it can), checked before any unit runs and before any row or
    ``RUN_INPUTS.json`` is written. For every organism of the runnable units (``truth.organism`` of
    their committed rows; units the plan excluded need nothing): ``--refs`` must hold
    ``collectri_<organism>.tsv``, and GRNBENCH_CISTARGET_DIR must be a folder with at least one
    rankings database and a motif table for that organism. A proxy by organism: scoring picks the
    species with ``guess_species(genes)`` (registered behaviour), which is why the data root must
    still pass ``data_manifest.py --check ... --require-all`` (record section 2)."""
    from grnbench.infer.cistarget import ENV_DIR, SPECIES_PREFIX, cistarget_dir, resource_files
    orgs = set()
    for u in units:
        if any(f.startswith("MISSING") for f in u.get("plan_flags", [])):
            continue
        try:
            row = json.loads((_resolve(u["committed"]) / "rows" / f"{u['task']}.json")
                             .read_text(encoding="utf-8"))
            orgs.add(json.loads(row["truth"]).get("organism"))
        except (OSError, ValueError, KeyError, TypeError, AttributeError):
            continue                    # unreadable committed row: the unit itself records the error
    problems = [f"--refs {_portable(Path(refs))}: no collectri_{org}.tsv"
                for org in sorted(orgs, key=str) if not (Path(refs) / f"collectri_{org}.tsv").is_file()]
    root = cistarget_dir()
    if root is None:
        problems.append(f"{ENV_DIR} is unset or not a folder")
    else:
        for org in sorted(o for o in orgs if o in SPECIES_PREFIX):
            dbs, ann = resource_files(org, root)
            if not dbs or ann is None:
                problems.append(f"{ENV_DIR}={_portable(root)}: no cisTarget rankings database and motif "
                                f"table for {org}")
    return problems


def code_commit() -> str:
    """The commit of the code that runs: ``GIT_COMMIT`` next to the code when it exists (written by
    ``hpc.sh sync``; the cluster path, no git call), else ``git rev-parse HEAD`` plus ``-dirty`` for
    a tree with changes (as ``hpc.sh sync`` marks it), else ``unknown``."""
    f = REPO_ROOT / "GIT_COMMIT"
    if f.is_file():
        return f.read_text(encoding="utf-8").strip() or "unknown"
    import subprocess
    try:
        git = ["git", "-C", str(REPO_ROOT)]
        head = subprocess.run(git + ["rev-parse", "HEAD"], capture_output=True, text=True, check=True,
                              timeout=60).stdout.strip()
        dirty = subprocess.run(git + ["status", "--porcelain"], capture_output=True, text=True, check=True,
                               timeout=60).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return "unknown"
    return head + ("-dirty" if dirty else "") if head else "unknown"


RERUN_RULE = ("Rule (c) resubmits with the same commit (re-sync {old}, do not run hpc.sh sync between the job "
              "and its resubmission). After a code fix, rerun every unit at the fixed commit into a new, empty "
              "--out and record the dated deviation naming both commits")
CODE_COMMIT_NOTE = (" hpc/engine_bench.sbatch has already copied this job's GIT_COMMIT to OUT/CODE_COMMIT before "
                    "this refusal: restore it with 'echo {old} > OUT/CODE_COMMIT' (or resubmit at {old}, which copies "
                    "it again), otherwise the next report shows MIXED CODE COMMITS")
CHECK_FAILED_RULE = ("a check_failed unit is never resubmitted (record section 1, code checks): after the code "
                     "fix, rerun every unit at the fixed commit into a new, empty --out and record the dated "
                     "deviation naming both commits")


def resume_problems(outdir: Path, commit: str, n_perm: int) -> list[str]:
    """Why ``run`` must not write into an OUT that already holds rows (empty: it may). Every row
    records the commit that wrote it; one OUT never mixes commits (rule c resubmits with the same
    commit; a code fix reruns every unit into a new OUT). A row that is skipped on resubmission
    (completed, excluded, check_failed) and records another ``n_perm`` would mix permutation seed
    counts in one set: refused too, nothing overwritten. Unreadable rows are rerun, as before."""
    others: dict[str, list[str]] = {}
    seeds: dict[int, list[str]] = {}
    for p in sorted((Path(outdir) / "rows").glob("*.json")):
        try:
            row = json.loads(p.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        old = row.get("code_commit")
        if old != commit:
            others.setdefault(old or "none recorded", []).append(row.get("id", p.stem))
        if row.get("status") in ("completed", "excluded", "check_failed") and "n_perm" in row \
                and row["n_perm"] != n_perm:
            seeds.setdefault(row["n_perm"], []).append(row.get("id", p.stem))
    return [f"{_portable(outdir)} holds {len(ids)} row(s) from commit {old} (e.g. {ids[0]}), this job runs "
            f"{commit}. " + RERUN_RULE.format(old=old) + ("" if old == "none recorded" else CODE_COMMIT_NOTE.format(old=old))
            for old, ids in sorted(others.items())] + \
        [f"{ids[0]}" + (f" and {len(ids) - 1} more row(s)" if len(ids) > 1 else "") +
         f": existing row has n_perm={old}, this run asks {n_perm}; use a fresh --out or rerun them with one "
         f"--n-perm (registered: {N_PERM})" for old, ids in sorted(seeds.items())]


def run_inputs_payload(path: Path, fingerprints: dict, invocation: dict) -> tuple[dict | None, str]:
    """(RUN_INPUTS.json to write, '') or (None, why not). The first submission writes the
    fingerprints, its platform and ``invocations`` [{platform, code_commit}]; a resubmission keeps
    the file (the top-level platform stays the first one's) and appends its invocation, and refuses
    when the reference or cisTarget fingerprints differ (never a blind overwrite)."""
    if not Path(path).exists():
        return {**fingerprints, "platform": invocation["platform"], "invocations": [invocation]}, ""
    old = json.loads(Path(path).read_text(encoding="utf-8"))
    changed = []
    for kind in ("references", "cistarget"):
        a, b = old.get(kind) or {}, fingerprints.get(kind) or {}
        changed += [f"{kind}/{n}" for n in sorted(set(a) | set(b)) if a.get(n) != b.get(n)]
    if changed:
        return None, (f"inputs changed between submissions: {', '.join(changed)} differ from "
                      f"{_portable(path)} (missing, added or another sha256). A resubmission reads the "
                      "same inputs; other inputs go into a new, empty --out")
    old.setdefault("invocations", []).append(invocation)
    return old, ""


def platform_info() -> dict:
    """OS, CPU model and BLAS of the machine that ran gate 1 (exactness across platforms is not
    guaranteed). No host name, user name or path (AGENTS.md section 11)."""
    cpu = platform.processor() or ""
    try:
        with open("/proc/cpuinfo", encoding="utf-8") as fh:
            cpu = next((ln.split(":", 1)[1].strip() for ln in fh if ln.startswith("model name")), cpu)
    except OSError:
        pass
    blas = {}
    try:
        blas = {k: v for k, v in np.show_config(mode="dicts")["Build Dependencies"]["blas"].items()
                if k in ("name", "version", "openblas configuration")}
    except Exception:  # noqa: BLE001 - older numpy: the version below still pins the wheel
        pass
    return {"os": platform.platform(), "machine": platform.machine(), "cpu": cpu,
            "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
            "blas": blas, "omp_num_threads": os.environ.get("OMP_NUM_THREADS")}


# ----------------------------------------------------------------------------- prior + variants
def rank_product(E: pd.DataFrame, m: np.ndarray) -> pd.DataFrame:
    """cistarget.motif_rank_product with the motif ranks passed in (same operations, same order)."""
    out = E.copy()
    base_rank = out["weight"].rank(ascending=False, method="average").to_numpy()
    mot_rank = pd.Series(m).rank(ascending=True, method="average").to_numpy()
    out["weight"] = -(np.log(base_rank) + np.log(mot_rank)) / 2
    return out.sort_values("weight", ascending=False, kind="stable").reset_index(drop=True)


def rank_product_weights(base_rank: np.ndarray, m: np.ndarray) -> np.ndarray:
    """The weights of ``rank_product`` in the input edge order (no sort: the scorer re-ranks)."""
    mot_rank = pd.Series(m).rank(ascending=True, method="average").to_numpy()
    return -(np.log(base_rank) + np.log(mot_rank)) / 2


def motif_parts(tfs_e: list[str], genes: tuple[str, ...]):
    """M (TFs x genes, the registered matrix), gene-in-database mask, K motifs per TF, worst rank."""
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
        K += np.array([sum(m in present for m in tf_motifs.get(t.casefold(), ())) for t in tfs_e], dtype=int)
        worst = n_total if worst is None else min(worst, n_total)
    return M, in_db, K, worst


def motif_sets(tfs_e: list[str], genes: tuple[str, ...]) -> list[frozenset]:
    """Per TF, the annotated motifs present in the ranking databases: the motifs ``motif_ranks`` takes
    its best rank over (K > 0 exactly when the set is not empty). Reads only the ``motifs`` column."""
    import pyarrow.feather as pf
    from grnbench.infer.cistarget import guess_species, load_annotation, resource_files
    dbs, ann = resource_files(guess_species(genes))
    tf_motifs = load_annotation(str(ann))
    present: set = set()
    for db in dbs:
        present |= set(pf.read_table(str(db), columns=["motifs"]).column("motifs").to_pylist())
    return [frozenset(m for m in tf_motifs.get(t.casefold(), ()) if m in present) for t in tfs_e]


def shared_matrix(sets: list[frozenset]) -> np.ndarray:
    """TFs x TFs, True where two different TFs share at least one motif (diagonal False)."""
    n = len(sets)
    out = np.zeros((n, n), dtype=bool)
    by_motif: dict = {}
    for i, st in enumerate(sets):
        for m in st:
            by_motif.setdefault(m, []).append(i)
    for members in by_motif.values():
        if len(members) > 1:
            out[np.ix_(members, members)] = True
    np.fill_diagonal(out, False)
    return out


_U64 = np.uint64


def _mix(x):
    """splitmix64 finalizer (uint64 in, uint64 out, wrapping arithmetic)."""
    with np.errstate(over="ignore"):
        x = x + _U64(0x9E3779B97F4A7C15)
        x = (x ^ (x >> _U64(30))) * _U64(0xBF58476D1CE4E5B9)
        x = (x ^ (x >> _U64(27))) * _U64(0x94D049BB133111EB)
        return x ^ (x >> _U64(31))


def name_keys(names) -> np.ndarray:
    """A fixed 64-bit key per TF or gene name (case-insensitive), the same in every dataset."""
    return np.array([int.from_bytes(hashlib.blake2b(str(n).upper().encode(), digest_size=8).digest(), "little")
                     for n in names], dtype=np.uint64)


def keyed_order(tf_keys: np.ndarray, gene_keys: np.ndarray, seed: int) -> np.ndarray:
    """Per TF (rows), the genes (columns) sorted by u = hash(seed, TF, gene). u does not depend on
    the dataset, so two contexts order the genes they share the same way."""
    tk = _mix(tf_keys ^ _mix(_U64(seed)))
    return np.argsort(_mix(gene_keys[None, :] ^ tk[:, None]), axis=1, kind="stable")


def variant_matrix(name: str, M, in_db, K, worst, seed: int | None = None, keys=None,
                   shared=None) -> np.ndarray:
    """One control of the motif matrix M (TFs x genes). Permutations need ``seed`` and
    ``keys`` = (TF keys, gene keys) from ``name_keys``; ``gene_tf_np`` needs ``shared`` =
    ``shared_matrix`` of the TFs' motif sets."""
    W = float(worst)
    has = K > 0
    if name == "real":
        return M
    if name in PERMS:
        # each TF's own values (in-database genes for perm_db, all genes for perm_all) re-assigned
        # in keyed random order: the j-th best value goes to the gene with the j-th smallest u
        if seed is None or keys is None:
            raise ValueError(f"{name} needs a seed and the TF / gene keys")
        tf_keys, gene_keys = keys
        out = M.astype(float, copy=True)
        cols = np.flatnonzero(in_db) if name == "perm_db" else np.arange(M.shape[1])
        if len(cols):
            sub = np.empty((M.shape[0], len(cols)))
            np.put_along_axis(sub, keyed_order(tf_keys, gene_keys[cols], seed),
                              np.sort(M[:, cols].astype(float), axis=1, kind="stable"), axis=1)
            out[:, cols] = sub
        return out
    if name in ("gene_mean", "gene_tf"):
        # TF-generic gene propensity: for TF i, the mean of M over the OTHER TFs with a motif on
        # each in-database gene. gene_mean uses it as the row; gene_tf keeps TF i's own sorted
        # in-database values (TF-level distribution, mask and ties of the real prior) and assigns
        # them in that order. Everything but TF-specific target identity.
        out = np.full(M.shape, W)
        rows, cols = np.flatnonzero(has), np.flatnonzero(in_db)
        if not len(rows) or not len(cols):
            return out
        sub = M[np.ix_(rows, cols)].astype(float)
        tot = sub.sum(0)
        for k, i in enumerate(rows):
            loo = (tot - sub[k]) / (len(rows) - 1) if len(rows) > 1 else sub[k]   # one TF: its own row
            if name == "gene_mean":
                out[i, cols] = loo
            else:
                out[i, cols[np.argsort(loo, kind="stable")]] = np.sort(sub[k], kind="stable")
        return out
    if name == "gene_tf_np":
        # descriptive: gene_tf with the propensity also leaving out every TF that shares a motif with
        # TF i (their rows carry part of TF i's own target information). A TF left with no other TF
        # keeps its real row (as gene_tf does with a single TF with a motif).
        if shared is None:
            raise ValueError("gene_tf_np needs the shared-motif matrix")
        out = np.full(M.shape, W)
        rows, cols = np.flatnonzero(has), np.flatnonzero(in_db)
        if not len(rows) or not len(cols):
            return out
        sub = M[np.ix_(rows, cols)].astype(float)
        tot = sub.sum(0)
        share = np.asarray(shared)[np.ix_(rows, rows)]
        for k, i in enumerate(rows):
            n_keep = len(rows) - 1 - int(share[k].sum())
            if n_keep >= 1:
                loo = (tot - sub[k] - sub[share[k]].sum(0)) / n_keep
                out[i, cols[np.argsort(loo, kind="stable")]] = np.sort(sub[k], kind="stable")
            else:
                out[i, cols] = sub[k]
        return out
    tf_mean = np.array([M[i, in_db].mean() if has[i] and in_db.any() else W for i in range(len(K))])
    if name == "tf_mean":
        return np.repeat(tf_mean[:, None], M.shape[1], 1)
    if name == "tf_kexp":
        return np.repeat(np.where(has, (W + 1) / (K + 1), W)[:, None], M.shape[1], 1)
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
    """Fast exact metrics for many weightings of one fixed edge list on one universe.

    ``pooled`` equals ``evaluate_ranking(...)["auprc_ratio"]`` (checked against it for every fixed
    variant and the first permutation seed); the macro and restricted metrics use the same
    tie-aware average precision within each TF / mask.
    """

    def __init__(self, universe: pd.DataFrame, gi, masks: dict[str, np.ndarray], edges: pd.DataFrame):
        from grnbench.eval.universe import canonical_pairs
        cu = canonical_pairs(universe, name="universe")
        key = pd.Series(np.arange(len(cu)), index=pd.MultiIndex.from_frame(cu))
        e = edges[["source", "target"]].astype(str).apply(lambda c: c.str.upper())
        pos = key.reindex(pd.MultiIndex.from_frame(e)).to_numpy()
        if np.isnan(pos.astype(float)).any():
            raise ValueError("edges outside the universe")
        self.pos = pos.astype(int)
        self.y = gi.labels(cu).astype(int)
        src = cu["source"].to_numpy()
        self.tf_names = list(pd.unique(src))
        self.blocks = [np.flatnonzero(src == s) for s in self.tf_names]
        self.masks = masks            # {'T': DB & labelable targets (U'), 'R': same x TFs with a motif (U'')}

    def universe_weights(self, w_edges: np.ndarray) -> np.ndarray:
        w = np.zeros(len(self.y))
        w[self.pos] = pd.Series(w_edges).rank(method="dense").to_numpy()   # as align_ranking
        return w

    @staticmethod
    def _ap_ratio(y, w) -> float:
        from grnbench.eval.metrics import average_precision
        if len(y) == 0 or y.sum() == 0:
            return np.nan
        o = np.argsort(-w, kind="stable")
        return average_precision(y[o], w[o]) / (y.sum() / len(y))

    def macro_blocks(self, mask=None) -> list[np.ndarray]:
        """Each TF's row of U (of the mask when given), kept when it holds at least MIN_POS_TF
        evaluation positives: the one definition of the TFs entering a macro metric."""
        rows = self.blocks if mask is None else [b[mask[b]] for b in self.blocks]
        return [b for b in rows if len(b) and self.y[b].sum() >= MIN_POS_TF]

    def macro(self, w, mask=None) -> float:
        vals = [self._ap_ratio(self.y[b], w[b]) for b in self.macro_blocks(mask)]
        return float(np.mean(vals)) if vals else np.nan

    def score(self, w_edges: np.ndarray) -> dict:
        w = self.universe_weights(w_edges)
        out = {"pooled": self._ap_ratio(self.y, w), "macro": self.macro(w)}
        for tag, m in self.masks.items():
            out[f"pooled_{tag}"] = self._ap_ratio(self.y[m], w[m])
            out[f"macro_{tag}"] = self.macro(w, m)
        return out


def _registered(edges: pd.DataFrame, weights: np.ndarray, gi, universe: pd.DataFrame, fast: dict) -> dict:
    """Replace the fast pooled value by the registered code path after checking they agree (both
    undefined, NaN, is agreement: a context with no evaluation pair in its universe)."""
    from grnbench.eval.metrics import evaluate_ranking
    r = evaluate_ranking(edges.assign(weight=weights), gi, with_curve=False, universe=universe)
    a, b = float(r["auprc_ratio"]), float(fast["pooled"])
    if not (abs(a - b) <= TOL or (np.isnan(a) and np.isnan(b))):
        raise CheckFailed(f"fast scorer {fast['pooled']} != evaluate_ranking {r['auprc_ratio']}")
    return {**fast, "pooled": float(r["auprc_ratio"]), "epr": float(r["epr"])}


def decompose(E: pd.DataFrame, expr: pd.DataFrame, universe: pd.DataFrame, gi, labelable: set[str],
              n_perm: int) -> tuple[dict, dict]:
    """Every metric for base, real prior and controls on one context; per-seed values apart."""
    from grnbench.infer.cistarget import motif_rank_product

    genes = tuple(map(str, expr.columns))
    tfs_e = list(dict.fromkeys(E["source"].astype(str)))           # = _motif_rank_of_edges
    M, in_db, K, worst = motif_parts(tfs_e, genes)
    sets = motif_sets(tfs_e, genes)
    shared = shared_matrix(sets)
    ti = E["source"].astype(str).map({t: i for i, t in enumerate(tfs_e)}).to_numpy()
    gix = E["target"].astype(str).map({g: i for i, g in enumerate(genes)}).to_numpy()
    if pd.isna(gix).any():
        raise ValueError("edge targets outside the expression genes")
    ti, gix = ti.astype(int), gix.astype(int)

    # self-consistency: the local rank product equals the registered one for the unchanged prior
    ref = motif_rank_product(E, expr)
    mine = rank_product(E, M[ti, gix].astype(float))
    if not (mine[["source", "target"]].equals(ref[["source", "target"]])
            and np.array_equal(mine["weight"].to_numpy(), ref["weight"].to_numpy())):
        raise CheckFailed("self-consistency failed: rank_product != cistarget.motif_rank_product")

    up = {g.upper(): i for i, g in enumerate(genes)}
    tf_has = {t.upper(): k > 0 for t, k in zip(tfs_e, K)}
    src_u, tgt_u = universe["source"].str.upper().to_numpy(), universe["target"].str.upper().to_numpy()
    tgt_ok = np.array([in_db[up[t]] and t in labelable for t in tgt_u], dtype=bool)
    restrict = tgt_ok & np.array([tf_has.get(s, False) for s in src_u], dtype=bool)
    sc = Scorer(universe, gi, {"T": tgt_ok, "R": restrict}, E)
    base_w = E["weight"].to_numpy(dtype=float)
    base_rank = E["weight"].rank(ascending=False, method="average").to_numpy()

    rec = {f"base_{k}": v for k, v in _registered(E, base_w, gi, universe, sc.score(base_w)).items()}
    for v in FIXED:
        w = rank_product_weights(base_rank, variant_matrix(v, M, in_db, K, worst, shared=shared)[ti, gix].astype(float))
        rec.update({f"{v}_{k}": x for k, x in _registered(E, w, gi, universe, sc.score(w)).items()})
    worst_tie_applies = bool((M[np.ix_(K > 0, in_db)] < worst).all())
    violations = builtin_check_violations(rec, worst_tie_applies)
    if violations:
        raise CheckFailed("built-in checks: " + "; ".join(violations))
    seeds = {}
    keys = (name_keys(tfs_e), name_keys(genes))
    for v in PERMS:
        per = []
        for s in range(n_perm):
            w = rank_product_weights(base_rank, variant_matrix(v, M, in_db, K, worst, s, keys)[ti, gix].astype(float))
            res = sc.score(w)
            if s == 0:
                _registered(E, w, gi, universe, res)                   # the check only
            per.append(res)
        for k in METRICS:
            vals = np.array([p[k] for p in per], dtype=float)
            rec[f"{v}_{k}"] = float(np.mean(vals[np.isfinite(vals)])) if np.isfinite(vals).any() else np.nan
            seeds[f"{v}_{k}"] = vals.tolist()

    gold_tf = {t: int(sc.y[b].sum()) for t, b in zip(sc.tf_names, sc.blocks)}
    worst_u = np.array([not tf_has.get(s, False) or not in_db[up[t]] for s, t in zip(src_u, tgt_u)])
    k_s = pd.Series(K, dtype=float)
    g_s = pd.Series([gold_tf.get(t.upper(), 0) for t in tfs_e], dtype=float)
    rho = float(k_s.corr(g_s, method="spearman")) if k_s.nunique() > 1 and g_s.nunique() > 1 else np.nan
    has_motif = np.array([len(st) > 0 for st in sets], dtype=bool)
    same = np.array([has_motif[i] and any(j != i and sets[j] == sets[i] for j in np.flatnonzero(has_motif))
                     for i in range(len(sets))], dtype=bool)
    n_motif = int(has_motif.sum())
    rec.update(
        n_tfs_motif=n_motif,
        frac_tfs_same_motifs=float(same.sum() / n_motif) if n_motif else np.nan,
        frac_tfs_shared_motifs=float((shared[has_motif].any(1)).sum() / n_motif) if n_motif else np.nan,
        n_universe=int(len(universe)), n_tfs=len(tfs_e), n_tfs_no_motif=int((K == 0).sum()),
        n_tfs_macro=len(sc.macro_blocks()), n_tfs_macro_R=len(sc.macro_blocks(restrict)),
        n_genes=len(genes), genes_not_in_db=int((~in_db).sum()),
        genes_ensembl=int(pd.Index(genes).str.match(r"^ENS[A-Z]{0,4}G\d+").sum()),
        frac_universe_worst_tie=float(worst_u.mean()), gold_in_worst_tie=int(sc.y[worst_u].sum()),
        gold_in_universe=int(sc.y.sum()), frac_universe_T=float(tgt_ok.mean()), gold_in_T=int(sc.y[tgt_ok].sum()),
        frac_universe_R=float(restrict.mean()), gold_in_R=int(sc.y[restrict].sum()), spearman_K_gold=rho,
        builtin_checks="held" if worst_tie_applies else "held (worst_tie on U'' not applicable: an "
                       "in-database gene has the worst rank for a TF with a motif)")
    return rec, seeds


def builtin_check_violations(rec: dict, worst_tie_applies: bool = True) -> list[str]:
    """PLAN.md section 2: controls that cannot move a metric by construction must leave it equal to
    the base model's, exactly (NaN equals NaN). A violation is a code error, never data."""
    expect = [(v, m) for v in ("tf_mean", "tf_kexp", "tf_presence") for m in ("macro", "macro_T", "macro_R")]
    expect += [(v, m) for v in ("tf_presence", "db_only") + (("worst_tie",) if worst_tie_applies else ())
               for m in ("pooled_R", "macro_R")]
    expect += [("tfmean_db", "macro_R")]
    out = []
    for v, m in expect:
        a, b = rec.get(f"{v}_{m}", np.nan), rec.get(f"base_{m}", np.nan)
        if not (a == b or (np.isnan(a) and np.isnan(b))):
            out.append(f"{v}_{m} {a!r} != base_{m} {b!r}")
    return out


# ----------------------------------------------------------------------------- one dataset
class _KeptPlan:
    """The one attribute score_models reads from a plan: the strategy (from decisions.json)."""

    def __init__(self, strategy: str):
        self.strategy = strategy

    def value(self, key, default=None):
        return self.strategy if key == "strategy" else default


def registered_rescore(tdir: Path, models: list[str], tfs: list[str], eval_pairs: set,
                       strategies: dict, work: Path) -> tuple[dict, dict]:
    """real_bench.score_models on the kept processed.h5ad of each arm (artifacts re-written to
    ``work``), and the re-written top edges compared with the kept files."""
    records, plans = {}, {}
    for arm in ARMS:
        recs: list[dict] = []
        for cdir in ctx_dirs(tdir, arm):
            i = int(cdir.name[3:])
            recs += [{"processed": None}] * (i - len(recs))
            info = json.loads((cdir / "context.json").read_text(encoding="utf-8")) \
                if (cdir / "context.json").exists() else {}
            recs.append({"processed": str(cdir / "processed.h5ad"), "context": info.get("context")})
        records[arm], plans[arm] = recs, _KeptPlan(strategies[arm])
    scores = real_bench.score_models(records, plans, {a: a for a in ARMS}, list(models), tfs,
                                     eval_pairs, artifacts_dir=work)
    edges = {f"{arm}/{c.name}/{m}": compare_edge_files(c / f"{m}.tsv.gz", work / arm / c.name / f"{m}.tsv.gz")
             for arm in ARMS for c in ctx_dirs(tdir, arm) for m in models}
    return scores, edges


def _finite_mean(values) -> float:
    x = [float(v) for v in values if v is not None and np.isfinite(v)]
    return float(np.mean(x)) if x else np.nan          # score_models: np.mean of the finite contexts


def run_unit(unit: dict, refs: str | None, n_perm: int = N_PERM,
             workroot: str | None = None) -> tuple[dict, dict]:
    """One dataset of one set -> (row, per-context permutation seeds)."""
    import anndata as ad
    from grnbench.autoprep import resolve_strategy
    from grnbench.eval.metrics import GoldIndex
    from grnbench.infer import prepare_expression, run_method
    from grnbench.tfs import resolve_tf_names

    pairs = [tuple(p) for p in unit["pairs"]]
    check_stage1_pairs(pairs)
    if refs is None:
        raise ValueError("--refs is required (labelable targets come from the full CollecTRI)")
    out = {"id": unit["id"], "exp": unit["exp"], "task": unit["task"], "experiment": EXPERIMENT,
           "status": "completed", "flags": list(unit.get("plan_flags", [])),
           "notes": list(unit.get("plan_notes", [])), "n_perm": int(n_perm),
           "pairs": [list(p) for p in pairs], "recompute": list(unit["recompute"]),
           "reference_sha256": real_bench.reference_fingerprints(refs)}
    if any(f.startswith("MISSING") for f in out["flags"]):
        out["status"] = "excluded"
        return out, {}
    tdir = Path(unit["run"]).expanduser() / "artifacts" / unit["task"]
    committed = json.loads((_resolve(unit["committed"]) / "rows" / f"{unit['task']}.json")
                           .read_text(encoding="utf-8"))
    sel, ev = read_pairs(tdir / "reference_selection.tsv"), read_pairs(tdir / "reference_evaluation.tsv")
    tfs = sorted(set(sel["source"]) | set(ev["source"]))                  # = run_real_task's tfs
    eval_pairs = set(zip(ev["source"], ev["target"]))
    gi = GoldIndex.coerce(eval_pairs)
    organism = json.loads(committed["truth"]).get("organism")
    if committed.get("auto_candidate") is not None:
        out["auto_is_baseline"] = committed.get("auto_candidate") == committed.get("baseline_candidate")
    ref_path = Path(refs) / f"collectri_{organism}.tsv"
    labelable = set(read_pairs(ref_path)["target"].str.upper())
    decisions = json.loads((tdir / "decisions.json").read_text(encoding="utf-8"))
    strategies = {arm: resolve_strategy(decisions["chosen_plans"][arm]["choices"])[0] for arm in ARMS}
    out["strategies"] = strategies

    # 1. calibration (gate 1): the registered scoring code on the kept inputs
    work = Path(tempfile.mkdtemp(prefix="sens_", dir=workroot))
    try:
        reg, edges = registered_rescore(tdir, unit["recompute"], tfs, eval_pairs, strategies, work)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    out["registered"] = {k: v for k, v in reg.items() if k.endswith(("_auprc_ratio", "_epr"))}
    out["edges"] = edges
    out["gate"] = {}
    errors = [f"score_models {k[2:-6]}: {v}" for k, v in reg.items() if k.endswith("_error")]
    if errors:          # an exception swallowed inside score_models is a runtime error, rule (c)
        raise ScoreModelsError("; ".join(errors))
    for m in unit["recompute"]:
        for arm in ARMS:
            for metric in ("auprc_ratio", "epr"):
                key = f"m_{m}_{arm}_{metric}"
                d = abs(float(reg.get(key, np.nan)) - float(committed.get(key, np.nan)))
                out["gate"][f"{m}|{arm}|{metric}"] = d
                if metric == "auprc_ratio" and not d <= TOL:
                    out["flags"].append(f"GATE score {m}/{arm}: |recomputed - committed| = {d}")
    for k, e in edges.items():
        if not e["match"]:
            out["flags"].append(f"GATE edges {k}: {e}")

    # 2. decomposition, per arm and context
    out["ctx"], seeds = [], {}
    for arm in ARMS:
        assume = strategies[arm] != "raw_counts"
        for cdir in ctx_dirs(tdir, arm):
            adata = ad.read_h5ad(cdir / "processed.h5ad")
            genes = pd.Index(adata.var_names).astype(str)
            sources = resolve_tf_names(genes, tfs)
            universe = pd.DataFrame([(s, g) for s in sources for g in genes if s != g],
                                    columns=["source", "target"])
            expr, _ = prepare_expression(adata, assume_normalized=assume)
            for combo, base in pairs:
                E, _ = run_method(base, adata, tf_list=tfs, assume_normalized=assume)
                rec, sd = decompose(E, expr, universe, gi, labelable, n_perm)
                out["ctx"].append({"arm": arm, "ctx": cdir.name, "pair": combo, **rec})
                seeds[f"{arm}|{cdir.name}|{combo}"] = sd

    # 3. dataset level: mean over finite contexts, as score_models; internal consistency
    out["means"], out["seed_means"] = {}, {}
    for arm in ARMS:
        for combo, base in pairs:
            recs = [r for r in out["ctx"] if r["arm"] == arm and r["pair"] == combo]
            for col in [c for c in (recs[0] if recs else {}) if c.endswith(METRICS + ("epr",))]:
                out["means"][f"{arm}|{combo}|{col}"] = _finite_mean(r[col] for r in recs)
            for v, b in (("real", combo), ("base", base)):
                mine = out["means"].get(f"{arm}|{combo}|{v}_pooled", np.nan)
                theirs = float(reg.get(f"m_{b}_{arm}_auprc_ratio", np.nan))
                d = abs(mine - theirs)
                if not (d <= TOL or (np.isnan(mine) and np.isnan(theirs))):    # both undefined: agree
                    raise CheckFailed(f"internal consistency: {arm}/{v} pooled {mine} differs from "
                                      f"score_models by {d}")
            for k in SEED_KEYS:
                per_ctx = [seeds[f"{arm}|{r['ctx']}|{combo}"][f"perm_db_{k}"] for r in recs]
                out["seed_means"][f"{arm}|{combo}|perm_db_{k}"] = [
                    _finite_mean(vals[s] for vals in per_ctx) for s in range(n_perm)]
    return out, seeds


def _scrub(text: str) -> str:
    """No personal home path in a committed row (as real_bench.record_screening_failures)."""
    for home in {str(Path.home()), str(Path.home().resolve())}:
        if home and home not in ("/", "."):
            text = text.replace(home, "~")
    return text


def _run_and_save(args) -> tuple[str, str]:
    unit, refs, outdir, n_perm, workroot, commit = args
    path = Path(outdir) / "rows" / f"{unit['id']}.json"
    if path.exists():
        try:
            status = json.loads(path.read_text(encoding="utf-8")).get("status")
            if status in ("completed", "excluded", "check_failed"):   # check_failed: fix the code first
                return unit["id"], f"skipped ({status})"
        except (OSError, ValueError):
            pass
    try:
        row, seeds = run_unit(unit, refs, n_perm, workroot)
    except Stage2Refused:
        raise
    except Exception as exc:  # noqa: BLE001 - recorded in the row
        # a failed code check is not transient (no resubmission); anything else, ScoreModelsError
        # included, is an error that every resubmission with the same commit retries (record,
        # exclusion rule c: the "once" is the operator's one resubmission)
        row, seeds = {"id": unit["id"], "exp": unit["exp"], "task": unit["task"],
                      "status": "check_failed" if isinstance(exc, CheckFailed) else "error",
                      "flags": [f"{'CHECK' if isinstance(exc, CheckFailed) else 'ERROR'} "
                                f"{type(exc).__name__}: {exc}"[:500]],
                      "traceback": traceback.format_exc()[-3000:]}, {}
    row["code_commit"] = commit                 # every status: one OUT never mixes commits (run audits it)
    row["flags"] = [_scrub(f) for f in row.get("flags", [])]
    if "traceback" in row:
        row["traceback"] = _scrub(row["traceback"])
    if seeds:
        _atomic_json(Path(outdir) / "perm_seeds" / f"{unit['id']}.json", seeds)
    _atomic_json(path, row)
    return unit["id"], row["status"] + ("" if not row["flags"] else ": " + "; ".join(row["flags"])[:300])


# ----------------------------------------------------------------------------- plan
def plan_units(runs: list[str], results_root: Path, pairs=STAGE1_PAIRS) -> list[dict]:
    """One unit per committed row of each set, so every registered dataset is accounted for.
    Missing artifacts or missing committed stage-1 scores exclude a unit (``plan_flags``); a
    run-folder row that is missing or differs from the committed row is a note (``plan_notes``):
    gate 1 compares with the committed row anyway."""
    pairs = [tuple(p) for p in pairs]
    check_stage1_pairs(pairs)
    need = {m for p in pairs for m in p}
    units = []
    for spec in runs:
        exp, _, run = spec.partition("=")
        if not exp or not run:
            raise ValueError(f"--run expects EXP=RUN_DIR, got {spec!r}")
        committed = Path(results_root) / exp
        rows = sorted((committed / "rows").glob("*.json"))
        if not rows or not (committed / "REPORT_REAL.md").exists():
            raise ValueError(f"{_portable(committed)}: no committed rows/ or REPORT_REAL.md "
                             "(the calibration reference)")
        run_dir = Path(run).expanduser()
        for c in rows:
            row = json.loads(c.read_text(encoding="utf-8"))
            flags, notes = [], []
            kept = run_dir / "rows" / c.name
            if not kept.exists():
                notes.append("no row in the run folder; gate 1 compares with the committed row of "
                             "the repository")
            elif _sha256(kept) != _sha256(c):
                notes.append("ROW of the run folder differs from the committed row (sha256); "
                             "gate 1 compares with the committed row of the repository")
            tdir = run_dir / "artifacts" / c.stem
            missing = [n for n in ("reference_selection.tsv", "reference_evaluation.tsv", "decisions.json")
                       if not (tdir / n).exists()]
            missing += [f"{arm}/ctx*/processed.h5ad" for arm in ARMS if not ctx_dirs(tdir, arm)]
            if missing:
                flags.append(f"MISSING artifacts: {', '.join(missing)}")
            recompute = [m for m in RECOMPUTABLE if f"m_{m}_auto_auprc_ratio" in row]
            if not need <= set(recompute):
                flags.append(f"MISSING committed scores for {sorted(need - set(recompute))}")
            units.append({"id": f"{exp}__{c.stem}", "exp": exp, "task": c.stem,
                          "run": _portable(run_dir), "committed": _portable(committed),
                          "pairs": [list(p) for p in pairs], "recompute": recompute,
                          "committed_row_sha256": _sha256(c), "plan_flags": flags, "plan_notes": notes})
    return units


# ----------------------------------------------------------------------------- report
def declared_lines(lines) -> dict:
    """{(heading, A, B): line} for the declared-comparison tables of a REPORT_REAL.md."""
    out, head = {}, None
    for line in lines:
        if line.startswith("#### Declared comparisons, "):
            head = line
        elif line.startswith("#"):
            head = None
        elif head and line.startswith("| ") and not line.startswith(("| A |", "|---")):
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            out[(head, cells[0], cells[1])] = line
    return out


def gate2(rows: list[dict], committed_dir: Path) -> tuple[list[dict], int, int]:
    """models_section on the committed rows, with the recomputed values substituted for every
    dataset that has them, must print the committed declared lines (AUPRC ratio, both arms) of
    every pair with a recomputed member. Other family members (e.g. C6's GRNBoost2 and
    RegDiffusion) keep their committed values, so the Holm family is the registered one.
    Returns (checks, datasets substituted, datasets in the committed set)."""
    committed = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((committed_dir / "rows").glob("*.json"))]
    df = pd.DataFrame(committed)
    recomputed = sorted({m for r in rows for m in r.get("recompute", [])})
    by_task = {r["task"]: r["registered"] for r in rows if r.get("registered")}
    for m in recomputed:
        for arm in ARMS:
            for metric in ("auprc_ratio", "epr"):
                key = f"m_{m}_{arm}_{metric}"
                if key in df:
                    df[key] = [float(by_task[t].get(key, np.nan)) if t in by_task else v
                               for t, v in zip(df["task"], df[key])]
    want = declared_lines((committed_dir / "REPORT_REAL.md").read_text(encoding="utf-8").splitlines())
    got = declared_lines(real_bench.models_section(df))
    checks = []
    for key, line in want.items():
        head, a, b = key
        if "AUPRC ratio" not in head or not ({a, b} & set(recomputed)):
            continue
        checks.append({"heading": head, "a": a, "b": b, "committed": line, "recomputed": got.get(key),
                       "equal": got.get(key) == line})
    return checks, int(df["task"].isin(by_task).sum()), int(len(df))


def _series(rows, arm, combo, variant, metric) -> pd.Series:
    return pd.Series({r["task"]: r["means"].get(f"{arm}|{combo}|{variant}_{metric}", np.nan) for r in rows},
                     dtype=float)


def _test(frame: pd.DataFrame) -> dict | None:
    from grnbench.eval.comparison import compare_units
    frame = frame.dropna()
    if len(frame) < 3:
        return None
    res = compare_units(frame)["pairs"][0]
    return {**res, "n": int(len(frame))}


def declared_tests(rows: list[dict], arm: str, combo: str) -> list[dict]:
    """T1, T2 and T3 for one set and arm; Holm over the directed hypotheses actually tested (6, or
    2 fewer per test left untested with fewer than 3 analysed datasets)."""
    from grnbench.eval.comparison import holm_adjust
    frames = [("T1", "pooled: real > perm_db", "perm_db", "pooled"),
              ("T2", "macro_R: real > base", "base", "macro_R"),
              ("T3", "macro_R: real > gene_tf", "gene_tf", "macro_R")]
    recs, raw = [], []
    for name, label, other, metric in frames:
        res = _test(pd.DataFrame({"real": _series(rows, arm, combo, "real", metric),
                                  other: _series(rows, arm, combo, other, metric)}))
        recs.append({"test": name, "label": label, "pair": combo, "res": res})
        if res is not None:
            raw += [res["p_a_gt_b"], res["p_b_gt_a"]]
    adj = iter(holm_adjust(raw) if raw else [])
    for rec in recs:
        if rec["res"] is not None:
            rec["p_a_gt_b_holm"], rec["p_b_gt_a_holm"] = next(adj), next(adj)
    for rec in recs:
        rec["holm_family"] = len(raw)
    return recs


def _lift(rows, arm, combo, variant, metric, other: str = "base") -> float:
    d = (_series(rows, arm, combo, variant, metric) - _series(rows, arm, combo, other, metric)).dropna()
    return float(d.mean()) if len(d) else np.nan


def ratio_interval(rows, arm, combo, num: tuple[str, str], seed: int = 0) -> tuple[float, float, float]:
    """mean(num[0] - num[1]) / L over datasets (pooled), with a bootstrap 95% interval over datasets
    (N_BOOT draws). Ratios are point estimates: they can exceed 1 or be negative."""
    f = pd.DataFrame({"n": _series(rows, arm, combo, num[0], "pooled") - _series(rows, arm, combo, num[1], "pooled"),
                      "L": _series(rows, arm, combo, "real", "pooled") - _series(rows, arm, combo, "base", "pooled")}).dropna()
    if len(f) < 2 or f["L"].mean() == 0:
        return np.nan, np.nan, np.nan
    idx = np.random.default_rng(seed).integers(0, len(f), (N_BOOT, len(f)))
    with np.errstate(divide="ignore", invalid="ignore"):
        bs = f["n"].to_numpy()[idx].mean(1) / f["L"].to_numpy()[idx].mean(1)
    return float(f["n"].mean() / f["L"].mean()), float(np.nanpercentile(bs, 2.5)), float(np.nanpercentile(bs, 97.5))


def permutation_rank_p(rows, arm, combo) -> tuple[float, int, int]:
    """T1 robustness to the one shared database: rank of the real prior's set mean (pooled) among
    the set means of the permuted priors, seed by seed (a seed is one keyed random database shared
    by every dataset). Minimum 1 / (n_perm + 1). Returns (p, datasets, seeds)."""
    ok = [r for r in rows if np.isfinite(r["means"].get(f"{arm}|{combo}|real_pooled", np.nan))
          and np.isfinite(r["seed_means"].get(f"{arm}|{combo}|perm_db_pooled", [np.nan])).all()]
    if not ok:
        return np.nan, 0, 0
    lens: dict[int, list[str]] = {}
    for r in ok:
        lens.setdefault(len(r["seed_means"][f"{arm}|{combo}|perm_db_pooled"]), []).append(str(r.get("id", "?")))
    if len(lens) > 1:                    # never truncated: that would change p_perm's resolution silently
        ids = "; ".join(f"{n} seeds: {len(v)} row(s), e.g. {v[0]}" for n, v in sorted(lens.items()))
        raise ValueError(f"{arm}|{combo}: rows mix permutation seed counts {sorted(lens)} ({ids}); rerun them "
                         "with one --n-perm")
    real = np.mean([r["means"][f"{arm}|{combo}|real_pooled"] for r in ok])
    perm = np.mean([r["seed_means"][f"{arm}|{combo}|perm_db_pooled"] for r in ok], axis=0)
    return float((1 + np.sum(perm >= real)) / (1 + len(perm))), len(ok), len(perm)


def _supported(tests: list[dict], name: str, side: str = "p_a_gt_b_holm") -> bool:
    r = next((t for t in tests if t["test"] == name), None)
    return bool(r and r.get("res") is not None and r[side] < ALPHA)


def target_level_shortfalls(tests: list[dict], L: float, d_tgt: float, l_art: float, l_gene: float,
                            *, t1_robust: bool) -> list[str]:
    """The conditions of "supports target-level evidence" beyond T1 and T2 that fail (empty: none)."""
    half = 0.5 * L
    out = []
    if not _supported(tests, "T3"):
        out.append("T3 not supported: the within-TF lift is not shown to exceed a TF-generic gene "
                   "propensity (gene_tf)")
    if not (np.isfinite(d_tgt) and d_tgt >= half):
        out.append("D < 0.5 L: the gap to the random assignment is less than half of the lift")
    if not (np.isfinite(l_art) and l_art < half):
        out.append("L_art >= 0.5 L: the design-effect control reaches at least half of the lift")
    if not (np.isfinite(l_gene) and l_gene < half):
        out.append("L_gene >= 0.5 L: a TF-generic gene propensity reaches at least half of the lift")
    if not t1_robust:
        out.append(T1_NOT_ROBUST)
    return out


def stage1_verdict(tests: list[dict], L: float, d_tgt: float, l_art: float, l_gene: float, *,
                   t1_robust: bool) -> str:
    """PLAN.md section 5 for one pair, with an explicit precedence (the more cautious row wins) and
    the record's additions: T2's reverse joins T1's; "target-level" also needs T3, L_gene < 0.5 L
    and T1's permutation rank p < ALPHA (record EXP-20261003-GB-1, decision rule)."""
    t1, t1_rev = _supported(tests, "T1"), _supported(tests, "T1", "p_b_gt_a_holm")
    t2, t2_rev = _supported(tests, "T2"), _supported(tests, "T2", "p_b_gt_a_holm")
    if not np.isfinite(L) or L <= 0:
        return "inconclusive"
    art = np.isfinite(l_art) and l_art >= 0.5 * L
    if t1_rev or t2_rev:
        return "design effects only"
    if t1 and t2:
        return "mixed" if target_level_shortfalls(tests, L, d_tgt, l_art, l_gene, t1_robust=t1_robust) \
            else "target-level"
    if t1 or t2:
        return "mixed" if art else "inconclusive"
    return "design effects only" if art else "inconclusive"


def t2_reverse_only(tests: list[dict], L: float, d_tgt: float, l_art: float, l_gene: float, *,
                    t1_robust: bool) -> tuple[bool, str | None]:
    """Reported, never decides (the registered rule 2 stands): did rule 2 fire through T2's reverse
    alone while T1 or T3 is supported? A prior that carries target identity but ranks within TF
    worse than the base also lowers macro_R. Returns (flag, the verdict of the registered rules with
    T2's reverse masked), or (False, None) when the flag is not set."""
    verdict = stage1_verdict(tests, L, d_tgt, l_art, l_gene, t1_robust=t1_robust)
    if not (verdict == "design effects only" and _supported(tests, "T2", "p_b_gt_a_holm")
            and not _supported(tests, "T1", "p_b_gt_a_holm")
            and (_supported(tests, "T1") or _supported(tests, "T3"))):
        return False, None
    masked = [dict(t, p_b_gt_a_holm=1.0) if t["test"] == "T2" else t for t in tests]
    return True, stage1_verdict(masked, L, d_tgt, l_art, l_gene, t1_robust=t1_robust)


def _test_summary(rec) -> dict:
    out = {k: v for k, v in rec.items() if k != "res"}
    res = rec["res"]
    if res is not None:
        out.update(n=res["n"], difference=res["difference"], interval=res["simultaneous_ci"],
                   p_a_gt_b=res["p_a_gt_b"], p_b_gt_a=res["p_b_gt_a"], status=res["status"])
    return out


def _fmt_test(rec) -> str:
    res = rec["res"]
    if res is None:
        return f"| {rec['test']} | {rec['label']} | {rec['pair']} | <3 | - | - | not tested | not tested |"
    lo, hi = res["simultaneous_ci"]
    ci = "unbounded" if lo is None else f"[{lo:+.4f}, {hi:+.4f}]"
    return (f"| {rec['test']} | {rec['label']} | {rec['pair']} | {res['n']} | {res['difference']:+.4f} | "
            f"{ci} | {rec['p_a_gt_b_holm']:.3g} | {rec['p_b_gt_a_holm']:.3g} |")


def _registered_hashes() -> dict[str, str]:
    path = REPO_ROOT / "experiments" / "datasets.tsv"
    if not path.exists():
        return {}
    tab = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    return {Path(p).name: h for p, h in zip(tab["path"], tab["sha256"])
            if p.startswith(("cistarget/", "references/"))}


def _snapshot_inputs() -> list[str] | None:
    """The ``cistarget/`` and ``references/`` paths of this record's D-020 snapshot manifest (None
    when the manifest is not in the checkout). The report names each one RUN_INPUTS.json lacks."""
    path = REPO_ROOT / "experiments" / f"{EXPERIMENT}.snapshot.tsv"
    if not path.exists():
        return None
    tab = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    return [p for p in tab["path"] if p.startswith(("cistarget/", "references/"))]


COMPONENTS = ("tf_mean", "tf_kexp", "tf_presence", "db_only", "worst_tie", "tfmean_db", "gene_mean", "gene_tf",
              "gene_tf_np")


def _n_lift(rows, arm, combo, variant, metric, other: str = "base") -> int:
    """Datasets with a finite value for both sides of ``_lift`` (its n)."""
    return int((_series(rows, arm, combo, variant, metric) - _series(rows, arm, combo, other, metric)).notna().sum())


def arm_analysis(ok: list[dict], arm: str, combo: str) -> dict:
    """Declared tests, decision quantities and verdict of one set and arm (datasets passing gate 1)."""
    tests = declared_tests(ok, arm, combo)
    L = _lift(ok, arm, combo, "real", "pooled")
    d_tgt = _lift(ok, arm, combo, "real", "pooled", other="perm_db")
    l_art = _lift(ok, arm, combo, "tfmean_db", "pooled")
    l_gene = _lift(ok, arm, combo, "gene_tf", "pooled")
    p_perm, n_sets, n_seeds = permutation_rank_p(ok, arm, combo)
    assessable = n_seeds > 0 and 1 / (n_seeds + 1) < ALPHA
    t1_robust = bool(assessable and p_perm < ALPHA)
    verdict = stage1_verdict(tests, L, d_tgt, l_art, l_gene, t1_robust=t1_robust)
    t1 = _supported(tests, "T1")
    both = t1 and _supported(tests, "T2")
    rule4 = []                      # "mixed" without both T1 and T2 is record rule 4; reported, never decides
    if verdict == "mixed" and not both:
        rule4 = [f"{'T2' if t1 else 'T1'} not supported",
                 "L_art >= 0.5 L: the design-effect control reaches at least half of the lift"]
        if t1 and not t1_robust:
            rule4.append(T1_NOT_ROBUST)
    via_t2, without_t2 = t2_reverse_only(tests, L, d_tgt, l_art, l_gene, t1_robust=t1_robust)
    return {
        "tests": tests, "L": L, "delta_target": d_tgt, "L_art": l_art, "L_gene": l_gene,
        "n_quantities": {"L": _n_lift(ok, arm, combo, "real", "pooled"),
                         "D": _n_lift(ok, arm, combo, "real", "pooled", other="perm_db"),
                         "L_art": _n_lift(ok, arm, combo, "tfmean_db", "pooled"),
                         "L_gene": _n_lift(ok, arm, combo, "gene_tf", "pooled")},
        "macro_R_lift": _lift(ok, arm, combo, "real", "macro_R"),
        "macro_R_over_gene_tf": _lift(ok, arm, combo, "real", "macro_R", other="gene_tf"),
        "epr_lift": _lift(ok, arm, combo, "real", "epr"),
        "components": {v: _lift(ok, arm, combo, v, "pooled") for v in COMPONENTS},
        "ratios": {"D/L": ratio_interval(ok, arm, combo, ("real", "perm_db")),
                   "L_art/L": ratio_interval(ok, arm, combo, ("tfmean_db", "base")),
                   "L_gene/L": ratio_interval(ok, arm, combo, ("gene_tf", "base"))},
        "p_perm": p_perm, "n_perm_datasets": n_sets, "n_seeds": n_seeds, "assessable": assessable,
        "t1_robust": t1_robust, "verdict": verdict,
        "shortfalls": target_level_shortfalls(tests, L, d_tgt, l_art, l_gene, t1_robust=t1_robust)
        if verdict == "mixed" and both else [],
        "rule4_reasons": rule4,
        "untested": [t["test"] for t in tests if t["res"] is None], "n_analysed": len(ok),
        "downgraded": verdict == "mixed" and both and not t1_robust and stage1_verdict(
            tests, L, d_tgt, l_art, l_gene, t1_robust=True) == "target-level",
        "rule2_via_t2_reverse_only": via_t2, "verdict_without_t2_reverse": without_t2}


def secondary_checks(primary: dict, secondary: dict | None) -> tuple[list[str], list[str]]:
    """Record: the baseline arm and EPR only have to agree in sign with the auto arm. A
    disagreement is named next to the verdict and never changes it. Returns (disagree,
    not_compared): a pair with a non-finite value (a test not run in one or both arms, no finite
    lift) has no sign, so it is not compared rather than called a disagreement."""
    def sign(x):
        return np.sign(x) if x is not None and np.isfinite(x) else None
    disagree, not_compared = [], []

    def check(text, x, ref):
        if sign(x) is None or sign(ref) is None:
            not_compared.append(text)
        elif sign(x) != sign(ref):
            disagree.append(text)
    check(f"EPR(real) - EPR(base) {primary['epr_lift']:+.4f} vs L {primary['L']:+.4f}",
          primary["epr_lift"], primary["L"])
    if secondary is not None:
        check(f"baseline arm L {secondary['L']:+.4f} vs {primary['L']:+.4f}", secondary["L"], primary["L"])
        for tp, ts in zip(primary["tests"], secondary["tests"]):
            dp = tp["res"]["difference"] if tp["res"] else np.nan
            ds = ts["res"]["difference"] if ts["res"] else np.nan
            check(f"baseline arm {ts['test']} mean difference {ds:+.4f} vs {dp:+.4f}", ds, dp)
    return disagree, not_compared


def secondary_text(disagree: list[str], not_compared: list[str]) -> str:
    """'all agree.' only when every pair was compared and none disagrees."""
    if not disagree and not not_compared:
        return "all agree."
    out = ("DISAGREE: " + "; ".join(disagree) + ".") if disagree else "no compared pair disagrees."
    if not_compared:
        out += " Not compared (no finite value in one or both arms): " + "; ".join(not_compared) + "."
    return out


def missing_row(unit: dict) -> dict:
    """The row of a planned unit that never wrote one (a killed worker, a job cut by its time limit,
    or a run still going): reported and excluded under exclusion rule (c), never dropped."""
    return {"id": unit["id"], "exp": unit["exp"], "task": unit["task"], "status": "missing",
            "flags": [f"MISSING row: no rows/{unit['id']}.json was written (worker killed or job cut "
                      "before the unit finished); rule (c): resubmit once with the same commit, then it "
                      "is excluded"]}


def report(outdir: Path) -> str:
    outdir = Path(outdir)
    rows = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((outdir / "rows").glob("*.json"))]
    planned = (outdir / "sens_tasks.json").exists()
    units = json.loads((outdir / "sens_tasks.json").read_text(encoding="utf-8")) if planned else []
    inputs = json.loads((outdir / "RUN_INPUTS.json").read_text(encoding="utf-8")) \
        if (outdir / "RUN_INPUTS.json").exists() else {}
    plat = inputs.pop("platform", None)
    inputs.pop("invocations", None)
    code_file = (outdir / "CODE_COMMIT").read_text(encoding="utf-8").strip() \
        if (outdir / "CODE_COMMIT").exists() else None
    registered = _registered_hashes()
    have = {r.get("id") for r in rows}
    missing = [u for u in units if u["id"] not in have]
    n_seeds = sorted({r.get("n_perm") for r in rows if "n_perm" in r})
    deviations = [f"DEVIATION: permutation seeds {x}, registered {N_PERM} (p_perm resolution 1/{x + 1})"
                  for x in n_seeds if x != N_PERM]          # reported, not refused (the mocks use 8 or 20)
    summary = {"experiment": EXPERIMENT, "type": "exploratory", "stage": 1, "platform": plat,
               "code_commit_file": code_file, "deviations": deviations,
               "missing_units": [u["id"] for u in missing] if planned else None, "sets": {}}
    lines = [f"# {EXPERIMENT}: motif-prior lift decomposition, stage 1 (Pearson, tier B)", ""]
    if missing:
        lines += [f"**INCOMPLETE: {len(missing)} planned unit(s) have no row (listed as excluded, status "
                  "missing, in their set below). Exclusion rule (c): resubmit once with the same commit; "
                  "a unit still without a row after that stays excluded and the tests use the other "
                  "datasets.**", ""]
    lines += ["EXPLORATORY. Unit = dataset. Generated by `benchmarks/motif_sensitivity.py report`; "
              "no number in this file is typed by hand.", "",
              f"- units planned: {len(units) if planned else 'unknown'}; rows: {len(rows)}; "
              f"permutation seeds: {n_seeds}"]
    lines += [f"- {d}" for d in deviations]
    if planned:
        lines.append(f"- units without a row: {len(missing)}"
                     + (f" ({', '.join(u['id'] for u in missing)})" if missing else ""))
    rows += [missing_row(u) for u in missing]          # after the real rows: their order is unchanged
    if plat:
        lines.append(f"- platform: {plat.get('os')}; CPU {plat.get('cpu') or '-'}; Python {plat.get('python')}; "
                     f"numpy {plat.get('numpy')}; BLAS {plat.get('blas') or '-'}")
    lines.append(f"- CODE_COMMIT (copied by hpc/engine_bench.sbatch at the last submission): {code_file}"
                 if code_file is not None else "- CODE_COMMIT: not in this folder")
    for kind, files in sorted(inputs.items()):
        for name, digest in sorted(files.items()):
            status = ("matches experiments/datasets.tsv" if registered.get(name) == digest
                      else "NOT the registered hash" if name in registered else "not in experiments/datasets.tsv")
            lines.append(f"- input {kind}/{name}: sha256 {digest[:16]}... ({status})")
    snapshot = _snapshot_inputs()
    if snapshot is None:
        lines.append(f"- input check: experiments/{EXPERIMENT}.snapshot.tsv not found")
    summary["missing_inputs"] = [p for p in snapshot or []
                                 if p.split("/", 1)[1] not in (inputs.get(p.split("/", 1)[0]) or {})]
    for p in summary["missing_inputs"]:
        lines.append(f"- input {p}: MISSING from RUN_INPUTS.json (registered in "
                     f"experiments/{EXPERIMENT}.snapshot.tsv; not read by this run or not recorded)")
    lines.append("")
    exps = [e for e in REGISTERED_SETS if any(r["exp"] == e for r in rows)]
    for exp in exps + sorted({r["exp"] for r in rows} - set(exps)):
        set_rows = [r for r in rows if r["exp"] == exp]
        ok = [r for r in set_rows if r.get("status") == "completed" and not r.get("flags")]
        bad = [r for r in set_rows if r not in ok]
        check_failed = [r for r in set_rows if r.get("status") == "check_failed"]
        committed_dir = next((_resolve(u["committed"]) for u in units if u["exp"] == exp), None) \
            or REPO_ROOT / "benchmarks" / "results" / exp
        role = REGISTERED_SETS.get(exp, "unregistered set")
        g2, n_sub, n_set = gate2(set_rows, committed_dir) if (committed_dir / "REPORT_REAL.md").exists() \
            else ([], 0, 0)
        g2_ok = bool(g2) and all(c["equal"] for c in g2) and n_sub == n_set and not bad
        not_interpretable = bool(check_failed) or (not g2_ok and not bad)
        same = [r.get("auto_is_baseline") for r in ok if r.get("auto_is_baseline") is not None]
        # the code commit of every written row (a missing row has none): one folder, one commit
        commits = sorted({r.get("code_commit") for r in set_rows if r.get("status") != "missing"}, key=str)
        mixed = len(commits) > 1 or (code_file is not None and any(c != code_file for c in commits))
        shown = ", ".join(c or "none recorded" for c in commits) or "none (no row written)"
        sset = summary["sets"][exp] = {
            "role": role, "units": len(set_rows), "gate_pass": len(ok),
            "excluded": [{"id": r["id"], "status": r.get("status"), "reasons": r.get("flags", [])} for r in bad],
            "gate2": {"pass": g2_ok, "substituted": n_sub, "committed_datasets": n_set, "checks": g2},
            "not_interpretable": not_interpretable, "auto_is_baseline": [int(sum(same)), len(same)],
            "code_commits": commits, "mixed_code_commits": mixed,
            "arms": {}}
        lines += [f"## {exp} ({role})", "",
                  f"- datasets: {len(set_rows)}; passed gate 1 and analysed: {len(ok)}; "
                  f"excluded (reported, never silently): {len(bad)}",
                  f"- datasets where the auto pipeline is the baseline pipeline (committed rows): "
                  f"{sum(same)}/{len(same)}; on these the two arms coincide",
                  f"- code commit of the rows: {shown}"]
        if mixed:
            lines += [f"- **MIXED CODE COMMITS**: the rows of this set come from {shown}"
                      + (f" and CODE_COMMIT is {code_file}" if code_file is not None else "")
                      + ". One results folder holds one commit: rule (c) resubmits with the same commit, and "
                      "after a code fix every unit is rerun at the fixed commit into a new results folder."]
        for r in bad:           # every reason in full: the size of every gate-1 failure (record, platform note)
            flags = r.get("flags", [])
            lines.append(f"  - {r['task']} ({r.get('status')}): {len(flags)} reason(s)")
            for f in flags:
                if not f.startswith("GATE ") and len(f) > 400:        # free text only (errors, checks)
                    f = f[:400] + f" [... {len(f) - 400} characters cut; full text in SENSITIVITY.json]"
                lines.append(f"    - {f}")
        for r in set_rows:
            for note in r.get("notes", []):
                lines.append(f"- note, {r['task']}: {note}")
        held = [c.get("builtin_checks", "") for r in ok for c in r.get("ctx", [])]
        lines.append(f"- built-in checks (constant-row and tie controls leave their metrics equal to the "
                     f"base): held on {sum(h.startswith('held') for h in held)}/{len(held)} contexts"
                     + (f"; worst_tie on U'' not applicable on {sum('not applicable' in h for h in held)}"
                        if any("not applicable" in h for h in held) else ""))
        if check_failed:
            lines += ["", f"**A code check failed on {len(check_failed)} dataset(s): this set is not "
                      "interpretable until the code is fixed (check_failed units are not resubmitted; after "
                      "the fix every unit is rerun at the fixed commit into a new results folder and the fix "
                      "is a dated deviation naming both commits).**"]
        lines += ["", "### Calibration", "",
                  "Gate 1, per dataset: every recomputed model's AUPRC ratio vs the committed row "
                  "(max |diff| over datasets and arms; tolerance 1e-12) and the re-written top edges vs "
                  "the kept files.", "", "| model | arm | datasets | max abs diff | edge files matching |",
                  "|---|---|---|---|---|"]
        for m in sorted({m for r in set_rows for m in r.get("recompute", [])}):
            for arm in ARMS:
                key = f"{m}|{arm}|auprc_ratio"
                ds = [r["gate"][key] for r in set_rows if key in r.get("gate", {})]
                ed = [e["match"] for r in set_rows for k, e in r.get("edges", {}).items()
                      if k.startswith(f"{arm}/") and k.endswith(f"/{m}")]
                lines.append(f"| {m} | {arm} | {len(ds)} | {max(ds) if ds else float('nan'):.3g} | "
                             f"{sum(ed)}/{len(ed)} |")
        lines += ["", f"Gate 2, per set: declared-comparison lines (AUPRC ratio) recomputed with "
                  f"`real_bench.models_section` vs the committed REPORT_REAL.md, recomputed values "
                  f"substituted for {n_sub}/{n_set} datasets, every dataset through gate 1: "
                  f"**{'pass' if g2_ok else 'FAIL'}** ({sum(c['equal'] for c in g2)}/{len(g2)} lines equal).", ""]
        for c in g2:
            lines.append(f"- {'equal' if c['equal'] else 'DIFFERENT'}: `{c['heading'][5:]}`")
            lines.append(f"  - committed:  `{c['committed']}`")
            if not c["equal"]:
                lines.append(f"  - recomputed: `{c['recomputed']}`")
        lines.append("")
        if not g2_ok and not bad:
            lines += ["**Gate 2 failed although every dataset passed gate 1: the registered report code "
                      "path does not reproduce. The tests below are not interpretable until the cause is "
                      "found and recorded.**", ""]
        elif not g2_ok:
            lines += [f"Gate 2 cannot pass with {len(bad)} dataset(s) excluded or failing gate 1; the "
                      f"tests below use the {len(ok)} datasets that passed gate 1.", ""]
        for combo, base in STAGE1_PAIRS:
            res = {arm: arm_analysis(ok, arm, combo) for arm in ARMS}
            for arm in (PRIMARY_ARM,) + tuple(a for a in ARMS if a != PRIMARY_ARM):
                a = res[arm]
                tests = a["tests"]
                verdict = "not interpretable" if not_interpretable else a["verdict"]
                label, meaning = VERDICTS[verdict]
                if not ok:
                    meaning = "Not tested: no dataset of this set passed gate 1."
                extra = ""
                if verdict == "mixed" and a["shortfalls"]:
                    extra += " Conditions of the strongest verdict that fail: " + "; ".join(a["shortfalls"]) + "."
                if verdict == "mixed" and any(s.startswith(("T3", "L_gene")) for s in a["shortfalls"]):
                    extra += " The lift is consistent with a TF-generic gene propensity."
                if verdict == "mixed" and a["rule4_reasons"]:      # the rule's own L_art and L, not the bootstrap
                    missing = a["rule4_reasons"][0].split()[0]
                    extra += (f" Record rule 4: only one of T1 and T2 is supported ({missing} is not) and the "
                              f"design-effect control reaches at least half of the lift (L_art/L "
                              f"{a['L_art'] / a['L']:.2f} >= 0.5).")
                    if T1_NOT_ROBUST in a["rule4_reasons"]:
                        extra += f" {T1_NOT_ROBUST}: reported only, it does not change this verdict (record)."
                if a["downgraded"]:
                    extra += (" Verdict downgraded from 'supports target-level evidence' to 'mixed' by the "
                              "record's robustness condition on T1 only.")
                if verdict == "design effects only" and a["rule2_via_t2_reverse_only"]:   # reported, never decides
                    sup = [n for n in ("T1", "T3") if _supported(tests, n)]
                    extra += (" Rule 2 fired through T2's reverse only: the rank-product combination lowers macro_R "
                              f"below the base within TF. {' and '.join(sup)} {'are' if len(sup) > 1 else 'is'} "
                              "supported. Without T2's reverse the registered rules would give "
                              f"'{VERDICTS[a['verdict_without_t2_reverse']][0]}'. The registered verdict and its "
                              "consequence stand; quote this together with D/L, L_art/L and L_gene/L.")
                    if np.isfinite(a["L_art"]) and a["L_art"] < 0.5 * a["L"]:
                        extra += " The design-effect control reaches less than half of the lift."
                if verdict != "not interpretable" and ok and a["untested"]:
                    none_run = {"T1", "T2"} <= set(a["untested"])
                    extra += (f" {', '.join(a['untested'])} not tested (fewer than 3 analysed datasets with a "
                              "value; the family shrinks, record): 'not supported' here means not tested"
                              + (", and the verdict rests on the point estimates L and L_art/L only."
                                 if none_run else "."))
                    if none_run:
                        meaning = meaning.replace("neither T1 nor T2 is supported", "neither T1 nor T2 was run")
                secondary = secondary_checks(a, res.get("baseline")) \
                    if arm == PRIMARY_ARM and verdict != "not interpretable" else None
                if secondary is not None:
                    extra += (" Secondary checks (sign agreement; they never change the verdict): "
                              + secondary_text(*secondary))
                kind = "declared tests" if arm == PRIMARY_ARM else "secondary (same tests, own Holm family)"
                ratio = ", ".join(f"{k} {v[0]:.2f} [{v[1]:.2f}, {v[2]:.2f}]" for k, v in a["ratios"].items())
                lines += [f"### {combo} vs {base}, {arm} pipeline: {kind}", "",
                          f"Paired t over datasets, one one-sided test per direction (compare_units), Holm over the "
                          f"{tests[0]['holm_family']} directed hypotheses of this table (alpha {ALPHA}).",
                          SHARED_DATABASE, "",
                          "| test | hypothesis (A > B) | pair | n | mean A-B | 95% interval | p(A>B) Holm | "
                          "p(B>A) Holm |",
                          "|---|---|---|---|---|---|---|---|"] + [_fmt_test(t) for t in tests] + [""]
                nq = a["n_quantities"]
                lines += [f"- registered lift L = pooled(real) - pooled(base): {a['L']:+.4f} (n = {nq['L']})",
                          f"- target-identity gap (T1) D = pooled(real) - pooled(perm_db): {a['delta_target']:+.4f} "
                          f"(n = {nq['D']})",
                          f"- within-TF lift on U'' (T2) = macro_R(real) - macro_R(base): {a['macro_R_lift']:+.4f}",
                          f"- within-TF lift over the gene propensity (T3) = macro_R(real) - macro_R(gene_tf): "
                          f"{a['macro_R_over_gene_tf']:+.4f}",
                          f"- design-effect control L_art = pooled(tfmean_db) - pooled(base): {a['L_art']:+.4f} "
                          f"(n = {nq['L_art']})",
                          f"- gene-propensity control L_gene = pooled(gene_tf) - pooled(base): {a['L_gene']:+.4f} "
                          f"(n = {nq['L_gene']})",
                          f"- ratios to L (point estimate [bootstrap 95% over datasets]; they can exceed 1 or be "
                          f"negative, the 0.5 cuts use the point estimates): {ratio}",
                          "- components (pooled lift of each control): "
                          + ", ".join(f"{k} {v:+.4f}" for k, v in a["components"].items()),
                          f"- EPR(real) - EPR(base): {a['epr_lift']:+.4f}",
                          f"- T1 permutation rank (set mean of the real prior among the set means of "
                          f"the keyed permuted priors, {a['n_perm_datasets']} datasets, {a['n_seeds']} seeds): "
                          f"p = {a['p_perm']:.3g}" + ("" if a["assessable"] else " (not assessable: too few seeds)"),
                          f"- **verdict ({arm}): {label}**. {meaning}{extra}", ""]
                t1_supported = _supported(tests, "T1")
                sset["arms"][arm] = {
                    "pair": [combo, base], "role": "declared" if arm == PRIMARY_ARM else "secondary",
                    "tests": [_test_summary(t) for t in tests],
                    **{k: a[k] for k in ("L", "delta_target", "macro_R_lift", "macro_R_over_gene_tf", "L_art",
                                         "L_gene", "epr_lift", "components", "ratios", "shortfalls",
                                         "rule4_reasons", "untested", "n_analysed", "n_quantities")},
                    "p_perm_T1": a["p_perm"], "verdict": verdict, "verdict_label": label,
                    "t1_robust_to_database": (a["p_perm"] < ALPHA) if (t1_supported and a["assessable"]) else None,
                    "downgraded_by_t1_robustness": a["downgraded"],
                    "rule2_via_t2_reverse_only": a["rule2_via_t2_reverse_only"],
                    "verdict_without_t2_reverse": a["verdict_without_t2_reverse"]}
                if secondary is not None:
                    sset["arms"][arm]["secondary_disagreements"], sset["arms"][arm]["secondary_not_compared"] = \
                        secondary
            lines += descriptive_tables(ok, combo) + diagnostics_table(ok, combo)
    text = "\n".join(lines) + "\n"
    (outdir / "SENSITIVITY.md").write_text(text, encoding="utf-8")
    _atomic_json(outdir / "SENSITIVITY.json", summary)
    return text


def descriptive_tables(rows, combo) -> list[str]:
    """Every variant against the base model, every metric (descriptive; Holm within each row's
    two directions only, not the declared family)."""
    out = []
    boot = np.random.default_rng(0)
    for arm in ARMS:
        out += [f"### Descriptive, {arm} pipeline: variant - base (not the declared family)", "",
                "| metric | variant | n | mean lift | 95% interval | p(>0) Holm 2 | p(<0) Holm 2 | "
                "share of the real lift [bootstrap 95%] |", "|---|---|---|---|---|---|---|---|"]
        for metric in METRICS:
            for v in VARIANTS:
                pair = pd.DataFrame({"variant": _series(rows, arm, combo, v, metric),
                                     "base": _series(rows, arm, combo, "base", metric),
                                     "real": _series(rows, arm, combo, "real", metric)}).dropna()
                res = _test(pair[["variant", "base"]])
                if res is None:
                    out.append(f"| {metric} | {v} | {len(pair)} | - | - | - | - | - |")
                    continue
                lo, hi = res["simultaneous_ci"]
                ci = "unbounded" if lo is None else f"[{lo:+.3f}, {hi:+.3f}]"
                lv = (pair["variant"] - pair["base"]).to_numpy()
                lr = (pair["real"] - pair["base"]).to_numpy()
                idx = boot.integers(0, len(lv), (N_BOOT, len(lv)))
                with np.errstate(divide="ignore", invalid="ignore"):
                    bs = lv[idx].mean(1) / lr[idx].mean(1)
                share = (f"{lv.mean() / lr.mean():.2f} [{np.nanpercentile(bs, 2.5):.2f}, "
                         f"{np.nanpercentile(bs, 97.5):.2f}]") if lr.mean() != 0 else "-"
                out.append(f"| {metric} | {v} | {res['n']} | {res['difference']:+.3f} | {ci} | "
                           f"{res['p_a_gt_b_adjusted']:.3g} | {res['p_b_gt_a_adjusted']:.3g} | {share} |")
        out.append("")
    return out


def diagnostics_table(rows, combo) -> list[str]:
    keys = ("n_universe", "n_tfs", "n_tfs_no_motif", "n_tfs_motif", "frac_tfs_same_motifs",
            "frac_tfs_shared_motifs", "n_tfs_macro", "n_tfs_macro_R", "n_genes",
            "genes_not_in_db", "genes_ensembl",
            "frac_universe_worst_tie", "gold_in_worst_tie", "gold_in_universe", "frac_universe_T",
            "gold_in_T", "frac_universe_R", "gold_in_R", "spearman_K_gold")
    out = ["### Diagnostics (per context, mean over contexts then over datasets)", "",
           "| quantity | " + " | ".join(ARMS) + " |", "|---|" + "---|" * len(ARMS)]
    for k in keys:
        vals = []
        for arm in ARMS:
            per = [_finite_mean(c.get(k, np.nan) for c in r["ctx"] if c["arm"] == arm and c["pair"] == combo)
                   for r in rows]
            vals.append(_finite_mean(per))
        out.append(f"| {k} | " + " | ".join(f"{v:.4g}" for v in vals) + " |")
    return out + [""]


# ----------------------------------------------------------------------------- CLI
def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("plan", help="list the units (datasets) of each set; check rows and artifacts")
    a.add_argument("--run", action="append", required=True, help="EXP-ID=RUN_DIR (kept artifacts of that run)")
    a.add_argument("--out", required=True)
    a.add_argument("--results-root", default=str(REPO_ROOT / "benchmarks" / "results"),
                   help="committed results (rows/, REPORT_REAL.md): the calibration reference")
    a.add_argument("--pairs", default=",".join(f"{c}:{b}" for c, b in STAGE1_PAIRS),
                   help="COMBO:BASE,... (stage 1: pearson_motif:pearson_signed only)")
    b = sub.add_parser("run")
    b.add_argument("--tasks", required=True)
    b.add_argument("--out", required=True)
    b.add_argument("--refs", default=None, help="folder with collectri_{human,mouse}.tsv (required)")
    b.add_argument("--jobs", type=int, default=1)
    b.add_argument("--index", type=int, default=None)
    b.add_argument("--chunk", type=int, default=1)
    b.add_argument("--workroot", default=None)
    b.add_argument("--n-perm", type=int, default=N_PERM, help="permutation seeds 0..n-1 (registered: 100)")
    b.add_argument("--models", default=None, help="accepted from engine_bench.sbatch; refit models are refused")
    b.add_argument("--keep-artifacts", action="store_true", help="accepted from engine_bench.sbatch; ignored")
    b.add_argument("--max-attempts", type=int, default=1, help="accepted from engine_bench.sbatch and "
                   "not used: a unit ending with a runtime error is a row with status error that every "
                   "resubmission reruns (check_failed rows never); the 'once' of exclusion rule (c) is "
                   "the operator's one resubmission with the same commit")
    c = sub.add_parser("report")
    c.add_argument("--out", required=True)
    args = p.parse_args(argv)
    try:
        if args.cmd == "plan":
            pairs = [tuple(x.split(":", 1)) for x in args.pairs.split(",") if x.strip()]
            units = plan_units(args.run, Path(args.results_root), pairs)
            Path(args.out).parent.mkdir(parents=True, exist_ok=True)
            Path(args.out).write_text(json.dumps(units, indent=1), encoding="utf-8")
            flagged = [u for u in units if u["plan_flags"] or u["plan_notes"]]
            print(f"{len(units)} units -> {args.out}; excluded by plan: "
                  f"{sum(bool(u['plan_flags']) for u in units)}; with notes: {sum(bool(u['plan_notes']) for u in units)}")
            for u in flagged:
                print(f"  {u['id']}: {'; '.join(u['plan_flags'] + u['plan_notes'])}")
            return 0
        if args.cmd == "report":
            print(report(Path(args.out)))
            return 0
        for m in (args.models or "").replace("+", ",").split(","):
            if any(m.strip().startswith(base) for base in REFIT_BASES):
                raise Stage2Refused(STAGE2_REFUSAL.format(model=m.strip()))
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        units = json.loads(Path(args.tasks).read_text(encoding="utf-8"))
        if args.index is not None:
            units = units[args.index * args.chunk:(args.index + 1) * args.chunk]
        for u in units:
            check_stage1_pairs([tuple(x) for x in u["pairs"]])
        if args.refs is None:
            p.error("--refs is required")
        problems = preflight_problems(units, args.refs)
        if problems:
            print("preflight failed (no unit was run, no row written):\n  " + "\n  ".join(problems),
                  file=sys.stderr)
            return 2
        commit = code_commit()
        problems = resume_problems(Path(args.out), commit, args.n_perm)
        if problems:
            print("refused (no unit was run, no row written):\n  " + "\n  ".join(problems), file=sys.stderr)
            return 2
        payload, why = run_inputs_payload(Path(args.out) / "RUN_INPUTS.json", input_fingerprints(args.refs),
                                          {"platform": platform_info(), "code_commit": commit})
        if payload is None:
            print(f"refused (no unit was run, no row written): {why}", file=sys.stderr)
            return 2
        Path(args.out).mkdir(parents=True, exist_ok=True)
        if args.workroot:                       # engine_bench.sbatch passes a path it does not create
            Path(args.workroot).mkdir(parents=True, exist_ok=True)
        _atomic_json(Path(args.out) / "RUN_INPUTS.json", payload)
        jobs = [(u, args.refs, args.out, args.n_perm, args.workroot, commit) for u in units]
        failed = skipped_check_failed = False
        results = pool_results(_run_and_save, jobs, args.jobs) if args.jobs > 1 else map(_run_and_save, jobs)
        try:
            for uid, status in results:
                failed |= status.startswith(("error", "check_failed"))
                skipped_check_failed |= status == "skipped (check_failed)"
                print(f"{uid}: {status}", flush=True)
        except BrokenProcessPool as exc:
            print(f"worker died ({exc}); finished rows are saved, resubmit to resume", flush=True)
            return 2
        if skipped_check_failed:
            print(f"skipped check_failed unit(s): {CHECK_FAILED_RULE}", file=sys.stderr)
        return 2 if failed or skipped_check_failed else 0
    except Stage2Refused as exc:
        print(f"refused: {exc}", file=sys.stderr)
        return 2
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
