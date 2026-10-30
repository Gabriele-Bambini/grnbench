"""Automatic preprocessing plan for GRN inference.

``diagnose`` loads the data and reads its context (matrix, metadata, README
files, download pages). ``plan_preprocessing`` profiles the data
(:mod:`grnbench.profile`), lets the knowledge-base selector choose one option
per decision point (:mod:`grnbench.selector`, :mod:`grnbench.kb`) and turns
the choices into a :class:`PrepPlan`: one :class:`Decision` per executable
parameter, each with the value, the reason, the evidence source and the
alternatives. ``execute_plan`` runs it (one ``run_prep`` per cellular
context).

The plan maps 1:1 onto :func:`grnbench.pipeline.run_prep` and the
``grnbench prep`` command line (``PrepPlan.cli_args``), so a plan edited in the
interface is exactly reproducible from a terminal. Numeric cell thresholds are
MAD outlier bounds on the log scale (scater convention): lower bound on
detected genes at 3 MADs, upper bounds on genes and counts at 5 MADs, and
mitochondrial fraction at median + 3 MADs, never below 5 %.
"""
from __future__ import annotations

import json
import math
import re
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Any

import numpy as np

from grnbench.context import DatasetContext, build_context, fetch_url_text, organism_note
from grnbench.qc import QCReport, Thresholds, run_qc

LOWER_NMADS = 3.0
UPPER_NMADS = 5.0
MITO_NMADS = 3.0
MITO_FLOOR = 5.0
MIN_GENES_FLOOR = 50
SUBSAMPLE_ABOVE = 20_000
SUBSAMPLE_TO = 10_000
DEFAULT_N_HVG = 1000
DOUBLET_MIN_CELLS = 1000
DROPLET_PLATFORMS = {"10x Chromium", "Drop-seq", "inDrop", "Seq-Well", "SPLiT-seq", "sci-RNA-seq"}

THRESHOLD_KEYS = ("min_genes_per_cell", "min_cells_per_gene", "max_genes_per_cell",
                  "max_counts_per_cell", "max_pct_mito_cell", "mito_warn", "mito_fail")
# Tissues whose cells are naturally rich in mitochondrial RNA (Osorio & Cai 2021): a
# dataset-wide median above the default 25 % gate is biology there, not a failed sample.
# Whole words only: as substrings, 'renal', 'atri', 'ventric', 'septum' and 'cardi' also matched adrenal
# gland, pediatric/psychiatric, brain ventricles, the lateral septum and the cardia of the stomach.
# Matched on the normalised text (`is_high_mito_tissue`), so snake_case, hyphenated, camelCase, digit-suffixed and
# plural spellings (skeletal_muscle_tissue, left-ventricle, SkeletalMuscle, heart1, hearts, endocardiums) match like
# the spaced labels. A lower-case run of words (skeletalmuscle) does not, on purpose: finding a word inside a run of
# letters is what relaxed liverwort, deliver and heartwood. The hepatic flexure of the colon is not liver, in any
# spelling ('hepatic flexure', 'hepaticflexure').
HIGH_MITO_TISSUES = re.compile(
    r"\b(?:livers?|\w*hepat(?!\w*flexure)\w*(?! flexure)|hearts?|cardiac|cardiomyocytes?|\w*cardi(?:ums?|al)"
    r"|atri(?:ums?|a|al)|atrioventricular|sinoatrial|interventricular|(?:left|right) ventricles?|kidneys?|renal"
    r"|\w*nephr\w*|muscles?|brown adipose)\b")
# Word breaks inside a label: lower -> Upper followed by a lower-case letter ('SkeletalMuscle'; 'iPSC' stays one
# word), and letter <-> digit ('heart1', 'Liver2Kidney').
_WORD_BREAKS = re.compile(r"(?<=[a-z])(?=[A-Z][a-z])|(?<=[A-Za-z])(?=\d)|(?<=\d)(?=[A-Za-z])")


def tissue_words(tissue) -> str:
    """The tissue text the high-mito matcher reads: words split at ``_WORD_BREAKS``, lower case, '_' and '-' read
    as spaces, runs of white space as one space. Only the matcher reads it: the context stores and shows the
    label in lower case, as written otherwise."""
    text = _WORD_BREAKS.sub(" ", str(tissue or ""))
    return " ".join(re.sub(r"[_-]", " ", text.lower()).split())


def is_high_mito_tissue(tissue) -> bool:
    """Whether a tissue text names a high-mitochondrial tissue: ``HIGH_MITO_TISSUES`` on ``tissue_words``."""
    return bool(HIGH_MITO_TISSUES.search(tissue_words(tissue)))
RUN_PREP_KEYS = ("strategy", "n_hvg", "batch_key", "group_key", "stratify_key",
                 "subsample_cells", "subsample_mode", "drop_doublets", "orientation", "exclude_qc_genes",
                 "per_batch_bounds", "cells_per_metacell", "cell_calling", "force")
MAX_CONTEXTS = 4   # largest cell types inferred separately; the interface can pick others
DEFAULT_SEED = 42  # execute_plan's seed when none is given (the interface cleans with it); = `grnbench prep --seed`
_CLI_FLAGS = {
    "strategy": "--strategy", "n_hvg": "--n-hvg", "batch_key": "--batch-key",
    "group_key": "--group-key", "stratify_key": "--stratify-key",
    "subsample_cells": "--subsample-cells", "subsample_mode": "--subsample-mode",
    "orientation": "--orientation",
    "min_genes_per_cell": "--min-genes-per-cell", "min_cells_per_gene": "--min-cells-per-gene",
    "max_genes_per_cell": "--max-genes-per-cell", "max_counts_per_cell": "--max-counts-per-cell",
    "max_pct_mito_cell": "--max-pct-mito", "cells_per_metacell": "--cells-per-metacell",
    "cell_calling": "--cell-calling", "mito_warn": "--mito-warn", "mito_fail": "--mito-fail",
}


# The knowledge-base decisions that set each plan decision: a user choice of one of them is
# named next to the value it changed ("chosen: fixed").
KB_DECISIONS_OF = {
    "cell_calling": ("cell_calling",),
    "min_genes_per_cell": ("cell_qc_thresholds",), "max_genes_per_cell": ("cell_qc_thresholds",),
    "max_counts_per_cell": ("cell_qc_thresholds",), "max_pct_mito_cell": ("cell_qc_thresholds",),
    "per_batch_bounds": ("cell_qc_thresholds",),
    "drop_doublets": ("doublets",), "exclude_qc_genes": ("gene_filter",),
    "strategy": ("normalization", "batch", "sparsity"), "cells_per_metacell": ("sparsity",),
    "n_hvg": ("hvg",), "subsample_cells": ("subsample",), "subsample_mode": ("subsample",),
    "tf_source": ("regulators",),
    "context_scope": ("context_scope",),
}
NOT_IN_ENGINE_PLAN = "not in the engine's own plan"


@dataclass
class Decision:
    """One executable parameter of the plan.

    After a user override (a typed value or a knowledge-base choice other than the
    selector's), ``auto_value``, ``auto_reason`` and ``auto_source`` keep the executable value
    the engine would have used, its reason and its source, and ``source`` is ``"user"``.
    Without an override they are ``None`` and ``plan.json`` leaves them out.
    """
    key: str
    label: str
    value: Any
    reason: str
    source: str = "default"
    options: list | None = None
    auto_value: Any = None
    auto_reason: str | None = None
    auto_source: str | None = None

    def __post_init__(self) -> None:
        _check_trio(self)


def _check_trio(d: Decision) -> None:
    """The engine's value, reason and source come as a whole, and only on a user value."""
    if d.auto_reason is None:
        if d.auto_value is not None or d.auto_source is not None:
            raise ValueError(f"decision {d.key!r}: auto_value/auto_source without auto_reason")
    elif d.auto_source is None:
        raise ValueError(f"decision {d.key!r}: auto_reason without auto_source")
    elif d.source != "user":
        raise ValueError(f"decision {d.key!r}: the engine's value is kept only on a user value, "
                         f"but source is {d.source!r}")


@dataclass
class PrepPlan:
    decisions: dict[str, Decision] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    methods: list[dict] = field(default_factory=list)
    context: dict = field(default_factory=dict)
    profile: dict = field(default_factory=dict)
    selection: dict = field(default_factory=dict)
    choices: dict = field(default_factory=dict)
    contexts: list = field(default_factory=lambda: [{"name": "all_cells", "subset_key": None,
                                                     "subset_value": None}])
    # Never in ``warnings`` (n_warnings is a column of benchmark rows): for each best option this
    # machine cannot run, why and what it would take ("blocked": {decision: Availability + option}).
    option_notes: dict = field(default_factory=lambda: {"blocked": {}})

    def add(self, key, label, value, reason, source="default", options=None) -> None:
        self.decisions[key] = Decision(key, label, value, reason, source, options)

    def value(self, key: str, default=None):
        d = self.decisions.get(key)
        return default if d is None else d.value

    def chosen_options(self, key: str) -> list[str]:
        """The knowledge-base options chosen instead of the selector's that set decision ``key``."""
        results = self.selection.get("results", {})
        return [self.choices[k] for k in KB_DECISIONS_OF.get(key, ())
                if k in self.choices and k in results and self.choices[k] != results[k].get("chosen")]

    def set(self, key: str, value, reason: str = "set by the user") -> None:
        """A user value. The first override keeps the engine's value, reason and source next
        to it; setting the engine's value again makes the decision the engine's again."""
        d = self.decisions[key]
        if d.value == value:
            return
        if d.auto_reason is not None and value == d.auto_value:
            d.value, d.reason, d.source = d.auto_value, d.auto_reason, d.auto_source
            d.auto_value = d.auto_reason = d.auto_source = None
            return
        if d.source != "user":
            d.auto_value, d.auto_reason, d.auto_source = d.value, d.reason, d.source
        d.value, d.reason, d.source = value, reason, "user"

    def thresholds(self) -> Thresholds:
        return Thresholds(**{k: self.value(k) for k in THRESHOLD_KEYS
                             if k in self.decisions and self.value(k) is not None})

    def run_prep_kwargs(self, tf_list: list[str] | None = None, *, gene_mapping_path=None) -> dict:
        kw = {k: self.value(k) for k in RUN_PREP_KEYS
              if k in self.decisions and self.value(k) is not None}
        kw["drop_doublets"] = bool(self.value("drop_doublets", False))
        kw["force"] = bool(self.value("force", False))
        kw["thresholds"] = self.thresholds()
        if tf_list is not None:
            kw["tf_list"] = list(tf_list)
        if gene_mapping_path is not None:
            kw["gene_mapping_path"] = gene_mapping_path
        return kw

    def cli_args(self, input_path, out_dir, tf_list_path=None, *, context: dict | None = None,
                 bounds_path=None, gene_mapping_path=None, seed: int,
                 force: bool | None = None) -> list[str]:
        """The ``grnbench prep`` argv that reproduces this plan (one context) exactly.

        ``per_batch_bounds`` travels as a JSON file: pass the path it was written
        to as ``bounds_path`` (``execute_plan`` does this). ``seed`` is the seed the
        context ran with, always written: without ``--seed`` the command would
        replay at the CLI default, other cells for a screened context (seed 0).
        ``force`` is whether the run went past failed QC gates (None: the plan's value).
        ``--n-hvg`` is always written, the number or ``all`` (what ``run_prep`` runs when the plan's value is
        None or absent), so a recorded command never depends on the command-line default (no ``--n-hvg``
        adds a note that ``execute_plan``'s run does not have; LD's review of #162).
        """
        argv = [str(input_path), "--out", str(out_dir)]
        for key, flag in _CLI_FLAGS.items():
            v = self.value(key)
            if key == "n_hvg":
                argv += [flag, "all" if v is None else _fmt(v)]
            elif key in self.decisions and v is not None:
                argv += [flag, _fmt(v)]
        if self.value("drop_doublets"):
            argv.append("--drop-doublets")
        if self.value("force") if force is None else force:
            argv.append("--force")
        if self.value("exclude_qc_genes"):
            argv.append("--exclude-qc-genes")
        if self.value("per_batch_bounds") and bounds_path is not None:
            argv += ["--per-batch-bounds", str(bounds_path)]
        if context and context.get("subset_key"):
            argv += ["--subset-key", str(context["subset_key"]),
                     f"--subset-value={context['subset_value']}"]
        if tf_list_path is not None:
            argv += ["--tf-list", str(tf_list_path)]
        if gene_mapping_path is not None:
            argv += ["--gene-mapping", str(gene_mapping_path)]
        return argv + ["--seed", str(seed)]

    def to_dict(self) -> dict:
        return {"decisions": {k: _decision_dict(d) for k, d in self.decisions.items()},
                "warnings": list(self.warnings), "methods": list(self.methods),
                "context": self.context, "profile": self.profile,
                "selection": self.selection, "choices": self.choices,
                "contexts": list(self.contexts), "option_notes": self.option_notes}

    @classmethod
    def from_dict(cls, data: dict) -> "PrepPlan":
        plan = cls(decisions={k: _decision_from_dict(k, d) for k, d in data["decisions"].items()},
                   warnings=list(data.get("warnings", [])),
                   methods=list(data.get("methods", [])),
                   context=dict(data.get("context", {})),
                   profile=dict(data.get("profile", {})),
                   selection=dict(data.get("selection", {})),
                   choices=dict(data.get("choices", {})),
                   option_notes=dict(data.get("option_notes") or {"blocked": {}}))
        if data.get("contexts"):
            plan.contexts = list(data["contexts"])
        return plan

    def save(self, path) -> Path:
        path = Path(path)
        path.write_text(json.dumps(self.to_dict(), indent=2, default=str), encoding="utf-8")
        return path

    @classmethod
    def load(cls, path) -> "PrepPlan":
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))


_AUTO_FIELDS = ("auto_value", "auto_reason", "auto_source")


def _decision_dict(d: Decision) -> dict:
    """``asdict`` of a decision; the engine's value, reason and source only after an override,
    so a plan nobody edited is written exactly as before overrides kept them."""
    _check_trio(d)
    out = asdict(d)
    if d.auto_reason is None:
        for k in _AUTO_FIELDS:
            del out[k]
    return out


def _decision_from_dict(name: str, data: dict) -> Decision:
    unknown = sorted(set(data) - {f.name for f in fields(Decision)})
    if unknown:
        raise ValueError(f"decision {name!r}: unknown key(s) {', '.join(unknown)}")
    return Decision(**data)


def _keep_engine_values(plan: PrepPlan, engine: PrepPlan) -> None:
    """Next to every decision the choices changed, the engine's own value, reason and source."""
    for key, d in plan.decisions.items():
        e = engine.decisions.get(key)
        own = (e.value, e.reason, e.source) if e else (None, NOT_IN_ENGINE_PLAN, "knowledge base")
        if (d.value, d.reason, d.source) != own:
            d.auto_value, d.auto_reason, d.auto_source = own
            d.source = "user"


@dataclass
class Diagnosis:
    """Everything ``plan_preprocessing`` looked at; the interface plots it.

    ``knee_no_cliff``: the matrix looks unfiltered but the barcode-rank knee found no cell/empty-droplet
    cliff, so every barcode was kept (issue #116); the plan says so in a warning.
    """
    adata: Any
    ingest: dict
    qc: QCReport
    context: DatasetContext
    id_verdict: dict
    raw_droplets: bool = False
    cell_calling: dict | None = None
    knee_no_cliff: bool = False


def diagnose(path, *, orientation: str = "auto", text_files=(), urls=(),
             notes: str | None = None, fetch=fetch_url_text,
             thresholds: Thresholds | None = None,
             gene_mapping_path=None) -> Diagnosis:
    """Load ``path`` counts-verified, read its context and run the QC gates."""
    from grnbench import hygiene as hyg
    from grnbench.ingest import load_any

    adata, ing = load_any(
        path, orientation=orientation, gene_mapping_path=gene_mapping_path)
    hyg.sparsify_if_dense(adata)
    from grnbench.profile import looks_unfiltered
    cell_calling = None
    knee_no_cliff = False
    totals = np.asarray(adata.X.sum(axis=1)).ravel()
    if looks_unfiltered(totals):
        # Profile and thresholds must describe cells, not empty droplets.
        cell_calling = hyg.knee_cell_threshold(totals)
        if cell_calling is not None:
            adata = adata[totals >= cell_calling["threshold"]].copy()
        else:
            knee_no_cliff = True       # every barcode is kept: the plan warns (issue #116)
    id_verdict = hyg.detect_id_type(adata.var_names)
    ctx = build_context(path, adata, text_files=text_files, urls=urls, notes=notes,
                        fetch=fetch)
    from grnbench.profile import technology_from_counts
    hint, ones, vmax = technology_from_counts(adata.X)
    if hint != "unknown":
        # Count values are stronger evidence than prose: a README that says
        # Smart-seq2 over a matrix of UMI counts loses the vote (and the conflict
        # is still reported).
        ctx.add("technology", hint, "matrix",
                f"{ones:.0%} of non-zero values are 1, max value {vmax:g}", 4)
    batch_col, _ = hyg.detect_batch_column(adata)
    qc = run_qc(adata, thresholds or Thresholds(), batch_key=batch_col,
                id_verdict=id_verdict, n_dup_barcodes=hyg.duplicate_barcodes_report(adata))
    return Diagnosis(adata, ing.to_dict(), qc, ctx, id_verdict, raw_droplets=looks_unfiltered(totals),
                     cell_calling=cell_calling, knee_no_cliff=knee_no_cliff)


def mad_bounds(values, nmads: float, log: bool = True) -> tuple[float, float, float]:
    """``(lower, upper, median)`` at ``nmads`` scaled MADs; NaN bounds when MAD is 0."""
    x = np.log1p(np.asarray(values, dtype=float)) if log else np.asarray(values, dtype=float)
    med = float(np.median(x))
    mad = 1.4826 * float(np.median(np.abs(x - med)))
    back = np.expm1 if log else (lambda v: v)
    if mad == 0:
        return float("nan"), float("nan"), float(back(med))
    return float(back(med - nmads * mad)), float(back(med + nmads * mad)), float(back(med))


def _pct(mask) -> str:
    return f"{int(np.sum(mask))} cells ({np.mean(mask):.1%})"


def _mad_thresholds(genes, counts, mito, n_mito_genes: int) -> dict:
    """MAD bounds for one group of cells: ``{key: (value, reason)}``."""
    default = Thresholds()
    out = {}
    lo, _, med = mad_bounds(genes, LOWER_NMADS)
    if math.isnan(lo):
        out["min_genes_per_cell"] = (default.min_genes_per_cell,
                                     "no spread in detected genes: fixed default kept")
    else:
        v = int(max(MIN_GENES_FLOOR, round(lo)))
        out["min_genes_per_cell"] = (v, f"median {med:.0f} genes/cell; lower {LOWER_NMADS:g}-MAD bound "
                                        f"(log) = {lo:.0f}, floor {MIN_GENES_FLOOR}: removes {_pct(genes < v)} "
                                        f"(fixed default {default.min_genes_per_cell} would remove "
                                        f"{_pct(genes < default.min_genes_per_cell)})")
    _, hi, _ = mad_bounds(genes, UPPER_NMADS)
    out["max_genes_per_cell"] = ((None, "no spread in detected genes: no upper bound") if math.isnan(hi)
                                 else (int(math.ceil(hi)), f"upper {UPPER_NMADS:g}-MAD bound (gross "
                                       f"multiplets only) flags {_pct(genes > hi)}"))
    _, hi_c, med_c = mad_bounds(counts, UPPER_NMADS)
    out["max_counts_per_cell"] = ((None, "no spread in library size: no upper bound") if math.isnan(hi_c)
                                  else (float(math.ceil(hi_c)), f"median {med_c:.0f} counts/cell; upper "
                                        f"{UPPER_NMADS:g}-MAD bound flags {_pct(counts > hi_c)}"))
    if n_mito_genes == 0:
        out["max_pct_mito_cell"] = (None, "no mitochondrial genes in the matrix: filter off")
    else:
        _, hi_m, med_m = mad_bounds(mito, MITO_NMADS, log=False)
        bound = MITO_FLOOR if math.isnan(hi_m) else max(hi_m, MITO_FLOOR)
        out["max_pct_mito_cell"] = (round(bound, 1), f"median {med_m:.1f}% mito; median + {MITO_NMADS:g} "
                                    f"MADs = {hi_m:.1f}% (floor {MITO_FLOOR:g}%) flags {_pct(mito > bound)}")
    return out


_STRATEGY_OF = {"cpm_log1p": "cpm_log1p", "pearson_residuals": "pearson_residuals",
                "raw_counts": "raw_counts"}
_BATCH_STRATEGY = {"combat": "cpm_log1p_combat", "center_per_batch": "cpm_log1p_center_batch"}
_N_HVG = {"tf_plus_500": 500, "tf_plus_1000": 1000, "tf_plus_2000": 2000, "all_genes": None}
# Recorded only: no step reads tf_source, the TF list that runs is the caller's (issue #131).
_TF_SOURCE = {"user_file": "file", "dorothea": "dorothea:{organism}", "all_genes": "none"}
_QC_LABELS = {"min_genes_per_cell": "Min detected genes per cell",
              "max_genes_per_cell": "Max detected genes per cell",
              "max_counts_per_cell": "Max counts per cell",
              "max_pct_mito_cell": "Max % mitochondrial counts"}


def resolve_strategy(choices: dict) -> tuple[str, str, list[str]]:
    """One executable strategy from the normalization/batch/sparsity choices.

    Returns ``(strategy, reason, compatibility notes)``. Priority: sparsity
    (metacells/pseudobulk change the samples themselves) > batch correction >
    normalization, because each of the first two already includes a
    normalization step. An option with no plan step
    (``executors.PLAN_STEPS``) is an error, never a silent fallback.
    """
    from grnbench import executors

    notes = []
    norm = choices.get("normalization", "cpm_log1p")
    batch = choices.get("batch", "none")
    sparsity = choices.get("sparsity", "single_cells")
    for decision, option in (("normalization", norm), ("batch", batch), ("sparsity", sparsity)):
        if option not in executors.PLAN_STEPS[decision]:
            raise ValueError(f"no plan step for {decision}={option}: grnbench cannot execute it")
    if sparsity == "metacells":
        if batch in _BATCH_STRATEGY:
            notes.append(f"batch option '{batch}' not combined with metacells: metacells are built "
                         "inside each batch instead")
        if norm != "cpm_log1p":
            notes.append(f"metacells are depth-normalized with cpm_log1p (not '{norm}')")
        return "metacells", "sparse data aggregated into metacells (built within batch)", notes
    if sparsity == "pseudobulk":
        reason_notes = []
        if batch != "none":
            reason_notes.append(
                f"batch option '{batch}' is not applied after pseudobulk aggregation")
        if norm != "cpm_log1p":
            reason_notes.append(
                f"pseudobulk groups are normalized with cpm_log1p, not '{norm}'")
        reason = "counts summed per group and normalized with cpm_log1p"
        if reason_notes:
            reason += "; " + "; ".join(reason_notes)
        return "pseudobulk", reason, notes
    if batch in _BATCH_STRATEGY:
        if norm != "cpm_log1p":
            notes.append(f"batch option '{batch}' is applied on cpm_log1p, not '{norm}'")
        return _BATCH_STRATEGY[batch], f"cpm_log1p + batch option '{batch}'", notes
    return _STRATEGY_OF[norm], f"normalization '{norm}'", notes


def _user_forced(choices: dict | None, selection, decision: str) -> bool:
    """Whether ``choices`` sets ``decision`` to an option other than the selector's."""
    return bool(choices) and decision in choices and choices[decision] != selection.results[decision].chosen


def _check_runnable(selection, chosen: dict, choices: dict | None) -> None:
    """Fail loudly on any chosen option grnbench cannot run (D-004, D-006; issue #128): missing
    dependency, no executor, planned, not applicable or not offered (D-002). The selector never
    chooses one; a user-forced one is an error, never a silent substitution."""
    from grnbench.kb import DECISION_INDEX

    for decision, option in chosen.items():
        if decision not in DECISION_INDEX:
            continue
        a = DECISION_INDEX[decision].availability_full(option)
        if a.ok:
            continue
        who = "was set by the user" if _user_forced(choices, selection, decision) else "was chosen by the selector"
        more = f" Install: {a.install_hint}." if a.install_hint else (f" {a.detail}." if a.detail else "")
        raise ValueError(f"{decision}={option} {who}, but grnbench cannot run it: {a.reason}.{more}")


def plan_preprocessing(diag: Diagnosis, *, methods: list[dict] | None = None,
                       tf_list: list[str] | None = None, selection=None,
                       choices: dict | None = None) -> PrepPlan:
    """Executable plan from the knowledge-base selection (or explicit ``choices``).

    ``choices`` ({decision_id: option_id}) overrides the selector - that is how
    an alternative candidate pipeline or a user edit becomes a plan. Every decision the
    choices change keeps the engine's own value, reason and source (``Decision.auto_*``), read
    from the plan built with the selector's choices. An option
    grnbench cannot run is never chosen by the selector, and forcing one fails
    loudly (``ValueError`` naming the option, the reason and the install hint).
    """
    from grnbench.kb import DECISION_INDEX
    from grnbench.profile import build_profile
    from grnbench.selector import select

    profile = build_profile(diag.adata, diag.qc.metrics, diag.context, id_verdict=diag.id_verdict,
                            tf_list=tf_list, raw_droplets=diag.raw_droplets)
    selection = selection or select(profile)
    unknown = sorted(set(choices or {}) - set(DECISION_INDEX))
    if unknown:
        raise KeyError(f"unknown preprocessing decision(s): {unknown}")
    plan = _build_plan(diag, profile, selection, choices, methods=methods)
    own = selection.chosen()
    if choices and any(own.get(k) != v for k, v in choices.items()):
        # the engine's own plan is only compared: methods=[] skips listing the methods again
        _keep_engine_values(plan, _build_plan(diag, profile, selection, None, methods=[]))
    return plan


def _build_plan(diag: Diagnosis, profile, selection, choices: dict | None, *,
                methods: list[dict] | None) -> PrepPlan:
    import pandas as pd

    from grnbench import hygiene as hyg
    from grnbench.kb import DECISION_INDEX
    from grnbench.profile import MIN_CELLS_PER_CONTEXT

    adata, qc, ctx = diag.adata, diag.qc, diag.context
    chosen = {**selection.chosen(), **(choices or {})}
    _check_runnable(selection, chosen, choices)
    plan = PrepPlan(context=ctx.to_dict())
    plan.profile = profile.to_dict()
    plan.selection = selection.to_dict()
    plan.choices = dict(chosen)
    plan.option_notes = {"blocked": {
        r.decision: {"option": r.best_overall,
                     **asdict(DECISION_INDEX[r.decision].availability_full(r.best_overall))}
        for r in selection.results.values() if r.blocked_by_dependency}}

    def why(decision: str) -> tuple[str, str]:
        r = selection.results[decision]
        opt = chosen[decision]
        so = next(o for o in r.options if o.option == opt)
        src = "user" if _user_forced(choices, selection, decision) else "knowledge base"
        text = ("; ".join(x["text"] for x in so.reasons if x["delta"] > 0)
                or DECISION_INDEX[decision].option(opt).summary)
        return f"[{decision} = {opt}] {text}", src

    obs = adata.obs
    genes = obs["_n_genes"].to_numpy()
    counts = obs["_total_counts"].to_numpy()
    mito = obs["_pct_mito"].to_numpy()
    n_mito = qc.metrics.get("n_mito_genes", 0)

    # -- cell calling -------------------------------------------------------------- #
    if chosen.get("cell_calling") == "knee":
        call = diag.cell_calling
        if call is not None:
            called = (f"barcode-rank knee keeps {call.get('n_called', '?')} of {call.get('n_barcodes', '?')} "
                      f"barcodes (>= {call.get('threshold', float('nan')):.0f} counts)")
        elif diag.knee_no_cliff:      # issue #116: no "? of ? barcodes (>= nan counts)"
            called = ("the barcode-rank curve shows no clear cell/empty-droplet cliff, so the knee called no "
                      "cells; running this plan stops at cell calling")
        else:
            called = ("the matrix does not look unfiltered, so no cells were called at diagnosis; the knee "
                      "runs when the plan is executed")
        plan.add("cell_calling", "Cell calling", "knee", f"{why('cell_calling')[0]}; {called}",
                 "data", ["knee", None])

    # -- high-mitochondrial tissues: relax the dataset-level gate, keep per-cell bounds -- #
    tissue = (ctx.get("tissue") or "").lower()
    if is_high_mito_tissue(ctx.spelling("tissue", ctx.get("tissue"))):   # as written: camelCase breaks are words
        med = qc.metrics.get("median_pct_mito", 0.0)
        plan.add("mito_warn", "Dataset mito gate: warn above median %", 35.0,
                 f"tissue '{tissue}' is naturally rich in mitochondrial RNA (Osorio & Cai 2021); "
                 f"median here {med:.1f}%", _source(ctx, "tissue"))
        plan.add("mito_fail", "Dataset mito gate: fail above median %", 60.0,
                 "the default 25% would refuse healthy tissue; damaged cells are still removed by "
                 "the per-cell mitochondrial bound", _source(ctx, "tissue"))

    # -- QC gates: fail loudly unless the user decides otherwise --------------------- #
    relaxed = plan.value("mito_fail")
    failing = [g for g in qc.gates if g.status == "fail"
               and not (relaxed and g.name == "mito_fraction"
                        and qc.metrics.get("median_pct_mito", 0.0) <= relaxed)]
    plan.add("force", "Proceed past failed QC gates", False,
             ("gates failing: " + "; ".join(f"{g.name} ({g.detail})" for g in failing) +
              ". Default refuses (a whole-dataset failure usually means a bad sample); set True "
              "only if it is expected, e.g. a high-mitochondrial tissue such as heart or kidney")
             if failing else "no QC gate fails", "knowledge base", [False, True])

    # -- cell QC ------------------------------------------------------------------ #
    qc_choice = chosen["cell_qc_thresholds"]
    qc_why, qc_src = why("cell_qc_thresholds")
    if qc_choice == "fixed":
        d = Thresholds()
        plan.add("min_genes_per_cell", _QC_LABELS["min_genes_per_cell"], d.min_genes_per_cell, qc_why, qc_src)
        plan.add("max_pct_mito_cell", _QC_LABELS["max_pct_mito_cell"], 20.0 if n_mito else None,
                 qc_why, qc_src)
    elif qc_choice == "mad_per_batch" and profile.batch_key:
        b = obs[profile.batch_key].astype(str).to_numpy()
        per_batch = {}
        for level in pd_unique(b):
            rows = b == level
            per_batch[level] = {k: v for k, (v, _) in
                                _mad_thresholds(genes[rows], counts[rows], mito[rows], n_mito).items()}
        plan.add("per_batch_bounds", "Per-batch cell bounds", per_batch,
                 qc_why + f"; MAD bounds computed inside each '{profile.batch_key}' level", qc_src)
        plan.add("min_genes_per_cell", _QC_LABELS["min_genes_per_cell"], MIN_GENES_FLOOR,
                 f"global floor {MIN_GENES_FLOOR}; per-batch bounds apply", qc_src)
    else:
        fallback = (f"[cell_qc_thresholds = mad_per_batch] no batch column was detected; "
                    "using global MAD bounds instead; "
                    if qc_choice == "mad_per_batch" else "")
        for k, (v, r) in _mad_thresholds(genes, counts, mito, n_mito).items():
            plan.add(k, _QC_LABELS[k], v, fallback + r, "data")
    if hyg.ensembl_without_mito(diag.id_verdict, n_mito):      # issue #149 item 1: say why mito QC is off
        mito_key = "max_pct_mito_cell" if "max_pct_mito_cell" in plan.decisions else "per_batch_bounds"
        plan.decisions[mito_key].reason += "; " + hyg.ensembl_mito_note(
            diag.id_verdict, "the per-cell mito filter is off")
    v = max(Thresholds().min_cells_per_gene, int(math.ceil(0.001 * profile.n_cells)))
    plan.add("min_cells_per_gene", "Min cells expressing a gene", v,
             f"genes seen in fewer than max(3, 0.1% of {profile.n_cells} cells) = {v} cells", "data")

    # -- doublets, genes ------------------------------------------------------------ #
    plan.add("drop_doublets", "Remove doublets (scrublet)", chosen["doublets"] == "scrublet",
             *why("doublets"), [True, False])
    plan.add("exclude_qc_genes", "Drop mitochondrial/ribosomal genes from features",
             chosen["gene_filter"] == "min_cells_drop_qc_genes", *why("gene_filter"), [True, False])

    # -- strategy (normalization x batch x sparsity) ------------------------------- #
    strategy, s_reason, notes = resolve_strategy(chosen)
    scope = chosen["context_scope"]
    separates_batch_into_contexts = (
        scope in ("per_cell_type", "selected_type")
        and profile.cell_type_key is not None
        and profile.batch_key == profile.cell_type_key
        and bool((obs[profile.cell_type_key].dropna().astype(str).value_counts()
                  >= MIN_CELLS_PER_CONTEXT).any())
    )
    if separates_batch_into_contexts and strategy in _BATCH_STRATEGY.values():
        # Each context contains exactly one value of this column, so ComBat or
        # per-batch centring cannot be defined after context subsetting.  Keep
        # the requested KB choice in ``plan.choices`` for provenance, but make
        # the executable strategy truthful and explain the compatibility
        # fallback in the plan and its warnings.
        requested_batch = chosen["batch"]
        effective = {**chosen, "batch": "none"}
        strategy, s_reason, fallback_notes = resolve_strategy(effective)
        fallback = (
            f"batch option '{requested_batch}' is not applied: batch key "
            f"'{profile.batch_key}' is also the cell-type context key, so every "
            "context contains only one batch"
        )
        s_reason += f"; {fallback}"
        notes = [*fallback_notes, fallback]
    plan.add("strategy", "Normalization / aggregation", strategy,
             f"{s_reason}. " + " ".join(why(k)[0] for k in ("normalization", "batch", "sparsity")),
             "knowledge base", _available_strategies())
    plan.warnings += notes
    cols = [None] + [str(c) for c in obs.columns if not str(c).startswith("_")]
    plan.add("batch_key", "Batch / donor column", profile.batch_key,
             f"'{profile.batch_key}' detected" if profile.batch_key else "no batch-like column",
             _source(ctx, "batch_key"), cols)
    plan.add("group_key", "Pseudobulk group column", None, "only used by 'pseudobulk'", "default", cols)
    if strategy == "metacells":
        plan.add("cells_per_metacell", "Cells per metacell", 20,
                 "~20 cells per metacell (Bilous et al. 2024: graining 10-50)", "knowledge base",
                 [10, 20, 30, 50])

    # -- features, cells ------------------------------------------------------------ #
    requested_hvg = _N_HVG[chosen["hvg"]]
    n_hvg = requested_hvg
    hvg_reason, hvg_source = why("hvg")
    if requested_hvg is not None and profile.n_genes <= 1.5 * requested_hvg:
        n_hvg = None
        hvg_reason += (f"; the matrix has {profile.n_genes:,} genes (<= 1.5 x "
                       f"{requested_hvg:,}), so all genes are kept instead")
    plan.add("n_hvg", "Highly variable genes (TFs always kept)", n_hvg,
             hvg_reason, hvg_source,
             [500, 1000, 2000, 5000, None])
    sub = chosen["subsample"]
    sub_reason, sub_source = why("subsample")
    subsample_mode = None
    if sub == "random":
        subsample_mode = "random"
    elif sub == "stratified":
        if profile.cell_type_key:
            subsample_mode = "stratified"
        else:
            subsample_mode = "random"
            sub_reason += "; no cell-type column was detected, so uniform random subsampling runs instead"
    plan.add("subsample_cells", "Subsample cells", None if sub == "none" else SUBSAMPLE_TO,
             sub_reason, sub_source)
    if subsample_mode is not None:
        plan.add("subsample_mode", "Subsample mode", subsample_mode,
                 sub_reason, sub_source, ["random", "stratified"])
    plan.add("stratify_key", "Stratify subsample / metacell purity by", profile.cell_type_key,
             f"cell-type column '{profile.cell_type_key}'" if profile.cell_type_key
             else "no cell-type column", _source(ctx, "cell_type_key"), cols)
    if diag.ingest.get("fmt") == "table":
        plan.add("orientation", "Table orientation",
                 "genes_x_cells" if diag.ingest.get("transposed") else "cells_x_genes",
                 "detected at load time; stated explicitly so the run is reproducible", "data",
                 ["cells_x_genes", "genes_x_cells"])

    # -- regulators ----------------------------------------------------------------- #
    organism = profile.organism or "human"
    if profile.organism not in (None, "human", "mouse"):
        plan.warnings.append(f"Organism is {profile.organism}: no curated TF list or reference network "
                             "for this species ships with grnbench - provide a TF list (human "
                             "orthologs are an approximation, regulatory edges are not always conserved).")
    tf_source = _TF_SOURCE[chosen["regulators"]].format(organism=organism)
    tf_reason, tf_reason_source = why("regulators")
    declared = organism_note(ctx, adata)      # issue #149: the organism is guessed, uns says another
    if declared:
        tf_reason += "; " + declared
    plan.add("tf_source", "Transcription-factor list", tf_source, tf_reason, tf_reason_source,
             ["file", "dorothea:human", "dorothea:mouse", "none"])

    # -- cellular context -------------------------------------------------------------- #
    plan.contexts = [{"name": "all_cells", "subset_key": None, "subset_value": None}]
    context_note = ""
    if scope in ("per_cell_type", "selected_type") and profile.cell_type_key:
        labels = obs[profile.cell_type_key]
        sizes = labels.dropna().astype(str).value_counts()      # a missing label is not a cell type
        if labels.isna().any():
            plan.warnings.append(f"{int(labels.isna().sum()):,} cells have no label in '{profile.cell_type_key}': "
                                 "they are in no per-cell-type context (the pooled context keeps them)")
        big = [str(k) for k, n in sizes.items() if n >= MIN_CELLS_PER_CONTEXT]
        if scope == "selected_type":
            big = big[:1]
            if big:
                context_note = (f"; selected the largest eligible type '{big[0]}' "
                                f"({int(sizes.loc[big[0]]):,} cells)")
        elif len(big) > MAX_CONTEXTS:
            plan.warnings.append(f"{len(big)} cell types have >= {MIN_CELLS_PER_CONTEXT} cells; the "
                                 f"{MAX_CONTEXTS} largest are inferred separately (pick others in the "
                                 "interface or with --subset-key/--subset-value)")
            big = big[:MAX_CONTEXTS]
            context_note = f"; running the {MAX_CONTEXTS} largest eligible cell types"
        if big:
            plan.contexts = [{"name": f"{profile.cell_type_key}={k}",
                              "subset_key": profile.cell_type_key, "subset_value": k} for k in big]
            # issue #149 item 3: "unknown" & co. are not NaN, so dropna() keeps them; named here, not excluded
            missing = hyg.missing_label_mask(pd.Series(big, dtype=object))
            for k, is_missing in zip(big, missing):
                if is_missing:
                    context_note += (f"; the group '{k}' ({int(sizes.loc[k]):,} cells) is a missing-label "
                                     "value, but it is inferred like a cell type")
        else:
            context_note = (f"; no cell type has at least {MIN_CELLS_PER_CONTEXT} cells, "
                            "so the pooled all-cells context runs instead")
    elif scope in ("per_cell_type", "selected_type"):
        context_note = "; no cell-type column was detected, so the pooled all-cells context runs instead"
    context_reason, context_source = why("context_scope")
    plan.add("context_scope", "Cellular context", scope,
             context_reason + context_note, context_source,
             ["pooled", "per_cell_type", "selected_type"])
    if strategy == "metacells":
        sizes = ([int((obs[c["subset_key"]].astype(str) == str(c["subset_value"])).sum())
                  for c in plan.contexts] if plan.contexts[0]["subset_key"] else [profile.n_cells])
        n_meta = min(sizes) // plan.value("cells_per_metacell")
        if n_meta < 100:
            plan.warnings.append(
                f"metacells: the smallest context has {min(sizes)} cells -> ~{n_meta} metacells of "
                f"{plan.value('cells_per_metacell')}; fewer than ~100 samples makes gene-gene "
                "statistics unstable (Bilous et al. 2024). Use single cells, pooled context or "
                "smaller metacells.")

    # -- flags carried to the user ------------------------------------------------------ #
    if diag.knee_no_cliff:          # issue #116 (never on the registered datasets: filtered, or a knee call)
        plan.warnings.append(f"The matrix looks unfiltered ({adata.n_obs:,} barcodes, most with < 100 counts) "
                             "but the barcode-rank knee found no clear cell/empty-droplet cliff: every barcode "
                             "was kept, so the profile, the QC and the cell bounds describe droplets, not cells.")
    for r in selection.results.values():
        if r.blocked_by_dependency:
            best = next(o for o in r.options if o.option == r.best_overall)
            plan.warnings.append(f"{r.decision}: best option '{best.name}' is unavailable "
                                 f"({best.unavailable_reason}); using '{chosen[r.decision]}'.")
    for key, votes in ctx.conflicts().items():
        plan.warnings.append(f"Sources disagree on {key}: {votes} - using '{ctx.get(key)}'.")
    if ctx.get("value_type") == "normalized":
        plan.warnings.append("A text source describes the values as normalized (TPM/log...). The matrix "
                             "passed the raw-count check; verify you loaded the count matrix.")
    plan.warnings += diag.ingest.get("warnings", [])  # integer but non-raw-looking counts (#42)
    for g in failing:
        plan.warnings.append(f"QC gate '{g.name}' fails: {g.detail}")
    plan.warnings += profile.notes
    plan.methods = recommend_methods(n_cells=min(profile.n_cells, SUBSAMPLE_TO),
                                     n_features=n_hvg or profile.n_genes, methods=methods)
    return plan


def pd_unique(values) -> list:
    seen, out = set(), []
    for v in values:
        if v not in seen:
            seen.add(v)
            out.append(v)
    return out


def _available_strategies() -> list[str]:
    from grnbench.strategies import list_strategies
    return [s["name"] for s in list_strategies() if s["available"]]


def safe_folder_name(name: str, limit: int = 80) -> str:
    """A folder name valid on every OS (cell-type labels contain ':', '/', '|', ...).

    It never ends in '.', not even where the length cut falls: Windows drops a trailing dot, so 'X.' and 'X'
    would be one folder there.
    """
    out = "".join(c if c.isalnum() or c in "-_." else ("-" if c == "=" else "_") for c in name)
    return out.strip("._")[:limit].rstrip(".") or "context"


def _folder_key(name: str) -> str:
    """When two folder names are one folder on disk: Windows and macOS ignore case, macOS also Unicode
    normalization (a precomposed and a decomposed spelling of one letter)."""
    import unicodedata
    return unicodedata.normalize("NFD", unicodedata.normalize("NFD", name).casefold())


def unique_folder_names(names, limit: int = 80) -> list[str]:
    """``safe_folder_name`` of each label, made unique.

    'CD4+T', 'CD4 T' and 'CD4:T' are three cell types but one safe name; the later ones get a short hash of
    their own label so no context overwrites, or is refused because of, another. A label can already be
    spelled like that hashed name, so every candidate is checked: then a longer hash, then a counter.
    Names are compared as the file systems compare them (``_folder_key``): 'T' and 't' are one folder.
    """
    import hashlib
    import itertools
    taken: set[str] = set()
    out = []
    for raw in names:
        sha = hashlib.sha1(raw.encode("utf-8")).hexdigest()
        tags = itertools.chain([sha[:6], sha[:12], sha], (f"{sha[:12]}-{i}" for i in itertools.count(2)))
        name = safe_folder_name(raw, limit)
        while _folder_key(name) in taken:
            tag = next(tags)
            name = f"{safe_folder_name(raw, limit - len(tag) - 1)}-{tag}"
        taken.add(_folder_key(name))
        out.append(name)
    return out


def execute_plan(plan: PrepPlan, input_path, out_dir, *, tf_list: list[str] | None = None,
                 gene_mapping_path=None, seed: int = DEFAULT_SEED,
                 force: bool = False) -> list[dict]:
    """Run ``run_prep`` once per cellular context of ``plan``.

    Writes ``plan.json`` (and ``per_batch_bounds.json`` when used) into
    ``out_dir`` and each context into ``out_dir/<context>``. Returns one record
    per context: ``{"context", "out_dir", "processed", "seed", "cli", "error"}`` - a
    context that fails (e.g. a QC gate) is reported, not raised, so the other
    contexts still run.
    """
    from grnbench.pipeline import PrepError, run_prep

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    plan.save(out / "plan.json")
    bounds_path = None
    if plan.value("per_batch_bounds"):
        bounds_path = out / "per_batch_bounds.json"
        bounds_path.write_text(json.dumps(plan.value("per_batch_bounds"), indent=2), encoding="utf-8")
    tf_path = None
    if tf_list:
        tf_path = out / "tf_list.txt"
        tf_path.write_text("\n".join(tf_list) + "\n", encoding="utf-8")
    records = []
    forced = bool(force or plan.value("force", False))      # what runs is what the command records
    folders = unique_folder_names([str(c["name"]) for c in plan.contexts])
    for ctx, name in zip(plan.contexts, folders):
        ctx_dir = out / name
        kwargs = plan.run_prep_kwargs(tf_list, gene_mapping_path=gene_mapping_path)
        kwargs.update(subset_key=ctx.get("subset_key"), subset_value=ctx.get("subset_value"),
                      seed=seed, force=forced)
        rec = {"context": ctx["name"], "out_dir": str(ctx_dir), "processed": None, "error": None, "seed": seed,
               "cli": ["grnbench", "prep"] + plan.cli_args(input_path, ctx_dir, tf_path, context=ctx,
                                                           bounds_path=bounds_path,
                                                           gene_mapping_path=gene_mapping_path,
                                                           seed=seed, force=forced)}
        try:
            res = run_prep(input_path, ctx_dir, **kwargs)
            rec["processed"] = str(res.processed_path)
            rec["manifest"] = str(res.manifest_path)
        except (PrepError, ValueError, KeyError) as exc:
            rec["error"] = f"{type(exc).__name__}: {exc}"
        records.append(rec)
    (out / "contexts.json").write_text(json.dumps(records, indent=2), encoding="utf-8")
    return records


def _fmt(v) -> str:
    if isinstance(v, float):
        return str(int(v)) if v.is_integer() else repr(v)
    return str(v)


def resolve_tf_list(source: str | None, file_path=None) -> list[str] | None:
    """TF names for a plan's ``tf_source`` ("dorothea:<organism>", "file", "none")."""
    from grnbench.tfs import fetch_dorothea, load_tf_list

    if source in (None, "none"):
        return None
    if source == "file":
        if file_path is None:
            raise ValueError("tf_source is 'file' but no TF list file was given")
        return load_tf_list(file_path)
    if source.startswith("dorothea:"):
        reg = fetch_dorothea(organism=source.split(":", 1)[1])
        return sorted(set(map(str, reg["source"])))
    raise ValueError(f"unknown tf_source {source!r}")


def _source(ctx: DatasetContext, key: str) -> str:
    facts = sorted(ctx.evidence(key), key=lambda f: -f.weight)
    return facts[0].source if facts else "default"


# Cheap baselines always; tree ensembles when installed (best BEELINE
# performers); partial correlation only when cells outnumber features.
_BASELINES = {"spearman_signed", "pearson_signed", "clr"}
_TREES = {"grnboost2", "genie3"}


def recommend_methods(n_cells: int, n_features: int, methods: list[dict] | None = None) -> list[dict]:
    if methods is None:
        from grnbench.infer import list_methods
        methods = list_methods()
    out = []
    for m in methods:
        name = m["name"]
        rec, why = False, m.get("description", "")
        if not m["available"]:
            why = f"unavailable: {m['unavailable_reason']}"
        elif name in _BASELINES:
            rec, why = True, "fast baseline (seconds); every method should beat it"
        elif name in _TREES:
            cost = n_cells * n_features
            rec = name == "grnboost2" or cost <= 5_000_000
            why = ("tree ensemble, top performer in BEELINE" if rec else
                   f"{n_cells} cells x {n_features} genes: GENIE3 would take hours - use GRNBoost2")
        elif name == "ppcor":
            rec = n_cells > n_features + 50
            why = ("partial correlation: removes indirect edges" if rec else
                   f"partial correlation needs more cells ({n_cells}) than genes ({n_features})")
        out.append({"name": name, "available": bool(m["available"]), "recommended": bool(rec),
                    "reason": why, "origin": m.get("origin", "")})
    return out
