"""Rule-based pipeline selector: knowledge base x dataset profile -> pipeline.

For every decision point in :mod:`grnbench.kb` each option gets a score from
explicit rules over the :class:`grnbench.profile.Profile`. Every rule that
fires leaves a reason (with its references) on the option, so the interface
can show why an option won or lost.

* The chosen option is the best **available** one (dependencies installed AND
  a grnbench executor, ``kb.DecisionPoint.availability``): an option grnbench
  cannot run is never chosen, and forcing one fails loudly in the plan. When a
  better option is unavailable, the selection says so ("best: scDblFinder - needs R").
* A decision is **uncertain** when the runner-up is within ``UNCERTAIN_MARGIN``
  of the winner. Uncertain high-impact decisions are where alternative
  pipelines are generated and tested on the data (``candidates``).
* Candidates change one decision at a time relative to the best pipeline, so
  each alternative's pros and cons stay attributable to one choice.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field

from grnbench.kb import DECISION_INDEX, DECISIONS, KB_VERSION
from grnbench.profile import MIN_CELLS_PER_CONTEXT, Profile

UNCERTAIN_MARGIN = 1.0
DEFAULT_BONUS = 0.5
IMPACT_ORDER = {"high": 0, "medium": 1, "low": 2}
DROPLET_PLATFORMS = {"10x Chromium", "Drop-seq", "inDrop", "Seq-Well", "SPLiT-seq", "sci-RNA-seq"}


@dataclass
class ScoredOption:
    option: str
    name: str
    score: float
    available: bool
    unavailable_reason: str
    reasons: list = field(default_factory=list)   # [{"delta", "text", "refs"}]


@dataclass
class DecisionResult:
    decision: str
    question: str
    impact: str
    chosen: str
    best_overall: str
    uncertain: bool
    runner_up: str | None
    options: list[ScoredOption]

    @property
    def blocked_by_dependency(self) -> bool:
        return self.best_overall != self.chosen


@dataclass
class Selection:
    kb_version: str
    results: dict[str, DecisionResult]

    def chosen(self) -> dict[str, str]:
        return {k: r.chosen for k, r in self.results.items()}

    def uncertain(self) -> list[DecisionResult]:
        return sorted((r for r in self.results.values() if r.uncertain and r.runner_up),
                      key=lambda r: IMPACT_ORDER[r.impact])

    def candidates(self, max_candidates: int = 4, always: tuple = ()) -> list[dict]:
        """Best pipeline first, then one-change alternatives.

        ``always`` lists decisions whose runner-up is always tried (model-dependent
        choices such as the cellular context: per-cell-type networks help
        correlation methods and starve mutual-information ones); the remaining
        slots go to uncertain decisions, high impact first.
        """
        base = self.chosen()
        out = [{"name": "recommended", "changes": {}, "choices": dict(base)}]
        forced = [self.results[d] for d in always if d in self.results and self.results[d].runner_up]
        rest = [r for r in self.uncertain() if r.decision not in always]
        for r in forced + rest:
            if len(out) >= max_candidates:
                break
            alt = dict(base)
            alt[r.decision] = r.runner_up
            out.append({"name": f"{r.decision}={r.runner_up}",
                        "changes": {r.decision: (r.chosen, r.runner_up)}, "choices": alt})
        return out

    def to_dict(self) -> dict:
        return {"kb_version": self.kb_version,
                "results": {k: {**asdict(r), "blocked_by_dependency": r.blocked_by_dependency}
                            for k, r in self.results.items()}}


class _Scores:
    def __init__(self, decision_id: str):
        self.d = DECISION_INDEX[decision_id]
        self.opts = {}
        for o in self.d.options:
            ok, why = self.d.availability(o.id)
            self.opts[o.id] = ScoredOption(o.id, o.name, 0.0, ok, why)
        self.add(self.d.default, DEFAULT_BONUS, "knowledge-base default")

    def add(self, option: str, delta: float, text: str, refs=()) -> None:
        so = self.opts[option]
        so.score += delta
        so.reasons.append({"delta": delta, "text": text, "refs": list(refs)})

    def result(self) -> DecisionResult:
        ranked = sorted(self.opts.values(), key=lambda s: -s.score)
        avail = [s for s in ranked if s.available]
        chosen = avail[0] if avail else ranked[0]
        runner = avail[1] if len(avail) > 1 else None
        uncertain = runner is not None and chosen.score - runner.score < UNCERTAIN_MARGIN
        return DecisionResult(self.d.id, self.d.question, self.d.impact, chosen.option,
                              ranked[0].option, uncertain, runner.option if runner else None, ranked)


def _droplet(p: Profile) -> bool:
    return p.platform in DROPLET_PLATFORMS or (p.technology == "umi" and p.platform is None)


def _identifier_mapping(p, s):
    if p.id_verdict in ("ensembl_only", "mixed"):
        s.add("ensembl_to_symbol", 3, f"gene identifiers are {p.id_verdict}: TF lists and gold standards "
              "use symbols", ["heumos2023"])
        s.add("keep", -2, "Ensembl IDs would match no TF")
    else:
        s.add("keep", 2, "identifiers already look like symbols")


def _cell_calling(p, s):
    if p.raw_droplets_suspected:
        s.add("emptydrops", 3, f"{p.n_cells} barcodes, most with < 100 counts: unfiltered droplet matrix",
              ["emptydrops2019"])
        s.add("knee", 2, "unfiltered droplet matrix; dependency-free fallback", ["emptydrops2019"])
        s.add("as_provided", -3, "would treat empty droplets as cells")
    else:
        s.add("as_provided", 3, "barcode counts look already filtered to cells")


def _ambient(p, s):
    if p.technology == "full_length":
        s.add("none", 3, "plate-based protocol: no droplet soup", ["ambient2026"])
    elif _droplet(p):
        s.add("decontx", 2, "droplet data; DecontX works without the raw matrix", ["ambient2026", "decontx2020"])
        s.add("soupx", 1.5 if p.raw_droplets_suspected else 0.5,
              "least over-correction; needs empty droplets from the raw matrix", ["ambient2026", "soupx2020"])
        s.add("cellbender", 1.5 if p.raw_droplets_suspected else 0.5,
              "most precise background estimate; needs raw matrix and GPU", ["janssen2023"])
        s.add("none", 1, "no GRN-specific evidence that correction helps; avoids over-correction",
              ["ambient2026"])
    else:
        s.add("none", 1, "protocol unknown: no correction by default")


def _cell_qc(p, s):
    if p.n_cells < 300:
        s.add("fixed", 2, f"only {p.n_cells} cells: distributions too thin for MAD", ["luecken2019"])
    if p.n_batches >= 2 and p.batch_depth_ratio > 1.5 and p.min_cells_per_batch >= 100:
        s.add("mad_per_batch", 3, f"batch medians differ {p.batch_depth_ratio:.1f}x: a global threshold "
              "would cut the shallow batch", ["amezquita2020", "heumos2023"])
    if p.median_pct_mito > 10:
        s.add("miqc", 2, f"median {p.median_pct_mito:.1f}% mitochondrial: high-mito tissue", ["miqc2021", "osorio2021"])
    s.add("mad", 2, "adaptive outlier bounds per dataset", ["scater2017", "pipecomp2020"])


def _doublets(p, s):
    if p.technology == "full_length":
        s.add("none", 3, "plate-based protocol: one cell per well")
    elif _droplet(p) and p.n_cells >= 1000:
        s.add("scdblfinder", 2.5, "fewest false positives on image-verified doublets", ["imagedoubler2024", "scdblfinder2021"])
        s.add("scrublet", 2, "droplet data with enough cells; per-batch calling", ["scrublet2019", "xi2021"])
        s.add("doubletfinder", 1.5, "best accuracy in Xi & Li 2021", ["xi2021"])
        s.add("solo", 1, "deep classifier; more false positives on verified data", ["imagedoubler2024"])
    elif _droplet(p):
        s.add("none", 2, f"droplet data but only {p.n_cells} cells: thresholds unreliable")
    else:
        s.add("none", 1, "protocol unknown: removal off (enable for droplet data)")
        s.add("scrublet", 0.5, "cheap to try if the data are droplet-based")


def _gene_filter(p, s):
    s.add("min_cells", 1, "standard minimum-detection filter", ["luecken2019"])
    if p.n_mito_genes:
        s.add("min_cells_drop_qc_genes", 1.5, f"{p.n_mito_genes} mitochondrial genes track cell quality, "
              "not regulation")


def _context(p, s):
    if p.time_course:
        s.add("pooled", 2, "time course: one continuous process, pooled cells follow the trajectory",
              ["beeline2020"])
    if p.n_cell_types_with_enough_cells >= 2:
        s.add("per_cell_type", 3, f"{p.n_cell_types_with_enough_cells} cell types with >= "
              f"{MIN_CELLS_PER_CONTEXT} cells: pooled co-expression would reflect composition",
              ["farahbod2020", "retana2026"])
        s.add("pooled", -1, "cell-type mixing is a known failure mode of GRN inference", ["retana2026"])
    elif p.cell_type_key is None:
        s.add("pooled", 1.5, "no cell-type annotation: annotate first to infer per type")
    else:
        s.add("pooled", 1.5, "cell types too small to infer separately")


def _normalization(p, s):
    if p.technology == "full_length":
        s.add("cpm_log1p", 3, "read counts: shifted log; residual methods model UMI noise", ["ahlmann2023"])
        s.add("pearson_residuals", -2, "assumes UMI (Poisson/NB) counts", ["lause2021"])
        s.add("sctransform", -2, "assumes UMI counts", ["hafemeister2019"])
    else:
        s.add("cpm_log1p", 2, "as good as more complex transforms; the input BEELINE used",
              ["ahlmann2023", "beeline2020"])
        if p.library_cv > 1.0 or p.batch_depth_ratio > 2:
            s.add("pearson_residuals", 2, f"very variable depth (CV {p.library_cv:.2f}): variance-stabilizing "
                  "residuals", ["lause2021"])
        s.add("scran", 1, "robust pooled size factors", ["lun2016", "pipecomp2020"])
        s.add("sctransform", 0.5, "regularized NB; weaker with heavy dropout", ["hafemeister2019"])


def _batch(p, s):
    if p.n_batches < 2:
        s.add("none", 5, "single batch")
        return
    v = p.batch_celltype_cramers_v
    if v is not None and v > 0.7:
        s.add("none", 3, f"batch and cell type are confounded (Cramer's V {v:.2f}): correction would "
              "remove biology", ["heumos2023"])
        s.add("combat", -1, "removes biology aligned with batch", ["cobra2024"])
    else:
        strong = p.batch_silhouette is not None and p.batch_silhouette > 0.1
        s.add("per_batch_consensus", 2 + strong, "batch-reproducible edges without correction assumptions")
        s.add("center_per_batch", 1.5 + strong, "removes between-batch shifts in co-expression", ["cobra2024"])
        s.add("combat", 1 + 0.5 * strong, "location/scale per gene; covariance effects persist", ["cobra2024", "tran2020"])
        s.add("none", 1.5 - 2 * strong, "batches mix well" if not strong else
              f"batches separate on PCA (silhouette {p.batch_silhouette:.2f})")
    s.add("scvi_denoised", 0.5, "model-denoised expression can inflate correlations", ["scvi2018"])
    s.add("harmony", -10, "embedding-only: no corrected gene matrix for GRN inference", ["harmony2019"])


def _confounders(p, s):
    s.add("none", 2, "regression would also remove real regulation (e.g. cell-cycle TFs)", ["luecken2019"])


def _sparsity(p, s):
    s.add("imputation", -3, "did not improve network inference and could be detrimental",
          ["mccalla2023", "andrews2018"])
    if p.sparsity > 0.95 and p.n_cells >= 2000:
        s.add("metacells", 3, f"{p.sparsity:.0%} zeros with {p.n_cells} cells: aggregate similar cells",
              ["bilous2024", "seacells2023"])
    else:
        s.add("single_cells", 2, f"{p.sparsity:.0%} zeros: single cells carry enough signal", ["beeline2020"])


def _hvg(p, s):
    if p.n_genes < 1500:
        s.add("all_genes", 3, f"only {p.n_genes} genes")
    elif p.n_cells < 500:
        s.add("tf_plus_500", 2, f"{p.n_cells} cells: smaller target space", ["beeline2020"])
    else:
        s.add("tf_plus_1000", 2, "BEELINE protocol", ["beeline2020"])
        if p.n_cells >= 2000:
            s.add("tf_plus_2000", 1.5, "enough cells for a broader target space")


def _regulators(p, s):
    if p.tf_requested:
        s.add("user_file", 4, f"TF list given: {p.tf_detected}/{p.tf_requested} detected in the data")
    elif p.organism in ("human", "mouse"):
        s.add("dorothea", 2, f"organism {p.organism}: curated TFs", ["dorothea2019"])
        s.add("collectri", 1.5, "larger curated collection", ["collectri2023"])
    else:
        s.add("all_genes", 0.5, "no TF list and organism unknown")


def _subsample(p, s):
    if p.n_cells > 20_000:
        if p.cell_type_key:
            s.add("stratified", 3, f"{p.n_cells} cells: keep every cell type's share")
        else:
            s.add("random", 2, f"{p.n_cells} cells: tractable runtime")
    else:
        s.add("none", 3, f"{p.n_cells} cells: no subsampling needed")


def _pseudotime(p, s):
    if p.time_course:
        s.add("slingshot", 1.5, "lineage-aware pseudotime used by BEELINE", ["beeline2020"])
        s.add("dpt", 1.5, "diffusion pseudotime", ["dpt2016"])
    else:
        s.add("none", 3, "no time course: static methods")


RULES = {
    "identifier_mapping": _identifier_mapping, "cell_calling": _cell_calling,
    "ambient_rna": _ambient, "cell_qc_thresholds": _cell_qc, "doublets": _doublets,
    "gene_filter": _gene_filter, "context_scope": _context, "normalization": _normalization,
    "batch": _batch, "confounders": _confounders, "sparsity": _sparsity, "hvg": _hvg,
    "regulators": _regulators, "subsample": _subsample, "pseudotime": _pseudotime,
}


def select(profile: Profile) -> Selection:
    results = {}
    for d in DECISIONS:
        s = _Scores(d.id)
        RULES[d.id](profile, s)
        results[d.id] = s.result()
    return Selection(KB_VERSION, results)
