"""Knowledge base: how to clean single-cell data for GRN inference.

Structure: pipeline **stages** -> **decision points** (the places where a
choice changes the result) -> **options** (every known solution). Each option
carries when to use / avoid it, pros and cons, what it does to regulatory
inference specifically, the level of evidence behind that advice, references,
and how grnbench implements it.

Evidence levels say what the advice rests on - most preprocessing literature
optimizes clustering or differential expression, which is not the same task:

* ``grn``       benchmarked for GRN inference itself;
* ``coexpr``    benchmarked for gene-gene association / co-expression;
* ``scrna``     benchmarked for general scRNA-seq tasks (clustering, DE, integration);
* ``consensus`` best-practice guideline or expert consensus;
* ``heuristic`` grnbench heuristic, not yet confirmed by a confirmatory experiment on the grnbench bench.

The selector (:mod:`grnbench.selector`) reads this module; the Markdown
document is generated from it (``render_markdown``), so code and documentation
cannot drift apart. ``KB_VERSION`` is recorded in every run manifest.
"""
from __future__ import annotations

import importlib.util
import shutil
from dataclasses import dataclass, field

KB_VERSION = "2026.10.2"

EVIDENCE_LEVELS = {
    "grn": "benchmarked for GRN inference",
    "coexpr": "benchmarked for gene-gene association / co-expression",
    "scrna": "benchmarked for general scRNA-seq tasks (clustering, DE, integration)",
    "consensus": "best-practice guideline / expert consensus",
    "heuristic": "grnbench heuristic, not yet confirmed by a confirmatory experiment on the grnbench bench",
}


@dataclass(frozen=True)
class Reference:
    key: str
    citation: str
    link: str  # DOI or URL


REFERENCES = {r.key: r for r in [
    Reference("luecken2019", "Luecken & Theis 2019, Mol Syst Biol - current best practices in scRNA-seq analysis", "https://doi.org/10.15252/msb.20188746"),
    Reference("heumos2023", "Heumos et al. 2023, Nat Rev Genet - best practices for single-cell analysis across modalities", "https://doi.org/10.1038/s41576-023-00586-w"),
    Reference("amezquita2020", "Amezquita et al. 2020, Nat Methods - Orchestrating single-cell analysis with Bioconductor (OSCA)", "https://doi.org/10.1038/s41592-019-0654-x"),
    Reference("beeline2020", "Pratapa et al. 2020, Nat Methods - BEELINE: benchmarking GRN inference from scRNA-seq", "https://doi.org/10.1038/s41592-019-0690-6"),
    Reference("mccalla2023", "McCalla et al. 2023, G3 - strengths and weaknesses of network inference from scRNA-seq (imputation does not help)", "https://doi.org/10.1093/g3journal/jkad004"),
    Reference("retana2026", "Fernandez-de-Retana et al. 2026, arXiv - when does GRN inference break (dropout, confounders, cell-type mixing)", "https://arxiv.org/abs/2605.04930"),
    Reference("dazzle2025", "DAZZLE 2025, PLoS Comput Biol - GRN inference with dropout augmentation", "https://doi.org/10.1371/journal.pcbi.1013603"),
    Reference("farahbod2020", "Farahbod & Pavlidis 2020, Genome Res - cellular composition drives coexpression", "https://doi.org/10.1101/gr.256735.119"),
    Reference("skinnider2019", "Skinnider et al. 2019, Nat Methods - measures of association for single-cell transcriptomics", "https://doi.org/10.1038/s41592-019-0372-4"),
    Reference("sccobench2025", "scCoBench 2025, bioRxiv - co-expression benchmark with promoter-reporter lines", "https://doi.org/10.1101/2025.05.26.656221"),
    Reference("ahlmann2023", "Ahlmann-Eltze & Huber 2023, Nat Methods - comparison of transformations for scRNA-seq", "https://doi.org/10.1038/s41592-023-01814-1"),
    Reference("lause2021", "Lause et al. 2021, Genome Biol - analytic Pearson residuals", "https://doi.org/10.1186/s13059-021-02451-7"),
    Reference("hafemeister2019", "Hafemeister & Satija 2019, Genome Biol - sctransform", "https://doi.org/10.1186/s13059-019-1874-1"),
    Reference("lun2016", "Lun et al. 2016, Genome Biol - scran pooled size factors", "https://doi.org/10.1186/s13059-016-0947-7"),
    Reference("emptydrops2019", "Lun et al. 2019, Genome Biol - EmptyDrops", "https://doi.org/10.1186/s13059-019-1662-y"),
    Reference("janssen2023", "Janssen et al. 2023, Genome Biol - background noise and its removal", "https://doi.org/10.1186/s13059-023-02978-x"),
    Reference("ambient2026", "Benchmarking computational decontamination of ambient RNA, 2026, Front Bioinform", "https://doi.org/10.3389/fbinf.2026.1844838"),
    Reference("soupx2020", "Young & Behjati 2020, GigaScience - SoupX", "https://doi.org/10.1093/gigascience/giaa151"),
    Reference("cellbender2023", "Fleming et al. 2023, Nat Methods - CellBender", "https://doi.org/10.1038/s41592-023-01943-7"),
    Reference("decontx2020", "Yang et al. 2020, Genome Biol - DecontX", "https://doi.org/10.1186/s13059-020-1950-6"),
    Reference("scater2017", "McCarthy et al. 2017, Bioinformatics - scater (MAD outlier QC)", "https://doi.org/10.1093/bioinformatics/btw777"),
    Reference("osorio2021", "Osorio & Cai 2021, Bioinformatics - mitochondrial thresholds by species/tissue", "https://doi.org/10.1093/bioinformatics/btaa751"),
    Reference("miqc2021", "Hippen et al. 2021, PLoS Comput Biol - miQC", "https://doi.org/10.1371/journal.pcbi.1009290"),
    Reference("pipecomp2020", "Germain et al. 2020, Genome Biol - pipeComp: filtering and normalization choices", "https://doi.org/10.1186/s13059-020-02136-7"),
    Reference("xi2021", "Xi & Li 2021, Cell Syst - benchmarking doublet detection", "https://doi.org/10.1016/j.cels.2020.11.008"),
    Reference("scdblfinder2021", "Germain et al. 2021, F1000Research - scDblFinder", "https://doi.org/10.12688/f1000research.73600.2"),
    Reference("scrublet2019", "Wolock et al. 2019, Cell Syst - Scrublet", "https://doi.org/10.1016/j.cels.2018.11.005"),
    Reference("doubletfinder2019", "McGinnis et al. 2019, Cell Syst - DoubletFinder", "https://doi.org/10.1016/j.cels.2019.03.003"),
    Reference("solo2020", "Bernstein et al. 2020, Cell Syst - Solo", "https://doi.org/10.1016/j.cels.2020.05.010"),
    Reference("imagedoubler2024", "ImageDoubler 2025, Nat Commun - image-based doublet ground truth", "https://doi.org/10.1038/s41467-024-55434-0"),
    Reference("seacells2023", "Persad et al. 2023, Nat Biotechnol - SEACells", "https://doi.org/10.1038/s41587-023-01716-9"),
    Reference("bilous2024", "Bilous et al. 2024, Mol Syst Biol - building and analyzing metacells", "https://doi.org/10.1038/s44320-024-00045-6"),
    Reference("mcrigor2025", "mcRigor 2025, Nat Commun - rigor of metacell partitioning", "https://doi.org/10.1038/s41467-025-63626-5"),
    Reference("andrews2018", "Andrews & Hemberg 2018, F1000Research - false signals induced by imputation", "https://doi.org/10.12688/f1000research.16613.2"),
    Reference("hou2020", "Hou et al. 2020, Genome Biol - imputation benchmark", "https://doi.org/10.1186/s13059-020-02132-x"),
    Reference("magic2018", "van Dijk et al. 2018, Cell - MAGIC", "https://doi.org/10.1016/j.cell.2018.05.061"),
    Reference("cobra2024", "COBRA 2024, Bioinformatics - batch effects persist in correlation networks", "https://doi.org/10.1093/bioinformatics/btae531"),
    Reference("scib2022", "Luecken et al. 2022, Nat Methods - benchmarking atlas-level integration (scIB)", "https://doi.org/10.1038/s41592-021-01336-8"),
    Reference("tran2020", "Tran et al. 2020, Genome Biol - benchmark of batch-effect correction", "https://doi.org/10.1186/s13059-019-1850-9"),
    Reference("combat2007", "Johnson et al. 2007, Biostatistics - ComBat", "https://doi.org/10.1093/biostatistics/kxj037"),
    Reference("harmony2019", "Korsunsky et al. 2019, Nat Methods - Harmony", "https://doi.org/10.1038/s41592-019-0619-0"),
    Reference("scvi2018", "Lopez et al. 2018, Nat Methods - scVI", "https://doi.org/10.1038/s41592-018-0229-2"),
    Reference("yip2019", "Yip et al. 2019, Brief Bioinform - evaluation of HVG selection", "https://doi.org/10.1093/bib/bby011"),
    Reference("stuart2019", "Stuart et al. 2019, Cell - Seurat v3 (VST HVG)", "https://doi.org/10.1016/j.cell.2019.05.031"),
    Reference("lambert2018", "Lambert et al. 2018, Cell - the human transcription factors", "https://doi.org/10.1016/j.cell.2018.01.029"),
    Reference("dorothea2019", "Garcia-Alonso et al. 2019, Genome Res - DoRothEA", "https://doi.org/10.1101/gr.240663.118"),
    Reference("collectri2023", "Mueller-Dott et al. 2023, Nucleic Acids Res - CollecTRI", "https://doi.org/10.1093/nar/gkad841"),
    Reference("dpt2016", "Haghverdi et al. 2016, Nat Methods - diffusion pseudotime", "https://doi.org/10.1038/nmeth.3971"),
    Reference("slingshot2018", "Street et al. 2018, BMC Genomics - Slingshot", "https://doi.org/10.1186/s12864-018-4772-0"),
    Reference("tirosh2016", "Tirosh et al. 2016, Science - cell-cycle scoring", "https://doi.org/10.1126/science.aad0501"),
    Reference("scenic2017", "Aibar et al. 2017, Nat Methods - SCENIC", "https://doi.org/10.1038/nmeth.4463"),
]}


@dataclass(frozen=True)
class Option:
    id: str
    name: str
    summary: str
    use_when: tuple[str, ...] = ()
    avoid_when: tuple[str, ...] = ()
    pros: tuple[str, ...] = ()
    cons: tuple[str, ...] = ()
    grn_note: str = ""
    evidence: str = "consensus"
    refs: tuple[str, ...] = ()
    implementation: str = "grnbench"   # grnbench | external | planned | not-applicable
    requires: tuple[str, ...] = ()     # "py:module", "r:package", "cli:tool"

    def availability(self) -> tuple[bool, str]:
        """Whether this option can run here, with the reason when it cannot."""
        if self.implementation == "not-applicable":
            return False, "not applicable to GRN inference (see note)"
        if self.implementation == "planned":
            return False, "planned in grnbench, not implemented yet"
        missing = [r for r in self.requires if not _have(r)]
        if missing:
            return False, "missing: " + ", ".join(missing)
        return True, ""


def _have(req: str) -> bool:
    kind, _, name = req.partition(":")
    if kind == "py":
        return importlib.util.find_spec(name) is not None
    if kind == "cli":
        return shutil.which(name) is not None
    if kind == "r":
        if importlib.util.find_spec("rpy2") is None:
            return False
        try:  # pragma: no cover - needs R
            from rpy2.robjects.packages import isinstalled
            return bool(isinstalled(name))
        except Exception:  # noqa: BLE001
            return False
    return False


@dataclass(frozen=True)
class DecisionPoint:
    id: str
    stage: str
    question: str
    why_critical: str
    impact: str               # high | medium | low
    options: tuple[Option, ...]
    default: str
    applies_when: str = "always"

    def option(self, option_id: str) -> Option:
        for o in self.options:
            if o.id == option_id:
                return o
        raise KeyError(f"{self.id}: unknown option {option_id!r}")


@dataclass(frozen=True)
class Stage:
    id: str
    order: int
    name: str
    goal: str
    pitfalls: tuple[str, ...] = field(default_factory=tuple)


STAGES = (
    Stage("ingest", 1, "Load and verify",
          "Get a counts-backed matrix in the right orientation, with provenance.",
          ("normalized/TPM matrices passed off as counts", "genes x cells tables read transposed",
           "duplicated barcodes or gene symbols", "processed data re-processed")),
    Stage("identifiers", 2, "Gene identifiers and organism",
          "Make gene names match the TF list and the gold standards.",
          ("Ensembl IDs vs symbols: TF lists silently match nothing", "version suffixes (ENSG...12)",
           "mouse vs human casing", "duplicated symbols after mapping")),
    Stage("droplets", 3, "Empty droplets and ambient RNA",
          "Keep real cells; remove the cell-free RNA soup every droplet contains.",
          ("raw (unfiltered) droplet matrices loaded as if they were cells",
           "ambient marker genes appearing in every cell type")),
    Stage("cell_qc", 4, "Cell quality control",
          "Remove damaged, empty-like and outlier cells without removing biology.",
          ("fixed thresholds that do not fit the tissue or protocol",
           "high-mito cell types (muscle, kidney, hepatocytes) removed wholesale",
           "one threshold across batches of different depth")),
    Stage("doublets", 5, "Doublets", "Remove droplets holding two cells.",
          ("doublets form fake intermediate states and fake co-expression of two programs",
           "calling on data without cluster structure over-calls")),
    Stage("gene_filter", 6, "Gene filtering", "Drop genes that carry no usable signal.",
          ("filtering away lowly expressed TFs", "mitochondrial/ribosomal genes becoming QC-driven hubs")),
    Stage("context", 7, "Cellular context",
          "Decide which cells one network is inferred from.",
          ("pooled cell types: composition creates correlations that are not regulation",
           "too few cells per context for the chosen method")),
    Stage("normalization", 8, "Normalization and transformation",
          "Remove depth differences without inventing or destroying gene-gene association.",
          ("scaling/z-scoring changes co-expression", "residual methods designed for UMI counts used on reads")),
    Stage("batch", 9, "Batch and confounders",
          "Keep technical and nuisance structure from becoming edges.",
          ("embedding-only integration used as if it corrected genes",
           "batch effects persist in covariance after mean correction",
           "regressing out biology correlated with the confounder")),
    Stage("sparsity", 10, "Sparsity and dropout",
          "Get enough signal per gene pair without smoothing in fake correlation.",
          ("imputation/smoothing inflating correlations", "impure metacells mixing states")),
    Stage("features", 11, "Features and regulators",
          "Choose candidate regulators and target genes (the candidate universe).",
          ("HVG selection dropping TFs", "TF list for the wrong organism", "universe differs between methods")),
    Stage("size", 12, "Cell budget", "Keep inference tractable without biasing the sample.",
          ("random subsampling losing rare states",)),
    Stage("dynamics", 13, "Pseudotime", "Order cells in time for methods that need it.",
          ("pseudotime computed on a branching process as if linear",)),
    Stage("handoff", 14, "Model-specific input", "Give each model the input it was designed for.",
          ("one input for every model: tree methods, discretizing MI methods and count models want different forms",)),
)


def _o(**kw) -> Option:
    for k in ("use_when", "avoid_when", "pros", "cons", "refs", "requires"):
        if k in kw and isinstance(kw[k], str):
            kw[k] = (kw[k],)
        elif k in kw:
            kw[k] = tuple(kw[k])
    return Option(**kw)


DECISIONS = (
    DecisionPoint(
        "identifier_mapping", "identifiers", "How are gene identifiers made to match TF lists and gold standards?",
        "TF lists and gold standards use symbols; a matrix of Ensembl IDs matches none of them and every "
        "regulator silently disappears.", "high",
        (
            _o(id="keep", name="Keep as loaded", summary="Use identifiers as they are.",
               use_when="identifiers are already symbols of the right organism",
               avoid_when="Ensembl IDs or mixed identifiers", evidence="consensus", refs="heumos2023"),
            _o(id="ensembl_to_symbol", name="Map Ensembl IDs to symbols",
               summary="Strip version suffixes, map with a gene annotation service, sum counts of duplicated symbols.",
               use_when="Ensembl or mixed identifiers", pros="TF lists and gold standards match",
               cons=("needs network access (mygene) or an annotation table", "some IDs map to nothing or to several symbols"),
               grn_note="Summing counts of IDs that map to one symbol keeps counts integer.",
               evidence="consensus", refs="heumos2023", implementation="planned", requires="py:mygene"),
            _o(id="orthologs", name="Map to another organism's orthologs",
               summary="Translate mouse symbols to human (or back) to use a TF list / gold of the other species.",
               use_when="no TF list or gold standard exists for the data's organism",
               cons="one-to-many orthologs; regulatory edges are not always conserved",
               evidence="consensus", implementation="planned"),
        ), default="keep"),
    DecisionPoint(
        "cell_calling", "droplets", "Which barcodes are cells?",
        "An unfiltered droplet matrix is mostly empty droplets; treating them as cells floods every "
        "statistic with ambient RNA.", "high",
        (
            _o(id="as_provided", name="Use the provided (filtered) matrix",
               summary="Trust the cell calling already done (e.g. Cell Ranger filtered_feature_bc_matrix).",
               use_when="the matrix is already filtered (typical GEO/CELLxGENE downloads)",
               evidence="consensus", refs="amezquita2020"),
            _o(id="knee", name="Knee / inflection point",
               summary="Rank barcodes by total counts and keep those whose total is at or above the inflection "
                       "(the steepest drop of the curve after the knee). Barcodes under 100 counts are assumed "
                       "empty. Stops with an error when it finds no call (fewer than 200 barcodes reach 100 "
                       "counts, or under 50 or over 90% of barcodes would be called). For raw droplet matrices; "
                       "on a filtered matrix it can cut real cells.",
               use_when="raw droplet matrix and no R available", pros="fast, dependency-free",
               cons=("misses small or low-RNA cells",
                     "the steepest drop can be the wrong cliff (heart_10k_v3 with a floor of 10), "
                     "so check how many cells it calls"),
               evidence="consensus", refs="emptydrops2019"),
            _o(id="emptydrops", name="EmptyDrops",
               summary="Test each barcode's profile against the ambient profile.",
               use_when="raw droplet matrix", pros="keeps low-RNA cells that differ from ambient",
               cons="needs R (DropletUtils)", evidence="scrna", refs="emptydrops2019",
               implementation="external", requires="r:DropletUtils"),
        ), default="as_provided", applies_when="raw droplet matrix"),
    DecisionPoint(
        "ambient_rna", "droplets", "Is cell-free (ambient) RNA removed?",
        "Ambient RNA adds the same background profile to every cell: markers of abundant cell types appear "
        "everywhere, blurring cell-type-specific expression.", "medium",
        (
            _o(id="none", name="No correction", summary="Leave counts as they are.",
               use_when=("plate-based protocols", "no raw matrix available", "low contamination"),
               pros="no risk of over-correction", cons="background inflates low-level expression everywhere",
               evidence="scrna", refs="ambient2026"),
            _o(id="soupx", name="SoupX", summary="Estimate the soup from empty droplets and subtract it per cell.",
               use_when="droplet data with the raw matrix; CPU only",
               pros="least over-correction when contamination is low", cons="needs raw matrix and clusters; R",
               grn_note="No GRN-specific evidence; expected to reduce spurious marker co-occurrence.",
               evidence="scrna", refs=("soupx2020", "ambient2026"), implementation="external", requires="r:SoupX"),
            _o(id="decontx", name="DecontX", summary="Bayesian mixture of native and contaminating counts.",
               use_when="droplet data; works without the raw matrix", cons="can over-correct when contamination is absent",
               evidence="scrna", refs=("decontx2020", "ambient2026"), implementation="external", requires="r:celda"),
            _o(id="cellbender", name="CellBender",
               summary="Deep generative model of ambient RNA and barcode swapping on the raw matrix.",
               use_when="droplet data, raw matrix and a GPU available",
               pros="most precise background estimates (Janssen 2023)",
               cons=("needs GPU and long runtime", "can over-correct when contamination is absent"),
               evidence="scrna", refs=("cellbender2023", "janssen2023", "ambient2026"),
               implementation="external", requires="cli:cellbender"),
        ), default="none"),
    DecisionPoint(
        "cell_qc_thresholds", "cell_qc", "How are low-quality cells defined?",
        "Thresholds decide which cells exist for inference: too loose keeps damaged cells whose "
        "mitochondrial/stress program becomes a hub; too strict removes whole cell types.", "high",
        (
            _o(id="fixed", name="Fixed thresholds",
               summary="Classic cut-offs (e.g. >=200 genes, <=20% mitochondrial).",
               use_when="very small datasets where distributions cannot be estimated",
               cons=("ignores depth, protocol and tissue", "high-mito tissues lose whole cell types"),
               evidence="consensus", refs=("luecken2019", "osorio2021")),
            _o(id="mad", name="MAD outliers (global)",
               summary="Outliers at k median absolute deviations on log scale (genes, counts) and on % mito.",
               use_when="one batch or batches of similar depth",
               pros=("adapts to each dataset", "standard in OSCA / scater / sc-best-practices"),
               cons="assumes most cells are good", evidence="scrna",
               refs=("scater2017", "pipecomp2020", "heumos2023")),
            _o(id="mad_per_batch", name="MAD outliers per batch",
               summary="Same rule computed inside each batch/sample.",
               use_when="several batches or samples with different depth",
               pros="a shallow batch is not removed wholesale", cons="needs >= ~100 cells per batch",
               evidence="scrna", refs=("amezquita2020", "heumos2023"), implementation="grnbench"),
            _o(id="miqc", name="miQC mixture model",
               summary="Joint model of % mito and gene count separates damaged from intact cells.",
               use_when="tissues with naturally high mitochondrial content",
               pros="probabilistic, tissue-adaptive", cons="assumes a two-population structure",
               evidence="scrna", refs="miqc2021", implementation="planned"),
        ), default="mad"),
    DecisionPoint(
        "doublets", "doublets", "Are doublets detected and removed?",
        "A doublet expresses two programs at once, creating co-expression between genes of different cell "
        "types - exactly the signal a GRN method reads as regulation.", "medium",
        (
            _o(id="none", name="No removal", summary="Keep every droplet.",
               use_when=("plate-based protocols", "already doublet-filtered data", "few cells"),
               evidence="consensus", refs="heumos2023"),
            _o(id="scrublet", name="Scrublet",
               summary="Simulate doublets from observed cells; kNN score; run per sample. Batches "
                       "under 100 cells are pooled and scored together; a pool still under 100 "
                       "cells is kept unscored, and the manifest records both.",
               use_when="droplet data with cluster structure", pros="Python, fast",
               cons=("threshold fails without cluster structure (grnbench refuses >20% calls)",
                     "pooled small batches mix samples when simulating doublets"),
               evidence="scrna", refs=("scrublet2019", "xi2021"), requires="py:scanpy"),
            _o(id="scdblfinder", name="scDblFinder",
               summary="Gradient-boosted classifier on simulated doublets.",
               use_when="droplet data; R available",
               pros="top accuracy in its benchmark; fewest false positives in image-verified data",
               cons=("R", "misses part of the doublets on some data"), evidence="scrna",
               refs=("scdblfinder2021", "imagedoubler2024"), implementation="external", requires="r:scDblFinder"),
            _o(id="doubletfinder", name="DoubletFinder",
               summary="Artificial nearest-neighbour doublet proportion.", pros="best accuracy in Xi & Li 2021",
               cons=("R/Seurat", "parameter sweep needed"), evidence="scrna",
               refs=("doubletfinder2019", "xi2021"), implementation="external", requires="r:DoubletFinder"),
            _o(id="solo", name="Solo", summary="Semi-supervised deep classifier on scVI latent space.",
               cons="GPU-friendly deep model; more false positives on image-verified data",
               evidence="scrna", refs=("solo2020", "imagedoubler2024"), implementation="external",
               requires="py:scvi"),
        ), default="none"),
    DecisionPoint(
        "gene_filter", "gene_filter", "Which genes are removed before inference?",
        "Genes seen in a handful of cells give unstable correlations; mitochondrial and ribosomal genes "
        "track cell quality and depth, not regulation.", "low",
        (
            _o(id="min_cells", name="Minimum cells expressing",
               summary="Drop genes detected in fewer than max(3, 0.1% of cells).", evidence="consensus",
               refs="luecken2019"),
            _o(id="min_cells_drop_qc_genes", name="Minimum cells + drop MT-/ribosomal genes as targets",
               summary="Also remove mitochondrial and ribosomal protein genes from the features.",
               pros="prevents QC-driven hubs", cons="loses genuine regulation of those genes (rarely the goal)",
               grn_note="MT-/RP genes co-vary with cell quality and depth across all cells.",
               evidence="heuristic", implementation="grnbench"),
        ), default="min_cells"),
    DecisionPoint(
        "context_scope", "context", "Which cells does one network come from?",
        "In mixed tissue most gene-gene co-expression reflects cell-type composition, not regulation "
        "inside a cell type; cell-type mixing is one of the stressors that break GRN inference.", "high",
        (
            _o(id="pooled", name="All cells pooled", summary="One network from every cell.",
               use_when=("a single cell type or one continuous trajectory", "few cells per type"),
               cons="composition-driven correlations become edges", evidence="coexpr",
               refs=("farahbod2020", "retana2026")),
            _o(id="per_cell_type", name="One network per cell type",
               summary="Split by the cell-type column and infer each context separately.",
               use_when="annotated, distinct cell types with enough cells each (>= ~200)",
               pros="edges reflect regulation within a state", cons="fewer cells per network; needs annotation",
               evidence="coexpr", refs=("farahbod2020", "retana2026"), implementation="grnbench"),
            _o(id="selected_type", name="One selected cell type",
               summary="Restrict the analysis to the cell type of interest.",
               use_when="the question is about one population", evidence="coexpr", refs="farahbod2020",
               implementation="grnbench"),
        ), default="pooled"),
    DecisionPoint(
        "normalization", "normalization", "How are depth differences removed?",
        "Normalization changes gene-gene association directly; scaling can lower co-expression recovery "
        "and residual methods assume UMI count noise.", "high",
        (
            _o(id="cpm_log1p", name="Shifted log (CPM + log1p)",
               summary="Scale to 10k counts per cell, log(x+1).",
               use_when="default for UMI and full-length data",
               pros=("as good as more complex transforms in a large comparison", "the input BEELINE used for all methods"),
               cons=("depth effects can remain in lowly expressed genes", "any scaling step can lower co-expression recovery"),
               grn_note="Keep it unscaled (no z-scoring) before correlation-based methods.",
               evidence="scrna", refs=("ahlmann2023", "beeline2020", "sccobench2025", "skinnider2019")),
            _o(id="pearson_residuals", name="Analytic Pearson residuals",
               summary="Residuals of a Poisson/NB null with fixed overdispersion.",
               use_when="UMI data with very variable depth",
               pros="variance-stable across depth", cons="designed for UMI counts, not read counts; negative values",
               evidence="scrna", refs=("lause2021", "ahlmann2023")),
            _o(id="sctransform", name="sctransform (regularized NB)",
               summary="Regularized negative-binomial regression on depth.",
               cons=("R", "weaker on data with heavy dropout"), evidence="scrna",
               refs=("hafemeister2019",), implementation="external", requires="r:sctransform"),
            _o(id="scran", name="scran pooled size factors + log",
               summary="Size factors pooled across similar cells, then log.",
               pros="robust size factors with many zeros", cons="R", evidence="scrna",
               refs=("lun2016", "pipecomp2020"), implementation="external", requires="r:scran"),
            _o(id="raw_counts", name="Raw counts",
               summary="No transformation, for models that model counts themselves.",
               use_when="the downstream model declares count input", evidence="consensus", refs="scenic2017"),
        ), default="cpm_log1p"),
    DecisionPoint(
        "batch", "batch", "How are batches handled?",
        "Batch effects create shared shifts in many genes - perfect fake co-expression - and they persist "
        "in the covariance even after mean correction.", "high",
        (
            _o(id="none", name="No correction", summary="Leave batch structure in.",
               use_when="a single batch, or batch fully confounded with biology (correction would remove it)",
               evidence="consensus", refs="heumos2023"),
            _o(id="center_per_batch", name="Center genes per batch (regress out batch)",
               summary="Subtract each gene's batch mean on the log scale.",
               use_when="batches share cell-type composition",
               pros="removes location shifts; transparent", cons="does not remove variance/covariance differences",
               grn_note="Removes the between-batch component of correlation.",
               evidence="coexpr", refs="cobra2024", implementation="grnbench"),
            _o(id="combat", name="ComBat", summary="Empirical-Bayes location/scale adjustment per gene.",
               use_when="batches share cell-type composition",
               cons="covariance-level batch effects persist; removes biology aligned with batch",
               evidence="coexpr", refs=("combat2007", "cobra2024", "tran2020")),
            _o(id="per_batch_consensus", name="Infer per batch, combine edges",
               summary="Run inference inside each batch and aggregate the edge ranks.",
               use_when="several batches with enough cells each",
               pros="no correction assumptions; batch-reproducible edges rank higher",
               cons="fewer cells per run; more compute", evidence="heuristic", implementation="planned"),
            _o(id="scvi_denoised", name="scVI denoised expression",
               summary="Batch-conditioned generative model; decoded expression at a reference batch.",
               cons=("model-smoothed expression can inflate correlations like imputation", "GPU recommended"),
               evidence="scrna", refs=("scvi2018", "scib2022"), implementation="external", requires="py:scvi"),
            _o(id="harmony", name="Harmony / embedding-only integration",
               summary="Corrects a low-dimensional embedding, not gene expression.",
               grn_note="Gives no corrected gene matrix, so it cannot feed gene-level GRN inference.",
               evidence="scrna", refs=("harmony2019", "scib2022"), implementation="not-applicable"),
        ), default="none", applies_when="batch column with >= 2 levels"),
    DecisionPoint(
        "confounders", "batch", "Are nuisance signals regressed out?",
        "Regression removes correlation with the confounder - including real regulation of cell-cycle or "
        "stress genes.", "medium",
        (
            _o(id="none", name="No regression", summary="Keep cell cycle, depth and stress signal.",
               use_when="default", pros="does not remove real regulation (cell-cycle TFs are real regulators)",
               evidence="consensus", refs=("luecken2019", "heumos2023")),
            _o(id="cell_cycle", name="Regress cell-cycle scores",
               summary="Score S/G2M (Tirosh genes) and regress out.",
               use_when="cycling is a nuisance for the question asked",
               cons="removes cell-cycle regulation from the network", evidence="consensus",
               refs=("tirosh2016", "luecken2019"), implementation="planned"),
        ), default="none"),
    DecisionPoint(
        "sparsity", "sparsity", "How is sparsity handled?",
        "Dropout and latent confounders erase the advantage of better methods; imputation does not help "
        "network inference and can create false correlations.", "high",
        (
            _o(id="single_cells", name="Single cells", summary="Infer from single-cell profiles as they are.",
               use_when="moderate sparsity and enough cells",
               grn_note="Dropout robustness belongs in the model (e.g. dropout augmentation), not in imputed data.",
               evidence="grn", refs=("beeline2020", "mccalla2023", "dazzle2025")),
            _o(id="metacells", name="Metacells",
               summary="Aggregate counts of small groups of transcriptionally similar cells (kNN partitions).",
               use_when="very sparse data with many cells (>= ~2000)",
               pros="less dropout; better gene-gene correlation", cons=("fewer samples", "impure metacells mix states"),
               grn_note="Check purity against the cell-type column; keep >= ~100 metacells.",
               evidence="coexpr", refs=("bilous2024", "seacells2023", "mcrigor2025"), implementation="grnbench"),
            _o(id="pseudobulk", name="Pseudobulk per group",
               summary="Sum counts per annotated group (sample x cell type).",
               use_when="many groups (tens) available", cons="needs many groups to have enough samples",
               evidence="consensus", refs="heumos2023"),
            _o(id="imputation", name="Imputation / smoothing (MAGIC, kNN)",
               summary="Share information between neighbouring cells to fill zeros.",
               avoid_when="network inference in general",
               cons=("did not improve network inference and could be detrimental", "creates false correlations"),
               evidence="grn", refs=("mccalla2023", "andrews2018", "hou2020", "magic2018"),
               implementation="external", requires="py:magic"),
        ), default="single_cells"),
    DecisionPoint(
        "hvg", "features", "Which target genes enter the network?",
        "Defines the candidate universe; runtime of most methods grows with genes x TFs.", "medium",
        (
            _o(id="tf_plus_1000", name="TFs + top 1000 HVGs",
               summary="Highly variable genes union every TF detected in the data.",
               pros="BEELINE protocol; tractable", evidence="grn", refs=("beeline2020", "yip2019", "stuart2019")),
            _o(id="tf_plus_500", name="TFs + top 500 HVGs", summary="Smaller BEELINE setting.",
               use_when="slow methods or few cells", evidence="grn", refs="beeline2020"),
            _o(id="tf_plus_2000", name="TFs + top 2000 HVGs", summary="Broader target space.",
               use_when="fast methods and >= 2000 cells", evidence="heuristic", refs="beeline2020"),
            _o(id="all_genes", name="All expressed genes", summary="No selection.",
               use_when="small gene space (< ~1500 genes)", cons="runtime and noise grow with genes",
               evidence="consensus"),
        ), default="tf_plus_1000"),
    DecisionPoint(
        "regulators", "features", "Where does the list of candidate regulators come from?",
        "Only listed TFs can be sources; a list for the wrong organism or naming matches nothing.", "high",
        (
            _o(id="user_file", name="User TF list", summary="A curated list supplied with the data (e.g. the human TF census).",
               evidence="consensus", refs="lambert2018"),
            _o(id="dorothea", name="DoRothEA TFs", summary="TFs with curated regulons (confidence levels A-E).",
               pros="regulons double as a silver standard", cons="~ few hundred TFs; circular if also used as gold",
               evidence="consensus", refs="dorothea2019", requires="py:decoupler"),
            _o(id="collectri", name="CollecTRI TFs", summary="Larger curated TF-target collection.",
               cons="circular if also used as gold", evidence="consensus", refs="collectri2023",
               requires="py:decoupler", implementation="planned"),
            _o(id="all_genes", name="Every gene", summary="No regulator list.",
               cons="universe grows quadratically; edges not restricted to TFs", evidence="consensus"),
        ), default="user_file"),
    DecisionPoint(
        "subsample", "size", "Are cells subsampled?",
        "Tree ensembles scale with cells x genes x TFs; subsampling must not drop rare states.", "low",
        (
            _o(id="none", name="All cells", summary="Use every cell.", evidence="consensus"),
            _o(id="stratified", name="Stratified subsample",
               summary="Subsample keeping each cell type's share.", use_when="> 20k cells with annotation",
               evidence="heuristic"),
            _o(id="random", name="Random subsample", summary="Uniform subsample.",
               cons="can lose rare states", evidence="heuristic"),
        ), default="none"),
    DecisionPoint(
        "pseudotime", "dynamics", "Is pseudotime computed?",
        "Dynamic methods (SCODE, SINCERITIES, SINGE, ...) require cells ordered in time.", "medium",
        (
            _o(id="none", name="No pseudotime", summary="Static methods only.", evidence="consensus"),
            _o(id="provided", name="Use provided pseudotime", summary="A column already in the metadata.",
               evidence="consensus", refs="beeline2020"),
            _o(id="dpt", name="Diffusion pseudotime", summary="Scanpy DPT from a root cell.",
               cons="needs a root; assumes connected trajectory", evidence="scrna", refs="dpt2016",
               implementation="planned"),
            _o(id="slingshot", name="Slingshot", summary="Lineage-aware pseudotime (used by BEELINE).",
               cons="R", evidence="grn", refs=("slingshot2018", "beeline2020"), implementation="external",
               requires="r:slingshot"),
        ), default="none", applies_when="dynamic methods selected or time-course data"),
)

DECISION_INDEX = {d.id: d for d in DECISIONS}


def validate_kb() -> list[str]:
    """Internal consistency problems (empty list = consistent)."""
    problems = []
    stages = {s.id for s in STAGES}
    for d in DECISIONS:
        if d.stage not in stages:
            problems.append(f"{d.id}: unknown stage {d.stage}")
        ids = [o.id for o in d.options]
        if len(ids) != len(set(ids)):
            problems.append(f"{d.id}: duplicated option ids")
        if d.default not in ids:
            problems.append(f"{d.id}: default {d.default!r} not an option")
        if d.impact not in {"high", "medium", "low"}:
            problems.append(f"{d.id}: impact {d.impact!r}")
        for o in d.options:
            if o.evidence not in EVIDENCE_LEVELS:
                problems.append(f"{d.id}.{o.id}: evidence {o.evidence!r}")
            problems += [f"{d.id}.{o.id}: unknown reference {r}" for r in o.refs if r not in REFERENCES]
            if o.implementation not in {"grnbench", "external", "planned", "not-applicable"}:
                problems.append(f"{d.id}.{o.id}: implementation {o.implementation!r}")
    return problems


def option_table() -> list[dict]:
    """One row per option with live availability - what the interface shows."""
    rows = []
    for d in DECISIONS:
        for o in d.options:
            ok, why = o.availability()
            rows.append({"stage": d.stage, "decision": d.id, "option": o.id, "name": o.name,
                         "evidence": o.evidence, "implementation": o.implementation,
                         "available": ok, "unavailable_reason": why,
                         "default": o.id == d.default})
    return rows


def render_markdown() -> str:
    """The knowledge base as a readable document (docs/knowledge/preprocessing_kb.md)."""
    out = [f"# grnbench preprocessing knowledge base (v{KB_VERSION})", "",
           "Generated from `src/grnbench/kb.py` - edit the code, not this file.", "",
           "## Evidence levels", ""]
    out += [f"- **{k}**: {v}" for k, v in EVIDENCE_LEVELS.items()]
    out += ["", "## Stages", ""]
    for s in STAGES:
        out.append(f"{s.order}. **{s.name}** - {s.goal}")
        out += [f"   - pitfall: {p}" for p in s.pitfalls]
    for s in STAGES:
        ds = [d for d in DECISIONS if d.stage == s.id]
        if not ds:
            continue
        out += ["", f"## {s.order}. {s.name}"]
        for d in ds:
            out += ["", f"### {d.question}", "",
                    f"*Impact: {d.impact}. Applies: {d.applies_when}. Default: `{d.default}`.*", "",
                    d.why_critical, ""]
            for o in d.options:
                out.append(f"- **{o.name}** (`{o.id}`; evidence: {o.evidence}; {o.implementation}) - {o.summary}")
                for label, items in (("use when", o.use_when), ("avoid when", o.avoid_when),
                                     ("pros", o.pros), ("cons", o.cons)):
                    if items:
                        out.append(f"  - {label}: " + "; ".join(items))
                if o.grn_note:
                    out.append(f"  - GRN: {o.grn_note}")
                if o.refs:
                    out.append("  - refs: " + "; ".join(f"[{r}]({REFERENCES[r].link})" for r in o.refs))
    out += ["", "## References", ""]
    out += [f"- **{r.key}**: {r.citation}. {r.link}" for r in REFERENCES.values()]
    return "\n".join(out) + "\n"
