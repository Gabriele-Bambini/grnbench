---
status: mapped
mapped_at: 2026-09-23T16:37:00+02:00
mapped_by: haiku
confidence: high
---
# Map - grnbench

## One-liner
Plugin-first GRN benchmark with bulletproof preprocessing, interchangeable inference, null-calibrated evaluation, and transparent comparison reporting.

## How to run

**One command — full pipeline:**
```bash
grnbench run --data d.h5ad --out runs/today --tf-list tfs.txt \
             --gold dorotheaABC --plugin my_model.json --stability
```

**Stage by stage:**
```bash
grnbench prep  data.h5ad --out runs/mydata --tf-list tfs.txt --n-hvg 3000
grnbench infer runs/mydata/processed.h5ad --out runs/infer --tf-list tfs.txt
grnbench eval  runs/infer/*.edges.tsv.gz --gold gold_pairs.tsv --out runs/eval --stability-data runs/mydata/processed.h5ad
grnbench methods       # list available/unavailable methods
```

**Python API:** `from grnbench import run_prep`; `from grnbench.infer import run_method`; `from grnbench.eval import evaluate_ranking`; `from grnbench.html_report import render_html`

## Architecture

Four-stage pipeline: **load-verify → preprocess → infer → evaluate → report**

1. **Ingest** (any format) → counts-verified AnnData with hyperlinked provenance
2. **QC gates** (fail loudly on cell/gene thresholds, mito%, ribo%, batch effects)
3. **Normalization** (raw/CPM/Pearson/log2; auto-recommended) + HVG ∪ TF feature selection
4. **Inference** (built-in: Pearson/MI/CLR/GENIE3/GRNBoost2/ARACNe3; user plugins via entry point or file ref)
5. **Evaluation** (null-calibrated metrics: AUPRC, AUROC, precision@k; regulator/degree-preserving randomization; subsample stability; gold reachability)
6. **Reporting** (self-contained HTML: heatmaps, PR curves, ID coverage, null diagnostics, ranking contradictions; JSON manifest with reproducibility metadata)

**Key design commitments:**
- Transparent null calibration: every metric includes a `(hits+1)/(n+1)` p-value against fixed-universe randomization
- Resolved comparisons only: paired t-test with Holm correction over independent experimental units (datasets, not edges/cells)
- Ordering instability reported: method rankings differ across gold-standard references; harness says so
- Counts-verified throughout: `.X` = raw counts, `.layers["counts"]` preserved across preprocessing
- Validated against BEELINE: `benchmarks/beeline_backtest.py` replicates BEELINE protocol; pass condition is regime match (AUPRC-ratio ~1.10, EPR ~1.20)

## Modules

### src - Core grnbench pipeline + adapters (6646 LOC, 34 files)

**Role**: Four-stage GRN benchmarking: (1) preprocess raw counts → normalized + features, (2) infer edges via plugins, (3) evaluate against gold standard, (4) report.

**Entry points**: CLI dispatch (`cli_main` → `cli_*.py`); programmatic API via `grnbench.pipeline:run_prep`, `infer:run_method`, `eval:evaluate_ranking`, `html_report:render_html`.

**Stage 1 — Preprocessing** (`pipeline.py`, `ingest.py`, `qc.py`, `hygiene.py`, `strategies.py`, `features.py`, `tfs.py`, `recommend.py`)
- `load_any()` — unified multi-format ingestion (h5ad, 10x h5/MTX, loom, CSV/TSV), counts integrity verified, gene-name hygiene logged
- `run_qc()` — loud-fail gates on cell/gene thresholds, mito%, ribo%, batch effects; `QCReport` tracks what failed
- `run_prep()` — full pipeline: load → QC → filter genes/cells → apply strategy (raw/CPM/Pearson/log-transform) → TF-union feature selection → outputs (h5ad + manifest + report)
- `Thresholds`, `Gate`, `QCReport` — QC configuration and results
- Strategies registry: list and apply preprocessing methods via `get_strategy()`, `apply_strategy()`
- Hygiene module: normalize labels (case, synonyms), missing-label masking, organism hints

**Stage 2 — Inference adapters** (`infer/`)
- Base contract (`base.py`): `MethodSpec`, `prepare_adata()`, `run_method()`, normalized edge schema `[source, target, weight]`
- Built-in methods (`classical.py`): Pearson, mutual info (MI), CLR; optional tree methods (GENIE3/GRNBoost2 via arboreto)
- Correlation baseline (`correlation.py`): signed-absolute Spearman
- External methods: ARACNe3 (`aracne3.py` — jar executor with timeout/availability probing)
- Plugin system (`plugins.py`): discover Python objects or external commands via entry points; edge-alias normalization (`source`→`tf`/`regulator`, `target`→`target_gene`, `weight`→`score`/`mi`)

**Stage 3 — Evaluation** (`eval/`)
- Gold standards (`gold.py`): load/validate gold edge sets; `coverage_report()` shows identifier overlap with predictions
- Metrics (`metrics.py`): `GoldIndex` (canonicalize + rank edge predictions), `average_precision`, `precision_at_k`, `pr_curve`, `auc_roc`, `top_k_overlap`, `label_edges()`; identifies genes missing from gold
- Nulls (`nulls.py`): diagnostic randomization (random, regulator-stratified, degree-preserving) to establish empirical significance
- Stability (`stability.py`): subsample-based stability of top-k edges across resampled cells
- Universe (`universe.py`): strict candidate-pair validation, alignment to fixed TF/target sets
- Significance (`significance.py`): gold thinning sensitivity (how ranking breaks when gold is incomplete), reference ordering diagnostics
- Comparison module (`comparison.py`): paired-t / Holm statistical tests across independent experimental units

**Stage 4 — Reports** (`html_report.py`, `report.py`)
- `render_html()` — self-contained comparison report: heatmaps, bars, PR curves, ID coverage, null diagnostics, ranking contradictions
- `build_manifest()`, `render_report()` — JSON provenance (steps, params, env, input fingerprint) + Markdown human summary
- `file_fingerprint()`, `sequence_identity()` — reproducibility tracking

**Calibration pilot** (`calibration.py`)
- `calibration_gate()` — Monte Carlo binomial bounds on false-superiority risk across paired t-tests
- Used by `cli_calibrate` for score-level validation

**Cross-stage dependencies**: All stages work on counts-backed AnnData (`.X` = raw counts, `.layers["counts"]` preserved); TF list resolved at ingest time; eval normalizes edges to canonical schema; gold standards are gene-name case-insensitive.

**Key symbols**: `run_prep()`, `run_method()`, `evaluate_ranking()`, `render_html()`, `MethodSpec`, `GoldIndex`, `QCReport`, `PrepResult`, `IngestReport`

**Depends on**: scanpy, AnnData, numpy, pandas, scipy (core); arboreto (tree methods); dask/distributed (large datasets); click (CLI)

### tests - Verification suite for preprocessing, inference, evaluation, and CLI workflows (4023 LOC, 21 files)

- `conftest.py` - Shared test fixtures: synthetic_adata (400 cells, known properties), extreme_mito_adata (high mitochondrial load), cascade_adata (planted TF→G cascade); depends on anndata, numpy, scipy
- `test_calibration.py` - Calibration gate contracts; binomial bounds, power/FPR pilots, config validation; depends on `grnbench.calibration`
- `test_cli.py` - Umbrella command, `grnbench run`, eval with stability; end-to-end workflows; depends on all CLI entry points
- `test_cli_prep.py` - Preprocessing CLI argument marshalling (thresholds, gates, strategy selection); depends on `grnbench.cli`
- `test_classical.py` - Classical methods (Pearson, MI, CLR, correlation); baseline recovery and cascade detection; depends on `grnbench.infer.classical`
- `test_comparison.py` - Reference inference on independent experimental units; Holm adjustment, paired-t; depends on `grnbench.eval.comparison`
- `test_eval.py` - Evaluation metrics (AP, AUROC, precision@k), null hypothesis testing (random/degree-preserving), gold standard loading; depends on `grnbench.eval.*`
- `test_eval_cli_provenance.py` - Replay-safe evaluation contracts; regression invariants from dated audits; depends on `grnbench.cli_eval`
- `test_eval_contract.py` - Evaluation contract invariants (truncation, invalid scores, edge validation); depends on `grnbench.eval`
- `test_features_pipeline.py` - TF-union rule for HVG selection; end-to-end pipeline from raw counts; depends on `grnbench.features`, `grnbench.pipeline`
- `test_hygiene.py` - Data-hygiene utilities (label normalization, ID matching, organism detection, stratification); depends on `grnbench.hygiene`
- `test_infer.py` - Inference adapter contracts; list_methods availability reporting, cascade recovery; depends on `grnbench.infer.*`
- `test_inference_reliability.py` - Regressions: MI/CLR invariance, TF identifier failures; depends on `grnbench.infer.classical`
- `test_ingest.py` - Format loading (h5ad, CSV roundtrip), orientation detection, counts integrity; depends on `grnbench.ingest`
- `test_plugins.py` - Plugin contract (Python objects, entry points, external commands); isolation fixture; depends on `grnbench.plugins`
- `test_prep_gaps.py`, `test_provenance_reliable.py`, `test_qc.py`, `test_report_html.py`, `test_significance.py`, `test_strategies.py` - Specialized tests for QC gates, report rendering, significance annotation, and preprocessing strategies

**Cross-test dependencies:** Tests depend primarily on src/ (core modules).

### audit - Regression fixtures and validation against published references (562 LOC, 3 files) (from scan)

- `audit/2026-09-13/implementation/beeline-reference-generateExpInputs.py` - Gene filtering, selection, feature parity with BEELINE: validates preprocessing alignment with canonical reference implementation
- `audit/2026-09-13/reproduce_eval_findings.py` - Read-only audit fixtures; regression invariants from dated scientific audits; not for repair/regression to prior behavior
- `audit/2026-09-13/reproduce_mi_findings.py` - Invariant: MI/CLR must not infer dependence from row independence; source-level validation

**Role:** Replay-safe regression fixtures and reproducibility audits; depends on `grnbench.eval`; referenced by test suite for invariant verification.

### benchmarks - Protocol validation and performance characterization (424 LOC, 2 files) (from scan)

- `benchmarks/beeline_backtest.py` - BEELINE-derived selection and grnbench fixed-universe comparison protocol; regime check (AUPRC-ratio ~1.10, EPR ~1.20 on cell-type-specific ChIP-seq); 83 of 84 cell-type results within 3× regime (125/126 with tree methods); found 8 defects in grnbench (each with regression test in CHANGELOG.md)
- `benchmarks/plot_calibration.py` - Plot recorded calibration JSON without rerunning experiments

**Role:** Validate harness against BEELINE reference; publish regime stability; depends on `grnbench.eval`, `grnbench.plugins`; outputs to `benchmarks/results/`.

### hpc - Environment and data provisioning scripts (105 LOC, 2 files) (from scan)

- `hpc/fetch_beeline.sh` - Download BEELINE benchmark datasets
- `hpc/setup_env.sh` - Environment setup for HPC runs

**Role:** Provisioning and setup for large-scale benchmarking runs; shell-level utilities.

### examples - Demonstration plugins and models (72 LOC, 1 file) (from scan)

- `examples/plugins/ridge_model.py` - Example Python plugin: multivariate ridge regression of each gene on TF expression; demonstrates plugin contract (`name` attribute, `infer(adata, tf_list, **params)` method returning `[source, target, weight]` DataFrame)

**Role:** Plugin template and proof-of-concept; no external dependencies beyond numpy/scipy.

## Data & artifacts

- **Input**: Counts matrices (h5ad, 10x HDF5, Loom, CSV/TSV); gold-standard edge lists (TSV); TF name lists (TXT)
- **Output**: Preprocessed AnnData (`processed.h5ad`), inferred edges (TSV.GZ), evaluation reports (HTML/JSON/MD), run manifest (JSON)
- **Results dirs (not scanned)**: `OUT/`, `benchmarks/results/`
- **Version**: v0.7.0 (ranking retraction noted in CHANGELOG.md; prior versions reported pseudo-replication; 0.7.0 uses paired t-test over experimental units)

## Conventions & gotchas

- Counts integrity: preprocessing preserves `.X` = raw counts; normalized data in `.obsm` layers or explicit `.layers["counts"]`
- TF list resolution: TFs matched case-insensitively at ingest time; unmatchable names logged; used by feature selection and all inference methods
- Gold standard matching: case-insensitive; identifier overlap reported in coverage_report
- Plugin discovery: lazy on first `list_methods()` call; entry-point methods auto-discovered; file refs require explicit path
- Edge schema: canonical `[source, target, weight]`; plugins alias (`tf`→`source`, `target_gene`→`target`, `score`/`mi`→`weight`)
- Comparison units: paired t-test uses *datasets* as independent units, not edges or cells; unresolved pairs not called "equivalent"
- Hard rules: QC gates fail loud; preprocessing is not optional; null metrics required in evaluation reports

## Open questions

- (unverified) Plugin performance bottlenecks: Dask parallelization scope and tuning for 10M+ edges
- (unverified) Cross-species gold-standard coverage: how many references per organism are needed for stable rankings
- (unverified) Hyperparameter sensitivity: methods like GENIE3/GRNBoost2 tuning regime and defaults in arboreto wrapper
