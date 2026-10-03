# Changelog

## Unreleased (0.8.0-dev) — interface and context-aware cleaning

### Project
- EXP-20261003-LM-1: without a reference the fixed baseline beats the engine's former choice (claim C7);
  D-021 makes the baseline the engine's choice when no reference is given.
- `engine_bench.py` / `real_bench.py`: a worker killed by the OOM killer now ends the run (retryable,
  resubmission resumes) instead of hanging until the Slurm time limit (`pool_results`).
- Engine robustness (found by EXP-20261002-LM-2): screening fails one candidate, never all, when a
  context cannot be screened; strategies may not emit NaN/inf (`PrepError`); Pearson residuals give a
  zero-count gene 0 instead of NaN; ComBat drops one-cell batches with a warning.
- D-020 pins reference networks per registered experiment: historical reproductions keep their
  hashes; new runs prefer versioned upstream releases or use a dated, exact snapshot manifest.
  `data_manifest.py` argument parsing is now order-independent for reliable strict checks.
- EXP-20261001-LD-1: completed 740/740 synthetic and 40/40 real task accounting, generated stress
  curves, real-data visualizations and LaTeX, retained every refusal/failure, and added portable
  task/data manifests. Real auto cleaning usually ties the baseline; the run remains exploratory.
- `real_bench.py report --out <RUN>` falls back to `tasks.json` when `real_tasks.json`
  is absent, so runs planned before D-019 (EXP-20260926-LM-1) regenerate unchanged.
- EXP-20261002-LM-1: the motif prior lift (C4, C5) replicates on 20 never-used datasets (claim C6).
- `benchmarks/candidate_study.py`: every engine candidate, the baseline and one-change variants scored
  on real datasets under controlled degradations (EXP-20261002-LM-2).
- Registered Lorenzo (`@lorenzo-de-santis`, initials `LD`) in the team metadata.
- Canonical repo is now the private GitHub repo (D-016). Rules for people and agents in
  `AGENTS.md`; onboarding prompt in `docs/START_HERE_AGENT_PROMPT.md`.
- Experiment registry `experiments/` with pre-registration (D-017), data provenance manifest
  `experiments/datasets.tsv` (`benchmarks/data_manifest.py`), `paper/` (claims ledger, reporting
  checklist, CRediT log, availability drafts, `paper/make_figures.py`), `CITATION.cff`.
- Pinned environment `requirements/lock-py311.txt`; `hpc/setup_env.sh` installs it. CI runs on
  Ubuntu / Python 3.11 with that lock, on PRs and `main` only.
- HPC: `hpc/hpc.sh` (check/sync/run/submit/status/log/cancel/fetch; per-user `.env.hpc`, never a
  username in the repo) replaces `hpc/sync_to_cluster.sh`; `docs/HPC.md`. `grnbench.sbatch` and
  `beeline_backtest.sbatch` now ask for 8 CPUs (the `stud` QOS cap), not 16.
- `.gitattributes` keeps shell scripts LF.
- Day-one team setup: `grnbench board` (runs in flight = issues labelled `running`, open PRs, STATE
  next steps); session start/end routines in AGENTS.md; short MAP.md and STATE.md.
- `grnbench viz` (`grnbench.viz`, one style file, Okabe-Ito palette): QC and embeddings before/after,
  decisions, score distributions, paired lifts, PR curves, top-k network, TF recovery and error
  analysis, motif lift by TF, calibration, stability, attribution helpers for learned models.
  `real_bench.py run --keep-artifacts` keeps edges, processed data and decisions for it.
- `grnbench report` (`grnbench.texreport`): LaTeX report per experiment and `--paper` draft with a
  bibliography from the KB; numbers only from result files (templates are digit-free, tested);
  `real_bench.py report` now also writes `viz/` and `latex/`. Structure follows the thesis
  pipeline's end-of-run report.
- Shared HPC skill `.claude/skills/grnbench-hpc/SKILL.md` (pointer in `.agents/skills/`), Antigravity
  rule moved to `.agents/rules/`; onboarding guide extended (from zero, how to work, look inside).
- Shared-work workflow: owner merge after an independent agent review when the other person is
  away (labelled as such; can be reopened); experiments are claimed with an `experiment` issue before
  the record, and the board lists claimed experiments; a recipe for stacked PRs after a squash;
  long local runs from a `git worktree`; claim status `reproduced-computational` (agent rerun from
  the record, human reproduction pending).

### Added
- Context-specific synthetic validation (D-019): the simulator writes shared, union and per-cell-type
  truth; selection/evaluation splits happen before context expansion; pooled and split predictions
  are macro-scored against each context. `engine_bench.py --stress-grid` fixes a 740-task design
  across four network-sharing regimes, eight progressive corruptions and mixed worst cases, with
  generated statistical/relative/operational breaking points. Benchmark task lists and
  score-affecting run configurations are hashed; atomic rows have explicit retry/terminal statuses;
  final reports require complete accounting (`--partial` is watermarked). HPC jobs can retain the
  12 pre-declared synthetic artifact sentinels and now propagate `MAX_ATTEMPTS`.
- `data_manifest.py --check --require-all`: optional exact-snapshot validation reports registered
  paths missing from the data root; the default subset-compatible check is unchanged. HPC and
  onboarding docs now distinguish the dynamic setup cache from registered experiment roots.
- `grnbench ui`: Streamlit interface (data -> profile -> pipeline -> compare pipelines ->
  clean -> models incl. user models -> benchmark). `grnbench demo`, `grnbench kb`.
- `grnbench.kb`: preprocessing knowledge base - 14 stages, 15 decision points, 56 options,
  51 references, evidence level per option; `docs/knowledge/preprocessing_kb.md` generated
  from it (a test keeps them in sync).
- `grnbench.context`: dataset context from the matrix, AnnData metadata, README files and
  download pages (GEO text records), rule-based with snippets and conflict reporting.
- `grnbench.profile` / `grnbench.selector`: dataset profile (protocol from count values,
  unfiltered-droplet check, batch x cell-type Cramer's V, batch silhouette) and rule-based
  selection with reasons, uncertainty and one-change candidate pipelines.
- `grnbench.autoprep`: executable plans (MAD cell bounds, per-batch bounds, strategy
  resolution, one run per cellular context) with the exact CLI; `grnbench.screen`: candidate
  pipelines run and scored on the data, reference split into selection/evaluation halves.
- Strategies `metacells` (k-means on PCA within batch, counts summed, purity reported) and
  `cpm_log1p_center_batch`.
- `run_prep`/`grnbench prep`: `--subset-key/--subset-value` (applied after doublet calling on
  all cells), `--exclude-qc-genes`, `--per-batch-bounds`, `--max-genes-per-cell`,
  `--max-counts-per-cell`, `--max-pct-mito`, `--cells-per-metacell`.

### Results
- EXP-20261001-LD-1 (exploratory): 740/740 synthetic tasks completed and 40/40 real tasks reached a
  terminal state (33 completed, seven expected BEELINE refusals). On 23 paired scored real datasets,
  auto usually ties baseline and shows no consistent recovery gain. Two completed datasets with a
  usable reference (`paul15`, `pbmc3k_raw`) produced no model score and remain explicit findings.
- EXP-20260925-LM-2 (confirmatory, pre-registered): the motif prior lifts Pearson's AUPRC ratio by
  +0.182 [+0.138, +0.227] on 20 fresh CELLxGENE datasets (p_Holm 1.4e-7); EPR lift unresolved.
  Reproduced computationally (independent agent rerun from the record, identical numbers);
  human reproduction pending.
- EXP-20260926-LM-1 (confirmatory, pre-registered, first HPC experiment): on the same 20 datasets
  the motif prior lifts the AUPRC ratio for RegDiffusion by +0.184 [+0.132, +0.237] (n = 19,
  p_Holm 1.2e-6) and for GRNBoost2 by +0.130 [+0.108, +0.153] (n = 20, p_Holm 3.6e-10). EPR lift
  resolved for GRNBoost2, not for RegDiffusion. Reproduced computationally (independent agent
  rerun on the HPC, identical numbers); human reproduction pending.
- `paper/tables/T_models_real_EXP-20260925-LM-1.md` regenerated: it still had the 45-character
  truncated dataset ids fixed in #7 (numbers unchanged).

### Added (HPC experiment tooling, EXP-20260926-LM-1)
- `hpc/engine_bench.sbatch`: `MODELS` (models joined by `+`, since `sbatch --export` splits on
  commas) and `KEEP_ARTIFACTS=1` reach `real_bench.py run`; the job copies `GIT_COMMIT` to
  `$OUT/CODE_COMMIT`. `hpc/setup_env.sh`: `WITH_TORCH=1` installs RegDiffusion + CPU torch from the
  pinned line in the lock header. `hpc/fetch_cistarget.sh`: the six cisTarget v10 files.
- `fetch_cellxgene.py --from-manifest`: fetch exactly the datasets of an earlier manifest (no
  listing, no re-selection). `data_manifest.py --check`: verify a data root against
  `experiments/datasets.tsv` (exit 1 on any mismatch or unlisted file; LF/CRLF-only differences
  in text files are reported as `ok-eol`).
- `grnbench report`: several `;`-separated co-primary comparisons in a record each get a verdict.

### Added (VPN tooling for collaborators)
- `hpc/vpn.ps1` (Windows GlobalProtect): `status`, `reconnect` (restart PanGPS, wait up to 90 s,
  no elevation after a one-time grant), `setup` (prints the admin command with your own SID) and
  `grant <SID>` (idempotent `sc sdset` of one query/start/stop entry). `hpc/hpc.sh`: when the login
  node does not resolve or answer, the first error line is the VPN command for your OS. Docs:
  `docs/HPC.md` section 1b per OS (Linux `globalprotect` CLI, macOS none, both cited from Palo
  Alto docs); agents reconnect the VPN themselves and ask only for SSO, the grant, or on macOS.
  Windows ssh-agent: `StartupType Automatic` once, no per-reboot step.

### Fixed
- Candidate screening now carries the diagnosis' raw-droplet flag into its rebuilt profile. Raw
  10x inputs therefore keep knee cell calling instead of reclassifying called cells as already
  filtered and then applying QC to hundreds of thousands of empty barcodes (issue #30).
- Ingest aggregates gene identifiers that collide under grnbench's case-insensitive matching
  contract, retaining the first spelling and recording case-fold merges separately. This prevents
  one duplicate symbol from blocking every TF-resolved pipeline while preserving strict rejection
  for ambiguous matrices passed directly to inference (issue #30).
- Real-benchmark rows retain a portable structured summary of candidate/context failures and report
  scoreless rows with screening failures explicitly instead of hiding the cause behind `completed`.
- `grnbench report`: record keys containing commas or backticks were read as continuation text and
  fenced code was appended to the key above, garbling the Methods section; a note after the
  primary comparison containing `;` was read as a second hypothesis. LaTeX of EXP-20260925-LM-2
  and EXP-20260926-LM-1 regenerated from unchanged rows (numbers identical).
- `hpc/engine_bench.sbatch` imports grnbench from the checkout it runs (`PYTHONPATH`), so
  `CODE_COMMIT` always names the code that ran, not the conda env's editable install.
- RegDiffusion adapter: genes that are constant in the scored cells (possible after the 3,000-cell
  subsample) are dropped before training instead of making RegDiffusion refuse the whole context
  (`1 genes have 0 variance`; 1/20 datasets in EXP-20260926-LM-1, issue #13). They get no edges.
- `grnbench report --paper`: every experiment section defined the same macro names, so a
  second confirmed experiment would redefine the first one's numbers. Macros are now prefixed per
  experiment (test with two experiments).
- `real_bench._model_matrix` truncated dataset ids to 45 characters; two datasets with a common
  prefix collided and `compare_units` refused the table. Full ids now.
- Doublet calling silently fell back to a crude heuristic because scanpy's scrublet needs
  scikit-image, which was not installed (issue #2). scikit-image is now a dependency (pyproject +
  lock). A fallback now raises a `RuntimeWarning`, is recorded in the manifest step
  (`fallback`, `fallback_reason`) and in real_bench rows (`<arm>_doublets`) and reports.
- Per-batch scrublet refused batches of 1-30 cells and the whole dataset fell back to the
  heuristic (4/20 datasets in EXP-20260925-LM-2, issue #6). Batches under 100 cells
  (`hygiene.SCRUBLET_MIN_CELLS`) are now pooled and scored together; a pool still under 100 cells
  is kept unscored; missing batch labels are one level. The manifest step records
  `scrublet_groups`. On those four datasets every cell is now scored by scrublet.
- `test_list_methods_reports_unavailability` assumed java on PATH; aracne3 has two
  documented prerequisites.

## 0.7.0 — the ranking is not real, and neither was the test that said so (2026-09-15)

### Retracted

The 0.6.0 "significance" layer treated **gold edges** as independent samples.
With ~4,500 gold edges per dataset that made 80% of pairwise method comparisons
look "resolved", and 47% of the resolved orderings appeared to contradict each
other across ground truths. Both numbers were pseudo-replication: the only
independent replicates in a BEELINE-style design are the **datasets** (n = 7),
and a paired test that knows this resolves almost nothing (see below). The
0.6.0 entry is kept for the record and marked retracted. Neither number should
be quoted.

### Removed

- `eval.significance.paired_bootstrap`, `significant_tiers`,
  `contradiction_report` and the `grnbench-eval --significance*` flags.
  Removed rather than deprecated: a deprecated function that still prints
  "significant" is a trap, not a courtesy.

### Added

- **`eval.comparison.compare_units`** and `grnbench compare`: paired Student t
  on declared independent units, Holm over both directions of every pair,
  Bonferroni simultaneous intervals. Refuses fewer than three units, never
  reads "unresolved" as "equivalent", and requires the caller to state what a
  unit is. This is the reference procedure, not a new estimator.
- `benchmarks/beeline_backtest.py` now ends with **`unit_comparison`** —
  `compare_units` with datasets as units, one family per ground truth × metric
  (`epr`, `auprc_ratio`) — and **`regime_check`**, the specific-chip scores
  against the published best (AUPRC-ratio 1.10 / EPR 1.20, factor 3,
  descriptive). Both land in `beeline_backtest.json` and `summary.md`;
  `--summarize existing.json` recomputes them on a finished run.
- `eval.universe`: an explicit, fixed candidate universe. Missing predictions
  get one bottom tie and stay in every denominator, so a method that emits
  fewer edges is not rewarded for it.
- `eval.significance.gold_thinning_sensitivity`: what remains of the old
  layer, renamed for what it is — annotation-removal *scenarios* with fixed
  `k`, reported as quantiles, not confidence intervals, and never a winner.
  `reference_ordering_differences` reports order reversals between references
  as differences: ChIP binding, STRING association and lof/gof response are
  different relations, and disagreeing about them is not a contradiction.
- `eval.metrics.metric_from_labels(metric, weights, y, k)`.

### Changed

- **`null_metrics` is ~30× faster** on real edge lists. A permutation of the
  scores over fixed pairs is the same draw as a permutation of the labels over
  the fixed, already-ranked score vector — same distribution, same tie blocks —
  so the `random` and `regulator` modes no longer re-validate 370k identifier
  strings on every draw (6.5 s → 0.2 s per draw on hESC). Pinned by a test that
  compares both routes on a heavily tied ranking. Draws for a given seed differ
  from 0.6.0; `degree` mode is unchanged because it does move the pairs.
- `tfs.resolve_tf_names` raises `TFListError` when a TF list matches nothing,
  instead of silently promoting every gene to regulator.
- `infer.classical` equal-frequency binning assigns a bin per distinct value
  (first rank), so tied expression values are no longer split by cell order.

### What the corrected inference says

9 methods (`pearson_signed`, `spearman_signed`, `mutual_info`, `clr`, `ppcor`,
`lasso`, `random_forest`, `grnboost2`, `genie3`) × 7 BEELINE datasets, pinned
BEELINE selection (significant TFs + top-500 significant genes), datasets as
units, α = 0.05, Holm:

| metric | resolved pairs (3 references × 36) |
|---|---|
| EPR | **0 / 108** |
| AUPRC-ratio | **0 / 108** |

GENIE3 and GRNBoost2 are not distinguishable from Pearson correlation on any
reference. Meanwhile the point-estimate winner changes with the reference on
7 of 7 datasets. Every ordering the point estimates suggest — including every
"winner flips between ground truths" — is within what seven paired
observations cannot distinguish. That is the finding: not that the methods
are equal, but that this design cannot rank them, and a leaderboard that does
is reporting noise. Archived in `benchmarks/results/beeline-2026-09-17-9methods/`.

The only resolved comparisons ever seen (4 / 45 on AUPRC-ratio, `pearson` /
`spearman` > `clr` / `mutual_info`) came from an earlier run that kept every
TF instead of the Bonferroni-significant ones, which BEELINE's own script does
not do. Fixing the selection removed them: the gene-selection rule alone moved
the only "significant" result the design had.

### Fixed

- `benchmarks/beeline_backtest.py` kept every TF present in the expression
  matrix; BEELINE's `generateExpInputs.py` keeps only TFs passing the
  Bonferroni p-value filter. hESC's candidate universe was 2.9 M pairs instead
  of 373 k. The selection is pinned to the official script's SHA-256 and tested
  against its output on a fixture (`test_subset_agrees_with_pinned_official_fixture`).
- `benchmarks/beeline_backtest.py` `regime_check` prints, for every EPR, the
  hit count behind it against the count expected at random. mDC has 756
  reachable ChIP edges in 263 k candidates: `clr`'s EPR of 4.55 there is 8 hits
  against 1.8 expected, and `pearson`'s 0.00 is 0 hits. Both are now labelled
  small-count instead of being read as a leak or a loss.
- **`run_prep` could write a 0-cell `processed.h5ad` and report success.**
  `--strategy pseudobulk --min-group-cells N` above every group's size (or
  `--force` past a QC gate on a gene-poor fixture) filtered every cell out
  *after* the gate check, and nothing looked at the strategy's output shape
  before writing it and exiting 0. Preprocessing now refuses an empty result
  the same way ingest refuses to guess an orientation: loudly, as a
  `PrepError`, before any output is written.
- **A strategy's `warnings.warn()` (e.g. pseudobulk dropping tiny groups) went
  only to stderr.** Nothing that only reads `manifest.json`/`report.md` — the
  documented run artifacts — could see it. Captured into `manifest["warnings"]`
  like every other pipeline warning.
- **`grnbench-prep --min-genes-per-cell 0` / `--min-cells-per-gene 0` silently
  reapplied the default (200 / 3).** `args.x or Thresholds.x` reads an explicit
  `0` as falsy. This CLI entry point had zero direct test coverage (every
  other preprocessing test calls `run_prep` itself, bypassing argument
  marshalling) — added `tests/test_cli_prep.py`.
- **`--group-key` and `--label-key` together audited only the group key.**
  `qc_group = group_key or label_key`: a 500/90/6 label split came back
  "balanced" because the four donors were. Now `class_imbalance` covers the
  labels and `group_imbalance` the group key, each naming its column.
- **Missing `--label-key` labels became a literal `"nan"` stratum** in
  stratified subsampling, with a proportional quota like any real class. They
  are now excluded from stratification, counted in the step and warned about;
  the cells themselves stay in the matrix (only `--group-key` drops rows, and
  the two help texts now say which does what).
- **Duplicated cell barcodes were reported and then carried through.** They
  are made unique with numeric suffixes; the original barcode is kept in
  `obs['barcode_as_loaded']` and the gate still warns (pooled libraries).

### Added — preprocessing

- **`--drop-doublets`**: scanpy's built-in scrublet (per batch when a batch
  column exists), scores kept in `obs['doublet_score']`, removal logged as a
  step. Falls back to the joint-extreme heuristic *with a manifest warning*
  when scrublet cannot run. Never on by default. **Refuses when the call rate
  exceeds `--max-doublet-rate` (0.20)**: on a homogeneous test set scrublet
  called 497 of 596 cells (83%) doublets, which is a threshold landing inside
  the singlet mode, not a doublet population — dropping it would have deleted
  the dataset and reported success.
- **`--strategy cpm_log1p_combat`**: ComBat on the detected/declared batch
  column after `cpm_log1p`. Never auto-recommended: it removes every
  batch-aligned covariance, including biology that was run one batch per
  condition. Verified on a synthetic two-batch set with gene-wise shifts: the
  mean per-gene batch gap drops below 20% of the uncorrected one. Raw counts
  survive in `layers['counts']`.
- `grnbench-prep` prints `[PREP REFUSED]` (no `--force` hint) for refusals
  that `--force` cannot fix — an empty result, a missing batch column, an
  implausible doublet call rate — and `[QC GATE FAILURE]` only for gates.

## 0.6.0 — how much of the ranking is real? (2026-09-11)

### Added

- **`eval.significance` — the layer a leaderboard cannot provide about itself.**
  *Retracted in 0.7.0: gold edges are not independent units; see above.*
  `grnbench-eval --significance` resamples the gold standard (subsampling without
  replacement, since the metrics are defined over a *set* of gold pairs and a
  with-replacement bootstrap would be silently collapsed), re-scores every method
  on the *same* draw so the comparison is paired, and reports:
  - `P(A > B)` and a 95% interval on the paired difference, per pair;
  - `fraction_resolved` — of all the pairwise comparisons a ranking implicitly
    claims, how many the data actually supports;
  - **tiers**: methods the gold standard cannot separate are grouped, instead of
    being printed in an order the reader will believe.

  Resampling the gold rather than the predictions is deliberate: confidence that
  A beats B is limited above all by which true edges happen to be in the
  ChIP-seq network, and those networks are partial, cell-line-specific and biased
  toward well-studied TFs.

  Cost: each method's gold-hit identity and tie-block structure are computed once
  and a draw is then a boolean mask over them, so 200 draws over a million-edge
  list take seconds rather than the tens of minutes a naive re-scoring needs.

- Six inference methods in `infer/classical.py`, numpy/scipy/scikit-learn only:
  `pearson_signed`, `mutual_info`, `clr`, `ppcor`, `lasso`, `random_forest`.
  The registry goes from 4 methods (of which 2 needed arboreto **and** a live
  dask cluster, and 1 needed a Java jar) to 10, of which 9 run in a plain
  `pip install grnbench[infer]` environment. A benchmark whose only tree-based
  method is arboreto-backed cannot run at all on a login node or in a
  reviewer's fresh venv.
  - `random_forest` computes what GENIE3 computes — per-target feature
    importance over the regulator set — with scikit-learn alone. On a 200-gene
    problem it took 24s against arboreto GENIE3's 171s, and the two agree on
    the top-20 edges (Jaccard > 0.4, pinned by a test that skips when arboreto
    is absent).
  - `clr` is Faith et al., *PLoS Biology* 5(1):e8 (2007); `ppcor` uses
    Ledoit-Wolf shrinkage because scRNA-seq is always in the p > n regime where
    R's `ppcor` inverts a singular covariance.
- `tests/test_classical.py` — 29 tests. Beyond the contract checks (canonical
  schema, descending sort, no self-loops, regulator whitelist honoured,
  determinism), each method is tested against **the thing it exists to do**:
  `mutual_info` must rank a U-shaped target that Pearson misses; `clr` must
  demote a broadcast hub that raw MI promotes; `ppcor` must break an indirect
  A->B->C edge that correlation keeps. A method that returns well-formed
  garbage passes a smoke test and fails the user.

### Changed

- The HTML report gained an **"Is the ranking real?"** section, and the README
  now positions grnbench against geneRNIB honestly: if you want a leaderboard,
  use geneRNIB.

### Audit of published harnesses

The tie-order defect fixed in 0.5.1 was checked against the field rather than
assumed to be general. Reading the actual sources:

- **BEELINE** (`BLEval/EarlyPrecision.py`) is **tie-safe**: it expands the
  selection to every edge tied at the cutoff and divides by the count actually
  selected. `BLEval/AUPRC.py` delegates to `precision_recall_curve` + `auc`.
  Its trapezoidal integration differs from average precision by at most ~8% in
  the GRN regime, and in the *conservative* direction — not a defect.
- **BEAR-GRN** (`R/reproduce_early_metrics.R`) does `arrange(desc(score))` then
  `head(n_edges)`. `arrange` is a stable sort, so where the cutoff falls inside
  a tied block the result depends on the input file's row order. Real, but
  measured at ~1.05x under uniformly distributed gold — a caveat, not a scandal.

The conclusion is the honest one: the published harnesses are largely fine, and
grnbench was the one with the severe version of the bug.

## 0.5.1 — what the BEELINE replication found (2026-09-07)

Running `benchmarks/beeline_backtest.py` against the real BEELINE corpus
(Zenodo 3378975, md5-verified) rather than self-authored synthetic data. Every
entry below is a defect the real data exposed and synthetic data had hidden.

### Fixed

- **Gene symbols were compared with exact string equality**, so an expression
  matrix and a ground truth that disagree on case shared *nothing*. BEELINE
  ships mESC expression in MGI case (`Pou5f1`) and every mouse network in upper
  case (`POU5F1`): 64,557 genuinely reachable gold edges scored as **0 hits**,
  every metric returned `NaN`, and not one line was printed. A wrong number
  wearing the costume of an honest null result — the exact failure the package's
  "unavailability is data, not a crash" rule exists to prevent. Six of the seven
  datasets were unaffected only because BEELINE happened to ship them already
  upper-cased, which is why no synthetic test had ever caught it.
  - `eval.metrics.label_edges` / `evaluate_ranking` / `eval.gold.coverage_report`
    now fold case by default; `grnbench-eval --case-sensitive` opts out.
  - New `eval.metrics.identifier_overlap()` classifies every edge-list/gold pair
    as `ok`, `case-mismatch` or `no-overlap`, and `grnbench-eval` prints
    `[GOLD CASE]` / `[GOLD MISMATCH]` with example symbols from **both** sides.
    A gold standard in a different namespace (Ensembl IDs vs symbols) or for the
    wrong organism now says so instead of silently reporting `NaN`.
- `benchmarks/beeline_backtest.py` applied the published cell-type-specific
  ChIP-seq ceiling (AUPRC-ratio 1.10 / EPR 1.20) to **every** ground truth, and
  so flagged six correct mHSC results as "above the published regime". The
  ceiling is now keyed by ground truth, only `specific-chip` is a falsification
  target, and the generic networks are reported with the reason their numbers
  legitimately run higher — including that **STRING contains co-expression
  channels**, which makes scoring a correlation method against it partly
  circular.

### Fixed — reporting defects found auditing the fix above

Chasing the identifier bug through the report surfaced four more ways a number
could come out wrong and *look* fine. Each was reproduced by running the code
before it was touched, and each now has a regression test.

- **A method that matched no gold edge at all was reported as statistically
  significant.** `null_metrics` computed `hits = sum(null_scores >= observed)`,
  and `>= NaN` is elementwise `False` in NumPy, so an undefined observation
  collected zero hits and got the *smallest* p-value the estimator can emit —
  `1/(n+1) = 0.022` — which `eval_report.html` then painted as a green **pass**
  badge beside a row of empty metric cells. The report simultaneously said "this
  method scored nothing" and "this method beat the null at p = 0.022".
  `null_metrics` now reports `observed_undefined: true` and a `NaN` p-value, and
  the HTML renders `n/a` in red.
- **`eval_report.html` never rendered the identifier diagnostic.** The verdict,
  the `--case-sensitive` state and the method-name collisions all existed in
  `eval_report.json` and the Markdown, but the HTML — the file `grnbench run`
  advertises and the one people actually forward — dropped them. A vacuous run
  now opens with a red banner *above* the verdict table, plus a
  **Gene-identifier matching** section; a renamed method says so in its own note.
- **`top_k_overlap` compared edges case-sensitively**, so two byte-identical
  networks differing only in `Pou5f1` vs `POU5F1` scored Jaccard **0.000** under
  a heading asserting that methods sharing no edges are not confirming each
  other. It now folds case like gold matching, with `case_sensitive=True` to opt
  out.
- **The HTML provenance block could not display the preprocessing environment.**
  `report["environment"]` is always set at scoring time, so the `or` fallback to
  the prep manifest was unreachable: a report could state `arboreto: not
  installed` directly beside the `genie3` edge list arboreto had produced on a
  cluster weeks earlier. Both environments now render, each labelled.

### Fixed — ranking metrics were measuring their own sort order

- **AUPRC, precision@k and EPR depended on the arbitrary order of edges that
  scored *identically*.** Both walked the stored row order, so for any method
  emitting coarse or repeated scores (thresholded output, integer scores, a
  large block of zeros) the headline number was partly noise. Six orderings of
  the **same** scores gave AUPRC 0.319–0.478 and EPR 3.0–5.0, while the
  tie-safe AUROC stayed at 0.9167 throughout — proof the ranking was identical
  and only the tie-breaking moved.
  - `average_precision` now reads precision/recall at distinct-score boundaries
    (the standard tie-safe definition; identical to the old result when scores
    are unique), and `pr_curve` draws the same points instead of the sawtooth
    the within-block order happened to trace.
  - `precision_at_k` returns the expectation under random tie-breaking when the
    cut at `k` lands inside a tied block.
  - `evaluate_ranking` additionally reports `tie_block_at_k` and `epr_range`,
    the worst and best EPR reachable by reordering tied edges alone. A wide
    interval means the headline EPR is a coin flip, and now says so.
  - Re-running the BEELINE backtest, this moves AUPRC-ratio by at most 0.027 and
    EPR not at all for `spearman_signed` and `ridge`, which emit continuous
    scores — the defect was real but those two methods never triggered it.

### Fixed — ingest and preprocessing

- **A genes×cells table was silently loaded transposed.** Orientation was
  decided by gene-likeness alone, but a 10x barcode (`TGGCCAAAATGTGGTG-1`)
  matches the gene-symbol pattern perfectly, so on the standard GEO/BEELINE
  layout both axes scored ~1.0, the 0.1 margin never cleared, and a
  4000×600 matrix was ingested as 4000 cells × 600 genes — recording
  `transposed: false` as though that were a determination. Every QC metric was
  then computed on the wrong axis: median genes/cell 331 instead of 2204, and
  `n_mito_genes` **0 instead of 2**, which silently disables the mitochondrial
  gate. Detection now keys on barcode shape and on the presence of mito/ribo
  markers, and **refuses to guess** when neither decides;
  `--orientation genes_x_cells|cells_x_genes` states the layout outright.
- **Stratified subsampling deleted whole strata.** Every zero-quota class was
  bumped to one cell, overshooting the budget, and the `picked[:n_keep]` slice
  that repaired it cut in alphabetical class order. With 200 singleton classes
  and one class holding half the cells, the majority class kept **zero** cells
  while the manifest still read `stratified by 'ct'`. Allocation is now purely
  largest-remainder (that class now gets its proportional 10 of 20), and any
  stratum that cannot be represented is counted in the manifest and raised as a
  run warning.
- **`filter_cells` removal counts did not add up.** Zero-count cells are
  necessarily also low-gene cells, but both were reported as separate figures,
  so "removed 25 cells with < 200 genes, 10 cells with zero counts" described 25
  actual removals. The step now reports the union plus `n_cells_removed`.

### Added

- `eval.metrics.GoldIndex` — the gold standard normalised once and reused. Null
  models re-score the same set `n_nulls` times per method; on a ~10^6-pair
  ChIP-seq network, re-normalising per draw was the entire runtime.
- The backtest records an `id_match` verdict per row and a `ceiling_verdict`
  (`in-published-regime` / `above-published-regime` / `no-published-ceiling`).

## 0.5.0 — plugin release (2026-09-05)

Closes the v0.4 roadmap: custom-model plugin contract, HTML comparison reports,
stability integrated into evaluation.

### Added

- **Plugin contract** (`grnbench.plugins`). Bring your own GRN method three ways,
  all landing in the same registry, edge schema, manifests and evaluation as the
  built-ins:
  - a `grnbench.methods` **entry point** (auto-discovered, no flags);
  - `--plugin path/to/model.py:MyModel` for a one-off script (no packaging);
  - `--plugin spec.json` describing a **Docker/CLI** command that reads
    `expr.parquet` + `tfs.csv` + `params.json` and writes `edges.parquet`
    (`"format": "tsv"` for pyarrow-less environments).
  A plugin cannot silently shadow a built-in name; a plugin whose import,
  binary or image is missing stays listed with its reason instead of vanishing.
- **`grnbench run`** — preprocess + infer + evaluate in one command, with a
  `run_manifest.json` recording every stage, its exit code and why a stage was
  skipped. Also a `grnbench` umbrella command (`run|prep|infer|eval|methods|
  strategies|version`).
- **HTML comparison report** (`grnbench.html_report`) — self-contained, no CDN,
  no JavaScript, no plotting dependency: PR curves and bar charts are inline SVG,
  so it opens on a compute node with no network. Sections: verdict table, PR
  curves, nulls (degenerate draws surfaced, not hidden), cross-method top-k
  Jaccard heatmap, stability, gold reachability, preprocessing QC gates,
  provenance.
- **Stability in `grnbench-eval`** — `--stability-data processed.h5ad` re-runs
  each scored method on cell subsamples; `--stability-{n,frac,top-k,methods,tf-list}`.
- **Cross-method overlap** — `eval.metrics.top_k_overlap`, reported in Markdown
  and as an HTML heatmap.
- `eval.metrics.pr_curve`; `infer.prepare_adata`; `infer.register_method`;
  `MethodSpec.probe` (availability re-evaluated per listing, so a jar appearing
  on disk is picked up without a restart), `MethodSpec.takes`, `MethodSpec.origin`.
- `examples/plugins/` — `ridge_model.py`, a real scikit-learn method usable as a
  plugin, plus a Docker spec template.
- `hpc/` — SLURM job scripts and notes.
- `benchmarks/beeline_backtest.py` — replicates the BEELINE protocol (Pratapa
  et al., Nature Methods 2020) on BEELINE's own data and ground truths, and
  checks the numbers against the published regime. Because the literature
  reports a *ceiling* of AUPRC-ratio 1.10 / EPR 1.20 on real data, this is a
  falsification test: a result far above that means the harness is leaking, not
  that the method is good. Ships with `hpc/fetch_beeline.sh` and
  `hpc/beeline_backtest.sbatch`.
- 48 new tests (115 total).

### Fixed

- **10x `.h5` loader silently split duplicated genes.** The `.h5` path called
  `var_names_make_unique()` before `_clean_gene_names`, renaming CellRanger's own
  duplicated symbols (e.g. `TBCE`) to `TBCE-1` and leaving the aggregation branch
  dead. The same biological data therefore produced a different matrix than the
  MTX path, the manifest recorded `n_duplicates_merged: 0`, and the renamed copy
  was invisible to TF matching — so a duplicated TF could be dropped by HVG
  selection, violating the "HVG never drops TFs" rule. Counts are now aggregated
  by sum on every loader path. Regression-tested against a real CellRanger-v3 file.
- **`pip install '.[infer]'` was unsatisfiable**: the pin was `arboreto>=1.1`, but
  arboreto's latest PyPI release is 0.1.6. Now `>=0.1.6`, with `dask`/`distributed`
  pinned alongside.
- `environment_block()` no longer reads the deprecated `scanpy.__version__` /
  `anndata.__version__` (FutureWarning on every run); it uses
  `importlib.metadata` and also records the arboreto version.
- `register_plugins` is idempotent for the same `--plugin` argument, so
  `grnbench run` passing it to both the inference and the evaluation stage is no
  longer mistaken for a name collision.
- Null models no longer compute a PR curve per draw (`evaluate_ranking(...,
  with_curve=False)`).
- `run_method` rejects a method returning a non-DataFrame or fewer than three
  columns with a clear contract error instead of an opaque `IndexError`.
- ARACNe3 availability is now a `probe` on its spec rather than a hardcoded
  special case in `list_methods()`.

Found by replicating the BEELINE benchmark on real published data
(`benchmarks/beeline_backtest.py`):

- **Every BEELINE ground-truth network was unreadable.** BEELINE ships its
  ChIP-seq / STRING / lof-gof networks as `Gene1,Gene2`, and that convention
  propagated to most published GRN ground truths. The alias table in
  `eval/gold.py` accepted `gene2` as a target but not `gene1` as a source, so
  all ten shipped networks raised "need source/regulator + target columns" —
  grnbench could not read the reference benchmark of its own field. `gene1`,
  `g1`, `g2` and the symbol variants are now accepted, and the error message
  lists what *is* accepted. Regression-tested in the real format.
- `load_gold` read `net.csv.gz` as TSV, because the separator was chosen from
  the `.gz` suffix rather than the name underneath it, silently yielding one
  wide column. The suffix is now resolved with `.gz` stripped.

Found by the pre-release audit:

- **`grnbench-eval` silently dropped an edge list when two inputs derived the
  same method name.** The name came from `path.name.split(".")[0]`, so scoring
  the same method across two runs (`runs/a/genie3.edges.tsv.gz` and
  `runs/b/genie3.edges.tsv.gz`) fed two files in and reported one method, with
  no warning — a provenance lie. Colliding names are now qualified with the
  parent directory (`runB/genie3`), announced on stderr, and recorded in the
  report under `name_collisions`.
- `canonical_edge_frame` resolved columns by iteration order, so a plugin
  emitting BOTH `source` and `regulator` (or `weight` and `score`) could have
  its regulator column chosen by luck. An exact canonical name now always wins
  over an alias.
- `hpc/grnbench.sbatch` split the `PLUGINS` variable on `:` — but every plugin
  reference grnbench accepts contains a colon (`pkg.mod:Class`,
  `file.py:Class`), so each reference was torn in half. The separator is now
  `;`, and the documented example in `hpc/README.md` was itself broken.
- `hpc/grnbench.sbatch` hardcoded `--stability`, an N-fold runtime multiplier
  with no way to switch it off on exactly the machine where runs cost hours.
  Now `STABILITY=0` / `STABILITY_N=...`.
- `pyyaml` was imported for YAML plugin specs but declared in no extra: a
  documented capability with no supported way to install it. Added to
  `[plugins]`, `[dev]` and `[all]`, and covered by a test.
- **A misspelled method name was reported as an unavailable method.**
  `--methods geni3` printed `[SKIP] geni3: unknown` — the same channel that says
  "arboreto not installed" — so a typo looked like a missing dependency. Unknown
  names now fail the run with `[UNKNOWN METHOD]`, a close-match suggestion
  (`Did you mean: genie3?`) and the list of registered methods. A known method
  that simply cannot run here still skips loudly, and a run where *nothing* can
  run now exits 2 instead of reporting success over zero methods.
- **`grnbench-infer` crashed writing its manifest when no method ran.** The
  output directory was only created as a side effect of `save_edges`, so a run
  with zero successful methods died with `FileNotFoundError:
  inference_manifest.json` instead of failing cleanly. The directory is now
  created up front.

### Changed

- `list_methods()` reports an `origin` (`builtin` / `entry_point` / `python` /
  `command`) and triggers entry-point discovery.
- `grnbench-eval` writes `eval_report.html` by default (`--no-html` to skip) and
  accepts `--prep-manifest` to embed the preprocessing QC gates in it.
- `inference_manifest.json` records registered plugins and the reason each
  skipped method was skipped.

## 0.4.0 — data-hygiene release

Preprocessing (`grnbench-prep`), inference adapters (`grnbench-infer`),
evaluation (`grnbench-eval`). Counts-verified ingest, loud QC gates, strategy
registry with auto-recommendation, HVG ∪ TF feature selection, canonical edge
schema, AUPRC/EPR/AUROC with degree-preserving nulls, full provenance manifests.
