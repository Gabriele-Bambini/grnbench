# Attempts - grnbench
Append-only. One entry per approach/experiment. Never delete; update status instead.
Format (keep the heading exactly, `pos attempt` parses it):

    ### A-001 | YYYY-MM-DD | agent | STATUS | short title
    - goal: what it should achieve
    - approach: what was done (config, params, commit)
    - result: outcome + metric
    - evidence: log/file/commit path
    - lesson: what to do or avoid next time
    - tags: kw1, kw2

STATUS: RUNNING | SUCCESS | PARTIAL | FAILED | ABANDONED | INCONCLUSIVE


### A-001 | 2026-09-07 | backfill | SUCCESS | Gene case matching failure in BEELINE replication
- goal: Fix silent mismatch between upper-case reference networks and mixed-case data
- approach: Added case-insensitive folding to label_edges, identifier_overlap, and top_k_overlap; new diag module identifier_overlap() returns case_sensitive flag
- result: hESC case mismatch identified; 64k reachable edges now scored instead of NaN
- evidence: eval/metrics.py:238-252; eval/gold.py
- lesson: Always test with real data that has case mismatches (BEELINE ships mESC in MGI case, networks in upper case)
- tags: beeline,case-sensitivity,identifier

### A-002 | 2026-09-07 | backfill | SUCCESS | Stratified subsampling deleted whole strata
- goal: Fix stratified sampling where zero-quota classes were bumped to 1, causing overshot and class deletion
- approach: Replaced one-bump-per-zero with largest-remainder allocation; dropped strata now counted and warned
- result: Majority class now receives correct quota instead of zero cells
- evidence: pipeline.py:subsample_balanced; regression test
- lesson: Stratification must use deterministic allocation (largest-remainder) not ad-hoc bumps; test with intentional imbalance
- tags: preprocessing,stratification,qc-gate

### A-003 | 2026-09-11 | backfill | FAILED | Gold-edge resampling pseudo-replication in significance layer
- goal: Bootstrap gold edges to claim method superiority with confidence intervals
- approach: Resampled gold edges 200x without replacement; re-scored all methods on same draw; reported P(A>B) and CI
- result: Claimed 80% resolved comparisons and ordering contradictions 47%; retracted in v0.7.0
- evidence: eval/significance.py (lines 69-80,134-145 in audit); CHANGELOG 0.7.0 retracted
- lesson: Gold edges are not independent replicates; datasets are the independent unit in BEELINE-style design. Confusing the two inflates false positives dramatically. Never treat measurement units as equivalent to sample units.
- tags: stats,pseudo-replication,significance-testing

### A-004 | 2026-09-07 | backfill | SUCCESS | Ranking metrics sensitive to arbitrary tie-break order
- goal: Fix AUPRC, EPR, precision@k that changed with row order when scores were identical
- approach: Changed to tie-safe definitions: AUPRC at distinct-score boundaries; precision@k as expectation under random tie-breaking; EPR range reported
- result: Six orderings of same scores gave AUPRC 0.319-0.478 (was non-deterministic); now fixed at tie-safe value
- evidence: eval/metrics.py:average_precision, precision_at_k, evaluate_ranking; regression test compares routines
- lesson: Sorting matters; test with coarse/repeated scores; always report tie-break uncertainty range
- tags: metrics,determinism,ties

### A-005 | 2026-09-07 | backfill | SUCCESS | Matrix orientation detection broken on 10x barcodes
- goal: Fix silent transpose of genes×cells matrices when barcodes matched gene-symbol patterns
- approach: Added barcode shape detection (TGGCCAAAATGTGGTG-1 pattern) and mitochondrial/ribosomal marker presence check; refuse to guess with --orientation flag
- result: Median genes/cell now 2204 instead of 331; mito gates no longer disabled
- evidence: ingest.py:detect_orientation; regression test vs real CellRanger-v3 file
- lesson: Never use simple pattern matching on identifiers (barcodes look like gene symbols); detect on multiple independent signals; refuse when ambiguous
- tags: ingest,data-hygiene,10x

### A-006 | 2026-09-15 | backfill | SUCCESS | Paired t-test with Holm correction over datasets
- goal: Resolve method comparisons by treating datasets as independent units, not gold edges
- approach: Paired Student t on declared independent units, Holm correction over both directions of every pair, Bonferroni simultaneous intervals
- result: On BEELINE 7 datasets: 0/216 pairs resolved (was 80% under edge pseudo-replication); identifies that design cannot rank
- evidence: eval/comparison.py:compare_units; benchmarks/beeline_backtest.py:unit_comparison
- lesson: Independent unit identification is the crux. Proper design reveals negative results (no resolved pairs). Honest reporting matters more than a publishable-looking number.
- tags: stats,comparison,resolved-vs-unresolved
- update 2026-09-23: 0.7.0 real run 9 methods x 7 BEELINE datasets: 0/108 EPR, 0/108 AUPRC-ratio resolved; archived benchmarks/results/beeline-2026-09-17-9methods/

### A-007 | 2026-09-13 | backfill | FAILED | Default cutoff k changes during gold subsampling, reversing rankings
- goal: Fix evaluation cutoff that recalculates k=n_pos after subsampling, changing what is measured
- approach: Audit fixture: A EPR=4 (k=100), B EPR=2.5; default resample k→50; resampled CI for B-A=[−6.8,−3.59] (B better); fixed k=100 gives [+0.8,+2.2] (A better)
- result: Default k is not fixed cutoff, causing methodological confusion about estimand identity
- evidence: eval/significance.py:69-80,134-145; audit/2026-09-13/eval_findings.json case:default_cutoff_inversion
- lesson: Cutoff must be frozen before predictions or subsampling. Do not recalculate k from resampled gold; it changes the estimand.
- tags: audit,evaluation,cutoff-invariance
- update 2026-09-23: fixed in 0.7.0: gold_thinning_sensitivity uses fixed k, reported as scenarios/quantiles, never a winner

### A-008 | 2026-09-13 | backfill | FAILED | Degree-preserving null creates self-loops and changes node degrees
- goal: Implement degree-preserving edge randomization as a diagnostic null
- approach: Degree null: sample from uniform edge-swap proposals respecting out-degree; aim to preserve row/col sums
- result: Audit fixture: [A→x,A→y,B→x,B→y] becomes [A→x,A→x,B→y,B→y]; self-loops introduced; node degree 2→1 for simple graph
- evidence: eval/nulls.py:26-38; audit/2026-09-13/eval_findings.json case:degree_null_duplicates
- lesson: Preserving multiset row/col sums is not the same as preserving simple-graph node degrees. Self-loops are forbidden in benchmarks but allowed in output. Never use degree null as statistical p-value; it is diagnostic only for topology mismatch.
- tags: audit,null-models,degree-preservation
- update 2026-09-23: 0.7.0 CHANGELOG: degree mode unchanged - verify whether still open
- update 2026-09-23: verified 2026-09-23: fixed in code (eval/nulls.py degree_preserving_nulls: simple double-edge swaps, rejects loops/duplicates, immobile graph raises 'degenerate'; degree mode never yields p-value). Both audit inputs (degree_null_duplicates, null_support_and_random_margins) now raise degenerate. Regressions: tests/test_eval_contract.py, tests/test_audit_fixtures.py. The 0.7.0 CHANGELOG line 'degree mode is unchanged' refers only to the null_metrics speedup, not to this bug (misread during backfill). Closed.

### A-009 | 2026-09-13 | backfill | FAILED | Methods evaluated on separate candidate universes, not fixed universe
- goal: Ensure all methods are scored against the same fixed TF×gene universe
- approach: Metrics currently use only edges present in method output; no penalty for missing predictions
- result: Method A outputs 100 edges with 5 gold hits (AP=.51); truncate to 1 edge, AP=1 (omitted hit disappears from denominator). Not a comparison over same universe.
- evidence: eval/metrics.py:326-340; cli_eval.py:275-293; audit evaluation
- lesson: Fix: Score all methods on explicit fixed candidate universe (TF list × target genes); missing predictions get bottom tie. Comparison universe must be identical across methods.
- tags: audit,universe-alignment,candidate-pairs
- update 2026-09-23: fixed in 0.7.0: eval.universe fixed candidate universe, missing predictions get one bottom tie

### A-010 | 2026-09-13 | backfill | FAILED | MI binning causes order-dependent scores on identically-binned genes
- goal: Implement mutual information with robust tie-breaking for constant-valued genes
- approach: MI computed on binned (equal-frequency) expression; binning assigns bin per distinct value; shared cell index breaks ties
- result: Audit fixture: 16×3 matrix all zeros; MI(A,B)=2.079 (expected 0); fixture permutation yields different scores despite identical data
- evidence: infer/classical.py:79-82; reproduce_mi_findings.json
- lesson: Do not include cell/sample indices in tie-breaking logic; they are not part of the statistical model. Either use rank-based binning or deprecate MI on scRNA-seq.
- tags: audit,classical-inference,mutual-information
- update 2026-09-23: fixed in 0.7.0: infer.classical equal-frequency binning assigns bin per distinct value (first rank)

### A-011 | 2026-09-13 | backfill | SUCCESS | Score-level calibration of paired t + Holm comparison (pilot 200 reps, confirmation 1000 reps)
- goal: Show paired t + Holm keeps nominal FWER and interval coverage when units are independent, and flags pseudo-replication
- approach: Monte Carlo over iid global/partial null, margin-boundary and pseudoreplication scenarios; alpha 0.05, margin 0.1; gate = FWER upper <= 0.075 and coverage lower >= 0.925 (configs/calibration-pilot.json 200 reps, seed 13092026; configs/calibration-confirmation.json 1000 reps, seed 13092027)
- result: Pilot: 3 iid scenarios pass, 2 inconclusive (iid_partial_null_n32, iid_margin_boundary_n32). Confirmation (1000 reps): all 5 iid scenarios pass. Pseudoreplication_null_8x4 fails as designed (pilot FWER 0.470, coverage 0.375).
- evidence: audit/2026-09-13/implementation/calibration-pilot/calibration.md; audit/2026-09-13/implementation/calibration-confirmation/calibration.md
- lesson: 200 reps leave Monte Carlo bounds too wide to pass the gate near the boundary; use >=1000. This is score-level validity only - not biological validation or evidence of novelty (stated in the report).
- tags: calibration,monte-carlo,validity,fwer

### A-012 | 2026-09-13 | backfill | SUCCESS | Provenance tracking with SHA256 fingerprints and recursive inventories
- goal: Record reproducible identifiers for all inputs, parameters, versions, and candidate universes
- approach: Full streamed SHA256 for files; deterministic recursive inventories for directories; ordered gene/TF lists hashed; environment versions captured; candidate universe.tsv.gz persisted before inference
- result: Prep/infer/eval/run manifests now include file hashes, TF matches, effective gene universe, seed, versions and edge-list provenance. Backtest persists expression, candidate U, gold, run_method records.
- evidence: report.py, pipeline.py, cli_*.py, benchmarks/beeline_backtest.py; task3-provenance-status.md
- lesson: Provenance must be recorded *before* outputs, not afterward. Candidate universe must be saved before method predictions. Always hash multi-file inputs deterministically (avoid sorting on iteration order).
- tags: provenance,reproducibility,manifest

### A-013 | 2026-09-15 | backfill | FAILED | BEELINE backtest kept every TF instead of Bonferroni-significant TFs
- goal: Replicate BEELINE input selection in beeline_backtest.py
- approach: Kept every TF present in the expression matrix
- result: hESC candidate universe 2.9M pairs instead of 373k; produced the only resolved comparisons ever seen (4/45 AUPRC-ratio), which vanished once selection was fixed
- evidence: CHANGELOG 0.7.0 Fixed; test_subset_agrees_with_pinned_official_fixture
- lesson: Pin external protocol steps to the official script (SHA-256) and test against its output; a gene-selection rule alone can create or erase significance
- tags: beeline,tf-selection,replication,protocol-fidelity

### A-014 | 2026-09-23 | claude | FAILED | Context subsetting before doublet calling
- goal: Per-cell-type networks via run_prep subset
- approach: subset cells to one cell type right after ingest, then QC + scrublet
- result: scrublet called 22-28% doublets inside single demo cell types (no cluster structure); the 20% ceiling refused the run
- evidence: session 2026-09-23 demo run; fix in pipeline.py (subset after drop_doublets, genes re-filtered in context); tests/test_engine.py::test_per_context_execution_calls_doublets_before_subsetting
- lesson: Order matters: cell QC and doublet calling on ALL cells (per batch), then context subsetting, then gene filtering inside the context.
- tags: engine,doublets,context,ordering

### A-015 | 2026-09-23 | claude | SUCCESS | Candidate-pipeline screening on the demo dataset
- goal: Check that rules + data screening pick sensible pipelines
- approach: screen_candidates on demo (1500 cells, 3 types, 2 batches, planted GRN): per-cell-type vs pooled, center_per_batch, metacells, pearson residuals; proxy pearson_signed vs selection half of the planted network
- result: per-cell-type AUPRC ratio 50-64 vs pooled 30.5; center_per_batch removed batch leakage (0.084->0.002) and was picked by the data over the rules' close call; metacells on 1.3k cells (66 metacells) worst (18.6); pearson residuals 48.9 with higher depth leakage
- evidence: session 2026-09-23 screening run; GUI run-20260923-222135
- lesson: Demo is built with cell-type-specific programs, so per-type wins by construction: confirm on the validation bench with varied scenarios before claiming it. Screening can overturn close rule decisions - keep it as the final judge.
- tags: engine,screening,demo,batch,metacells

### A-016 | 2026-09-24 | claude | PARTIAL | Engine mass test v1 (80 synthetic datasets, pre-fix)
- goal: Does the automatic cleaning beat a standard pipeline on data with a known network?
- approach: benchmarks/engine_bench.py: 19 presets x3 + 40 random; arms baseline/rules/auto/auto_noref/oracle; 3 models scored on the evaluation half of the planted network
- result: 72 completed: auto vs baseline win 83%, median +41%, Wilcoxon p=1.6e-11; per model Pearson +78% / Spearman +88% (86% wins), CLR ties (15% wins); rules-only hurt CLR (48->33); auto without reference = baseline (14% wins, p=0.2). 8 errors from engine bugs (raw droplets, dataset-wide mito gate, ranking crash).
- evidence: benchmarks/results/engine-2026-09-24-v1-prefix/REPORT.md
- lesson: Best cleaning is model-dependent (per-cell-type helps correlation, starves MI): choose per model. No-reference screening does not beat baseline - say so. Rerun after fixes (v2).
- tags: engine,mass-test,validation

### A-017 | 2026-09-24 | claude | SUCCESS | Knee cell calling (dependency-free) for unfiltered droplet matrices
- goal: Handle raw 10x matrices without R/EmptyDrops
- approach: hygiene.knee_cell_threshold: log-rank vs log-total, knee = max distance above chord, threshold at steepest descent after it; diagnose profiles called cells; run_prep --cell-calling knee
- result: Real PBMC3k raw (737,280 barcodes): 2,816 cells called vs ~2,700 by Cell Ranger; simulated raw_droplets: 0/22,000 empty barcodes kept
- evidence: session 2026-09-24; ~/grnbench_data/10x_pbmc3k/raw
- lesson: Profile and MAD thresholds must be computed on called cells, never on all barcodes.
- tags: engine,cell-calling,droplets

### A-018 | 2026-09-24 | claude | SUCCESS | decoupler 2.x removed dc.get_dorothea / get_collectri
- goal: Keep DoRothEA/CollecTRI options working
- approach: tfs.decoupler_network supports dc.op.* (2.x) and dc.get_* (1.x); gold and TF fetch use it
- result: CollecTRI human 42,990 / mouse 43,226, DoRothEA ABC human 32,286 / mouse 32,529 edges saved
- evidence: ~/grnbench_data/references
- lesson: Optional dependencies change APIs silently; the old code path only failed when real data needed it.
- tags: dependencies,decoupler,references

### A-019 | 2026-09-24 | claude | FAILED | Overnight mass test run from the Google Drive repo
- goal: Finish 97 synthetic + 40 real datasets overnight
- approach: benches started from <drive>/grnbench with PYTHONPATH=src, multiprocessing spawn workers
- result: G: (Google Drive) became unreachable during the night: 58 workers FileNotFoundError on engine_bench.py, 2 ModuleNotFoundError; 2 datasets hung ~10 h (sleep). 37/97 synthetic, 11/40 real completed
- evidence: scratchpad bench_local2/run.log, bench_real1/run.log (2026-09-24)
- lesson: Run long benches from a frozen local code snapshot (~/grnbench_bench/code) - also prevents mixing code versions mid-run; hold a per-process no-sleep request (SetThreadExecutionState); the runners are resumable, so resume instead of restarting.
- tags: infrastructure,mass-test,drive,sleep

### A-020 | 2026-09-24 | claude | SUCCESS | Real data: high-mitochondrial tissue refused by the dataset mito gate
- goal: Clean CELLxGENE human liver (565 cells) without refusing healthy tissue
- approach: planner relaxes the dataset-level mito gate (warn 35 / fail 60) when context tissue is liver/heart/kidney/muscle (Osorio & Cai 2021); per-cell MAD mito bound unchanged; --mito-warn/--mito-fail CLI flags
- result: liver median mito 29.7% now passes (was PrepError); per-cell bound 76.3%
- evidence: tests/test_engine.py::test_high_mito_tissue_relaxes_the_dataset_gate_only
- lesson: Dataset-level gates need tissue context; fail-loudly stays the default elsewhere.
- tags: engine,qc,mito,real-data

### A-021 | 2026-09-24 | claude | SUCCESS | Real data: non-model species detected as human/mouse
- goal: Report the true species of CELLxGENE datasets
- approach: context: obs['organism'] curated label = strongest fact (weight 12), species_label keeps binomials; planner warns that no TF list/reference ships for that species
- result: Before: microcebus->mouse, macaque/marmoset->human (4/4 wrong). After: curated binomial wins (tests).
- evidence: tests/test_context.py::test_curated_species_beats_gene_casing
- lesson: Gene-symbol casing only separates human-like from mouse-like; curated metadata must outrank it.
- tags: engine,context,organism,real-data
