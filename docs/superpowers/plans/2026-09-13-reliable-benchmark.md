# Reliable GRN benchmark implementation plan

> **For agentic workers:** use superpowers:subagent-driven-development. Track acceptance in `docs/superpowers/plans/2026-09-13-progress.md`.

**Goal:** repair the audited scientific defects, make benchmark runs replayable, and provide an executable calibration pilot with honest statistical scope.

**Architecture:** a fixed candidate-universe contract separates ranking recovery from missing output; label thinning is descriptive sensitivity. Valid paired inference consumes independent experimental-unit scores. A separate pilot measures that inference procedure on known data-generating processes before biological claims.

**Tech Stack:** Python >=3.10, NumPy/SciPy/pandas, AnnData/Scanpy, pytest; optional arboreto. No new foundation model or hosted service.

**Spec:** `audit/2026-09-13/AUDIT_E_PIANO_PUBBLICAZIONE.md`, especially sections 2, 5 and 6. User explicitly authorized planning and implementation on 13 September 2026.

## Global constraints

- Preserve the dated audit and its original JSON findings. New verification artifacts go in `audit/2026-09-13/implementation/`.
- This local copy has no Git repository. Make a source backup and hashes, work here, do not initialize/push/commit a repository on the user's behalf.
- Do not imply a paper is validated by software tests or a synthetic calibration pilot.
- Do not equate unlabelled edges with certified biological negatives, nonsignificance with equivalence, or different reference modalities with conflicting truths.
- Reuse the external audit runtime, pin actual dependency resolution, record skipped optional checks.
- Before fixes, add meaningful regression tests that fail on current behavior; retain scientific invariants rather than implementation-shaped assertions.

## Task 1 — evaluation contract and diagnostic semantics

**Owner/files:** evaluation worker: `eval/metrics.py`, new `eval/universe.py`, `eval/nulls.py`, `eval/significance.py`, `eval/stability.py`, `eval/__init__.py`, `cli_eval.py`, `html_report.py`, related evaluation/HTML tests.

**Consumes:** edge frames source,target,weight and GoldIndex. **Produces:** explicit fixed-universe alignment, safe evaluation CLI, descriptive thinning and well-labelled null diagnostics.

- [ ] Add regressions for the exact audited EPR inversion, list truncation with explicit U, different candidate sets, null duplicates/diagonal, and misleading tiers.
- [ ] Align score lists on supplied U; fill missing output as a tied bottom block. Without supplied U, comparisons require identical supports. Reject duplicate canonical pairs/colliding mappings and invalid input scores. Record support origin, missing predictions and gold reachability.
- [ ] Freeze k before thinning; validate frac/n_draws. Rename it gold-thinning sensitivity; remove inferential significant/equivalence/contradiction claims from CLI/HTML. Keep old flag/function only as clearly deprecated compatibility paths.
- [ ] Random-ranking null permutes scores on unchanged U. Regulator-stratified score null is a separately named hypothesis. Any retained sparse degree-rewiring generator must reject duplicates/self-loops, preserve simple degrees and report/raise degeneracy; do not use uncalibrated rewiring as default significance.
- [ ] For overlap, include the complete tied boundary rather than arbitrary row order; separate cell resampling from algorithm seed changes.
- [ ] Run evaluation, significance, HTML and affected CLI tests; report commands and API changes to coordinator.

## Task 2 — inference and identifiers

**Owner/files:** inference worker: `infer/classical.py`, `infer/correlation.py`, `infer/base.py` if needed, `tfs.py`, `examples/plugins/ridge_model.py`, targeted inference/plugin tests. Do not modify pipeline/report/CLI files owned elsewhere.

- [ ] Add tests asserting MI(constant,constant)=0, dependence-score invariance under common cell permutation and finite CLR on zero-variance backgrounds.
- [ ] Replace row-index tie splitting with a value-based discretization; record zero-information behavior.
- [ ] Fit Ridge per target excluding the target predictor, retain sign, and test against an independently fitted sklearn Ridge.
- [ ] Match TF identifiers consistently preserving original expression names; reject ambiguous case collisions, verify mouse symbols and zero-match behavior.
- [ ] Run classical/infer/plugin and affected TF tests; provide a short change manifest.

## Task 3 — provenance, backtest and failure status

**Owner/files:** provenance worker: `report.py`, `pipeline.py`, `cli_infer.py`, `cli_run.py`, `benchmarks/beeline_backtest.py`, `hpc/`, new provenance/backtest tests. Do not modify evaluation worker files. Coordinate candidate-universe forwarding with evaluation worker.

- [ ] Add tests for complete file hashing, nested input directories, same-size content changes, provenance retention and nonzero status on actual failures.
- [ ] Stream full SHA256, hash recursive directory inventories, record effective genes/TFs/gold, all runtime versions and execution params. Avoid scanning outside the named input.
- [ ] Persist backtest edge lists and inference provenance per method; export input/feature/gold identity, args, seeds and dataset source info. Refuse stale output reuse across incompatible manifests.
- [ ] Implement documented BEELINE selection against a pinned official reference and fixture. Remove universal empirical-ceiling/pass claims. Clearly label protocol/version, validate numerical metrics on the same edge input.
- [ ] HPC setup must fail on installation/test errors. Partial failures must be visible in manifest and exit code while retaining useful diagnostic outputs.
- [ ] Generate the inference candidate universe from effective TFs and genes and pass it to evaluation; coordinate the exact `--universe` interface.
- [ ] Run targeted tests; supply verified commands and limitations.

## Task 4 — independent-unit inference and calibration pilot

**Owner/files:** coordinator: new `eval/comparison.py`, `calibration.py`, `cli_calibrate.py`, `cli_compare.py`, `cli_main.py`, new comparison/calibration tests, `configs/calibration-pilot.json`.

**Interfaces:** `compare_units(scores, alpha=0.05, margin=0.0)` receives a wide finite DataFrame (rows independent units; columns methods). Paired Student t differences with Holm adjustment are the reference procedure, not a novel estimator. A compare CLI accepts a long table unit,method,score and requires an explicit independent-unit declaration. No automatic inference from cells or edges.

- [ ] Test paired differences against scipy t-test; identical methods unresolved; constant nonzero sample differences conservatively unresolved with zero-variance warning; Holm monotonicity and simultaneous Bonferroni t intervals.
- [ ] Return signed differences, unadjusted/adjusted p-values, simultaneous intervals and directed supported superiority at fixed margin; no equivalence tiers. Validate duplicate/missing units and n>=3.
- [ ] Pilot with known Gaussian method-score means, shared unit noise, correlated methods, null and partial-null cases. Also include deliberate dependence violations to expose assumption limits. Infer biological meaning only when real data supply appropriate units.
- [ ] Count familywise false superiority, simultaneous coverage, power and unresolved fraction; return one-sided exact binomial limits. Gate uses prespecified alpha/epsilon and bounds, with pass/fail/inconclusive separated.
- [ ] Record configuration, seeds, versions and per-experiment results; run the 200-replication pilot and save artifacts. This validates only the specified score-level model, not GRN biology or MNAR identifiability.

## Task 5 — integration, replay and documentation

**Owner/files:** coordinator: README, CHANGELOG, pyproject version, `requirements/`, CI definition, protocol docs and integration tests.

- [ ] Replace unverified novelty/ceiling/significance claims with current behavior; preserve old numerical results only in the dated audit.
- [ ] Document fixed U, sparse-output policy, sensitivity vs inference, negative controls, finite-sample AP baseline and unresolved vs equivalent.
- [ ] Pin a reproducible environment and run a clean install/test check within available runtime limits.
- [ ] Execute a small end-to-end count-data → two-method inference → explicit-U evaluation run, save edges, manifests and report; independently check AP/ROC against sklearn on identical data.
- [ ] Run full test suite once integrated; run optional arboreto checks if installation is feasible; inspect failures rather than hiding skips.

## Task 6 — review and completion

- [ ] Independent scientific/code review of evaluation and pilot; address material findings and rerun affected tests.
- [ ] Verify all changed source paths against backup and no unintended modifications; preserve red/green and final logs.
- [ ] Refresh Graphify and record acts ledger. Complete compact post-mortem.
- [ ] Report implemented behavior, test results, pilot outcomes, and remaining publication gates (external protocol audit, mechanistic simulators, independent biological validation). Those research gates are deliberately not claimed complete by this software plan.
