# Decisions - grnbench
Append-only ADR-lite. Supersede instead of editing history.

    ### D-001 | YYYY-MM-DD | agent | ACTIVE | short title
    - decision: what was chosen
    - why: reasons / evidence
    - alternatives: what was rejected and why
    - supersedes: D-000 (optional)

STATUS: PROPOSED | ACTIVE | SUPERSEDED by D-xxx | REVERTED


### D-001 | 2026-09-05 | backfill | ACTIVE | Counts or nothing: raw integer counts only
- decision: Ingest rejects anything that does not look like raw integer counts; normalized/log inputs raise IngestError loudly
- why: Preprocessing expects raw integer counts and handles depth normalization explicitly. Feeding normalized data breaks QC gates (mito/ribo %), TF discovery, and strategy reversibility.
- alternatives: Could auto-detect and auto-normalize: loses transparency and couples steps unnecessarily. Could accept but warn: silently wrong answers worse than loud failure.

### D-002 | 2026-09-05 | backfill | ACTIVE | No imputation: strategies cannot smooth or impute data
- decision: The strategy registry is explicitly prohibited from containing smoothing/imputation methods; a test enforces this.
- why: Imputation hides data loss and changes the inference problem. Counts carry more information than smoothed/imputed values. Making the choice visible in preprocessing options is better than burying it.
- alternatives: Could allow imputation as an optional strategy: too easy to use without noticing its assumptions. Could document it as a separate tool: loses the force of design constraints.

### D-003 | 2026-09-05 | backfill | ACTIVE | HVG never drops TFs: kept features = HVGs ∪ TFs detected in data
- decision: HVG filtering is followed by union with TF list; no TF is removed even if not highly variable
- why: A TF can be present in the data but show low variance (housekeeping, stochastic dropout in scRNA-seq). Removing it violates the model contract and silently breaks inference on that regulator. TF list defines the biological task scope.
- alternatives: Could require users to filter TF list: places burden on user and error is silent. Could autoselect TFs by variance: biologically wrong; causal regulators can be lowly variable.

### D-004 | 2026-09-05 | backfill | ACTIVE | QC gates fail loudly: extreme thresholds abort execution
- decision: Thresholds on cell/gene counts, mitochondrial %, ribosomal %, doublet rate and batch imbalance fail with exit code 2 unless --force is given and recorded.
- why: Silent filtering produces garbage outputs with correct-looking manifests. Loud failures force users to choose: accept risk explicitly or tune thresholds. The force flag is recorded in the manifest.
- alternatives: Could auto-tune thresholds: loses transparency and shifts blame to the tool. Could warn silently: the worst outcome (wrong result, confident user).

### D-005 | 2026-09-05 | backfill | ACTIVE | Full provenance: every step logged with parameters, versions, thresholds, input fingerprints
- decision: Manifests record parameter values, effective gene/TF lists, SHA256 input hashes, dependency versions, candidate universe and seed.
- why: Reproducibility is the baseline for any benchmark claim. Provenance must include enough to recreate the same measurement even if no code is available. Hidden choices (default thresholds, autosort) are the enemy of reproducibility.
- alternatives: Could log only major parameters: misses the edge cases that usually break reproduction. Could assume code version suffices: code can change its defaults between runs.

### D-006 | 2026-09-05 | backfill | ACTIVE | Unavailability is data, not a crash: missing deps stay listed with reason
- decision: Methods with missing external dependencies (arboreto, java/jar, Docker image, broken plugin import) are skipped loudly and stay in the available list with their reason (e.g. "arboreto not installed").
- why: Discovering which methods can run is a runtime property that varies per environment. Silently omitting methods loses information. Explicitly stating "method X unavailable: reason Y" is honest and helps debugging.
- alternatives: Could crash if any requested method is missing: forces the user to install everything. Could silently drop unavailable methods: lies about what was computed.

### D-007 | 2026-09-15 | backfill | ACTIVE | Null-calibrated scores: every metric includes (hits+1)/(n+1) p-value vs fixed-universe permutation
- decision: AUPRC, AUROC, EPR and precision@k are reported alongside p-values from score permutation on the *same* fixed candidate universe. A method can lead the table and still be indistinguishable from its own null.
- why: Leaderboards assert a total order over methods. That order is untested. Null calibration reveals whether observed differences exceed what random ranking produces. Conservative (hits+1)/(n+1) p-values avoid overstating significance.
- alternatives: Could report only observed metrics: no indication of whether differences matter. Could use resampling of the reference instead of the scores: conflates ranking uncertainty with reference incompleteness.

### D-008 | 2026-09-15 | backfill | ACTIVE | Resolved comparisons only: paired t-test with Holm correction over independent experimental units (datasets)
- decision: grnbench compare runs paired Student t on declared independent units, Holm-corrects both directions of every pair, Bonferroni simultaneous intervals, and refuses to call an unresolved pair "equivalent".
- why: BEELINE-style benchmarks have independent units (datasets/organisms), not edges (millions) or cells (billions). Confusing these inflates false positives from 0 resolved to 80%. Holm is valid under dependence if individual p-values are valid.
- alternatives: Could use Wilcoxon signed-rank: different estimand (median vs mean difference). Could use false-discovery rate: less conservative but requires more resolved pairs to publish. Could skip testing: leaderboard asserting a total order anyway.

### D-009 | 2026-09-05 | backfill | ACTIVE | Plugin-first architecture: custom models discoverable via entry point, CLI spec file, or Python file reference
- decision: User methods are registered via three pluggable routes: [project.entry-points."grnbench.methods"] in pyproject.toml, --plugin path/to/model.py:ClassName for scripts, --plugin spec.json for Docker/CLI with env variable placeholders.
- why: No fixed method list means users can add their own without forking the code. All three paths go through the same edge-schema registry, manifests and evaluation, so plugins are first-class citizens. Lowering barriers to benchmarking increases adoption.
- alternatives: Could require packaging as a PyPI package: too much friction. Could only support one plugin path: loses flexibility. Could run plugins in different evaluation branches: violates consistency.

### D-010 | 2026-09-05 | backfill | ACTIVE | Not a leaderboard: geneRNIB handles rankings; grnbench audits ranking validity
- decision: grnbench is positioned as an audit tool, not a ranking system. geneRNIB is the live benchmark with a leaderboard. grnbench asks: "Is this ranking real?" not "Who wins?"
- why: Leaderboards assert a total order. That order is untested and unstable across incomplete gold standards. Benchmarking infrastructure should be pluggable (plugins, strategies, null models); benchmarking philosophy should be honest (admit when design cannot rank).
- alternatives: Could build a leaderboard anyway: then own the false positives. Could try to scale to geneRNIB: different goals and tooling (Python vs modular CI). Could deprecate the tool: leaves a gap for auditing harnesses.

### D-011 | 2026-09-13 | backfill | SUPERSEDED by D-012 | Direction: metascience study on when GRN benchmarks support superiority claims
- decision: Pursue the audit-recommended conditional GO: restricted methodological pilot E0-E4 with gates G0-G5; software/portfolio as fallback
- why: Audit 2026-09-13: current benchmark paper NO-GO (overlap + flawed statistics); metascience direction is the recommended conditional GO
- alternatives: Publish current plugin/null benchmark (NO-GO); software-only / geneRNIB component (engineering GO, lower scientific value)

### D-012 | 2026-09-23 | claude | ACTIVE | Direction: universal GRN data-cleaning tool with usable interface, pre-mounted models and side-by-side comparison
- decision: grnbench is a product: any scRNA-seq input -> intelligent, explained cleaning -> GRN-ready data -> run the built-in models -> compare them, all from a usable interface (not CLI only). Metascience paper is NOT the goal.
- why: User direction 2026-09-23 (explicit correction of the audit-driven D-011 proposal).
- alternatives: D-011 metascience study (rejected by user); CLI-only software (insufficient: interface required)
- supersedes: D-011
- update 2026-09-23: user: 'deve diventare un tool da pubblicare' -> publication target = software tool paper + released package. Interface default chosen by agent (user dismissed the question): local web app `grnbench ui` (Streamlit), revisable.

### D-013 | 2026-09-23 | claude | ACTIVE | Cleaning engine = knowledge base of stages/decisions/options + dataset profile + rule selector + tested candidate pip...
- decision: Method set by user: (1) enumerate cleaning stages, (2) critical/decision points, (3) all solutions per problem, implemented (missing deps listed with reason), (4) dataset profile with visual data analysis + context, (5) literature guidelines per option distilled into a machine-readable knowledge base (best-practice combos from top publications), (6) rule-based selector assembles best pipeline + alternatives, runs and scores them, shows pros/cons, user can edit. Must be practical and integrable.
- why: User direction 2026-09-23

### D-014 | 2026-09-23 | claude | ACTIVE | Interface = local Streamlit app (grnbench ui), thin layer over the library
- decision: 7-step Streamlit app; every screen reproducible via plan.json + exact grnbench prep commands; logic in grnbench.ui.logic (tested), app smoke-tested with AppTest
- why: User dismissed the interface question; agent default. Publishable tools commonly ship Streamlit; pip-installable; reuses pipeline directly.
- alternatives: FastAPI+React web app; desktop app; terminal wizard

### D-015 | 2026-09-24 | claude | ACTIVE | Validate the engine on real data before trusting synthetic gains
- decision: Headline claims about the automatic cleaning must come from real datasets with reference networks (CELLxGENE + CollecTRI halves) and be paired with the synthetic bench; synthetic-only gains are reported as upper bounds of the simulator's assumptions.
- why: Interim v2: synthetic +42% vs real ~ties; rules alone hurt on real data

### D-016 | 2026-09-25 | claude (for LM) | ACTIVE | Canonical repo = private GitHub elle13eth/grnbench; Drive folder frozen
- decision: all code changes go through the GitHub repo (feature branches <initials>/<topic>, PRs, squash-merge, one review by another person). The Drive folder is read-only history; HANDOFF.md is replaced by PR descriptions + STATE.md.
- why: three people and their agents edited one Drive folder; Drive disconnects killed runs; no history, no review.
- alternatives: keep Drive + local git copies (diverges; no review); public repo now (not ready to publish).
- note: pushed as a squashed, sanitized snapshot: the old history contained personal absolute paths and host names (no keys or tokens found).

### D-017 | 2026-09-25 | claude (for LM) | ACTIVE | One experiment registry with pre-registration; paper claims only from reproduced confirmatory runs
- decision: every run gets experiments/EXP-YYYYMMDD-<initials>-<n>.md committed before it starts. Confirmatory runs fix hypothesis, endpoint, unit = dataset, datasets, exclusions and the comparison family; tested with compare_units (paired, Holm). paper/CLAIMS.md maps claims to EXP IDs; a claim is usable only after a second person reproduced it. Tables/figures come from paper/make_figures.py. ATTEMPTS.md keeps non-run attempts only.
- why: the goal is a publishable paper (team request 2026-09-25); D-008 (datasets as units) and the 0.6.0 retraction show what happens without declared families.
- alternatives: ATTEMPTS.md for runs (one file, merge conflicts between three people, no pre-registration); a lab notebook outside git (not traceable to commits).

### D-018 | 2026-09-25 | claude (for LM) | ACTIVE | Licence: keep MIT ("The grnbench authors")
- decision: keep the existing MIT LICENSE; confirm the copyright holder wording and the author list (CITATION.cff) before any public release.
- why: MIT is OSI-approved and accepted by journal code policies; nothing in the repo conflicts with it. Redistributed data keep their own licences (CELLxGENE: CC BY 4.0).
- alternatives: Apache-2.0 (patent grant), BSD-3. Needs team confirmation.
- status: ACTIVE since 2026-10-16 (LM, on #73). MIT is confirmed; the copyright holder is "The grnbench authors"
  (`LICENSE` and `pyproject.toml` changed from "GRN-Bench contributors"). The author list is not decided here: it goes into `CITATION.cff` later.

### D-019 | 2026-10-01 | codex (for LD) | ACTIVE | Context-specific synthetic truth and controlled breaking-point validation
- decision: Synthetic engine validation uses partially shared cell-type GRNs (shared-edge regimes
  100%, 75%, 50%, 25%), standardises planted TF activity inside each cell type, splits the union
  truth before expanding it by context, and macro-scores pooled or context-specific predictions
  against each context's held-out truth. Stress curves change one factor at a time at four declared
  levels (depth, cells, doublets, ambient RNA, damaged cells, batch confounding, batch depth and
  type imbalance), plus four mixed levels and five seeds. Exploratory breaking points are: first
  level whose 95% t interval no longer lies wholly above random, first level at or below half the
  clean mean, and first level with any missing model score. Expensive models share the `auto`
  pipeline selected by the fast proxy family on the selection half; exploratory runs use 100
  GRNBoost2 estimators and 100 RegDiffusion steps and report that choice. Runners pin tasks and
  score-affecting configuration by hash, retry transient rows, retain terminal errors, write rows
  atomically and refuse an incomplete final report unless it is explicitly watermarked partial.
- why: The old simulator assigned one network to all types while giving types different mean
  expression programs, so splitting by type could look beneficial without testing biological
  network heterogeneity. Progressive stress curves answer the owner's requested worst-case
  question, while explicit completeness prevents the partial 17/40 real-data run from being
  mistaken for a finished experiment. Approved by the owner on 2026-10-01 (issues #22 and #23).
- alternatives: keep one shared truth (cannot test context specificity); score against the union
  network (credits a context for edges true only elsewhere); run every heavy model on every
  rejected candidate (prohibitive and lets the evaluation design dominate compute); accept partial
  reports silently (repeats the earlier failure mode).
- compatibility: D-019 intentionally changes the default simulator composition, regulatory-activity
  standardisation and random streams. Synthetic runs from before D-019 (including
  `engine-2026-09-24` and `engine-2026-09-25`) therefore do not regenerate from current `main`;
  their committed outputs remain historical records. A checksum regression test pins the new
  default so future drift is explicit.

### D-020 | 2026-10-02 | codex (for LD) | ACTIVE | Immutable reference snapshots per registered experiment
- decision: Reproductions of EXP-20260925-LM-2 and EXP-20260926-LM-1 use their existing registered
  reference hashes. New registered experiments prefer an upstream versioned release when it
  provides the exact required resource and organism; otherwise they create a new experiment-specific
  snapshot manifest with retrieval date and SHA256 values. The data root must pass
  `data_manifest.py --check <snapshot> --require-all` before submission, the runner records the
  actual reference hashes in `RUN_CONFIG.json`, and an existing registered hash is never replaced.
  Runs using different snapshots are not directly comparable unless rerun on the same snapshot.
- why: Pinning `decoupler` does not pin the live OmniPath response. The read-only audit in #19 found
  different bytes for all four regenerated CollecTRI/DoRothEA files. The owner approved this policy
  in #21 on 2026-10-02.
- alternatives: Always reuse the 2026-09-25 files (best historical comparability, but no durable
  public fetch is established); always take the current live response (silently changes the
  biological reference); replace rows in the cumulative registry (destroys provenance).

### D-021 | 2026-10-03 | claude (for LM) | ACTIVE | Without a reference, the engine's choice is the fixed baseline pipeline
- decision: `screen_candidates` always screens a candidate named `baseline` (`screen.BASELINE_CHOICES`,
  the pipeline the benchmarks already use). `rank_candidates` without reference columns ranks
  `baseline` first when it has no failed context; if it is refused or fails, the previous rule
  (edge stability - depth leakage - batch leakage) picks among the remaining candidates. With a
  reference, ranking is unchanged. ACTIVE since the confirmatory EXP-20261003-LM-1 supported H1
  (baseline +0.089 [+0.014, +0.163] over the former rule, n = 20 unused datasets, p_Holm 0.044;
  claim C7). Its H2 (gene-set-preserving screening with a reference) was inconclusive, so the
  reference-based rule is unchanged.
- why: EXP-20261002-LM-2 (section 3, `benchmarks/results/EXP-20261002-LM-2/SUMMARY.md`): on real data
  the unreferenced screened choice lost to the fixed baseline in every condition. The baseline is
  sometimes refused by QC gates at very low depth (7/20 datasets at 10% depth) where adaptive
  candidates still run, hence the fallback.
- alternatives: keep the screening-only rule (loses to the baseline in the evidence); drop the
  fallback (no output when the baseline is refused); change the screening metric (a larger redesign,
  not supported by the evidence yet).

### D-022 | 2026-10-09 | claude (for LM) | ACTIVE | Next experiments test networks against knockdowns; a single-screen experiment may take the TF as its unit
- decision: (1) direction: the next registered experiments ask whether inferred edges agree with
  causal evidence rather than only with CollecTRI: CRISPRi knockdowns as the truth (Replogle 2022, RPE1
  as design screen and K562 genome-wide as replication). The cleaning engine stays the product (D-012)
  and D-010 is unchanged. TP53 in breast cancer is parked (GB, #79, 2026-10-22): it is an exploratory
  example of use, run after the engine plan, and not an aim or a headline claim. Its own record will
  declare, before any scoring, which TP53 alteration classes count as loss of canonical p53 function
  (and how missense, VUS, multi-hit and LOH states are handled) and take the unique patient as the unit,
  with one eligible primary tumour per patient or a declared deduplication rule (LD, #79, 2026-10-22).
  (2) unit of analysis: when an experiment rests on one screen, the dataset cannot be the unit
  (`experiments/TEMPLATE.md`, AGENTS.md section 7). Such an experiment may take the TF as its unit, if
  its record (a) uses a null that keeps the dependence between units (shared controls; TFs whose motif
  rows or targets are correlated, by permuting labels jointly across units rather than per unit), and
  (b) names an independent replication set analysed as its own family. A headline claim from such an
  experiment needs both sets. Experiments with several datasets keep the dataset as the unit.
- why: our claims (C4-C7) measure agreement with CollecTRI, which is built from the literature, and the
  team's earlier projects (TS-GNN, grnprep) never had a truth from an intervention. Feasibility, descriptive
  only: RPE1 has 83 of 154 perturbed TFs at knockdown >= 0.70 with more than 10 responsive genes (#67,
  job 687167). On the
  dependence point: on the kept C6 artifacts the TFs' cisTarget rows are correlated (mean pairwise
  Spearman r 0.30-0.38 in four datasets), so a per-TF permutation understates the spread of a
  shared-database null (review of #55, 2026-10-09); a TF-level test has to keep that structure.
- alternatives: keep the dataset as the only unit (one screen gives n = 1; pooling cells or edges as
  units is worse); pool TFs across screens into one test (mixes two cell lines and hides a failed
  replication); stay with CollecTRI-only validation (cannot tell recovering the literature from being
  right about regulation).
- accepted 2026-10-27: GB agreed on #79 (2026-10-22: the knockdown part goes ahead, TP53 parked as an
  exploratory example of use run after the engine plan, the wording above) and LD on #79 (review of
  2026-10-25 at db5a270: agrees with the contract and asks for ACTIVE). Records may now cite it.

### D-024 | 2026-10-10 | claude (for LM) | PROPOSED | Cells are subsampled by a hash of their barcode, not by position
- decision: screening (the 3,000-cell sample and the two halves of edge stability) and benchmark
  scoring (`MAX_SCORE_CELLS` in `real_bench.py` and `engine_bench.py`) pick cells by a 64-bit blake2b
  hash of `f"{seed}|{obs_name}"` (`screen.stable_cell_order`, `screen.stable_subsample`): first K
  of that order, returned sorted; halves = first and second half of the order. n_obs <= K keeps all
  cells. `pipeline.subsample_cells` (user-facing, stratified) and `profile.py` are not changed.
- why: EXP-20260925-LM-2, dataset lung__7b3368a5, same commit on Windows and Linux (#82): QC metrics are
  identical, but scrublet run per batch removed 82 predicted doublets on Windows and 87 on Linux
  (4,047 vs 4,042 cells left), every position after them shifted, and the positional 3,000-cell sample
  was mostly different cells. The pooled candidate screened 1.4932 on Windows and 1.5164 on Linux against
  the baseline's 1.5045 on both, and the screened choice flipped. With the hash, removing 3 of 5,000
  cells changes at most 3 sampled barcodes (positional rule: 2,046 of 3,000 shared). Checked on the same
  dataset with this change plus the symmetric-weights fix (#86): the pooled candidate's platform gap
  drops from 0.023 to 0.009 (1.5310 vs 1.5403), the baseline agrees to 4e-8 (1.559887390 vs
  1.559887452) and both platforms choose the baseline. One dataset; the per-batch scrublet difference
  itself remains.
- alternatives: keep positional sampling (one changed cell reshuffles the sample); stratified
  sampling only (does not fix the reshuffle).
- consequence: registered results keep their numbers and their kept `processed.h5ad` files. Two
  effects change what new runs sample. (a) The cap (`MAX_SCREEN_CELLS`, `MAX_SCORE_CELLS`, K = 3,000)
  changes only above K: at or below K all cells are kept. (b) The two fixed-size halves of
  `edge_stability` (150 cells each) are drawn from the same hash order, so they change below K too,
  on every matrix large enough to have two halves. Example, 800 cells named `AAAC000000-1`..`AAAC000799-1`, 150-cell halves, seed 0: the
  old and new first halves share 26 cells, the second halves 22 (chance level is about 28). So stability
  values, and anything that depends on them (the screened choice), differ from old rows on
  datasets of any size, not only above K. New runs do not reproduce old rows exactly. Do not merge into a branch that
  must reproduce registered rows (EXP-20261003-GB-1 gate 1).
- pin: the rule is named `screen.CELL_SAMPLING = "blake2b64-seed-barcode-v1"`. Both runners record it
  in `RUN_CONFIG.json` (`real_bench` also records `max_score_cells`), and `ensure_run_config` refuses to
  resume an output written without it or with another value. `report` does not check it, so committed
  runs can still be reported. A golden vector in `tests/test_stable_subsample.py` fixes the algorithm
  across processes (`PYTHONHASHSEED`).

### D-023 | 2026-10-09 | claude (for LM) | PROPOSED | With a reference, screening ranks a candidate scored on a strict subset of a viable candidate's contexts after it
- decision: in `rank_candidates` with a reference, a candidate is failed/unrankable when it failed on
  any context or any of its scores is non-finite (NaN or inf); such candidates are excluded from the
  comparison below and rank last (best stays None when no candidate is rankable). Among the other
  candidates, one whose set of scored contexts (rows without error, by the `context` column) is a
  strict subset of a non-failed candidate's set is ranked after every unflagged candidate and marked
  `fewer_contexts`; a failed superset never flags anyone (a failed pipeline is not a viable
  alternative), so a lone subset whose only superset failed is chosen and not flagged. The returned
  rule says so. The pooled `baseline` (context `all`) is not a subset of per-cell-type sets, so it is unaffected. The
  no-reference rule (D-021) and the handling of failed candidates are unchanged. PROPOSED until GB reviews.
  Amends D-021's "with a reference, ranking is unchanged".
- why: issue #75. The score is the mean over a candidate's own contexts, so `context_scope=selected_type`,
  which cleans and scores only the largest cell type, can beat `recommended` on the easiest context
  with bit-identical values there, and the cleaning that follows keeps one cell type.
- impact on the registered experiments (`python benchmarks/screen_subset_impact.py`, committed rows,
  `real_bench.models_section`; the first line of each block counts the datasets where the auto arm was
  that candidate; the rows are the declared motif pairs with those datasets dropped, then given the
  baseline arm's values):

  ```text
  EXP-20260925-LM-2: 20 rows, auto = context_scope=selected_type in 3: auto_n_contexts [1, 1, 1] vs rules_n_contexts [2, 2, 4]
  -- affected rows dropped
  | pearson_motif | pearson_signed | 17 | +0.197 | [+0.148, +0.247] | 5.16e-07 | 1 | pearson_motif |
  | pearson_motif | motif_prior | 17 | +0.171 | [+0.080, +0.262] | 0.00162 | 1 | pearson_motif |
  -- affected rows -> baseline arm values
  | pearson_motif | pearson_signed | 20 | +0.186 | [+0.142, +0.229] | 6.45e-08 | 1 | pearson_motif |
  | pearson_motif | motif_prior | 20 | +0.154 | [+0.075, +0.233] | 0.00097 | 1 | pearson_motif |
  EXP-20260926-LM-1: 20 rows, auto = context_scope=selected_type in 3: auto_n_contexts [1, 1, 1] vs rules_n_contexts [2, 2, 4]
  -- affected rows dropped
  | grnboost2_motif | grnboost2 | 17 | +0.140 | [+0.119, +0.162] | 5.34e-10 | 1 | grnboost2_motif |
  | regdiffusion_motif | regdiffusion | 16 | +0.204 | [+0.147, +0.260] | 2.02e-06 | 1 | regdiffusion_motif |
  -- affected rows -> baseline arm values
  | grnboost2_motif | grnboost2 | 20 | +0.132 | [+0.111, +0.153] | 1.03e-10 | 1 | grnboost2_motif |
  | regdiffusion_motif | regdiffusion | 19 | +0.189 | [+0.139, +0.240] | 4.49e-07 | 1 | regdiffusion_motif |
  EXP-20261002-LM-1: 20 rows, auto = context_scope=selected_type in 4: auto_n_contexts [1, 1, 1, 1] vs rules_n_contexts [4, 4, 2, 4]
  -- affected rows dropped
  | pearson_motif | pearson_signed | 16 | +0.203 | [+0.140, +0.266] | 2.48e-05 | 1 | pearson_motif |
  | grnboost2_motif | grnboost2 | 16 | +0.149 | [+0.118, +0.180] | 1.73e-07 | 1 | grnboost2_motif |
  | regdiffusion_motif | regdiffusion | 16 | +0.231 | [+0.159, +0.304] | 2.48e-05 | 1 | regdiffusion_motif |
  -- affected rows -> baseline arm values
  | pearson_motif | pearson_signed | 20 | +0.193 | [+0.143, +0.244] | 7.59e-07 | 1 | pearson_motif |
  | grnboost2_motif | grnboost2 | 20 | +0.141 | [+0.115, +0.167] | 3.07e-09 | 1 | grnboost2_motif |
  | regdiffusion_motif | regdiffusion | 20 | +0.219 | [+0.160, +0.279] | 1.05e-06 | 1 | regdiffusion_motif |
  ```
  The registered rows were produced under the old rule and stay as registered; each record carries a dated note.
  Consequence: future auto-arm choices can differ from the registered ones; reproductions of registered rows use each record's CODE_COMMIT.
- alternatives: stop forcing `selected_type` as the `context_scope` runner-up (a user-added one would stay
  unprotected); leave the rule and only disclose it (#75 options 2 and 3).

### D-026 | 2026-10-16 | LM (accepted on #73) | ACTIVE | Reproduction verdicts have two levels: exact on the registered platform class, or equivalent on another platform
- decision: Level 1, exact, on the registered platform class (OS family, Python minor version, the lock;
  CPU vendor and BLAS written in the record): all numeric row fields within rtol 1e-9. Level 2, another
  platform: the same task set and the same n per declared comparison after the registered exclusion rules;
  the same direction and the same Holm decision for every declared comparison; and the mean paired difference
  between the rerun's and the committed per-dataset primary values within a margin delta declared in section 1
  of the record (default 0.02 AUPRC-ratio units), by two one-sided t tests at alpha 0.05. Planner choices that
  flipped and the number of per-dataset primary inputs outside rtol 1e-6 are listed in section 4 and do not
  decide. Either level, directed and signed by a second person, gives `reproduced`; a rerun from the owner's
  side gives `reproduced-computational` as now. Records registered from 2026-10-16 on declare the platform class
  and delta in section 1. Section 1 of a record registered earlier is never edited: when such a record has no
  delta, a reproduction of it launched after 2026-10-16 uses the default 0.02. Level 1 is available only when every
  field of the platform-class definition is present in the record or its run artifacts; when any is missing, the
  reproduction is judged under Level 2, and missing fields are never inferred or added afterwards (LD, #103,
  2026-10-22). Scope: applies only to reproductions launched after 2026-10-16. A reproduction whose result is already
  known is never judged again under it: #82 (EXP-20260925-LM-2 rerun on another platform) stays DIFFERENT under
  the rule of #24, and C4's route to `reproduced` is a new second-person rerun judged under Level 2, since
  EXP-20260925-LM-2 does not name its CPU vendor or BLAS. The `--equivalence` option of `compare_runs.py` is a separate PR after #43 (not in the PR that
  records this decision). Approved by GB on #103 (2026-10-21) with a working condition: whenever a reproduction is
  run by an agent, GB is told before it is recorded. A reproduction is human
  (`reproduced`, "reproduced by a second person" in the manuscript) only when a second person directed and
  checked each step and signed section 4 (AGENTS.md section 8); an agent-only rerun stays
  `reproduced-computational` and the manuscript labels it so.
- why: rtol 1e-9 held within one platform (the 26/9 rerun of C4: 0 of 1,180 fields differ) and cannot hold
  across platforms (#82: 40 of 40 primary inputs differ, relative 2.2e-07 to 2.8e-02, while the headline moved
  by 0.001). A reproduction should test the finding and not only the determinism of one machine. Text proposed
  in #89, accepted by LM on #73 on 2026-10-16.
- alternatives: exactness only on a pinned container (tests determinism); applying a new rule to known results
  (a rule chosen after the result).

### D-029 | 2026-10-16 | claude (for LM) | PROPOSED | A preprint may carry `confirmed` and `reproduced-computational` claims, each labelled with its status
- decision: A preprint may carry `confirmed` and `reproduced-computational` claims; every claim shows its
  status where it appears (text, tables, figure captions); the abstract says which findings a second person
  reproduced; exploratory results are labelled exploratory everywhere; a journal submission still needs
  `reproduced`. Agreed by LD (#73, 2026-10-18) and GB (#103, 2026-10-21). GB's condition: the step to a
  journal submission is enforced by the tool, not only stated. Today `grnbench report --paper` accepts
  `confirmed`, `reproduced-computational` and `reproduced` alike (`texreport.py`), so the PR that activates this
  decision adds a submission build (`grnbench report --paper --submission`) that fails and names every claim
  that is not `reproduced`, with a test, and amends rule 3 of `paper/README.md` and section 8 of `AGENTS.md`.
  GB's side writes that PR after this one merges; until it merges, D-029 stays PROPOSED and both texts stand
  as written.
- why: issue #73. Rule 3 came with D-017 after the 0.6.0 retraction, and its wording ("submitted draft")
  predates any preprint plan. bioRxiv cannot remove a posted paper (withdrawal leaves the original version
  online), so a preprint needs the same care in what it claims; the status labels are how a reader sees which
  claims a second person has checked.
- alternatives: a preprint with only `reproduced` claims; any status, unlabelled.
### D-027 | 2026-10-10 | codex (for LD) | SUPERSEDED by D-028 | Trace identifier transforms and refuse measured unsafe collision sets
- decision: v1 maps gene identifiers only from a symbol column carried by the input or an explicit
  two-column `id -> symbol` table whose SHA256 is recorded; it never performs a live MyGene lookup.
  Supported human and mouse Ensembl version suffixes are stripped. Exact duplicates are summed and
  recorded but do not count as collision groups. Every distinct-spelling case-fold group keeps the
  first loaded spelling as representative and is aggregated by sum. Ingest refuses when such groups
  exceed 1% of input features, or when groups spanning spelling conventions (upper, mixed, lower,
  classified after removing the measured `orf` token) exceed 0.1%. An input with more than 90%
  Ensembl identifiers is refused if an input column or supplied table cannot reduce that share to
  90% or below. The complete source-to-representative audit is written to the manifest and portable
  diagnostics without paths, barcodes or expression values.
- why: #34 measured zero cross-convention groups in the separate approved human, mouse and reference
  lists, but 18.7% in their full join. A rule on all duplicate positions would falsely reject the
  bronchus dataset's 733 repeated `Y_RNA` entries; a rule on any single case group would falsely
  reject paul15's `LPCAT4/Lpcat4`. The boundary accepts two cross-convention groups in 2,000 features
  and refuses three. The mapping boundary makes runs reproducible and prevents silent zero overlap
  with symbol-based TF lists and references.
- alternatives: live MyGene mapping (network-dependent, does not resolve the measured unmatched
  reference genes); always merge every case-fold collision (can silently combine mixed-species
  annotations); refuse every duplicate (rejects valid repeated features); drop colliding genes
  (changes counts and task scope without a traceable representative).

### D-028 | 2026-10-11 | codex (for LD) | ACTIVE | Case-only identifier merging is human/mouse-specific
- decision: D-027's mapping, Ensembl, exact-duplicate and provenance rules remain. Case-insensitive
  merging and its 1% / 0.1% refusal guards apply only when the input is confidently identified as
  human or mouse, using curated `obs['organism']` metadata first, human/mouse or other recognised
  Ensembl gene identifiers second and the existing symbol-casing hint last. A non-human/mouse
  Ensembl prefix establishes the `other` scope. For a named third organism or an organism that
  cannot be established, exact duplicate names are still summed, but distinct case-only spellings
  remain separate features; every group and convention is recorded and a warning explains why it
  was preserved. Human+mouse metadata remains inside the guarded scope.
- why: exact-snapshot calibration EXP-20261010-LD-2 accounted for 33/33 approved real inputs but
  D-027 falsely refused the `Danio rerio` input (91 cross-convention groups among 32,031 genes).
  Upper/mixed/lower casing is species nomenclature, not a universal equivalence relation. Merging
  zebrafish spellings under a human/mouse rule can combine distinct identifiers; preserving them is
  the reversible safe action. Unknown inputs are treated conservatively for the same reason.
- alternatives: merge but merely disable the refusal outside human/mouse (irreversibly combines
  potentially distinct genes); invent thresholds for every organism from one zebrafish dataset
  (uncalibrated); drop all case variants (loses counts); refuse every non-human/mouse input (blocks
  approved data without evidence of corruption).
- compatibility: Inputs that remain more than 90% Ensembl identifiers without a usable symbol source
  stop before QC, because symbol-based references cannot score them. This supersedes #149's
  pure-Ensembl mitochondrial-QC warning; that warning remains applicable and tested for mixed inputs
  with 5-90% Ensembl identifiers.
- status: ACTIVE since 2026-10-28 (GB, on #98); the hard refusal remains after review against the
  loud-QC behaviour, while the proposed double-lock override remains a separate design preference.
- supersedes: D-027.

### D-031 | 2026-10-21 | claude (for GB) | ACTIVE | An option is available only if its dependencies are installed and grnbench can run it
- decision: availability is decided per (decision, option) (`DecisionPoint.availability`). An option is
  available only when its dependencies are installed AND `executors.PLAN_STEPS` maps it onto the
  pipeline (or `RECORDED_ONLY` says it runs nothing). The selector never chooses an unavailable option;
  a chosen or forced option grnbench cannot run raises a `ValueError` that names the decision, the
  option, who set it and why, instead of being replaced in silence. The executor registry refuses
  `sparsity=imputation` and `batch=scvi_denoised` (D-002), and executors that return a matrix pass
  `check_integer_counts` (D-001). The reason text `missing: <req>` is unchanged byte for byte.
- why: the 2026-09-27 audit found that an option whose dependency is installed but that grnbench
  cannot run (for example scDblFinder) could be chosen by the selector and then dropped by the plan
  without a word, so the recorded choice was not what ran (#128). In the registered environments none
  of those dependencies is installed, so plans and rows there are unchanged (golden frozen from main
  6917c88: option table, 40 whole plans, screening candidates and 451 committed rows, 0 differences).
- alternatives: keep dependency-only availability (the recorded plan can differ from what ran); fall
  back to the default with a warning (still a silent change of the user's choice, and `n_warnings` is
  a field of committed rows); hide external options from the knowledge base (loses the evidence and the
  install guidance).
### D-032 | 2026-10-21 | claude (for GB) | ACTIVE | Optional R bridge through an Rscript subprocess
- decision: the R-only cleaning options of the knowledge base (EmptyDrops, SoupX, DecontX, scDblFinder, ...)
  run through `grnbench.rbridge`: `Rscript --vanilla` in a child process, no rpy2. R comes from a dedicated
  conda env `grnbench-r` built from `requirements/` (`r-bridge.yml` on Linux, `conda create` + `r-bridge.R`
  on Windows; `docs/R_BRIDGE.md`); grnbench never installs anything. The wrapper `r/_wrapper.R` fixes
  `RNGkind("Mersenne-Twister", "Inversion", "Rejection")` and the seed and records them in `session.json`
  with the R and package versions. A missing requirement keeps the knowledge base's text `missing: <req>`.
  grnbench never switches quietly to another R: an Rscript that was named (`--rscript`, `GRNBENCH_RSCRIPT`,
  `GRNBENCH_R_ENV`) or a conda env that exists but does not work is an error, not a fall-through to PATH.
  R options become choosable only after step 2c (issue #135), and opt-in under the current preset.
- why: rpy2 is in no registered environment, ties R to the Python process (an R crash takes the run down,
  no timeout) and breaks with conda R on Windows (its DLLs are found only through PATH); an isolated env
  keeps the Python lock and every registered row unchanged.
- alternatives: rpy2 (above); the system R (versions differ per machine and the user library leaks in);
  Docker only (one more runtime on every laptop, not assumed on the cluster).
- note: D-031 is proposed in PR #133; this PR takes D-032 (renumber if main takes it first).
- update 2026-10-23 (step 2c, issue #135, PR #144): an R option also needs a grnbench executor (D-031)
  before it can be chosen. The knowledge base asks R about `r:` packages only when `GRNBENCH_R_BRIDGE=1`;
  off (the default), they count as missing without starting R, so planning never starts R by surprise and
  registered reruns keep their texts.
### D-057 | 2026-10-22 | claude (for GB) | ACTIVE | A pool of 200 never-used CELLxGENE datasets, one per study, split by a seed before any look
- decision: new real-data experiments may draw from a fixed pool of 200 CELLxGENE Discover datasets
  from studies (collections) that no earlier record drew from: human or mouse, single-cell or single-nucleus RNA of any assay except
  targeted panels and ATAC, H5AD of at most 2 GB with 500 to 500,000 cells, and one dataset per
  collection (the smallest H5AD of each study; the first 200 of the selector's round-robin over assay x
  organism). The pool comes from `fetch_cellxgene.py --one-per-collection --exclude-collections` on the
  public listing of 2026-10-22 and is committed as
  `experiments/EXP-20261022-GB.pool.tsv`. Before any dataset is downloaded or scored,
  `split_pool.py --seed 42` puts 120 datasets in a confirmatory part and 80 in an exploratory part
  (`EXP-20261022-GB-1.datasets.tsv`, `EXP-20261022-GB-2.datasets.tsv`). Confirmatory records use only
  the confirmatory part; tuning, stress tests and internal checks use only the exploratory part; a
  dataset never changes part. A confirmatory part already scored for one question may serve a
  confirmatory test of another question only if nothing was developed or tuned on it (as C7 used C6's
  datasets).
- why: the owner asked for it ("usane anche 200 di dataset, anche grandi, basta che abbia senso",
  2026-10-18; "prendi 200 datasets di tutti i tipi e fai le prove", 2026-10-22). C4's EPR question
  needs about 100 to 120 never-used datasets (EXP-20261020-GB-1); D-015 and D-017 want confirmation on
  unused data; one dataset per study, and no study used before, avoid the shared-study pairs found
  in C4 (erratum of 2026-10-09). The smallest file per study keeps compute bounded and is not chosen by content.
- alternatives: keep drawing 20 datasets at a time (too few for EPR); a random dataset within each
  study (larger files, a seed per study); no split (exploratory looks would contaminate the
  confirmatory data).
- applies: D-015, D-017 and D-020 (references of each record are its own snapshot manifest); the new
  data have an issue first (AGENTS.md section 9).
- approved: GB, 2026-10-23 ("go" on the plan of the registration PR: download, pilot, then the two jobs).
- update 2026-10-23 (LM's review, before anything was scored): "used before" means any CELLxGENE id in
  `experiments/` or `benchmarks/results/` on any branch (`benchmarks/used_datasets.py`; 71 datasets in 43
  studies); the pool is redrawn from the saved listing of 2026-10-22 (`fetch_cellxgene.py --listing`),
  dropping 6 datasets EXP-20261001-LD-1 had scored. Each record's data provenance is its own snapshot
  (`EXP-...data.tsv`, checked with `--require-all`); data and large outputs live in a data folder outside the
  code checkout (the cluster's scratch area is not writable for this account; home has no quota).
