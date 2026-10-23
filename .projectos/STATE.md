# State - grnbench (where things left off)

Keep it short: every PR that changes focus, in-progress work or next steps updates this file.
Live runs are NOT tracked here: they are GitHub issues labelled `running` (see `python -m grnbench.board`).
Updated: 2026-10-21 (R bridge, step 2b, #135, GB's agent); 2026-10-15 (#87, #71, #69 merged at GB's request, GB's agent); 2026-10-10 (interface branch merged with main after #55, GB's agent); 2026-10-09 (interface redesign; tool end to end, GB's agent); 2026-10-07 (interface explorer, GB's agent); 2026-10-06 (EXP-20261003-GB-1: review fixes and record note 2.1 on PR #55, agent); 2026-10-03 (issue #42, GB's agent); 2026-10-02 by LD's and LM's agents. PR #3, PR #7 and the shared-workflow PR were merged by the owner
after independent agent reviews, not human reviews (GB was away).

## Focus
Product (D-012/D-013): context-aware GRN data-cleaning engine + interface + mounted models +
comparison, published as a tool with a paper. Headline claims from real data (D-015), through
pre-registered experiments (D-017).

## In progress
- Preprocessing engine step 2b, optional R bridge (GB, `gb/r-bridge`, issue #135, D-032 proposed):
  `grnbench.rbridge` finds an Rscript (conda env `grnbench-r` by default), probes its packages once per
  cache, runs R scripts with a fixed seed and RNG kind and records the session; `grnbench env` and the
  interface's Environment popover show what is there and the install command for the rest. Checked with a
  fake Rscript and against the real `grnbench-r` env on GB's Windows PC; not on Linux or the cluster yet
  (installing R there needs the owner's OK). Planner behaviour unchanged. Next: step 2c (after PR #133)
  wires `requirement_status` into the knowledge base's availability, R options opt-in.
- More models (GB, `gb/more-models`): eleven literature methods registered as *extra* (`bicor`,
  `wgcna_tom`, `propr_rho`, `aracne_dpi`, `mrnet`, `tigress`, `bbsr`, `pcnet`, `extra_trees`,
  `gradient_boosting`, `community`); `--methods available` still runs exactly the 17 older ones.
  Checked on unit properties and simulated planted networks only: no real-data score exists, so none
  is a result. `docs/models_catalog.md` lists differences from the originals and what is not built in.
  Merged on 2026-10-15 (LD approval; LM's independent agent review, not a second-human review). Next: an
  exploratory record to score them on real data (HPC, GB's own account).
- Interface (GB, `gb/ui-explorer`): quick run page, model picker (presets, family, speed, tier), data explorer
  (cell maps before/after cleaning, genes and TFs) and results explorer (edges, network, model agreement,
  precision-recall, per-TF recovery) added to `grnbench ui`. Logic is in `ui/logic.py` and tested without
  Streamlit; the pages are rendered by Streamlit's AppTest and were checked in a browser on the demo dataset. The
  picker lists whatever the registry holds (the extra models of #69 appear once it lands). Reviewed against the other open PRs
  (2026-10-07): fixes for ties, detection, TF denominators, reference choice and edge validation are in; still open and
  not in this PR: the Fast preset ticks default-tier models only (extras via "Fast, including extras"), thread budget and
  the run provenance file wait for #40, `pr_figure`/`tf_recovery` duplicate `viz.py` (one definition per quantity, follow-up).
  Walked through end to end as a newcomer on 2026-10-09 (Chromium + AppTest): the one blocker and every major
  interface defect found are fixed or mitigated on this branch, with a test per behaviour change
  (`tests/test_ui_regressions.py`); what is left is the last section of #76, and the screening ranking rule
  behind one of them is a decision for LM (#75).
  Redesigned on 2026-10-09 (screen list `ui/steps.py`, stylesheet `ui/style.py`, Streamlit options `ui/theme.py`):
  three menu groups, a header and numbered cards on every screen, empty states with a way forward, results first after
  a quick run, plain-language quality verdict and pipeline overview. Looked at in headless Chromium (light, dark, phone
  width), by AppTest, and by one independent agent session reading the screenshots as a newcomer (its interface-only
  findings are fixed; open: the compare step still preselects a winner that may leave out cell types, see #75, and the
  menu ticks Profile and Pipeline only once they were on screen). Nobody but GB's agent has used it. Next: LM or LD click through it on a real dataset and say what
  is confusing or missing.
  Not yet tried on a real dataset or by someone else.
  LD's re-review of 8321f74 (2026-10-10): Explore data no longer makes a whole sparse matrix dense (B2: sparse
  per-gene moments, only sampled cells x chosen genes converted; the precision-recall universe is counted, not
  listed), the network is directed with one arrow per edge, drawn over the whole edge so no head is dropped
  between markers that do not overlap (B4), the two-dataset session test for B1 is in, and
  the before/after caption no longer implies aligned cells. Still open: `benchmark()` lists the full universe too.
- Tool end to end (GB, `gb/tool-fixes`, issue #74): the integrated tool was driven like a newcomer on
  2026-10-09 (CLI and interface, hard inputs). The small, local defects in `main` are fixed with a test each
  (input guards, option validation, messages, `python -m grnbench`, ...). A screening rule that lets a
  candidate scored on fewer cell types win is a decision for LM (#75; it touches the `auto` arm behind
  C4-C6, no claim is said to be wrong). Everything else is the checklist in #76.
- EXP-20261003-GB-1 (exploratory, issue #52): **done 2026-10-10** (HPC job 692696, commit 42c6fb2,
  40/40 units, every dataset through gate 1 and gate 2). Verdicts from `SENSITIVITY.md`: C6 (primary)
  inconclusive, C4 mixed (T3 not supported; a TF-generic gene propensity reaches the lift), combined
  inconclusive: the caveat on C4/C6 stays as worded. Record section 3; results in
  `benchmarks/results/EXP-20261003-GB-1/`. Next: the caveat wording (#41, #45); #86 and #88 can now land.
  History: registered, stage 1
  Pearson (C6 primary, C4 secondary), tests T1-T3 (T3 = beyond a TF-generic gene propensity). Needs the kept artifacts copied from the login node's `/tmp` first (#41). The
  review fixes on PR #55 (2026-10-06) changed the script and its tests only. The independent agent
  review for LM (same day) asked for dated notes in section 2 before launch: record section 2.1 now
  has them (provenance, simulations as rationale only, cross-set order, motif-sharing TFs with
  descriptive diagnostics `gene_tf_np` and motif-sharing shares, Windows rerun procedure). LM answered
  on 2026-10-10 (option c: not the second human before Monday 2026-10-12): the review of `661514b` stands as
  the independent agent review (AGENTS.md section 5), GB merges and launches, and section 2.2 of the
  record adds LM's note 7 (C6 and C4 combine only with T1-T3 run on at least 15 analysed datasets each).
  #86 and #88 must not be in the commit the job runs (LM holds them until it has run).
- Issue #56 (LD, `ld/portable-diagnostics`): schema-v1 `diagnostics.json` gives prep, screening and
  model outcomes one portable structure for reports/UI, with bounded messages, before/after
  dimensions, effective-vs-requested strategy, real subset context labels and recursive privacy
  validation. Identifier mapping/collision policy and goal-aware selection remain separate
  follow-ups to #34.
- Issue #42 (GB, `gb/raw-count-check`): ingest warns, never refuses, on integer input that looks
  non-raw (rounded normalized, 0/1) and gives specific advice when it refuses. Calibration (LM, HPC
  job 681752, review of PR #53): 56 real files (CELLxGENE + 10x raw/filtered), would warn 0; the
  tight-CV signal is safe there only with its "no 1s" condition (lowest CV 0.164 vs 0.15). A
  refusal still needs its own record and PR. Item 3 (count recovery) comes later, own D-xxx.
- EXP-20261004-LM-2 (exploratory, 8 public 10x raw matrices, #51): from raw droplets the engine's
  knee calling matches Cell Ranger's cell count and the vendor-filtered recovery on 7/8 datasets;
  heart_10k_v3 exposed a knee-floor bug (10 -> 100 UMIs, DropletUtils default), fixed; then 8/8 cell
  counts within 0.93-1.03x of Cell Ranger (recovery rerun for heart_10k only).
- EXP-20261004-LM-1 (exploratory, job 681527): scored on the genes a candidate shares with the
  baseline, no cleaning change beats the fixed baseline by more than ~0.01 (hvg=2000: +0.037 own
  universe -> +0.008 shared). HVG default unchanged; no confirmatory test of gene-set changers.
- EXP-20261003-LM-1 (confirmatory, job 681213): without a reference the fixed baseline beats the
  engine's former reference-free choice by +0.089 [+0.014, +0.163] (n = 20 unused datasets, p_Holm
  0.044; also under doublets and damaged cells). Claim C7 `confirmed`: the rerun directed by GB and typed by an
  agent session (job 702612, record section 4b, #95; signature pending) has both primaries identical,
  but the registered row check says DIFFERENT because of per-candidate timings inside lists. A corrected
  check (`benchmarks/diff_rows_nested.py`) would apply only to a new rerun registered after it. D-021 ACTIVE: the engine now
  returns the baseline when no reference is given. Gene-set-preserving reference screening: +0.023,
  inconclusive, not shipped.
- EXP-20261002-LM-1 (confirmatory, fresh-set rerun of C4/C5, HPC job 679411, merged in #31): the
  motif prior lifts the AUPRC ratio on 20 never-used datasets for Pearson +0.203, RegDiffusion
  +0.234, GRNBoost2 +0.146 (n = 20 each, p_Holm <= 1.6e-6); all supported. Claim C6 `reproduced`
  (GB's second-person rerun, job 693269, record section 4b signed 2026-10-17, #94): the three
  primaries and the report (timing aside) are identical; one planner-warning count differs (4 -> 3,
  one dataset, not a score; most likely the dataset README, which the data check does not hash).
- EXP-20261002-LM-2 (exploratory, HPC job 679448, issue #28, branch `lm/candidate-study`): every
  cleaning candidate on 20 real datasets x 7 controlled degradations.
- Issue #30: the two scoreless EXP-20261001-LD-1 rows were reproduced without changing that
  experiment. `pbmc3k_raw` lost the raw-droplet flag when screening rebuilt its profile, so all
  candidates treated empty barcodes as cells; `paul15` has a case-fold gene collision
  (`LPCAT4`/`Lpcat4`) that made TF matching refuse affected contexts. Focused regressions and fixes
  are on `ld/scoreless-real-rows`; isolated HPC validation at commit `57e737a` produced zero failed
  contexts and finite Pearson/CLR scores for both datasets. Post-merge full suite: 427 passed.
- Issue #19: a read-only HPC audit found no random corruption; the setup cache is broader and more
  dynamic than the cumulative provenance registry. A small follow-up adds optional exact-snapshot
  checking and clarifies the data-root workflow; no dataset or past result is changed. Issue #21
  was decided by the owner and recorded as D-020 for future CollecTRI/DoRothEA snapshots.
- EXP-20260926-LM-1 (confirmatory, HPC job 670902, results PR on `lm/exp-hpc-motif-models`):
  the motif prior lifts the AUPRC ratio for RegDiffusion by +0.184 [+0.132, +0.237] (n = 19,
  p_Holm 1.2e-6) and for GRNBoost2 by +0.130 [+0.108, +0.153] (n = 20, p_Holm 3.6e-10); both
  supported. EPR: GRNBoost2 lift resolved (+0.61), RegDiffusion unresolved (about 0). Claim C5
  `reproduced` (GB's second-person rerun on the HPC, job 693267, every non-timing row field, and the
  report apart from its timing columns, identical; record section 4b signed 2026-10-17, #93; earlier agent rerun identical, #15). One dataset lost to a RegDiffusion zero-variance error (issue #13).
  Edges and processed data (190 MB) stay on the cluster in `~/gb_runs/EXP-20260926-LM-1`.
- Issue #6 merged (#11): scrublet pools batches under 100 cells.
- EXP-20260925-LM-2 (confirmatory, merged in #7): the motif prior lifts Pearson's AUPRC ratio by
  +0.182 [+0.138, +0.227], n = 20, p_Holm 1.4e-7 (H1 supported; the EPR lift is unresolved).
  An independent agent reran it from the record on 2026-09-26 and got identical numbers (record
  section 4, issue #8). Paper claim C4 is therefore `reproduced-computational`. The human
  reproduction is still pending, so C4 is not `reproduced`.
- The scrublet fix (#3) is merged. Rebuild your venv from the lock (`scikit-image` is new).

## Next steps
1. (done 2026-10-08, #24, record section 4b) GB's rerun of EXP-20260925-LM-2 on the Bocconi HPC
   (Linux, job 690538): `compare_runs.py` verdict DIFFERENT, exact and `--rtol 1e-9` (H1 +0.183
   [+0.138, +0.229] vs +0.182 [+0.138, +0.227]; numeric drift across platforms, planner choices
   changed in 3 datasets). By #24's rule C4 stays `reproduced-computational`. D-026 (#89) applies only
   to reproductions launched after it; C4's route to `reproduced` is a rerun directed by GB on the
   registered platform class (Windows, venv from the lock, Python 3.11). Still open: review, or challenge, the
   agent-reviewed merges #3, #7 and the workflow PR (AGENTS.md section 5: open an issue citing the PR).
2. (done, #11) Issue #6 small-batch scrublet rule.
3. (done 2026-10-15, #93, record section 4b) GB's rerun of EXP-20260926-LM-1 on the Bocconi HPC (job
   693267): `compare_runs.py` verdict MATCH, exact and `--rtol 1e-9`; 4b signed, C5 `reproduced`.
   C6 (#94) and C7 (#95) reruns: their own PRs.
4. Review and merge issue #30's validated fix; then use EXP-20261001-LD-1's real ties plus synthetic
   breakpoints to register a separate confirmatory auto-vs-baseline run. Scheduled monitoring stays
   disabled.
5. Recalibrate `context_scope` in `selector.py` with bench evidence (add a "goal" input); regenerate the KB doc.
6. After review of EXP-20261001-LD-1: register the confirmatory synthetic + real engine run under D-020.
7. GUI: per-model best pipeline, knee cell calling, mito-gate decision, motif models, viz figures.

## Blockers / open questions
- Publication route (JOSS / Bioinformatics App Note / Nature-family journal) - team decision.
- Licence and copyright holder (D-018, proposed: keep MIT).

## Notes for the next agent
- Never quote the 0.6.0 numbers (retracted). Unresolved != equivalent.
- `reproduced-computational` is not `reproduced`. Never call C4 "reproduced by a second person".
- Heavy work (suites, runs) goes off laptops: HPC for experiments, a workstation for suites (AGENTS.md section 4).
- README top still carries the 0.7.0 framing; a README rewrite for the product direction needs team approval.
- The long "improved plan" list of 2026-09-23 is in the previous STATE.md (`git show 4de687d:.projectos/STATE.md`).
