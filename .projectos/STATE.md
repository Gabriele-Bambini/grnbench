# State - grnbench (where things left off)

Keep it short: every PR that changes focus, in-progress work or next steps updates this file.
Live runs are NOT tracked here: they are GitHub issues labelled `running` (see `python -m grnbench.board`).
Updated: 2026-10-06 (EXP-20261003-GB-1: review fixes and record note 2.1 on PR #55, agent);
2026-10-03 (issue #42, GB's agent); 2026-10-02 by LD's and LM's agents. PR #3, PR #7 and the shared-workflow PR were merged by the owner
after independent agent reviews, not human reviews (GB was away).

## Focus
Product (D-012/D-013): context-aware GRN data-cleaning engine + interface + mounted models +
comparison, published as a tool with a paper. Headline claims from real data (D-015), through
pre-registered experiments (D-017).

## In progress
- EXP-20261003-GB-1 (exploratory, issue #52, branch `gb/motif-sensitivity`): registered, stage 1
  Pearson (C6 primary, C4 secondary), tests T1-T3 (T3 = beyond a TF-generic gene propensity),
  not launched. Needs the kept artifacts copied from the login node's `/tmp` first (#41). The
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
  0.044; also under doublets and damaged cells). Claim C7 `confirmed`; D-021 ACTIVE: the engine now
  returns the baseline when no reference is given. Gene-set-preserving reference screening: +0.023,
  inconclusive, not shipped.
- EXP-20261002-LM-1 (confirmatory, fresh-set rerun of C4/C5, HPC job 679411, merged in #31): the
  motif prior lifts the AUPRC ratio on 20 never-used datasets for Pearson +0.203, RegDiffusion
  +0.234, GRNBoost2 +0.146 (n = 20 each, p_Holm <= 1.6e-6); all supported. Claim C6 `confirmed`.
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
  `reproduced-computational` (agent rerun on the HPC, identical numbers, issue #15); human pending. One dataset lost to a RegDiffusion zero-variance error (issue #13).
  Edges and processed data (190 MB) stay on the cluster in `~/gb_runs/EXP-20260926-LM-1`.
- Issue #6 merged (#11): scrublet pools batches under 100 cells.
- EXP-20260925-LM-2 (confirmatory, merged in #7): the motif prior lifts Pearson's AUPRC ratio by
  +0.182 [+0.138, +0.227], n = 20, p_Holm 1.4e-7 (H1 supported; the EPR lift is unresolved).
  An independent agent reran it from the record on 2026-09-26 and got identical numbers (record
  section 4, issue #8). Paper claim C4 is therefore `reproduced-computational`. The human
  reproduction is still pending, so C4 is not `reproduced`.
- The scrublet fix (#3) is merged. Rebuild your venv from the lock (`scikit-image` is new).

## Next steps
1. GB: reproduce EXP-20260925-LM-2 from the record (#24: GB directs and checks on his own account; an
   agent may type the commands; on another OS or machine if possible),
   then fill section 4 and set C4 to `reproduced`. Review, or challenge, the agent-reviewed merges
   #3, #7 and the workflow PR (AGENTS.md section 5: open an issue citing the PR).
2. (done, #11) Issue #6 small-batch scrublet rule.
3. GB: reproduce EXP-20260926-LM-1 under the same rule (C5 is only `reproduced-computational`). Issue #13 is
   fixed (#16). A fresh-set rerun would remove the "seen in LM-2" caveat of C5.
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
