# State - grnbench (where things left off)

Keep it short: every PR that changes focus, in-progress work or next steps updates this file.
Live runs are NOT tracked here: they are GitHub issues labelled `running` (see `python -m grnbench.board`).
Updated: 2026-09-27 by LM's agent. PR #3, PR #7 and the shared-workflow PR were merged by the owner
after independent agent reviews, not human reviews (GB was away).

## Focus
Product (D-012/D-013): context-aware GRN data-cleaning engine + interface + mounted models +
comparison, published as a tool with a paper. Headline claims from real data (D-015), through
pre-registered experiments (D-017).

## In progress
- Nothing running. No open `running` issues.
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
1. GB: reproduce EXP-20260925-LM-2 by hand from the record (on another OS or machine if possible),
   then fill section 4 and set C4 to `reproduced`. Review, or challenge, the agent-reviewed merges
   #3, #7 and the workflow PR (AGENTS.md section 5: open an issue citing the PR).
2. (done, #11) Issue #6 small-batch scrublet rule.
3. GB: reproduce EXP-20260926-LM-1 by hand (C5 is only `reproduced-computational`). Issue #13 is
   fixed (#16). A fresh-set rerun would remove the "seen in LM-2" caveat of C5.
4. First confirmatory synthetic + real engine experiment on HPC (auto vs baseline cleaning).
5. Recalibrate `context_scope` in `selector.py` with bench evidence (add a "goal" input); regenerate the KB doc.
6. Simulator: partially cell-type-specific networks, then rerun synthetic.
7. GUI: per-model best pipeline, knee cell calling, mito-gate decision, motif models, viz figures.

## Blockers / open questions
- Publication route (JOSS / Bioinformatics App Note / Nature-family journal) - team decision.
- Licence and copyright holder (D-018, proposed: keep MIT). Third member's GitHub username.
- No TeX engine on Luca's laptop: `.tex` is generated, PDFs need latexmk/tectonic (START_HERE, from zero).

## Notes for the next agent
- Never quote the 0.6.0 numbers (retracted). Unresolved != equivalent.
- `reproduced-computational` is not `reproduced`. Never call C4 "reproduced by a second person".
- Heavy work (suites, runs) goes off laptops: HPC for experiments, a workstation for suites (AGENTS.md section 4).
- README top still carries the 0.7.0 framing; a README rewrite for the product direction needs team approval.
- The long "improved plan" list of 2026-09-23 is in the previous STATE.md (`git show 4de687d:.projectos/STATE.md`).
