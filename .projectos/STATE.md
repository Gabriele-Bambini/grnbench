# State - grnbench (where things left off)

Keep it short: every PR that changes focus, in-progress work or next steps updates this file.
Live runs are NOT tracked here: they are GitHub issues labelled `running` (see `python -m grnbench.board`).
Updated: 2026-09-25 by LM (setup PR: board, viz, LaTeX report, HPC skill, onboarding guide)

## Focus
Product (D-012/D-013): context-aware GRN data-cleaning engine + interface + mounted models +
comparison, published as a tool with a paper. Headline claims from real data (D-015), through
pre-registered experiments (D-017).

## In progress
- PR #3 scrublet import fix (issue #2): scikit-image in the lock, doublet fallback recorded per arm. Review: Gabriele. Merged? rebuild your venv from the lock.
- EXP-20260925-LM-2 (confirmatory) done: motif prior lifts Pearson AUPRC ratio +0.182 [+0.138, +0.227],
  n = 20, p_Holm 1.4e-7 (H1 supported; EPR lift unresolved). PR open on `lm/confirmatory-motif-prior`;
  Gabriele to review and reproduce (paper/CLAIMS.md C4 stays `confirmed` until then).

## Next steps
1. Review and merge the two PRs above; rebuild your venv from the lock; reproduce EXP-20260925-LM-2.
2. Small-batch rule for per-batch scrublet (issue #6).
3. RegDiffusion / GRNBoost2 (+ motif) on HPC as a pre-registered experiment (CPU time too long locally).
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
- README top still carries the 0.7.0 framing; a README rewrite for the product direction needs team approval.
- The long "improved plan" list of 2026-09-23 is in the previous STATE.md (`git show 4de687d:.projectos/STATE.md`).
