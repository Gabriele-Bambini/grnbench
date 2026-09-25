# State - grnbench (where things left off)

Keep it short: every PR that changes focus, in-progress work or next steps updates this file.
Live runs are NOT tracked here: they are GitHub issues labelled `running` (see `python -m grnbench.board`).
Updated: 2026-09-25 by LM (setup PR: board, viz, LaTeX report, HPC skill, onboarding guide)

## Focus
Product (D-012/D-013): context-aware GRN data-cleaning engine + interface + mounted models +
comparison, published as a tool with a paper. Headline claims from real data (D-015), through
pre-registered experiments (D-017).

## In progress
- PR "scrublet import fix" (#2 / PR #3): scikit-image in the lock, doublet fallback recorded. Review: Gabriele.
- EXP-20260925-LM-2 confirmatory motif-prior experiment on 20 fresh datasets (branch
  `lm/confirmatory-motif-prior`, pre-registered). Review + reproduction: Gabriele.

## Next steps
1. Review and merge the two PRs above; rebuild your venv from the lock.
2. RegDiffusion / GRNBoost2 (+ motif) on HPC as a pre-registered experiment (CPU time too long locally).
3. First confirmatory synthetic + real engine experiment on HPC (auto vs baseline cleaning).
4. Recalibrate `context_scope` in `selector.py` with bench evidence (add a "goal" input); regenerate the KB doc.
5. Simulator: partially cell-type-specific networks, then rerun synthetic.
6. GUI: per-model best pipeline, knee cell calling, mito-gate decision, motif models, viz figures.

## Blockers / open questions
- Publication route (JOSS / Bioinformatics App Note / Nature-family journal) - team decision.
- Licence and copyright holder (D-018, proposed: keep MIT). Third member's GitHub username.
- No TeX engine on Luca's laptop: `.tex` is generated, PDFs need latexmk/tectonic (START_HERE, from zero).

## Notes for the next agent
- Never quote the 0.6.0 numbers (retracted). Unresolved != equivalent.
- README top still carries the 0.7.0 framing; a README rewrite for the product direction needs team approval.
- The long "improved plan" list of 2026-09-23 is in the previous STATE.md (`git show 4de687d:.projectos/STATE.md`).
