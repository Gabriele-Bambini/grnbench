# State - grnbench
Updated: 2026-09-25 by claude for LM (repo moved to GitHub; registry and paper standard added)

## Focus
Product (D-012/D-013): universal, context-aware GRN data-cleaning engine + interface + mounted models + comparison, to be published as a tool with a paper. Headline claims from real data only (D-015), through the experiment registry (D-017).

## Done recently
- 2026-09-25: canonical repo is now github.com/elle13eth/grnbench (private), squashed sanitized snapshot; Drive folder frozen (D-016). Experiment registry `experiments/`, `paper/` standard, `docs/HPC.md`, `hpc/hpc.sh`, pinned `requirements/lock-py311.txt`, `experiments/datasets.tsv`.
- 2026-09-25 (EXP-20260925-LM-1, exploratory): motif prior (cisTarget) lifts every model on 12 real CELLxGENE datasets; pearson_motif vs pearson +0.23 AUPRC ratio [0.16, 0.30], p_Holm 1e-4; regdiffusion_motif best (1.65); GRNBoost2 at random on EPR; auto vs baseline cleaning still ties. `docs/models_that_work.md`.
- 2026-09-24: engine mass tests v1 and v2-partial (EXP-20260924-PRE-1/2): synthetic auto >> baseline, real ~ties, rules alone worse on real.
- 0.7.0 and earlier: see CHANGELOG.md and EXP-20260917-PRE-1 (BEELINE: 0/108 comparisons resolved).

## In progress
- Nothing running. The v2 runs on the collaborator PC were abandoned (partial results archived).
- Known defect: scikit-image missing from the lock, so scrublet falls back to the heuristic doublet caller (test_drop_doublets_uses_scrublet... fails; CI deselects it). Likely affected EXP-20260925-LM-1.

## Next steps (engine validation)
1. Add scikit-image to pyproject + lock; confirm scrublet runs; check the 2026-09-25 manifests for the fallback.
2. First confirmatory experiment on HPC (docs/HPC.md): synthetic ~490 + real datasets not used before, pre-registered; one job per run (QOS: 1 running job, 8 CPUs).
3. Recalibrate the `context_scope` rule in selector.py with bench evidence (pooled for generic regulons; per_cell_type when the goal is a cell-type network; add a "goal" input); regenerate the KB doc.
4. Simulator realism: partially cell-type-specific networks, then rerun synthetic.
5. GUI: per-model best pipeline from screening, knee cell calling, mito-gate decision, motif models.

## Next steps (improved plan, proposed to user 2026-09-23)
0. (NEXT) Define "best": measurable objective + validation bench (known-truth data with injected corruptions: doublets, ambient, damaged cells, batch, dropout) - the selector must beat defaults; this is also the paper's main figure.
1. [done v1] Stage map extended: empty droplets/ambient, gene IDs/orthologs, samples/donors, cell context (per cell type vs pooled), confounders, pseudotime, per-model last-mile inputs.
2. Decision points ranked by impact x uncertainty.
3. [partial: 35 native options; still planned: ensembl_to_symbol (mygene), orthologs, knee cell calling, miQC, per_batch_consensus, cell-cycle regression, DPT, CollecTRI; external R/GPU tools need adapters] Options implemented (missing deps listed with reason).
4. Dataset profile incl. structure (cell types, batch-biology confounding, composition).
5. [done v1] Knowledge base with evidence level per option (GRN-benchmarked / general scRNA / consensus / none), versioned, cited.
6. [selector done; screening of candidates TODO] Selector: best pipeline + alternatives only on uncertain high-impact points; screen with cheap models, confirm with full models.
7. [done in screen.py + GUI benchmark] Anti-circularity: prior network used to pick the pipeline != network used to compare models (split).
8. Model contract declares its input (counts/log/discretized/pseudotime); user models same treatment, containerized.
9. Publication: positioning vs BEAR-GRN / geneRNIB / BEELINE, docs, pip/conda/Docker, "why this pipeline" report.

## Blockers / open questions
- Publication route (JOSS first vs Bioinformatics App Note vs a Nature-family journal) - team decision.
- Licence and copyright holder (D-018, proposed: keep MIT).
- Third team member's GitHub username (to invite).

## Notes for the next agent
- Rules: AGENTS.md. Runs: experiments/. Paper: paper/. HPC: docs/HPC.md. Onboarding: docs/START_HERE_AGENT_PROMPT.md.
- MAP.md predates the engine/GUI modules (context, kb, profile, selector, autoprep, screen, demo, ui/): remap before relying on it.
- README top still carries the 0.7.0 "ranking is not real" framing below the 0.8 section; a README rewrite for the product direction needs team approval.
- Never quote the 0.6.0 numbers (80% resolved, 47% contradictions): retracted. Unresolved != equivalent.
