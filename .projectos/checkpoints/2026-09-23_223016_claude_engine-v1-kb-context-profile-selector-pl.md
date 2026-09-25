---
time: 2026-09-23T22:30:16+02:00
agent: claude
kind: manual
reason: milestone
summary: Engine v1 (KB, context, profile, selector, planner, screening, new strategies) + Streamlit GUI v0 verified end-to-end on demo
usage: ctx 510k tok
machine: <host>
---
# Checkpoint - Engine v1 (KB, context, profile, selector, planner, screening, new strategies) + Streamlit GUI v0 verified end-to-end...

## Files edited since previous checkpoint (10)
- src/grnbench/kb.py (Write)
- src/grnbench/profile.py (Write)
- src/grnbench/selector.py (Write)
- src/grnbench/autoprep.py (Edit)
- tests/test_engine.py (Write)
- src/grnbench/screen.py (Edit)
- src/grnbench/ui/logic.py (Write)
- src/grnbench/ui/app.py (Write)
- tests/test_ui.py (Write)
- .claude/launch.json (Write)

## STATE.md snapshot

#### State - grnbench
Updated: 2026-09-23 by claude (G0 audit-fixture regressions; A-008 closed)

#### Focus
Product (D-012/D-013): universal, context-aware GRN data-cleaning engine + interface + mounted models + comparison, to be published as a tool. Engine = stages/decision points -> all options implemented -> dataset profile -> literature KB with evidence levels -> rule selector -> tested candidate pipelines.

#### Done recently
- G0 (2026-09-23): every audit/2026-09-13 eval_findings case now mapped to a regression on its exact input; new tests/test_audit_fixtures.py (coverage guard + null_support and finite_random_ap_baseline gaps), 3/3 pass. MI fixtures already exact in tests/test_inference_reliability.py. A-008 degree null verified fixed and closed.
- G0 suite (2026-09-23, this machine, Python 3.12, `PYTHONPATH=src`): 290 passed, 3 skipped, 1 env-dependent failure fixed - tests/test_infer.py::test_list_methods_reports_unavailability assumed java on PATH; now accepts either documented aracne3 prerequisite (java / ARACNE3_JAR).
- 0.7.0 (2026-09-15): retracted the 0.6.0 edge-level significance layer (pseudo-replication, A-003); added `compare_units` / `grnbench compare` (paired t, Holm, Bonferroni intervals) and backtest `unit_comparison` + `regime_check`; added `eval.universe` fixed candidate universe (fixes A-009); gold thinning now fixed-k scenarios (fixes A-007); MI binning by first rank (fixes A-010); BEELINE TF selection pinned to official script SHA-256 (A-013).
- Real run: 9 methods x 7 BEELINE datasets, datasets as units, Holm: **0/108 EPR, 0/108 AUPRC-ratio resolved**; point-estimate winner changes with the reference on 7/7 datasets. Archived `benchmarks/results/beeline-2026-09-17-9methods/`.
- Score-level calibration: confirmation run (1000 reps) passes all iid scenarios; pseudoreplication fails as designed (A-011) -> G3 met at score level.
- Preprocessing fixes in 0.7.0 (0-cell output refused, strategy warnings captured in manifest, explicit 0 thresholds honored, label/group imbalance split, NaN label stratum, duplicate barcodes).

#### In progress
- Nothing running. Engine v1 + GUI v0 done and verified (2026-09-23):
  - Engine: kb.py (14 stages / 15 decisions / 56 options / 51 refs, evidence levels; docs/knowledge/preprocessing_kb.md generated, `grnbench kb --write-docs`), context.py, profile.py, selector.py, autoprep.py (plan + execute_plan per context), screen.py (candidates run + scored; reference split selection/evaluation halves), demo.py; strategies metacells + cpm_log1p_center_batch; run_prep subset (after doublets), exclude_qc_genes, per_batch_bounds, upper per-cell filters; CLI flags; `grnbench kb|demo|ui`.
  - GUI: src/grnbench/ui/app.py (Streamlit, 7 steps) + ui/logic.py; walked end-to-end in the browser on the demo (analyse -> pipeline -> compare -> adopt center_per_batch -> clean 3 contexts -> 4 models -> benchmark on evaluation half). Launch: `.claude/launch.json` config grnbench-ui (PYTHONPATH=src) or `grnbench ui` once installed.
  - Tests: full suite 325 pass / 3 skip before GUI; tests/test_ui.py + test_engine.py + test_context.py + CLI tests 58 pass after.
  - README "New in 0.8" section, CHANGELOG Unreleased, pyproject extra [ui].

#### Next steps (improved plan, proposed to user 2026-09-23)
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

#### Blockers / open questions
- User to confirm the improved plan order (0 -> 9) and the venue (Bioinformatics App Note vs JOSS).

#### Notes for the next agent
- MAP.md predates the engine/GUI modules (context, kb, profile, selector, autoprep, screen, demo, ui/): run an incremental remap (pos-mapper, haiku) before relying on it.
- README top still carries the 0.7.0 "ranking is not real" framing below the new 0.8 section; full README rewrite for the product direction needs user approval.
- This machine (collab-pc): grnbench not installed; run tests with `PYTHONPATH=src` (python 3.12 has numpy/pandas/scipy/anndata/scanpy). Suite is slow on Google Drive.
- Read `.projectos/MAP.md` for code layout. DIRECTION: product/tool (D-012, D-013). Do NOT steer toward the audit's metascience plan (D-011 rejected; user was angry about it).
- Never quote the 0.6.0 numbers (80% resolved, 47% contradictions): retracted.
- Unresolved != equivalent: the finding is that the design cannot rank, not that methods are equal.
- ProjectOS kit is vendored in `.projectos/kit/` (`python .projectos/kit/install.py` on a new machine).

