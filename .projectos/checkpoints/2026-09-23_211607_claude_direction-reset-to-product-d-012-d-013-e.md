---
time: 2026-09-23T21:16:07+02:00
agent: claude
kind: manual
reason: milestone
summary: Direction reset to product (D-012/D-013); engine v0 (context, autoprep, demo, per-cell filters); improved plan in STATE
usage: ctx 274k tok
machine: <host>
---
# Checkpoint - Direction reset to product (D-012/D-013); engine v0 (context, autoprep, demo, per-cell filters); improved plan in STATE

## Files edited since previous checkpoint (5)
- .projectos/PURPOSE.md (Write)
- src/grnbench/context.py (Edit)
- tests/test_context.py (Write)
- src/grnbench/autoprep.py (Write)
- src/grnbench/demo.py (Write)

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
- Engine v0 written this session (untested beyond smoke run): src/grnbench/context.py (+tests/test_context.py, 11 pass), autoprep.py (MAD thresholds, planner, CLI argv), demo.py (synthetic data with planted GRN), per-cell upper filters in qc/pipeline/cli (max genes/counts/%mito).

#### Next steps (improved plan, proposed to user 2026-09-23)
0. Define "best": measurable objective + validation bench (known-truth data with injected corruptions: doublets, ambient, damaged cells, batch, dropout) - the selector must beat defaults; this is also the paper's main figure.
1. Stage map extended: empty droplets/ambient, gene IDs/orthologs, samples/donors, cell context (per cell type vs pooled), confounders, pseudotime, per-model last-mile inputs.
2. Decision points ranked by impact x uncertainty.
3. Options implemented (missing deps listed with reason).
4. Dataset profile incl. structure (cell types, batch-biology confounding, composition).
5. Knowledge base with evidence level per option (GRN-benchmarked / general scRNA / consensus / none), versioned, cited.
6. Selector: best pipeline + alternatives only on uncertain high-impact points; screen with cheap models, confirm with full models.
7. Anti-circularity: prior network used to pick the pipeline != network used to compare models (split).
8. Model contract declares its input (counts/log/discretized/pseudotime); user models same treatment, containerized.
9. Publication: positioning vs BEAR-GRN / geneRNIB / BEELINE, docs, pip/conda/Docker, "why this pipeline" report.

#### Blockers / open questions
- User to confirm the improved plan order (0 -> 9) and the venue (Bioinformatics App Note vs JOSS).

#### Notes for the next agent
- This machine (collab-pc): grnbench not installed; run tests with `PYTHONPATH=src` (python 3.12 has numpy/pandas/scipy/anndata/scanpy). Suite is slow on Google Drive.
- Read `.projectos/MAP.md` for code layout. DIRECTION: product/tool (D-012, D-013). Do NOT steer toward the audit's metascience plan (D-011 rejected; user was angry about it).
- Never quote the 0.6.0 numbers (80% resolved, 47% contradictions): retracted.
- Unresolved != equivalent: the finding is that the design cannot rank, not that methods are equal.
- ProjectOS kit is vendored in `.projectos/kit/` (`python .projectos/kit/install.py` on a new machine).

