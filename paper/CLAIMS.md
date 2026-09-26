# Claims ledger

One row per sentence-level claim. Status, in order:
- `draft`: exploratory evidence.
- `confirmed`: a confirmatory (pre-registered) experiment supports it.
- `reproduced-computational`: a different agent reran the experiment from its record alone and got
  the same numbers (record section 4). The computation is reproducible, but a human reproduction is
  still pending. This is **not** a second-person reproduction.
- `reproduced`: a second person reran it from the record.

Only `reproduced` rows go into a submission.

| # | claim (as it would read in the paper) | experiments | result files / generated table | status | reproduced by |
|---|---|---|---|---|---|
| C1 | (exploratory, superseded by C4) Adding a cisTarget motif prior to Pearson co-expression raised the AUPRC ratio on the CollecTRI evaluation half by 0.23 (95% CI 0.16-0.30; n = 12 real datasets; paired two-sided t, Holm over 14 declared hypotheses, p = 1e-4). | EXP-20260925-LM-1 (exploratory) | `paper/tables/T_models_real_EXP-20260925-LM-1.md` | draft | - |
| C4 | On 20 fresh real datasets, pre-registered: adding the cisTarget motif prior to Pearson co-expression raises the AUPRC ratio on the CollecTRI evaluation half by 0.18 (95% CI 0.14-0.23; n = 20 datasets; paired two-sided t, Holm over 4 declared hypotheses, p = 1.4e-7); the gain in top-of-list precision (EPR) is not resolved. | EXP-20260925-LM-2 (confirmatory) | `benchmarks/results/EXP-20260925-LM-2/latex/tables/comparisons_auprc_ratio_auto.tex`, `REPORT_REAL.md` | reproduced-computational | independent agent rerun 2026-09-26, identical numbers (record section 4); human: pending (GB) |
| C2 | A BEELINE-style design with 7 datasets resolves none of 108 pairwise method comparisons after Holm correction. | EXP-20260917-PRE-1 | `benchmarks/results/beeline-2026-09-17-9methods/summary.md` | draft | - |
| C3 | On real data the engine's automatic cleaning ties with a fixed baseline pipeline for every model tested. | EXP-20260925-LM-1, EXP-20260924-PRE-2 | `paper/tables/T_models_real_EXP-20260925-LM-1.md` | draft | - |
