# Claims ledger

One row per sentence-level claim. Status: draft (exploratory evidence) -> confirmed (confirmatory
experiment) -> reproduced (second person reran it). Only `reproduced` rows go into a submission.

| # | claim (as it would read in the paper) | experiments | result files / generated table | status | reproduced by |
|---|---|---|---|---|---|
| C1 | Adding a cisTarget motif prior to Pearson co-expression raises the AUPRC ratio on the CollecTRI evaluation half by 0.23 (95% simultaneous CI 0.16-0.30; n = 12 real datasets; paired two-sided t, Holm over 14 declared hypotheses, p = 1e-4). | EXP-20260925-LM-1 | `paper/tables/T_models_real_EXP-20260925-LM-1.md` | draft | - |
| C2 | A BEELINE-style design with 7 datasets resolves none of 108 pairwise method comparisons after Holm correction. | EXP-20260917-PRE-1 | `benchmarks/results/beeline-2026-09-17-9methods/summary.md` | draft | - |
| C3 | On real data the engine's automatic cleaning ties with a fixed baseline pipeline for every model tested. | EXP-20260925-LM-1, EXP-20260924-PRE-2 | `paper/tables/T_models_real_EXP-20260925-LM-1.md` | draft | - |
