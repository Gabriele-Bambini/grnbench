# EXP-20261022-GB-1: report (text fallback, no TeX engine found)

- experiment: EXP-20261022-GB-1
- type: confirmatory (pre-registration = the pushed commit that adds this file, before any dataset of the list below is scored)
- owner: GB (agent-prepared)   issue: #139 (data: #141)   registered: 2026-10-22 and LD's reviews of the registration PR: (1) six pool datasets had been scored by EXP-20261001-LD-1; the exclusion now covers `benchmarks/results/` and every branch (`used_datasets.py`), and the pool was redrawn from the same saved listing (194 datasets kept, those 6 and their studies replaced by 6 others); (2) the scope is the registered pool, not every scRNA-seq dataset, and subgroup breakdowns are declared below; (3) data provenance is an experiment-specific snapshot checked with `--require-all`; (4) data and large outputs live in a dedicated data folder outside the code checkout; (5) T3 of EXP-20261003-GB-1 is added as a secondary confirmatory family (LM's suggestion). A partial download made on 2026-10-23 under the first version (4.5 GB of exploratory files, nothing scored) is kept; files not in the redrawn pool are deleted before the data check.
- datasets (rows): 120
- models: pearson_motif, pearson_signed
- arms: auto, baseline
- code: branch `gb/pool-200`, synced to a separate cluster checkout `~/grnbench-pool200` with `hpc/hpc.sh sync`; the job copies `GIT_COMMIT` into the results as `CODE_COMMIT`. The commit that runs is the one that adds the data snapshot (below), pushed before any model runs; `real_bench.py` and the engine are `main`'s (6917c88), the branch adds only the selection tools and records.
- code commit (run folder): 7e4bcb69f1d1c31bc6fbe85820a454cae84ddd68
- environment: conda env `grnbench` on the cluster, Python 3.11, `requirements/lock-py311.txt`; no extra packages (no torch).
- machine: Bocconi HPC, partition/QOS `stud` (1 running + 1 queued job per user), 8 CPUs, 24 h per job.
- seeds: reference split seed 0 (`screen.split_reference`), cell subsample seed 0 (at most 3,000 cells per context, `engine_bench.MAX_SCORE_CELLS`), scrublet seed 0; regulators = the 150 TFs with most reference targets in the data (`MAX_TFS`). The same settings as C4.
- run folder: EXP-20261022-GB-1
- generated: 2026-10-28 18:26

## AUPRC ratio, auto pipeline (paired t, Holm over 2 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 119 | +0.206 | [+0.183, +0.229] | 1.6e-35 | 1.000 | pearson_motif |

## AUPRC ratio, baseline pipeline (paired t, Holm over 2 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 119 | +0.202 | [+0.181, +0.223] | 9.2e-38 | 1.000 | pearson_motif |

## EPR, auto pipeline (paired t, Holm over 2 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 119 | +0.149 | [+0.067, +0.230] | 0.00042 | 1.000 | pearson_motif |

## EPR, baseline pipeline (paired t, Holm over 2 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 119 | +0.199 | [+0.131, +0.267] | 5.3e-08 | 1.000 | pearson_motif |
