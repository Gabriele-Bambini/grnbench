# EXP-20261001-LD-1: report (text fallback, no TeX engine found)

- experiment: EXP-20261001-LD-1
- type: exploratory. This run validates new simulator/runner machinery and estimates stress curves; it cannot support a headline or paper claim without a separately registered confirmatory rerun
- owner: LD (agent-run)   issue: #23   registered: 2026-10-01
- datasets (rows): 40
- models: clr, grnboost2, pearson_signed, regdiffusion
- arms: auto, baseline
- code: clean commit on branch `ld/contextual-stress-bench`; `hpc/hpc.sh sync` writes its hash to `GIT_COMMIT`, and each job copies that into `CODE_COMMIT`. The task SHA256 and score-affecting options are pinned separately in `tasks.json` and `RUN_CONFIG.json`.
- code commit (run folder): 9743907d659510c4b63d53473adae7cc7bfc79b0
- environment: Bocconi HPC conda env `grnbench`, Python 3.11, `requirements/lock-py311.txt` plus its pinned RegDiffusion/CPU-torch header line via `WITH_TORCH=1 bash hpc/setup_env.sh`.
- machine: Bocconi HPC, partition/QOS `stud`, one job at a time, 8 CPUs, 32 GB, 24 h. Work files use `~/gb_tmp`, not compute-node `/tmp`; resubmission uses the identical command and output folder.
- seeds: scenario seeds 0-4; union-truth split seed = scenario seed; score-cell subsample seed 0; GRNBoost2 seed 42; RegDiffusion seed 0. Shared edges are split before context expansion, so the same shared edge can never enter selection for one type and evaluation for another.
- run folder: benchmarks/results/EXP-20261001-LD-1/real
- generated: 2026-10-02 14:31

## AUPRC ratio, auto pipeline (paired t, Holm over 4 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| grnboost2 | pearson_signed | 24 | -0.225 | [-0.308, -0.143] | 1.000 | 1.9e-05 | pearson_signed |
| regdiffusion | pearson_signed | 24 | +0.068 | [+0.010, +0.126] | 0.036 | 1.000 | regdiffusion |

## AUPRC ratio, baseline pipeline (paired t, Holm over 4 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| grnboost2 | pearson_signed | 23 | -0.220 | [-0.291, -0.149] | 1.000 | 3.5e-06 | pearson_signed |
| regdiffusion | pearson_signed | 23 | +0.071 | [+0.010, +0.133] | 0.038 | 1.000 | regdiffusion |

## EPR, auto pipeline (paired t, Holm over 4 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| grnboost2 | pearson_signed | 24 | -0.802 | [-1.134, -0.471] | 1.000 | 9.3e-05 | pearson_signed |
| regdiffusion | pearson_signed | 24 | +0.582 | [+0.340, +0.825] | 9.3e-05 | 1.000 | regdiffusion |

## EPR, baseline pipeline (paired t, Holm over 4 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| grnboost2 | pearson_signed | 23 | -0.836 | [-1.079, -0.594] | 1.000 | 7.2e-07 | pearson_signed |
| regdiffusion | pearson_signed | 23 | +0.577 | [+0.297, +0.858] | 0.00047 | 1.000 | regdiffusion |
