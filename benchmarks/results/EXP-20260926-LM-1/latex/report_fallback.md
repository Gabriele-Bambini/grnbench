# EXP-20260926-LM-1: report (text fallback, no TeX engine found)

- experiment: EXP-20260926-LM-1
- type: **confirmatory** (pre-registration = the commit that adds this file, pushed before any of the four models below was run on any dataset of this list)
- owner: LM (agent-run)   issue: #10   registered: 2026-09-26
- datasets (rows): 20
- models: grnboost2, grnboost2_motif, regdiffusion, regdiffusion_motif
- arms: auto, baseline
- code: branch `lm/exp-hpc-motif-models` after merging `main` with PR #11 (issue #6: batches under 100 cells are pooled for scrublet, so the 4 LM-2 datasets that fell back to the heuristic now get scrublet; this changes their processed data relative to LM-2, which is one more reason the comparison with LM-2 is descriptive only). The job copies `GIT_COMMIT` (written by `hpc/hpc.sh sync`) to the results folder as `CODE_COMMIT`; that commit must be clean.
- code commit (run folder): 84f28e1edfdee819cce4e5323c719cc59a32b1fc
- environment: conda env `grnbench` on the cluster, Python 3.11, `requirements/lock-py311.txt` plus the RegDiffusion line of its header (regdiffusion 0.2.3, torch 2.14.0 CPU) via `WITH_TORCH=1 bash hpc/setup_env.sh`.
- machine: Bocconi HPC, partition/QOS `stud`, 1 job, 8 CPUs (8 dataset workers, one thread each), 32 GB. Expected cost: RegDiffusion took 808 s per dataset and arm on the laptop in LM-1; with `--keep-artifacts` both arms are scored, so about 40 RegDiffusion fits plus screening, roughly 3-6 h on 8 single-threaded workers.
- seeds: reference split seed 0 (`screen.split_reference`), cell subsample seed 0 (at most 3,000 cells per context, `MAX_SCORE_CELLS`), scrublet seed 0, RegDiffusion seed 0 (1,000 steps, CPU), GRNBoost2 seed 42 (the adapter defaults); regulators = the 150 TFs with most reference targets in the data (`MAX_TFS`).
- run folder: benchmarks/results/EXP-20260926-LM-1
- generated: 2026-09-27 08:34

## AUPRC ratio, auto pipeline (paired t, Holm over 4 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| grnboost2_motif | grnboost2 | 20 | +0.130 | [+0.108, +0.153] | 3.6e-10 | 1.000 | grnboost2_motif |
| regdiffusion_motif | regdiffusion | 19 | +0.184 | [+0.132, +0.237] | 1.2e-06 | 1.000 | regdiffusion_motif |

## AUPRC ratio, baseline pipeline (paired t, Holm over 4 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| grnboost2_motif | grnboost2 | 20 | +0.129 | [+0.110, +0.149] | 4.8e-11 | 1.000 | grnboost2_motif |
| regdiffusion_motif | regdiffusion | 19 | +0.187 | [+0.138, +0.236] | 3.7e-07 | 1.000 | regdiffusion_motif |

## EPR, auto pipeline (paired t, Holm over 4 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| grnboost2_motif | grnboost2 | 20 | +0.606 | [+0.473, +0.739] | 2.2e-08 | 1.000 | grnboost2_motif |
| regdiffusion_motif | regdiffusion | 19 | -0.003 | [-0.197, +0.191] | 1.000 | 1.000 | unresolved |

## EPR, baseline pipeline (paired t, Holm over 4 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| grnboost2_motif | grnboost2 | 20 | +0.548 | [+0.422, +0.673] | 4.4e-08 | 1.000 | grnboost2_motif |
| regdiffusion_motif | regdiffusion | 19 | +0.006 | [-0.201, +0.214] | 1.000 | 1.000 | unresolved |
