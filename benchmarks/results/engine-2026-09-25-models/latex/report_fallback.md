# EXP-20260925-LM-1: report (text fallback, no TeX engine found)

- experiment: EXP-20260925-LM-1
- type: exploratory. The comparison family was declared in code (`real_bench.DECLARED_PAIRS`) before the run, but there was no committed registration record. Needs a confirmatory rerun before it feeds a headline claim.
- owner: LM (agent run on Luca's laptop)
- datasets (rows): 12
- models: clr, grnboost2, grnboost2_cistarget, grnboost2_motif, motif_prior, pearson_motif, pearson_signed, regdiffusion, regdiffusion_motif, spearman_signed
- arms: auto, baseline
- code: commit 942cdea of Luca's offline history (its tree is part of the first commit of this repo); `CODE_COMMIT` in the results folder
- code commit (run folder): 942cdea61936cbc96bac0a83c344e5cd03f94bcd
- environment: Python 3.11.9 (Windows), requirements/lock-py311.txt + regdiffusion 0.2.3, torch 2.14.0 CPU
- machine: local laptop, 3 workers, 03:28-05:18, 13 datasets, 0 errors
- seeds: reference split seed 0; at most 3,000 cells per context
- run folder: benchmarks/results/engine-2026-09-25-models/test12
- generated: 2026-09-25 18:40

## AUPRC ratio, auto pipeline (paired t, Holm over 14 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 12 | +0.229 | [+0.160, +0.297] | 9.9e-05 | 1.000 | pearson_motif |
| grnboost2_motif | grnboost2 | 12 | +0.143 | [+0.110, +0.177] | 9.6e-06 | 1.000 | grnboost2_motif |
| regdiffusion_motif | regdiffusion | 12 | +0.222 | [+0.136, +0.308] | 0.00085 | 1.000 | regdiffusion_motif |
| grnboost2_cistarget | grnboost2 | 12 | +0.076 | [+0.047, +0.106] | 0.00085 | 1.000 | grnboost2_cistarget |
| pearson_motif | motif_prior | 12 | +0.114 | [-0.029, +0.257] | 0.429 | 1.000 | unresolved |
| grnboost2 | pearson_signed | 12 | -0.221 | [-0.375, -0.067] | 1.000 | 0.046 | pearson_signed |
| regdiffusion | pearson_signed | 12 | +0.113 | [+0.025, +0.201] | 0.074 | 1.000 | unresolved |

## AUPRC ratio, baseline pipeline (paired t, Holm over 14 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 12 | +0.221 | [+0.149, +0.292] | 0.0002 | 1.000 | pearson_motif |
| grnboost2_motif | grnboost2 | 12 | +0.145 | [+0.111, +0.178] | 8.3e-06 | 1.000 | grnboost2_motif |
| regdiffusion_motif | regdiffusion | 12 | +0.226 | [+0.142, +0.310] | 0.00061 | 1.000 | regdiffusion_motif |
| grnboost2_cistarget | grnboost2 | 12 | +0.074 | [+0.045, +0.102] | 0.00074 | 1.000 | grnboost2_cistarget |
| pearson_motif | motif_prior | 12 | +0.105 | [-0.019, +0.230] | 0.356 | 1.000 | unresolved |
| grnboost2 | pearson_signed | 12 | -0.210 | [-0.331, -0.088] | 1.000 | 0.015 | pearson_signed |
| regdiffusion | pearson_signed | 12 | +0.114 | [+0.006, +0.221] | 0.179 | 1.000 | unresolved |

## EPR, auto pipeline (paired t, Holm over 14 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 12 | +0.372 | [+0.063, +0.681] | 0.101 | 1.000 | unresolved |
| grnboost2_motif | grnboost2 | 12 | +0.674 | [+0.511, +0.837] | 1.3e-05 | 1.000 | grnboost2_motif |
| regdiffusion_motif | regdiffusion | 12 | +0.043 | [-0.307, +0.392] | 1.000 | 1.000 | unresolved |
| grnboost2_cistarget | grnboost2 | 12 | +0.941 | [+0.668, +1.213] | 6.9e-05 | 1.000 | grnboost2_cistarget |
| pearson_motif | motif_prior | 12 | +0.648 | [+0.294, +1.002] | 0.010 | 1.000 | pearson_motif |
| grnboost2 | pearson_signed | 12 | -0.934 | [-1.380, -0.489] | 1.000 | 0.004 | pearson_signed |
| regdiffusion | pearson_signed | 12 | +0.772 | [+0.422, +1.122] | 0.003 | 1.000 | regdiffusion |

## EPR, baseline pipeline (paired t, Holm over 14 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 12 | +0.350 | [+0.040, +0.660] | 0.137 | 1.000 | unresolved |
| grnboost2_motif | grnboost2 | 12 | +0.699 | [+0.520, +0.878] | 2.3e-05 | 1.000 | grnboost2_motif |
| regdiffusion_motif | regdiffusion | 12 | -0.137 | [-0.504, +0.231] | 1.000 | 1.000 | unresolved |
| grnboost2_cistarget | grnboost2 | 12 | +1.037 | [+0.768, +1.306] | 2.4e-05 | 1.000 | grnboost2_cistarget |
| pearson_motif | motif_prior | 12 | +0.586 | [+0.263, +0.909] | 0.011 | 1.000 | pearson_motif |
| grnboost2 | pearson_signed | 12 | -0.942 | [-1.351, -0.533] | 1.000 | 0.002 | pearson_signed |
| regdiffusion | pearson_signed | 12 | +0.899 | [+0.543, +1.255] | 0.001 | 1.000 | regdiffusion |
