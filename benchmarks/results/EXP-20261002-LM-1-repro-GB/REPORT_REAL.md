# Engine mass test on real datasets - 20 datasets

- ran end to end: 20
- refused as expected (not raw counts): 0
- errors: 0
- doublet calling fell back to the heuristic (scrublet failed) in 0/20 datasets (baseline or auto arm)

## Detection vs curated metadata

| check | n | correct |
|---|---|---|
| organism | 20 | 100% |
| technology | 20 | 100% |

## Reference recovery (CollecTRI evaluation half, AUPRC ratio; 1 = random)

| model | n | auto > baseline | auto = baseline | mean auto | mean baseline | mean rules |
|---|---|---|---|---|---|---|
| pearson_signed | 20 | 25% | 60% | 1.39 | 1.36 | 1.21 |
| spearman_signed | 20 | 35% | 50% | 1.40 | 1.37 | 1.24 |
| clr | 20 | 25% | 70% | 1.34 | 1.33 | 1.29 |

## Models (CollecTRI evaluation half; unit = dataset)

### AUPRC ratio (1 = random)

| model | n | mean auto | median auto | mean baseline | auto > baseline | s per dataset (auto) |
|---|---|---|---|---|---|---|
| regdiffusion_motif | 20 | 1.664 | 1.645 | 1.637 | 5/20 | 2 |
| pearson_motif | 20 | 1.587 | 1.575 | 1.556 | 4/20 | 2 |
| regdiffusion | 20 | 1.430 | 1.402 | 1.419 | 5/20 | 405 |
| pearson_signed | 20 | 1.384 | 1.378 | 1.359 | 4/20 | 1 |
| grnboost2_motif | 20 | 1.269 | 1.280 | 1.269 | 3/20 | 2 |
| grnboost2 | 20 | 1.122 | 1.130 | 1.129 | 2/20 | 56 |

#### Declared comparisons, AUPRC ratio, auto pipeline (paired t, Holm over 10 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 20 | +0.203 | [+0.147, +0.259] | 1.62e-06 | 1 | pearson_motif |
| grnboost2_motif | grnboost2 | 20 | +0.146 | [+0.119, +0.173] | 3.25e-09 | 1 | grnboost2_motif |
| regdiffusion_motif | regdiffusion | 20 | +0.234 | [+0.172, +0.296] | 9.06e-07 | 1 | regdiffusion_motif |
| grnboost2 | pearson_signed | 20 | -0.262 | [-0.363, -0.161] | 1 | 0.000109 | pearson_signed |
| regdiffusion | pearson_signed | 20 | +0.046 | [-0.041, +0.133] | 0.854 | 1 | unresolved |

#### Declared comparisons, AUPRC ratio, baseline pipeline (paired t, Holm over 10 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 20 | +0.197 | [+0.149, +0.245] | 2.63e-07 | 1 | pearson_motif |
| grnboost2_motif | grnboost2 | 20 | +0.140 | [+0.114, +0.167] | 4.08e-09 | 1 | grnboost2_motif |
| regdiffusion_motif | regdiffusion | 20 | +0.218 | [+0.158, +0.278] | 1.3e-06 | 1 | regdiffusion_motif |
| grnboost2 | pearson_signed | 20 | -0.230 | [-0.330, -0.130] | 1 | 0.000411 | pearson_signed |
| regdiffusion | pearson_signed | 20 | +0.061 | [-0.023, +0.144] | 0.439 | 1 | unresolved |

### EPR (1 = random)

| model | n | mean auto | median auto | mean baseline | auto > baseline | s per dataset (auto) |
|---|---|---|---|---|---|---|
| regdiffusion_motif | 20 | 2.711 | 2.786 | 2.707 | 3/20 | 2 |
| regdiffusion | 20 | 2.643 | 2.662 | 2.625 | 4/20 | 405 |
| pearson_motif | 20 | 2.388 | 2.346 | 2.352 | 5/20 | 2 |
| pearson_signed | 20 | 2.284 | 2.221 | 2.196 | 5/20 | 1 |
| grnboost2_motif | 20 | 1.848 | 1.794 | 1.875 | 4/20 | 2 |
| grnboost2 | 20 | 1.303 | 1.347 | 1.352 | 2/20 | 56 |

#### Declared comparisons, EPR, auto pipeline (paired t, Holm over 10 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 20 | +0.104 | [-0.095, +0.303] | 1 | 1 | unresolved |
| grnboost2_motif | grnboost2 | 20 | +0.544 | [+0.405, +0.684] | 6e-07 | 1 | grnboost2_motif |
| regdiffusion_motif | regdiffusion | 20 | +0.068 | [-0.162, +0.297] | 1 | 1 | unresolved |
| grnboost2 | pearson_signed | 20 | -0.981 | [-1.297, -0.665] | 1 | 1.44e-05 | pearson_signed |
| regdiffusion | pearson_signed | 20 | +0.359 | [+0.034, +0.683] | 0.128 | 1 | unresolved |

#### Declared comparisons, EPR, baseline pipeline (paired t, Holm over 10 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 20 | +0.157 | [-0.008, +0.322] | 0.216 | 1 | unresolved |
| grnboost2_motif | grnboost2 | 20 | +0.523 | [+0.407, +0.639] | 6.74e-08 | 1 | grnboost2_motif |
| regdiffusion_motif | regdiffusion | 20 | +0.082 | [-0.158, +0.322] | 1 | 1 | unresolved |
| grnboost2 | pearson_signed | 20 | -0.843 | [-1.134, -0.553] | 1 | 3.46e-05 | pearson_signed |
| regdiffusion | pearson_signed | 20 | +0.430 | [+0.095, +0.764] | 0.0582 | 1 | unresolved |

### Per dataset, AUPRC ratio (auto pipeline)

| dataset | pearson_signed | pearson_motif | regdiffusion | regdiffusion_motif | grnboost2 | grnboost2_motif |
|---|---|---|---|---|---|---|
| 10x-3-v2__homo-sapiens__breast__44941fdb | 1.66 | 1.89 | 1.86 | 2.12 | 1.11 | 1.28 |
| 10x-3-v2__homo-sapiens__cortex-of-kidney__2fc9c59f | 1.38 | 1.56 | 1.40 | 1.62 | 1.07 | 1.21 |
| 10x-3-v2__homo-sapiens__duodenum__ee195b7d | 1.37 | 1.77 | 1.41 | 1.84 | 1.19 | 1.39 |
| 10x-3-v2__homo-sapiens__entorhinal-cortex__f9ad5649 | 1.05 | 1.29 | 1.04 | 1.32 | 1.02 | 1.11 |
| 10x-3-v2__homo-sapiens__kidney__4c6f9f26 | 1.38 | 1.59 | 1.39 | 1.60 | 1.18 | 1.29 |
| 10x-3-v2__homo-sapiens__peripheral-zone-of-prostate__ea426edb | 1.35 | 1.44 | 1.29 | 1.42 | 1.10 | 1.18 |
| 10x-3-v2__homo-sapiens__superior-frontal-gyrus__1492eb6b | 1.28 | 1.66 | 1.21 | 1.64 | 1.05 | 1.21 |
| 10x-3-v3__homo-sapiens__chest-wall__743473de | 1.38 | 1.55 | 1.40 | 1.57 | 1.18 | 1.35 |
| 10x-3-v3__homo-sapiens__colon__4269074c | 1.38 | 1.52 | 1.40 | 1.56 | 1.18 | 1.31 |
| 10x-3-v3__homo-sapiens__liver__ef0d813e | 1.85 | 1.86 | 1.70 | 1.82 | 1.16 | 1.26 |
| 10x-3-v3__homo-sapiens__respiratory-airway__b25f3834 | 1.54 | 1.76 | 1.80 | 1.94 | 1.15 | 1.29 |
| 10x-3-v3__homo-sapiens__retina__d95ab381 | 1.06 | 1.52 | 1.17 | 1.70 | 1.07 | 1.28 |
| 10x-3-v3__homo-sapiens__telencephalon__6fc37633 | 1.35 | 1.53 | 1.15 | 1.45 | 1.10 | 1.25 |
| 10x-3-v3__mus-musculus__autopod-skin__1398f3a9 | 1.38 | 1.65 | 1.43 | 1.74 | 1.18 | 1.40 |
| 10x-3-v3__mus-musculus__autopod-skin__5bc26721 | 1.81 | 2.13 | 1.68 | 2.07 | 1.16 | 1.47 |
| 10x-3-v3__mus-musculus__ovary__c32b7e3d | 1.80 | 1.80 | 1.65 | 1.73 | 1.21 | 1.32 |
| 10x-3-v3__mus-musculus__primary-motor-cortex__4ec57620 | 1.03 | 1.25 | 1.24 | 1.45 | 1.04 | 1.22 |
| 10x-multiome__homo-sapiens__brain__00f15b43 | 1.13 | 1.22 | 1.19 | 1.35 | 1.06 | 1.14 |
| 10x-multiome__homo-sapiens__brain__5662e542 | 1.03 | 1.15 | 1.68 | 1.70 | 1.05 | 1.12 |
| smart-seq2__mus-musculus__brown-adipose-tissue__8f1bc86b | 1.47 | 1.60 | 1.50 | 1.65 | 1.18 | 1.30 |

## Per dataset

| dataset | truth | detected | ref edges | outcome | Pearson auto/base | CLR auto/base | s |
|---|---|---|---|---|---|---|---|
| cxg__10x-3-v2__homo-sapiens__breast__44941fdb | 10x 3' v2 human | umi human | 20396 | ok, scored | 1.77/1.77 | 1.59/1.58 | 1182.4 |
| cxg__10x-3-v2__homo-sapiens__cortex-of-kidney__2fc9c59f | 10x 3' v2 human | umi human | 27041 | ok, scored | 1.38/1.38 | 1.43/1.33 | 1378.2 |
| cxg__10x-3-v2__homo-sapiens__duodenum__ee195b7d | 10x 3' v2 human | umi human | 26785 | ok, scored | 1.37/1.28 | 1.45/1.45 | 1474.2 |
| cxg__10x-3-v2__homo-sapiens__entorhinal-cortex__f9ad564 | 10x 3' v2 human | umi human | 27002 | ok, scored | 1.05/1.05 | 1.03/1.03 | 1059.1 |
| cxg__10x-3-v2__homo-sapiens__kidney__4c6f9f26 | 10x 3' v2 human | umi human | 27041 | ok, scored | 1.38/1.19 | 1.37/1.28 | 1087.1 |
| cxg__10x-3-v2__homo-sapiens__peripheral-zone-of-prostat | 10x 3' v2 human | umi human | 24590 | ok, scored | 1.35/1.28 | 1.37/1.37 | 1163.0 |
| cxg__10x-3-v2__homo-sapiens__superior-frontal-gyrus__14 | 10x 3' v2 human | umi human | 27002 | ok, scored | 1.38/1.28 | 1.19/1.19 | 1698.6 |
| cxg__10x-3-v3__homo-sapiens__chest-wall__743473de | 10x 3' v3 human | umi human | 26986 | ok, scored | 1.38/1.38 | 1.29/1.23 | 1001.3 |
| cxg__10x-3-v3__homo-sapiens__colon__4269074c | 10x 3' v3 human | umi human | 26981 | ok, scored | 1.38/1.38 | 1.42/1.42 | 1112.7 |
| cxg__10x-3-v3__homo-sapiens__liver__ef0d813e | 10x 3' v3 human | umi human | 26986 | ok, scored | 1.85/1.51 | 1.36/1.28 | 1050.4 |
| cxg__10x-3-v3__homo-sapiens__respiratory-airway__b25f38 | 10x 3' v3 human | umi human | 27058 | ok, scored | 1.54/1.54 | 1.56/1.56 | 1121.9 |
| cxg__10x-3-v3__homo-sapiens__retina__d95ab381 | 10x 3' v3 human | umi human | 27102 | ok, scored | 1.03/1.06 | 1.02/1.15 | 1034.3 |
| cxg__10x-3-v3__homo-sapiens__telencephalon__6fc37633 | 10x 3' v3 human | umi human | 27063 | ok, scored | 1.35/1.35 | 1.11/1.11 | 1299.8 |
| cxg__10x-3-v3__mus-musculus__autopod-skin__1398f3a9 | 10x 3' v3 mouse | umi mouse | 28284 | ok, scored | 1.38/1.38 | 1.28/1.28 | 1248.3 |
| cxg__10x-3-v3__mus-musculus__autopod-skin__5bc26721 | 10x 3' v3 mouse | umi mouse | 28284 | ok, scored | 1.81/1.81 | 1.50/1.50 | 1260.3 |
| cxg__10x-3-v3__mus-musculus__ovary__c32b7e3d | 10x 3' v3 mouse | umi mouse | 23679 | ok, scored | 1.80/1.80 | 1.41/1.41 | 1219.9 |
| cxg__10x-3-v3__mus-musculus__primary-motor-cortex__4ec5 | 10x 3' v3 mouse | umi mouse | 28157 | ok, scored | 1.03/1.05 | 1.18/1.18 | 2451.2 |
| cxg__10x-multiome__homo-sapiens__brain__00f15b43 | 10x multiome human | umi human | 26893 | ok, scored | 1.13/1.18 | 1.17/1.17 | 1214.0 |
| cxg__10x-multiome__homo-sapiens__brain__5662e542 | 10x multiome human | umi human | 26893 | ok, scored | 1.03/1.03 | 1.63/1.63 | 1064.4 |
| cxg__smart-seq2__mus-musculus__brown-adipose-tissue__8f | Smart-seq2 mouse | full_length mouse | 27064 | ok, scored | 1.47/1.47 | 1.45/1.45 | 1073.3 |
