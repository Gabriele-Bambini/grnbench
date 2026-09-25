# Engine mass test on real datasets - 13 datasets

- ran end to end: 13
- refused as expected (not raw counts): 0
- errors: 0

## Detection vs curated metadata

| check | n | correct |
|---|---|---|
| organism | 13 | 100% |
| technology | 13 | 92% |

## Reference recovery (CollecTRI evaluation half, AUPRC ratio; 1 = random)

| model | n | auto > baseline | auto = baseline | mean auto | mean baseline | mean rules |
|---|---|---|---|---|---|---|
| pearson_signed | 13 | 23% | 69% | 1.34 | 1.33 | 1.22 |
| spearman_signed | 13 | 31% | 54% | 1.39 | 1.36 | 1.21 |
| clr | 13 | 23% | 62% | 1.33 | 1.33 | 1.25 |

## Models (CollecTRI evaluation half; unit = dataset)

### AUPRC ratio (1 = random)

| model | n | mean auto | median auto | mean baseline | auto > baseline | s per dataset (auto) |
|---|---|---|---|---|---|---|
| regdiffusion_motif | 13 | 1.661 | 1.626 | 1.670 | 3/13 | 3 |
| pearson_motif | 13 | 1.556 | 1.519 | 1.556 | 2/13 | 2 |
| regdiffusion | 13 | 1.442 | 1.399 | 1.441 | 3/13 | 778 |
| motif_prior | 13 | 1.437 | 1.441 | 1.443 | 3/13 | 3 |
| spearman_signed | 13 | 1.371 | 1.414 | 1.358 | 3/13 | 3 |
| pearson_signed | 13 | 1.325 | 1.295 | 1.330 | 3/13 | 2 |
| clr | 13 | 1.306 | 1.304 | 1.332 | 1/13 | 18 |
| grnboost2_motif | 13 | 1.248 | 1.220 | 1.261 | 1/13 | 2 |
| grnboost2_cistarget | 13 | 1.178 | 1.172 | 1.187 | 1/13 | 59 |
| grnboost2 | 13 | 1.099 | 1.100 | 1.110 | 0/13 | 117 |

#### Declared comparisons, AUPRC ratio, auto pipeline (paired t, Holm over 14 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 13 | +0.230 | [+0.168, +0.293] | 2.45e-05 | 1 | pearson_motif |
| grnboost2_motif | grnboost2 | 13 | +0.149 | [+0.116, +0.182] | 3.03e-06 | 1 | grnboost2_motif |
| regdiffusion_motif | regdiffusion | 13 | +0.218 | [+0.140, +0.297] | 0.000321 | 1 | regdiffusion_motif |
| grnboost2_cistarget | grnboost2 | 13 | +0.079 | [+0.051, +0.107] | 0.000277 | 1 | grnboost2_cistarget |
| pearson_motif | motif_prior | 13 | +0.119 | [-0.012, +0.250] | 0.285 | 1 | unresolved |
| grnboost2 | pearson_signed | 13 | -0.226 | [-0.367, -0.085] | 1 | 0.0219 | pearson_signed |
| regdiffusion | pearson_signed | 13 | +0.117 | [+0.037, +0.198] | 0.0365 | 1 | regdiffusion |

#### Declared comparisons, AUPRC ratio, baseline pipeline (paired t, Holm over 14 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 13 | +0.226 | [+0.160, +0.292] | 5.13e-05 | 1 | pearson_motif |
| grnboost2_motif | grnboost2 | 13 | +0.151 | [+0.118, +0.185] | 2.97e-06 | 1 | grnboost2_motif |
| regdiffusion_motif | regdiffusion | 13 | +0.229 | [+0.152, +0.306] | 0.000179 | 1 | regdiffusion_motif |
| grnboost2_cistarget | grnboost2 | 13 | +0.077 | [+0.050, +0.104] | 0.000242 | 1 | grnboost2_cistarget |
| pearson_motif | motif_prior | 13 | +0.112 | [-0.002, +0.227] | 0.213 | 1 | unresolved |
| grnboost2 | pearson_signed | 13 | -0.220 | [-0.333, -0.107] | 1 | 0.00568 | pearson_signed |
| regdiffusion | pearson_signed | 13 | +0.111 | [+0.013, +0.209] | 0.135 | 1 | unresolved |

### EPR (1 = random)

| model | n | mean auto | median auto | mean baseline | auto > baseline | s per dataset (auto) |
|---|---|---|---|---|---|---|
| regdiffusion_motif | 13 | 2.786 | 2.698 | 2.760 | 3/13 | 3 |
| regdiffusion | 13 | 2.735 | 2.860 | 2.831 | 1/13 | 778 |
| pearson_motif | 13 | 2.421 | 2.284 | 2.387 | 2/13 | 2 |
| spearman_signed | 13 | 2.181 | 2.208 | 2.190 | 2/13 | 3 |
| pearson_signed | 13 | 2.053 | 1.941 | 2.040 | 2/13 | 2 |
| grnboost2_cistarget | 13 | 2.035 | 2.006 | 2.057 | 2/13 | 59 |
| clr | 13 | 1.906 | 1.963 | 1.942 | 2/13 | 18 |
| motif_prior | 13 | 1.728 | 1.741 | 1.768 | 2/13 | 3 |
| grnboost2_motif | 13 | 1.709 | 1.702 | 1.769 | 2/13 | 2 |
| grnboost2 | 13 | 1.044 | 0.964 | 1.085 | 4/13 | 117 |

#### Declared comparisons, EPR, auto pipeline (paired t, Holm over 14 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 13 | +0.368 | [+0.087, +0.649] | 0.0659 | 1 | unresolved |
| grnboost2_motif | grnboost2 | 13 | +0.665 | [+0.515, +0.815] | 3.6e-06 | 1 | grnboost2_motif |
| regdiffusion_motif | regdiffusion | 13 | +0.052 | [-0.268, +0.371] | 1 | 1 | unresolved |
| grnboost2_cistarget | grnboost2 | 13 | +0.991 | [+0.720, +1.262] | 2.55e-05 | 1 | grnboost2_cistarget |
| pearson_motif | motif_prior | 13 | +0.693 | [+0.357, +1.030] | 0.00409 | 1 | pearson_motif |
| grnboost2 | pearson_signed | 13 | -1.009 | [-1.446, -0.572] | 1 | 0.00177 | pearson_signed |
| regdiffusion | pearson_signed | 13 | +0.682 | [+0.308, +1.056] | 0.00927 | 1 | regdiffusion |

#### Declared comparisons, EPR, baseline pipeline (paired t, Holm over 14 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 13 | +0.346 | [+0.064, +0.629] | 0.0917 | 1 | unresolved |
| grnboost2_motif | grnboost2 | 13 | +0.684 | [+0.518, +0.851] | 8.12e-06 | 1 | grnboost2_motif |
| regdiffusion_motif | regdiffusion | 13 | -0.072 | [-0.435, +0.292] | 1 | 1 | unresolved |
| grnboost2_cistarget | grnboost2 | 13 | +0.973 | [+0.691, +1.255] | 4.6e-05 | 1 | grnboost2_cistarget |
| pearson_motif | motif_prior | 13 | +0.619 | [+0.316, +0.921] | 0.00436 | 1 | pearson_motif |
| grnboost2 | pearson_signed | 13 | -0.956 | [-1.329, -0.582] | 1 | 0.000727 | pearson_signed |
| regdiffusion | pearson_signed | 13 | +0.791 | [+0.390, +1.191] | 0.00516 | 1 | regdiffusion |

### Per dataset, AUPRC ratio (auto pipeline)

| dataset | pearson_signed | spearman_signed | clr | grnboost2 | regdiffusion | motif_prior | pearson_motif | grnboost2_motif | regdiffusion_motif | grnboost2_cistarget |
|---|---|---|---|---|---|---|---|---|---|---|
| 10x-3-v2__homo-sapiens__superior-frontal-gyru | 1.14 | 1.27 | 1.15 | 1.04 | 1.27 | 1.54 | 1.48 | 1.19 | 1.63 | 1.10 |
| 10x-3-v2__mus-musculus__large-intestine__de4e | 1.26 | 1.43 | 1.40 | 1.10 | 1.54 | 1.28 | 1.45 | 1.22 | 1.67 | 1.19 |
| 10x-3-v3__homo-sapiens__liver__d7476ae2 | 1.97 | 1.94 | 1.62 | 1.07 | 1.94 | 1.44 | 2.05 | 1.18 | 2.07 | 1.12 |
| 10x-3-v3__mus-musculus__cervical-spinal-cord_ | 0.97 | 1.05 | 1.23 | 1.10 | 1.20 | 1.64 | 1.36 | 1.29 | 1.62 | 1.14 |
| 10x-5-v1__homo-sapiens__lung__703f00e6 | 1.30 | 1.44 | 1.45 | 1.11 | 1.65 | 1.26 | 1.39 | 1.21 | 1.62 | 1.26 |
| 10x-5-v2__mus-musculus__caecum__0a1ad9d5 | 1.30 | 1.30 | 1.26 | 1.13 | 1.30 | 1.64 | 1.68 | 1.40 | 1.71 | 1.28 |
| 10x-multiome__homo-sapiens__atrioventricular- | 1.25 | 1.58 | 1.30 | 1.06 | 1.47 | 1.50 | 1.52 | 1.21 | 1.71 | 1.10 |
| cel-seq2__homo-sapiens__left-main-bronchus__0 | 1.40 | 1.41 | 1.26 | 1.10 | 1.40 | 1.36 | 1.61 | 1.22 | 1.58 | 1.21 |
| cel-seq__homo-sapiens__epithelium-of-main-bro | 1.16 | 1.17 | 1.32 | 1.05 | 1.35 | 1.43 | 1.41 | 1.21 | 1.58 | 1.17 |
| patch-seq__mus-musculus__primary-motor-cortex | 1.22 | 1.16 | 1.14 | 1.14 | 1.31 | 1.27 | 1.37 | 1.20 | 1.43 | 1.17 |
| smart-seq2__homo-sapiens__blood__a5d95a42 | 1.30 | 1.11 | 1.14 | 1.08 | 1.30 | 1.53 | 1.59 | 1.22 | 1.63 | 1.12 |
| smart-seq2__mus-musculus__aorta__0380ddce | 1.55 | 1.43 | 1.37 | 1.20 | 1.45 | 1.30 | 1.67 | 1.33 | 1.62 | 1.23 |
| pbmc3k_filtered | 1.40 | 1.52 | 1.35 | 1.11 | 1.57 | 1.48 | 1.65 | 1.33 | 1.74 | 1.23 |

## Per dataset

| dataset | truth | detected | ref edges | outcome | Pearson auto/base | CLR auto/base | s |
|---|---|---|---|---|---|---|---|
| cxg__10x-3-v2__homo-sapiens__superior-frontal-gyrus__bd | 10x 3' v2 human | umi human | 27002 | ok, scored | 1.39/1.39 | 1.17/1.17 | 1438.8 |
| cxg__10x-3-v2__mus-musculus__large-intestine__de4e7a0c | 10x 3' v2 mouse | umi mouse | 26662 | ok, scored | 1.26/1.21 | 1.44/1.44 | 3541.5 |
| cxg__10x-3-v3__homo-sapiens__liver__d7476ae2 | 10x 3' v3 human | umi human | 16898 | ok, scored | 1.97/1.77 | 1.62/1.66 | 1139.2 |
| cxg__10x-3-v3__mus-musculus__cervical-spinal-cord__9fcd | 10x 3' v3 mouse | umi mouse | 28159 | ok, scored | 0.97/0.97 | 1.23/1.23 | 1167.4 |
| cxg__10x-5-v1__homo-sapiens__lung__703f00e6 | 10x 5' v1 human | umi human | 18058 | ok, scored | 1.30/1.30 | 1.45/1.45 | 1495.3 |
| cxg__10x-5-v2__mus-musculus__caecum__0a1ad9d5 | 10x 5' v2 mouse | umi mouse | 23848 | ok, scored | 1.30/1.30 | 1.35/1.26 | 1304.4 |
| cxg__10x-multiome__homo-sapiens__atrioventricular-node_ | 10x multiome human | umi human | 27046 | ok, scored | 1.25/1.25 | 1.30/1.30 | 1475.6 |
| cxg__cel-seq2__homo-sapiens__left-main-bronchus__03da97 | CEL-seq2 human | umi human | 19450 | ok, scored | 1.40/1.40 | 1.26/1.26 | 1260.7 |
| cxg__cel-seq__homo-sapiens__epithelium-of-main-bronchus | CEL-seq human | umi human | 27025 | ok, scored | 1.16/1.13 | 1.32/1.32 | 1471.4 |
| cxg__patch-seq__mus-musculus__primary-motor-cortex__8f9 | Patch-seq mouse | None mouse | 27590 | ok, scored | 1.22/1.22 | 1.14/1.14 | 696.5 |
| cxg__smart-seq2__homo-sapiens__blood__a5d95a42 | Smart-seq2 human | full_length human | 21606 | ok, scored | 1.30/1.30 | 1.16/1.14 | 1386.8 |
| cxg__smart-seq2__mus-musculus__aorta__0380ddce | Smart-seq2 mouse | full_length mouse | 27064 | ok, scored | 1.55/1.55 | 1.27/1.37 | 1491.2 |
| pbmc3k_filtered | umi human | umi human | 26471 | ok, scored | 1.40/1.49 | 1.59/1.59 | 1311.2 |
