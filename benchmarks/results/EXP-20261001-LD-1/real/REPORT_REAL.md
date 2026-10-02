# Engine mass test on real datasets - 40 planned datasets

- completed without a top-level error: 33
- with a requested-model reference score: 24
- completed without a requested-model score: 9
- refused as expected (not raw counts): 7
- permanent errors: 0
- missing or retryable: 0
- doublet calling fell back to the heuristic (scrublet failed) in 0/40 datasets (baseline or auto arm)

## Detection vs curated metadata

| check | n | correct |
|---|---|---|
| cell_calling_best | 1 | 100% |
| organism | 33 | 79% |
| technology | 33 | 97% |
| time_course | 1 | 100% |

## Reference recovery (CollecTRI evaluation half, AUPRC ratio; 1 = random)

| model | n | auto > baseline | auto = baseline | mean auto | mean baseline | mean rules |
|---|---|---|---|---|---|---|
| pearson_signed | 23 | 13% | 78% | 1.33 | 1.32 | 1.20 |
| spearman_signed | 23 | 22% | 70% | 1.40 | 1.38 | 1.22 |
| clr | 23 | 26% | 52% | 1.32 | 1.32 | 1.26 |

## Models (CollecTRI evaluation half; unit = dataset)

### AUPRC ratio (1 = random)

| model | n | mean auto | median auto | mean baseline | auto > baseline | s per dataset (auto) |
|---|---|---|---|---|---|---|
| regdiffusion | 24 | 1.395 | 1.350 | 1.395 | 4/23 | 414 |
| pearson_signed | 24 | 1.326 | 1.278 | 1.324 | 3/23 | 1 |
| clr | 24 | 1.302 | 1.305 | 1.317 | 3/23 | 15 |
| grnboost2 | 24 | 1.101 | 1.100 | 1.104 | 2/23 | 62 |

#### Declared comparisons, AUPRC ratio, auto pipeline (paired t, Holm over 4 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| grnboost2 | pearson_signed | 24 | -0.225 | [-0.308, -0.143] | 1 | 1.92e-05 | pearson_signed |
| regdiffusion | pearson_signed | 24 | +0.068 | [+0.010, +0.126] | 0.0357 | 1 | regdiffusion |

#### Declared comparisons, AUPRC ratio, baseline pipeline (paired t, Holm over 4 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| grnboost2 | pearson_signed | 23 | -0.220 | [-0.291, -0.149] | 1 | 3.48e-06 | pearson_signed |
| regdiffusion | pearson_signed | 23 | +0.071 | [+0.010, +0.133] | 0.0382 | 1 | regdiffusion |

### EPR (1 = random)

| model | n | mean auto | median auto | mean baseline | auto > baseline | s per dataset (auto) |
|---|---|---|---|---|---|---|
| regdiffusion | 24 | 2.543 | 2.682 | 2.623 | 3/23 | 414 |
| pearson_signed | 24 | 1.961 | 2.030 | 2.046 | 2/23 | 1 |
| clr | 24 | 1.822 | 1.978 | 1.899 | 3/23 | 15 |
| grnboost2 | 24 | 1.158 | 1.217 | 1.210 | 3/23 | 62 |

#### Declared comparisons, EPR, auto pipeline (paired t, Holm over 4 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| grnboost2 | pearson_signed | 24 | -0.802 | [-1.134, -0.471] | 1 | 9.27e-05 | pearson_signed |
| regdiffusion | pearson_signed | 24 | +0.582 | [+0.340, +0.825] | 9.27e-05 | 1 | regdiffusion |

#### Declared comparisons, EPR, baseline pipeline (paired t, Holm over 4 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| grnboost2 | pearson_signed | 23 | -0.836 | [-1.079, -0.594] | 1 | 7.19e-07 | pearson_signed |
| regdiffusion | pearson_signed | 23 | +0.577 | [+0.297, +0.858] | 0.000469 | 1 | regdiffusion |

### Per dataset, AUPRC ratio (auto pipeline)

| dataset | pearson_signed | clr | grnboost2 | regdiffusion |
|---|---|---|---|---|
| beeline__GSD-2000-1 | - | - | - | - |
| beeline__GSD-2000-10 | - | - | - | - |
| beeline__GSD-2000-2 | - | - | - | - |
| beeline__GSD-2000-3 | - | - | - | - |
| beeline__GSD-2000-4 | - | - | - | - |
| beeline__GSD-2000-5 | - | - | - | - |
| beeline__GSD-2000-6 | - | - | - | - |
| 10x-3-v2__homo-sapiens__entorhinal-cortex__9f1049ac | 1.21 | 1.13 | 1.04 | 1.21 |
| 10x-3-v2__homo-sapiens__superior-frontal-gyrus__bdacc907 | 1.39 | 1.17 | 1.07 | 1.35 |
| 10x-3-v2__microcebus-murinus__hypothalamus-pituitary-axis__62ec8c21 | - | - | - | - |
| 10x-3-v2__microcebus-murinus__uterus__f69adac8 | - | - | - | - |
| 10x-3-v2__mus-musculus__large-intestine__de4e7a0c | 1.26 | 1.37 | 1.13 | 1.54 |
| 10x-3-v2__mus-musculus__spinal-cord__989a76f9 | 1.25 | 1.12 | 1.07 | 1.15 |
| 10x-3-v3__callithrix-jacchus__hippocampal-formation__19c78a83 | - | - | - | - |
| 10x-3-v3__callithrix-jacchus__primary-motor-cortex__787138a6 | - | - | - | - |
| 10x-3-v3__homo-sapiens__liver__d7476ae2 | 2.00 | 1.62 | 1.08 | 1.94 |
| 10x-3-v3__homo-sapiens__prostate-gland__182f8181 | 1.18 | 1.36 | 1.08 | 1.32 |
| 10x-3-v3__macaca-mulatta__primary-motor-cortex__1f59341c | - | - | - | - |
| 10x-3-v3__mus-musculus__cervical-spinal-cord__9fcd3812 | 0.97 | 1.23 | 1.10 | 1.20 |
| 10x-3-v3__mus-musculus__diencephalon__6de332e1 | 1.21 | 1.20 | 1.15 | 1.16 |
| 10x-5-v1__homo-sapiens__lung__703f00e6 | 1.30 | 1.45 | 1.11 | 1.65 |
| 10x-5-v1__homo-sapiens__lung__cdefb878 | 1.23 | 1.33 | 1.15 | 1.40 |
| 10x-5-v2__homo-sapiens__blood__ee782eb7 | 1.66 | 1.68 | 1.13 | 1.66 |
| 10x-5-v2__mus-musculus__caecum__0a1ad9d5 | 1.30 | 1.26 | 1.13 | 1.30 |
| 10x-5-v2__mus-musculus__caecum__95dfddd0 | 1.36 | 1.27 | 1.16 | 1.35 |
| 10x-gene-expression-flex__homo-sapiens__blood__95498b77 | 1.47 | 1.31 | 1.11 | 1.50 |
| 10x-multiome__homo-sapiens__atrioventricular-node__680029ee | 1.25 | 1.30 | 1.07 | 1.47 |
| bd-rhapsody-targeted-mrna__homo-sapiens__blood__d3566d6a | 1.37 | 1.10 | 1.17 | 1.21 |
| cel-seq2__homo-sapiens__left-main-bronchus__03da97e9 | 1.40 | 1.26 | 1.10 | 1.40 |
| cel-seq__homo-sapiens__epithelium-of-main-bronchus__61d327d1 | 1.16 | 1.32 | 1.05 | 1.35 |
| modified-strt-seq__homo-sapiens__interventricular-septum__f3ee7613 | 1.32 | 1.26 | 1.04 | 1.34 |
| patch-seq__mus-musculus__primary-motor-cortex__8f98c236 | 1.22 | 1.14 | 1.13 | 1.31 |
| sci-rna-seq3__danio-rerio__endocrine-pancreas__e2309842 | - | - | - | - |
| sci-rna-seq3__macaca-mulatta__brodmann-1909-area-19__885d325f | - | - | - | - |
| sci-rna-seq3__mus-musculus__embryo__0bcda669 | 1.24 | 1.31 | 1.06 | 1.31 |
| smart-seq2__homo-sapiens__islet-of-langerhans__66d15835 | 1.21 | 1.31 | 1.06 | 1.27 |
| smart-seq2__mus-musculus__aorta__0380ddce | 1.55 | 1.37 | 1.20 | 1.45 |
| paul15 | - | - | - | - |
| pbmc3k_filtered | 1.31 | 1.36 | 1.07 | 1.60 |
| pbmc3k_raw | - | - | - | - |

## Per dataset

| dataset | truth | detected | ref edges | outcome | Pearson auto/base | CLR auto/base | s |
|---|---|---|---|---|---|---|---|
| beeline__GSD-2000-1 |   | nan nan | nan | refused: not counts | - | - | 0.2 |
| beeline__GSD-2000-10 |   | nan nan | nan | refused: not counts | - | - | 0.2 |
| beeline__GSD-2000-2 |   | nan nan | nan | refused: not counts | - | - | 0.2 |
| beeline__GSD-2000-3 |   | nan nan | nan | refused: not counts | - | - | 0.2 |
| beeline__GSD-2000-4 |   | nan nan | nan | refused: not counts | - | - | 0.2 |
| beeline__GSD-2000-5 |   | nan nan | nan | refused: not counts | - | - | 0.2 |
| beeline__GSD-2000-6 |   | nan nan | nan | refused: not counts | - | - | 0.2 |
| cxg__10x-3-v2__homo-sapiens__entorhinal-cortex__9f1049a | 10x 3' v2 human | umi human | 27002.0 | ok, scored | 1.21/1.21 | 1.13/1.13 | 1146.9 |
| cxg__10x-3-v2__homo-sapiens__superior-frontal-gyrus__bd | 10x 3' v2 human | umi human | 27002.0 | ok, scored | 1.39/1.39 | 1.17/1.17 | 1062.0 |
| cxg__10x-3-v2__microcebus-murinus__hypothalamus-pituita | 10x 3' v2 Microcebus murinus | umi mouse | 0.0 | completed, unscored (2 contexts; no reference edges) | - | - | 21.0 |
| cxg__10x-3-v2__microcebus-murinus__uterus__f69adac8 | 10x 3' v2 Microcebus murinus | umi mouse | 0.0 | completed, unscored (3 contexts; no reference edges) | - | - | 25.3 |
| cxg__10x-3-v2__mus-musculus__large-intestine__de4e7a0c | 10x 3' v2 mouse | umi mouse | 26662.0 | ok, scored | 1.26/1.21 | 1.44/1.44 | 1818.6 |
| cxg__10x-3-v2__mus-musculus__spinal-cord__989a76f9 | 10x 3' v2 mouse | umi mouse | 24615.0 | ok, scored | 1.25/1.25 | 1.09/1.12 | 1165.7 |
| cxg__10x-3-v3__callithrix-jacchus__hippocampal-formatio | 10x 3' v3 Callithrix jacchus | umi None | 0.0 | completed, unscored (4 contexts; no reference edges) | - | - | 45.1 |
| cxg__10x-3-v3__callithrix-jacchus__primary-motor-cortex | 10x 3' v3 Callithrix jacchus | umi human | 0.0 | completed, unscored (4 contexts; no reference edges) | - | - | 47.3 |
| cxg__10x-3-v3__homo-sapiens__liver__d7476ae2 | 10x 3' v3 human | umi human | 16898.0 | ok, scored | 2.00/1.77 | 1.62/1.66 | 955.3 |
| cxg__10x-3-v3__homo-sapiens__prostate-gland__182f8181 | 10x 3' v3 human | umi human | 27063.0 | ok, scored | 1.18/1.18 | 1.36/1.36 | 979.5 |
| cxg__10x-3-v3__macaca-mulatta__primary-motor-cortex__1f | 10x 3' v3 Macaca mulatta | umi human | 0.0 | completed, unscored (4 contexts; no reference edges) | - | - | 72.0 |
| cxg__10x-3-v3__mus-musculus__cervical-spinal-cord__9fcd | 10x 3' v3 mouse | umi mouse | 28159.0 | ok, scored | 0.97/0.97 | 1.23/1.23 | 1014.8 |
| cxg__10x-3-v3__mus-musculus__diencephalon__6de332e1 | 10x 3' v3 mouse | umi mouse | 28191.0 | ok, scored | 1.21/1.21 | 1.20/1.22 | 1096.1 |
| cxg__10x-5-v1__homo-sapiens__lung__703f00e6 | 10x 5' v1 human | umi human | 18058.0 | ok, scored | 1.30/1.30 | 1.45/1.45 | 947.3 |
| cxg__10x-5-v1__homo-sapiens__lung__cdefb878 | 10x 5' v1 human | umi human | 21785.0 | ok, scored | 1.24/1.24 | 1.33/1.29 | 1171.0 |
| cxg__10x-5-v2__homo-sapiens__blood__ee782eb7 | 10x 5' v2 human | umi human | 22729.0 | ok, scored | 1.66/1.66 | 1.80/1.68 | 1540.5 |
| cxg__10x-5-v2__mus-musculus__caecum__0a1ad9d5 | 10x 5' v2 mouse | umi mouse | 23848.0 | ok, scored | 1.30/1.30 | 1.38/1.26 | 1132.3 |
| cxg__10x-5-v2__mus-musculus__caecum__95dfddd0 | 10x 5' v2 mouse | umi mouse | 23848.0 | ok, scored | 1.36/1.36 | 1.38/1.27 | 1156.8 |
| cxg__10x-gene-expression-flex__homo-sapiens__blood__954 | 10x gene expression flex human | umi human | 25951.0 | ok, scored | 1.47/1.47 | 1.31/1.31 | 1676.0 |
| cxg__10x-multiome__homo-sapiens__atrioventricular-node_ | 10x multiome human | umi human | 27046.0 | ok, scored | 1.25/1.25 | 1.30/1.30 | 963.4 |
| cxg__bd-rhapsody-targeted-mrna__homo-sapiens__blood__d3 | BD Rhapsody Targeted mRNA human | umi human | 250.0 | ok, scored | 1.37/nan | 1.10/nan | 569.2 |
| cxg__cel-seq2__homo-sapiens__left-main-bronchus__03da97 | CEL-seq2 human | umi human | 19450.0 | ok, scored | 1.40/1.40 | 1.26/1.26 | 1030.7 |
| cxg__cel-seq__homo-sapiens__epithelium-of-main-bronchus | CEL-seq human | umi human | 27025.0 | ok, scored | 1.16/1.14 | 1.32/1.32 | 1305.8 |
| cxg__modified-strt-seq__homo-sapiens__interventricular- | modified STRT-seq human | umi human | 27046.0 | ok, scored | 1.32/1.32 | 1.26/1.26 | 1210.7 |
| cxg__patch-seq__mus-musculus__primary-motor-cortex__8f9 | Patch-seq mouse | None mouse | 27590.0 | ok, scored | 1.22/1.22 | 1.14/1.14 | 1048.5 |
| cxg__sci-rna-seq3__danio-rerio__endocrine-pancreas__e23 | sci-RNA-seq3 Danio rerio | umi mouse | 0.0 | completed, unscored (3 contexts; no reference edges) | - | - | 26.9 |
| cxg__sci-rna-seq3__macaca-mulatta__brodmann-1909-area-1 | sci-RNA-seq3 Macaca mulatta | umi None | 0.0 | completed, unscored (2 contexts; no reference edges) | - | - | 26.4 |
| cxg__sci-rna-seq3__mus-musculus__embryo__0bcda669 | sci-RNA-seq3 mouse | umi mouse | 28262.0 | ok, scored | 1.24/1.24 | 1.31/1.31 | 1169.0 |
| cxg__smart-seq2__homo-sapiens__islet-of-langerhans__66d | Smart-seq2 human | full_length human | 26582.0 | ok, scored | 1.30/1.30 | 1.31/1.15 | 2196.7 |
| cxg__smart-seq2__mus-musculus__aorta__0380ddce | Smart-seq2 mouse | full_length mouse | 27064.0 | ok, scored | 1.55/1.55 | 1.27/1.37 | 979.2 |
| paul15 | umi mouse | umi mouse | 3778.0 | completed, no model score | - | - | 64.6 |
| pbmc3k_filtered | umi human | umi human | 26471.0 | ok, scored | 1.31/1.49 | 1.36/1.59 | 907.4 |
| pbmc3k_raw | umi human | umi human | 26471.0 | completed, no model score | - | - | 13.8 |
