# Engine mass test on real datasets - 20 datasets

- ran end to end: 20
- refused as expected (not raw counts): 0
- errors: 0
- doublet calling fell back to the heuristic (scrublet failed) in 0/20 datasets (baseline or auto arm)

## Detection vs curated metadata

| check | n | correct |
|---|---|---|
| organism | 20 | 100% |
| technology | 20 | 90% |

## Reference recovery (CollecTRI evaluation half, AUPRC ratio; 1 = random)

| model | n | auto > baseline | auto = baseline | mean auto | mean baseline | mean rules |
|---|---|---|---|---|---|---|
| pearson_signed | 20 | 20% | 60% | 1.35 | 1.34 | 1.22 |
| spearman_signed | 20 | 25% | 55% | 1.39 | 1.39 | 1.25 |
| clr | 20 | 40% | 55% | 1.34 | 1.32 | 1.29 |

## Models (CollecTRI evaluation half; unit = dataset)

### AUPRC ratio (1 = random)

| model | n | mean auto | median auto | mean baseline | auto > baseline | s per dataset (auto) |
|---|---|---|---|---|---|---|
| regdiffusion_motif | 19 | 1.600 | 1.566 | 1.599 | 3/19 | 3 |
| regdiffusion | 19 | 1.416 | 1.330 | 1.412 | 6/19 | 432 |
| grnboost2_motif | 20 | 1.232 | 1.223 | 1.237 | 2/20 | 3 |
| grnboost2 | 20 | 1.102 | 1.087 | 1.108 | 3/20 | 64 |

#### Declared comparisons, AUPRC ratio, auto pipeline (paired t, Holm over 4 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| grnboost2_motif | grnboost2 | 20 | +0.130 | [+0.108, +0.153] | 3.63e-10 | 1 | grnboost2_motif |
| regdiffusion_motif | regdiffusion | 19 | +0.184 | [+0.132, +0.237] | 1.18e-06 | 1 | regdiffusion_motif |

#### Declared comparisons, AUPRC ratio, baseline pipeline (paired t, Holm over 4 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| grnboost2_motif | grnboost2 | 20 | +0.129 | [+0.110, +0.149] | 4.8e-11 | 1 | grnboost2_motif |
| regdiffusion_motif | regdiffusion | 19 | +0.187 | [+0.138, +0.236] | 3.7e-07 | 1 | regdiffusion_motif |

### EPR (1 = random)

| model | n | mean auto | median auto | mean baseline | auto > baseline | s per dataset (auto) |
|---|---|---|---|---|---|---|
| regdiffusion | 19 | 2.797 | 2.577 | 2.725 | 5/19 | 432 |
| regdiffusion_motif | 19 | 2.794 | 2.587 | 2.732 | 5/19 | 3 |
| grnboost2_motif | 20 | 1.911 | 1.869 | 1.885 | 6/20 | 3 |
| grnboost2 | 20 | 1.306 | 1.214 | 1.338 | 3/20 | 64 |

#### Declared comparisons, EPR, auto pipeline (paired t, Holm over 4 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| grnboost2_motif | grnboost2 | 20 | +0.606 | [+0.473, +0.739] | 2.17e-08 | 1 | grnboost2_motif |
| regdiffusion_motif | regdiffusion | 19 | -0.003 | [-0.197, +0.191] | 1 | 1 | unresolved |

#### Declared comparisons, EPR, baseline pipeline (paired t, Holm over 4 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| grnboost2_motif | grnboost2 | 20 | +0.548 | [+0.422, +0.673] | 4.41e-08 | 1 | grnboost2_motif |
| regdiffusion_motif | regdiffusion | 19 | +0.006 | [-0.201, +0.214] | 1 | 1 | unresolved |

### Per dataset, AUPRC ratio (auto pipeline)

| dataset | regdiffusion | regdiffusion_motif | grnboost2 | grnboost2_motif |
|---|---|---|---|---|
| 10x-3-v2__homo-sapiens__entorhinal-cortex__06b91002 | - | - | 1.05 | 1.20 |
| 10x-3-v2__homo-sapiens__peripheral-zone-of-prostate__df287f8d | 1.22 | 1.34 | 1.07 | 1.17 |
| 10x-3-v2__homo-sapiens__superior-frontal-gyrus__b94e3bdf | 1.54 | 1.81 | 1.05 | 1.16 |
| 10x-3-v2__mus-musculus__spinal-cord__8ef59750 | 1.17 | 1.52 | 1.07 | 1.20 |
| 10x-3-v3__homo-sapiens__primary-motor-cortex__9b686bb6 | 1.12 | 1.31 | 1.03 | 1.13 |
| 10x-3-v3__homo-sapiens__prostate-gland__20de8570 | 1.46 | 1.67 | 1.14 | 1.28 |
| 10x-3-v3__homo-sapiens__respiratory-airway__eec804b9 | 1.77 | 1.92 | 1.16 | 1.33 |
| 10x-3-v3__mus-musculus__autopod-skin__389b1fd4 | 1.50 | 1.67 | 1.12 | 1.29 |
| 10x-3-v3__mus-musculus__epididymal-fat-pad__e2529f66 | 1.23 | 1.45 | 1.06 | 1.16 |
| 10x-5-v1__homo-sapiens__lung__62315937 | 1.28 | 1.61 | 1.09 | 1.31 |
| 10x-5-v1__homo-sapiens__lung__7b3368a5 | 1.55 | 1.95 | 1.13 | 1.33 |
| 10x-multiome__homo-sapiens__brain__bdc7e498 | 1.36 | 1.57 | 1.13 | 1.23 |
| 10x-multiome__homo-sapiens__cranial-neural-crest__805832dc | 1.33 | 1.33 | 1.07 | 1.10 |
| cel-seq2__homo-sapiens__bronchial-mucosa__d2bc2703 | 1.79 | 2.02 | 1.22 | 1.41 |
| cel-seq2__homo-sapiens__pancreas__b07e5164 | 1.26 | 1.38 | 1.07 | 1.17 |
| patch-seq__mus-musculus__primary-motor-cortex__341a0702 | 1.31 | 1.45 | 1.13 | 1.21 |
| patch-seq__mus-musculus__primary-motor-cortex__bf6a5c78 | 1.25 | 1.49 | 1.08 | 1.27 |
| smart-seq2__homo-sapiens__cerebral-cortex__774c18c5 | 1.33 | 1.38 | 1.07 | 1.14 |
| smart-seq2__mus-musculus__diaphragm__cbd62079 | 1.62 | 1.62 | 1.15 | 1.28 |
| smart-seq2__mus-musculus__kidney__524179b0 | 1.80 | 1.89 | 1.13 | 1.27 |

### Model errors

- cxg__10x-3-v2__homo-sapiens__entorhinal-cortex__06b91002 / regdiffusion_baseline: ValueError: 1 genes have 0 variance. Please remove these genes from your data.
- cxg__10x-3-v2__homo-sapiens__entorhinal-cortex__06b91002 / regdiffusion_motif_baseline: ValueError: 1 genes have 0 variance. Please remove these genes from your data.
- cxg__10x-3-v2__homo-sapiens__entorhinal-cortex__06b91002 / regdiffusion_auto: ValueError: 1 genes have 0 variance. Please remove these genes from your data.
- cxg__10x-3-v2__homo-sapiens__entorhinal-cortex__06b91002 / regdiffusion_motif_auto: ValueError: 1 genes have 0 variance. Please remove these genes from your data.

## Per dataset

| dataset | truth | detected | ref edges | outcome | Pearson auto/base | CLR auto/base | s |
|---|---|---|---|---|---|---|---|
| cxg__10x-3-v2__homo-sapiens__entorhinal-cortex__06b9100 | 10x 3' v2 human | umi human | 27002 | ok, scored | 1.31/1.31 | 1.23/1.23 | 518.8 |
| cxg__10x-3-v2__homo-sapiens__peripheral-zone-of-prostat | 10x 3' v2 human | umi human | 24465 | ok, scored | 1.18/1.15 | 1.33/1.33 | 1051.1 |
| cxg__10x-3-v2__homo-sapiens__superior-frontal-gyrus__b9 | 10x 3' v2 human | umi human | 27002 | ok, scored | 1.39/1.39 | 1.26/1.28 | 1174.7 |
| cxg__10x-3-v2__mus-musculus__spinal-cord__8ef59750 | 10x 3' v2 mouse | umi mouse | 24615 | ok, scored | 1.26/1.26 | 1.14/1.14 | 1072.5 |
| cxg__10x-3-v3__homo-sapiens__primary-motor-cortex__9b68 | 10x 3' v3 human | umi human | 26968 | ok, scored | 1.01/1.02 | 1.00/1.00 | 1779.5 |
| cxg__10x-3-v3__homo-sapiens__prostate-gland__20de8570 | 10x 3' v3 human | umi human | 27063 | ok, scored | 1.58/1.59 | 1.42/1.41 | 1000.4 |
| cxg__10x-3-v3__homo-sapiens__respiratory-airway__eec804 | 10x 3' v3 human | umi human | 27058 | ok, scored | 1.80/1.75 | 1.61/1.50 | 1049.3 |
| cxg__10x-3-v3__mus-musculus__autopod-skin__389b1fd4 | 10x 3' v3 mouse | umi mouse | 28284 | ok, scored | 1.35/1.35 | 1.34/1.34 | 1109.5 |
| cxg__10x-3-v3__mus-musculus__epididymal-fat-pad__e2529f | 10x 3' v3 mouse | umi mouse | 26448 | ok, scored | 1.15/1.15 | 1.35/1.35 | 1326.3 |
| cxg__10x-5-v1__homo-sapiens__lung__62315937 | 10x 5' v1 human | umi human | 20834 | ok, scored | 1.21/1.24 | 1.29/1.25 | 2456.0 |
| cxg__10x-5-v1__homo-sapiens__lung__7b3368a5 | 10x 5' v1 human | umi human | 19445 | ok, scored | 1.40/1.40 | 1.59/1.59 | 1348.8 |
| cxg__10x-multiome__homo-sapiens__brain__bdc7e498 | 10x multiome human | umi human | 26893 | ok, scored | 1.30/1.30 | 1.30/1.30 | 1228.9 |
| cxg__10x-multiome__homo-sapiens__cranial-neural-crest__ | 10x multiome human | umi human | 26866 | ok, scored | 1.07/1.10 | 1.28/1.28 | 1251.4 |
| cxg__cel-seq2__homo-sapiens__bronchial-mucosa__d2bc2703 | CEL-seq2 human | umi human | 21042 | ok, scored | 1.36/1.36 | 1.57/1.57 | 1176.1 |
| cxg__cel-seq2__homo-sapiens__pancreas__b07e5164 | CEL-seq2 human | umi human | 21909 | ok, scored | 1.25/1.24 | 1.22/1.20 | 1180.2 |
| cxg__patch-seq__mus-musculus__primary-motor-cortex__341 | Patch-seq mouse | None mouse | 27590 | ok, scored | 1.23/1.23 | 1.14/1.14 | 1020.3 |
| cxg__patch-seq__mus-musculus__primary-motor-cortex__bf6 | Patch-seq mouse | None mouse | 27590 | ok, scored | 1.31/1.31 | 1.12/1.12 | 1023.0 |
| cxg__smart-seq2__homo-sapiens__cerebral-cortex__774c18c | Smart-seq2 human | full_length human | 27029 | ok, scored | 1.39/1.34 | 1.35/1.31 | 1806.4 |
| cxg__smart-seq2__mus-musculus__diaphragm__cbd62079 | Smart-seq2 mouse | full_length mouse | 27064 | ok, scored | 1.51/1.51 | 1.60/1.54 | 1010.7 |
| cxg__smart-seq2__mus-musculus__kidney__524179b0 | Smart-seq2 mouse | full_length mouse | 27064 | ok, scored | 1.85/1.85 | 1.60/1.60 | 974.4 |
