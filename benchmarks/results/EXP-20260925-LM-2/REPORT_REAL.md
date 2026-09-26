# Engine mass test on real datasets - 20 datasets

- ran end to end: 20
- refused as expected (not raw counts): 0
- errors: 0
- doublet calling fell back to the heuristic (scrublet failed) in 4/20 datasets (baseline or auto arm)

## Detection vs curated metadata

| check | n | correct |
|---|---|---|
| organism | 20 | 100% |
| technology | 20 | 90% |

## Reference recovery (CollecTRI evaluation half, AUPRC ratio; 1 = random)

| model | n | auto > baseline | auto = baseline | mean auto | mean baseline | mean rules |
|---|---|---|---|---|---|---|
| pearson_signed | 20 | 30% | 55% | 1.35 | 1.34 | 1.22 |
| spearman_signed | 20 | 30% | 50% | 1.40 | 1.39 | 1.26 |
| clr | 20 | 35% | 55% | 1.34 | 1.32 | 1.29 |

## Models (CollecTRI evaluation half; unit = dataset)

### AUPRC ratio (1 = random)

| model | n | mean auto | median auto | mean baseline | auto > baseline | s per dataset (auto) |
|---|---|---|---|---|---|---|
| pearson_motif | 20 | 1.529 | 1.539 | 1.525 | 6/20 | 3 |
| motif_prior | 20 | 1.373 | 1.386 | 1.373 | 5/20 | 3 |
| pearson_signed | 20 | 1.346 | 1.307 | 1.341 | 5/20 | 2 |

#### Declared comparisons, AUPRC ratio, auto pipeline (paired t, Holm over 4 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 20 | +0.182 | [+0.138, +0.227] | 1.36e-07 | 1 | pearson_motif |
| pearson_motif | motif_prior | 20 | +0.155 | [+0.077, +0.234] | 0.000826 | 1 | pearson_motif |

#### Declared comparisons, AUPRC ratio, baseline pipeline (paired t, Holm over 4 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 20 | +0.184 | [+0.143, +0.225] | 2.91e-08 | 1 | pearson_motif |
| pearson_motif | motif_prior | 20 | +0.152 | [+0.076, +0.228] | 0.000733 | 1 | pearson_motif |

### EPR (1 = random)

| model | n | mean auto | median auto | mean baseline | auto > baseline | s per dataset (auto) |
|---|---|---|---|---|---|---|
| pearson_motif | 20 | 2.448 | 2.482 | 2.393 | 5/20 | 3 |
| pearson_signed | 20 | 2.328 | 2.209 | 2.241 | 7/20 | 2 |
| motif_prior | 20 | 1.608 | 1.581 | 1.571 | 5/20 | 3 |

#### Declared comparisons, EPR, auto pipeline (paired t, Holm over 4 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 20 | +0.120 | [-0.046, +0.286] | 0.219 | 1 | unresolved |
| pearson_motif | motif_prior | 20 | +0.839 | [+0.495, +1.184] | 0.000127 | 1 | pearson_motif |

#### Declared comparisons, EPR, baseline pipeline (paired t, Holm over 4 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 20 | +0.152 | [-0.014, +0.318] | 0.105 | 1 | unresolved |
| pearson_motif | motif_prior | 20 | +0.823 | [+0.491, +1.154] | 0.000105 | 1 | pearson_motif |

### Per dataset, AUPRC ratio (auto pipeline)

| dataset | pearson_signed | pearson_motif | motif_prior |
|---|---|---|---|
| 10x-3-v2__homo-sapiens__entorhinal-cortex__06b91002 | 1.31 | 1.57 | 1.51 |
| 10x-3-v2__homo-sapiens__peripheral-zone-of-prostate__df287f8d | 1.18 | 1.31 | 1.23 |
| 10x-3-v2__homo-sapiens__superior-frontal-gyrus__b94e3bdf | 1.38 | 1.64 | 1.50 |
| 10x-3-v2__mus-musculus__spinal-cord__8ef59750 | 1.26 | 1.53 | 1.45 |
| 10x-3-v3__homo-sapiens__primary-motor-cortex__9b686bb6 | 1.01 | 1.20 | 1.32 |
| 10x-3-v3__homo-sapiens__prostate-gland__20de8570 | 1.60 | 1.78 | 1.40 |
| 10x-3-v3__homo-sapiens__respiratory-airway__eec804b9 | 1.75 | 1.89 | 1.50 |
| 10x-3-v3__mus-musculus__autopod-skin__389b1fd4 | 1.35 | 1.53 | 1.35 |
| 10x-3-v3__mus-musculus__epididymal-fat-pad__e2529f66 | 1.15 | 1.35 | 1.38 |
| 10x-5-v1__homo-sapiens__lung__62315937 | 1.22 | 1.55 | 1.53 |
| 10x-5-v1__homo-sapiens__lung__7b3368a5 | 1.44 | 1.83 | 1.63 |
| 10x-multiome__homo-sapiens__brain__bdc7e498 | 1.30 | 1.46 | 1.39 |
| 10x-multiome__homo-sapiens__cranial-neural-crest__805832dc | 1.08 | 1.13 | 1.12 |
| cel-seq2__homo-sapiens__bronchial-mucosa__d2bc2703 | 1.36 | 1.62 | 1.48 |
| cel-seq2__homo-sapiens__pancreas__b07e5164 | 1.25 | 1.37 | 1.26 |
| patch-seq__mus-musculus__primary-motor-cortex__341a0702 | 1.23 | 1.37 | 1.28 |
| patch-seq__mus-musculus__primary-motor-cortex__bf6a5c78 | 1.31 | 1.57 | 1.41 |
| smart-seq2__homo-sapiens__cerebral-cortex__774c18c5 | 1.39 | 1.43 | 1.21 |
| smart-seq2__mus-musculus__diaphragm__cbd62079 | 1.51 | 1.57 | 1.25 |
| smart-seq2__mus-musculus__kidney__524179b0 | 1.85 | 1.89 | 1.27 |

## Per dataset

| dataset | truth | detected | ref edges | outcome | Pearson auto/base | CLR auto/base | s |
|---|---|---|---|---|---|---|---|
| cxg__10x-3-v2__homo-sapiens__entorhinal-cortex__06b9100 | 10x 3' v2 human | umi human | 27002 | ok, scored | 1.31/1.31 | 1.23/1.23 | 314.1 |
| cxg__10x-3-v2__homo-sapiens__peripheral-zone-of-prostat | 10x 3' v2 human | umi human | 24465 | ok, scored | 1.18/1.14 | 1.33/1.33 | 213.9 |
| cxg__10x-3-v2__homo-sapiens__superior-frontal-gyrus__b9 | 10x 3' v2 human | umi human | 27002 | ok, scored | 1.38/1.38 | 1.27/1.28 | 270.1 |
| cxg__10x-3-v2__mus-musculus__spinal-cord__8ef59750 | 10x 3' v2 mouse | umi mouse | 24615 | ok, scored | 1.26/1.26 | 1.14/1.14 | 179.7 |
| cxg__10x-3-v3__homo-sapiens__primary-motor-cortex__9b68 | 10x 3' v3 human | umi human | 26968 | ok, scored | 1.01/1.02 | 1.01/1.00 | 185.3 |
| cxg__10x-3-v3__homo-sapiens__prostate-gland__20de8570 | 10x 3' v3 human | umi human | 27063 | ok, scored | 1.60/1.59 | 1.40/1.41 | 137.3 |
| cxg__10x-3-v3__homo-sapiens__respiratory-airway__eec804 | 10x 3' v3 human | umi human | 27058 | ok, scored | 1.80/1.75 | 1.61/1.50 | 175.5 |
| cxg__10x-3-v3__mus-musculus__autopod-skin__389b1fd4 | 10x 3' v3 mouse | umi mouse | 28284 | ok, scored | 1.35/1.35 | 1.34/1.34 | 248.3 |
| cxg__10x-3-v3__mus-musculus__epididymal-fat-pad__e2529f | 10x 3' v3 mouse | umi mouse | 26448 | ok, scored | 1.15/1.15 | 1.35/1.35 | 456.1 |
| cxg__10x-5-v1__homo-sapiens__lung__62315937 | 10x 5' v1 human | umi human | 20834 | ok, scored | 1.22/1.24 | 1.30/1.25 | 683.5 |
| cxg__10x-5-v1__homo-sapiens__lung__7b3368a5 | 10x 5' v1 human | umi human | 19445 | ok, scored | 1.44/1.40 | 1.59/1.59 | 603.7 |
| cxg__10x-multiome__homo-sapiens__brain__bdc7e498 | 10x multiome human | umi human | 26893 | ok, scored | 1.30/1.30 | 1.30/1.30 | 470.4 |
| cxg__10x-multiome__homo-sapiens__cranial-neural-crest__ | 10x multiome human | umi human | 26866 | ok, scored | 1.08/1.10 | 1.29/1.28 | 471.5 |
| cxg__cel-seq2__homo-sapiens__bronchial-mucosa__d2bc2703 | CEL-seq2 human | umi human | 21042 | ok, scored | 1.36/1.36 | 1.57/1.57 | 319.1 |
| cxg__cel-seq2__homo-sapiens__pancreas__b07e5164 | CEL-seq2 human | umi human | 21909 | ok, scored | 1.25/1.25 | 1.22/1.20 | 414.6 |
| cxg__patch-seq__mus-musculus__primary-motor-cortex__341 | Patch-seq mouse | None mouse | 27590 | ok, scored | 1.23/1.23 | 1.14/1.14 | 211.4 |
| cxg__patch-seq__mus-musculus__primary-motor-cortex__bf6 | Patch-seq mouse | None mouse | 27590 | ok, scored | 1.31/1.31 | 1.12/1.12 | 281.6 |
| cxg__smart-seq2__homo-sapiens__cerebral-cortex__774c18c | Smart-seq2 human | full_length human | 27029 | ok, scored | 1.39/1.34 | 1.35/1.31 | 505.3 |
| cxg__smart-seq2__mus-musculus__diaphragm__cbd62079 | Smart-seq2 mouse | full_length mouse | 27064 | ok, scored | 1.51/1.51 | 1.60/1.54 | 335.3 |
| cxg__smart-seq2__mus-musculus__kidney__524179b0 | Smart-seq2 mouse | full_length mouse | 27064 | ok, scored | 1.85/1.85 | 1.60/1.60 | 229.8 |
