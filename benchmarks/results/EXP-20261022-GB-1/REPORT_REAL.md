# PARTIAL - Engine mass test on real datasets - 120 planned datasets

- completed without a top-level error: 119
- with a requested-model reference score: 119
- completed without a requested-model score: 0
- refused as expected (not raw counts): 0
- permanent errors: 0
- missing or retryable: 1
- doublet calling fell back to the heuristic (scrublet failed) in 0/120 datasets (baseline or auto arm)

## Detection vs curated metadata

| check | n | correct |
|---|---|---|
| organism | 120 | 100% |
| technology | 120 | 100% |

## Reference recovery (CollecTRI evaluation half, AUPRC ratio; 1 = random)

| model | n | auto > baseline | auto = baseline | mean auto | mean baseline | mean rules |
|---|---|---|---|---|---|---|
| pearson_signed | 119 | 45% | 50% | 1.39 | 1.31 | 1.26 |
| spearman_signed | 119 | 44% | 48% | 1.47 | 1.39 | 1.30 |
| clr | 119 | 29% | 58% | 1.38 | 1.36 | 1.25 |

## Models (CollecTRI evaluation half; unit = dataset)

### AUPRC ratio (1 = random)

| model | n | mean auto | median auto | mean baseline | auto > baseline | s per dataset (auto) |
|---|---|---|---|---|---|---|
| pearson_motif | 119 | 1.579 | 1.550 | 1.514 | 45/119 | 3 |
| pearson_signed | 119 | 1.373 | 1.349 | 1.313 | 44/119 | 2 |

#### Declared comparisons, AUPRC ratio, auto pipeline (paired t, Holm over 2 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 119 | +0.206 | [+0.183, +0.229] | 1.58e-35 | 1 | pearson_motif |

#### Declared comparisons, AUPRC ratio, baseline pipeline (paired t, Holm over 2 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 119 | +0.202 | [+0.181, +0.223] | 9.21e-38 | 1 | pearson_motif |

### EPR (1 = random)

| model | n | mean auto | median auto | mean baseline | auto > baseline | s per dataset (auto) |
|---|---|---|---|---|---|---|
| pearson_motif | 119 | 2.490 | 2.431 | 2.246 | 50/119 | 3 |
| pearson_signed | 119 | 2.341 | 2.297 | 2.047 | 47/119 | 2 |

#### Declared comparisons, EPR, auto pipeline (paired t, Holm over 2 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 119 | +0.149 | [+0.067, +0.230] | 0.000424 | 1 | pearson_motif |

#### Declared comparisons, EPR, baseline pipeline (paired t, Holm over 2 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 119 | +0.199 | [+0.131, +0.267] | 5.29e-08 | 1 | pearson_motif |

### Per dataset, AUPRC ratio (auto pipeline)

| dataset | pearson_signed | pearson_motif |
|---|---|---|
| 10x-3-v2__homo-sapiens__anterior-wall-of-left-ventricl__4ed927e9 | 1.24 | 1.31 |
| 10x-3-v2__homo-sapiens__blood__ae5341b8 | 1.54 | 1.72 |
| 10x-3-v2__homo-sapiens__brain__dc30c3ec | 1.11 | 1.29 |
| 10x-3-v2__homo-sapiens__breast__9a64bf99 | 1.49 | 1.58 |
| 10x-3-v2__homo-sapiens__caudate-lobe-of-liver__a43aa46b | 1.27 | 1.49 |
| 10x-3-v2__homo-sapiens__cortex-of-kidney__105c7dad | 1.78 | 1.95 |
| 10x-3-v2__homo-sapiens__dental-pulp__bdd109f1 | 1.98 | 2.24 |
| 10x-3-v2__homo-sapiens__epithelium-of-rectum__80d7111d | 1.22 | 1.45 |
| 10x-3-v2__homo-sapiens__fovea-centralis__babbf9f3 | 1.15 | 1.42 |
| 10x-3-v2__homo-sapiens__kidney__8dc29642 | 1.55 | 1.61 |
| 10x-3-v2__homo-sapiens__lung__6e00ccf7 | 1.23 | 1.51 |
| 10x-3-v2__homo-sapiens__mammary-gland-connective-tissu__68b6114f | 1.85 | 2.18 |
| 10x-3-v2__homo-sapiens__ovary__1f1c5c14 | 1.45 | 1.55 |
| 10x-3-v2__homo-sapiens__respiratory-airway__8fcf0ccc | 1.00 | 1.20 |
| 10x-3-v2__homo-sapiens__skin-of-prepuce-of-penis__f512b8b6 | 1.40 | 1.49 |
| 10x-3-v2__homo-sapiens__thymus__11ff73e8 | 1.50 | 1.65 |
| 10x-3-v2__homo-sapiens__zone-of-skin__124744b8 | 1.41 | 1.53 |
| 10x-3-v2__mus-musculus__hindlimb-muscle__bedee313 | 1.61 | 1.62 |
| 10x-3-v3__homo-sapiens__abdomen__7b95971d | 1.15 | 1.22 |
| 10x-3-v3__homo-sapiens__abdomen__97d9238c | 1.37 | 1.50 |
| 10x-3-v3__homo-sapiens__apex-of-heart__2e9d2f32 | 1.39 | 1.47 |
| 10x-3-v3__homo-sapiens__blood__19e46756 | 1.65 | 1.82 |
| 10x-3-v3__homo-sapiens__blood__20fae564 | 1.35 | 1.50 |
| 10x-3-v3__homo-sapiens__blood__738942eb | 1.42 | 1.85 |
| 10x-3-v3__homo-sapiens__blood__a199ca73 | 1.69 | 1.91 |
| 10x-3-v3__homo-sapiens__blood__d9f3a82c | 1.22 | 1.47 |
| 10x-3-v3__homo-sapiens__blood__de2c780c | 1.65 | 1.73 |
| 10x-3-v3__homo-sapiens__bone-marrow__343ff97c | 1.56 | 1.74 |
| 10x-3-v3__homo-sapiens__breast__24ec2dc5 | 1.39 | 1.58 |
| 10x-3-v3__homo-sapiens__brodmann-1909-area-13__f16f4108 | 1.02 | 1.30 |
| 10x-3-v3__homo-sapiens__brodmann-1909-area-4__c8f83821 | 1.03 | 1.28 |
| 10x-3-v3__homo-sapiens__brown-preadipocyte__137ca42b | 1.17 | 1.33 |
| 10x-3-v3__homo-sapiens__calcaneal-tendon__ef0fc6ac | 1.44 | 1.81 |
| 10x-3-v3__homo-sapiens__calcaneal-tendon__f641d68a | 1.24 | 1.63 |
| 10x-3-v3__homo-sapiens__caudate-nucleus__c893ddc3 | 1.35 | 1.59 |
| 10x-3-v3__homo-sapiens__cerebellum__be884a28 | 1.26 | 1.50 |
| 10x-3-v3__homo-sapiens__cornea__8058d760 | 1.19 | 1.44 |
| 10x-3-v3__homo-sapiens__cortex-of-kidney__07854d9c | 1.30 | 1.46 |
| 10x-3-v3__homo-sapiens__cortex-of-kidney__0bae7ebf | 1.38 | 1.57 |
| 10x-3-v3__homo-sapiens__cranial-nerve-ii__374edd6c | 1.33 | 1.48 |
| 10x-3-v3__homo-sapiens__cranial-nerve-ii__9f1a899f | 1.20 | 1.43 |
| 10x-3-v3__homo-sapiens__dome-of-urinary-bladder__601c7454 | 1.22 | 1.41 |
| 10x-3-v3__homo-sapiens__endometrium__76150f40 | 1.46 | 1.70 |
| 10x-3-v3__homo-sapiens__endometrium__7a419302 | 1.48 | 1.63 |
| 10x-3-v3__homo-sapiens__entorhinal-cortex__24584be9 | 1.08 | 1.36 |
| 10x-3-v3__homo-sapiens__fallopian-tube__30f5e171 | 2.09 | 1.97 |
| 10x-3-v3__homo-sapiens__fovea-centralis__11ef37ee | 1.15 | 1.51 |
| 10x-3-v3__homo-sapiens__fovea-centralis__af8b241a | 1.96 | 2.01 |
| 10x-3-v3__homo-sapiens__fovea-centralis__e404d417 | 1.08 | 1.51 |
| 10x-3-v3__homo-sapiens__fovea-centralis__f5be9ed2 | 1.50 | 1.61 |
| 10x-3-v3__homo-sapiens__heart-left-ventricle__f15e263b | 1.25 | 1.37 |
| 10x-3-v3__homo-sapiens__kidney__ea01c125 | 1.70 | 1.74 |
| 10x-3-v3__homo-sapiens__left-atrium-auricular-region__87a06a2a | 1.67 | 1.73 |
| 10x-3-v3__homo-sapiens__lens-of-camera-type-eye__489318a0 | 1.02 | 1.45 |
| 10x-3-v3__homo-sapiens__liver__0bbf93aa | 1.44 | 1.80 |
| 10x-3-v3__homo-sapiens__liver__ec6ceff8 | 1.29 | 1.32 |
| 10x-3-v3__homo-sapiens__lung-parenchyma__01ff5cf0 | 1.35 | 1.45 |
| 10x-3-v3__homo-sapiens__lung-parenchyma__8fbed309 | 1.50 | 1.68 |
| 10x-3-v3__homo-sapiens__lung__6ba45a2b | 1.43 | 2.21 |
| 10x-3-v3__homo-sapiens__lung__d8da613f | 1.21 | 1.42 |
| 10x-3-v3__homo-sapiens__mesenteric-artery__42b6a476 | 1.61 | 1.74 |
| 10x-3-v3__homo-sapiens__minor-salivary-gland__9aef5aa0 | 1.63 | 1.81 |
| 10x-3-v3__homo-sapiens__neocortex__3a8aec06 | 1.58 | 1.90 |
| 10x-3-v3__homo-sapiens__neocortex__4859642c | 1.25 | 1.39 |
| 10x-3-v3__homo-sapiens__occipital-cortex__e6361237 | 1.31 | 1.55 |
| 10x-3-v3__homo-sapiens__oropharynx__b6b5ea88 | 1.83 | 1.91 |
| 10x-3-v3__homo-sapiens__ovary__3a702020 | 1.19 | 1.33 |
| 10x-3-v3__homo-sapiens__ovary__a9bedd04 | 1.23 | 1.35 |
| 10x-3-v3__homo-sapiens__placenta__b7efc0ba | 1.44 | 1.69 |
| 10x-3-v3__homo-sapiens__prefrontal-cortex__fe2eecbc | 1.27 | 1.54 |
| 10x-3-v3__homo-sapiens__right-atrium-auricular-region__8f4f8502 | 1.41 | 1.61 |
| 10x-3-v3__homo-sapiens__right-temporal-lobe__145dcf6a | 1.59 | 2.01 |
| 10x-3-v3__homo-sapiens__skin-of-body__b796f27e | 1.56 | 1.71 |
| 10x-3-v3__homo-sapiens__tendon-of-semitendinosus__06ef6b36 | 1.13 | 1.45 |
| 10x-3-v3__homo-sapiens__upper-lobe-of-left-lung__2ac76f1b | 1.81 | 1.86 |
| 10x-3-v3__homo-sapiens__visceral-fat__1405dd3e | 1.39 | 1.74 |
| 10x-3-v3__mus-musculus__aorta__f7417545 | 1.36 | 1.42 |
| 10x-3-v3__mus-musculus__brain-gray-matter__6347cc90 | 1.23 | 1.40 |
| 10x-3-v3__mus-musculus__brain__e056b20b | 1.69 | 1.72 |
| 10x-3-v3__mus-musculus__cerebral-cortex__0e1c7067 | 1.21 | 1.55 |
| 10x-3-v3__mus-musculus__cerebral-cortex__12990215 | 1.00 | 1.10 |
| 10x-3-v3__mus-musculus__gastrocnemius__2ab892f9 | 1.06 | 1.18 |
| 10x-3-v3__mus-musculus__gingiva__e04cb7d9 | 1.47 | 1.54 |
| 10x-3-v3__mus-musculus__liver__e9de826e | 1.25 | 1.34 |
| 10x-3-v3__mus-musculus__middle-lobe-of-right-lung__89840cec | 1.40 | 1.67 |
| 10x-3-v3__mus-musculus__ovary__55d66f2a | 1.37 | 1.57 |
| 10x-3-v3__mus-musculus__tracheal-epithelial-cell__732a691e | 1.32 | 1.42 |
| 10x-5-transcription-profiling__homo-sapiens__adrenal-tissue__5af90777 | 1.35 | 1.66 |
| 10x-5-transcription-profiling__homo-sapiens__lung__8dc858b0 | 1.60 | 1.87 |
| 10x-5-v1__homo-sapiens__blood__5b25a46f | 1.42 | 1.81 |
| 10x-5-v1__homo-sapiens__blood__c2a461b1 | 1.51 | 1.93 |
| 10x-5-v1__homo-sapiens__blood__fa8605cf | 1.45 | 1.89 |
| 10x-5-v1__homo-sapiens__bone-marrow__965386e9 | - | - |
| 10x-5-v1__homo-sapiens__lymph-node__99950e99 | 1.30 | 1.54 |
| 10x-5-v1__homo-sapiens__t-cell__0c774045 | 1.31 | 1.64 |
| 10x-5-v1__homo-sapiens__t-cell__460dfe21 | 1.23 | 1.63 |
| 10x-5-v1__homo-sapiens__thymus__8d73847b | 1.15 | 1.44 |
| 10x-5-v1__mus-musculus__gonad__a13bda79 | 1.40 | 1.48 |
| 10x-5-v2__homo-sapiens__blood__30c2a6fd | 1.28 | 1.61 |
| 10x-5-v2__homo-sapiens__blood__bd65a70f | 1.48 | 1.60 |
| 10x-multiome__homo-sapiens__amygdala__835fa5c4 | 1.27 | 1.44 |
| 10x-multiome__homo-sapiens__anterior-cingulate-cortex__56a4bc14 | 1.16 | 1.40 |
| 10x-multiome__homo-sapiens__frontal-lobe__01b84709 | 1.23 | 1.44 |
| 10x-multiome__homo-sapiens__kidney__867757c1 | 1.16 | 1.35 |
| 10x-multiome__homo-sapiens__retina__322c816f | 1.10 | 1.28 |
| 10x-multiome__mus-musculus__kidney__25818bf7 | 1.10 | 1.29 |
| 10x-next-gem-flex-v1__homo-sapiens__endometrium__fb90c70d | 1.47 | 1.46 |
| 10x-next-gem-flex-v1__homo-sapiens__parenchyma-of-mammary-gland__fbdd8c17 | 1.50 | 1.56 |
| asteria-scrna-seq-kit__homo-sapiens__blood__27d3952f | 1.29 | 1.68 |
| bd-rhapsody-whole-transcriptome-analysis__homo-sapiens__liver__e3ed2ba4 | 1.53 | 1.65 |
| mars-seq__homo-sapiens__bone-marrow__e22482ee | 1.33 | 1.60 |
| mct-seq__homo-sapiens__brodmann-1909-area-46__43aa19d2 | 1.11 | 1.61 |
| sci-plex__homo-sapiens__k-562__38a78e91 | 1.19 | 1.44 |
| seq-well-s3__homo-sapiens__nasopharynx__13b61a7d | 2.14 | 2.33 |
| seq-well__homo-sapiens__blood__456e8b9b | 1.32 | 1.56 |
| smart-seq-v4__mus-musculus__primary-motor-cortex__1304e107 | 1.03 | 1.22 |
| smart-seq2__homo-sapiens__blood__e04daea4 | 1.51 | 1.67 |
| smart-seq2__homo-sapiens__heart-left-ventricle__07e6773a | 1.21 | 1.37 |
| smart-seq__mus-musculus__primary-visual-cortex__3a15ab1c | 1.06 | 1.32 |
| split-seq__homo-sapiens__cerebellum__b2dda353 | 1.14 | 1.32 |

## Per dataset

| dataset | truth | detected | ref edges | outcome | Pearson auto/base | CLR auto/base | s |
|---|---|---|---|---|---|---|---|
| cxg__10x-3-v2__homo-sapiens__anterior-wall-of-left-vent | 10x 3' v2 human | umi human | 27041 | ok, scored | 1.24/1.18 | 1.40/1.22 | 14580.4 |
| cxg__10x-3-v2__homo-sapiens__blood__ae5341b8 | 10x 3' v2 human | umi human | 18372 | ok, scored | 1.54/1.54 | 1.40/1.40 | 1275.9 |
| cxg__10x-3-v2__homo-sapiens__brain__dc30c3ec | 10x 3' v2 human | umi human | 21808 | ok, scored | 1.11/1.16 | 1.10/1.08 | 1051.4 |
| cxg__10x-3-v2__homo-sapiens__breast__9a64bf99 | 10x 3' v2 human | umi human | 23024 | ok, scored | 1.37/1.37 | 1.54/1.54 | 3137.4 |
| cxg__10x-3-v2__homo-sapiens__caudate-lobe-of-liver__a43 | 10x 3' v2 human | umi human | 27041 | ok, scored | 1.27/1.27 | 1.48/1.48 | 429.8 |
| cxg__10x-3-v2__homo-sapiens__cortex-of-kidney__105c7dad | 10x 3' v2 human | umi human | 27009 | ok, scored | 1.78/1.78 | 1.52/1.52 | 2121.6 |
| cxg__10x-3-v2__homo-sapiens__dental-pulp__bdd109f1 | 10x 3' v2 human | umi human | 22123 | ok, scored | 1.98/1.98 | 1.54/1.54 | 526.8 |
| cxg__10x-3-v2__homo-sapiens__epithelium-of-rectum__80d7 | 10x 3' v2 human | umi human | 18813 | ok, scored | 1.42/1.22 | 1.44/1.44 | 427.8 |
| cxg__10x-3-v2__homo-sapiens__fovea-centralis__babbf9f3 | 10x 3' v2 human | umi human | 27000 | ok, scored | 1.20/1.20 | 1.28/1.22 | 1230.2 |
| cxg__10x-3-v2__homo-sapiens__kidney__8dc29642 | 10x 3' v2 human | umi human | 25253 | ok, scored | 1.55/1.26 | 1.39/1.32 | 1099.4 |
| cxg__10x-3-v2__homo-sapiens__lung__6e00ccf7 | 10x 3' v2 human | umi human | 24245 | ok, scored | 1.27/1.23 | 1.46/1.46 | 518.8 |
| cxg__10x-3-v2__homo-sapiens__mammary-gland-connective-t | 10x 3' v2 human | umi human | 26209 | ok, scored | 1.85/1.85 | 1.71/1.71 | 1938.7 |
| cxg__10x-3-v2__homo-sapiens__ovary__1f1c5c14 | 10x 3' v2 human | umi human | 27041 | ok, scored | 1.58/1.45 | 1.40/1.40 | 611.0 |
| cxg__10x-3-v2__homo-sapiens__respiratory-airway__8fcf0c | 10x 3' v2 human | umi human | 23392 | ok, scored | 1.20/1.00 | 1.33/1.33 | 1437.8 |
| cxg__10x-3-v2__homo-sapiens__skin-of-prepuce-of-penis__ | 10x 3' v2 human | umi human | 22061 | ok, scored | 1.40/1.06 | 1.34/1.34 | 4010.8 |
| cxg__10x-3-v2__homo-sapiens__thymus__11ff73e8 | 10x 3' v2 human | umi human | 27041 | ok, scored | 1.50/1.21 | 1.36/1.21 | 5162.2 |
| cxg__10x-3-v2__homo-sapiens__zone-of-skin__124744b8 | 10x 3' v2 human | umi human | 26785 | ok, scored | 1.41/1.41 | 1.40/1.40 | 692.8 |
| cxg__10x-3-v2__mus-musculus__hindlimb-muscle__bedee313 | 10x 3' v2 mouse | umi mouse | 21227 | ok, scored | 1.61/1.30 | 1.66/1.66 | 6492.5 |
| cxg__10x-3-v3__homo-sapiens__abdomen__7b95971d | 10x 3' v3 human | umi human | 27056 | ok, scored | 1.15/1.15 | 1.11/1.11 | 1169.3 |
| cxg__10x-3-v3__homo-sapiens__abdomen__97d9238c | 10x 3' v3 human | umi human | 26997 | ok, scored | 1.37/1.37 | 1.34/1.34 | 1011.0 |
| cxg__10x-3-v3__homo-sapiens__apex-of-heart__2e9d2f32 | 10x 3' v3 human | umi human | 27107 | ok, scored | 1.41/1.19 | 1.30/1.21 | 3696.6 |
| cxg__10x-3-v3__homo-sapiens__blood__19e46756 | 10x 3' v3 human | umi human | 26964 | ok, scored | 1.65/1.65 | 1.48/1.48 | 2175.8 |
| cxg__10x-3-v3__homo-sapiens__blood__20fae564 | 10x 3' v3 human | umi human | 26992 | ok, scored | 1.35/1.35 | 1.58/1.36 | 10533.3 |
| cxg__10x-3-v3__homo-sapiens__blood__738942eb | 10x 3' v3 human | umi human | 21171 | ok, scored | 1.42/1.24 | 1.29/1.27 | 1780.1 |
| cxg__10x-3-v3__homo-sapiens__blood__a199ca73 | 10x 3' v3 human | umi human | 27057 | ok, scored | 1.69/1.69 | 1.35/1.38 | 6272.5 |
| cxg__10x-3-v3__homo-sapiens__blood__d9f3a82c | 10x 3' v3 human | umi human | 26339 | ok, scored | 1.28/1.22 | 1.32/1.32 | 19291.4 |
| cxg__10x-3-v3__homo-sapiens__blood__de2c780c | 10x 3' v3 human | umi human | 27058 | ok, scored | 1.52/1.65 | 1.64/1.64 | 3013.3 |
| cxg__10x-3-v3__homo-sapiens__bone-marrow__343ff97c | 10x 3' v3 human | umi human | 27057 | ok, scored | 1.56/1.56 | 1.34/1.39 | 1001.1 |
| cxg__10x-3-v3__homo-sapiens__breast__24ec2dc5 | 10x 3' v3 human | umi human | 26863 | ok, scored | 1.39/1.39 | 1.30/1.30 | 656.2 |
| cxg__10x-3-v3__homo-sapiens__brodmann-1909-area-13__f16 | 10x 3' v3 human | umi human | 27020 | ok, scored | 1.00/1.00 | 1.10/1.10 | 434.8 |
| cxg__10x-3-v3__homo-sapiens__brodmann-1909-area-4__c8f8 | 10x 3' v3 human | umi human | 27058 | ok, scored | 1.09/1.09 | 1.24/1.22 | 159.1 |
| cxg__10x-3-v3__homo-sapiens__brown-preadipocyte__137ca4 | 10x 3' v3 human | umi human | 27058 | ok, scored | 1.17/1.17 | 1.30/1.30 | 543.6 |
| cxg__10x-3-v3__homo-sapiens__calcaneal-tendon__ef0fc6ac | 10x 3' v3 human | umi human | 27107 | ok, scored | 1.44/1.31 | 1.60/1.57 | 5005.7 |
| cxg__10x-3-v3__homo-sapiens__calcaneal-tendon__f641d68a | 10x 3' v3 human | umi human | 27062 | ok, scored | 1.24/1.24 | 1.34/1.34 | 5729.0 |
| cxg__10x-3-v3__homo-sapiens__caudate-nucleus__c893ddc3 | 10x 3' v3 human | umi human | 25652 | ok, scored | 1.37/1.00 | 1.39/0.97 | 7813.6 |
| cxg__10x-3-v3__homo-sapiens__cerebellum__be884a28 | 10x 3' v3 human | umi human | 24484 | ok, scored | 1.26/1.21 | 1.14/1.35 | 1630.2 |
| cxg__10x-3-v3__homo-sapiens__cornea__8058d760 | 10x 3' v3 human | umi human | 27062 | ok, scored | 1.19/1.10 | 1.45/1.37 | 9979.5 |
| cxg__10x-3-v3__homo-sapiens__cortex-of-kidney__07854d9c | 10x 3' v3 human | umi human | 25958 | ok, scored | 1.30/1.20 | 1.48/1.30 | 10545.5 |
| cxg__10x-3-v3__homo-sapiens__cortex-of-kidney__0bae7ebf | 10x 3' v3 human | umi human | 25885 | ok, scored | 1.38/1.28 | 1.37/1.37 | 1366.3 |
| cxg__10x-3-v3__homo-sapiens__cranial-nerve-ii__374edd6c | 10x 3' v3 human | umi human | 27062 | ok, scored | 1.33/1.33 | 1.38/1.38 | 1236.3 |
| cxg__10x-3-v3__homo-sapiens__cranial-nerve-ii__9f1a899f | 10x 3' v3 human | umi human | 25689 | ok, scored | 1.23/1.21 | 1.24/1.29 | 10901.6 |
| cxg__10x-3-v3__homo-sapiens__dome-of-urinary-bladder__6 | 10x 3' v3 human | umi human | 27061 | ok, scored | 1.22/1.22 | 1.52/1.52 | 237.2 |
| cxg__10x-3-v3__homo-sapiens__endometrium__76150f40 | 10x 3' v3 human | umi human | 25509 | ok, scored | 1.47/1.03 | 1.63/1.63 | 12607.7 |
| cxg__10x-3-v3__homo-sapiens__endometrium__7a419302 | 10x 3' v3 human | umi human | 27062 | ok, scored | 1.48/1.48 | 1.42/1.42 | 758.9 |
| cxg__10x-3-v3__homo-sapiens__entorhinal-cortex__24584be | 10x 3' v3 human | umi human | 25241 | ok, scored | 1.08/1.00 | 1.18/1.20 | 683.0 |
| cxg__10x-3-v3__homo-sapiens__fallopian-tube__30f5e171 | 10x 3' v3 human | umi human | 27063 | ok, scored | 2.09/2.09 | 1.79/1.79 | 213.0 |
| cxg__10x-3-v3__homo-sapiens__fovea-centralis__11ef37ee | 10x 3' v3 human | umi human | 26785 | ok, scored | 1.17/1.14 | 1.12/1.13 | 369.8 |
| cxg__10x-3-v3__homo-sapiens__fovea-centralis__af8b241a | 10x 3' v3 human | umi human | 23322 | ok, scored | 1.96/1.79 | 1.67/1.67 | 456.4 |
| cxg__10x-3-v3__homo-sapiens__fovea-centralis__e404d417 | 10x 3' v3 human | umi human | 27062 | ok, scored | 1.10/0.98 | 1.30/1.27 | 1307.9 |
| cxg__10x-3-v3__homo-sapiens__fovea-centralis__f5be9ed2 | 10x 3' v3 human | umi human | 26993 | ok, scored | 1.50/1.50 | 1.41/1.41 | 2265.7 |
| cxg__10x-3-v3__homo-sapiens__heart-left-ventricle__f15e | 10x 3' v3 human | umi human | 21232 | ok, scored | 1.25/1.17 | 1.32/1.32 | 1951.4 |
| cxg__10x-3-v3__homo-sapiens__kidney__ea01c125 | 10x 3' v3 human | umi human | 27041 | ok, scored | 1.70/1.61 | 1.44/1.31 | 632.6 |
| cxg__10x-3-v3__homo-sapiens__left-atrium-auricular-regi | 10x 3' v3 human | umi human | 23100 | ok, scored | 1.67/1.30 | 1.71/1.71 | 2630.8 |
| cxg__10x-3-v3__homo-sapiens__lens-of-camera-type-eye__4 | 10x 3' v3 human | umi human | 27102 | ok, scored | 1.02/1.02 | 1.13/1.12 | 1106.9 |
| cxg__10x-3-v3__homo-sapiens__liver__0bbf93aa | 10x 3' v3 human | umi human | 27059 | ok, scored | 1.50/1.26 | 1.39/1.49 | 2450.5 |
| cxg__10x-3-v3__homo-sapiens__liver__ec6ceff8 | 10x 3' v3 human | umi human | 27058 | ok, scored | 1.29/1.29 | 1.28/1.28 | 6030.1 |
| cxg__10x-3-v3__homo-sapiens__lung-parenchyma__01ff5cf0 | 10x 3' v3 human | umi human | 27108 | ok, scored | 1.35/1.35 | 1.54/1.54 | 6413.2 |
| cxg__10x-3-v3__homo-sapiens__lung-parenchyma__8fbed309 | 10x 3' v3 human | umi human | 27007 | ok, scored | 1.52/1.50 | 1.84/1.84 | 3798.1 |
| cxg__10x-3-v3__homo-sapiens__lung__6ba45a2b | 10x 3' v3 human | umi human | 21335 | ok, scored | 1.43/1.43 | 1.49/1.49 | 902.2 |
| cxg__10x-3-v3__homo-sapiens__lung__d8da613f | 10x 3' v3 human | umi human | 26660 | ok, scored | 1.21/1.21 | 1.41/1.41 | 5290.7 |
| cxg__10x-3-v3__homo-sapiens__mesenteric-artery__42b6a47 | 10x 3' v3 human | umi human | 27058 | ok, scored | 1.65/1.59 | 1.75/1.59 | 772.7 |
| cxg__10x-3-v3__homo-sapiens__minor-salivary-gland__9aef | 10x 3' v3 human | umi human | 26661 | ok, scored | 1.63/1.63 | 1.22/1.30 | 328.9 |
| cxg__10x-3-v3__homo-sapiens__neocortex__3a8aec06 | 10x 3' v3 human | umi human | 23808 | ok, scored | 1.58/1.58 | 1.42/1.42 | 2877.2 |
| cxg__10x-3-v3__homo-sapiens__neocortex__4859642c | 10x 3' v3 human | umi human | 19978 | ok, scored | 1.30/1.25 | 1.14/1.11 | 118.2 |
| cxg__10x-3-v3__homo-sapiens__occipital-cortex__e6361237 | 10x 3' v3 human | umi human | 25934 | ok, scored | 1.31/1.06 | 1.46/1.16 | 5740.0 |
| cxg__10x-3-v3__homo-sapiens__oropharynx__b6b5ea88 | 10x 3' v3 human | umi human | 27061 | ok, scored | 1.83/1.57 | 1.71/1.46 | 4825.5 |
| cxg__10x-3-v3__homo-sapiens__ovary__3a702020 | 10x 3' v3 human | umi human | 25009 | ok, scored | 1.16/1.19 | 1.36/1.36 | 2756.7 |
| cxg__10x-3-v3__homo-sapiens__ovary__a9bedd04 | 10x 3' v3 human | umi human | 27062 | ok, scored | 1.31/1.31 | 1.45/1.35 | 254.1 |
| cxg__10x-3-v3__homo-sapiens__placenta__b7efc0ba | 10x 3' v3 human | umi human | 27057 | ok, scored | 1.44/1.31 | 1.32/1.30 | 4709.1 |
| cxg__10x-3-v3__homo-sapiens__prefrontal-cortex__fe2eecb | 10x 3' v3 human | umi human | 27107 | ok, scored | 1.27/1.27 | 1.19/1.19 | 286.9 |
| cxg__10x-3-v3__homo-sapiens__right-atrium-auricular-reg | 10x 3' v3 human | umi human | 27062 | ok, scored | 1.41/1.13 | 1.31/1.30 | 11588.9 |
| cxg__10x-3-v3__homo-sapiens__right-temporal-lobe__145dc | 10x 3' v3 human | umi human | 27020 | ok, scored | 1.59/1.23 | 1.50/1.37 | 10226.4 |
| cxg__10x-3-v3__homo-sapiens__skin-of-body__b796f27e | 10x 3' v3 human | umi human | 24870 | ok, scored | 1.52/1.36 | 1.54/1.54 | 2126.6 |
| cxg__10x-3-v3__homo-sapiens__tendon-of-semitendinosus__ | 10x 3' v3 human | umi human | 21240 | ok, scored | 1.13/1.13 | 1.51/1.51 | 1108.5 |
| cxg__10x-3-v3__homo-sapiens__upper-lobe-of-left-lung__2 | 10x 3' v3 human | umi human | 27093 | ok, scored | 1.81/1.27 | 1.53/1.53 | 3331.4 |
| cxg__10x-3-v3__homo-sapiens__visceral-fat__1405dd3e | 10x 3' v3 human | umi human | 27062 | ok, scored | 1.39/1.39 | 1.47/1.47 | 925.2 |
| cxg__10x-3-v3__mus-musculus__aorta__f7417545 | 10x 3' v3 mouse | umi mouse | 26539 | ok, scored | 1.36/1.36 | 1.30/1.30 | 4464.8 |
| cxg__10x-3-v3__mus-musculus__brain-gray-matter__6347cc9 | 10x 3' v3 mouse | umi mouse | 21843 | ok, scored | 1.23/1.23 | 1.13/1.13 | 515.6 |
| cxg__10x-3-v3__mus-musculus__brain__e056b20b | 10x 3' v3 mouse | umi mouse | 28097 | ok, scored | 1.69/1.32 | 1.42/1.36 | 6464.1 |
| cxg__10x-3-v3__mus-musculus__cerebral-cortex__0e1c7067 | 10x 3' v3 mouse | umi mouse | 28284 | ok, scored | 1.21/1.20 | 1.15/1.15 | 1211.0 |
| cxg__10x-3-v3__mus-musculus__cerebral-cortex__12990215 | 10x 3' v3 mouse | umi mouse | 28059 | ok, scored | 1.00/1.00 | 1.18/1.18 | 9232.5 |
| cxg__10x-3-v3__mus-musculus__gastrocnemius__2ab892f9 | 10x 3' v3 mouse | umi mouse | 25525 | ok, scored | 1.06/1.06 | 1.40/1.40 | 2901.8 |
| cxg__10x-3-v3__mus-musculus__gingiva__e04cb7d9 | 10x 3' v3 mouse | umi mouse | 28191 | ok, scored | 1.47/1.47 | 1.69/1.69 | 1042.9 |
| cxg__10x-3-v3__mus-musculus__liver__e9de826e | 10x 3' v3 mouse | umi mouse | 28194 | ok, scored | 1.25/1.20 | 1.22/1.37 | 1698.8 |
| cxg__10x-3-v3__mus-musculus__middle-lobe-of-right-lung_ | 10x 3' v3 mouse | umi mouse | 27837 | ok, scored | 1.40/1.35 | 1.46/1.46 | 6720.0 |
| cxg__10x-3-v3__mus-musculus__ovary__55d66f2a | 10x 3' v3 mouse | umi mouse | 28172 | ok, scored | 1.51/1.37 | 1.78/1.78 | 3170.7 |
| cxg__10x-3-v3__mus-musculus__tracheal-epithelial-cell__ | 10x 3' v3 mouse | umi mouse | 28157 | ok, scored | 1.32/1.32 | 1.35/1.23 | 440.2 |
| cxg__10x-5-transcription-profiling__homo-sapiens__adren | 10x 5' transcription profiling human | umi human | 23509 | ok, scored | 1.35/1.24 | 1.25/1.25 | 20092.4 |
| cxg__10x-5-transcription-profiling__homo-sapiens__lung_ | 10x 5' transcription profiling human | umi human | 27030 | ok, scored | 1.60/1.60 | 1.39/1.39 | 507.3 |
| cxg__10x-5-v1__homo-sapiens__blood__5b25a46f | 10x 5' v1 human | umi human | 27062 | ok, scored | 1.42/1.42 | 1.33/1.33 | 791.5 |
| cxg__10x-5-v1__homo-sapiens__blood__c2a461b1 | 10x 5' v1 human | umi human | 22480 | ok, scored | 1.51/1.70 | 1.38/1.55 | 5544.3 |
| cxg__10x-5-v1__homo-sapiens__blood__fa8605cf | 10x 5' v1 human | umi human | 23167 | ok, scored | 1.45/1.45 | 1.51/1.51 | 2956.1 |
| cxg__10x-5-v1__homo-sapiens__bone-marrow__965386e9 | 10x 5' v1 human | umi human | 24394 | ERROR: ZeroDivisionError: division by zero | - | - | 347.3 |
| cxg__10x-5-v1__homo-sapiens__lymph-node__99950e99 | 10x 5' v1 human | umi human | 26785 | ok, scored | 1.30/1.34 | 1.27/1.27 | 6120.3 |
| cxg__10x-5-v1__homo-sapiens__t-cell__0c774045 | 10x 5' v1 human | umi human | 21671 | ok, scored | 1.31/1.22 | 1.30/1.15 | 7009.4 |
| cxg__10x-5-v1__homo-sapiens__t-cell__460dfe21 | 10x 5' v1 human | umi human | 26995 | ok, scored | 1.23/1.23 | 1.20/1.20 | 1426.5 |
| cxg__10x-5-v1__homo-sapiens__thymus__8d73847b | 10x 5' v1 human | umi human | 27058 | ok, scored | 1.30/1.15 | 1.28/1.28 | 717.9 |
| cxg__10x-5-v1__mus-musculus__gonad__a13bda79 | 10x 5' v1 mouse | umi mouse | 26562 | ok, scored | 1.40/1.40 | 1.48/1.48 | 2791.0 |
| cxg__10x-5-v2__homo-sapiens__blood__30c2a6fd | 10x 5' v2 human | umi human | 22186 | ok, scored | 1.28/1.28 | 1.16/1.36 | 16044.6 |
| cxg__10x-5-v2__homo-sapiens__blood__bd65a70f | 10x 5' v2 human | umi human | 21973 | ok, scored | 1.48/1.18 | 1.35/1.33 | 13681.3 |
| cxg__10x-multiome__homo-sapiens__amygdala__835fa5c4 | 10x multiome human | umi human | 27062 | ok, scored | 1.26/1.04 | 1.24/1.05 | 6902.7 |
| cxg__10x-multiome__homo-sapiens__anterior-cingulate-cor | 10x multiome human | umi human | 27056 | ok, scored | 1.27/1.27 | 1.24/1.21 | 8182.6 |
| cxg__10x-multiome__homo-sapiens__frontal-lobe__01b84709 | 10x multiome human | umi human | 27062 | ok, scored | 1.23/1.23 | 1.18/1.18 | 5037.7 |
| cxg__10x-multiome__homo-sapiens__kidney__867757c1 | 10x multiome human | umi human | 25431 | ok, scored | 1.16/1.16 | 1.26/1.26 | 2452.0 |
| cxg__10x-multiome__homo-sapiens__retina__322c816f | 10x multiome human | umi human | 27056 | ok, scored | 1.10/1.10 | 1.18/1.18 | 3242.9 |
| cxg__10x-multiome__mus-musculus__kidney__25818bf7 | 10x multiome mouse | umi mouse | 28097 | ok, scored | 1.10/1.10 | 1.31/1.31 | 11465.9 |
| cxg__10x-next-gem-flex-v1__homo-sapiens__endometrium__f | 10x Next GEM Flex v1 human | umi human | 26012 | ok, scored | 1.47/1.15 | 1.35/1.35 | 5333.3 |
| cxg__10x-next-gem-flex-v1__homo-sapiens__parenchyma-of- | 10x Next GEM Flex v1 human | umi human | 25951 | ok, scored | 1.56/1.50 | 1.62/1.62 | 921.8 |
| cxg__asteria-scrna-seq-kit__homo-sapiens__blood__27d395 | Asteria scRNA-seq kit human | umi human | 16045 | ok, scored | 1.47/1.47 | 1.47/1.35 | 347.8 |
| cxg__bd-rhapsody-whole-transcriptome-analysis__homo-sap | BD Rhapsody Whole Transcriptome Analysis human | umi human | 24570 | ok, scored | 1.47/1.33 | 1.50/1.41 | 6057.9 |
| cxg__mars-seq__homo-sapiens__bone-marrow__e22482ee | MARS-seq human | umi human | 25940 | ok, scored | 1.33/1.33 | 1.24/1.24 | 353.7 |
| cxg__mct-seq__homo-sapiens__brodmann-1909-area-46__43aa | mCT-seq human | umi human | 27102 | ok, scored | 1.27/1.15 | 1.16/1.17 | 5030.5 |
| cxg__sci-plex__homo-sapiens__k-562__38a78e91 | sci-Plex human | umi human | 27082 | ok, scored | 1.19/1.19 | 1.12/1.12 | 23477.6 |
| cxg__seq-well-s3__homo-sapiens__nasopharynx__13b61a7d | Seq-Well S3 human | umi human | 27027 | ok, scored | 2.16/2.07 | 1.52/1.44 | 1532.3 |
| cxg__seq-well__homo-sapiens__blood__456e8b9b | Seq-Well human | umi human | 24972 | ok, scored | 1.38/1.26 | 1.20/1.27 | 1782.9 |
| cxg__smart-seq-v4__mus-musculus__primary-motor-cortex__ | Smart-seq v4 mouse | full_length mouse | 28258 | ok, scored | 1.06/1.06 | 0.99/1.07 | 491.9 |
| cxg__smart-seq2__homo-sapiens__blood__e04daea4 | Smart-seq2 human | full_length human | 25668 | ok, scored | 1.51/1.51 | 1.37/1.37 | 299.3 |
| cxg__smart-seq2__homo-sapiens__heart-left-ventricle__07 | Smart-seq2 human | full_length human | 27048 | ok, scored | 1.21/1.21 | 1.23/1.23 | 241.9 |
| cxg__smart-seq__mus-musculus__primary-visual-cortex__3a | Smart-seq mouse | full_length mouse | 27432 | ok, scored | 1.09/1.09 | 1.15/1.15 | 167.5 |
| cxg__split-seq__homo-sapiens__cerebellum__b2dda353 | SPLiT-seq human | umi human | 25962 | ok, scored | 1.19/1.19 | 1.08/1.14 | 8874.7 |
