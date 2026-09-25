# Engine mass test on real datasets - 17 datasets

- ran end to end: 17
- refused as expected (not raw counts): 0
- errors: 0

## Detection vs curated metadata

| check | n | correct |
|---|---|---|
| organism | 17 | 71% |
| technology | 17 | 100% |

## Reference recovery (CollecTRI evaluation half, AUPRC ratio; 1 = random)

| model | n | auto > baseline | auto = baseline | mean auto | mean baseline | mean rules |
|---|---|---|---|---|---|---|
| pearson_signed | 11 | 18% | 82% | 1.27 | 1.27 | 1.15 |
| spearman_signed | 11 | 9% | 82% | 1.34 | 1.33 | 1.22 |
| clr | 11 | 36% | 45% | 1.33 | 1.30 | 1.27 |

## Per dataset

| dataset | truth | detected | ref edges | outcome | Pearson auto/base | CLR auto/base | s |
|---|---|---|---|---|---|---|---|
| cxg__10x-3-v2__homo-sapiens__entorhinal-cortex__9f1049a | 10x 3' v2 human | umi human | 27002 | ok, scored | 1.21/1.21 | 1.24/1.13 | 692.0 |
| cxg__10x-3-v2__homo-sapiens__superior-frontal-gyrus__bd | 10x 3' v2 human | umi human | 27002 | ok, scored | 1.39/1.39 | 1.17/1.17 | 504.2 |
| cxg__10x-3-v2__microcebus-murinus__hypothalamus-pituita | 10x 3' v2 Microcebus murinus | umi mouse | 0 | ok (2 contexts) | - | - | 47.4 |
| cxg__10x-3-v2__microcebus-murinus__uterus__f69adac8 | 10x 3' v2 Microcebus murinus | umi mouse | 0 | ok (3 contexts) | - | - | 55.8 |
| cxg__10x-3-v2__mus-musculus__large-intestine__de4e7a0c | 10x 3' v2 mouse | umi mouse | 26662 | ok, scored | 1.26/1.21 | 1.44/1.44 | 734.4 |
| cxg__10x-3-v2__mus-musculus__spinal-cord__989a76f9 | 10x 3' v2 mouse | umi mouse | 24615 | ok, scored | 1.25/1.25 | 1.17/1.12 | 789.8 |
| cxg__10x-3-v3__callithrix-jacchus__hippocampal-formatio | 10x 3' v3 Callithrix jacchus | umi None | 0 | ok (4 contexts) | - | - | 106.4 |
| cxg__10x-3-v3__callithrix-jacchus__primary-motor-cortex | 10x 3' v3 Callithrix jacchus | umi human | 0 | ok (4 contexts) | - | - | 116.5 |
| cxg__10x-3-v3__homo-sapiens__liver__d7476ae2 | 10x 3' v3 human | umi human | 16898 | ok, scored | - | - | 10.4 |
| cxg__10x-3-v3__homo-sapiens__prostate-gland__182f8181 | 10x 3' v3 human | umi human | 27063 | ok, scored | 1.18/1.18 | 1.36/1.36 | 37590.1 |
| cxg__10x-3-v3__macaca-mulatta__primary-motor-cortex__1f | 10x 3' v3 Macaca mulatta | umi human | 0 | ok (4 contexts) | - | - | 200.3 |
| cxg__10x-3-v3__mus-musculus__cervical-spinal-cord__9fcd | 10x 3' v3 mouse | umi mouse | 28159 | ok, scored | 0.97/0.97 | 1.23/1.23 | 607.2 |
| cxg__10x-3-v3__mus-musculus__diencephalon__6de332e1 | 10x 3' v3 mouse | umi mouse | 28191 | ok, scored | 1.21/1.21 | 1.20/1.22 | 833.4 |
| cxg__10x-5-v1__homo-sapiens__lung__703f00e6 | 10x 5' v1 human | umi human | 18058 | ok, scored | 1.30/1.30 | 1.45/1.45 | 370.9 |
| cxg__10x-5-v1__homo-sapiens__lung__cdefb878 | 10x 5' v1 human | umi human | 21785 | ok, scored | 1.27/1.24 | 1.26/1.29 | 830.5 |
| cxg__10x-5-v2__homo-sapiens__blood__ee782eb7 | 10x 5' v2 human | umi human | 22729 | ok, scored | 1.66/1.66 | 1.78/1.68 | 2806.4 |
| cxg__10x-5-v2__mus-musculus__caecum__0a1ad9d5 | 10x 5' v2 mouse | umi mouse | 23848 | ok, scored | 1.30/1.30 | 1.38/1.26 | 937.9 |
