# Engine mass test - 80 datasets, 80 completed, 0 errors

## Overall (paired over datasets = independent units)

| comparison | n | win rate | median relative gain | mean A | mean B | Wilcoxon p |
|---|---|---|---|---|---|---|
| auto vs baseline | 76 | 79% | +38.7% | 56.81 | 43.29 | 8e-12 |
| rules vs baseline | 74 | 78% | +22.9% | 51.37 | 43.50 | 7e-08 |
| auto_noref vs baseline | 76 | 14% | +0.0% | 45.17 | 43.29 | 0.16 |
| auto vs rules | 75 | 80% | +9.5% | 56.81 | 50.86 | 1.5e-11 |
| oracle vs auto | 77 | 17% | +0.0% | 56.73 | 56.26 | 0.0015 |

## Per model: pipeline chosen for that model vs baseline (evaluation half)

| model | n | auto wins | median relative gain | mean auto | mean baseline | mean rules | mean oracle | Wilcoxon p |
|---|---|---|---|---|---|---|---|---|
| pearson_signed | 76 | 82% | +76.6% | 64.99 | 41.50 | 59.84 | 65.40 | 4.8e-12 |
| spearman_signed | 76 | 80% | +80.5% | 66.23 | 40.73 | 61.40 | 66.57 | 5.1e-12 |
| clr | 76 | 12% | +0.0% | 48.50 | 47.66 | 32.88 | 49.16 | 0.65 |

## Per scenario and model: auto-for-model vs baseline (mean AUPRC ratio)

| scenario | n | pearson_signed auto / base | spearman_signed auto / base | clr auto / base |
|---|---|---|---|---|
| ambient | 3 | 65.7 / 36.1 | 69.8 / 34.2 | 43.2 / 43.2 |
| batch_balanced | 3 | 79.3 / 49.6 | 78.3 / 48.7 | 62.7 / 62.7 |
| batch_confounded | 3 | 46.1 / 29.6 | 44.8 / 27.9 | 46.8 / 46.1 |
| clean | 3 | 76.9 / 40.9 | 77.5 / 39.1 | 56.7 / 56.7 |
| damaged | 3 | 76.6 / 39.7 | 77.4 / 38.5 | 56.3 / 56.3 |
| doublets | 3 | 76.2 / 41.6 | 76.9 / 39.5 | 56.9 / 56.9 |
| empty_droplets | 3 | 77.8 / 40.5 | 78.3 / 38.9 | 56.9 / 56.9 |
| ensembl_ids | 3 | 79.5 / 42.5 | 79.7 / 41.4 | 59.0 / 60.8 |
| full_length | 3 | 73.7 / 47.8 | 75.1 / 44.8 | 68.6 / 68.6 |
| low_depth | 3 | 48.4 / 24.8 | 54.4 / 24.0 | 29.4 / 29.4 |
| misleading_readme | 3 | 76.9 / 40.9 | 77.5 / 39.1 | 56.7 / 56.7 |
| mixed_problems | 3 | 49.7 / 39.2 | 57.4 / 39.1 | 38.6 / 38.6 |
| mouse_ids | 3 | 76.9 / 40.9 | 77.5 / 39.1 | 56.7 / 56.7 |
| no_readme | 3 | 76.9 / 40.9 | 77.5 / 39.1 | 56.7 / 56.7 |
| random | 23 | 55.7 / 42.8 | 56.5 / 43.4 | 37.8 / 36.0 |
| raw_droplets | 3 | nan / nan | nan / nan | nan / nan |
| single_type | 3 | 86.3 / 86.5 | 86.1 / 86.1 | 83.1 / 84.0 |
| small | 3 | 52.7 / 53.4 | 50.2 / 50.5 | 45.0 / 46.1 |
| strong_composition | 3 | 60.4 / 19.3 | 61.7 / 21.1 | 39.4 / 39.4 |
| very_sparse | 3 | 45.5 / 23.0 | 51.8 / 22.3 | 27.2 / 27.2 |

## Per scenario: mean AUPRC ratio on the evaluation half (mean over models)

| scenario | n paired/total | baseline | rules | auto | auto_noref | oracle | auto wins vs baseline |
|---|---|---|---|---|---|---|---|
| ambient | 3/3 | 37.84 | 47.95 | 56.33 | 37.84 | 56.33 | 100% |
| batch_balanced | 3/3 | 53.67 | 63.26 | 69.49 | 62.81 | 69.49 | 100% |
| batch_confounded | 1/3 | 32.69 | 58.41 | 66.81 | 66.81 | 66.81 | 33% |
| clean | 3/3 | 45.56 | 62.27 | 67.41 | 45.56 | 67.41 | 100% |
| damaged | 3/3 | 44.83 | 60.75 | 66.44 | 44.83 | 66.44 | 100% |
| doublets | 3/3 | 46.03 | 62.29 | 66.54 | 46.03 | 66.54 | 100% |
| empty_droplets | 3/3 | 45.47 | 61.59 | 68.13 | 45.47 | 68.13 | 100% |
| ensembl_ids | 3/3 | 48.25 | 65.61 | 70.45 | 48.25 | 70.45 | 100% |
| full_length | 3/3 | 53.73 | 59.83 | 63.08 | 53.73 | 63.95 | 100% |
| low_depth | 3/3 | 26.05 | 30.53 | 37.27 | 26.05 | 37.27 | 100% |
| misleading_readme | 3/3 | 45.56 | 62.27 | 67.41 | 45.56 | 67.41 | 100% |
| mixed_problems | 3/3 | 38.95 | 35.39 | 39.71 | 38.95 | 42.98 | 67% |
| mouse_ids | 3/3 | 45.56 | 62.27 | 67.41 | 45.56 | 67.41 | 100% |
| no_readme | 3/3 | 45.56 | 62.27 | 67.41 | 45.56 | 67.41 | 100% |
| random | 22/23 | 40.75 | 42.47 | 50.25 | 44.47 | 51.20 | 68% |
| raw_droplets | 0/3 | nan | nan | nan | nan | nan | nan% |
| single_type | 3/3 | 85.53 | 85.67 | 85.36 | 85.75 | 85.83 | 0% |
| small | 3/3 | 50.00 | 49.55 | 49.46 | 49.55 | 50.15 | 0% |
| strong_composition | 3/3 | 26.60 | 37.97 | 48.26 | 26.60 | 48.26 | 100% |
| very_sparse | 3/3 | 24.18 | 26.89 | 34.55 | 24.18 | 34.55 | 100% |

## Detection accuracy (share of datasets where the engine got it right)

| check | n | accuracy | scenarios where it failed |
|---|---|---|---|
| batch | 5 | 100% | - |
| batch_detected | 80 | 100% | - |
| cell_calling_best | 3 | 100% | - |
| conflict | 8 | 100% | - |
| context_scope | 63 | 100% | - |
| doublets | 13 | 100% | - |
| doublets_removed | 57 | 100% | - |
| identifier_mapping_best | 3 | 100% | - |
| normalization | 13 | 100% | - |
| organism | 80 | 100% | - |
| technology | 80 | 100% | - |

## Cleaning accuracy (mean over datasets)

| measure | baseline | auto |
|---|---|---|
| bad_removed | 69.2% | 71.2% |
| ok_retained | 91.1% | 96.2% |
| removed_damaged | 100.0% | 100.0% |
| removed_doublet | 13.6% | 13.3% |
| removed_empty | 99.3% | 100.0% |

Median runtime per dataset: 334 s.
