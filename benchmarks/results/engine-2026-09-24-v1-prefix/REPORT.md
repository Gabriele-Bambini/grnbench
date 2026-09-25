# Engine mass test - 80 datasets, 72 completed, 8 errors

## Overall (paired over datasets = independent units)

| comparison | n | win rate | median relative gain | mean A | mean B | Wilcoxon p |
|---|---|---|---|---|---|---|
| auto vs baseline | 71 | 83% | +41.3% | 57.88 | 43.27 | 1.6e-11 |
| rules vs baseline | 69 | 81% | +24.3% | 52.02 | 43.50 | 6.8e-08 |
| auto_noref vs baseline | 71 | 14% | +0.0% | 43.96 | 43.27 | 0.2 |
| auto vs rules | 70 | 83% | +10.5% | 57.88 | 51.43 | 1.7e-11 |
| oracle vs auto | 72 | 14% | +0.0% | 57.53 | 57.26 | 0.0051 |

## Per model: pipeline chosen for that model vs baseline (evaluation half)

| model | n | auto wins | median relative gain | mean auto | mean baseline | mean rules | mean oracle | Wilcoxon p |
|---|---|---|---|---|---|---|---|---|
| pearson_signed | 71 | 86% | +78.2% | 66.11 | 41.29 | 60.75 | 66.46 | 7.6e-12 |
| spearman_signed | 71 | 86% | +88.1% | 67.38 | 40.42 | 62.43 | 67.74 | 1.1e-11 |
| clr | 71 | 15% | +0.0% | 49.53 | 48.11 | 32.87 | 50.00 | 0.061 |

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
| misleading_readme | 3 | 77.5 / 40.9 | 78.0 / 39.1 | 56.7 / 56.7 |
| mixed_problems | 3 | 49.7 / 39.2 | 57.4 / 39.1 | 38.6 / 38.6 |
| mouse_ids | 3 | 76.9 / 40.9 | 77.5 / 39.1 | 56.7 / 56.7 |
| no_readme | 3 | 76.9 / 40.9 | 77.5 / 39.1 | 56.7 / 56.7 |
| random | 21 | 56.7 / 44.0 | 56.9 / 44.1 | 39.8 / 36.2 |
| single_type | 3 | 86.3 / 86.5 | 86.1 / 86.1 | 83.1 / 84.0 |
| strong_composition | 3 | 60.4 / 19.3 | 61.7 / 21.1 | 39.4 / 39.4 |
| very_sparse | 3 | 45.5 / 23.0 | 51.8 / 22.3 | 27.2 / 27.2 |

## Per scenario: mean AUPRC ratio on the evaluation half (mean over models)

| scenario | n | baseline | rules | auto | auto_noref | oracle | auto wins vs baseline |
|---|---|---|---|---|---|---|---|
| ambient | 3 | 37.84 | 47.95 | 56.33 | 37.84 | 56.33 | 100% |
| batch_balanced | 3 | 53.67 | 63.26 | 69.49 | 62.81 | 69.49 | 100% |
| batch_confounded | 3 | 34.55 | 58.41 | 45.92 | 45.92 | 45.92 | 33% |
| clean | 3 | 45.56 | 62.27 | 67.41 | 45.56 | 67.41 | 100% |
| damaged | 3 | 44.83 | 60.75 | 66.44 | 44.83 | 66.44 | 100% |
| doublets | 3 | 46.03 | 62.29 | 66.54 | 46.03 | 66.54 | 100% |
| empty_droplets | 3 | 45.47 | 61.59 | 68.13 | 45.47 | 68.13 | 100% |
| ensembl_ids | 3 | 48.25 | 65.61 | 70.45 | 48.25 | 70.45 | 100% |
| full_length | 3 | 53.73 | 59.83 | 63.08 | 53.73 | 63.95 | 100% |
| low_depth | 3 | 26.05 | 30.53 | 37.27 | 26.05 | 37.27 | 100% |
| misleading_readme | 3 | 45.56 | 63.20 | 68.15 | 45.56 | 68.15 | 100% |
| mixed_problems | 3 | 38.95 | 35.39 | 39.71 | 38.95 | 42.98 | 67% |
| mouse_ids | 3 | 45.56 | 62.27 | 67.41 | 45.56 | 67.41 | 100% |
| no_readme | 3 | 45.56 | 62.27 | 67.41 | 45.56 | 67.41 | 100% |
| random | 21 | 41.43 | 41.84 | 50.33 | 39.33 | 50.58 | 70% |
| single_type | 3 | 85.53 | 85.67 | 85.36 | 85.75 | 85.83 | 0% |
| strong_composition | 3 | 26.60 | 37.97 | 48.26 | 26.60 | 48.26 | 100% |
| very_sparse | 3 | 24.18 | 26.89 | 34.55 | 24.18 | 34.55 | 100% |

## Detection accuracy (share of datasets where the engine got it right)

| check | n | accuracy | scenarios where it failed |
|---|---|---|---|
| batch | 5 | 100% | - |
| batch_detected | 72 | 100% | - |
| cell_calling_best | 0 | nan% | - |
| conflict | 7 | 100% | - |
| context_scope | 59 | 100% | - |
| doublets | 12 | 100% | - |
| doublets_removed | 54 | 93% | misleading_readme, random |
| identifier_mapping_best | 3 | 100% | - |
| normalization | 12 | 100% | - |
| organism | 72 | 100% | - |
| technology | 72 | 100% | - |

## Cleaning accuracy (mean over datasets)

| measure | baseline | auto |
|---|---|---|
| bad_removed | 74.7% | 73.4% |
| ok_retained | 89.1% | 96.3% |
| removed_damaged | 100.0% | 100.0% |
| removed_doublet | 18.3% | 14.8% |
| removed_empty | 98.3% | 99.6% |

## Errors

- `random__s10`: KeyError: 'edge_stability'
- `random__s2`: KeyError: 'edge_stability'
- `raw_droplets__s0`: KeyError: 'edge_stability'
- `raw_droplets__s1`: KeyError: 'edge_stability'
- `raw_droplets__s2`: KeyError: 'edge_stability'
- `small__s0`: KeyError: 'edge_stability'
- `small__s1`: KeyError: 'edge_stability'
- `small__s2`: KeyError: 'edge_stability'

Median runtime per dataset: 273 s.
