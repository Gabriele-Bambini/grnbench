# Engine mass test - 740 planned datasets

- completed: 740
- expected refusals: 0
- permanent errors: 0
- missing or retryable: 0

## Overall (paired over datasets = independent units)

| comparison | n | win rate | median relative gain | mean A | mean B | Wilcoxon p |
|---|---|---|---|---|---|---|
| auto vs baseline | 700 | 96% | +252.3% | 68.64 | 21.65 | 4.8e-114 |
| rules vs baseline | 700 | 97% | +241.5% | 67.86 | 21.65 | 3.7e-115 |
| auto_noref vs baseline | 700 | 45% | +0.0% | 43.13 | 21.65 | 4.1e-55 |
| auto vs rules | 740 | 28% | +0.0% | 65.56 | 64.81 | 0.00069 |
| oracle vs auto | 740 | 61% | +1.2% | 67.02 | 65.56 | 2.8e-76 |

## Per model: pipeline chosen for that model vs baseline (evaluation half)

| model | n | auto wins | median relative gain | mean auto | mean baseline | mean rules | mean oracle | Wilcoxon p |
|---|---|---|---|---|---|---|---|---|
| pearson_signed | 700 | 96% | +231.9% | 73.82 | 24.24 | 73.57 | 75.60 | 5.8e-114 |
| spearman_signed | 700 | 96% | +256.3% | 74.15 | 23.06 | 73.88 | 75.82 | 1.8e-113 |
| clr | 700 | 95% | +243.7% | 57.89 | 17.66 | 56.12 | 59.35 | 7.1e-113 |

## Per scenario and model: auto-for-model vs baseline (mean AUPRC ratio)

| scenario | n | pearson_signed auto / base | spearman_signed auto / base | clr auto / base |
|---|---|---|---|---|
| sf025_ambient_extreme | 5 | 70.8 / 11.3 | 72.8 / 10.8 | 36.1 / 5.3 |
| sf025_ambient_mild | 5 | 81.2 / 14.7 | 81.1 / 13.9 | 64.5 / 8.3 |
| sf025_ambient_moderate | 5 | 80.1 / 14.0 | 81.1 / 13.3 | 61.5 / 7.9 |
| sf025_ambient_severe | 5 | 76.9 / 12.9 | 78.1 / 12.2 | 49.4 / 6.3 |
| sf025_batch_confounding_extreme | 5 | 82.4 / 3.4 | 82.3 / 3.0 | 70.0 / 8.0 |
| sf025_batch_confounding_mild | 5 | 82.8 / 10.5 | 81.8 / 9.1 | 64.6 / 10.7 |
| sf025_batch_confounding_moderate | 5 | 85.7 / 6.1 | 85.3 / 5.2 | 69.4 / 8.1 |
| sf025_batch_confounding_severe | 5 | 81.7 / 4.2 | 81.9 / 3.7 | 69.6 / 7.8 |
| sf025_batch_depth_extreme | 5 | 81.8 / 17.4 | 82.1 / 16.3 | 75.9 / 12.7 |
| sf025_batch_depth_mild | 5 | 81.4 / 16.0 | 81.6 / 15.0 | 75.1 / 11.3 |
| sf025_batch_depth_moderate | 5 | 82.6 / 16.5 | 82.6 / 15.4 | 76.9 / 11.9 |
| sf025_batch_depth_severe | 5 | 82.0 / 17.0 | 82.2 / 15.8 | 77.5 / 12.8 |
| sf025_cells_extreme | 5 | 6.1 / 6.1 | 5.7 / 5.7 | 2.2 / 2.2 |
| sf025_cells_mild | 5 | 82.0 / 14.4 | 82.0 / 13.4 | 61.0 / 8.2 |
| sf025_cells_moderate | 5 | 78.4 / 12.1 | 78.5 / 11.4 | 42.6 / 5.4 |
| sf025_cells_severe | 5 | 8.5 / 8.4 | 8.1 / 8.1 | 3.4 / 3.3 |
| sf025_clean_baseline | 5 | 82.3 / 15.5 | 82.2 / 14.7 | 72.2 / 9.3 |
| sf025_damaged_extreme | 5 | 79.9 / 14.8 | 80.1 / 14.0 | 67.4 / 8.8 |
| sf025_damaged_mild | 5 | 81.7 / 15.4 | 81.9 / 14.6 | 73.2 / 9.3 |
| sf025_damaged_moderate | 5 | 81.7 / 15.4 | 82.0 / 14.5 | 71.2 / 9.1 |
| sf025_damaged_severe | 5 | 81.9 / 15.4 | 81.6 / 14.5 | 69.8 / 9.2 |
| sf025_depth_extreme | 5 | 27.7 / nan | 30.4 / nan | 4.5 / nan |
| sf025_depth_mild | 5 | 78.1 / 12.7 | 78.9 / 12.1 | 49.8 / 5.4 |
| sf025_depth_moderate | 5 | 71.4 / 10.3 | 73.2 / 10.1 | 24.7 / 3.6 |
| sf025_depth_severe | 5 | 51.6 / 6.5 | 56.3 / 6.4 | 9.2 / 2.6 |
| sf025_doublets_extreme | 5 | 81.6 / 15.5 | 81.5 / 14.7 | 68.5 / 10.1 |
| sf025_doublets_mild | 5 | 82.2 / 15.6 | 82.1 / 14.8 | 71.7 / 9.4 |
| sf025_doublets_moderate | 5 | 82.3 / 15.6 | 82.1 / 14.8 | 70.7 / 9.5 |
| sf025_doublets_severe | 5 | 82.1 / 15.6 | 82.0 / 14.8 | 69.2 / 9.7 |
| sf025_mixed_mixed_1 | 5 | 76.0 / 8.0 | 70.6 / 7.4 | 31.7 / 3.9 |
| sf025_mixed_mixed_2 | 5 | 44.8 / 7.4 | 41.2 / 7.6 | 7.5 / 2.2 |
| sf025_mixed_mixed_3 | 5 | 4.3 / 1.7 | 6.5 / 1.8 | 1.3 / 1.1 |
| sf025_mixed_mixed_4 | 5 | 1.3 / nan | 1.2 / nan | 1.0 / nan |
| sf025_type_imbalance_extreme | 5 | 80.4 / 23.9 | 81.2 / 27.2 | 81.4 / 21.0 |
| sf025_type_imbalance_mild | 5 | 81.5 / 16.0 | 81.7 / 15.1 | 78.0 / 10.2 |
| sf025_type_imbalance_moderate | 5 | 82.4 / 18.0 | 82.7 / 17.8 | 79.6 / 13.4 |
| sf025_type_imbalance_severe | 5 | 82.7 / 20.7 | 82.6 / 21.8 | 81.2 / 17.6 |
| sf050_ambient_extreme | 5 | 72.8 / 17.0 | 74.7 / 15.9 | 34.5 / 8.4 |
| sf050_ambient_mild | 5 | 81.0 / 20.5 | 80.9 / 19.0 | 66.3 / 12.3 |
| sf050_ambient_moderate | 5 | 79.7 / 20.0 | 80.6 / 18.6 | 60.0 / 11.2 |
| sf050_ambient_severe | 5 | 77.7 / 18.6 | 78.5 / 17.3 | 46.4 / 9.9 |
| sf050_batch_confounding_extreme | 5 | 83.0 / 6.0 | 83.3 / 5.2 | 72.5 / 13.2 |
| sf050_batch_confounding_mild | 5 | 82.7 / 14.0 | 80.7 / 12.1 | 64.1 / 13.4 |
| sf050_batch_confounding_moderate | 5 | 80.3 / 12.8 | 79.1 / 10.9 | 61.6 / 13.9 |
| sf050_batch_confounding_severe | 5 | 80.9 / 7.7 | 80.6 / 6.5 | 69.0 / 11.2 |
| sf050_batch_depth_extreme | 5 | 82.0 / 22.2 | 83.2 / 20.8 | 76.7 / 16.0 |
| sf050_batch_depth_mild | 5 | 83.1 / 20.9 | 83.1 / 19.6 | 76.3 / 15.0 |
| sf050_batch_depth_moderate | 5 | 83.8 / 21.3 | 83.7 / 19.9 | 78.5 / 16.0 |
| sf050_batch_depth_severe | 5 | 81.7 / 21.3 | 82.8 / 20.0 | 78.4 / 16.2 |
| sf050_cells_extreme | 5 | 16.7 / 16.7 | 16.5 / 16.5 | 4.7 / 4.8 |
| sf050_cells_mild | 5 | 79.7 / 23.7 | 81.3 / 22.1 | 61.0 / 15.4 |
| sf050_cells_moderate | 5 | 79.4 / 19.5 | 80.1 / 18.2 | 45.3 / 10.6 |
| sf050_cells_severe | 5 | 18.9 / 18.9 | 17.8 / 17.8 | 8.1 / 7.9 |
| sf050_clean_baseline | 5 | 80.8 / 21.1 | 80.9 / 19.7 | 72.1 / 13.6 |
| sf050_damaged_extreme | 5 | 81.0 / 20.7 | 81.0 / 19.1 | 68.4 / 12.8 |
| sf050_damaged_mild | 5 | 81.5 / 21.1 | 81.2 / 19.7 | 72.4 / 13.5 |
| sf050_damaged_moderate | 5 | 80.8 / 20.8 | 80.9 / 19.4 | 70.1 / 13.2 |
| sf050_damaged_severe | 5 | 80.3 / 20.8 | 80.5 / 19.4 | 70.4 / 13.3 |
| sf050_depth_extreme | 5 | 29.6 / nan | 31.2 / nan | 4.7 / nan |
| sf050_depth_mild | 5 | 78.5 / 17.8 | 78.9 / 16.7 | 49.1 / 9.4 |
| sf050_depth_moderate | 5 | 68.4 / 14.7 | 71.7 / 14.4 | 27.1 / 6.5 |
| sf050_depth_severe | 5 | 52.3 / 10.5 | 57.1 / 10.4 | 12.5 / 3.8 |
| sf050_doublets_extreme | 5 | 80.4 / 21.3 | 80.4 / 19.8 | 68.2 / 14.5 |
| sf050_doublets_mild | 5 | 80.9 / 21.2 | 81.0 / 19.7 | 70.9 / 13.7 |
| sf050_doublets_moderate | 5 | 81.1 / 21.2 | 81.0 / 19.8 | 70.2 / 13.9 |
| sf050_doublets_severe | 5 | 81.1 / 21.2 | 81.0 / 19.8 | 69.6 / 14.0 |
| sf050_mixed_mixed_1 | 5 | 71.4 / 14.3 | 68.0 / 13.3 | 31.2 / 6.9 |
| sf050_mixed_mixed_2 | 5 | 42.6 / 9.5 | 43.9 / 10.1 | 5.7 / 3.0 |
| sf050_mixed_mixed_3 | 5 | 6.6 / 2.9 | 7.5 / 3.3 | 1.4 / 1.2 |
| sf050_mixed_mixed_4 | 5 | 1.5 / nan | 1.3 / nan | 1.2 / nan |
| sf050_type_imbalance_extreme | 5 | 80.0 / 32.8 | 80.0 / 36.3 | 78.4 / 26.4 |
| sf050_type_imbalance_mild | 5 | 81.2 / 21.9 | 80.5 / 20.6 | 75.3 / 14.5 |
| sf050_type_imbalance_moderate | 5 | 80.2 / 24.6 | 79.9 / 23.9 | 77.5 / 17.9 |
| sf050_type_imbalance_severe | 5 | 80.7 / 27.9 | 80.7 / 28.8 | 78.6 / 22.4 |
| sf075_ambient_extreme | 5 | 73.8 / 25.8 | 75.3 / 24.0 | 42.3 / 14.5 |
| sf075_ambient_mild | 5 | 83.7 / 29.0 | 83.9 / 27.1 | 71.1 / 20.4 |
| sf075_ambient_moderate | 5 | 83.0 / 28.2 | 83.2 / 26.3 | 62.9 / 19.1 |
| sf075_ambient_severe | 5 | 80.7 / 27.4 | 81.0 / 25.5 | 52.8 / 17.0 |
| sf075_batch_confounding_extreme | 5 | 81.4 / 12.3 | 81.4 / 9.8 | 71.6 / 18.2 |
| sf075_batch_confounding_mild | 5 | 83.9 / 29.8 | 82.3 / 26.7 | 65.8 / 28.5 |
| sf075_batch_confounding_moderate | 5 | 80.7 / 21.6 | 79.4 / 18.6 | 61.0 / 21.3 |
| sf075_batch_confounding_severe | 5 | 81.3 / 15.0 | 81.2 / 13.0 | 71.5 / 20.6 |
| sf075_batch_depth_extreme | 5 | 83.8 / 37.1 | 84.9 / 34.2 | 79.6 / 28.2 |
| sf075_batch_depth_mild | 5 | 85.8 / 34.9 | 85.8 / 32.1 | 79.9 / 27.7 |
| sf075_batch_depth_moderate | 5 | 84.8 / 35.2 | 85.4 / 32.5 | 81.1 / 28.5 |
| sf075_batch_depth_severe | 5 | 84.6 / 36.3 | 85.8 / 33.4 | 81.6 / 29.0 |
| sf075_cells_extreme | 5 | 17.1 / 16.7 | 16.1 / 15.7 | 4.3 / 3.9 |
| sf075_cells_mild | 5 | 81.0 / 31.0 | 81.1 / 28.9 | 60.1 / 20.6 |
| sf075_cells_moderate | 5 | 77.1 / 32.5 | 77.6 / 30.3 | 42.8 / 20.1 |
| sf075_cells_severe | 5 | 29.5 / 29.5 | 28.0 / 28.0 | 14.5 / 14.4 |
| sf075_clean_baseline | 5 | 84.6 / 28.9 | 84.8 / 27.4 | 76.7 / 22.0 |
| sf075_damaged_extreme | 5 | 84.1 / 28.8 | 84.1 / 27.2 | 71.0 / 21.4 |
| sf075_damaged_mild | 5 | 83.9 / 28.8 | 84.3 / 27.2 | 74.7 / 22.1 |
| sf075_damaged_moderate | 5 | 84.0 / 28.9 | 84.2 / 27.3 | 73.7 / 21.9 |
| sf075_damaged_severe | 5 | 83.9 / 28.9 | 84.3 / 27.3 | 75.3 / 21.9 |
| sf075_depth_extreme | 5 | 30.4 / nan | 34.1 / nan | 4.4 / nan |
| sf075_depth_mild | 5 | 81.5 / 25.8 | 83.1 / 24.6 | 57.5 / 17.4 |
| sf075_depth_moderate | 5 | 72.6 / 22.1 | 75.8 / 21.8 | 31.6 / 12.0 |
| sf075_depth_severe | 5 | 53.4 / 17.2 | 58.5 / 17.3 | 10.0 / 5.9 |
| sf075_doublets_extreme | 5 | 84.1 / 29.3 | 84.2 / 27.8 | 72.9 / 22.5 |
| sf075_doublets_mild | 5 | 84.5 / 29.0 | 84.8 / 27.5 | 75.7 / 22.0 |
| sf075_doublets_moderate | 5 | 84.4 / 29.0 | 84.7 / 27.5 | 74.5 / 22.0 |
| sf075_doublets_severe | 5 | 84.1 / 29.1 | 84.4 / 27.6 | 73.3 / 22.0 |
| sf075_mixed_mixed_1 | 5 | 75.9 / 22.6 | 71.2 / 20.6 | 31.8 / 11.8 |
| sf075_mixed_mixed_2 | 5 | 43.8 / 20.2 | 50.7 / 20.2 | 6.8 / 5.7 |
| sf075_mixed_mixed_3 | 5 | 7.0 / 3.3 | 7.7 / 4.3 | 1.4 / 1.5 |
| sf075_mixed_mixed_4 | 5 | 2.3 / nan | 1.4 / nan | 1.1 / nan |
| sf075_type_imbalance_extreme | 5 | 83.3 / 45.6 | 83.2 / 51.4 | 84.5 / 42.1 |
| sf075_type_imbalance_mild | 5 | 83.3 / 30.9 | 83.1 / 29.8 | 78.3 / 23.9 |
| sf075_type_imbalance_moderate | 5 | 83.6 / 35.2 | 83.4 / 35.7 | 81.6 / 29.8 |
| sf075_type_imbalance_severe | 5 | 83.9 / 40.3 | 84.3 / 43.0 | 83.2 / 36.0 |
| sf100_ambient_extreme | 5 | 75.0 / 37.4 | 77.5 / 34.5 | 38.7 / 22.2 |
| sf100_ambient_mild | 5 | 84.4 / 40.3 | 84.5 / 37.4 | 71.9 / 31.6 |
| sf100_ambient_moderate | 5 | 83.1 / 39.9 | 83.5 / 36.7 | 63.1 / 29.3 |
| sf100_ambient_severe | 5 | 80.4 / 38.4 | 81.7 / 35.2 | 52.5 / 25.1 |
| sf100_batch_confounding_extreme | 5 | 84.9 / 11.1 | 85.2 / 9.2 | 75.2 / 18.5 |
| sf100_batch_confounding_mild | 5 | 79.8 / 35.2 | 77.8 / 31.7 | 60.6 / 29.4 |
| sf100_batch_confounding_moderate | 5 | 84.6 / 28.4 | 83.6 / 24.8 | 66.2 / 27.6 |
| sf100_batch_confounding_severe | 5 | 83.3 / 24.8 | 83.0 / 20.9 | 72.0 / 30.2 |
| sf100_batch_depth_extreme | 5 | 83.8 / 46.1 | 84.5 / 43.8 | 78.6 / 38.3 |
| sf100_batch_depth_mild | 5 | 84.6 / 43.2 | 84.5 / 40.9 | 78.5 / 37.3 |
| sf100_batch_depth_moderate | 5 | 84.7 / 44.4 | 85.0 / 41.6 | 78.8 / 38.7 |
| sf100_batch_depth_severe | 5 | 84.2 / 45.1 | 84.7 / 42.4 | 79.1 / 39.1 |
| sf100_cells_extreme | 5 | 29.1 / 28.5 | 27.8 / 27.2 | 8.2 / 7.8 |
| sf100_cells_mild | 5 | 81.7 / 37.5 | 82.7 / 35.4 | 59.6 / 27.2 |
| sf100_cells_moderate | 5 | 80.2 / 44.1 | 80.4 / 41.0 | 43.2 / 28.9 |
| sf100_cells_severe | 5 | 38.1 / 37.9 | 35.2 / 35.1 | 18.1 / 17.9 |
| sf100_clean_baseline | 5 | 84.2 / 41.1 | 84.3 / 38.2 | 76.7 / 34.2 |
| sf100_damaged_extreme | 5 | 82.6 / 40.9 | 83.2 / 38.2 | 71.2 / 33.2 |
| sf100_damaged_mild | 5 | 84.5 / 41.2 | 84.3 / 38.3 | 75.8 / 34.0 |
| sf100_damaged_moderate | 5 | 84.1 / 41.0 | 84.5 / 38.2 | 75.0 / 33.7 |
| sf100_damaged_severe | 5 | 84.7 / 41.2 | 84.8 / 38.4 | 74.8 / 33.4 |
| sf100_depth_extreme | 5 | 28.5 / nan | 32.6 / nan | 4.9 / nan |
| sf100_depth_mild | 5 | 82.3 / 36.4 | 83.2 / 33.7 | 53.9 / 25.7 |
| sf100_depth_moderate | 5 | 73.0 / 32.2 | 76.0 / 30.6 | 27.8 / 17.6 |
| sf100_depth_severe | 5 | 54.1 / 27.0 | 58.6 / 26.3 | 9.5 / 8.1 |
| sf100_doublets_extreme | 5 | 83.5 / 41.8 | 83.9 / 38.8 | 72.5 / 34.2 |
| sf100_doublets_mild | 5 | 84.1 / 41.2 | 84.3 / 38.3 | 75.5 / 34.1 |
| sf100_doublets_moderate | 5 | 83.7 / 41.3 | 84.1 / 38.4 | 74.0 / 34.0 |
| sf100_doublets_severe | 5 | 83.8 / 41.6 | 84.2 / 38.7 | 74.1 / 34.1 |
| sf100_mixed_mixed_1 | 5 | 76.7 / 33.8 | 75.2 / 31.2 | 32.7 / 18.4 |
| sf100_mixed_mixed_2 | 5 | 40.9 / 19.6 | 46.9 / 19.5 | 6.3 / 5.4 |
| sf100_mixed_mixed_3 | 5 | 6.1 / 4.5 | 6.7 / 5.6 | 1.3 / 1.4 |
| sf100_mixed_mixed_4 | 5 | 1.8 / nan | 1.5 / nan | 1.1 / nan |
| sf100_type_imbalance_extreme | 5 | 84.5 / 67.6 | 85.0 / 73.0 | 84.0 / 62.2 |
| sf100_type_imbalance_mild | 5 | 85.1 / 44.4 | 85.2 / 42.4 | 81.5 / 36.7 |
| sf100_type_imbalance_moderate | 5 | 85.2 / 51.8 | 85.9 / 52.1 | 82.6 / 46.4 |
| sf100_type_imbalance_severe | 5 | 85.3 / 59.8 | 85.5 / 63.1 | 83.3 / 54.8 |

## Per scenario: mean AUPRC ratio on the evaluation half (mean over models)

| scenario | n paired/total | baseline | rules | auto | auto_noref | oracle | auto wins vs baseline |
|---|---|---|---|---|---|---|---|
| sf025_ambient_extreme | 5/5 | 9.15 | 60.58 | 59.65 | 20.14 | 61.67 | 100% |
| sf025_ambient_mild | 5/5 | 12.30 | 76.38 | 75.44 | 50.05 | 76.87 | 100% |
| sf025_ambient_moderate | 5/5 | 11.78 | 74.08 | 73.73 | 36.13 | 75.10 | 100% |
| sf025_ambient_severe | 5/5 | 10.49 | 68.35 | 68.11 | 34.39 | 69.76 | 100% |
| sf025_batch_confounding_extreme | 5/5 | 4.80 | 81.02 | 78.22 | 79.49 | 81.11 | 100% |
| sf025_batch_confounding_mild | 5/5 | 10.14 | 77.82 | 76.40 | 78.53 | 78.53 | 100% |
| sf025_batch_confounding_moderate | 5/5 | 6.45 | 80.87 | 80.14 | 80.97 | 81.84 | 100% |
| sf025_batch_confounding_severe | 5/5 | 5.23 | 78.02 | 77.62 | 79.32 | 80.56 | 100% |
| sf025_batch_depth_extreme | 5/5 | 15.45 | 81.71 | 79.97 | 81.95 | 82.32 | 100% |
| sf025_batch_depth_mild | 5/5 | 14.08 | 81.09 | 79.37 | 51.95 | 81.28 | 100% |
| sf025_batch_depth_moderate | 5/5 | 14.60 | 82.00 | 80.80 | 79.53 | 82.30 | 100% |
| sf025_batch_depth_severe | 5/5 | 15.22 | 82.58 | 80.37 | 68.28 | 83.36 | 100% |
| sf025_cells_extreme | 5/5 | 4.65 | 4.69 | 4.69 | 4.71 | 4.71 | 60% |
| sf025_cells_mild | 5/5 | 11.97 | 75.52 | 75.41 | 76.56 | 78.12 | 100% |
| sf025_cells_moderate | 5/5 | 9.63 | 65.74 | 67.34 | 44.62 | 68.77 | 100% |
| sf025_cells_severe | 5/5 | 6.60 | 6.65 | 6.66 | 6.57 | 6.67 | 40% |
| sf025_clean_baseline | 5/5 | 13.18 | 79.43 | 78.99 | 79.87 | 80.19 | 100% |
| sf025_damaged_extreme | 5/5 | 12.55 | 77.15 | 75.38 | 37.41 | 78.10 | 100% |
| sf025_damaged_mild | 5/5 | 13.11 | 79.04 | 78.64 | 51.94 | 79.80 | 100% |
| sf025_damaged_moderate | 5/5 | 13.01 | 78.80 | 78.33 | 52.24 | 79.53 | 100% |
| sf025_damaged_severe | 5/5 | 13.03 | 78.45 | 77.90 | 38.69 | 78.73 | 100% |
| sf025_depth_extreme | 0/5 | nan | nan | nan | nan | nan | nan% |
| sf025_depth_mild | 5/5 | 10.05 | 69.53 | 68.98 | 33.61 | 70.96 | 100% |
| sf025_depth_moderate | 5/5 | 8.00 | 56.66 | 56.15 | 17.97 | 57.32 | 100% |
| sf025_depth_severe | 5/5 | 5.14 | 39.03 | 39.03 | 5.14 | 40.00 | 100% |
| sf025_doublets_extreme | 5/5 | 13.44 | 77.78 | 77.22 | 51.21 | 78.62 | 100% |
| sf025_doublets_mild | 5/5 | 13.26 | 79.16 | 78.72 | 52.01 | 80.05 | 100% |
| sf025_doublets_moderate | 5/5 | 13.31 | 78.90 | 78.41 | 51.98 | 79.88 | 100% |
| sf025_doublets_severe | 5/5 | 13.35 | 78.30 | 77.71 | 51.62 | 79.21 | 100% |
| sf025_mixed_mixed_1 | 5/5 | 6.46 | 50.39 | 59.44 | 47.14 | 59.44 | 100% |
| sf025_mixed_mixed_2 | 5/5 | 5.71 | 17.70 | 31.38 | 20.38 | 31.38 | 100% |
| sf025_mixed_mixed_3 | 5/5 | 1.52 | 1.68 | 4.02 | 2.10 | 4.02 | 100% |
| sf025_mixed_mixed_4 | 0/5 | nan | nan | nan | nan | nan | nan% |
| sf025_type_imbalance_extreme | 5/5 | 24.05 | 77.36 | 81.86 | 67.71 | 81.86 | 100% |
| sf025_type_imbalance_mild | 5/5 | 13.79 | 78.23 | 80.22 | 52.40 | 81.06 | 100% |
| sf025_type_imbalance_moderate | 5/5 | 16.39 | 75.93 | 81.54 | 52.29 | 81.54 | 100% |
| sf025_type_imbalance_severe | 5/5 | 20.02 | 73.26 | 82.19 | 68.54 | 82.21 | 100% |
| sf050_ambient_extreme | 5/5 | 13.78 | 61.72 | 60.92 | 22.50 | 61.87 | 100% |
| sf050_ambient_mild | 5/5 | 17.27 | 77.96 | 76.07 | 28.30 | 78.77 | 100% |
| sf050_ambient_moderate | 5/5 | 16.59 | 75.11 | 73.44 | 39.23 | 75.68 | 100% |
| sf050_ambient_severe | 5/5 | 15.30 | 69.26 | 67.54 | 24.85 | 69.68 | 100% |
| sf050_batch_confounding_extreme | 5/5 | 8.13 | 79.40 | 79.29 | 80.88 | 81.59 | 100% |
| sf050_batch_confounding_mild | 5/5 | 13.18 | 77.45 | 75.84 | 79.07 | 79.80 | 100% |
| sf050_batch_confounding_moderate | 5/5 | 12.52 | 74.28 | 73.48 | 76.43 | 76.43 | 100% |
| sf050_batch_confounding_severe | 5/5 | 8.47 | 75.83 | 76.85 | 62.64 | 77.59 | 100% |
| sf050_batch_depth_extreme | 5/5 | 19.65 | 82.30 | 80.67 | 69.37 | 84.13 | 100% |
| sf050_batch_depth_mild | 5/5 | 18.52 | 81.78 | 80.84 | 42.54 | 83.20 | 100% |
| sf050_batch_depth_moderate | 5/5 | 19.07 | 82.79 | 82.01 | 55.19 | 84.34 | 100% |
| sf050_batch_depth_severe | 5/5 | 19.18 | 82.41 | 80.81 | 42.85 | 84.14 | 100% |
| sf050_cells_extreme | 5/5 | 12.66 | 12.58 | 12.57 | 12.57 | 12.88 | 40% |
| sf050_cells_mild | 5/5 | 20.40 | 74.50 | 73.94 | 50.75 | 74.87 | 100% |
| sf050_cells_moderate | 5/5 | 16.08 | 67.30 | 68.06 | 25.10 | 70.33 | 100% |
| sf050_cells_severe | 5/5 | 14.84 | 14.97 | 14.88 | 14.93 | 14.99 | 20% |
| sf050_clean_baseline | 5/5 | 18.14 | 80.07 | 77.94 | 52.91 | 80.80 | 100% |
| sf050_damaged_extreme | 5/5 | 17.51 | 78.81 | 76.79 | 28.96 | 79.87 | 100% |
| sf050_damaged_mild | 5/5 | 18.08 | 80.33 | 78.38 | 40.79 | 80.83 | 100% |
| sf050_damaged_moderate | 5/5 | 17.79 | 79.73 | 77.28 | 28.98 | 80.58 | 100% |
| sf050_damaged_severe | 5/5 | 17.84 | 79.18 | 77.13 | 17.84 | 80.21 | 100% |
| sf050_depth_extreme | 0/5 | nan | nan | nan | nan | nan | nan% |
| sf050_depth_mild | 5/5 | 14.65 | 70.98 | 68.87 | 25.16 | 71.93 | 100% |
| sf050_depth_moderate | 5/5 | 11.86 | 57.96 | 55.20 | 11.86 | 58.54 | 100% |
| sf050_depth_severe | 5/5 | 8.24 | 41.37 | 40.70 | 8.24 | 42.15 | 100% |
| sf050_doublets_extreme | 5/5 | 18.56 | 78.27 | 76.34 | 29.44 | 78.73 | 100% |
| sf050_doublets_mild | 5/5 | 18.20 | 79.62 | 77.60 | 30.31 | 80.40 | 100% |
| sf050_doublets_moderate | 5/5 | 18.29 | 79.38 | 77.44 | 41.27 | 80.34 | 100% |
| sf050_doublets_severe | 5/5 | 18.34 | 79.11 | 77.21 | 40.56 | 79.82 | 100% |
| sf050_mixed_mixed_1 | 5/5 | 11.53 | 49.48 | 56.86 | 54.45 | 56.86 | 100% |
| sf050_mixed_mixed_2 | 5/5 | 7.54 | 18.45 | 30.75 | 21.34 | 30.75 | 100% |
| sf050_mixed_mixed_3 | 5/5 | 2.47 | 2.43 | 4.84 | 4.00 | 4.84 | 100% |
| sf050_mixed_mixed_4 | 0/5 | nan | nan | nan | nan | nan | nan% |
| sf050_type_imbalance_extreme | 5/5 | 31.85 | 78.04 | 80.10 | 80.10 | 80.72 | 100% |
| sf050_type_imbalance_mild | 5/5 | 18.99 | 79.36 | 79.29 | 41.24 | 80.15 | 100% |
| sf050_type_imbalance_moderate | 5/5 | 22.15 | 75.89 | 79.21 | 54.17 | 79.39 | 100% |
| sf050_type_imbalance_severe | 5/5 | 26.38 | 71.42 | 80.01 | 67.91 | 80.01 | 100% |
| sf075_ambient_extreme | 5/5 | 21.44 | 62.59 | 63.70 | 38.99 | 66.12 | 100% |
| sf075_ambient_mild | 5/5 | 25.48 | 78.80 | 79.34 | 34.96 | 82.01 | 100% |
| sf075_ambient_moderate | 5/5 | 24.53 | 75.83 | 76.38 | 33.71 | 78.71 | 100% |
| sf075_ambient_severe | 5/5 | 23.31 | 70.70 | 71.44 | 31.92 | 73.72 | 100% |
| sf075_batch_confounding_extreme | 5/5 | 13.45 | 77.69 | 77.73 | 78.21 | 78.81 | 100% |
| sf075_batch_confounding_mild | 5/5 | 28.34 | 78.05 | 77.20 | 78.28 | 79.26 | 100% |
| sf075_batch_confounding_moderate | 5/5 | 20.47 | 74.29 | 73.71 | 74.40 | 75.25 | 100% |
| sf075_batch_confounding_severe | 5/5 | 16.21 | 77.30 | 77.91 | 78.89 | 79.09 | 100% |
| sf075_batch_depth_extreme | 5/5 | 33.18 | 84.13 | 83.20 | 52.15 | 86.66 | 100% |
| sf075_batch_depth_mild | 5/5 | 31.55 | 84.01 | 83.91 | 70.64 | 85.91 | 100% |
| sf075_batch_depth_moderate | 5/5 | 32.04 | 84.23 | 83.92 | 69.79 | 86.56 | 100% |
| sf075_batch_depth_severe | 5/5 | 32.86 | 84.47 | 84.46 | 51.14 | 86.49 | 100% |
| sf075_cells_extreme | 5/5 | 12.10 | 12.51 | 12.42 | 12.56 | 12.58 | 60% |
| sf075_cells_mild | 5/5 | 26.84 | 74.07 | 73.81 | 34.12 | 75.33 | 100% |
| sf075_cells_moderate | 5/5 | 27.63 | 64.78 | 66.09 | 27.63 | 66.46 | 100% |
| sf075_cells_severe | 5/5 | 23.97 | 24.03 | 24.03 | 23.99 | 24.13 | 40% |
| sf075_clean_baseline | 5/5 | 26.12 | 81.70 | 82.19 | 36.55 | 83.88 | 100% |
| sf075_damaged_extreme | 5/5 | 25.83 | 79.78 | 79.93 | 25.83 | 82.68 | 100% |
| sf075_damaged_mild | 5/5 | 26.01 | 80.77 | 80.96 | 26.01 | 83.03 | 100% |
| sf075_damaged_moderate | 5/5 | 26.04 | 80.87 | 80.68 | 26.04 | 82.97 | 100% |
| sf075_damaged_severe | 5/5 | 26.02 | 80.41 | 81.04 | 35.72 | 83.00 | 100% |
| sf075_depth_extreme | 0/5 | nan | nan | nan | nan | nan | nan% |
| sf075_depth_mild | 5/5 | 22.59 | 73.07 | 74.01 | 22.59 | 76.80 | 100% |
| sf075_depth_moderate | 5/5 | 18.62 | 58.40 | 59.30 | 18.62 | 61.96 | 100% |
| sf075_depth_severe | 5/5 | 13.48 | 40.45 | 40.90 | 13.48 | 42.81 | 100% |
| sf075_doublets_extreme | 5/5 | 26.57 | 80.04 | 80.38 | 35.74 | 82.18 | 100% |
| sf075_doublets_mild | 5/5 | 26.15 | 81.32 | 81.81 | 26.15 | 83.61 | 100% |
| sf075_doublets_moderate | 5/5 | 26.18 | 80.89 | 81.32 | 26.18 | 83.07 | 100% |
| sf075_doublets_severe | 5/5 | 26.23 | 80.35 | 80.68 | 26.23 | 82.52 | 100% |
| sf075_mixed_mixed_1 | 5/5 | 18.35 | 52.06 | 60.06 | 52.94 | 60.06 | 100% |
| sf075_mixed_mixed_2 | 5/5 | 15.36 | 23.21 | 33.88 | 27.34 | 33.88 | 100% |
| sf075_mixed_mixed_3 | 5/5 | 3.03 | 3.22 | 4.65 | 4.49 | 5.02 | 80% |
| sf075_mixed_mixed_4 | 0/5 | nan | nan | nan | nan | nan | nan% |
| sf075_type_imbalance_extreme | 5/5 | 46.39 | 78.46 | 82.96 | 57.56 | 84.92 | 100% |
| sf075_type_imbalance_mild | 5/5 | 28.22 | 80.28 | 81.77 | 37.00 | 82.41 | 100% |
| sf075_type_imbalance_moderate | 5/5 | 33.59 | 77.21 | 83.67 | 39.71 | 83.67 | 100% |
| sf075_type_imbalance_severe | 5/5 | 39.77 | 74.44 | 82.97 | 62.31 | 84.45 | 100% |
| sf100_ambient_extreme | 5/5 | 31.35 | 63.81 | 63.34 | 42.61 | 65.86 | 100% |
| sf100_ambient_mild | 5/5 | 36.43 | 80.07 | 80.26 | 36.43 | 81.35 | 100% |
| sf100_ambient_moderate | 5/5 | 35.30 | 77.36 | 76.66 | 41.39 | 79.19 | 100% |
| sf100_ambient_severe | 5/5 | 32.88 | 71.99 | 71.34 | 32.88 | 73.30 | 100% |
| sf100_batch_confounding_extreme | 5/5 | 12.94 | 81.87 | 81.91 | 81.94 | 83.24 | 100% |
| sf100_batch_confounding_mild | 5/5 | 32.12 | 73.73 | 72.55 | 74.27 | 74.46 | 100% |
| sf100_batch_confounding_moderate | 5/5 | 26.94 | 78.31 | 78.21 | 79.98 | 79.98 | 100% |
| sf100_batch_confounding_severe | 5/5 | 25.28 | 78.67 | 79.70 | 79.19 | 80.14 | 100% |
| sf100_batch_depth_extreme | 5/5 | 42.73 | 82.87 | 82.26 | 42.73 | 83.89 | 100% |
| sf100_batch_depth_mild | 5/5 | 40.49 | 82.32 | 82.33 | 40.49 | 82.98 | 100% |
| sf100_batch_depth_moderate | 5/5 | 41.55 | 82.91 | 82.91 | 55.04 | 84.26 | 100% |
| sf100_batch_depth_severe | 5/5 | 42.21 | 83.31 | 83.22 | 56.01 | 83.98 | 100% |
| sf100_cells_extreme | 5/5 | 21.16 | 21.61 | 21.74 | 21.74 | 21.87 | 60% |
| sf100_cells_mild | 5/5 | 33.37 | 75.23 | 74.67 | 33.37 | 76.55 | 100% |
| sf100_cells_moderate | 5/5 | 37.98 | 67.18 | 68.35 | 37.98 | 70.03 | 100% |
| sf100_cells_severe | 5/5 | 30.32 | 30.39 | 30.41 | 30.34 | 30.43 | 40% |
| sf100_clean_baseline | 5/5 | 37.83 | 81.99 | 81.88 | 37.83 | 82.98 | 100% |
| sf100_damaged_extreme | 5/5 | 37.46 | 79.36 | 78.96 | 37.46 | 80.15 | 100% |
| sf100_damaged_mild | 5/5 | 37.81 | 81.92 | 81.47 | 37.81 | 83.27 | 100% |
| sf100_damaged_moderate | 5/5 | 37.62 | 81.42 | 81.18 | 37.62 | 83.17 | 100% |
| sf100_damaged_severe | 5/5 | 37.66 | 81.54 | 81.44 | 37.66 | 82.86 | 100% |
| sf100_depth_extreme | 0/5 | nan | nan | nan | nan | nan | nan% |
| sf100_depth_mild | 5/5 | 31.93 | 73.62 | 73.01 | 31.93 | 74.69 | 100% |
| sf100_depth_moderate | 5/5 | 26.77 | 59.28 | 58.93 | 26.77 | 60.07 | 100% |
| sf100_depth_severe | 5/5 | 20.45 | 40.28 | 40.79 | 20.45 | 41.65 | 100% |
| sf100_doublets_extreme | 5/5 | 38.27 | 80.51 | 79.97 | 38.27 | 81.23 | 100% |
| sf100_doublets_mild | 5/5 | 37.88 | 81.59 | 81.41 | 37.88 | 82.46 | 100% |
| sf100_doublets_moderate | 5/5 | 37.90 | 81.15 | 80.46 | 37.90 | 82.22 | 100% |
| sf100_doublets_severe | 5/5 | 38.15 | 80.86 | 80.53 | 38.15 | 81.77 | 100% |
| sf100_mixed_mixed_1 | 5/5 | 27.81 | 52.10 | 61.54 | 53.65 | 61.54 | 100% |
| sf100_mixed_mixed_2 | 5/5 | 14.84 | 25.45 | 31.79 | 26.55 | 31.79 | 100% |
| sf100_mixed_mixed_3 | 5/5 | 3.81 | 4.16 | 4.40 | 4.82 | 5.10 | 40% |
| sf100_mixed_mixed_4 | 0/5 | nan | nan | nan | nan | nan | nan% |
| sf100_type_imbalance_extreme | 5/5 | 67.63 | 80.00 | 85.30 | 77.88 | 85.30 | 100% |
| sf100_type_imbalance_mild | 5/5 | 41.18 | 81.37 | 85.03 | 41.18 | 85.03 | 100% |
| sf100_type_imbalance_moderate | 5/5 | 50.14 | 78.54 | 84.77 | 60.50 | 84.77 | 100% |
| sf100_type_imbalance_severe | 5/5 | 59.22 | 75.46 | 84.96 | 61.93 | 84.96 | 100% |

## Detection accuracy (share of datasets where the engine got it right)

| check | n | accuracy | scenarios where it failed |
|---|---|---|---|
| batch | 40 | 100% | - |
| batch_detected | 740 | 100% | - |
| context_scope | 660 | 99% | sf025_mixed_mixed_2, sf050_mixed_mixed_2, sf075_mixed_mixed_2 |
| doublets_removed | 620 | 100% | - |
| organism | 740 | 100% | - |
| technology | 740 | 100% | - |

## Cleaning accuracy (mean over datasets)

| measure | baseline | auto |
|---|---|---|
| bad_removed | 52.3% | 53.0% |
| ok_retained | 93.4% | 97.3% |
| removed_damaged | 100.0% | 99.9% |
| removed_doublet | 5.5% | 6.1% |

## Requested inference models

The same held-out context truths are used for all arms. The common `auto` pipeline is chosen by the fast proxy methods on the selection half; expensive models never see the evaluation half during pipeline selection.

### AUPRC ratio (1 = random)

| model | n auto | mean baseline | mean rules | mean auto | auto > baseline |
|---|---|---|---|---|---|
| clr | 740 | 17.662 | 53.231 | 54.695 | 668/700 |
| grnboost2 | 740 | 54.222 | 65.719 | 66.411 | 548/700 |
| pearson_signed | 740 | 24.237 | 70.409 | 70.793 | 674/700 |
| regdiffusion | 740 | 2.065 | 3.790 | 3.714 | 642/700 |

### EPR (1 = random)

| model | n auto | mean baseline | mean rules | mean auto | auto > baseline |
|---|---|---|---|---|---|
| clr | 740 | 24.500 | 64.489 | 65.340 | 657/700 |
| grnboost2 | 740 | 63.286 | 72.087 | 72.231 | 562/700 |
| pearson_signed | 740 | 26.322 | 73.457 | 73.721 | 664/700 |
| regdiffusion | 740 | 5.333 | 9.701 | 9.467 | 621/700 |


## Controlled stress-test breaking points

Exploratory diagnostics: `statistical` is the first level whose 95% t interval no longer lies wholly above random (AUPRC ratio 1); `half-clean` is the first level at or below 50% of the clean mean; `operational` is the first level with any missing model score.

| shared edges | model | stressor | statistical | half-clean | operational |
|---|---|---|---|---|---|
| 25% | pearson_signed | ambient | not reached | not reached | not reached |
| 25% | pearson_signed | batch_confounding | not reached | not reached | not reached |
| 25% | pearson_signed | batch_depth | not reached | not reached | not reached |
| 25% | pearson_signed | cells | not reached | severe | not reached |
| 25% | pearson_signed | damaged | not reached | not reached | not reached |
| 25% | pearson_signed | depth | not reached | extreme | not reached |
| 25% | pearson_signed | doublets | not reached | not reached | not reached |
| 25% | pearson_signed | mixed | mixed_4 | mixed_3 | not reached |
| 25% | pearson_signed | type_imbalance | not reached | not reached | not reached |
| 25% | clr | ambient | not reached | extreme | not reached |
| 25% | clr | batch_confounding | not reached | not reached | not reached |
| 25% | clr | batch_depth | not reached | not reached | not reached |
| 25% | clr | cells | not reached | severe | not reached |
| 25% | clr | damaged | not reached | not reached | not reached |
| 25% | clr | depth | not reached | moderate | not reached |
| 25% | clr | doublets | not reached | not reached | not reached |
| 25% | clr | mixed | mixed_4 | mixed_1 | not reached |
| 25% | clr | type_imbalance | not reached | not reached | not reached |
| 25% | grnboost2 | ambient | not reached | not reached | not reached |
| 25% | grnboost2 | batch_confounding | not reached | not reached | not reached |
| 25% | grnboost2 | batch_depth | not reached | not reached | not reached |
| 25% | grnboost2 | cells | not reached | severe | not reached |
| 25% | grnboost2 | damaged | not reached | not reached | not reached |
| 25% | grnboost2 | depth | not reached | severe | not reached |
| 25% | grnboost2 | doublets | not reached | not reached | not reached |
| 25% | grnboost2 | mixed | not reached | mixed_2 | not reached |
| 25% | grnboost2 | type_imbalance | not reached | not reached | not reached |
| 25% | regdiffusion | ambient | not reached | extreme | not reached |
| 25% | regdiffusion | batch_confounding | not reached | not reached | not reached |
| 25% | regdiffusion | batch_depth | not reached | not reached | not reached |
| 25% | regdiffusion | cells | not reached | severe | not reached |
| 25% | regdiffusion | damaged | not reached | not reached | not reached |
| 25% | regdiffusion | depth | not reached | moderate | not reached |
| 25% | regdiffusion | doublets | not reached | not reached | not reached |
| 25% | regdiffusion | mixed | mixed_3 | mixed_1 | not reached |
| 25% | regdiffusion | type_imbalance | not reached | not reached | not reached |
| 50% | pearson_signed | ambient | not reached | not reached | not reached |
| 50% | pearson_signed | batch_confounding | not reached | not reached | not reached |
| 50% | pearson_signed | batch_depth | not reached | not reached | not reached |
| 50% | pearson_signed | cells | not reached | severe | not reached |
| 50% | pearson_signed | damaged | not reached | not reached | not reached |
| 50% | pearson_signed | depth | not reached | extreme | not reached |
| 50% | pearson_signed | doublets | not reached | not reached | not reached |
| 50% | pearson_signed | mixed | not reached | mixed_3 | not reached |
| 50% | pearson_signed | type_imbalance | not reached | not reached | not reached |
| 50% | clr | ambient | not reached | extreme | not reached |
| 50% | clr | batch_confounding | not reached | not reached | not reached |
| 50% | clr | batch_depth | not reached | not reached | not reached |
| 50% | clr | cells | not reached | severe | not reached |
| 50% | clr | damaged | not reached | not reached | not reached |
| 50% | clr | depth | not reached | moderate | not reached |
| 50% | clr | doublets | not reached | not reached | not reached |
| 50% | clr | mixed | mixed_3 | mixed_1 | not reached |
| 50% | clr | type_imbalance | not reached | not reached | not reached |
| 50% | grnboost2 | ambient | not reached | not reached | not reached |
| 50% | grnboost2 | batch_confounding | not reached | not reached | not reached |
| 50% | grnboost2 | batch_depth | not reached | not reached | not reached |
| 50% | grnboost2 | cells | not reached | severe | not reached |
| 50% | grnboost2 | damaged | not reached | not reached | not reached |
| 50% | grnboost2 | depth | not reached | severe | not reached |
| 50% | grnboost2 | doublets | not reached | not reached | not reached |
| 50% | grnboost2 | mixed | not reached | mixed_2 | not reached |
| 50% | grnboost2 | type_imbalance | not reached | not reached | not reached |
| 50% | regdiffusion | ambient | not reached | extreme | not reached |
| 50% | regdiffusion | batch_confounding | not reached | not reached | not reached |
| 50% | regdiffusion | batch_depth | not reached | not reached | not reached |
| 50% | regdiffusion | cells | not reached | severe | not reached |
| 50% | regdiffusion | damaged | not reached | not reached | not reached |
| 50% | regdiffusion | depth | not reached | moderate | not reached |
| 50% | regdiffusion | doublets | not reached | not reached | not reached |
| 50% | regdiffusion | mixed | mixed_2 | mixed_2 | not reached |
| 50% | regdiffusion | type_imbalance | not reached | not reached | not reached |
| 75% | pearson_signed | ambient | not reached | not reached | not reached |
| 75% | pearson_signed | batch_confounding | not reached | not reached | not reached |
| 75% | pearson_signed | batch_depth | not reached | not reached | not reached |
| 75% | pearson_signed | cells | not reached | severe | not reached |
| 75% | pearson_signed | damaged | not reached | not reached | not reached |
| 75% | pearson_signed | depth | not reached | extreme | not reached |
| 75% | pearson_signed | doublets | not reached | not reached | not reached |
| 75% | pearson_signed | mixed | not reached | mixed_3 | not reached |
| 75% | pearson_signed | type_imbalance | not reached | not reached | not reached |
| 75% | clr | ambient | not reached | not reached | not reached |
| 75% | clr | batch_confounding | not reached | not reached | not reached |
| 75% | clr | batch_depth | not reached | not reached | not reached |
| 75% | clr | cells | not reached | severe | not reached |
| 75% | clr | damaged | not reached | not reached | not reached |
| 75% | clr | depth | not reached | moderate | not reached |
| 75% | clr | doublets | not reached | not reached | not reached |
| 75% | clr | mixed | mixed_4 | mixed_1 | not reached |
| 75% | clr | type_imbalance | not reached | not reached | not reached |
| 75% | grnboost2 | ambient | not reached | not reached | not reached |
| 75% | grnboost2 | batch_confounding | not reached | not reached | not reached |
| 75% | grnboost2 | batch_depth | not reached | not reached | not reached |
| 75% | grnboost2 | cells | not reached | extreme | not reached |
| 75% | grnboost2 | damaged | not reached | not reached | not reached |
| 75% | grnboost2 | depth | not reached | severe | not reached |
| 75% | grnboost2 | doublets | not reached | not reached | not reached |
| 75% | grnboost2 | mixed | not reached | mixed_2 | not reached |
| 75% | grnboost2 | type_imbalance | not reached | not reached | not reached |
| 75% | regdiffusion | ambient | not reached | extreme | not reached |
| 75% | regdiffusion | batch_confounding | not reached | not reached | not reached |
| 75% | regdiffusion | batch_depth | not reached | not reached | not reached |
| 75% | regdiffusion | cells | not reached | extreme | not reached |
| 75% | regdiffusion | damaged | not reached | not reached | not reached |
| 75% | regdiffusion | depth | extreme | moderate | not reached |
| 75% | regdiffusion | doublets | not reached | not reached | not reached |
| 75% | regdiffusion | mixed | mixed_4 | mixed_2 | not reached |
| 75% | regdiffusion | type_imbalance | not reached | not reached | not reached |
| 100% | pearson_signed | ambient | not reached | not reached | not reached |
| 100% | pearson_signed | batch_confounding | not reached | not reached | not reached |
| 100% | pearson_signed | batch_depth | not reached | not reached | not reached |
| 100% | pearson_signed | cells | not reached | severe | not reached |
| 100% | pearson_signed | damaged | not reached | not reached | not reached |
| 100% | pearson_signed | depth | not reached | extreme | not reached |
| 100% | pearson_signed | doublets | not reached | not reached | not reached |
| 100% | pearson_signed | mixed | not reached | mixed_2 | not reached |
| 100% | pearson_signed | type_imbalance | not reached | not reached | not reached |
| 100% | clr | ambient | not reached | extreme | not reached |
| 100% | clr | batch_confounding | not reached | not reached | not reached |
| 100% | clr | batch_depth | not reached | not reached | not reached |
| 100% | clr | cells | not reached | severe | not reached |
| 100% | clr | damaged | not reached | not reached | not reached |
| 100% | clr | depth | not reached | moderate | not reached |
| 100% | clr | doublets | not reached | not reached | not reached |
| 100% | clr | mixed | mixed_4 | mixed_1 | not reached |
| 100% | clr | type_imbalance | not reached | not reached | not reached |
| 100% | grnboost2 | ambient | not reached | not reached | not reached |
| 100% | grnboost2 | batch_confounding | not reached | not reached | not reached |
| 100% | grnboost2 | batch_depth | not reached | not reached | not reached |
| 100% | grnboost2 | cells | not reached | extreme | not reached |
| 100% | grnboost2 | damaged | not reached | not reached | not reached |
| 100% | grnboost2 | depth | not reached | severe | not reached |
| 100% | grnboost2 | doublets | not reached | not reached | not reached |
| 100% | grnboost2 | mixed | not reached | mixed_2 | not reached |
| 100% | grnboost2 | type_imbalance | not reached | not reached | not reached |
| 100% | regdiffusion | ambient | not reached | extreme | not reached |
| 100% | regdiffusion | batch_confounding | not reached | not reached | not reached |
| 100% | regdiffusion | batch_depth | not reached | not reached | not reached |
| 100% | regdiffusion | cells | extreme | not reached | not reached |
| 100% | regdiffusion | damaged | not reached | not reached | not reached |
| 100% | regdiffusion | depth | moderate | moderate | not reached |
| 100% | regdiffusion | doublets | not reached | not reached | not reached |
| 100% | regdiffusion | mixed | mixed_3 | mixed_2 | not reached |
| 100% | regdiffusion | type_imbalance | not reached | not reached | not reached |


Median runtime per completed dataset: 670 s.
