# Candidate study - SUMMARY

EXPLORATORY. Descriptive means only; no p-values, no claims. Unit = task (dataset x condition); eval = AUPRC ratio on the CollecTRI evaluation half, mean of the evaluation methods (1 = random).

## Per condition

| condition | ok | skipped | error | baseline | engine_auto | engine_auto_noref | rules | oracle engine+baseline | oracle all | noref worse than baseline |
|---|---|---|---|---|---|---|---|---|---|---|
| clean | 20 | 0 | 0 | 1.351 | 1.357 | 1.227 | 1.253 | 1.361 | 1.407 | 11/20 |
| depth_30 | 20 | 0 | 0 | 1.371 | 1.376 | 1.206 | 1.266 | 1.384 | 1.436 | 11/19 |
| depth_10 | 20 | 0 | 0 | 1.339 | 1.388 | 1.277 | 1.295 | 1.398 | 1.423 | 5/13 |
| ambient_20 | 20 | 0 | 0 | 1.336 | 1.336 | 1.181 | 1.248 | 1.343 | 1.379 | 14/20 |
| doublets_20 | 20 | 0 | 0 | 1.361 | 1.371 | 1.206 | 1.282 | 1.376 | 1.424 | 13/20 |
| damaged_15 | 20 | 0 | 0 | 1.388 | 1.390 | 1.201 | 1.261 | 1.394 | 1.442 | 15/20 |
| mixed | 19 | 0 | 1 | 1.327 | 1.324 | 1.105 | 1.243 | 1.332 | 1.368 | 14/18 |

## Variants: mean (eval_mean - baseline eval_mean), n tasks in brackets

| variant | clean | depth_30 | depth_10 | ambient_20 | doublets_20 | damaged_15 | mixed |
|---|---|---|---|---|---|---|---|
| base+batch=center_per_batch | -0.255 [19] | -0.280 [18] | -0.240 [13] | -0.213 [19] | -0.267 [19] | -0.277 [19] | -0.220 [17] |
| base+batch=combat | -0.239 [17] | -0.278 [15] | -0.236 [11] | -0.202 [17] | -0.261 [17] | -0.263 [17] | -0.210 [15] |
| base+cell_qc_thresholds=mad | +0.003 [20] | +0.017 [19] | -0.001 [13] | -0.010 [20] | -0.002 [20] | -0.016 [20] | +0.004 [18] |
| base+cell_qc_thresholds=mad_per_batch | -0.002 [20] | +0.001 [19] | +0.001 [13] | -0.023 [20] | -0.005 [20] | -0.021 [20] | -0.002 [18] |
| base+context_scope=per_cell_type | -0.075 [20] | -0.082 [19] | -0.086 [13] | -0.097 [20] | -0.063 [20] | -0.096 [20] | -0.083 [18] |
| base+context_scope=selected_type | -0.066 [20] | -0.058 [19] | -0.076 [13] | -0.099 [20] | -0.047 [20] | -0.085 [20] | -0.084 [18] |
| base+doublets=scrublet | +0.002 [19] | -0.001 [19] | +0.005 [13] | -0.008 [18] | +0.000 [20] | +0.003 [20] | -0.003 [18] |
| base+gene_filter=min_cells_drop_qc_genes | +0.010 [20] | -0.001 [19] | +0.006 [13] | -0.001 [20] | +0.004 [20] | +0.004 [20] | -0.003 [18] |
| base+hvg=tf_plus_2000 | +0.037 [20] | +0.020 [19] | +0.010 [13] | +0.018 [20] | +0.033 [20] | +0.015 [20] | +0.010 [18] |
| base+hvg=tf_plus_500 | -0.041 [20] | -0.045 [19] | -0.027 [13] | -0.055 [20] | -0.060 [20] | -0.048 [20] | -0.048 [18] |
| base+normalization=pearson_residuals | -0.308 [20] | -0.318 [19] | -0.259 [13] | -0.280 [20] | -0.319 [20] | -0.334 [20] | -0.288 [18] |
| base+sparsity=metacells | -0.159 [20] | -0.126 [19] | -0.126 [13] | -0.143 [20] | -0.160 [20] | -0.190 [20] | -0.118 [18] |

## Selection rules: mean (chosen - baseline) eval, wins/losses vs baseline

Tasks where the baseline itself failed (QC gate refusal or error) are left out of the mean; `+k` counts those where the rule still produced a scored network.

Baseline failed: clean 0/20, depth_30 1/20, depth_10 7/20, ambient_20 0/20, doublets_20 0/20, damaged_15 0/20, mixed 1/19

| rule | clean | depth_30 | depth_10 | ambient_20 | doublets_20 | damaged_15 | mixed |
|---|---|---|---|---|---|---|---|
| with reference: selection half, engine pool | +0.006 (5/2) | +0.021 (2/1) +1 | +0.008 (2/1) +3 | -0.001 (2/2) | +0.010 (5/3) | +0.002 (3/3) | +0.011 (2/3) +1 |
| with reference: selection half, all candidates | +0.044 (16/2) | +0.069 (17/1) +1 | +0.023 (8/3) +3 | +0.029 (16/3) | +0.041 (17/3) | +0.039 (16/4) | +0.040 (13/5) +1 |
| no reference: shipped (stability - depth - batch), engine pool | -0.124 (3/11) | -0.159 (1/11) +1 | -0.085 (2/5) +3 | -0.156 (1/14) | -0.155 (0/13) | -0.187 (0/15) | -0.220 (1/14) +1 |
| no reference: stability - depth leakage, engine pool | -0.095 (4/10) | -0.138 (2/10) +1 | -0.088 (1/6) +3 | -0.131 (1/14) | -0.143 (0/12) | -0.153 (1/14) | -0.220 (1/14) +1 |
| no reference: stability - depth leakage, all candidates | -0.209 (1/19) | -0.157 (2/17) +1 | -0.114 (4/9) +3 | -0.191 (3/17) | -0.239 (0/20) | -0.224 (1/19) | -0.235 (0/18) +1 |
| no reference: stability only, all candidates | -0.224 (1/19) | -0.249 (1/18) +1 | -0.221 (1/12) +3 | -0.202 (2/18) | -0.244 (0/20) | -0.226 (1/19) | -0.265 (1/17) +1 |

Variants skipped (option unavailable): batch=harmony, normalization=scran

## Task errors

- cxg__10x-3-v3__homo-sapiens__prostate-gland__20de8570__mixed: ValueError: Input X contains NaN.
PCA does not accept missing values encoded as NaN natively. For supervised learning, you might want to consider sklearn.ensemble.HistGradientBoostingClassifier and Re
