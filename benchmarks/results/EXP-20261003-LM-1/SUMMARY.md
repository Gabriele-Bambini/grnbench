# Candidate study - SUMMARY

EXPLORATORY. Descriptive means only; no p-values, no claims. Unit = task (dataset x condition); eval = AUPRC ratio on the CollecTRI evaluation half, mean of the evaluation methods (1 = random).

## Per condition

| condition | ok | skipped | error | baseline | engine_auto | engine_auto_noref | rules | oracle engine+baseline | oracle all | noref worse than baseline |
|---|---|---|---|---|---|---|---|---|---|---|
| clean | 20 | 0 | 0 | 1.355 | 1.371 | 1.266 | 1.248 | 1.384 | 1.430 | 13/20 |
| doublets_20 | 20 | 0 | 0 | 1.367 | 1.368 | 1.255 | 1.245 | 1.379 | 1.426 | 13/20 |
| damaged_15 | 20 | 0 | 0 | 1.398 | 1.415 | 1.272 | 1.225 | 1.416 | 1.453 | 14/20 |

## Variants: mean (eval_mean - baseline eval_mean), n tasks in brackets

| variant | clean | doublets_20 | damaged_15 |
|---|---|---|---|
| base+batch=center_per_batch | -0.245 [15] | -0.264 [15] | -0.276 [15] |
| base+batch=combat | -0.251 [15] | -0.265 [15] | -0.285 [15] |
| base+cell_qc_thresholds=mad | -0.008 [20] | -0.013 [20] | -0.038 [20] |
| base+cell_qc_thresholds=mad_per_batch | -0.009 [20] | -0.019 [20] | -0.040 [20] |
| base+context_scope=per_cell_type | -0.108 [20] | -0.100 [20] | -0.128 [20] |
| base+context_scope=selected_type | -0.093 [20] | -0.082 [20] | -0.118 [20] |
| base+doublets=scrublet | -0.001 [18] | -0.004 [19] | -0.000 [19] |
| base+gene_filter=min_cells_drop_qc_genes | +0.004 [20] | +0.008 [20] | +0.009 [20] |
| base+hvg=tf_plus_2000 | +0.031 [20] | +0.041 [20] | +0.016 [20] |
| base+hvg=tf_plus_500 | -0.063 [20] | -0.056 [20] | -0.076 [20] |
| base+normalization=pearson_residuals | -0.305 [20] | -0.319 [20] | -0.335 [20] |
| base+sparsity=metacells | -0.142 [20] | -0.157 [20] | -0.183 [20] |

## Selection rules: mean (chosen - baseline) eval, wins/losses vs baseline

Tasks where the baseline itself failed (QC gate refusal or error) are left out of the mean; `+k` counts those where the rule still produced a scored network.

Baseline failed: clean 0/20, doublets_20 0/20, damaged_15 0/20

| rule | clean | doublets_20 | damaged_15 |
|---|---|---|---|
| with reference: selection half, engine pool | +0.016 (4/3) | +0.001 (4/3) | +0.017 (3/1) |
| with reference: selection half, all candidates | +0.060 (16/4) | +0.057 (17/2) | +0.044 (16/4) |
| no reference: shipped (stability - depth - batch), engine pool | -0.089 (2/13) | -0.112 (1/13) | -0.126 (2/14) |
| no reference: stability - depth leakage, engine pool | -0.056 (2/11) | -0.112 (1/13) | -0.126 (2/14) |
| no reference: stability - depth leakage, all candidates | -0.155 (2/18) | -0.149 (3/17) | -0.190 (3/17) |
| no reference: stability only, all candidates | -0.225 (2/18) | -0.186 (4/16) | -0.207 (2/18) |

Variants skipped (option unavailable): batch=harmony, normalization=scran
