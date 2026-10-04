# Candidate study - SUMMARY

EXPLORATORY. Descriptive means only; no p-values, no claims. Unit = task (dataset x condition); eval = AUPRC ratio on the CollecTRI evaluation half, mean of the evaluation methods (1 = random).

## Per condition

| condition | ok | skipped | error | baseline | engine_auto | engine_auto_noref | rules | oracle engine+baseline | oracle all | noref worse than baseline |
|---|---|---|---|---|---|---|---|---|---|---|
| clean | 1 | 0 | 0 | 1.323 | 1.342 | 1.323 | 1.343 | 1.343 | 1.434 | 0/1 |

## Variants: mean (eval_mean - baseline eval_mean), n tasks in brackets

| variant | clean |
|---|---|
| base+batch=center_per_batch | - |
| base+batch=combat | - |
| base+cell_qc_thresholds=mad | +0.027 [1] |
| base+cell_qc_thresholds=mad_per_batch | +0.027 [1] |
| base+context_scope=per_cell_type | +0.000 [1] |
| base+context_scope=selected_type | +0.000 [1] |
| base+doublets=scrublet | -0.001 [1] |
| base+gene_filter=min_cells_drop_qc_genes | -0.004 [1] |
| base+hvg=tf_plus_2000 | +0.111 [1] |
| base+hvg=tf_plus_500 | -0.061 [1] |
| base+normalization=pearson_residuals | -0.295 [1] |
| base+sparsity=metacells | -0.070 [1] |

## Variants vs baseline on shared genes: mean (pair_eval_mean - pair_base_eval_mean), both on TF x the genes the variant and the baseline share, n tasks in brackets

| variant | clean |
|---|---|
| base+batch=center_per_batch | - |
| base+batch=combat | - |
| base+cell_qc_thresholds=mad | +0.035 [1] |
| base+cell_qc_thresholds=mad_per_batch | +0.035 [1] |
| base+context_scope=per_cell_type | +0.000 [1] |
| base+context_scope=selected_type | +0.000 [1] |
| base+doublets=scrublet | -0.001 [1] |
| base+gene_filter=min_cells_drop_qc_genes | -0.008 [1] |
| base+hvg=tf_plus_2000 | +0.009 [1] |
| base+hvg=tf_plus_500 | +0.000 [1] |
| base+normalization=pearson_residuals | -0.295 [1] |
| base+sparsity=metacells | -0.074 [1] |

## Selection rules: mean (chosen - baseline) eval, wins/losses vs baseline

Tasks where the baseline itself failed (QC gate refusal or error) are left out of the mean; `+k` counts those where the rule still produced a scored network.

Baseline failed: clean 0/1

| rule | clean |
|---|---|
| with reference: selection half, engine pool | +0.018 (1/0) |
| with reference: selection half, all candidates | +0.111 (1/0) |
| with reference: selection half, all candidates, judged on shared genes | +0.009 (1/0) |
| no reference: shipped (stability - depth - batch), engine pool | +0.000 (0/0) |
| no reference: stability - depth leakage, engine pool | +0.000 (0/0) |
| no reference: stability - depth leakage, all candidates | -0.070 (0/1) |
| no reference: stability only, all candidates | -0.070 (0/1) |

Variants skipped (option unavailable): batch=harmony, normalization=scran
