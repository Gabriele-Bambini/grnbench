# Candidate study - SUMMARY

EXPLORATORY. Descriptive means only; no p-values, no claims. Unit = task (dataset x condition); eval = AUPRC ratio on the CollecTRI evaluation half, mean of the evaluation methods (1 = random).

## Per condition

| condition | ok | skipped | error | baseline | engine_auto | engine_auto_noref | rules | oracle engine+baseline | oracle all | noref worse than baseline |
|---|---|---|---|---|---|---|---|---|---|---|
| clean | 20 | 0 | 0 | 1.351 | 1.357 | 1.227 | 1.253 | 1.361 | 1.407 | 11/20 |

## Variants: mean (eval_mean - baseline eval_mean), n tasks in brackets

| variant | clean |
|---|---|
| base+batch=center_per_batch | -0.255 [19] |
| base+batch=combat | -0.249 [19] |
| base+cell_qc_thresholds=mad | +0.003 [20] |
| base+cell_qc_thresholds=mad_per_batch | -0.002 [20] |
| base+context_scope=per_cell_type | -0.075 [20] |
| base+context_scope=selected_type | -0.066 [20] |
| base+doublets=scrublet | +0.002 [19] |
| base+gene_filter=min_cells_drop_qc_genes | +0.010 [20] |
| base+hvg=tf_plus_2000 | +0.037 [20] |
| base+hvg=tf_plus_500 | -0.041 [20] |
| base+normalization=pearson_residuals | -0.308 [20] |
| base+sparsity=metacells | -0.159 [20] |

## Variants vs baseline on shared genes: mean (pair_eval_mean - pair_base_eval_mean), both on TF x the genes the variant and the baseline share, n tasks in brackets

| variant | clean |
|---|---|
| base+batch=center_per_batch | -0.255 [19] |
| base+batch=combat | -0.249 [19] |
| base+cell_qc_thresholds=mad | -0.003 [20] |
| base+cell_qc_thresholds=mad_per_batch | -0.003 [20] |
| base+context_scope=per_cell_type | -0.118 [20] |
| base+context_scope=selected_type | -0.107 [20] |
| base+doublets=scrublet | -0.004 [19] |
| base+gene_filter=min_cells_drop_qc_genes | +0.002 [20] |
| base+hvg=tf_plus_2000 | +0.008 [20] |
| base+hvg=tf_plus_500 | -0.008 [20] |
| base+normalization=pearson_residuals | -0.308 [20] |
| base+sparsity=metacells | -0.078 [20] |

## Selection rules: mean (chosen - baseline) eval, wins/losses vs baseline

Tasks where the baseline itself failed (QC gate refusal or error) are left out of the mean; `+k` counts those where the rule still produced a scored network.

Baseline failed: clean 0/20

| rule | clean |
|---|---|
| with reference: selection half, engine pool | +0.006 (5/2) |
| with reference: selection half, all candidates | +0.044 (16/2) |
| with reference: selection half, all candidates, judged on shared genes | +0.003 (14/4) |
| no reference: shipped (stability - depth - batch), engine pool | -0.124 (3/11) |
| no reference: stability - depth leakage, engine pool | -0.095 (4/10) |
| no reference: stability - depth leakage, all candidates | -0.209 (1/19) |
| no reference: stability only, all candidates | -0.224 (1/19) |

Variants skipped (option unavailable): batch=harmony, normalization=scran
