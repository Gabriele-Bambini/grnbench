# Candidate study - SUMMARY

EXPLORATORY. Descriptive means only; no p-values, no claims. Unit = task (dataset x condition); eval = AUPRC ratio on the CollecTRI evaluation half, mean of the evaluation methods (1 = random).

## Per condition

| condition | ok | skipped | error | baseline | engine_auto | engine_auto_noref | rules | oracle engine+baseline | oracle all | noref worse than baseline |
|---|---|---|---|---|---|---|---|---|---|---|
| clean | 16 | 0 | 0 | 1.377 | 1.382 | 1.347 | 1.353 | 1.399 | 1.457 | 5/15 |

## Variants: mean (eval_mean - baseline eval_mean), n tasks in brackets

| variant | clean |
|---|---|
| base+batch=center_per_batch | - |
| base+batch=combat | - |
| base+cell_qc_thresholds=mad | -0.004 [15] |
| base+cell_qc_thresholds=mad_per_batch | -0.004 [15] |
| base+context_scope=per_cell_type | +0.000 [15] |
| base+context_scope=selected_type | +0.000 [15] |
| base+doublets=scrublet | +0.001 [15] |
| base+gene_filter=min_cells_drop_qc_genes | +0.024 [15] |
| base+hvg=tf_plus_2000 | +0.049 [15] |
| base+hvg=tf_plus_500 | -0.003 [15] |
| base+normalization=pearson_residuals | -0.306 [15] |
| base+sparsity=metacells | -0.027 [15] |

## Variants vs baseline on shared genes: mean (pair_eval_mean - pair_base_eval_mean), both on TF x the genes the variant and the baseline share, n tasks in brackets

| variant | clean |
|---|---|
| base+batch=center_per_batch | - |
| base+batch=combat | - |
| base+cell_qc_thresholds=mad | +0.004 [15] |
| base+cell_qc_thresholds=mad_per_batch | +0.004 [15] |
| base+context_scope=per_cell_type | +0.000 [15] |
| base+context_scope=selected_type | +0.000 [15] |
| base+doublets=scrublet | -0.002 [15] |
| base+gene_filter=min_cells_drop_qc_genes | -0.006 [15] |
| base+hvg=tf_plus_2000 | +0.003 [15] |
| base+hvg=tf_plus_500 | +0.002 [15] |
| base+normalization=pearson_residuals | -0.306 [15] |
| base+sparsity=metacells | -0.016 [15] |

## Selection rules: mean (chosen - baseline) eval, wins/losses vs baseline

Tasks where the baseline itself failed (QC gate refusal or error) are left out of the mean; `+k` counts those where the rule still produced a scored network.

Baseline failed: clean 1/16

| rule | clean |
|---|---|
| with reference: selection half, engine pool | +0.022 (5/1) +1 |
| with reference: selection half, all candidates | +0.077 (14/0) +1 |
| with reference: selection half, all candidates, judged on shared genes | +0.025 (11/3) |
| no reference: shipped (stability - depth - batch), engine pool | -0.015 (3/5) +1 |
| no reference: stability - depth leakage, engine pool | -0.015 (3/5) +1 |
| no reference: stability - depth leakage, all candidates | -0.021 (6/8) +1 |
| no reference: stability only, all candidates | -0.068 (3/11) +1 |

Variants skipped (option unavailable): batch=harmony, normalization=scran

## Raw droplet matrix vs vendor-filtered matrix (per dataset)

| dataset | cells called from raw | cells in filtered | engine pick (raw) | engine (raw) | baseline (raw) | baseline (filtered) | engine (raw) - baseline (filtered) |
|---|---|---|---|---|---|---|---|
| 10x__5k_pbmc_protein_v3 | 5215 | 5247 | baseline | 1.510 | 1.510 | 1.496 | +0.014 |
| 10x__heart_10k_v3 | 92106 | 7713 | recommended | 1.110 | refused/failed | 1.396 | -0.286 |
| 10x__heart_1k_v3 | 983 | 1011 | gene_filter=min_cells_drop_qc_genes | 1.444 | 1.351 | 1.320 | +0.124 |
| 10x__neuron_10k_v3 | 11732 | 11843 | baseline | 1.216 | 1.216 | 1.238 | -0.022 |
| 10x__neuron_1k_v3 | 1232 | 1301 | baseline | 1.256 | 1.256 | 1.244 | +0.012 |
| 10x__pbmc_10k_v3 | 11811 | 11769 | baseline | 1.563 | 1.563 | 1.591 | -0.028 |
| 10x__pbmc_1k_v2 | 998 | 996 | recommended | 1.361 | 1.362 | 1.358 | +0.003 |
| 10x__pbmc_1k_v3 | 1224 | 1222 | baseline | 1.385 | 1.385 | 1.374 | +0.010 |
