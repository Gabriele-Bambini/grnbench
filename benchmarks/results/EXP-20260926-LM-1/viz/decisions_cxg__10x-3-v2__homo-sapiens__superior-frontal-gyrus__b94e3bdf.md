| decision | value | reason |
|---|---|---|
| identifier_mapping | keep | rules recommendation |
| cell_calling | as_provided | rules recommendation |
| ambient_rna | none | rules recommendation |
| cell_qc_thresholds | mad | rules recommendation |
| doublets | scrublet | rules recommendation |
| gene_filter | min_cells | rules recommendation |
| context_scope | per_cell_type | rules recommendation |
| normalization | cpm_log1p | rules recommendation |
| batch | none | rules recommendation |
| confounders | none | rules recommendation |
| sparsity | metacells | rules recommendation |
| hvg | tf_plus_1000 | rules recommendation |
| regulators | user_file | rules recommendation |
| subsample | none | rules recommendation |
| pseudotime | none | rules recommendation |
| candidate (baseline) | baseline | chosen for the baseline arm |
| warning (baseline) | - | identifier_mapping: best option 'Map Ensembl IDs to symbols' is unavailable (planned in grnbench, not implemented yet); using 'keep'. |
| warning (baseline) | - | ambient_rna: best option 'DecontX' is unavailable (missing: r:celda); using 'none'. |
| warning (baseline) | - | doublets: best option 'scDblFinder' is unavailable (missing: r:scDblFinder); using 'scrublet'. |
| candidate (rules) | recommended | chosen for the rules arm |
| warning (rules) | - | metacells: the smallest context has 1016 cells -> ~50 metacells of 20; fewer than ~100 samples makes gene-gene statistics unstable (Bilous et al. 2024). Use single cells, pooled context or smaller metacells. |
| warning (rules) | - | identifier_mapping: best option 'Map Ensembl IDs to symbols' is unavailable (planned in grnbench, not implemented yet); using 'keep'. |
| warning (rules) | - | ambient_rna: best option 'DecontX' is unavailable (missing: r:celda); using 'none'. |
| warning (rules) | - | doublets: best option 'scDblFinder' is unavailable (missing: r:scDblFinder); using 'scrublet'. |
| candidate (auto) | baseline | chosen for the auto arm |
| warning (auto) | - | identifier_mapping: best option 'Map Ensembl IDs to symbols' is unavailable (planned in grnbench, not implemented yet); using 'keep'. |
| warning (auto) | - | ambient_rna: best option 'DecontX' is unavailable (missing: r:celda); using 'none'. |
| warning (auto) | - | doublets: best option 'scDblFinder' is unavailable (missing: r:scDblFinder); using 'scrublet'. |
| candidate (auto_noref) | context_scope=selected_type | chosen for the auto_noref arm |
| warning (auto_noref) | - | identifier_mapping: best option 'Map Ensembl IDs to symbols' is unavailable (planned in grnbench, not implemented yet); using 'keep'. |
| warning (auto_noref) | - | ambient_rna: best option 'DecontX' is unavailable (missing: r:celda); using 'none'. |
| warning (auto_noref) | - | doublets: best option 'scDblFinder' is unavailable (missing: r:scDblFinder); using 'scrublet'. |
| candidate (oracle) | baseline | chosen for the oracle arm |
| warning (oracle) | - | identifier_mapping: best option 'Map Ensembl IDs to symbols' is unavailable (planned in grnbench, not implemented yet); using 'keep'. |
| warning (oracle) | - | ambient_rna: best option 'DecontX' is unavailable (missing: r:celda); using 'none'. |
| warning (oracle) | - | doublets: best option 'scDblFinder' is unavailable (missing: r:scDblFinder); using 'scrublet'. |
