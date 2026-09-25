# Scan - grnbench
generated 2026-09-22T20:27:59+02:00 by `pos scan` (tier-0, no LLM). files 132 | code files 63 | ~11832 LOC | truncated: no
size class **M** -> mapping plan: 2x pos-mapper (cheap tier) in parallel, one per group src | tests,benchmarks,audit,hpc,examples; then 1x pos-mapper merges parts into MAP.md

## Layout
| group | files | code | ~LOC | langs |
|---|---|---|---|---|
| src | 34 | 34 | 6646 | Python |
| tests | 21 | 21 | 4023 | Python |
| audit | 51 | 3 | 562 | Python |
| benchmarks | 3 | 2 | 424 | Python |
| hpc | 5 | 2 | 105 | Shell |
| examples | 4 | 1 | 72 | Python |
| . | 9 | 0 | 0 |  |
| configs | 2 | 0 | 0 |  |
| docs | 2 | 0 | 0 |  |
| requirements | 1 | 0 | 0 |  |

Data/artifact dirs (NOT scanned - never read wholesale): OUT/ (2 entries), benchmarks/results/ (2 entries)

## Languages
Python 61 files ~11727 LOC; Shell 2 files ~105 LOC

## README (README.md)
headings: grnbench — a plugin-first GRN benchmark with bulletproof preprocessing / Why this and not BEELINE, or geneRNIB / Install / Quickstart / One command / Or stage by stage / Bring your own model / Data-hygiene coverage (what breaks in real datasets → what we do) / Hard rules (non-negotiable by design) / Inference adapters (canonical edge schema) / Evaluation / Preprocessing strategies / Input formats / HPC / Roadmap
intro: **v0.7.0 — the ranking is not real, and neither was the test that said so.** One command takes raw scRNA-seq counts to a comparison report: `preprocess → infer → evaluate → report`, every step logged, every method interchangeable. Positioning and design notes: `docs/GRN_PLATFORM_PLAN.md` in the parent project.

## Dependencies
- requirements: ["anndata", "annotated-types", "arboreto", "array-api-compat", "bokeh", "certifi", "click", "cloudpickle", "colorama", "contourpy", "cycler", "dask", "distributed", "donfig", "fast-array-utils", "fonttools", "formulaic", "fsspec", "google-crc32c", "h5py", "iniconfig", "interface-meta", "jinja2", "joblib", "kiwisolver", "legacy-api-wrap", "llvmlite", "locket", "lz4", "markupsafe", "matplotlib", "msgpack", "narwhals", "natsort", "networkx", "numba", "numcodecs", "numpy", "packaging", "pandas"]
- pyproject: {"name": "grnbench", "description": "GRN benchmarking with explicit candidate universes, provenance and calibrated comparison tools", "deps": ["anndata", "scanpy", "numpy", "pandas", "scipy"], "scripts": {"grnbench": "grnbench.cli_main:main", "grnbench-prep": "grnbench.cli:main", "grnbench-infer": "grnbench.cli_infer:main", "grnbench-eval": "grnbench.cli_eval:main", "grnbench-run": "grnbench.cli_run:main"}, "optional": ["plot", "infer", "dorothea", "plugins", "dev", "all"]}

## Entry points
- audit/2026-09-13/implementation/beeline-reference-generateExpInputs.py (__main__)
- audit/2026-09-13/reproduce_eval_findings.py (__main__)
- audit/2026-09-13/reproduce_mi_findings.py (__main__)
- benchmarks/beeline_backtest.py (__main__)
- benchmarks/plot_calibration.py (__main__)
- src/grnbench/cli.py (__main__)
- src/grnbench/cli_calibrate.py (__main__)
- src/grnbench/cli_compare.py (__main__)
- src/grnbench/cli_eval.py (__main__)
- src/grnbench/cli_infer.py (__main__)
- src/grnbench/cli_main.py (__main__)
- src/grnbench/cli_run.py (__main__)
- pyproject script `grnbench` -> grnbench.cli_main:main
- pyproject script `grnbench-prep` -> grnbench.cli:main
- pyproject script `grnbench-infer` -> grnbench.cli_infer:main
- pyproject script `grnbench-eval` -> grnbench.cli_eval:main
- pyproject script `grnbench-run` -> grnbench.cli_run:main

## Tests
20 test files

## Existing context files
AGENTS.md, CLAUDE.md, graphify-out/GRAPH_REPORT.md, CHANGELOG.md, docs (graphify report: read its summary / god-nodes first)

## Symbol index (top-level defs; notebooks: headings/imports)
### src
- src/grnbench/__init__.py (33) - "GRN-Bench: explicit preprocessing and evaluation contracts for GRN..." |
- src/grnbench/calibration.py (194) - "Score-level calibration pilot for the established paired-t/Holm ref..." | binomial_bounds(events, trials, confidence); calibration_gate(false_superiority, covered, tr...
- src/grnbench/cli.py (141) - "Command-line interface." | main(argv)
- src/grnbench/cli_calibrate.py (30) - "Run a configured score-level calibration pilot." | main(argv)
- src/grnbench/cli_compare.py (50) - "Compare scores across explicitly declared independent experimental..." | main(argv)
- src/grnbench/cli_eval.py (480) - "grnbench-eval CLI: score one or more ranked edge lists against a go..." | main(argv)
- src/grnbench/cli_infer.py (181) - "grnbench-infer CLI: run GRN inference methods on a (preprocessed) d..." | main(argv)
- src/grnbench/cli_main.py (87) - "``grnbench`` — umbrella command dispatching to the stage CLIs." | main(argv)
- src/grnbench/cli_run.py (235) - "grnbench run — the whole loop in one command." | main(argv)
- src/grnbench/eval/__init__.py (42) - "Evaluation subpackage (Parte 3): metrics, gold standards, nulls, st..." |
- src/grnbench/eval/comparison.py (122) - "Reference inference on independent experimental-unit scores." | holm_adjust(pvalues); scores_from_long(table); compare_units(scores, alpha, margin); render_comparison...
- src/grnbench/eval/gold.py (117) - "Gold-standard loading for GRN evaluation." | class GoldStandard(); load_gold(spec); coverage_report(edges_files, gold, case_sensitive)
- src/grnbench/eval/metrics.py (389) - "Ranking metrics for GRN edge lists against a gold standard." | class GoldIndex()[2m]; label_edges(edges, gold); identifier_overlap(edges, gold, n_examples); average_p...
- src/grnbench/eval/nulls.py (159) - "Explicit conditional randomization diagnostics for ranked edge lists." | random_null(edges, n, seed); regulator_stratified_nulls(edges, n, seed); degree_preserving_nu...
- src/grnbench/eval/significance.py (141) - "Descriptive annotation-thinning scenarios, not statistical inference." | class _MethodState()[1m]; gold_thinning_sensitivity(edge_sets, gold); reference_ordering_diff...
- src/grnbench/eval/stability.py (76) - "Subsampling stability: do top-k edges survive resampling the cells?" | subsample_stability(adata, method, tf_list)
- src/grnbench/eval/universe.py (114) - "Fixed candidate support and strict identifier/score validation." | canonical_pairs(frame); validate_edges(edges); load_universe(path); align_ranking(edges, universe);...
- src/grnbench/features.py (72) - "Feature selection with the TF-union guarantee." | select_features(adata, n_hvg, tf_list)
- src/grnbench/html_report.py (646) - "Self-contained HTML comparison report for a benchmarking run." | render_html(report, prep_manifest)
- src/grnbench/hygiene.py (270) - "Data-hygiene utilities for the messy realities of scRNA-seq metadata." | missing_label_mask(series); normalize_labels(series, mapping, case); load_label_mapping(path)...
- src/grnbench/infer/__init__.py (33) - "GRN inference adapters (Parte 2). Public API mirrors grnbench.strat..." |
- src/grnbench/infer/aracne3.py (134) - "ARACNe3 adapter (Java jar) with honest availability probing." | aracne3_availability(); run_aracne3(expr, tf_names, threads, timeout_s); parse_aracne3_text(text)
- src/grnbench/infer/base.py (299) - "Unified GRN inference adapter contract." | class MethodError(ValueError); class MethodSpec(); prepare_adata(adata, assume_normalized); prepare_expression(adata, assum...
- src/grnbench/infer/classical.py (251) - "Classical GRN inference methods, in numpy/scipy/scikit-learn only." | pearson_signed(expr, tf_names); mutual_info(expr, tf_names, n_bins); clr(expr, tf_names, n_bins)...
- src/grnbench/infer/correlation.py (34) - "Deterministic correlation baseline: signed absolute Spearman correl..." | spearman_signed_abs(expr, tf_names)
- src/grnbench/infer/tree_methods.py (101) - "GENIE3 / GRNBoost2 adapters over the arboreto package." | genie3_runner(expr, tf_names, seed); grnboost2_runner(expr, tf_names, seed)
- src/grnbench/ingest.py (348) - "Dataset ingestion: common scRNA-seq formats -> counts-backed AnnData." | class IngestError(ValueError); class IngestReport()[1m]; load_any(path, orientation); looks_l...
- src/grnbench/pipeline.py (427) - "End-to-end preprocessing pipeline." | class PrepError(RuntimeError); class PrepResult(); run_prep(input_path, out_dir)
- src/grnbench/plugins.py (536) - "Custom-model plugin contract — bring your own GRN method." | class PluginError(MethodError); ensure_plugins_loaded(); load_entry_point_methods(); python_method_spec(o...
- src/grnbench/qc.py (272) - "Loud-failure QC gates for GRN-bound scRNA-seq datasets." | class Thresholds(); class Gate()[2m]; class QCReport()[2m]; run_qc(adata, thresholds); pd_index_upper(idx)
- src/grnbench/recommend.py (55) - "Auto-recommendation of a preprocessing strategy from dataset statis..." | recommend_strategy(metrics, group_key, n_groups, median_cells_per_group)
- src/grnbench/report.py (205) - "Markdown run report + JSON manifest generation." | file_fingerprint(path, full_below_mb); sequence_identity(values); require_empty_output(path); environment_block();...
- src/grnbench/strategies.py (245) - "Preprocessing strategy registry." | class StrategyError(ValueError); class StrategySpec(); raw_counts(adata); cpm_log1p(adata, target_sum); pearson_residuals(adata, t...
- src/grnbench/tfs.py (127) - "Transcription-factor list utilities." | class TFListError(ValueError); resolve_tf_names(gene_names, tf_names); load_tf_list(path); match_tfs_in_data(tf_list, var_name...
### tests
- tests/conftest.py (65) - "Shared fixtures: synthetic count datasets with known properties." | synthetic_adata(); extreme_mito_adata(); cascade_adata()
- tests/test_calibration.py (80) - "Calibration gates use finite Monte Carlo bounds and expose model vi..." | config(reps); test_exact_binomial_bound_zero_errors_is_not_zero_risk(); test_pilot_reproduci...
- tests/test_classical.py (226) - "The numpy/sklearn method family." | planted(); test_method_is_registered_and_described(method); test_edge_list_obeys_the_canonical_contract(method, planted); test_met...
- tests/test_cli.py (301) - "CLI tests: the umbrella command, the one-shot `run`, and eval-with-..." | dataset(tmp_path, synthetic_adata); test_umbrella_help_and_version(capsys); test_umbrella_li...
- tests/test_cli_prep.py (94) - "Tests for the grnbench-prep CLI wiring (grnbench.cli.main)." | test_explicit_zero_threshold_is_honoured_not_replaced_by_default(tmp_path); test_nonzero_threshold_stil...
- tests/test_comparison.py (90) - "Independent-unit inference: reference calculations, not gold resamp..." | test_pair_matches_scipy_and_two_directed_hypotheses(); test_holm_matches_known_stepdown_and_...
- tests/test_eval.py (393) - "Evaluation module tests: hand-computed metrics, nulls, gold loader,..." | test_label_edges(); test_hand_computed_metric_block(); test_average_precision_and_precision_...
- tests/test_eval_cli_provenance.py (64) - "Independent-review regressions for replay-safe standalone evaluation." | inputs(tmp_path, method); test_eval_refuses_nonempty_output_before_loading_or_mutating(tmp_pa...
- tests/test_eval_contract.py (160) - "Regression invariants from the dated scientific audit." | frame(n); test_truncation_does_not_remove_gold_from_denominator(); test_invalid_scores_rejected(bad); test_d...
- tests/test_features_pipeline.py (194) - "Feature selection (TF-union rule) and end-to-end pipeline tests." | test_tf_union_rule_rescues_weak_tf(synthetic_adata); test_select_features_without_tf_list(syntheti...
- tests/test_hygiene.py (288) - "Data-hygiene tests: labels, IDs, organism, gates, stratification, m..." | test_normalize_labels_strips_and_collapses(); test_normalize_labels_case_and_synonyms(); tes...
- tests/test_infer.py (121) - "Inference adapter tests: contract, baseline recovery, tree methods,..." | test_list_methods_reports_unavailability(); test_spearman_recovers_cascade_and_sign(cascade_...
- tests/test_inference_reliability.py (134) - "Regressions for audited inference and TF-identifier failures." | test_mutual_information_of_constant_variables_is_zero(); test_mutual_information_is_invariant_to_comm...
- tests/test_ingest.py (190) - "Ingest tests: format loading, orientation detection, counts integrity." | test_load_h5ad_roundtrip(tmp_path, synthetic_adata); test_csv_orientation_detection_both_way...
- tests/test_plugins.py (433) - "Plugin contract tests: Python objects, entry points, external comma..." | clean_registry(); class DemoModel()[1m]; demo_function(adata, tf_list); test_class_plugin_ru...
- tests/test_prep_gaps.py (289)
- tests/test_provenance_reliable.py (341)
- tests/test_qc.py (55)
- tests/test_report_html.py (335)
- tests/test_significance.py (92)
- tests/test_strategies.py (78)
### audit
- audit/2026-09-13/implementation/beeline-reference-generateExpInputs.py (308) - get_parser(); parse_arguments(); filter_genes_by_pvalue(gene_df, pval_cutoff, bf_corr); select_gene_set(gene_df, variable_tfs, num_genes, sort_by_variance); filter_net...
- audit/2026-09-13/reproduce_eval_findings.py (180) - "Reproduce audit findings; these are audit fixtures, not repair/regr..." | load_eval_modules(); edges(rows); clean_json(value); main()
- audit/2026-09-13/reproduce_mi_findings.py (74) - "Read-only source audit: MI/CLR must not infer dependence from row i..." | main()
### benchmarks
- benchmarks/beeline_backtest.py (375) - "Replayable BEELINE-derived selection and grnbench fixed-universe me..." | build_parser(); regime_check(rows); unit_comparison(rows, alpha, metrics); render_summary_md...
- benchmarks/plot_calibration.py (49) - "Plot recorded calibration JSON without rerunning or selecting exper..." | main()
### hpc
- hpc/fetch_beeline.sh (65)
- hpc/setup_env.sh (40)
### examples
- examples/plugins/ridge_model.py (72) - "Example Python plugin: multivariate ridge regression of each gene o..." | class RidgeGRN()[1m]
