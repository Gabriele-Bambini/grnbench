# Model catalog: what grnbench runs, what it does not, and why

Every model is one registered method behind the same contract (`grnbench.infer`): same cleaned
data in, a ranked `source, target, weight` edge list out, scored by the same evaluation. This page
lists what is built in, how each method differs from its original, and what is left out with the
reason. Bring your own model: `examples/plugins/README.md` (D-009).

**Status of a method**

| status | meaning |
|---|---|
| default | runs with `grnbench infer --methods available` (and `grnbench run`). Fixed since the first registered experiments; adding a method never changes it. |
| extra | registered, runs with `--methods all` or by name (`--methods bicor,mrnet`). Shown as `[extra]` in `grnbench methods`. |
| not included | not built in; the reason is in section 3. A plugin can still run it. |

A method that cannot run here (missing package, program or resource) is listed with its reason, never
hidden (D-006). No score from any extra method is part of a claim until an experiment record
declares it (section 4).

## 1. Built in, default (17)

| method | what | needs |
|---|---|---|
| `pearson_signed`, `spearman_signed` | \|correlation\| with sign | numpy |
| `mutual_info` | plug-in MI on equal-frequency bins (relevance networks, Butte & Kohane 2000) | numpy |
| `clr` | MI z-scored against each gene's background (Faith et al. 2007) | numpy |
| `ppcor` | partial correlation from a Ledoit-Wolf shrunk precision matrix | scikit-learn |
| `lasso` | per-target L1 regression over the regulators | scikit-learn |
| `random_forest` | per-target random-forest importance, GENIE3's estimator without dask | scikit-learn |
| `genie3`, `grnboost2` | arboreto (Huynh-Thu et al. 2010; Moerman et al. 2019) | arboreto, dask |
| `aracne3` | the ARACNe3 program (Califano lab; native C++, not a jar) | `GRNBENCH_ARACNE3_BIN` or `ARACNe3_app_release` on PATH |
| `regdiffusion` | diffusion-model adjacency (Zhu & Slonim 2024) | torch |
| `motif_prior`, `pearson_motif`, `grnboost2_motif`, `regdiffusion_motif`, `pearson_cistarget`, `grnboost2_cistarget` | cisTarget motif prior alone, rank-multiplied with a model, or as pySCENIC's pruning step | `GRNBENCH_CISTARGET_DIR` |

## 2. Built in, extra (11)

Native numpy / scipy / scikit-learn reimplementations written from the published description. They
are **not** the original packages: the third column says where they differ. All are deterministic
(`seed` where there is randomness) and pass the same contract tests as the default methods.

| method | what, source | differs from the original |
|---|---|---|
| `bicor` | biweight midcorrelation, robust to outlier cells (WGCNA; Langfelder & Horvath, *J Stat Softw* 46(11), 2012) | genes with MAD 0 (at least half the cells at one value, usual for sparse genes) use Pearson, as WGCNA's `pearsonFallback = "individual"`. On sparse data bicor therefore equals Pearson for most pairs. |
| `wgcna_tom` | soft-thresholded \|correlation\| and topological overlap (Zhang & Horvath, *SAGMB* 4(1):17, 2005), power from the scale-free fit (`pickSoftThreshold`) | Pearson adjacency, unsigned network; if no power reaches R² 0.85 the best one is used. TOM is symmetric: it ranks pairs, it does not orient them. Constant genes are dropped before the power is chosen, as WGCNA does. Refuses more than 6000 genes (genes × genes matrix). |
| `propr_rho` | proportionality ρp = 2 cov / (var x + var y) (Quinn et al., *Sci Rep* 7:16252, 2017; best-performing family in Skinnider et al., *Nat Methods* 16, 2019) | applied to the normalised log expression, not to clr-transformed counts as `propr` does. |
| `aracne_dpi` | MI, significance threshold, data-processing-inequality pruning (ARACNe, Margolin et al., *BMC Bioinformatics* 7(Suppl 1):S7, 2006) | equal-frequency plug-in MI (not adaptive partitioning or kernels); threshold from the asymptotic G-test, not permutations; no bootstrap consolidation; DPI only over TF-TF-target triangles; an edge tied with the indirect path is removed (ARACNe's `<=`). Pruned edges tie at 0 at the bottom. For the program itself use `aracne3`. |
| `mrnet` | maximum relevance, minimum redundancy per target (Meyer et al., *EURASIP J Bioinform Syst Biol* 2007:79879) | keeps the target-wise direction; `minet` symmetrises; stops when the best score is not positive; the MI is the equal-frequency 8-bin plug-in estimator, not `minet`'s default. Unselected regulators tie at 0. |
| `tigress` | LARS stability selection with random regulator weights (Haury et al., *BMC Syst Biol* 6:145, 2012) | 50 resamples by default (CPU time); cell subsets shared by all targets of a resample; columns standardised once; the "area" counts a regulator at step L only while its coefficient is non-zero. Never-selected regulators tie at 0. |
| `bbsr` | CLR preselection and BIC best-subset regression per target, in the style of the Inferelator (Greenfield et al., *Bioinformatics* 29:1060, 2013; Skok Gibbs et al., *Bioinformatics* 38:2519, 2022) | written from the description, not the Inferelator code: no prior network, no multitask mode, plain Gaussian BIC. Unchosen regulators tie at 0. Choosing the best 10 of many regulators and then testing them with a BIC that ignores that choice overstates inclusion on noise: with 100 regulators and 1000 cells about one pure-noise target in five still gets a non-empty model, more with fewer cells (measured in review, not a calibrated rate). Refuses `n_preselect` above 16 and fewer than 2 × `n_preselect` cells. |
| `pcnet` | principal-component regression of each gene on the others, the network step of scTenifoldNet (Osorio et al., *Patterns* 1(9), 2020) | only that step (no tensor decomposition or network comparison). With 3 components it only resolves networks governed by a few axes of variation; it is checked against an independent PCR, not against a planted network of many independent modules. One partial eigendecomposition of a genes × genes matrix per gene: cost grows with the cube of the number of genes (measured: about a minute at 1000 genes on a free 4-core machine; extrapolated, ten minutes at 2000 and over an hour at 4000), so it refuses more than 4000 genes by default. |
| `extra_trees` | per-target Extra-Trees importance, GENIE3's second option (Huynh-Thu et al. 2010) | scikit-learn, no dask; 100 trees by default where `genie3` (arboreto) uses 1000; scikit-learn's `feature_importances_` sum to 1 per target, so a target with no real regulator still gives its best regulator a large weight (the existing `random_forest` shares this). |
| `gradient_boosting` | per-target gradient-boosted trees importance, GRNBoost2's idea (Moerman et al., *Bioinformatics* 35:2159, 2019) | scikit-learn; 100 trees, learning rate 0.1 and `max_features` 0.3 where arboreto uses up to 5000 trees, 0.01 and 0.1 with early stopping; importances normalised to sum to 1 per target, as for `extra_trees`. Use `grnboost2` for arboreto's scores. |
| `community` | mean rank of `spearman_signed, clr, ppcor, lasso, mrnet, bicor` (Marbach et al., *Nat Methods* 9:796, 2012) | fixed member list; a member that cannot run is an error, never silently dropped (otherwise the community would mean something different on every machine). `members=a,b,c` overrides; `seed` is passed on to every member. |

Combining with the motif prior: `benchmarks/real_bench.py --models` resolves a `_motif` or
`_cistarget` suffix on any registered base method (for example `mrnet_motif`), so every extra
method can be tried with the sequence prior without a new registration.

## 3. Not included, and why

| method | reason |
|---|---|
| PIDC (Chan et al. 2017) | needs partial information decomposition over gene triplets; the original is Julia. A native version would have to sample the third gene, which changes the method. Run the original as a command plugin. |
| SINCERITIES, SCODE, LEAP, SCRIBE, SINGE, GRISLI, dynGENIE3 | need time stamps, a pseudotime or RNA velocity; grnbench inputs are cross-sectional snapshots. |
| SCENIC (AUCell step) | AUCell scores regulon activity per cell, not a network. Its first two steps are `grnboost2` and `grnboost2_cistarget`. |
| SCENIC+, CellOracle, Pando, FigR, GRaNIE, LINGER | need chromatin accessibility or multiome data (`docs/models_that_work.md`); the input here is scRNA-seq only. |
| DeepSEM, GENELink, CNNC, scGeneRAI, scPRINT and other deep or foundation models | GPU or heavy installs; DeepSEM is not on PyPI and is pinned to old Python and torch (`docs/models_that_work.md`). `regdiffusion` is the one deep model built in. A verified reimplementation is a separate piece of work; a plugin runs the original. |
| methods that read curated TF-target collections (DoRothEA, CollecTRI, TRRUST, decoupler activity) | circular against the evaluation reference (`docs/models_that_work.md` section 2). |
| graphical lasso | cost grows with the cube of the number of genes; `ppcor` (shrunk partial correlation) and `lasso` (neighbourhood selection) cover the same idea at the engine's sizes. |

## 4. What these methods may be used for

- **Exploratory until declared.** The extra methods were checked on simulated planted networks and
  on unit properties (tests in `tests/test_literature_methods.py`), not on real data. No score of
  theirs is a result. To compare them on real data, open an experiment record that names the
  methods, the primary metric and the comparison family before the run (`experiments/README.md`);
  the pre-declared pairs of EXP-20260925-LM-2, EXP-20260926-LM-1 and EXP-20261002-LM-1 are unchanged.
- **More methods, more comparisons.** With `k` methods there are `k(k-1)/2` pairs. A confirmatory run
  declares the pairs it tests; the Holm correction runs over that family (D-008), never over
  everything a run happened to produce.
- **No leaderboard** (D-010). The point of several models side by side is to see whether a cleaning
  decision or a prior changes the ranking, not to crown a method.
