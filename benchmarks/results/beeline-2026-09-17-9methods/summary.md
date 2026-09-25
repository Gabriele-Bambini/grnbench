# BEELINE backtest summary

Status: complete; 198 scored results; 0 errors.

## Published-regime check (specific-chip only)

1 of 126 scores exceed 3.0x the published best (AUPRC ratio 1.1, EPR 1.2); 1 of those are small-count EPRs.

Sanity check on the evaluation code, not a leaderboard. A score far above the published best on the same reference is either a scoring/matching defect or a small-count reference (hits_at_k against expected_random_hits below 5); read the counts before quoting the ratio.

| dataset | method | metric | value | ratio to published best | hits at k / expected at random |
|---|---|---|---:|---:|---:|
| hESC | pearson_signed | auprc_ratio | 1.456 | 1.32 |  |
| hESC | pearson_signed | epr | 2.734 | 2.28 | 152 / 55.6 |
| hESC | spearman_signed | auprc_ratio | 1.447 | 1.32 |  |
| hESC | spearman_signed | epr | 2.662 | 2.22 | 148 / 55.6 |
| hESC | mutual_info | auprc_ratio | 1.366 | 1.24 |  |
| hESC | mutual_info | epr | 2.698 | 2.25 | 150 / 55.6 |
| hESC | clr | auprc_ratio | 1.246 | 1.13 |  |
| hESC | clr | epr | 1.978 | 1.65 | 110 / 55.6 |
| hESC | ppcor | auprc_ratio | 1.052 | 0.96 |  |
| hESC | ppcor | epr | 1.277 | 1.06 | 71 / 55.6 |
| hESC | lasso | auprc_ratio | 1.110 | 1.01 |  |
| hESC | lasso | epr | 2.464 | 2.05 | 137 / 55.6 |
| hHEP | pearson_signed | auprc_ratio | 1.027 | 0.93 |  |
| hHEP | pearson_signed | epr | 1.074 | 0.89 | 250 / 232.8 |
| hHEP | spearman_signed | auprc_ratio | 1.007 | 0.92 |  |
| hHEP | spearman_signed | epr | 0.971 | 0.81 | 226 / 232.8 |
| hHEP | mutual_info | auprc_ratio | 1.040 | 0.95 |  |
| hHEP | mutual_info | epr | 1.314 | 1.10 | 306 / 232.8 |
| hHEP | clr | auprc_ratio | 1.056 | 0.96 |  |
| hHEP | clr | epr | 1.293 | 1.08 | 301 / 232.8 |
| hHEP | ppcor | auprc_ratio | 0.956 | 0.87 |  |
| hHEP | ppcor | epr | 0.820 | 0.68 | 191 / 232.8 |
| hHEP | lasso | auprc_ratio | 0.996 | 0.91 |  |
| hHEP | lasso | epr | 0.876 | 0.73 | 204 / 232.8 |
| mDC | pearson_signed | auprc_ratio | 1.167 | 1.06 |  |
| mDC | pearson_signed | epr | 0.000 | 0.00 | 0 / 1.8 (small count) |
| mDC | spearman_signed | auprc_ratio | 1.082 | 0.98 |  |
| mDC | spearman_signed | epr | 1.708 | 1.42 | 3 / 1.8 (small count) |
| mDC | mutual_info | auprc_ratio | 1.081 | 0.98 |  |
| mDC | mutual_info | epr | 0.000 | 0.00 | 0 / 1.8 (small count) |
| mDC | clr | auprc_ratio | 1.149 | 1.04 |  |
| mDC | clr | epr | 4.554 | 3.79 **above regime** | 8 / 1.8 (small count) |
| mDC | ppcor | auprc_ratio | 1.276 | 1.16 |  |
| mDC | ppcor | epr | 2.846 | 2.37 | 5 / 1.8 (small count) |
| mDC | lasso | auprc_ratio | 1.043 | 0.95 |  |
| mDC | lasso | epr | 1.138 | 0.95 | 2 / 1.8 (small count) |
| mESC | pearson_signed | auprc_ratio | 1.564 | 1.42 |  |
| mESC | pearson_signed | epr | 2.223 | 1.85 | 2223 / 1000.2 |
| mESC | spearman_signed | auprc_ratio | 1.465 | 1.33 |  |
| mESC | spearman_signed | epr | 2.060 | 1.72 | 2060 / 1000.2 |
| mESC | mutual_info | auprc_ratio | 1.537 | 1.40 |  |
| mESC | mutual_info | epr | 2.234 | 1.86 | 2234 / 1000.2 |
| mESC | clr | auprc_ratio | 1.324 | 1.20 |  |
| mESC | clr | epr | 1.907 | 1.59 | 1907 / 1000.2 |
| mESC | ppcor | auprc_ratio | 1.058 | 0.96 |  |
| mESC | ppcor | epr | 1.103 | 0.92 | 1103 / 1000.2 |
| mESC | lasso | auprc_ratio | 1.101 | 1.00 |  |
| mESC | lasso | epr | 1.563 | 1.30 | 1563 / 1000.2 |
| mHSC-E | pearson_signed | auprc_ratio | 1.037 | 0.94 |  |
| mHSC-E | pearson_signed | epr | 0.818 | 0.68 | 696 / 850.8 |
| mHSC-E | spearman_signed | auprc_ratio | 1.072 | 0.97 |  |
| mHSC-E | spearman_signed | epr | 0.811 | 0.68 | 690 / 850.8 |
| mHSC-E | mutual_info | auprc_ratio | 1.028 | 0.93 |  |
| mHSC-E | mutual_info | epr | 0.738 | 0.62 | 628 / 850.8 |
| mHSC-E | clr | auprc_ratio | 0.910 | 0.83 |  |
| mHSC-E | clr | epr | 0.711 | 0.59 | 605 / 850.8 |
| mHSC-E | ppcor | auprc_ratio | 1.007 | 0.92 |  |
| mHSC-E | ppcor | epr | 1.010 | 0.84 | 859 / 850.8 |
| mHSC-E | lasso | auprc_ratio | 1.027 | 0.93 |  |
| mHSC-E | lasso | epr | 1.202 | 1.00 | 1023 / 850.8 |
| mHSC-GM | pearson_signed | auprc_ratio | 1.059 | 0.96 |  |
| mHSC-GM | pearson_signed | epr | 1.009 | 0.84 | 578 / 573.1 |
| mHSC-GM | spearman_signed | auprc_ratio | 1.035 | 0.94 |  |
| mHSC-GM | spearman_signed | epr | 0.916 | 0.76 | 525 / 573.1 |
| mHSC-GM | mutual_info | auprc_ratio | 1.015 | 0.92 |  |
| mHSC-GM | mutual_info | epr | 0.776 | 0.65 | 445 / 573.1 |
| mHSC-GM | clr | auprc_ratio | 0.925 | 0.84 |  |
| mHSC-GM | clr | epr | 0.736 | 0.61 | 422 / 573.1 |
| mHSC-GM | ppcor | auprc_ratio | 1.051 | 0.96 |  |
| mHSC-GM | ppcor | epr | 1.091 | 0.91 | 625 / 573.1 |
| mHSC-GM | lasso | auprc_ratio | 1.121 | 1.02 |  |
| mHSC-GM | lasso | epr | 1.640 | 1.37 | 940 / 573.1 |
| mHSC-L | pearson_signed | auprc_ratio | 1.142 | 1.04 |  |
| mHSC-L | pearson_signed | epr | 1.301 | 1.08 | 617 / 474.2 |
| mHSC-L | spearman_signed | auprc_ratio | 1.174 | 1.07 |  |
| mHSC-L | spearman_signed | epr | 1.362 | 1.14 | 646 / 474.2 |
| mHSC-L | mutual_info | auprc_ratio | 1.311 | 1.19 |  |
| mHSC-L | mutual_info | epr | 1.451 | 1.21 | 688 / 474.2 |
| mHSC-L | clr | auprc_ratio | 1.120 | 1.02 |  |
| mHSC-L | clr | epr | 1.166 | 0.97 | 553 / 474.2 |
| mHSC-L | ppcor | auprc_ratio | 1.035 | 0.94 |  |
| mHSC-L | ppcor | epr | 1.061 | 0.88 | 503 / 474.2 |
| mHSC-L | lasso | auprc_ratio | 1.274 | 1.16 |  |
| mHSC-L | lasso | epr | 1.826 | 1.52 | 866 / 474.2 |
| hESC | random_forest | auprc_ratio | 1.483 | 1.35 |  |
| hESC | random_forest | epr | 3.399 | 2.83 | 189 / 55.6 |
| hESC | grnboost2 | auprc_ratio | 1.389 | 1.26 |  |
| hESC | grnboost2 | epr | 3.004 | 2.50 | 167 / 55.6 |
| hESC | genie3 | auprc_ratio | 1.434 | 1.30 |  |
| hESC | genie3 | epr | 2.698 | 2.25 | 150 / 55.6 |
| hHEP | random_forest | auprc_ratio | 1.004 | 0.91 |  |
| hHEP | random_forest | epr | 1.310 | 1.09 | 305 / 232.8 |
| hHEP | grnboost2 | auprc_ratio | 1.112 | 1.01 |  |
| hHEP | grnboost2 | epr | 1.357 | 1.13 | 316 / 232.8 |
| hHEP | genie3 | auprc_ratio | 1.021 | 0.93 |  |
| hHEP | genie3 | epr | 1.258 | 1.05 | 293 / 232.8 |
| mDC | random_forest | auprc_ratio | 1.097 | 1.00 |  |
| mDC | random_forest | epr | 1.138 | 0.95 | 2 / 1.8 (small count) |
| mDC | grnboost2 | auprc_ratio | 1.316 | 1.20 |  |
| mDC | grnboost2 | epr | 2.846 | 2.37 | 5 / 1.8 (small count) |
| mDC | genie3 | auprc_ratio | 1.133 | 1.03 |  |
| mDC | genie3 | epr | 1.708 | 1.42 | 3 / 1.8 (small count) |
| mESC | random_forest | auprc_ratio | 1.383 | 1.26 |  |
| mESC | random_forest | epr | 2.224 | 1.85 | 2224 / 1000.2 |
| mESC | grnboost2 | auprc_ratio | 1.328 | 1.21 |  |
| mESC | grnboost2 | epr | 2.390 | 1.99 | 2390 / 1000.2 |
| mESC | genie3 | auprc_ratio | 1.465 | 1.33 |  |
| mESC | genie3 | epr | 2.474 | 2.06 | 2474 / 1000.2 |
| mHSC-E | random_forest | auprc_ratio | 1.106 | 1.01 |  |
| mHSC-E | random_forest | epr | 0.649 | 0.54 | 552 / 850.8 |
| mHSC-E | grnboost2 | auprc_ratio | 1.011 | 0.92 |  |
| mHSC-E | grnboost2 | epr | 0.719 | 0.60 | 612 / 850.8 |
| mHSC-E | genie3 | auprc_ratio | 1.107 | 1.01 |  |
| mHSC-E | genie3 | epr | 0.582 | 0.48 | 495 / 850.8 |
| mHSC-GM | random_forest | auprc_ratio | 1.395 | 1.27 |  |
| mHSC-GM | random_forest | epr | 1.209 | 1.01 | 693 / 573.1 |
| mHSC-GM | grnboost2 | auprc_ratio | 1.114 | 1.01 |  |
| mHSC-GM | grnboost2 | epr | 1.026 | 0.86 | 588 / 573.1 |
| mHSC-GM | genie3 | auprc_ratio | 1.421 | 1.29 |  |
| mHSC-GM | genie3 | epr | 1.202 | 1.00 | 689 / 573.1 |
| mHSC-L | random_forest | auprc_ratio | 1.467 | 1.33 |  |
| mHSC-L | random_forest | epr | 1.565 | 1.30 | 742 / 474.2 |
| mHSC-L | grnboost2 | auprc_ratio | 1.345 | 1.22 |  |
| mHSC-L | grnboost2 | epr | 1.535 | 1.28 | 728 / 474.2 |
| mHSC-L | genie3 | auprc_ratio | 1.501 | 1.36 |  |
| mHSC-L | genie3 | epr | 1.632 | 1.36 | 774 / 474.2 |

## Dataset-unit comparison

Units: datasets. 0 of 216 method pairs resolved at alpha=0.05.

Datasets are the only independent replicates in this design; gold edges and cells are not units. Unresolved pairs are not shown equivalent.

### lofgof / epr - n=1 datasets

Too few complete units for inference (1 < 3).

### lofgof / auprc_ratio - n=1 datasets

Too few complete units for inference (1 < 3).

### nonspecific-chip / epr - n=7 datasets

Units: 7; alpha: 0.05; margin: 0.0.
Paired t tests with Holm adjustment over both directions of all pairs.
Intervals use simultaneous Bonferroni coverage; unresolved does not mean equivalent.

| A | B | mean A−B | simultaneous interval | supported winner |
|---|---|---:|---|---|
| clr | genie3 | 0.560596 | [-0.0801308, 1.20132] | unresolved |
| clr | grnboost2 | 2.40467 | [-3.37238, 8.18173] | unresolved |
| clr | lasso | 3.28728 | [-3.698, 10.2726] | unresolved |
| clr | mutual_info | 0.210277 | [-3.63959, 4.06015] | unresolved |
| clr | pearson_signed | 0.214647 | [-2.53266, 2.96196] | unresolved |
| clr | ppcor | 5.15398 | [-5.05722, 15.3652] | unresolved |
| clr | random_forest | 0.806582 | [-0.720074, 2.33324] | unresolved |
| clr | spearman_signed | 0.028574 | [-2.20791, 2.26505] | unresolved |
| genie3 | grnboost2 | 1.84408 | [-3.47387, 7.16202] | unresolved |
| genie3 | lasso | 2.72668 | [-3.95331, 9.40668] | unresolved |
| genie3 | mutual_info | -0.350319 | [-4.30947, 3.60883] | unresolved |
| genie3 | pearson_signed | -0.345949 | [-3.276, 2.5841] | unresolved |
| genie3 | ppcor | 4.59339 | [-5.25547, 14.4422] | unresolved |
| genie3 | random_forest | 0.245987 | [-0.91984, 1.41181] | unresolved |
| genie3 | spearman_signed | -0.532022 | [-3.15089, 2.08685] | unresolved |
| grnboost2 | lasso | 0.882606 | [-2.23295, 3.99816] | unresolved |
| grnboost2 | mutual_info | -2.1944 | [-10.8531, 6.46435] | unresolved |
| grnboost2 | pearson_signed | -2.19002 | [-9.76571, 5.38566] | unresolved |
| grnboost2 | ppcor | 2.74931 | [-3.44765, 8.94628] | unresolved |
| grnboost2 | random_forest | -1.59809 | [-7.00639, 3.81022] | unresolved |
| grnboost2 | spearman_signed | -2.3761 | [-9.86558, 5.11338] | unresolved |
| lasso | mutual_info | -3.077 | [-13.5873, 7.43331] | unresolved |
| lasso | pearson_signed | -3.07263 | [-12.4133, 6.26801] | unresolved |
| lasso | ppcor | 1.86671 | [-2.13679, 5.87021] | unresolved |
| lasso | random_forest | -2.4807 | [-9.01373, 4.05234] | unresolved |
| lasso | spearman_signed | -3.2587 | [-12.1856, 5.6682] | unresolved |
| mutual_info | pearson_signed | 0.00437006 | [-1.74849, 1.75723] | unresolved |
| mutual_info | ppcor | 4.94371 | [-8.59144, 18.4789] | unresolved |
| mutual_info | random_forest | 0.596306 | [-3.79323, 4.98585] | unresolved |
| mutual_info | spearman_signed | -0.181703 | [-2.72789, 2.36448] | unresolved |
| pearson_signed | ppcor | 4.93934 | [-7.72197, 17.6006] | unresolved |
| pearson_signed | random_forest | 0.591936 | [-3.14373, 4.3276] | unresolved |
| pearson_signed | spearman_signed | -0.186073 | [-1.34352, 0.971372] | unresolved |
| ppcor | random_forest | -4.3474 | [-13.6868, 4.99203] | unresolved |
| ppcor | spearman_signed | -5.12541 | [-17.4211, 7.17031] | unresolved |
| random_forest | spearman_signed | -0.778008 | [-4.2263, 2.67028] | unresolved |

### nonspecific-chip / auprc_ratio - n=7 datasets

Units: 7; alpha: 0.05; margin: 0.0.
Paired t tests with Holm adjustment over both directions of all pairs.
Intervals use simultaneous Bonferroni coverage; unresolved does not mean equivalent.

| A | B | mean A−B | simultaneous interval | supported winner |
|---|---|---:|---|---|
| clr | genie3 | 0.254788 | [-0.493476, 1.00305] | unresolved |
| clr | grnboost2 | 0.948656 | [-1.35238, 3.24969] | unresolved |
| clr | lasso | 1.20975 | [-1.59666, 4.01616] | unresolved |
| clr | mutual_info | 0.203198 | [-0.838013, 1.24441] | unresolved |
| clr | pearson_signed | 0.416575 | [-0.282898, 1.11605] | unresolved |
| clr | ppcor | 1.60705 | [-2.47646, 5.69055] | unresolved |
| clr | random_forest | 0.38131 | [-0.673103, 1.43572] | unresolved |
| clr | spearman_signed | 0.230894 | [-0.602115, 1.0639] | unresolved |
| genie3 | grnboost2 | 0.693868 | [-1.14486, 2.5326] | unresolved |
| genie3 | lasso | 0.954964 | [-1.27594, 3.18587] | unresolved |
| genie3 | mutual_info | -0.0515903 | [-1.61162, 1.50843] | unresolved |
| genie3 | pearson_signed | 0.161787 | [-0.735168, 1.05874] | unresolved |
| genie3 | ppcor | 1.35226 | [-2.09248, 4.79699] | unresolved |
| genie3 | random_forest | 0.126522 | [-0.352651, 0.605696] | unresolved |
| genie3 | spearman_signed | -0.0238938 | [-1.41503, 1.36724] | unresolved |
| grnboost2 | lasso | 0.261096 | [-0.440295, 0.962487] | unresolved |
| grnboost2 | mutual_info | -0.745458 | [-3.62063, 2.12971] | unresolved |
| grnboost2 | pearson_signed | -0.532081 | [-2.33318, 1.26902] | unresolved |
| grnboost2 | ppcor | 0.658389 | [-1.33061, 2.64738] | unresolved |
| grnboost2 | random_forest | -0.567346 | [-2.36372, 1.22902] | unresolved |
| grnboost2 | spearman_signed | -0.717762 | [-3.014, 1.57848] | unresolved |
| lasso | mutual_info | -1.00655 | [-4.42855, 2.41545] | unresolved |
| lasso | pearson_signed | -0.793177 | [-3.17491, 1.58855] | unresolved |
| lasso | ppcor | 0.397293 | [-0.934797, 1.72938] | unresolved |
| lasso | random_forest | -0.828442 | [-2.96984, 1.31296] | unresolved |
| lasso | spearman_signed | -0.978858 | [-3.89053, 1.93281] | unresolved |
| mutual_info | pearson_signed | 0.213378 | [-0.94433, 1.37109] | unresolved |
| mutual_info | ppcor | 1.40385 | [-3.2994, 6.1071] | unresolved |
| mutual_info | random_forest | 0.178112 | [-1.39057, 1.74679] | unresolved |
| mutual_info | spearman_signed | 0.0276965 | [-1.15857, 1.21397] | unresolved |
| pearson_signed | ppcor | 1.19047 | [-2.50741, 4.88835] | unresolved |
| pearson_signed | random_forest | -0.0352652 | [-1.06991, 0.999375] | unresolved |
| pearson_signed | spearman_signed | -0.185681 | [-0.882259, 0.510897] | unresolved |
| ppcor | random_forest | -1.22573 | [-4.54855, 2.09708] | unresolved |
| ppcor | spearman_signed | -1.37615 | [-5.6139, 2.8616] | unresolved |
| random_forest | spearman_signed | -0.150416 | [-1.79396, 1.49313] | unresolved |

### specific-chip / epr - n=7 datasets

Units: 7; alpha: 0.05; margin: 0.0.
Paired t tests with Holm adjustment over both directions of all pairs.
Intervals use simultaneous Bonferroni coverage; unresolved does not mean equivalent.

| A | B | mean A−B | simultaneous interval | supported winner |
|---|---|---:|---|---|
| clr | genie3 | 0.113082 | [-2.51949, 2.74565] | unresolved |
| clr | grnboost2 | -0.075965 | [-1.88331, 1.73138] | unresolved |
| clr | lasso | 0.233582 | [-2.9148, 3.38196] | unresolved |
| clr | mutual_info | 0.447749 | [-3.41649, 4.31199] | unresolved |
| clr | pearson_signed | 0.455358 | [-3.41522, 4.32593] | unresolved |
| clr | ppcor | 0.448312 | [-1.07056, 1.96718] | unresolved |
| clr | random_forest | 0.121612 | [-3.11577, 3.35899] | unresolved |
| clr | spearman_signed | 0.265162 | [-2.21908, 2.7494] | unresolved |
| genie3 | grnboost2 | -0.189047 | [-1.14072, 0.762621] | unresolved |
| genie3 | lasso | 0.120499 | [-1.05896, 1.29995] | unresolved |
| genie3 | mutual_info | 0.334666 | [-1.01128, 1.68061] | unresolved |
| genie3 | pearson_signed | 0.342276 | [-0.994219, 1.67877] | unresolved |
| genie3 | ppcor | 0.335229 | [-1.6181, 2.28856] | unresolved |
| genie3 | random_forest | 0.00852954 | [-0.80566, 0.822719] | unresolved |
| genie3 | spearman_signed | 0.15208 | [-0.321322, 0.625481] | unresolved |
| grnboost2 | lasso | 0.309547 | [-1.44816, 2.06725] | unresolved |
| grnboost2 | mutual_info | 0.523714 | [-1.65481, 2.70224] | unresolved |
| grnboost2 | pearson_signed | 0.531323 | [-1.64725, 2.7099] | unresolved |
| grnboost2 | ppcor | 0.524277 | [-1.04723, 2.09579] | unresolved |
| grnboost2 | random_forest | 0.197577 | [-1.26437, 1.65953] | unresolved |
| grnboost2 | spearman_signed | 0.341127 | [-0.480907, 1.16316] | unresolved |
| lasso | mutual_info | 0.214167 | [-1.22344, 1.65178] | unresolved |
| lasso | pearson_signed | 0.221777 | [-1.09333, 1.53688] | unresolved |
| lasso | ppcor | 0.21473 | [-1.74265, 2.17211] | unresolved |
| lasso | random_forest | -0.11197 | [-1.32418, 1.10024] | unresolved |
| lasso | spearman_signed | 0.0315802 | [-1.02723, 1.09039] | unresolved |
| mutual_info | pearson_signed | 0.0076095 | [-0.317091, 0.33231] | unresolved |
| mutual_info | ppcor | 0.000562662 | [-2.98533, 2.98645] | unresolved |
| mutual_info | random_forest | -0.326137 | [-1.29191, 0.639641] | unresolved |
| mutual_info | spearman_signed | -0.182587 | [-1.64361, 1.27844] | unresolved |
| pearson_signed | ppcor | -0.00704684 | [-2.95215, 2.93806] | unresolved |
| pearson_signed | random_forest | -0.333746 | [-1.26078, 0.593282] | unresolved |
| pearson_signed | spearman_signed | -0.190196 | [-1.61314, 1.23275] | unresolved |
| ppcor | random_forest | -0.3267 | [-2.85318, 2.19978] | unresolved |
| ppcor | spearman_signed | -0.18315 | [-1.92625, 1.55995] | unresolved |
| random_forest | spearman_signed | 0.14355 | [-0.727674, 1.01477] | unresolved |

### specific-chip / auprc_ratio - n=7 datasets

Units: 7; alpha: 0.05; margin: 0.0.
Paired t tests with Holm adjustment over both directions of all pairs.
Intervals use simultaneous Bonferroni coverage; unresolved does not mean equivalent.

| A | B | mean A−B | simultaneous interval | supported winner |
|---|---|---:|---|---|
| clr | genie3 | -0.193388 | [-0.603144, 0.216368] | unresolved |
| clr | grnboost2 | -0.126376 | [-0.290293, 0.0375398] | unresolved |
| clr | lasso | 0.00846581 | [-0.334831, 0.351763] | unresolved |
| clr | mutual_info | -0.0923445 | [-0.309606, 0.124917] | unresolved |
| clr | pearson_signed | -0.103243 | [-0.319866, 0.11338] | unresolved |
| clr | ppcor | 0.0421427 | [-0.297021, 0.381307] | unresolved |
| clr | random_forest | -0.17216 | [-0.592949, 0.248628] | unresolved |
| clr | spearman_signed | -0.0787784 | [-0.298678, 0.141121] | unresolved |
| genie3 | grnboost2 | 0.0670116 | [-0.277532, 0.411555] | unresolved |
| genie3 | lasso | 0.201854 | [-0.0847943, 0.488502] | unresolved |
| genie3 | mutual_info | 0.101043 | [-0.231759, 0.433846] | unresolved |
| genie3 | pearson_signed | 0.0901453 | [-0.313849, 0.494139] | unresolved |
| genie3 | ppcor | 0.235531 | [-0.247371, 0.718432] | unresolved |
| genie3 | random_forest | 0.0212275 | [-0.0625479, 0.105003] | unresolved |
| genie3 | spearman_signed | 0.11461 | [-0.23994, 0.46916] | unresolved |
| grnboost2 | lasso | 0.134842 | [-0.13229, 0.401975] | unresolved |
| grnboost2 | mutual_info | 0.0340318 | [-0.249507, 0.31757] | unresolved |
| grnboost2 | pearson_signed | 0.0231336 | [-0.28893, 0.335197] | unresolved |
| grnboost2 | ppcor | 0.168519 | [-0.122813, 0.459852] | unresolved |
| grnboost2 | random_forest | -0.0457841 | [-0.391111, 0.299543] | unresolved |
| grnboost2 | spearman_signed | 0.047598 | [-0.240511, 0.335707] | unresolved |
| lasso | mutual_info | -0.10081 | [-0.487106, 0.285485] | unresolved |
| lasso | pearson_signed | -0.111709 | [-0.57233, 0.348913] | unresolved |
| lasso | ppcor | 0.0336769 | [-0.259793, 0.327147] | unresolved |
| lasso | random_forest | -0.180626 | [-0.470272, 0.109019] | unresolved |
| lasso | spearman_signed | -0.0872442 | [-0.486828, 0.31234] | unresolved |
| mutual_info | pearson_signed | -0.0108982 | [-0.196377, 0.174581] | unresolved |
| mutual_info | ppcor | 0.134487 | [-0.356862, 0.625836] | unresolved |
| mutual_info | random_forest | -0.0798159 | [-0.434935, 0.275303] | unresolved |
| mutual_info | spearman_signed | 0.0135662 | [-0.1425, 0.169633] | unresolved |
| pearson_signed | ppcor | 0.145385 | [-0.327908, 0.618679] | unresolved |
| pearson_signed | random_forest | -0.0689178 | [-0.481971, 0.344135] | unresolved |
| pearson_signed | spearman_signed | 0.0244643 | [-0.0854954, 0.134424] | unresolved |
| ppcor | random_forest | -0.214303 | [-0.702217, 0.27361] | unresolved |
| ppcor | spearman_signed | -0.120921 | [-0.580429, 0.338587] | unresolved |
| random_forest | spearman_signed | 0.0933821 | [-0.256663, 0.443427] | unresolved |

### string / epr - n=7 datasets

Units: 7; alpha: 0.05; margin: 0.0.
Paired t tests with Holm adjustment over both directions of all pairs.
Intervals use simultaneous Bonferroni coverage; unresolved does not mean equivalent.

| A | B | mean A−B | simultaneous interval | supported winner |
|---|---|---:|---|---|
| clr | genie3 | 1.61783 | [-2.52195, 5.75761] | unresolved |
| clr | grnboost2 | 6.645 | [-16.469, 29.759] | unresolved |
| clr | lasso | 8.58558 | [-21.7168, 38.8879] | unresolved |
| clr | mutual_info | 0.845768 | [-5.00579, 6.69733] | unresolved |
| clr | pearson_signed | 1.40214 | [-5.41492, 8.21919] | unresolved |
| clr | ppcor | 13.1409 | [-22.1092, 48.391] | unresolved |
| clr | random_forest | 2.0689 | [-2.84094, 6.97873] | unresolved |
| clr | spearman_signed | 1.46238 | [-4.92724, 7.85199] | unresolved |
| genie3 | grnboost2 | 5.02717 | [-14.2889, 24.3432] | unresolved |
| genie3 | lasso | 6.96775 | [-19.721, 33.6565] | unresolved |
| genie3 | mutual_info | -0.77206 | [-7.15675, 5.61263] | unresolved |
| genie3 | pearson_signed | -0.215691 | [-5.20043, 4.76905] | unresolved |
| genie3 | ppcor | 11.5231 | [-20.0036, 43.0497] | unresolved |
| genie3 | random_forest | 0.451068 | [-1.2748, 2.17694] | unresolved |
| genie3 | spearman_signed | -0.155453 | [-4.66247, 4.35156] | unresolved |
| grnboost2 | lasso | 1.94058 | [-6.14572, 10.0269] | unresolved |
| grnboost2 | mutual_info | -5.79923 | [-29.5432, 17.9447] | unresolved |
| grnboost2 | pearson_signed | -5.24287 | [-24.5614, 14.0756] | unresolved |
| grnboost2 | ppcor | 6.49589 | [-7.10387, 20.0956] | unresolved |
| grnboost2 | random_forest | -4.57611 | [-24.0023, 14.85] | unresolved |
| grnboost2 | spearman_signed | -5.18263 | [-24.2928, 13.9275] | unresolved |
| lasso | mutual_info | -7.73981 | [-39.1286, 23.6489] | unresolved |
| lasso | pearson_signed | -7.18344 | [-33.9998, 19.6329] | unresolved |
| lasso | ppcor | 4.55531 | [-3.45788, 12.5685] | unresolved |
| lasso | random_forest | -6.51668 | [-33.3025, 20.2691] | unresolved |
| lasso | spearman_signed | -7.12321 | [-33.6007, 19.3543] | unresolved |
| mutual_info | pearson_signed | 0.556369 | [-4.32166, 5.4344] | unresolved |
| mutual_info | ppcor | 12.2951 | [-24.3115, 48.9017] | unresolved |
| mutual_info | random_forest | 1.22313 | [-5.99753, 8.44378] | unresolved |
| mutual_info | spearman_signed | 0.616607 | [-4.75658, 5.98979] | unresolved |
| pearson_signed | ppcor | 11.7388 | [-20.6638, 44.1413] | unresolved |
| pearson_signed | random_forest | 0.666759 | [-5.68069, 7.01421] | unresolved |
| pearson_signed | spearman_signed | 0.0602382 | [-1.10029, 1.22077] | unresolved |
| ppcor | random_forest | -11.072 | [-42.356, 20.212] | unresolved |
| ppcor | spearman_signed | -11.6785 | [-43.7604, 20.4034] | unresolved |
| random_forest | spearman_signed | -0.606521 | [-6.5371, 5.32406] | unresolved |

### string / auprc_ratio - n=7 datasets

Units: 7; alpha: 0.05; margin: 0.0.
Paired t tests with Holm adjustment over both directions of all pairs.
Intervals use simultaneous Bonferroni coverage; unresolved does not mean equivalent.

| A | B | mean A−B | simultaneous interval | supported winner |
|---|---|---:|---|---|
| clr | genie3 | 3.38911 | [-11.076, 17.8542] | unresolved |
| clr | grnboost2 | 6.98695 | [-21.8702, 35.8441] | unresolved |
| clr | lasso | 8.13054 | [-23.0633, 39.3244] | unresolved |
| clr | mutual_info | 0.573179 | [-4.02226, 5.16862] | unresolved |
| clr | pearson_signed | 3.52116 | [-14.3181, 21.3604] | unresolved |
| clr | ppcor | 9.6131 | [-24.8161, 44.0423] | unresolved |
| clr | random_forest | 3.58197 | [-11.2174, 18.3813] | unresolved |
| clr | spearman_signed | 3.61336 | [-15.7775, 23.0042] | unresolved |
| genie3 | grnboost2 | 3.59784 | [-10.8072, 18.0029] | unresolved |
| genie3 | lasso | 4.74143 | [-12.078, 21.5609] | unresolved |
| genie3 | mutual_info | -2.81593 | [-15.9374, 10.3055] | unresolved |
| genie3 | pearson_signed | 0.132049 | [-5.2217, 5.4858] | unresolved |
| genie3 | ppcor | 6.224 | [-13.9355, 26.3835] | unresolved |
| genie3 | random_forest | 0.192866 | [-0.476976, 0.862707] | unresolved |
| genie3 | spearman_signed | 0.224257 | [-6.43051, 6.87903] | unresolved |
| grnboost2 | lasso | 1.14359 | [-2.05658, 4.34375] | unresolved |
| grnboost2 | mutual_info | -6.41377 | [-33.6191, 20.7916] | unresolved |
| grnboost2 | pearson_signed | -3.46579 | [-15.9134, 8.98185] | unresolved |
| grnboost2 | ppcor | 2.62615 | [-4.24552, 9.49783] | unresolved |
| grnboost2 | random_forest | -3.40498 | [-17.5186, 10.7086] | unresolved |
| grnboost2 | spearman_signed | -3.37358 | [-14.6559, 7.90876] | unresolved |
| lasso | mutual_info | -7.55736 | [-37.3874, 22.2727] | unresolved |
| lasso | pearson_signed | -4.60938 | [-20.0202, 10.8014] | unresolved |
| lasso | ppcor | 1.48257 | [-2.40839, 5.37352] | unresolved |
| lasso | random_forest | -4.54856 | [-21.1454, 12.0482] | unresolved |
| lasso | spearman_signed | -4.51717 | [-18.8039, 9.76955] | unresolved |
| mutual_info | pearson_signed | 2.94798 | [-12.457, 18.353] | unresolved |
| mutual_info | ppcor | 9.03992 | [-24.1548, 42.2346] | unresolved |
| mutual_info | random_forest | 3.00879 | [-10.2605, 16.2781] | unresolved |
| mutual_info | spearman_signed | 3.04019 | [-13.9611, 20.0415] | unresolved |
| pearson_signed | ppcor | 6.09195 | [-12.9838, 25.1677] | unresolved |
| pearson_signed | random_forest | 0.060817 | [-4.79201, 4.91365] | unresolved |
| pearson_signed | spearman_signed | 0.0922085 | [-1.65927, 1.84369] | unresolved |
| ppcor | random_forest | -6.03113 | [-25.9857, 13.9234] | unresolved |
| ppcor | spearman_signed | -5.99974 | [-24.0335, 12.034] | unresolved |
| random_forest | spearman_signed | 0.0313915 | [-6.16591, 6.2287] | unresolved |

## Assumptions (all comparisons above)

- Rows are independent experimental units, paired across all methods.
- Student t inference is exact for normally distributed unit differences; otherwise approximate.
- Methods, score, comparison family and practical margin were specified before evaluation.
- This procedure does not correct incomplete labels or establish direct causal regulation.

