# EXP-20261018-GB-1: ARACNe3 smoke run (written by benchmarks/aracne3_smoke.py)

- code: `6927b9fd353541f1a9e006f7438f84bec63f0de5`; ARACNe3 source `3d8791a23e3bd8fd0d74f3b8d48f912e81d00f14`, program sha256 `29e6e77d09137dcaed428d24021c1f88708f986c7c9f606a0906873b7f3fca03`
- subset: 500 of 2700 cells, 300 genes, 30 regulators (subset seed 0); CPUs given to the job: 4
- relaxed rerun (`--alpha 0.2`): no; private-text scan: clean

| run | threads | status | exit | seconds | samples read | edges | consolidated rows | null model sha256 |
|---|---|---|---|---|---|---|---|---|
| t1_a | 1 | ok | 0 | 5.901 | 500 | 478 | 478 | 5a1048b4f1e9,6c3a7c896714 |
| t1_b | 1 | ok | 0 | 5.665 | 500 | 478 | 478 | 5a1048b4f1e9,6c3a7c896714 |
| t4_a | 4 | ok | 0 | 2.095 | 500 | 478 | 478 | f73ce545a0b4,2918fdd23168 |
| t4_b | 4 | ok | 0 | 2.101 | 500 | 478 | 478 | 6009dffbaae9,b70ec492a271 |

| a | b | same file | same null model | edges a / b | shared | Jaccard | max abs weight diff | identical edges |
|---|---|---|---|---|---|---|---|---|
| t1_a | t1_b | True | True | 478 / 478 | 478 | 1.0 | 0.0 | True |
| t1_a | t4_a | True | False | 478 / 478 | 478 | 1.0 | 0.0 | True |
| t4_a | t4_b | True | False | 478 / 478 | 478 | 1.0 | 0.0 | True |

| run | regulators with edges | MI min / median / max | rows with MI 0 | Spearman + / - / 0 | count.values | distinct log.p.values |
|---|---|---|---|---|---|---|
| t1_a | 30 | 0.0 / 0.056187249999999994 / 0.226802 | 6 | 279 / 199 / 0 | [1] | 1 |
| t1_b | 30 | 0.0 / 0.056187249999999994 / 0.226802 | 6 | 279 / 199 / 0 | [1] | 1 |
| t4_a | 30 | 0.0 / 0.056187249999999994 / 0.226802 | 6 | 279 / 199 / 0 | [1] | 1 |
| t4_b | 30 | 0.0 / 0.056187249999999994 / 0.226802 | 6 | 279 / 199 / 0 | [1] | 1 |
