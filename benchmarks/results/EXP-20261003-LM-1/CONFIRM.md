# Confirmatory comparisons (CONFIRM_FAMILY)

Primary condition: `clean`; other conditions are secondary. Paired t over datasets, one-sided per direction (each direction its own hypothesis), Holm over the directed hypotheses of each condition. Arms: baseline, noref_shipped, ref_geneset_preserving.

### clean

| A | B | n | excluded | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|---|
| baseline | noref_shipped | 20 | 0 | +0.089 | [+0.014, +0.163] | 0.0435 | 1 | baseline |
| ref_geneset_preserving | baseline | 20 | 0 | +0.023 | [-0.008, +0.053] | 0.201 | 1 | unresolved |

### doublets_20

| A | B | n | excluded | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|---|
| baseline | noref_shipped | 20 | 0 | +0.112 | [+0.044, +0.180] | 0.00544 | 1 | baseline |
| ref_geneset_preserving | baseline | 20 | 0 | +0.024 | [+0.008, +0.039] | 0.00613 | 1 | ref_geneset_preserving |

### damaged_15

| A | B | n | excluded | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|---|
| baseline | noref_shipped | 20 | 0 | +0.126 | [+0.046, +0.206] | 0.00778 | 1 | baseline |
| ref_geneset_preserving | baseline | 20 | 0 | +0.018 | [-0.010, +0.046] | 0.292 | 1 | unresolved |
