# Confirmatory comparisons (CONFIRM_FAMILY)

Primary condition: `clean`; other conditions are secondary. Paired t over datasets, two-sided per direction, Holm over the directed hypotheses of each condition. Arms: baseline, noref_shipped, ref_geneset_preserving.

### clean

| A | B | n | excluded | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|---|
| baseline | noref_shipped | 20 | 0 | +0.124 | [+0.035, +0.212] | 0.0174 | 1 | baseline |
| ref_geneset_preserving | baseline | 20 | 0 | +0.019 | [+0.002, +0.036] | 0.0465 | 1 | ref_geneset_preserving |

### depth_30

| A | B | n | excluded | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|---|
| baseline | noref_shipped | 19 | 1 | +0.159 | [+0.046, +0.272] | 0.0171 | 1 | baseline |
| ref_geneset_preserving | baseline | 19 | 1 | +0.054 | [+0.001, +0.107] | 0.0719 | 1 | unresolved |

### depth_10

| A | B | n | excluded | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|---|
| baseline | noref_shipped | 13 | 7 | +0.085 | [-0.051, +0.221] | 0.296 | 1 | unresolved |
| ref_geneset_preserving | baseline | 13 | 7 | +0.016 | [-0.003, +0.035] | 0.184 | 1 | unresolved |

### ambient_20

| A | B | n | excluded | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|---|
| baseline | noref_shipped | 20 | 0 | +0.156 | [+0.080, +0.232] | 0.00076 | 1 | baseline |
| ref_geneset_preserving | baseline | 20 | 0 | +0.010 | [-0.007, +0.028] | 0.339 | 1 | unresolved |

### doublets_20

| A | B | n | excluded | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|---|
| baseline | noref_shipped | 20 | 0 | +0.155 | [+0.065, +0.246] | 0.0037 | 1 | baseline |
| ref_geneset_preserving | baseline | 20 | 0 | +0.025 | [+0.004, +0.045] | 0.0288 | 1 | ref_geneset_preserving |

### damaged_15

| A | B | n | excluded | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|---|
| baseline | noref_shipped | 20 | 0 | +0.187 | [+0.092, +0.282] | 0.00118 | 1 | baseline |
| ref_geneset_preserving | baseline | 20 | 0 | +0.018 | [+0.004, +0.033] | 0.0259 | 1 | ref_geneset_preserving |

### mixed

| A | B | n | excluded | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|---|
| baseline | noref_shipped | 18 | 1 | +0.220 | [+0.123, +0.317] | 0.00034 | 1 | baseline |
| ref_geneset_preserving | baseline | 18 | 1 | +0.029 | [+0.001, +0.057] | 0.0697 | 1 | unresolved |
