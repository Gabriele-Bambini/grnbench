# Score-level calibration pilot

This is not biological validation or evidence of novelty.

| scenario | independent units / rows | FWER (upper bound) | coverage (lower bound) | power | gate |
|---|---:|---:|---:|---:|---|
| iid_global_null_n8 | 8 / 8 | 0.015 (0.038) | 0.960 (0.929) | n/a (no alternative) | pass |
| iid_global_null_n32 | 32 / 32 | 0.000 (0.015) | 0.980 (0.955) | n/a (no alternative) | pass |
| iid_partial_null_n16 | 16 / 16 | 0.000 (0.015) | 0.965 (0.935) | 0.079 | pass |
| iid_partial_null_n32 | 32 / 32 | 0.010 (0.031) | 0.945 (0.911) | 0.226 | inconclusive |
| iid_margin_boundary_n32 | 32 / 32 | 0.035 (0.065) | 0.955 (0.923) | n/a (no alternative) | inconclusive |
| pseudoreplication_null_8x4 | 8 / 32 | 0.470 (0.531) | 0.375 (0.318) | n/a (no alternative) | fail |
| pseudoreplication_partial_8x4 | 8 / 32 | 0.330 (0.389) | 0.345 (0.289) | 0.395 | fail |

Bounds: separate one-sided exact 95% Monte Carlo limits.
Gate: error upper <= 0.075 and coverage lower >= 0.925.
Finite precision may yield an inconclusive gate. A failure in repeated-unit controls is expected.

- No real expression matrix, GRN simulator, or incomplete-label model is tested.
- Means are known exactly in this score-level model; this does not hold for real algorithm rankings.
- Each gate is scenario-specific; no simultaneous confidence across gates is claimed.
- Repeated-unit scenarios intentionally violate independence.
- Power is averaged across independent experiments; dependent method pairs are not counted as independent trials.
