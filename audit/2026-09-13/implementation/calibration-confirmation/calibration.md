# Score-level calibration pilot

This is not biological validation or evidence of novelty.

| scenario | independent units / rows | FWER (upper bound) | coverage (lower bound) | power | gate |
|---|---:|---:|---:|---:|---|
| iid_global_null_n8 | 8 / 8 | 0.022 (0.031) | 0.958 (0.946) | n/a (no alternative) | pass |
| iid_global_null_n32 | 32 / 32 | 0.014 (0.022) | 0.963 (0.952) | n/a (no alternative) | pass |
| iid_partial_null_n16 | 16 / 16 | 0.005 (0.010) | 0.963 (0.952) | 0.088 | pass |
| iid_partial_null_n32 | 32 / 32 | 0.006 (0.012) | 0.951 (0.938) | 0.238 | pass |
| iid_margin_boundary_n32 | 32 / 32 | 0.025 (0.035) | 0.949 (0.936) | n/a (no alternative) | pass |
| pseudoreplication_null_8x4 | 8 / 32 | 0.484 (0.511) | 0.390 (0.364) | n/a (no alternative) | fail |
| pseudoreplication_partial_8x4 | 8 / 32 | 0.301 (0.326) | 0.363 (0.338) | 0.402 | fail |

Bounds: separate one-sided exact 95% Monte Carlo limits.
Gate: error upper <= 0.075 and coverage lower >= 0.925.
Finite precision may yield an inconclusive gate. A failure in repeated-unit controls is expected.

- No real expression matrix, GRN simulator, or incomplete-label model is tested.
- Means are known exactly in this score-level model; this does not hold for real algorithm rankings.
- Each gate is scenario-specific; no simultaneous confidence across gates is claimed.
- Repeated-unit scenarios intentionally violate independence.
- Power is averaged across independent experiments; dependent method pairs are not counted as independent trials.
