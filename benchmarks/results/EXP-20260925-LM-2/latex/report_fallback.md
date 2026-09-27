# EXP-20260925-LM-2: report (text fallback, no TeX engine found)

- experiment: EXP-20260925-LM-2
- type: **confirmatory** (pre-registration = the commit that adds this file, pushed before any dataset of this list was downloaded or scored)
- owner: LM   issue: board issue opened at launch   registered: 2026-09-25
- datasets (rows): 20
- models: motif_prior, pearson_motif, pearson_signed
- arms: auto, baseline
- code: branch `lm/confirmatory-motif-prior`, commit recorded in the results folder (`CODE_COMMIT`), which includes the scrublet fix (PR #3: scikit-image present, doublet fallback recorded per arm) and `real_bench --keep-artifacts` for the explainability figures.
- code commit (run folder): 67ca8809661baf5c4982a407805cccae2ebcf002
- environment: fresh venv from `requirements/lock-py311.txt` (Python 3.11.9, Windows), no torch.
- machine: local laptop (Windows, 3 workers)
- seeds: reference split seed 0 (`screen.split_reference`), cell subsample seed 0 (at most 3,000 cells per context, `MAX_SCORE_CELLS`), scrublet seed 0; regulators = the 150 TFs with most reference targets in the data (`MAX_TFS`).
- run folder: benchmarks/results/EXP-20260925-LM-2
- generated: 2026-09-27 08:34

## AUPRC ratio, auto pipeline (paired t, Holm over 4 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 20 | +0.182 | [+0.138, +0.227] | 1.4e-07 | 1.000 | pearson_motif |
| pearson_motif | motif_prior | 20 | +0.155 | [+0.077, +0.234] | 0.00083 | 1.000 | pearson_motif |

## AUPRC ratio, baseline pipeline (paired t, Holm over 4 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 20 | +0.184 | [+0.143, +0.225] | 2.9e-08 | 1.000 | pearson_motif |
| pearson_motif | motif_prior | 20 | +0.152 | [+0.076, +0.228] | 0.00073 | 1.000 | pearson_motif |

## EPR, auto pipeline (paired t, Holm over 4 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 20 | +0.120 | [-0.046, +0.286] | 0.219 | 1.000 | unresolved |
| pearson_motif | motif_prior | 20 | +0.839 | [+0.495, +1.184] | 0.00013 | 1.000 | pearson_motif |

## EPR, baseline pipeline (paired t, Holm over 4 directed hypotheses)

| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |
|---|---|---|---|---|---|---|---|
| pearson_motif | pearson_signed | 20 | +0.152 | [-0.014, +0.318] | 0.105 | 1.000 | unresolved |
| pearson_motif | motif_prior | 20 | +0.823 | [+0.491, +1.154] | 0.0001 | 1.000 | pearson_motif |
