# Experiment registry

One file per run, local or HPC. This is the only place runs are recorded (not
`.projectos/ATTEMPTS.md`, not chat, not Drive).

## Rules

1. **Record first.** Copy `TEMPLATE.md` to `EXP-YYYYMMDD-<initials>-<n>.md` (n counts your own runs
   that day, from 1). Fill sections 1-2, commit and push them on your branch **before** the run.
   That commit is the pre-registration. IDs are never reused or renamed.
2. **Confirmatory or exploratory.** Confirmatory = hypothesis, primary endpoint, unit (dataset),
   dataset list, exclusion rules and the full comparison family are fixed before the run, and the
   test is paired over datasets with Holm (`grnbench compare`, i.e. `compare_units`).
   Everything else is exploratory. Only confirmatory results can back a headline claim.
3. **Code and data are pinned.** Commit hash (clean tree, or `GIT_COMMIT` from `hpc/hpc.sh sync`),
   `requirements/lock-py311.txt` plus any extras with versions, seeds, and datasets cited by their
   `path` in `datasets.tsv` (source, version, sha256). New data: regenerate `datasets.tsv` with
   `python benchmarks/data_manifest.py DATA_ROOT > experiments/datasets.tsv` in the same PR.
4. **Results by script only.** Reports are written by the bench scripts; section 3 quotes them.
   Small outputs (reports, CSV/JSON summaries; keep a folder under ~5 MB) go to
   `benchmarks/results/<EXP-ID>/`; large outputs stay outside git, path in the record.
5. **Your folder only.** Never write into another person's `benchmarks/results/<EXP-ID>/` or rerun
   into it. A rerun is a new experiment ID that links the old one. A reproduction (rule 8) is not a new
   experiment: it goes in section 4 of the original record, with its small outputs in
   `benchmarks/results/<EXP-ID>-repro-<initials>/`.
6. **Report everything.** n and what n is, the exact test (paired t: two-sided interval, one-sided p in
   each direction, each direction its own hypothesis), effect size with CI, the correction and its
   family, every declared comparison including the null ones. Negative and inconclusive results stay
   in the registry. "Significant" only inside the declared family.
7. **Look inside.** Each record links its figures (`viz/index.md`: QC, embeddings, decisions,
   PR curves, networks, TF-level errors, motif lift, calibration, stability, attributions) and its
   generated report (`latex/`). Scores alone are not a result.
8. **Reproduction.** Before a result enters `paper/`, a second person reruns it from the record
   alone and fills section 4. They direct and check every step on their own account, with data they fetched or copied themselves and checked against the registered sha256 (D-020); an agent may
   type the commands. Reruns directed by the owner's side are `reproduced-computational`.

## Index

Add a row in the same PR that adds a record. Status: registered / running / done / reproduced-computational (rerun directed by the owner's side, second-person reproduction pending) / reproduced / abandoned.

| ID | question | type | status | verdict | results |
|---|---|---|---|---|---|
| [EXP-20260917-PRE-1](EXP-20260917-PRE-1.md) | Can a BEELINE-style design rank 9 methods? | exploratory | done | inconclusive (0/108 resolved) | `benchmarks/results/beeline-2026-09-17-9methods/` |
| [EXP-20260924-PRE-1](EXP-20260924-PRE-1.md) | Engine v1 vs baseline, synthetic | exploratory | done | supports (synthetic only) | `benchmarks/results/engine-2026-09-24-v1-prefix/` |
| [EXP-20260924-PRE-2](EXP-20260924-PRE-2.md) | Engine v2, synthetic + real (partial) | exploratory | abandoned | synthetic supports, real inconclusive | `benchmarks/results/engine-2026-09-24-v2-partial/` |
| [EXP-20260925-LM-1](EXP-20260925-LM-1.md) | Motif prior lifts recovery on real data? | exploratory | done | supports | `benchmarks/results/engine-2026-09-25-models/` |
| [EXP-20260925-LM-2](EXP-20260925-LM-2.md) | Confirmatory: motif prior lifts recovery on 20 fresh real datasets? | confirmatory | reproduced-computational | supports (+0.18, n = 20, p_Holm 1e-7) | `benchmarks/results/EXP-20260925-LM-2/` |
| [EXP-20260926-LM-1](EXP-20260926-LM-1.md) | Confirmatory: does the motif prior lift RegDiffusion and GRNBoost2 (LM-2's 20 datasets, HPC)? | confirmatory | reproduced | supports both (RegDiffusion +0.18, n = 19, p_Holm 1e-6; GRNBoost2 +0.13, n = 20, p_Holm 4e-10); second-person rerun: every non-timing field identical (section 4b, signed, #93) | `benchmarks/results/EXP-20260926-LM-1/` |
| [EXP-20261001-LD-1](EXP-20261001-LD-1.md) | Context-specific stress curves and complete approved real-data validation | exploratory | done | 740/740 synthetic and 40/40 real accounted; auto ties baseline on most paired real data | `benchmarks/results/EXP-20261001-LD-1/` |
| [EXP-20261002-LM-1](EXP-20261002-LM-1.md) | Confirmatory: does the motif prior lift Pearson, GRNBoost2 and RegDiffusion on 20 fresh datasets? | confirmatory | reproduced | supports all three (Pearson +0.20, RegDiffusion +0.23, GRNBoost2 +0.15; n = 20; p_Holm <= 1.6e-6); second-person rerun: primaries identical, one warning counter differs (section 4b, signed, #94) | `benchmarks/results/EXP-20261002-LM-1/` |
| [EXP-20261002-LM-2](EXP-20261002-LM-2.md) | Exploratory: where does cleaning help on real data under controlled degradations, and which selection rule finds it? | exploratory | done | little headroom; without a reference the engine loses to the baseline in all 7 conditions; Pearson residuals, batch correction, metacells, per-type contexts lower recovery | `benchmarks/results/EXP-20261002-LM-2/` |
| [EXP-20261003-LM-1](EXP-20261003-LM-1.md) | Confirmatory: baseline vs the reference-free choice; gene-set-preserving screening vs baseline (20 unused datasets) | confirmatory | done | H1 supported (baseline beats the reference-free choice by +0.089, p_Holm 0.044; also under doublets and damaged cells); H2 inconclusive; rerun directed by GB, typed by an agent: primaries identical, row check DIFFERENT on timings only, rule not met (section 4b, #95) | `benchmarks/results/EXP-20261003-LM-1/` |
| [EXP-20261003-GB-1](EXP-20261003-GB-1.md) | Exploratory: how much of the motif-prior lift for Pearson (C6 primary, C4 secondary) needs TF-specific target identity, beyond TF-level, database-membership and gene-propensity effects? (stage 1 of the 2026-09-27 audit's decomposition) | exploratory | done | C6 (primary) inconclusive, C4 mixed; combined inconclusive (n = 20 + 20) | `benchmarks/results/EXP-20261003-GB-1/` |
| [EXP-20261004-LM-1](EXP-20261004-LM-1.md) | Exploratory: on the genes a candidate shares with the baseline, does it still beat it? (HVG size, context, metacells) | exploratory | done | on shared genes no cleaning change beats the baseline by more than ~0.01; the hvg=2000 and selection gains were the larger edge universe | `benchmarks/results/EXP-20261004-LM-1/` |
| [EXP-20261004-LM-2](EXP-20261004-LM-2.md) | Exploratory: from a raw droplet matrix, does the engine reach the recovery of the vendor-filtered matrix? (8 public 10x datasets) | exploratory | done | from raw droplets the engine reaches vendor-filtered recovery on 7/8 (heart_10k failed: knee floor bug, fixed) | `benchmarks/results/EXP-20261004-LM-2/` |
| [EXP-20261010-LD-1](EXP-20261010-LD-1.md) | Exploratory product calibration: do D-027 identifier guards avoid false refusals on approved real data? | exploratory | failed | manifest preflight invalid (0 inputs inspected) | `benchmarks/results/EXP-20261010-LD-1/` |
| [EXP-20261010-LD-2](EXP-20261010-LD-2.md) | Exact-snapshot rerun of the D-027 approved-data refusal calibration | exploratory | done | fails shipping gate: 32/33 accepted, one zebrafish false refusal | `benchmarks/results/EXP-20261010-LD-2/` |
| [EXP-20261011-LD-1](EXP-20261011-LD-1.md) | Does D-028 remove the zebrafish false refusal without weakening human/mouse guards? | exploratory | done | fails shipping gate: 32/33 accepted; non-human Ensembl evidence was not recognised | `benchmarks/results/EXP-20261011-LD-1/` |
| [EXP-20261011-LD-2](EXP-20261011-LD-2.md) | Does non-human Ensembl evidence complete D-028 on approved inputs? | exploratory | done | passes product gate: 33/33 accepted, 0 refusals/errors; zebrafish preserved as other | `benchmarks/results/EXP-20261011-LD-2/` |
| [EXP-20261016-LM-1](EXP-20261016-LM-1.md) | Confirmatory: do networks inferred from control cells predict CRISPRi knockdown responders beyond gene propensity, and does the motif prior help? (RPE1 design, K562 replication, TF as the unit) | confirmatory | registered | pending | `benchmarks/results/EXP-20261016-LM-1/` |
| [EXP-20261018-LD-1](EXP-20261018-LD-1.md) | Does grnbench cisTarget pruning reproduce official pySCENIC 0.12.1 on identical inputs? | exploratory | done | primary endpoint not met: final motif-supported pairs diverge (Jaccard 0.43-0.58); secondary union of activating-module pairs matches exactly | `benchmarks/results/EXP-20261018-LD-1/` |
| [EXP-20261018-GB-1](EXP-20261018-GB-1.md) | Exploratory: does the `aracne3` adapter work against the real ARACNe3 program (first HPC build), and do its runs repeat with the same seed? | exploratory | done | supports: the adapter works as written (exit 0, 478 edges); same-seed runs repeat exactly with 1 thread; with 4 threads the null model differs, the network did not | `benchmarks/results/EXP-20261018-GB-1/` |
| [EXP-20261020-GB-1](EXP-20261020-GB-1.md) | Exploratory: why is C4's EPR lift of the motif prior not resolved, and how many never-used datasets would a new study need? (committed rows only) | exploratory | done | underpowered at n = 20, no heavy tails; leaving one dataset out changes the verdict in 0/20 (auto) and 2/20 (baseline) drops; +0.10 EPR needs 101 (one pair) to 123 (family of 4) datasets for 80% power; C4 unchanged | `benchmarks/results/EXP-20261020-GB-1/` |
| [EXP-20261022-GB-1](EXP-20261022-GB-1.md) | Confirmatory: does the motif prior lift Pearson's EPR on the 120 confirmatory datasets of the D-057 CELLxGENE pool? | confirmatory | done | H1 supported: EPR +0.149 [+0.067, +0.230], n = 119, p_Holm 0.000424; T3 not run; 1 dataset missing (rule e, #169 error); report `--partial`, deviations await GB | `benchmarks/results/EXP-20261022-GB-1/` |
| [EXP-20261022-GB-2](EXP-20261022-GB-2.md) | Exploratory stress test: where does cleaning beat the baseline under doublets and damaged cells, on the 80 exploratory datasets of the D-057 pool? | exploratory | registered | pending | `benchmarks/results/EXP-20261022-GB-2/` |

`PRE` = run before this registry existed (owner not recorded, not pre-registered).
