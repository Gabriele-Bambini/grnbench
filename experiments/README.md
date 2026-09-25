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
   into it. A rerun is a new experiment ID that links the old one.
6. **Report everything.** n and what n is, the exact test, two-sided, effect size with CI, the
   correction and its family, every declared comparison including the null ones. Negative and
   inconclusive results stay in the registry. "Significant" only inside the declared family.
7. **Reproduction.** Before a result enters `paper/`, a second person (or an agent run by another
   person) reruns it from the record alone and fills section 4.

## Index

Add a row in the same PR that adds a record. Status: registered / running / done / reproduced / abandoned.

| ID | question | type | status | verdict | results |
|---|---|---|---|---|---|
| [EXP-20260917-PRE-1](EXP-20260917-PRE-1.md) | Can a BEELINE-style design rank 9 methods? | exploratory | done | inconclusive (0/108 resolved) | `benchmarks/results/beeline-2026-09-17-9methods/` |
| [EXP-20260924-PRE-1](EXP-20260924-PRE-1.md) | Engine v1 vs baseline, synthetic | exploratory | done | supports (synthetic only) | `benchmarks/results/engine-2026-09-24-v1-prefix/` |
| [EXP-20260924-PRE-2](EXP-20260924-PRE-2.md) | Engine v2, synthetic + real (partial) | exploratory | abandoned | synthetic supports, real inconclusive | `benchmarks/results/engine-2026-09-24-v2-partial/` |
| [EXP-20260925-LM-1](EXP-20260925-LM-1.md) | Motif prior lifts recovery on real data? | exploratory | done | supports | `benchmarks/results/engine-2026-09-25-models/` |

`PRE` = run before this registry existed (owner not recorded, not pre-registered).
