# EXP-YYYYMMDD-XX-n: <one-line question>

<!-- Copy to experiments/EXP-YYYYMMDD-<initials>-<n>.md. Commit sections 1-2 BEFORE the run.
     Never edit section 1 after the run started; add a dated "Deviations" line instead. -->

## 1. Registration (before the run)

- type: confirmatory | exploratory
- owner: <initials>   issue: #<n>   registered: YYYY-MM-DD (commit of this file = the pre-registration)
- hypothesis: <one sentence, directional if confirmatory>
- primary endpoint: <metric, e.g. AUPRC ratio on the CollecTRI evaluation half>; secondary: <...>
- unit of analysis: dataset (n = <k>). Never edges, cells or genes.
- datasets: <paths from experiments/datasets.tsv, or the tasks.json that lists them>
- exclusion rules (decided now): <e.g. targeted panels < 1,000 genes; datasets refused as not counts>
- comparison family: <every pair that will be tested, e.g. pearson_motif vs pearson_signed, ...>
- test: paired, two-sided, Holm over the family above, alpha 0.05, via `grnbench compare`
  (`grnbench.eval.comparison.compare_units`). Effect size = mean paired difference with CI.
- primary comparison: <A> > <B> | <metric, e.g. auprc_ratio> | <arm, e.g. auto>   (machine-readable H1)
- smallest effect that matters (margin): <value or "none declared">
- seeds: <list>; reference split seed: <n>

## 2. How it runs

- code: git commit <hash> (clean tree; for HPC the `GIT_COMMIT` file written by `hpc/hpc.sh sync`)
- environment: requirements/lock-py311.txt, Python <x.y.z>, extra packages: <none | list+versions>
- command(s): <exact CLI, copy-pasteable, paths relative or ~/...>
- machine: <local OS / HPC job id(s), partition, CPUs, memory>
- results: `benchmarks/results/<EXP-ID>/` (small) and <large outputs location, outside git>
- report input: `benchmarks/results/<EXP-ID>` (folder with rows/*.json, read by `grnbench report`)

## 3. Results (filled by scripts, not by hand)

- status: running | done | failed | abandoned
- produced by: <command that wrote the report>, report file: <path>
- headline: <metric> = <value> [95% CI], n = <k> datasets, <test>, p_Holm = <value>
- all declared comparisons (including the null ones): <table from the report, or link>
- verdict: supports | refutes | inconclusive  (inconclusive is not "equal")
- figures (exploration + explainability, required): `benchmarks/results/<EXP-ID>/viz/index.md`
  (from `--keep-artifacts` + `grnbench viz`); list what was not produced and why
- generated report: `benchmarks/results/<EXP-ID>/latex/report.tex` (`grnbench report <EXP-ID>`)
- deviations from section 1: <none | dated list>

## 4. Reproduction (required before anything reaches paper/)

- reproduced by: <initials>, YYYY-MM-DD, commit <hash>, machine <...>
- directed and checked by: <initials> (signs this section); commands typed by: <the same person | an
  agent session run for <initials>>
- data: fetched or copied by <initials>, `data_manifest.py --check` output (D-020)
- match: <same numbers | differences and why>

## 5. Next

- what changes next, and which issue/experiment follows: <...>
