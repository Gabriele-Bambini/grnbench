# EXP-20261003-LM-1, second-person rerun by GB (issue #95)

The rerun of [EXP-20261003-LM-1](../../../experiments/EXP-20261003-LM-1.md) described in its
section 4b. It is not a new experiment: it uses that record, its data hashes and its code commit.
LM's results stay in `../EXP-20261003-LM-1/` (not touched).

| file | what |
|---|---|
| `rows/` | one JSON row per task (20 datasets x 3 conditions), written by `candidate_study.py run` at the code under test (`CODE_COMMIT`) on the Bocconi HPC, Slurm job 702612 |
| `SUMMARY.md`, `candidates.csv` | `candidate_study.py report` at the same commit, from these rows |
| `CONFIRM.md`, `CONFIRM.json` | `candidate_study.py confirm` at the same commit, from these rows |
| `CODE_COMMIT` | commit the job ran |
| `SHA256SUMS` | written on the cluster before the fetch; `sha256sum -c SHA256SUMS` checks the items above |
| `compare_rows_exact.txt`, `compare_rows_rtol1e-9.txt` | check (1) of #95: `benchmarks/compare_runs.py --rows-only` (PR #43 head `a39d0fd`, byte-identical to the `2cceed0` used for #24) against `../EXP-20261003-LM-1`, exact and `--rtol 1e-9`; last line = exit code |
| `compare_confirm_exact.txt`, `compare_confirm_rtol1e-9.txt` | check (2) of #95: `benchmarks/compare_confirm.py ../EXP-20261003-LM-1/CONFIRM.json CONFIRM.json --record experiments/EXP-20261003-LM-1.md`, exact and `--rtol 1e-9`; last line = exit code |
| `diff_rows_nested.txt` | diagnosis, not part of the rule: `python benchmarks/diff_rows_nested.py benchmarks/results/EXP-20261003-LM-1/rows benchmarks/results/EXP-20261003-LM-1-repro-GB/rows` (where the rows differ at any depth, timing fields apart) |

`candidates.csv` differs from the committed one only in its timing columns (`screen_seconds`,
`eval_seconds`, `seconds`): checked with pandas, `a.drop(columns=t).equals(b.drop(columns=t))` is
`True` for the 30 other columns of the 1,008 rows.

Not here: `tasks.json` (it holds data paths), the job log, the per-task work folders.
