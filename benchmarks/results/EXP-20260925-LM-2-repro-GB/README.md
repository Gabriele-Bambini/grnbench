# EXP-20260925-LM-2, second-person rerun by GB (issue #24)

The rerun of [EXP-20260925-LM-2](../../../experiments/EXP-20260925-LM-2.md) described in its
section 4b. It is not a new experiment: it uses that record, its data hashes and its code commit.
LM's results stay in `../EXP-20260925-LM-2/` (not touched).

| file | what |
|---|---|
| `rows/` | one JSON row per dataset, written by `real_bench.py run` at the code under test (`CODE_COMMIT`) on the Bocconi HPC, Slurm job 690538 |
| `REPORT_REAL.md`, `all_rows_real.csv` | `real_bench.py report --no-final-report` at the same commit, from these rows |
| `CODE_COMMIT` | commit the job ran |
| `SHA256SUMS` | written on the cluster before the fetch; `sha256sum -c SHA256SUMS` checks the four items above |
| `compare_exact.txt`, `compare_rtol1e-9.txt` | `benchmarks/compare_runs.py` (PR #43 head `2cceed0`) against `../EXP-20260925-LM-2`, with `--record experiments/EXP-20260925-LM-2.md`, exact and `--rtol 1e-9`, as fixed in #24; last line = exit code |
| `diagnose_rerun.txt` | `python benchmarks/diagnose_rerun.py benchmarks/results/EXP-20260925-LM-2/rows benchmarks/results/EXP-20260925-LM-2-repro-GB/rows` (exploratory, not part of the rule) |

Not here: `tasks.json` (it holds absolute data paths), the `--json` outputs of `compare_runs.py`,
the job log, edges and processed data. The large outputs stay on GB's cluster account in
`~/gb_runs/EXP-20260925-LM-2-repro-GB`.
