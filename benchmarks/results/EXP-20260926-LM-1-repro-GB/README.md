# EXP-20260926-LM-1, second-person rerun by GB (issue #93)

The rerun of [EXP-20260926-LM-1](../../../experiments/EXP-20260926-LM-1.md) described in its
section 4b. It is not a new experiment: it uses that record, its data hashes and its code commit.
LM's results stay in `../EXP-20260926-LM-1/` (not touched).

| file | what |
|---|---|
| `rows/` | one JSON row per dataset, written by `real_bench.py run` at the code under test (`CODE_COMMIT`) on the Bocconi HPC, Slurm job 693267 |
| `REPORT_REAL.md`, `all_rows_real.csv` | `real_bench.py report --no-final-report` at the same commit, from these rows |
| `CODE_COMMIT` | commit the job ran |
| `SHA256SUMS` | written on the cluster before the fetch; `sha256sum -c SHA256SUMS` checks the four items above |
| `compare_exact.txt`, `compare_rtol1e-9.txt` | `benchmarks/compare_runs.py` (PR #43 head `a39d0fd`; the file is byte-identical to the `2cceed0` used for #24) against `../EXP-20260926-LM-1`, with `--record experiments/EXP-20260926-LM-1.md`, exact and `--rtol 1e-9`, as fixed in #93; last line = exit code |

Not here: `tasks.json` (it holds data paths), the job log, edges and processed data. The large
outputs (190 MB of artifacts) stay on GB's cluster account in `~/gb_runs/EXP-20260926-LM-1-repro-GB`.
