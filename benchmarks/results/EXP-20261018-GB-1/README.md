# EXP-20261018-GB-1 results

Run folder of `benchmarks/aracne3_smoke.py`: smoke job 705969 on the Bocconi HPC, code `6927b9f`,
against the ARACNe3 program that build job 705955 compiled. Record: `experiments/EXP-20261018-GB-1.md`.

- `SUMMARY.md`, `summary.json`: written by the script. Both were regenerated after the run with
  `python benchmarks/aracne3_smoke.py --resummarize --out benchmarks/results/EXP-20261018-GB-1`
  (code `5d96311`): the run's parser did not read ARACNe3's whole-second log times (fixed in
  `aa16d71`), and the statistics of the consolidated columns were added (`5d96311`). Only
  `log_times_s` and `consolidated.columns` changed (record, section 3, deviations).
- `SHA256SUMS.cluster`: the checksums written on the cluster right after the run. Every file still
  matches them except the two regenerated summaries. `SHA256SUMS`: the folder as committed.
- `CODE_COMMIT`: the `GIT_COMMIT` of the cluster checkout the job ran.
- `subset_cells.txt`, `subset_genes.txt`, `regulators.txt`: the fixed subset.
- `runs/<label>/`: `edges.tsv` (what `run_method` returned), `messages.txt` (ARACNe3's exit status,
  stdout and stderr, and the adapter's answer) and `aracne3_out/` (ARACNe3's own files: the
  consolidated network, the subnetwork, its logs).
- `logs/`: the Slurm logs of the build job and of the smoke job.

Home paths are written as `~` and the adapter's temporary folders as `<tmp>`.
