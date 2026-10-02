---
name: grnbench-hpc
description: Run grnbench work on the Bocconi HPC Slurm cluster - connect, sync code, fetch data, submit, wait for, monitor, fetch and record jobs, with the per-user .env.hpc and hpc/hpc.sh. Use for "run on HPC", "submit to the cluster", "check my job", "fetch results from the cluster", "Not allowed at this time", "sbatch", "squeue", "MIG".
---

# grnbench on the Bocconi HPC

Everything goes through `hpc/hpc.sh` from the repo root (bash; Git Bash on Windows). Your host
alias comes from your own `~/.ssh/config`, your settings from the untracked `.env.hpc`. The full
human guide is `docs/HPC.md`; this skill is the operating procedure for agents.

## Hard rules

- Never type, store, echo or ask for a passphrase or password. The person loads the key into
  ssh-agent. `hpc.sh` runs ssh with `BatchMode=yes`, so a missing key fails at once.
- `Permission denied (publickey)`: stop and ask the person to run `ssh-add` (and check `User` in
  their `~/.ssh/config`). Try again only once, after they confirm.
- `Not allowed at this time`: lockout after failed logins. **Stop, wait 10 minutes, fix the cause,
  try once.** Never loop on authentication.
- `Could not resolve hostname` or a timeout: the VPN is down. **Reconnect it yourself first**
  (`hpc.sh` prints the command for the OS as the first error line). Windows:
  `powershell -NoProfile -ExecutionPolicy Bypass -File hpc/vpn.ps1 reconnect`. Linux:
  `globalprotect connect --portal vpnstudents.unibocconi.it`. Then retry once. Ask the person only when:
  `reconnect` exits 1 (90 s without a tunnel: the portal wants their SSO sign-in in the
  GlobalProtect window), exits 2 (the one-time rights grant is missing: they run
  `hpc/vpn.ps1 setup` and the admin command it prints, `docs/HPC.md` section 1b), or on macOS
  (no GlobalProtect CLI: they connect in the app).
- No compute on the login node: downloads, `sbatch`, `squeue`, `sacct`, small file operations only.
- Never write a username, student ID, home path or host name into the repo, a record or a PR.
  Write `~/...` for cluster paths.
- Before launching anything, run the board (`python -m grnbench.board`) and do not relaunch a
  run another person's agent has in flight (open `running` issue).

## 1. Connect

```bash
test -f .env.hpc || cp .env.hpc.example .env.hpc      # then the person sets GRNBENCH_HPC_HOST
bash hpc/hpc.sh check                                 # user@host, Slurm version, remote dir
```

## 2. Sync code (every time before a run)

```bash
git status --short            # commit first: a registered run cites a clean commit
bash hpc/hpc.sh sync          # tar over ssh to ~/$GRNBENCH_HPC_DIR, writes GIT_COMMIT
bash hpc/hpc.sh run 'bash hpc/setup_env.sh'   # first time, or when requirements/lock-py311.txt changed
```

## 3. Data (login node)

```bash
bash hpc/hpc.sh run 'mkdir -p logs && nohup bash hpc/fetch_engine_datasets.sh ~/grnbench_data 30 > logs/fetch.log 2>&1 &'
bash hpc/hpc.sh run 'tail -3 logs/fetch.log'          # poll; ends with "done: ..."
```

Only the approved datasets of that script. Anything else needs an issue first.

## 4. Submit

Limits of QOS `stud` (measured): 1 running + 1 queued job per user, 8 CPUs per user, 24 h wall
time, `--qos=stud` must be explicit. So one job runs all tasks of an experiment:

```bash
A='source /software/miniconda3/etc/profile.d/conda.sh && conda activate grnbench'
EXP=EXP-YYYYMMDD-XX-n                                  # the record must already be committed
bash hpc/hpc.sh run "$A && mkdir -p ~/gb_runs/$EXP && python benchmarks/real_bench.py plan --data-root ~/grnbench_data --out ~/gb_runs/$EXP/real_tasks.json"
bash hpc/hpc.sh submit "--time=24:00:00 --export=ALL,OUT=\$HOME/gb_runs/$EXP,SCRIPT=benchmarks/real_bench.py,TASKS=\$HOME/gb_runs/$EXP/real_tasks.json,REFS=\$HOME/grnbench_data/references,CHUNK=100000,MAX_ATTEMPTS=2 hpc/engine_bench.sbatch"
```

`submit` prints the job id. Then open the board issue (AGENTS.md "Live board"):
`gh issue create --label running --title "$EXP <short name>" --body "<who, HPC job id, started, expected finish, output ~/gb_runs/$EXP, check: bash hpc/hpc.sh status <id>>"`.

First time only, a 1-minute test job:
`bash hpc/hpc.sh submit "--partition=stud --qos=stud --time=2 --cpus-per-task=1 --mem=1G -o logs/hello_%j.out --wrap='hostname; nproc'"`.

## 5. Wait correctly

Poll the real condition with a growing interval; never a blind long sleep:

```bash
id=<jobid>; d=60
while :; do
  out=$(bash hpc/hpc.sh status $id) || { echo "ssh failed: stop and report"; break; }
  echo "$out" | grep -Eq 'COMPLETED|FAILED|CANCELLED|TIMEOUT|OUT_OF_ME|NODE_FAIL' && { echo "$out"; break; }
  echo "$out" | grep -q 'QOSMax' && { echo "$out"; echo "blocked by QOS limit: stop and report"; break; }
  sleep $d; d=$(( d < 900 ? d * 2 : 900 ))
done
```

- `PENDING` with reason `QOSMaxJobsPerUserLimit` / `QOSMaxSubmitJobPerUserLimit`: you already
  have a job running; stop and report, do not submit more.
- `PENDING (Resources)` / `(Priority)`: normal; the queue can take hours. Report the ETA
  (`squeue --start`) and comment it on the board issue instead of waiting in the session.
- A job that left `squeue` is not a success: `sacct` State and ExitCode decide, then the job's own output.
- Slurm snapshots the script at submit time: an edited `.sbatch` needs `cancel` + resubmit.

## 6. Monitor, cancel

```bash
bash hpc/hpc.sh status            # your queue
bash hpc/hpc.sh status <id>       # + sacct State, Elapsed, ExitCode
bash hpc/hpc.sh log <id> 100      # tail logs/*_<id>*.out
bash hpc/hpc.sh cancel <id>
```

## 7. Finished: fetch and record

```bash
bash hpc/hpc.sh run "$A && python benchmarks/real_bench.py report --out ~/gb_runs/$EXP"
bash hpc/hpc.sh run "cd ~/gb_runs/$EXP && ls rows | wc -l && find . -type f -name '*.json' -print0 | xargs -0 sha256sum > SHA256SUMS"
bash hpc/hpc.sh fetch "~/gb_runs/$EXP" <local folder outside the repo>
(cd <local folder>/$EXP && sha256sum -c SHA256SUMS | grep -v ': OK$' || echo all OK)
```

The report validates the task-manifest hash and refuses missing or retryable rows; only
`--partial` writes a watermarked incomplete report. `COMPLETE.json` must exist, the row count must
equal the task count and the report must show 0 unexpected errors. Copy the small
files to `benchmarks/results/$EXP/`, run `grnbench viz` and `grnbench report` on it, fill section 3
of the record from the generated report (never type numbers by hand), close the board issue with
label `done` or `failed` and a one-line result.

## 8. Failures and resubmits

- TIMEOUT or node failure: resubmit the exact same command. The runners skip terminal rows,
  retry `transient_error` rows and missing tasks, and turn a repeated row-level failure into
  `permanent_error` at `MAX_ATTEMPTS` (default 2). Never delete a row merely to make a report pass.
- A dataset erroring every time: it is a finding, record it; do not drop it silently.
- Out of memory: raise `--mem` (per job) rather than the CPU count; 8 CPUs is the cap.

## 9. GPU and MIG facts (grnbench is CPU-only; for learned models)

- Nodes: 32 CPU / 990 GB and 24 CPU / 500 GB; GRES `gpu:4g.40gb:4,gpu:3g.40gb:4` = 4 A100 80GB
  cards per node, each split into one `4g.40gb` and one `3g.40gb` MIG slice. A job gets at most 4 slices.
- `--gres=gpu:4` gives one slice on each of four cards: no P2P or NVLink, NCCL over PCIe ~14-16
  GB/s, GPU0+1 and GPU2+3 on different NUMA nodes (tensor parallel 2 inside a pair is cheapest).
- `4g.40gb`: 56 SMs, ~42 GB usable, ~150 TFLOPS bf16 after warm-up (unwarmed probes read 5x lower).
  `3g.40gb`: same memory, 42 SMs.
- `CUDA_VISIBLE_DEVICES` is ordinal (0,1,...), not MIG UUIDs.
- torch spawn: guard entry points with `if __name__ == "__main__":`; never feed a script on stdin.
- `module load cuda/12.4` gives nvcc but not `CUDA_HOME`; no containers on the cluster.
- Home is on BeeGFS; large scratch at `/mnt/beegfsstudents`.
