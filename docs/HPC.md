# HPC guide (Bocconi Slurm cluster)

Each person uses **their own** Bocconi HPC account. Nothing in the repo names a user: your host
alias and user live in your `~/.ssh/config`, your settings in the untracked `.env.hpc`.
All commands below run from the repo root in bash (Git Bash on Windows).

## 1. First-time setup (once per person, done by the person, not by an agent)

1. **Key.** If you have none: `ssh-keygen -t ed25519 -f ~/.ssh/id_bocconi` (set a passphrase).
   Install the public key on the cluster the way Bocconi IT indicates (`ssh-copy-id` if password
   login is allowed for your account).
2. **Alias** in your own `~/.ssh/config` (never in the repo):

   ```
   Host bocconi-hpc
       HostName slogin.hpc.unibocconi.it
       User <your-student-id>
       IdentityFile ~/.ssh/id_bocconi
       IdentitiesOnly yes
       ServerAliveInterval 60
   ```

3. **Passphrase into ssh-agent, by you.** Agents never type, store or ask for it.
   - Windows (PowerShell): once, as admin, `Set-Service ssh-agent -StartupType Manual`; then
     `Start-Service ssh-agent` after each reboot and `ssh-add $HOME\.ssh\id_bocconi`.
     The key stays registered; only the service needs starting.
   - macOS: `ssh-add --apple-use-keychain ~/.ssh/id_bocconi`. Linux: `eval "$(ssh-agent)"; ssh-add ~/.ssh/id_bocconi`.
4. **Settings:** `cp .env.hpc.example .env.hpc` and set `GRNBENCH_HPC_HOST=bocconi-hpc`
   (your alias). `.env.hpc` is gitignored; check with `git check-ignore .env.hpc`.
5. **Check:** `bash hpc/hpc.sh check` prints your user, the host and the Slurm version.

If it fails:

| message | meaning | do |
|---|---|---|
| `Could not resolve hostname` / timeout | off campus: the login node resolves to a private IP | connect the Bocconi VPN |
| `Permission denied (publickey)` | key not in the agent, or wrong `User` | `ssh-add`, check `~/.ssh/config`; do **not** retry in a loop |
| `Not allowed at this time` | lockout after failed logins (per-source penalty) | **stop, wait 10 minutes**, fix the cause, try once |

`hpc/hpc.sh` always runs ssh with `BatchMode=yes`, so it fails at once instead of waiting for a
passphrase. On Windows it uses `C:\Windows\System32\OpenSSH\ssh.exe`, because Git Bash's own ssh
cannot reach the Windows ssh-agent.

## 2. Code to the cluster

```bash
bash hpc/hpc.sh sync                          # tar over ssh into ~/$GRNBENCH_HPC_DIR, writes GIT_COMMIT
bash hpc/hpc.sh run 'bash hpc/setup_env.sh'   # once, and after the lock file changes
```

`sync` copies your working tree (uncommitted changes included, marked `-dirty` in `GIT_COMMIT`).
Commit before a registered run, so the record can cite a clean hash. `sync` does not delete remote
files you deleted locally. Everyone syncs into their own home: no shared remote checkout.

## 3. Data on the login node

Compute nodes have no outbound network: fetch on the login node (downloads only, no compute).

```bash
bash hpc/hpc.sh run 'nohup bash hpc/fetch_engine_datasets.sh ~/grnbench_data 30 > logs/fetch.log 2>&1 &'
bash hpc/hpc.sh run 'tail -5 logs/fetch.log'
```

About 1.9 GB (references, BEELINE, PBMC3k, paul15, 30 CELLxGENE datasets); `INDEX.md` lists them.
Motif models also need the cisTarget files (about 1.3 GB, resources.aertslab.org/cistarget) in a
folder named by `GRNBENCH_CISTARGET_DIR`. Record provenance for the experiment:

```bash
A='source /software/miniconda3/etc/profile.d/conda.sh && conda activate grnbench'
bash hpc/hpc.sh run "$A && python benchmarks/data_manifest.py ~/grnbench_data" > /tmp/datasets-hpc.tsv
```

## 4. Submit

Measured QOS limits (`sacctmgr show qos stud`, 2026-08-22; re-check with
`bash hpc/hpc.sh run "sacctmgr show qos stud format=Name,MaxWall,MaxTRESPU,MaxJobsPU,MaxSubmitPU"`):

| limit | value | consequence |
|---|---|---|
| partition / QOS | `stud` / `stud`, always pass `--qos=stud` (no default QOS) | scripts already do |
| MaxWall | 24 h | runners are resumable: resubmit after a TIMEOUT |
| MaxTRESPU | 8 CPUs per user | ask for `--cpus-per-task=8` at most |
| MaxJobsPU / MaxSubmitPU | 1 running, 2 submitted per user | **one job that runs many tasks**, not a big array |

Slurm counts array elements toward these per-user limits, so the `--array=0-N%20` pattern from the
old handoff is refused or runs one element at a time (not yet measured here; record what you see).
Run everything in one job instead (`CHUNK` larger than the task count, index 0):

```bash
bash hpc/hpc.sh run "$A && OUT=~/gb_runs/<EXP-ID> && mkdir -p \$OUT && \
  python benchmarks/engine_bench.py plan --presets-reps 10 --random 300 --out \$OUT/tasks.json"
bash hpc/hpc.sh submit "--time=24:00:00 --export=ALL,OUT=\$HOME/gb_runs/<EXP-ID>,CHUNK=100000 hpc/engine_bench.sbatch"

# real data
bash hpc/hpc.sh run "$A && python benchmarks/real_bench.py plan --data-root ~/grnbench_data --out ~/gb_runs/<EXP-ID>/real_tasks.json"
bash hpc/hpc.sh submit "--time=24:00:00 --export=ALL,OUT=\$HOME/gb_runs/<EXP-ID>,SCRIPT=benchmarks/real_bench.py,TASKS=\$HOME/gb_runs/<EXP-ID>/real_tasks.json,REFS=\$HOME/grnbench_data/references,CHUNK=100000 hpc/engine_bench.sbatch"
```

Cost: about 5 min per synthetic dataset and 10 min per real one per core; with 8 parallel
workers ~490 synthetic datasets fit in one 24 h job. Use `<EXP-ID>` from your record as the run
folder, so runs from different people never share a folder.

First time, try a 1-minute job: `bash hpc/hpc.sh submit "--partition=stud --qos=stud --time=2 --cpus-per-task=1 --mem=1G -o logs/hello_%j.out --wrap='hostname; nproc'"`,
then `bash hpc/hpc.sh status <jobid>` and `bash hpc/hpc.sh log <jobid>`.

## 5. Monitor

```bash
bash hpc/hpc.sh status            # squeue for you
bash hpc/hpc.sh status <jobid>    # plus sacct: State and ExitCode
bash hpc/hpc.sh log <jobid> 100   # tail of logs/*_<jobid>*.out
bash hpc/hpc.sh cancel <jobid>
```

A job that left `squeue` has not necessarily succeeded: read `sacct` State/ExitCode and the job's
own output. Slurm snapshots the job script at submit time: editing the `.sbatch` does not change
a pending job (cancel and resubmit).

## 6. Results back

```bash
bash hpc/hpc.sh run "$A && python benchmarks/engine_bench.py report --out ~/gb_runs/<EXP-ID>"
bash hpc/hpc.sh run "$A && python benchmarks/real_bench.py report --out ~/gb_runs/<EXP-ID>"
bash hpc/hpc.sh fetch '~/gb_runs/<EXP-ID>' /path/outside/the/repo
```

Copy only the small files (REPORT*.md, all_rows*.csv, rows/*.json, tasks.json) into
`benchmarks/results/<EXP-ID>/` and fill section 3 of the record.

## 7. Hardware facts (measured 2026-07/08; they matter only for GPU work, grnbench is CPU-only)

- Nodes `gnode01` (32 CPU, 990 GB) and `gnode02` (24 CPU, 500 GB). GRES per node
  `gpu:4g.40gb:4,gpu:3g.40gb:4`: 4 physical A100 80GB cards, each split into one `4g.40gb` and one
  `3g.40gb` MIG slice. QOS `stud` caps a job at 4 GPU slices.
- `--gres=gpu:4` gives one slice on each of four different cards. No peer-to-peer or NVLink;
  NCCL works over PCIe at about 14-16 GB/s. GPU0+1 and GPU2+3 sit on different NUMA nodes.
- A `4g.40gb` slice: 56 SMs, about 42 GB usable, about 150 TFLOPS bf16 after warm-up. A `3g.40gb`
  slice exposes the same memory with 42 SMs.
- Slurm gives ordinal `CUDA_VISIBLE_DEVICES` (0,1,...), not MIG UUIDs.
- torch spawn workers: guard scripts with `if __name__ == "__main__":` and never feed a script on
  stdin (spawned workers cannot re-import it).
- `module load cuda/12.4` gives nvcc but does not set `CUDA_HOME`. No containers on the cluster.
- Home is on BeeGFS; large scratch at `/mnt/beegfsstudents`.
