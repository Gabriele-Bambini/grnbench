# Start here: setup prompt for your coding agent

For the human (read this, then paste everything below the line into Claude Code, Codex or Antigravity):

1. Accept the GitHub invitation to `elle13eth/grnbench` first (email, or github.com/notifications).
2. Your agent will do the setup; it will stop and ask you for the few things only you can do.
3. Only you type passwords and the SSH key passphrase. Never paste them into the agent chat.
4. You need your own Bocconi HPC account (student ID) and, off campus, the Bocconi VPN.
5. Windows: install Git for Windows (Git Bash) and Python 3.11; macOS/Linux: Python 3.11.
6. Keep the repo and the data on a local disk, not in Google Drive.
7. Budget about 30 minutes; the data download (about 1.9 GB) runs on the cluster.
8. The agent ends with a checklist: read it, and fix anything marked BLOCKED.
9. Questions or access problems: ask Luca.
10. After setup, all work goes through GitHub issues and pull requests (see AGENTS.md).

---

You are setting up the `grnbench` project for a new team member on this machine. Work step by step,
run every command yourself unless it is marked PERSON, show real command output as evidence, and
stop at PERSON steps until the person confirms. Never type, store or ask for passwords or SSH
passphrases. Never change system settings or install services yourself: give the person the
command. Use bash (Git Bash on Windows) unless a step says PowerShell.

## 1. Access and clone

- PERSON: accept the invitation to `github.com/elle13eth/grnbench`.
- HTTPS (simplest): `gh auth login` (PERSON completes the browser login), then
  `gh repo clone elle13eth/grnbench ~/grnbench`. SSH alternative: `git clone git@github.com:elle13eth/grnbench.git ~/grnbench`
  if the person already has a GitHub SSH key.
- Clone to a local disk path (for example `~/grnbench`). Never into Google Drive or another synced
  folder: Drive disconnects killed long runs before, and sync tools rewrite files mid-run.
- Windows: if checkout fails on long paths, PERSON runs `git config --global core.longpaths true`.
- `.gitattributes` keeps `*.sh`, `*.sbatch` and `*.py` as LF. Do not convert them to CRLF.

## 2. Python environment (pinned)

```bash
cd ~/grnbench
python3.11 -m venv .venv            # Windows: py -3.11 -m venv .venv
. .venv/bin/activate                # Windows Git Bash: . .venv/Scripts/activate
python -m pip install -r requirements/lock-py311.txt
python -m pip install --no-deps -e .
```

Python must be 3.11 (the lock and the cluster use it). Optional RegDiffusion model: the torch
install line is at the top of `requirements/lock-py311.txt`. If you cannot install editable, prefix
every Python command with `PYTHONPATH=src`.

## 3. Sanity tests (expected result)

```bash
python -m pytest -q -p no:cacheprovider tests/test_engine.py tests/test_context.py tests/test_simulate.py tests/test_ui.py
python -m pytest -q -p no:cacheprovider
```

Expected on the full suite (Windows, 2026-09-25): `357 passed, 1 failed`; the only failure is
`tests/test_prep_gaps.py::test_drop_doublets_uses_scrublet_and_logs_what_it_removed` (needs the
scrublet backend). Another OS may differ slightly: report every failing test ID; do not "fix" code to make tests pass.

## 4. Read before any work

Read in this order and summarise each in two lines to the person: `AGENTS.md` (the rules),
`.projectos/STATE.md`, `.projectos/PURPOSE.md`, the ACTIVE entries of `.projectos/DECISIONS.md`
(D-012, D-013, D-015 set the direction: a product, not a metascience paper; honesty is a feature,
not the pitch), `experiments/README.md` (the registry), `paper/README.md`, `docs/HPC.md`.

## 5. HPC access (the person's own account)

1. PERSON: create or reuse an SSH key with a passphrase (`ssh-keygen -t ed25519 -f ~/.ssh/id_bocconi`)
   and install the public key on the cluster as Bocconi IT indicates.
2. PERSON: add a `Host bocconi-hpc` block with their own `User` and `IdentityFile` to their own
   `~/.ssh/config` (template in `docs/HPC.md` section 1). It never goes into the repo.
3. PERSON loads the passphrase into ssh-agent. Windows PowerShell: once as admin
   `Set-Service ssh-agent -StartupType Manual`, then `Start-Service ssh-agent` and
   `ssh-add $HOME\.ssh\id_bocconi`. macOS/Linux: `ssh-add ~/.ssh/id_bocconi`.
4. You: `cp .env.hpc.example .env.hpc`, set `GRNBENCH_HPC_HOST=bocconi-hpc`, confirm
   `git check-ignore .env.hpc` prints the file.
5. You: `bash hpc/hpc.sh check`. On `Could not resolve hostname` or a timeout: PERSON connects the
   Bocconi VPN (the login node resolves to a private IP). On `Permission denied`: PERSON re-runs
   `ssh-add`; try once more only after that. On `Not allowed at this time`: STOP, wait 10 minutes,
   try once. Never loop on authentication.
6. You: `bash hpc/hpc.sh sync`, then `bash hpc/hpc.sh run 'bash hpc/setup_env.sh'` (ends with `SETUP_DONE`).
7. Dry run: `bash hpc/hpc.sh submit "--partition=stud --qos=stud --time=2 --cpus-per-task=1 --mem=1G -o logs/hello_%j.out --wrap='hostname; nproc'"`,
   then `bash hpc/hpc.sh status <jobid>` until sacct says COMPLETED, and `bash hpc/hpc.sh log <jobid>`.
   Limits: 1 running + 1 queued job and 8 CPUs per user, 24 h wall time.

## 6. Data

- Only the datasets fetched by `hpc/fetch_engine_datasets.sh` are approved (references, BEELINE,
  PBMC3k, paul15, CELLxGENE). New sources need an issue first.
- On the cluster (login node, downloads only): see `docs/HPC.md` section 3.
- Locally, only if needed: `bash hpc/fetch_engine_datasets.sh ~/grnbench_data 30` (outside the repo,
  never on Google Drive). Motif models also need the cisTarget files and `GRNBENCH_CISTARGET_DIR`.
- Data is never committed. Provenance: `python benchmarks/data_manifest.py ~/grnbench_data` and
  compare with `experiments/datasets.tsv`.

## 7. Running things

- Long local runs: run from a frozen copy of `src/` + `benchmarks/` (not the working tree you are
  editing), keep the laptop plugged in; the runners hold a keep-awake request on Windows and are
  resumable (rerun the same command; finished datasets are skipped).
- Windows quirks: Git Bash heredocs mangle backslashes, so write paths with forward slashes;
  PowerShell and bash syntax differ (`$env:X` vs `$X`); the Windows ssh client is used for HPC
  because Git Bash's ssh cannot see the Windows agent (`hpc/hpc.sh` handles it).
- Known engine traps: decoupler 2.x removed `get_dorothea/get_collectri` (use
  `grnbench.tfs.decoupler_network`); BEELINE inputs are log-normalized and the engine refuses them;
  doublets are called on all cells before any subsetting; high-mito tissues (liver, heart, kidney,
  muscle) relax the 25% mito gate from context; cell-type labels with `:` or `|` go through
  `autoprep.safe_folder_name`.

## 8. How to work

- Pick or open a GitHub issue; assign yourself; branch `<initials>/<topic>` from `main`.
- Starting an experiment: copy `experiments/TEMPLATE.md` to `experiments/EXP-YYYYMMDD-<initials>-<n>.md`,
  fill registration + how-it-runs, commit and push BEFORE the run. Confirmatory runs fix hypothesis,
  endpoint, datasets, exclusions and the comparison family (Holm via `grnbench compare`) in advance.
- Results to `benchmarks/results/<EXP-ID>/` (small files) or outside git. Never write into someone
  else's EXP folder or reuse an ID; your initials keep your ID namespace separate.
- PR to `main` with the template filled, a CHANGELOG line under `## Unreleased`, tests with real
  output; one review by another person; squash-merge. No force-push.
- Never commit: secrets, passphrases, `.env.hpc`, HPC usernames or student IDs, emails, absolute
  personal paths (`C:\Users\<name>`, `/home/<id>`), host names, data.
- Handoff = the PR description (done / running where / next / do not) + `.projectos/STATE.md`.

## 9. Report back to your person

Reply with this checklist, each line DONE / BLOCKED (why) / SKIPPED (why), with the evidence:

- [ ] repo cloned at <path> (local disk, not Drive), branch `main`, `git status` clean
- [ ] Python 3.11 venv, lock installed (`python --version`, `pip check` output)
- [ ] quick tests: <summary line>; full suite: <summary line>, failing IDs = the known one only
- [ ] read AGENTS.md, STATE, PURPOSE, DECISIONS, experiments/README, paper/README, docs/HPC.md
- [ ] `~/.ssh/config` alias exists (person), key in ssh-agent (person), `.env.hpc` ignored by git
- [ ] `hpc/hpc.sh check` output (user and host redacted in your reply is fine)
- [ ] `hpc/hpc.sh sync` + `setup_env.sh` ended with SETUP_DONE
- [ ] dry-run job id, sacct State and ExitCode
- [ ] data location (cluster and/or local path outside the repo), manifest compared
- [ ] person's initials for IDs and branches: __
- [ ] anything the person must still do
