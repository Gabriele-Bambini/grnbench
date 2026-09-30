# Start here: setup prompt for your coding agent

For the human (read this, then paste everything below the line into Claude Code, Codex or Antigravity):

1. Accept the GitHub invitation to `elle13eth/grnbench` first (email, or github.com/notifications).
2. Your agent does the setup and stops at PERSON steps: only you can do those.
3. Only you type passwords and the SSH key passphrase. Never paste them into the agent chat.
4. You need your own Bocconi HPC account (student ID) and, off campus, the Bocconi VPN.
5. An empty laptop is fine: part A installs everything (Git, Python 3.11, gh, optional TeX).
6. Keep the repo and the data on a local disk, never in Google Drive.
7. Budget about one hour the first time; the data download (about 1.9 GB) runs on the cluster.
8. The agent ends with a checklist: read it, and fix anything marked BLOCKED.
9. Questions or access problems: ask Luca.
10. After setup, all work goes through GitHub issues and pull requests (part B and AGENTS.md).

---

You are setting up the `grnbench` project for a team member on this machine and then working on
it with them. Work step by step, run every command yourself unless it is marked PERSON, show
real command output as evidence, and stop at PERSON steps until the person confirms. Never type,
store or ask for passwords or SSH passphrases. Never change system settings or install services
yourself: give the person the exact command. Use bash (Git Bash on Windows) unless a step says
PowerShell.

# Part A. From zero to a first registered experiment

## A1. Tools (PERSON installs; you verify with the version commands)

Windows (PowerShell, as the person):
```powershell
winget install --id Git.Git -e            # Git + Git Bash
winget install --id Python.Python.3.11 -e
winget install --id GitHub.cli -e
# optional, for PDF reports: winget install --id MiKTeX.MiKTeX -e   (or: quarto install tinytex)
```
macOS: `xcode-select --install; brew install python@3.11 gh` (optional TeX: `brew install --cask mactex-no-gui`
or `brew install tectonic`). Linux: `sudo apt install git python3.11 python3.11-venv gh`
(optional: `sudo apt install latexmk texlive-latex-extra`).

Verify: `git --version`, `python3.11 --version` (Windows: `py -3.11 --version`), `gh --version`.
Windows: PERSON runs `git config --global core.longpaths true` once (deep result paths).

## A2. GitHub access and clone

- PERSON: accept the invitation to `github.com/elle13eth/grnbench`.
- PERSON: `gh auth login` (GitHub.com, HTTPS, browser login).
- You: `gh repo clone elle13eth/grnbench ~/grnbench && cd ~/grnbench`.
  SSH alternative if the person already has a GitHub SSH key: `git clone git@github.com:elle13eth/grnbench.git ~/grnbench`.
- Local disk only (`~/grnbench`), never a Google Drive or other synced folder.
- `.gitattributes` keeps `*.sh`, `*.sbatch`, `*.py` as LF. Do not convert them to CRLF.

## A3. Python environment (pinned)

```bash
cd ~/grnbench
python3.11 -m venv .venv            # Windows: py -3.11 -m venv .venv
. .venv/bin/activate                # Windows Git Bash: . .venv/Scripts/activate
python -m pip install -r requirements/lock-py311.txt
python -m pip install --no-deps -e .
python -m pip check
```
Optional RegDiffusion model: the torch install line is at the top of `requirements/lock-py311.txt`.
Not installed editable? Prefix every Python command with `PYTHONPATH=src`.

## A4. Sanity tests (expected result)

```bash
python -m pytest -q -p no:cacheprovider tests/test_engine.py tests/test_context.py tests/test_simulate.py tests/test_ui.py
python -m pytest -q -p no:cacheprovider
```
Expected: exactly the floor written in `AGENTS.md` section 4 (failing and skipped test IDs, not
only counts). Report every failing test ID; never change code just to make a test pass.

## A5. Read before any work

Read `AGENTS.md`, run `python -m grnbench.board`, then read `.projectos/MAP.md`, `.projectos/STATE.md`,
`.projectos/PURPOSE.md`, the ACTIVE entries of `.projectos/DECISIONS.md` (D-012, D-013, D-015, D-017:
a product, not a metascience paper; honesty is a feature, not the pitch; headline claims only from
pre-registered, reproduced experiments), `experiments/README.md`, `paper/README.md`, `docs/HPC.md`.
Summarise each in two lines to the person.

## A6. HPC account and SSH (the person's own account)

1. PERSON: request / activate the Bocconi HPC account as Bocconi IT indicates (student ID), and do
   the first login the way they describe (it may require a password once, typed by the person).
2. PERSON: create a key with a passphrase: `ssh-keygen -t ed25519 -f ~/.ssh/id_bocconi`, and
   install the public key on the cluster (`ssh-copy-id -i ~/.ssh/id_bocconi.pub <user>@<login host>`
   where allowed, or as IT indicates).
3. PERSON: load the passphrase into ssh-agent.
   - Windows: once, PowerShell as admin, `Set-Service ssh-agent -StartupType Automatic; Start-Service ssh-agent`;
     then once, normal PowerShell, `ssh-add $HOME\.ssh\id_bocconi`. The service starts at every boot and the
     key stays registered: no per-reboot step.
   - macOS: `ssh-add --apple-use-keychain ~/.ssh/id_bocconi` (Keychain remembers it).
   - Linux: `eval "$(ssh-agent)"; ssh-add ~/.ssh/id_bocconi`.
4. PERSON: add to their own `~/.ssh/config` (never to the repo):
   ```
   Host bocconi-hpc
       HostName slogin.hpc.unibocconi.it
       User <your-student-id>
       IdentityFile ~/.ssh/id_bocconi
       IdentitiesOnly yes
       ServerAliveInterval 60
   ```
4b. VPN, so you can reconnect it later without the person (`docs/HPC.md` section 1b):
   - Windows: you run `powershell -NoProfile -ExecutionPolicy Bypass -File hpc/vpn.ps1 setup` in the
     person's normal shell. It prints one command with their account SID. PERSON runs that
     command once in an administrator PowerShell. Then you check with `hpc/vpn.ps1 status`.
   - Linux: you check `globalprotect show --status`. SAML sign-in from the CLI needs
     GlobalProtect 6.2.1+ and Firefox; otherwise PERSON connects with the GUI app.
   - macOS: no GlobalProtect CLI. PERSON connects in the app when you ask.
5. You: `cp .env.hpc.example .env.hpc`, set `GRNBENCH_HPC_HOST=bocconi-hpc`, check
   `git check-ignore .env.hpc` prints the file name.
6. You: `bash hpc/hpc.sh check`.
   - `Could not resolve hostname` / timeout: the VPN is down (the login node resolves to a private IP).
     You reconnect it first, with the command `hpc.sh` prints on its first error line (Windows:
     `powershell -NoProfile -ExecutionPolicy Bypass -File hpc/vpn.ps1 reconnect`; Linux:
     `globalprotect connect --portal vpnstudents.unibocconi.it`), then retry once. PERSON only if
     it exits 1 (sign in once in the GlobalProtect window), exits 2 (PERSON does the one-time
     grant: `hpc/vpn.ps1 setup` prints the admin command, `docs/HPC.md` section 1b), or on macOS
     (no CLI: PERSON connects in the GlobalProtect app).
   - `Permission denied (publickey)`: PERSON re-runs `ssh-add`; try once more only after that.
   - `Not allowed at this time`: lockout after failed logins. STOP, wait 10 minutes, fix the cause, try once.
     Never loop on authentication.
7. You: `bash hpc/hpc.sh sync`, then `bash hpc/hpc.sh run 'bash hpc/setup_env.sh'` (ends with `SETUP_DONE`).
8. First tiny job end to end:
   `bash hpc/hpc.sh submit "--partition=stud --qos=stud --time=2 --cpus-per-task=1 --mem=1G -o logs/hello_%j.out --wrap='hostname; nproc'"`,
   then wait for it as in B3 and read it with `bash hpc/hpc.sh log <jobid>`.
   Limits: 1 running + 1 queued job and 8 CPUs per user, 24 h wall time.

## A7. Data

- Approved sources only: `hpc/fetch_engine_datasets.sh` (references, BEELINE, PBMC3k, paul15,
  CELLxGENE) and `benchmarks/fetch_cellxgene.py`. New sources need an issue first.
- On the cluster (login node, downloads only): `docs/HPC.md` section 3.
- Locally, only if needed: `bash hpc/fetch_engine_datasets.sh ~/grnbench_data 30` (outside the repo,
  never on Google Drive: Drive disconnects killed long runs before). Motif models also need the
  cisTarget files and `GRNBENCH_CISTARGET_DIR`.
- Provenance: `python benchmarks/data_manifest.py ~/grnbench_data` snapshots the current setup
  cache. That cache is intentionally broad and partly dynamic (CELLxGENE selection and live
  reference networks), so do not overwrite `experiments/datasets.tsv` or expect the two files to
  be identical. Registered experiments use a separate data root populated from their committed
  selection manifest and verify its files with `--check`. Use `--require-all` only when the
  supplied manifest is intended to describe that root exactly.

## A8. First local run, then first registered experiment

1. Tiny local run on synthetic data (no downloads):
   ```bash
   grnbench demo --out ~/gb_scratch/demo
   grnbench run --data ~/gb_scratch/demo/counts.h5ad --tf-list ~/gb_scratch/demo/tfs.txt \
       --gold ~/gb_scratch/demo/gold.tsv --methods pearson_signed,spearman_signed --out ~/gb_scratch/demo_run
   ```
   (paths outside the repo; about 15 s; without `--methods` it also runs GENIE3/GRNBoost2 and takes
   more than 10 minutes). Open `~/gb_scratch/demo_run/eval/eval_report.html`.
2. First experiment record: pick a small question with the person, open an issue, create a branch
   `<initials>/<topic>`, copy `experiments/TEMPLATE.md` to `experiments/EXP-YYYYMMDD-<initials>-1.md`,
   fill sections 1-2, commit and push it **before** running, then follow part B.

# Part B. How to work (agents and humans)

## B1. Daily loop

1. `git pull`, `python -m grnbench.board`: what runs, what is open, where things left off.
2. Pick or open an issue, assign yourself, branch `<initials>/<topic>` from `main`.
3. Experiment? Claim it first: its issue labelled `experiment`, assigned to you (the board shows
   it as claimed). Then create its record (A8.2) and push. Confirmatory = hypothesis, endpoint,
   datasets, exclusions and comparison family (Holm via `grnbench compare`) fixed in the record.
4. Launch, then add the label `running` to the same issue (AGENTS.md section 6).
5. Wait correctly (B3), process the results (B4), update the record.
6. PR to `main` with the template filled, a CHANGELOG line under `## Unreleased`, test output,
   STATE/MAP updated; one review by another person; squash-merge. No force-push.
   The other person is away (exam, leave)? The owner may merge after an independent agent review
   posted on the PR with no blocking finding left open, labelled as such (AGENTS.md section 5);
   they can reopen or challenge it later.

## B2. Not colliding with the others

- Never relaunch a run with an open `running` issue; never write into another person's
  `benchmarks/results/<EXP-ID>/`; never reuse an ID (your initials keep your namespace separate).
- Split work by issue and by file; if two PRs touch the same file, the second one merges `main`
  into its branch (no rebase, no force-push on shared branches).
- Avoid stacked PRs. If one was stacked and its base was squash-merged, follow the recipe in
  AGENTS.md section 5: merge `origin/main` and keep your side for the files the base PR touched.

## B3. Waiting for HPC outputs

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
- `PENDING (QOSMaxJobsPerUserLimit / QOSMaxSubmitJobPerUserLimit)`: you already have a job; stop and
  report. The 1 running + 1 queued limit means batching all tasks of an experiment into one job.
- `PENDING (Resources / Priority)`: normal, can take hours. Comment the ETA (`squeue --start`) on the
  board issue and end the session; do not keep a session alive just to wait.
- Leaving `squeue` is not success: `sacct` State and ExitCode decide, then the job's own output.

## B4. Processing a finished run

1. `real_bench.py report --out <run> [--tasks-for-viz <tasks.json>]` on the cluster or locally: it
   writes `REPORT_REAL.md`, then `viz/` and `latex/` (a failure there is logged, not fatal).
2. Fetch (`bash hpc/hpc.sh fetch '~/gb_runs/<EXP-ID>' <local folder outside the repo>`), verify
   `sha256sum -c SHA256SUMS` and that the number of `rows/*.json` equals the number of tasks.
3. Copy the small files (reports, rows, CSV, `viz/`, `latex/`) to `benchmarks/results/<EXP-ID>/`
   (keep it under about 5 MB), fill section 3 of the record from the generated report: copy the
   numbers, never retype or round them by hand. Regenerate `python paper/make_figures.py` if a
   paper table uses it.
4. Close the board issue (`done` / `failed`) with a one-line result.

## B5. Failures and resubmits

- Runners are resumable: resubmit the same command, finished datasets (`rows/<task>.json`) are
  skipped. Delete a row file only to recompute that dataset on purpose, and say so in the record.
- A dataset that errors every time is a finding: keep it in the record, never drop it silently.
- TIMEOUT: resubmit; OOM: more `--mem`, not more CPUs (8 is the cap).

## B6. Long local runs

- Heavy work does not belong on a laptop. Experiments and model runs go to the HPC; full suites
  and long benches go to a team member's workstation if there is one. The laptop coordinates:
  it dispatches, polls, and reviews. The rules below are for when a laptop has to run the job anyway.
- Run from a frozen checkout (`git worktree add ../grnbench-run-<EXP> <commit>`, installed editable
  in the run's venv), not the tree you are editing, and from a local disk, never from Drive. While
  it runs, do not `pip install -e` another tree into that venv: `--jobs` workers re-import from it.
  The runners hold a keep-awake request on Windows; keep the laptop plugged in.
- Analysis code fixed after the run: fix + test, regenerate the report from the unchanged rows
  (`real_bench.py report --out <run>`), and record it as a dated deviation in the record.
- Windows quirks: Git Bash heredocs can eat backslashes (write paths with forward slashes, or write
  files with your editor tool); PowerShell and bash differ (`$env:X` vs `$X`); the HPC helper uses
  the Windows ssh client because Git Bash's ssh cannot see the Windows agent.

## B7. Handoff

The PR description is the handoff: done / running where (board issue, job ids) / next / do not.
Also update `.projectos/STATE.md` (and `MAP.md` if things moved). Personal agent memory stays
personal; team facts live in the repo and on GitHub.

# Part C. Look inside, from day one (visualization, exploration, explainability)

Scores alone keep the models a black box; the figures turn it grey. Rule: **every experiment record
links its exploration and explainability figures, not only scores.**

- Run real_bench with `--keep-artifacts` (edges, processed data, decisions per dataset), then
  `grnbench viz <run>` (or `real_bench.py report`, which calls it). The minimum per run:
  QC per cell before/after cleaning (counts, genes, % mito); PCA/UMAP coloured by cell type and by
  batch, before/after; the pipeline decisions the engine took and why; per-model score distributions
  and paired lifts over datasets; PR curves vs random; top-k predicted network vs reference; per-TF
  recovered vs missed edges and the top false positives (error analysis); for the motif prior, the
  TFs that drive the lift; calibration of edge scores; stability of top edges across contexts/seeds;
  for learned models, attributions (`grnbench.viz.gradient_x_input`, `integrated_gradients`).
- Every figure is titled with the question it answers and shows n; one style file
  (`src/grnbench/viz.mplstyle`), colour-blind-safe palette.
- `viz/index.md` lists what was made and what was **not** produced and why; copy that list into
  the record's deviations if it matters.
- `grnbench report <run|EXP-ID>` turns the numbers and the figures into a LaTeX report (PDF when a TeX
  engine is installed, else `.tex` + `report_fallback.md`). Worked example:
  `benchmarks/results/engine-2026-09-25-models/viz/index.md` and `.../latex/report.tex`.

# Part D. What never goes into git

Secrets, passphrases, `.env.hpc`, HPC usernames or student IDs, emails, absolute personal paths
(`C:\Users\<name>`, `/home/<id>`), host names, data files. Found one? Tell Luca.

# Part E. Report back to your person

Reply with this checklist, each line DONE / BLOCKED (why) / SKIPPED (why), with the evidence:

- [ ] tools: git, Python 3.11, gh versions (TeX optional)
- [ ] repo cloned at <path> (local disk, not Drive), `git status` clean, `main`
- [ ] venv from the lock, `pip check` clean
- [ ] quick tests + full suite summary lines; failing IDs vs the floor in AGENTS.md
- [ ] board output read; AGENTS, MAP, STATE, PURPOSE, DECISIONS, experiments/README, paper/README, docs/HPC.md summarised
- [ ] `~/.ssh/config` alias (person), key in ssh-agent (person), `.env.hpc` ignored by git
- [ ] VPN: Windows `hpc/vpn.ps1 status` says up and `reconnect` works without admin (grant done by the person); Linux `globalprotect show --status`; macOS: person connects in the app
- [ ] `hpc/hpc.sh check` output (user and host may be redacted in your reply)
- [ ] `hpc/hpc.sh sync` + `setup_env.sh` ended with SETUP_DONE
- [ ] tiny job: id, sacct State and ExitCode
- [ ] data location (cluster and/or local path outside the repo), manifest compared
- [ ] demo run done; first experiment record committed on a branch (EXP id)
- [ ] person's initials for IDs and branches: __
- [ ] anything the person must still do
