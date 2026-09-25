# grnbench - rules for every agent and every person

This file is the single source of rules. Claude Code, Codex and Antigravity all read it
(`CLAUDE.md` only imports it). If something here conflicts with another file, this file wins;
fix the other file in the same PR.

## 1. What we are building (do not re-open)

- grnbench is a **product to publish**: a context-aware data-cleaning engine for gene regulatory
  network (GRN) inference, with an interface, GRN models mounted on top, and side-by-side model
  comparison, including the user's own model (D-012, D-013 in `.projectos/DECISIONS.md`).
- It is **not a metascience paper** (D-011 was rejected). Statistical honesty is a feature of the
  product, never its pitch.
- Headline claims come from **real data**; synthetic gains are upper bounds of the simulator (D-015).
- The target is a paper at the level of a Nature-family journal: every claim traceable to a
  registered experiment, reproduced by a second person (see `paper/README.md`).

## 2. People

| person | GitHub | initials | reviews |
|---|---|---|---|
| Luca | @elle13eth | LM | engine, statistics, paper/ |
| Gabriele | @Gabriele-Bambini | GB | engine, HPC, paper/ |
| (third colleague, TBD) | - | - | - |

Initials are used in experiment IDs and branch names.

## 3. Start of every session (in this order)

1. This file.
2. `.projectos/STATE.md` (current focus) and `.projectos/PURPOSE.md`.
3. `.projectos/DECISIONS.md`: every ACTIVE decision. Do not act against one; propose a new D-xxx.
4. `experiments/README.md`: the registry index. Check whether your idea already ran or failed.
5. The GitHub Issue you are working on (`gh issue view N`). No issue: create one first (section 6).

## 4. Environment and tests (exact commands)

Python 3.11. From the repo root:

```bash
python -m venv .venv                              # Windows: py -3.11 -m venv .venv
. .venv/bin/activate                              # Windows Git Bash: . .venv/Scripts/activate
python -m pip install -r requirements/lock-py311.txt
python -m pip install --no-deps -e .
# quick sanity (about 1 min)
python -m pytest -q -p no:cacheprovider tests/test_engine.py tests/test_context.py tests/test_simulate.py tests/test_ui.py
# full suite (about 3 min)
python -m pytest -q -p no:cacheprovider
```

Known floor (2026-09-25, Windows, Python 3.11.9): `357 passed, 1 failed`. The failure is
`tests/test_prep_gaps.py::test_drop_doublets_uses_scrublet_and_logs_what_it_removed`: scanpy's
scrublet needs scikit-image, which the lock lacks, so doublet calling falls back to `heuristic`
(issue #2; CI deselects this one test until the lock is fixed). Compare failing test IDs before and
after your change, not counts. Not installed editable? Prefix commands with `PYTHONPATH=src`.
Optional models: RegDiffusion needs torch (install line at the top of the lock file); the motif
models need `GRNBENCH_CISTARGET_DIR` pointing at the cisTarget files (`experiments/datasets.tsv`).

## 5. Branches, PRs, merges

- Never commit to `main`. Branch per issue: `<initials>/<topic>`, e.g. `gb/knee-calling-gui`.
- Small commits. **No force-push**, ever, on any shared branch.
- Every PR: fill `.github/pull_request_template.md` (what, why, experiment IDs, tests with real
  output, risks, handoff). Add a line under `## Unreleased` in `CHANGELOG.md`.
- Review: one approval from a person other than the author (a team rule: the free plan has no
  branch protection for private repos, so GitHub does not enforce it; only squash-merge is enabled). Changes to `paper/`, to
  `src/grnbench/eval/` or to statistical code need a reviewer from the "reviews" column above.
- Merge: **squash-merge** (one commit per PR on `main`; the PR text is the record). Delete the branch.
- A change of direction, scope, metric, data split or library = a decision `D-xxx` in
  `.projectos/DECISIONS.md` in the same PR (next free number; renumber if `main` took it first).

## 6. Tracking the work

- **Issues are the task list.** One issue per piece of work; assignee = whoever does it. Labels:
  `experiment`, `engine`, `hpc`, `ui`, `bug`, `decision-needed`.
- **Weekly status**: every Monday each person comments 3 lines on the pinned issue
  "Weekly status" (done / next / blocked). Nothing else to maintain.
- `.projectos/STATE.md` is the agents' view of the current focus: update it in the PR that
  changes the focus or the next steps.

## 7. Experiments (full rules: `experiments/README.md`)

- **No run without a record.** Create `experiments/EXP-YYYYMMDD-<initials>-<n>.md` from
  `experiments/TEMPLATE.md` and commit it before the run starts (local or HPC).
- Mark it `confirmatory` (pre-registered: hypothesis, primary metric, unit = dataset, dataset list,
  exclusions, declared comparison family, Holm via `grnbench compare` / `compare_units`) or
  `exploratory`. Exploratory results never feed a headline claim without a confirmatory rerun.
- Results go to `benchmarks/results/<EXP-ID>/` (small files only) or, if large, outside git with
  the path in the record. **Never write into another person's EXP folder**; never reuse an ID.
- Numbers come from scripts only. Never hand-edit a result or a table.
- Negative and null results are recorded with the same care as positive ones.
- Past attempts that are not runs (debugging hypotheses, libraries tried, reverted refactors) go
  to `.projectos/ATTEMPTS.md`. Runs never go there.

## 8. Paper (full rules: `paper/README.md`)

- Tables and figures are generated by `python paper/make_figures.py` from committed result files.
- `paper/CLAIMS.md` maps every sentence-level claim to experiment IDs and result files.
- No result enters `paper/` until a second person reproduced it from its record.
- Log contributions in `paper/CONTRIBUTIONS.md` (CRediT roles) in the PR that does the work.

## 9. Data

- Data never goes into git and never onto Google Drive (Drive disconnects killed long runs).
  Keep it in a local folder outside the repo (e.g. `~/grnbench_data`) or on the cluster.
- Only the datasets fetched by `hpc/fetch_engine_datasets.sh` are approved. Anything new: open an
  issue with source, licence and size first.
- Provenance lives in `experiments/datasets.tsv` (`python benchmarks/data_manifest.py DATA_ROOT`).
- Run code from a local clone, never from a Drive folder.

## 10. HPC (full guide: `docs/HPC.md`)

- Everyone uses their **own** Bocconi account through their own `~/.ssh/config` alias and an
  untracked `.env.hpc`. Scripts take the alias from `.env.hpc`; no username is ever committed.
- The SSH key passphrase is loaded **by the person** into ssh-agent. Agents never type, store or
  ask for passphrases or passwords; `hpc/hpc.sh` runs ssh with `BatchMode=yes`.
- `Not allowed at this time` = lockout after failed logins: **stop, wait 10 minutes**, fix the key.
  Never retry auth in a loop.
- QOS `stud` allows 1 running + 1 queued job per user and 8 CPUs: batch many tasks into one job.
- No compute on the login node (downloads and `sbatch` only).

## 11. Secrets and personal information

Never commit: keys, passphrases, tokens, passwords, `.env.hpc`, HPC usernames or student IDs,
emails, absolute personal paths (`C:\Users\<name>\...`, `/home/<id>/...`), machine host names,
data files. Write paths relative to the repo or as `~/...`. Found one already committed? Tell
Luca; do not just delete it in a new commit (it stays in history).

## 12. Handoff

The PR description is the handoff: what was done, what is running where (job IDs, EXP IDs),
what is next, what the next person must not do. Also update `.projectos/STATE.md`.
The old Drive `HANDOFF.md` habit is retired (`docs/history/` keeps the last one).

## 13. Known traps

- decoupler 2.x removed `dc.get_dorothea/get_collectri`: use `grnbench.tfs.decoupler_network`.
- BEELINE inputs are log-normalized, not counts: the engine refuses them by design.
- Doublet calling runs on all cells before any context subsetting (A-014).
- High-mitochondrial tissues (liver, heart, kidney, muscle): the planner relaxes the 25% gate from
  the tissue context (A-020); elsewhere QC gates fail loudly (D-004).
- Cell-type labels contain `:` or `|`: context folders go through `autoprep.safe_folder_name`.
- Do not edit engine code while a local run imports it: run long benches from a frozen copy of
  `src/` + `benchmarks/`. Runners are resumable (finished datasets are skipped).
- Long local runs on Windows: the runners hold a keep-awake request; keep the machine plugged in.
- Shell scripts must stay LF (`.gitattributes` enforces it); CRLF breaks bash on the cluster.
- ProjectOS: `pos` is optional (`python .projectos/kit/install.py` changes your global agent
  hooks, so only install it on your own machine if you want it). Without it, edit the markdown
  files directly in the formats they already use.
