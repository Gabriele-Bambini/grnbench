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

## 3. Session routine (every agent, every machine)

Start:
1. `git pull` on `main`, then read this file.
2. `python -m grnbench.board` (or `grnbench board`): runs in flight (issues labelled `running`),
   open PRs, and the in-progress / next-steps / blockers of `.projectos/STATE.md`.
   **Never relaunch a run that has an open `running` issue**: comment on the issue instead.
3. `.projectos/MAP.md` (where things are), `.projectos/PURPOSE.md`, the ACTIVE entries of
   `.projectos/DECISIONS.md` (do not act against one; propose a new D-xxx).
4. `experiments/README.md` and the records relevant to your task: did it already run or fail?
5. The GitHub issue you work on (`gh issue view N`). No issue: create one first (section 6).

End:
1. Update `.projectos/STATE.md` (and `MAP.md` if you added, moved or removed things) in your PR.
2. Update or close your board issues (section 6).
3. Write the handoff in the PR description (section 12).

Team memory lives only in the repo and on GitHub (MAP, STATE, DECISIONS, ATTEMPTS, `experiments/`,
issues, PRs). Your own tool's memory (Claude memory, Codex/Antigravity notes) stays personal:
never rely on it for team facts, and never copy personal paths or accounts from it into the repo.

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

Known floor (2026-09-25, fresh venv from the lock, Windows, Python 3.11.9): `368 passed, 1 skipped`
(the skip is the RegDiffusion test when torch/regdiffusion are not installed). Compare failing test
IDs before and after your change, not counts. Not installed editable? Prefix commands with `PYTHONPATH=src`.
CI (`.github/workflows/tests.yml`) is disabled: GitHub refuses to start Actions jobs on this
account for billing reasons, and we do not pay. The test output pasted in each PR is the check;
re-enable with `gh workflow enable tests` if that changes.
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
- `.projectos/STATE.md` = where things left off; `.projectos/MAP.md` = where things are. Both stay
  short and are updated in every PR that changes them (PR template checkbox).
- **Live board (what runs right now).** Never a shared file (concurrent edits conflict). Each launched
  run (local or HPC) gets one issue labelled `running`:
  - title `<EXP-ID> <short name>`; assignee = the person whose agent launched it;
  - body: who, machine or HPC, job ids, started at, expected finish, where outputs land, how to
    check (`bash hpc/hpc.sh status <id>`);
  - the agent comments when it checks the job (state, ETA), and when results are fetched and
    recorded it closes the issue with label `done` or `failed` and a one-line result.
  ```bash
  gh issue create --label running --assignee @me --title "EXP-... short name" --body "..."
  gh issue comment N --body "checked: RUNNING 3h10m, 9/20 rows"
  gh issue edit N --remove-label running --add-label done && gh issue close N --comment "result: ..."
  ```

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
- **Look inside, not only at scores** (from day one): every record links its exploration and
  explainability figures. Run real_bench with `--keep-artifacts`; `real_bench.py report` then
  writes `viz/` (`grnbench viz`) and `latex/` (`grnbench report`). Minimum per run: QC before/after,
  embeddings by batch and cell type before/after, the pipeline decisions and why, per-model score
  distributions and paired lifts, PR curves vs random, top-k network vs reference, TF-level
  recovered/missed edges (error analysis), which TFs drive a motif lift, calibration of edge scores,
  stability of top edges across contexts/seeds; for learned models, attributions
  (`viz.gradient_x_input`, `viz.integrated_gradients`). What could not be produced is listed in
  `viz/index.md` "Not produced": say why in the record.
- Reports are generated, never typed: `grnbench report <run|EXP>` writes LaTeX with every number
  from result files (a test forbids digits in templates); no TeX engine -> `.tex` + `report_fallback.md`.
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

## 10. HPC (full guide: `docs/HPC.md`; agent procedure: `.claude/skills/grnbench-hpc/SKILL.md`)

Codex and Antigravity: read that skill file before any cluster action (Antigravity also finds it
through `.agents/skills/grnbench-hpc/`). Claude Code loads it as the `grnbench-hpc` skill.

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

The PR description is the handoff: what was done, what is running where (board issue, job IDs,
EXP IDs), what is next, what the next person must not do. Also update `.projectos/STATE.md`.
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
