## What changed

## Why (issue #, decision D-xxx if direction changed)

## Experiments
<!-- EXP IDs registered, run or reproduced in this PR; "none" if none -->

## Tests (command + real output)
```
python -m pytest -q -p no:cacheprovider
```
<!-- paste the summary line and any failing test IDs; compare IDs with the known floor in AGENTS.md -->

## Risks / what a reviewer should check

## Handoff
<!-- done / still running where (job IDs) / next step / what the next person must not do -->

- [ ] CHANGELOG.md `## Unreleased` updated
- [ ] no secrets, HPC usernames, emails, personal paths or data in the diff
- [ ] `.projectos/STATE.md` (where things left off) and `.projectos/MAP.md` (where things are) updated if changed
- [ ] board issues updated or closed (label `done` / `failed`) for runs this PR records
- [ ] experiment records link their `viz/` figures and `latex/` report
- [ ] `paper/` touched? claims reproduced by a second person, `python paper/make_figures.py` leaves no diff
