---
name: pos-log
description: Use when an approach, experiment, training run, hypothesis test, debugging theory, refactor or tool evaluation starts or concludes (success OR failure), when a design/scope decision is made, or before trying something to check whether it was already tried. Maintains the append-only ATTEMPTS.md and DECISIONS.md so agents reason on old errors instead of repeating them.
---

# pos-log - attempts & decisions

## Before trying something
`pos attempt list --grep "kw1 kw2"` (add `--full` for details, `--status failed` for failures only).
If a match is FAILED/ABANDONED/INCONCLUSIVE: write one line on what is different this time, or choose another path.

## Log an attempt
```
pos attempt add -t "<short title>" -s <running|success|partial|failed|abandoned|inconclusive> \
  --goal "<target, ideally a metric>" --approach "<what exactly: params, config, commit>" \
  --result "<outcome + numbers>" --evidence "<log/metrics file, commit, notebook cell>" \
  --lesson "<what to do/avoid next time>" --tags "<kw1,kw2>"
```
- Long jobs: log `-s running` at launch; later `pos attempt update A-007 -s success --result "..." --lesson "..."`.
- FAILED/ABANDONED without `--lesson` is useless: always say why and what next.
- One attempt = one approach. Tuning sweeps: one entry with the grid + best/worst, not one per run.
- Numbers over adjectives ("AUPRC 0.41 vs 0.47 baseline", not "worse").

## Record a decision
```
pos decide -t "<title>" --decision "<what>" --why "<evidence/reasons>" --alternatives "<rejected + why>" [--supersedes D-002] [--proposed]
```
Decisions: architecture/library choice, data split, metric, scope cut, PURPOSE change (always `--proposed` until the user confirms). `pos decisions` lists those in force.

## No CLI
Append by hand with the exact heading format (parsed by tools):
`### A-012 | 2026-09-22 | codex | FAILED | title` followed by `- goal: / - approach: / - result: / - evidence: / - lesson: / - tags:` lines.
