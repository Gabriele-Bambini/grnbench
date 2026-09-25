---
name: projectos
description: Use when starting, resuming or planning work in any project folder, before choosing an approach (to check what was already tried), when the user asks about project goals/history/progress, or when a "ProjectOS" hook message appears. Shared project memory (purpose, map, state, attempts, decisions, checkpoints) in `.projectos/`, identical for Claude Code, Codex and Antigravity, driven by the `pos` CLI.
---

# ProjectOS - operating loop

Full rules: `.projectos/PROTOCOL.md` (read once per project if unsure). This is the working summary.

## Loop
1. **Orient (<=2k tokens)**: use the injected brief; if none, run `pos brief`. Then PURPOSE -> STATE -> only the MAP.md section you need.
   No `.projectos/`? Run `pos init` (free, no LLM), then skill `pos-map`.
2. **Check history before acting**: `pos attempt list --grep "kw1 kw2"` and `pos decisions`.
   Similar attempt FAILED/ABANDONED -> say what differs, or pick another path.
3. **Purpose check**: name the PURPOSE success criterion the work advances. None -> tell the user (scope creep).
4. **Work** in small verifiable units. Other agent active (brief says so)? `pos claim "<task>"`, stay off its files.
5. **Close each unit**: log (`pos-log` skill) -> update STATE.md -> `pos checkpoint -m "<summary>"` (`pos-checkpoint` skill).
6. **Hook says usage high** -> `pos-checkpoint` emergency mode immediately, then continue.

## Command cheat-sheet
| need | command |
|---|---|
| context digest | `pos brief` |
| health (map freshness, last checkpoint, sessions, usage) | `pos status` |
| what was tried | `pos attempt list [--grep kw] [--status failed] [--full]` |
| log attempt | `pos attempt add -t "..." -s failed --approach "..." --result "..." --lesson "..."` |
| decision | `pos decide -t "..." --decision "..." --why "..." [--proposed]` |
| checkpoint | `pos checkpoint -m "..." [--advances "criterion"] [--reason usage]` |
| usage now | `pos usage` |
| claim / release task | `pos claim "task"` / `pos release` |
| handoff bundle | `pos export --out HANDOFF.md` |
| all projects | `pos list` |

`pos` missing on PATH -> `python ~/.projectos/engine/pos.py ...`. No Python -> edit the markdown files by hand (PROTOCOL section 9).

## Sibling skills
- mapping / remapping / map all projects -> `pos-map`
- writing a checkpoint, usage warnings -> `pos-checkpoint`
- logging attempts and decisions -> `pos-log`
- resume, "where were we", handoff to another agent -> `pos-resume`

## Model economy
Mechanical work (scan, digest, checkpoint) = `pos` script, zero tokens. Mapping/summaries = cheapest tier (Claude: haiku subagent `pos-mapper`; Codex: mini model). Sonnet-tier only for large-map synthesis (`pos-architect`) or real code work. Top tier only for novel reasoning. Never re-read a whole codebase: MAP.md exists for that.
