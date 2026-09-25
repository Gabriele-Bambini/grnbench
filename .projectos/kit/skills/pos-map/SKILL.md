---
name: pos-map
description: Use when a project has no .projectos/MAP.md, when the ProjectOS brief says "MAP MISSING" or "MAP STALE", when the user asks to map/remap a project, or to map ALL projects under a folder. Produces the architecture map cheaply (tier-0 script scan + cheapest model), never with the top-tier model.
---

# pos-map - cheap, incremental project mapping

Principle: the script does the reading (free), the cheapest model does the writing, the main agent only dispatches.

## 1. Tier-0 (always first, zero tokens)
`pos init` (if no `.projectos/`) or `pos scan` -> `.projectos/scan.md`: layout, deps, entry points, tests, git, symbol index, **size class and mapping plan**.

## 2. Dispatch by size class (plan line in scan.md)
| class | plan |
|---|---|
| S | 1 mapper on whole project |
| M | 2-3 mappers in parallel (one per group listed in plan) + 1 mapper merge |
| L | up to 5 mappers in parallel + 1 `pos-architect` (sonnet tier) synthesis from parts only |

**Claude Code**: `Agent(subagent_type="pos-mapper", run_in_background=true, prompt=...)`. pos-mapper is pinned to haiku.
- S: prompt `Map project at <root>. Mode: full.`
- M/L: one agent per group: `Map project at <root>. Mode: part. Groups: <g1,g2>. Write .projectos/.runtime/map-parts/<first-group>.md`.
  When all parts are done: M -> one more pos-mapper `Mode: merge`; L -> `Agent(subagent_type="pos-architect", prompt="Synthesize MAP.md for <root> from .projectos/.runtime/map-parts/")`.
- Keep serving the user's request while mappers run; use scan.md meanwhile.

**Codex**: if sub-agents are available, spawn one on the smallest model (e.g. a `*-mini` model, low reasoning effort) with the same prompts; otherwise do it yourself following section 3 at low reasoning effort.
**Antigravity/others**: do it yourself following section 3; prefer a fast/cheap model setting.

## 3. Mapper procedure (what pos-mapper does; follow it when mapping yourself)
Budget: scan.md + graphify summary + README (<=150 lines) + at most 8 files (S) / 6 files per group (M/L). Files over 300 lines: first 120 lines + the symbol index. Never open data/output dirs, lock files, binaries, notebook outputs.
MAP.md skeleton (keep under ~250 lines):
```
# Map - <name>
## One-liner
## How to run          (setup, main commands, tests; copy from scan entry points, verify quickly)
## Architecture        (data/control flow in 3-8 bullets; optional small mermaid)
## Modules
### <group> - <role>
- `path` - what it does; key symbols; depends on
## Data & artifacts    (where inputs/outputs live; sizes from scan; "do not read" list)
## Conventions & gotchas
## Open questions      (things not verified - mark them, don't guess)
```
Finish with `pos map-done --model <haiku|sonnet|model-id> --confidence <high|medium|low> [--groups g1,g2]`.

## 4. Incremental remap (MAP STALE)
Only the groups listed as changed: one pos-mapper `Mode: incremental. Groups: <g...>` -> it rewrites only those `### <group>` subsections (and How to run / Architecture if affected), then `pos map-done --groups <g...>`.

## 5. Map ALL projects (`mode: all`)
1. `pos discover <root> --init --json` (tier-0 for every unmapped project, graphify bootstrap included, no LLM).
2. For each project with status `initialized`: dispatch mapping per its size class, max 4 mappers in flight, S projects first.
3. Report a table: project, size class, mapped by, confidence. `pos list` shows the result.
