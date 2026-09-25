---
name: pos-checkpoint
description: Use after finishing any unit of work (feature, fix, experiment, analysis), before risky operations, at the end of a session, and IMMEDIATELY when a ProjectOS message reports high usage (context or rate limit) or asks for a periodic checkpoint. Writes STATE.md + an immutable checkpoint so another agent or a later session can continue without loss.
---

# pos-checkpoint

## Normal mode (~1 min)
1. Concluded an approach/experiment? `pos attempt add ...` / `pos attempt update A-xxx ...` (skill `pos-log`). Made a design choice? `pos decide ...`.
2. Edit `.projectos/STATE.md` (re-read it first if a tool reports it changed):
   - `Updated:` line -> today + your agent name
   - **Focus**: 1-3 lines, current objective
   - **Done recently**: what this unit achieved (+ evidence path/commit)
   - **In progress**: anything half-done, with exact file/function and what remains
   - **Next steps**: numbered, concrete, first step executable without re-reading code
   - **Blockers / open questions**
   - **Notes for the next agent**: gotchas, commands, where results are
   Trim stale bullets; STATE is live state, not history (history = checkpoints + ATTEMPTS).
3. `pos checkpoint -m "<one-line outcome>" [--advances "<PURPOSE criterion>"]`

## Emergency mode (usage forced; must finish in <2 min)
Skip step 1 unless an attempt is RUNNING without status. In STATE.md update only Focus / In progress / Next steps / Notes (terse bullets). Then:
`pos checkpoint -m "<summary>" --reason usage`
Then tell the user in one line that a checkpoint was saved and where to resume.

## End of session
Normal mode with `--reason end`. If another agent will continue: `pos export --out HANDOFF.md` only when the next agent cannot read `.projectos/`.

## What NOT to do
- No prose essays in STATE.md; no copying code; no duplicating MAP.md content.
- Do not edit or delete files in `checkpoints/` (immutable). Auto checkpoints (precompact/session-end) are written by hooks; don't replicate them.
- No `pos`? Write `checkpoints/YYYY-MM-DD_HHMMSS_<agent>_<slug>.md` by hand (PROTOCOL section 9).
