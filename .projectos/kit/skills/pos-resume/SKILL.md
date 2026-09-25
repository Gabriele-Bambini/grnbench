---
name: pos-resume
description: Use when resuming a project after a break or after another agent (Claude, Codex, Antigravity) worked on it, when the user asks "where were we / what's the status / what's next", after context compaction, or when handing work over to another agent or machine.
---

# pos-resume - pick up cold, hand over clean

## Resume (budget ~3k tokens before touching code)
1. Brief: injected by hook, else `pos brief`. Note: last checkpoint, auto-checkpoints after it, RUNNING attempts, other active agents, map freshness.
2. If auto-checkpoints exist after the last manual one: open the newest; its "Files edited" list shows work STATE.md may not mention. Verify those files quickly (diff/git status) before trusting STATE.md.
3. RUNNING attempts: check their evidence path (logs/metrics) and close them with `pos attempt update`.
4. MAP STALE for groups you will touch: dispatch incremental remap (skill `pos-map`) in background.
5. Tell the user in <=5 lines: purpose, where things stand, the next step you propose. Then proceed.

## Handover (to another agent/machine/session)
1. Skill `pos-checkpoint`, normal mode, `--reason handoff`; "Notes for the next agent" must contain commands, paths, gotchas.
2. `pos release` your claim.
3. Target agent can read the project folder (same Drive/git): nothing else needed; its SessionStart hook or `pos brief` picks it up.
   Target cannot read files (web chat, other account): `pos export --out HANDOFF.md` and give that file.

## Cross-project view
`pos list` shows every ProjectOS project on both machines (via graphify hub) with last checkpoint and map status.
