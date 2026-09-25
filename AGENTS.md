<!-- projectos:begin v1 -->
## ProjectOS (read first)
This project keeps shared agent memory in `.projectos/` (purpose, map, state, attempts, decisions, checkpoints) for every agent: Claude, Codex, Antigravity, humans.
1. Start: read `.projectos/STATE.md` and `PURPOSE.md` (or `pos brief`). Hooks may have injected this already.
2. Code questions: `.projectos/MAP.md` first, raw files last.
3. Before a new approach: `pos attempt list --grep <kw>`. Never repeat a FAILED attempt without stating what differs.
4. Log outcomes (`pos attempt add`, `pos decide`). After each finished sub-task: update STATE.md, then `pos checkpoint -m "..."`.
Rules: `.projectos/PROTOCOL.md`. No `pos` on this machine? `python .projectos/kit/install.py` if `.projectos/kit/` exists (full setup), else edit the markdown files directly, same formats.
<!-- projectos:end -->
