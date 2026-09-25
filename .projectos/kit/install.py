#!/usr/bin/env python3
"""
ProjectOS installer - wires the `pos` engine, skills, subagents, hooks and global instruction
blocks into Claude Code, Codex and Antigravity/Gemini on THIS machine. Idempotent; re-run after
updating the package. Every touched config is backed up to ~/.projectos/backups/<timestamp>/.

  python install.py              install / update
  python install.py --dry-run    show what would change
  python install.py --uninstall  remove hooks, blocks, skills, shims (keeps project data)
  flags: --no-claude --no-codex --no-gemini --no-statusline
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import shutil
import sys
from pathlib import Path

PKG = Path(__file__).resolve().parent
HOME = Path.home()
POS_HOME = HOME / ".projectos"
ENGINE = POS_HOME / "engine"
BIN = HOME / ".local" / "bin"
SKILLS = ["projectos", "pos-map", "pos-checkpoint", "pos-log", "pos-resume"]
AGENTS = ["pos-mapper", "pos-architect"]
MARK = "projectos/engine/pos.py"
BEGIN, END = "<!-- projectos:begin", "<!-- projectos:end -->"
PY = "python" if os.name == "nt" else "python3"
STAMP = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
DRY = False
LOG = []


def say(msg):
    LOG.append(msg)
    print(("[dry] " if DRY else "") + msg)


def read(p):
    try:
        return Path(p).read_text(encoding="utf-8-sig")
    except Exception:
        return ""


def write(p, text):
    if DRY:
        return
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)


def backup(p):
    p = Path(p)
    if DRY or not p.exists():
        return
    dst = POS_HOME / "backups" / STAMP / (p.parent.name + "__" + p.name)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(p, dst)


def copytree(src, dst):
    if DRY:
        return
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)


def cmd(event, agent):
    engine = (ENGINE / "pos.py").as_posix()
    return f'{PY} "{engine}" hook {event} --agent {agent}'


# ------------------------------------------------------------------ instruction blocks

COMMON = """- Every project folder carries `.projectos/` = single source of truth: PURPOSE (why/done-when), MAP (architecture), STATE (live), ATTEMPTS (what was tried + lessons), DECISIONS, checkpoints. Same protocol for Claude, Codex, Antigravity: `.projectos/PROTOCOL.md` (master copy `~/.projectos/engine/PROTOCOL.md`). Skills: `projectos`, `pos-map`, `pos-checkpoint`, `pos-log`, `pos-resume`.
- Code questions: `.projectos/MAP.md` first, raw files last. Before a new approach: `pos attempt list --grep <kw>`; never repeat a FAILED attempt without stating what differs.
- Log outcomes (`pos attempt add`, `pos decide`). After each finished unit of work: update STATE.md + `pos checkpoint -m "..."`. Usage warning from ProjectOS = checkpoint immediately (skill pos-checkpoint, emergency mode).
- Another agent active on the project (brief shows it) = stay off its files, `pos claim "<task>"`, coordinate through STATE.md.
- Model economy: mechanical work = `pos` script (free); mapping/summaries = cheapest tier; top tier only for novel reasoning. Never map a project with the top model or by reading the whole repo.
- New projects are initialized automatically (`pos init`, which also runs the graphify hub bootstrap) - no manual bootstrap needed.
- `.projectos/` wins over other memory tools (claude-mem, cortex, openviking) when they disagree."""

START = {
    "claude": "- Session start: a SessionStart hook injects the ProjectOS brief (auto-initializes unmapped projects). If the brief says MAP MISSING/STALE, dispatch the haiku subagent `pos-mapper` in background as instructed; never map with Opus.",
    "codex": "- Session start: a SessionStart hook injects the ProjectOS brief (needs one-time trust via `/hooks`). If absent, run `pos brief` first. MAP MISSING/STALE: map with the cheapest model / a mini sub-agent (skill pos-map). `pos usage` shows context + rate limits.",
    "gemini": "- Session start (no hooks here): FIRST run `pos brief` in the project root; if it says no project, run `pos init`, then map via skill pos-map using the fastest/cheapest model. Checkpoint after EVERY unit of work (you cannot see usage limits).",
}


def block(agent):
    return f"""<!-- projectos:begin v1 -->
# ProjectOS - project memory, attempts, checkpoints (all agents, all projects)
{START[agent]}
{COMMON}
<!-- projectos:end -->"""


def upsert(path, text):
    path = Path(path)
    cur = read(path)
    if BEGIN in cur and END in cur:
        new = re.sub(re.escape(BEGIN) + r".*?" + re.escape(END), lambda _: text, cur, flags=re.S)
    else:
        new = (cur.rstrip() + "\n\n" + text + "\n") if cur.strip() else text + "\n"
    if new != cur:
        backup(path)
        write(path, new)
        say(f"block -> {path}")


def remove_block(path):
    path = Path(path)
    cur = read(path)
    if BEGIN in cur:
        backup(path)
        new = re.sub(r"\n*" + re.escape(BEGIN) + r".*?" + re.escape(END) + r"\n?", "\n", cur, flags=re.S)
        write(path, new.rstrip() + "\n")
        say(f"block removed <- {path}")


# ------------------------------------------------------------------ hooks

def strip_ours(hooks):
    out = {}
    for ev, groups in (hooks or {}).items():
        kept = []
        for g in groups:
            hs = [h for h in g.get("hooks", []) if MARK not in str(h.get("command", "")).replace("\\", "/")]
            if hs:
                kept.append({**g, "hooks": hs})
        if kept:
            out[ev] = kept
    return out


def claude_hooks():
    return {
        "SessionStart": [{"hooks": [{"type": "command", "command": cmd("session-start", "claude"), "timeout": 60,
                                     "statusMessage": "ProjectOS brief..."}]}],
        "UserPromptSubmit": [{"hooks": [{"type": "command", "command": cmd("prompt", "claude"), "timeout": 10}]}],
        "Stop": [{"hooks": [{"type": "command", "command": cmd("stop", "claude"), "timeout": 15}]}],
        "PreCompact": [{"hooks": [{"type": "command", "command": cmd("precompact", "claude"), "timeout": 30}]}],
        "SessionEnd": [{"hooks": [{"type": "command", "command": cmd("session-end", "claude"), "timeout": 5}]}],
    }


def codex_hooks():
    return {
        "SessionStart": [{"matcher": "startup|resume|clear|compact",
                          "hooks": [{"type": "command", "command": cmd("session-start", "codex"), "timeout": 60,
                                     "additionalContextLimit": 8000, "statusMessage": "ProjectOS brief"}]}],
        "UserPromptSubmit": [{"hooks": [{"type": "command", "command": cmd("prompt", "codex"), "timeout": 10}]}],
        "Stop": [{"hooks": [{"type": "command", "command": cmd("stop", "codex"), "timeout": 15}]}],
        "PreCompact": [{"hooks": [{"type": "command", "command": cmd("precompact", "codex"), "timeout": 30}]}],
        "SessionEnd": [{"hooks": [{"type": "command", "command": cmd("session-end", "codex"), "timeout": 3}]}],
    }


def merge_hooks(existing, ours):
    res = strip_ours(existing)
    for ev, groups in ours.items():
        res.setdefault(ev, []).extend(groups)
    return res


def load_json(p):
    t = read(p)
    if not t.strip():
        return {}
    return json.loads(t)


# ------------------------------------------------------------------ steps

def install_engine():
    for f in ("pos.py", "PROTOCOL.md"):
        src = PKG / "engine" / f
        if not DRY:
            ENGINE.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, ENGINE / f)
    write(ENGINE / "VERSION", read(PKG / "VERSION") or "1.0.0\n")
    write(ENGINE / "SOURCE", str(PKG) + "\n")
    say(f"engine -> {ENGINE}")
    if not (POS_HOME / "config.json").exists():
        write(POS_HOME / "config.json", json.dumps({
            "warn_pct": 75, "force_pct": 88, "periodic_checkpoint_min": 60, "auto_init": True,
            "graphify_bootstrap": True, "ignore_paths": []}, indent=2) + "\n")
        say(f"config -> {POS_HOME / 'config.json'} (edit thresholds / ignore_paths here)")
    engine = ENGINE / "pos.py"
    write(BIN / "pos.cmd", f'@echo off\r\n{PY} "{engine}" %*\r\n')
    write(BIN / "pos", f'#!/usr/bin/env sh\nexec {PY} "{engine.as_posix()}" "$@"\n')
    if not DRY and os.name != "nt":
        os.chmod(BIN / "pos", 0o755)
    say(f"CLI shims -> {BIN} (pos, pos.cmd)")
    if str(BIN).lower() not in os.environ.get("PATH", "").lower():
        say(f"WARNING: {BIN} is not on PATH; add it or call `{PY} {engine}`")


def install_skills(dest_root, label):
    dest_root = Path(dest_root)
    for s in SKILLS:
        copytree(PKG / "skills" / s, dest_root / s)
    say(f"skills ({label}) -> {dest_root}: {', '.join(SKILLS)}")


def install_claude(statusline=True):
    c = HOME / ".claude"
    install_skills(c / "skills", "claude")
    for a in AGENTS:
        if not DRY:
            (c / "agents").mkdir(parents=True, exist_ok=True)
            shutil.copy2(PKG / "agents" / "claude" / f"{a}.md", c / "agents" / f"{a}.md")
    say(f"subagents -> {c / 'agents'}: pos-mapper (haiku), pos-architect (sonnet)")
    sp = c / "settings.json"
    s = load_json(sp)
    s["hooks"] = merge_hooks(s.get("hooks"), claude_hooks())
    if statusline:
        cur = s.get("statusLine")
        if not cur or MARK in str(cur.get("command", "")).replace("\\", "/"):
            s["statusLine"] = {"type": "command", "command": f'{PY} "{(ENGINE / "pos.py").as_posix()}" statusline'}
            say("statusLine -> pos statusline (captures 5h/7d rate limits + context %; shows [CAVEMAN] badge)")
        else:
            say("statusLine: custom one present, left untouched (rate-limit capture for Claude CLI disabled)")
    backup(sp)
    write(sp, json.dumps(s, indent=2, ensure_ascii=False) + "\n")
    say(f"hooks -> {sp}: SessionStart, UserPromptSubmit, Stop, PreCompact, SessionEnd")
    upsert(c / "CLAUDE.md", block("claude"))


def install_codex():
    x = HOME / ".codex"
    if not x.exists():
        say("codex: ~/.codex not found, skipped")
        return
    install_skills(x / "skills", "codex")
    hp = x / "hooks.json"
    h = load_json(hp)
    h["hooks"] = merge_hooks(h.get("hooks"), codex_hooks())
    backup(hp)
    write(hp, json.dumps(h, indent=2) + "\n")
    say(f"hooks -> {hp} (open Codex and approve them once via /hooks)")
    upsert(x / "AGENTS.md", block("codex"))


def install_gemini():
    g = HOME / ".gemini"
    if not g.exists():
        say("gemini/antigravity: ~/.gemini not found, skipped")
        return
    install_skills(g / "skills", "antigravity/gemini")
    upsert(g / "GEMINI.md", block("gemini"))


def ledger():
    try:
        sys.path.insert(0, str(ENGINE))
        import pos  # noqa
        h = pos.hub_dir()
        if h and not DRY:
            anchor = pos.jload(h / "shared" / "anchor.json", {})
            entry = {"ts_utc": dt.datetime.now(dt.timezone.utc).isoformat(), "anchor_id": anchor.get("anchor_id", ""),
                     "project_id": "global-config", "project_path": str(PKG), "action": "projectos-install",
                     "details": f"ProjectOS {read(PKG / 'VERSION').strip()} installed on {pos.MACHINE}",
                     "user": os.environ.get("USERNAME", ""), "host": pos.MACHINE}
            pos.append(h / "shared" / "acts" / "ledger.ndjson", json.dumps(entry) + "\n")
    except Exception as e:
        say(f"(ledger skipped: {e})")


def uninstall(a):
    if not a.no_claude:
        c = HOME / ".claude"
        sp = c / "settings.json"
        s = load_json(sp)
        s["hooks"] = strip_ours(s.get("hooks"))
        if MARK in str((s.get("statusLine") or {}).get("command", "")).replace("\\", "/"):
            s.pop("statusLine")
        backup(sp)
        write(sp, json.dumps(s, indent=2, ensure_ascii=False) + "\n")
        say(f"hooks removed <- {sp}")
        for sk in SKILLS:
            if (c / "skills" / sk).exists() and not DRY:
                shutil.rmtree(c / "skills" / sk)
        for ag in AGENTS:
            if (c / "agents" / f"{ag}.md").exists() and not DRY:
                (c / "agents" / f"{ag}.md").unlink()
        remove_block(c / "CLAUDE.md")
    if not a.no_codex and (HOME / ".codex").exists():
        x = HOME / ".codex"
        hp = x / "hooks.json"
        if hp.exists():
            h = load_json(hp)
            h["hooks"] = strip_ours(h.get("hooks"))
            backup(hp)
            write(hp, json.dumps(h, indent=2) + "\n")
            say(f"hooks removed <- {hp}")
        for sk in SKILLS:
            if (x / "skills" / sk).exists() and not DRY:
                shutil.rmtree(x / "skills" / sk)
        remove_block(x / "AGENTS.md")
    if not a.no_gemini and (HOME / ".gemini").exists():
        g = HOME / ".gemini"
        for sk in SKILLS:
            if (g / "skills" / sk).exists() and not DRY:
                shutil.rmtree(g / "skills" / sk)
        remove_block(g / "GEMINI.md")
    for f in ("pos", "pos.cmd"):
        if (BIN / f).exists() and not DRY:
            (BIN / f).unlink()
    say("uninstalled. Kept: ~/.projectos (config, runtime, backups) and every project's .projectos/ folder.")


def main():
    global DRY
    ap = argparse.ArgumentParser()
    ap.add_argument("--uninstall", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--no-claude", action="store_true")
    ap.add_argument("--no-codex", action="store_true")
    ap.add_argument("--no-gemini", action="store_true")
    ap.add_argument("--no-statusline", action="store_true")
    a = ap.parse_args()
    DRY = a.dry_run
    if sys.version_info < (3, 9):
        sys.exit("ProjectOS needs Python >= 3.9")
    if a.uninstall:
        uninstall(a)
        return
    install_engine()
    if not a.no_claude:
        install_claude(statusline=not a.no_statusline)
    if not a.no_codex:
        install_codex()
    if not a.no_gemini:
        install_gemini()
    ledger()
    print("\nDone. Verify with: pos doctor")
    print("Codex: open it once and approve the new hooks with /hooks.")
    print("Existing projects are initialized automatically the next time an agent opens them;")
    print("to do all at once (tier-0, no LLM): pos discover \"<folder>\" --init")


if __name__ == "__main__":
    main()
