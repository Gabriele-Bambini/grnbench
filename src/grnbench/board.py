"""Team board: what is running, what is open, where things left off. Run it first in every session.

    python -m grnbench.board          # or: grnbench board

Live state lives in GitHub (issues labelled ``running``, claimed ``experiment`` issues, open PRs); durable state in
``.projectos/STATE.md``. Nothing here writes anything. Needs ``gh`` logged in; without it the
STATE.md part still prints.
"""
from __future__ import annotations

import json
import os
import shlex
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _gh(gh: list[str], args: list[str]) -> list[dict] | str:
    try:
        res = subprocess.run(gh + args, capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return f"gh not usable ({type(exc).__name__}: {exc})"
    if res.returncode != 0:
        return f"gh failed: {res.stderr.strip()[:200]}"
    return json.loads(res.stdout or "[]")


def _section(text: str, title: str, max_lines: int = 12) -> list[str]:
    out, on = [], False
    for line in text.splitlines():
        if line.startswith("## "):
            on = line[3:].strip().lower().startswith(title.lower())
            continue
        if on and line.strip():
            out.append(line)
    return out[:max_lines]


def _who(issue: dict) -> str:
    return ",".join(a.get("login", "?") for a in issue.get("assignees", [])) or "unassigned"


def board(gh: list[str] | None = None, root: Path | None = None) -> str:
    gh = gh or shlex.split(os.environ.get("GRNBENCH_GH", "gh"))
    root = root or (Path.cwd() if (Path.cwd() / ".projectos").is_dir() else ROOT)
    lines = ["# grnbench board", "", "## In flight (issues labelled `running`)"]
    running = _gh(gh, ["issue", "list", "--label", "running", "--state", "open",
                       "--json", "number,title,assignees,updatedAt,url"])
    if isinstance(running, str):
        lines.append(f"- ? {running}")
    elif not running:
        lines.append("- nothing running")
    for i in running if isinstance(running, list) else []:
        lines.append(f"- #{i['number']} {i['title']} ({_who(i)}, updated {i.get('updatedAt', '?')[:16]})")
    lines += ["", "## Claimed experiments (issues labelled `experiment`, not yet `running`)"]
    claimed = _gh(gh, ["issue", "list", "--label", "experiment", "--state", "open",
                       "--json", "number,title,assignees,labels"])
    if isinstance(claimed, str):
        lines.append(f"- ? {claimed}")
    else:
        claimed = [i for i in claimed if "running" not in {lb.get("name") for lb in i.get("labels", [])}]
        lines += [f"- #{i['number']} {i['title']} ({_who(i)})" for i in claimed] or ["- none"]
    lines += ["", "## Open pull requests"]
    prs = _gh(gh, ["pr", "list", "--state", "open", "--json", "number,title,author,headRefName,isDraft"])
    if isinstance(prs, str):
        lines.append(f"- ? {prs}")
    elif not prs:
        lines.append("- none")
    for p in prs if isinstance(prs, list) else []:
        draft = " [draft]" if p.get("isDraft") else ""
        lines.append(f"- #{p['number']} {p['title']} ({p.get('author', {}).get('login', '?')}, "
                     f"{p.get('headRefName', '?')}){draft}")
    state = root / ".projectos" / "STATE.md"
    text = state.read_text(encoding="utf-8") if state.exists() else ""
    for title in ("In progress", "Next steps", "Blockers"):
        sec = _section(text, title)
        if sec:
            lines += ["", f"## STATE.md: {title}"] + sec
    lines += ["", "Rule: never relaunch a run that has an open `running` issue; comment on it instead."]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    print(board())
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
