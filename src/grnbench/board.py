"""Team board: what is running, what is open, where things left off. Run it first in every session.

    python -m grnbench.board          # or: grnbench board

Live state lives in GitHub (issues labelled ``running``, claimed ``experiment`` issues, open PRs); durable state in
``.projectos/STATE.md``. Nothing here writes anything. Needs ``gh`` logged in; without it the
STATE.md part still prints.

    python -m grnbench.board --digest [--since YYYY-MM-DD] [--post ISSUE]   # weekly digest (default 7 days)
    python -m grnbench.board --claims                                        # claims status table

The digest reads git history and files from ``origin/main`` when it exists, else ``main``, else ``HEAD``, and
says which in its first line. ``--post ISSUE`` writes (a comment via ``gh issue comment``); without it, print only.
"""
from __future__ import annotations

import argparse
import datetime
import json
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path

from grnbench.project import find_project_root


def _split_command(value: str, *, posix: bool | None = None) -> list[str]:
    """Split an executable override without corrupting Windows paths.

    ``shlex`` uses POSIX escaping by default, so a value such as
    ``C:\\Users\\me\\bin\\gh.cmd`` loses every backslash.  Non-POSIX mode
    preserves the path but retains wrapping quotes; remove only matching
    outer quotes because ``subprocess`` receives an argument list.
    """
    use_posix = os.name != "nt" if posix is None else posix
    parts = shlex.split(value, posix=use_posix)
    if not use_posix:
        parts = [
            part[1:-1]
            if len(part) >= 2 and part[0] == part[-1] and part[0] in {"'", '"'}
            else part
            for part in parts
        ]
    return parts


def _gh(gh: list[str], args: list[str], cwd: Path | None = None) -> list[dict] | str:
    try:
        res = subprocess.run(
            gh + args,
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(cwd) if cwd is not None else None,
        )
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
    gh = gh or _split_command(os.environ.get("GRNBENCH_GH", "gh"))
    resolved_root = Path(root).resolve() if root is not None else find_project_root()
    if resolved_root is None:
        raise FileNotFoundError(
            "repository root not found (must be run within the grnbench repository)"
        )
    lines = ["# grnbench board", "", "## In flight (issues labelled `running`)"]
    running = _gh(
        gh,
        ["issue", "list", "--label", "running", "--state", "open", "--json", "number,title,assignees,updatedAt,url"],
        cwd=resolved_root,
    )
    if isinstance(running, str):
        lines.append(f"- ? {running}")
    elif not running:
        lines.append("- nothing running")
    for i in running if isinstance(running, list) else []:
        lines.append(f"- #{i['number']} {i['title']} ({_who(i)}, updated {i.get('updatedAt', '?')[:16]})")
    lines += ["", "## Claimed experiments (issues labelled `experiment`, not yet `running`)"]
    claimed = _gh(
        gh,
        ["issue", "list", "--label", "experiment", "--state", "open", "--json", "number,title,assignees,labels"],
        cwd=resolved_root,
    )
    if isinstance(claimed, str):
        lines.append(f"- ? {claimed}")
    else:
        claimed = [i for i in claimed if "running" not in {lb.get("name") for lb in i.get("labels", [])}]
        lines += [f"- #{i['number']} {i['title']} ({_who(i)})" for i in claimed] or ["- none"]
    lines += ["", "## Open pull requests"]
    prs = _gh(
        gh,
        ["pr", "list", "--state", "open", "--json", "number,title,author,headRefName,isDraft"],
        cwd=resolved_root,
    )
    if isinstance(prs, str):
        lines.append(f"- ? {prs}")
    elif not prs:
        lines.append("- none")
    for p in prs if isinstance(prs, list) else []:
        draft = " [draft]" if p.get("isDraft") else ""
        lines.append(f"- #{p['number']} {p['title']} ({p.get('author', {}).get('login', '?')}, "
                     f"{p.get('headRefName', '?')}){draft}")
    state = resolved_root / ".projectos" / "STATE.md"
    text = state.read_text(encoding="utf-8") if state.exists() else ""
    for title in ("In progress", "Next steps", "Blockers"):
        sec = _section(text, title)
        if sec:
            lines += ["", f"## STATE.md: {title}"] + sec
    lines += ["", "Rule: never relaunch a run that has an open `running` issue; comment on it instead."]
    return "\n".join(lines)


def _git(root: Path, args: list[str]) -> str:
    res = subprocess.run(["git", "-C", str(root)] + args, capture_output=True, text=True, encoding="utf-8")
    return res.stdout if res.returncode == 0 else ""


def _base(root: Path) -> str:
    """The ref the digest reads: origin/main when it exists, else main, else HEAD."""
    return next((r for r in ("origin/main", "main") if _git(root, ["rev-parse", "--verify", "-q", r]).strip()), "HEAD")


def _claims(text: str) -> dict[str, list[str]]:
    """Claim id -> [id, experiments, status, reproduced by, whole row]."""
    rows = {}
    for line in text.splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if re.fullmatch(r"C\d+", cells[0]):
            rows[cells[0]] = [cells[0], cells[2], cells[-2], cells[-1], line]
    return rows


def _decisions(text: str) -> dict[str, str]:
    return {m[1]: m[2] for m in re.finditer(r"^### (D-\d+) \| [^|]*\|[^|]*\|[^|]*\| (.*)$", text, re.M)}


def _lines(res: list[dict] | str, fmt) -> list[str]:
    return [f"- ? {res}"] if isinstance(res, str) else [fmt(r) for r in res] or ["- none"]


def digest(since: str | None = None, gh: list[str] | None = None, root: Path | None = None) -> str:
    gh = gh or _split_command(os.environ.get("GRNBENCH_GH", "gh"))
    root = Path(root).resolve() if root is not None else find_project_root()
    if root is None:
        raise FileNotFoundError(
            "repository root not found (must be run within the grnbench repository)"
        )
    since = since or (datetime.date.today() - datetime.timedelta(days=7)).isoformat()
    login = lambda p: p.get("author", {}).get("login", "?")
    base = _base(root)
    out = [f"# grnbench digest since {since} (git history and files from {base})", "", "## PRs merged"]
    out += _lines(_gh(gh, ["pr", "list", "--state", "merged", "--search", f"merged:>={since}",
                           "--json", "number,title,author"]),
                  lambda p: f"- #{p['number']} {p['title']} ({login(p)})")
    out += ["", "## PRs open"]
    out += _lines(_gh(gh, ["pr", "list", "--state", "open", "--json", "number,title,author,isDraft,reviewDecision"]),
                  lambda p: f"- #{p['number']} {p['title']} ({login(p)}, "
                            f"{'draft' if p.get('isDraft') else (p.get('reviewDecision') or 'no review yet').lower()})")
    out += ["", "## Experiment records added or changed"]
    files = sorted({f for f in _git(root, ["log", base, f"--since={since}", "--name-only", "--format=",
                                           "--", "experiments/EXP-*.md"]).split() if f.endswith(".md")})
    recs = []
    for f in files:
        txt = _git(root, ["show", f"{base}:{f}"])
        typ = re.search(r"^- type: \**(\w+)", txt, re.M)
        st = re.search(r"^- status: (.*)", txt, re.M)
        recs.append(f"- {Path(f).stem} ({typ[1] if typ else '?'})" + (f": status {st[1].strip()[:80]}" if st else ""))
    out += recs or ["- none"]
    rev = _git(root, ["rev-list", "-1", f"--before={since}", base]).strip()
    old = lambda path: _git(root, ["show", f"{rev}:{path}"]) if rev else ""
    now = lambda path: _git(root, ["show", f"{base}:{path}"])
    before, after = _claims(old("paper/CLAIMS.md")), _claims(now("paper/CLAIMS.md"))
    out += ["", "## Claims whose status changed"]
    out += [f"- {c}: {before[c][2] if c in before else '(new)'} -> {r[2]}"
            for c, r in after.items() if c not in before or before[c][2] != r[2]] or ["- none"]
    d_old, d_new = _decisions(old(".projectos/DECISIONS.md")), _decisions(now(".projectos/DECISIONS.md"))
    out += ["", "## Decisions added"]
    out += [f"- {d}: {t}" for d, t in d_new.items() if d not in d_old] or ["- none"]
    out += ["", "## Open issues labelled running/experiment"]
    issues, err = {}, None
    for label in ("running", "experiment"):
        res = _gh(gh, ["issue", "list", "--label", label, "--state", "open", "--json", "number,title,assignees"])
        if isinstance(res, str):
            err = res
        else:
            issues.update({i["number"]: i for i in res})
    out += _lines(err or list(issues.values()), lambda i: f"- #{i['number']} {i['title']} ({_who(i)})")
    return "\n".join(out)


def claims_table(root: Path | None = None) -> str:
    root = Path(root).resolve() if root is not None else find_project_root()
    if root is None:
        raise FileNotFoundError(
            "repository root not found (must be run within the grnbench repository)"
        )
    outputs = [p for d in ("tables", "figures") for p in sorted((root / "paper" / d).glob("*")) if p.is_file()]
    texts = {p: p.name + (p.read_text(encoding="utf-8", errors="ignore") if p.suffix in (".md", ".tex", ".csv") else "")
             for p in outputs}
    out = ["| claim | status | reproduced by | experiments | issues/PRs | generated files |", "|---|---|---|---|---|---|"]
    for cid, (_, exps, status, by, row) in _claims((root / "paper" / "CLAIMS.md").read_text(encoding="utf-8")).items():
        ids = sorted(set(re.findall(r"EXP-\d{8}-[A-Z]+-\d+", exps)))
        refs = sorted(set(re.findall(r"#\d+", row)), key=lambda r: int(r[1:]))
        gen = [p.relative_to(root).as_posix() for p, t in texts.items() if any(i in t for i in ids)]
        out.append(f"| {cid} | {status} | {by} | {', '.join(ids) or '-'} | {', '.join(refs) or '-'} | "
                   f"{', '.join(gen) or '-'} |")
    return "\n".join(out)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="grnbench board")
    ap.add_argument("--digest", action="store_true", help="weekly digest as Markdown")
    ap.add_argument("--since", help="digest start date YYYY-MM-DD (default: 7 days ago)")
    ap.add_argument("--claims", action="store_true", help="claims status table")
    ap.add_argument("--post", metavar="ISSUE", help="with --digest: also post it as a comment on this issue")
    args = ap.parse_args(argv)
    try:
        if args.claims:
            print(claims_table())
        elif args.digest or args.post:
            text = digest(args.since)
            print(text)
            if args.post:
                gh = _split_command(os.environ.get("GRNBENCH_GH", "gh"))
                res = subprocess.run(
                    gh + ["issue", "comment", args.post, "--body-file", "-"],
                    input=text,
                    capture_output=True,
                    text=True,
                    encoding="utf-8",
                )
                print(res.stdout.strip() or res.stderr.strip(), file=sys.stderr)
                return res.returncode
        else:
            print(board())
    except (OSError, ValueError) as exc:
        print(f"grnbench board: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
