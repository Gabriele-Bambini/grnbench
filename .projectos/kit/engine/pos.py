#!/usr/bin/env python3
"""
ProjectOS (`pos`) - agent-agnostic project memory for Claude Code, Codex, Antigravity/Gemini
and humans. Stdlib only, single file.

Per-project data  : <project>/.projectos/   (portable markdown, travels with the project)
Machine runtime   : ~/.projectos/            (config, usage snapshots, per-session hook state)

Run `pos --help` for commands. Hooks call `pos hook <event> --agent <claude|codex>`.
"""
from __future__ import annotations

import argparse
import ast
import datetime as dt
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time
import uuid
from pathlib import Path

VERSION = "1.0.0"
PROTOCOL_VERSION = 1
DIRNAME = ".projectos"

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

HOME = Path.home()
POS_HOME = Path(os.environ.get("PROJECTOS_HOME") or (HOME / ".projectos"))
RUNTIME = POS_HOME / "runtime"
ENGINE_DIR = Path(__file__).resolve().parent
MACHINE = platform.node() or "unknown"

DEFAULT_CONFIG = {
    "warn_pct": 75,                 # usage % that triggers a soft reminder
    "force_pct": 88,                # usage % that forces a checkpoint (second force at +7)
    "periodic_checkpoint_min": 60,  # Stop hook nudges a checkpoint if edits and none for this long
    "auto_init": True,              # SessionStart initializes unmapped project folders
    "graphify_bootstrap": True,     # run the graphify hub bootstrap for projects not in project-map.json
    "context_windows": {},          # model-id prefix -> window, for transcript-only usage (e.g. {"claude-sonnet-4": 200000})
    "brief_max_chars": 7000,
    "keep_checkpoints": 40,         # older checkpoints are condensed into checkpoints/_archive.md
    "session_active_min": 45,       # a session is "active" if its heartbeat is newer than this
    "hub_candidates": [
        "~/graphify-context",
        "~/Desktop/graphify total context",
    ],
    "ignore_paths": [],
}

# --------------------------------------------------------------------------- basics

_CFG = None


def cfg():
    global _CFG
    if _CFG is None:
        c = dict(DEFAULT_CONFIG)
        user = jload(POS_HOME / "config.json", {})
        if isinstance(user, dict):
            c.update(user)
        _CFG = c
    return _CFG


def now():
    return dt.datetime.now().astimezone()


def iso(t=None):
    return (t or now()).isoformat(timespec="seconds")


def epoch():
    return time.time()


def parse_ts(s):
    if not s:
        return None
    try:
        return dt.datetime.fromisoformat(str(s).replace("Z", "+00:00")).timestamp()
    except Exception:
        return None


def ago(sec):
    if sec is None:
        return "never"
    sec = max(0, int(sec))
    if sec < 90:
        return f"{sec}s ago"
    if sec < 5400:
        return f"{sec // 60}m ago"
    if sec < 172800:
        return f"{sec // 3600}h ago"
    return f"{sec // 86400}d ago"


def read(p, default=""):
    try:
        return Path(p).read_text(encoding="utf-8-sig", errors="replace")
    except Exception:
        return default


def write(p, text):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_name(p.name + ".tmp")
    try:
        with open(tmp, "w", encoding="utf-8", newline="\n") as f:
            f.write(text)
        os.replace(tmp, p)
    except Exception:
        try:
            tmp.unlink()
        except Exception:
            pass
        with open(p, "w", encoding="utf-8", newline="\n") as f:
            f.write(text)


def append(p, text):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "a", encoding="utf-8", newline="\n") as f:
        f.write(text)


def jload(p, default):
    t = read(p, "")
    if not t.strip():
        return default
    try:
        return json.loads(t)
    except Exception:
        return default


def jsave(p, obj):
    write(p, json.dumps(obj, indent=2, ensure_ascii=False))


def slugify(s, n=48):
    s = re.sub(r"[^a-z0-9]+", "-", str(s).lower()).strip("-")
    return s[:n].strip("-") or "x"


def clip(s, n):
    s = (s or "").strip()
    return s if len(s) <= n else s[: max(0, n - 3)].rstrip() + "..."


def one_line(s, n=160):
    return clip(re.sub(r"\s+", " ", s or ""), n)


def log_err(where, e):
    try:
        append(RUNTIME / "errors.log", f"{iso()} {where}: {type(e).__name__}: {e}\n")
    except Exception:
        pass


def detect_agent(default="agent"):
    e = os.environ
    if e.get("CLAUDECODE") or e.get("CLAUDE_CODE_ENTRYPOINT"):
        return "claude"
    if any(k.startswith("CODEX_") and k != "CODEX_HOME" for k in e):
        return "codex"
    if any("ANTIGRAVITY" in k for k in e):
        return "antigravity"
    if e.get("GEMINI_CLI"):
        return "gemini"
    return default


def env_session_id():
    e = os.environ
    return e.get("CLAUDE_CODE_SESSION_ID") or e.get("CODEX_THREAD_ID") or e.get("CODEX_SESSION_ID")


# --------------------------------------------------------------------------- project resolution

MARKERS = {
    ".git", "pyproject.toml", "setup.py", "setup.cfg", "requirements.txt", "package.json", "Cargo.toml",
    "go.mod", "environment.yml", "environment.yaml", "DESCRIPTION", "Makefile", "Snakefile", "main.nf",
    "main.tex", "pom.xml", "build.gradle", "pubspec.yaml", "Gemfile", "renv.lock", "CMakeLists.txt",
    "Dockerfile", "Project.toml", DIRNAME,
}
MARKER_SUFFIXES = (".Rproj", ".ipynb", ".sln", ".csproj", ".qmd", ".Rmd")
IGNORE_NAMES = {
    "desktop", "downloads", "documents", "my drive", "other computers", "my pc", "onedrive",
    "shared drives", "computers", "pictures", "music", "videos", "users",
}


def find_root(start):
    try:
        p = Path(start).resolve()
    except Exception:
        return None
    for d in [p, *p.parents]:
        if (d / DIRNAME / "project.json").exists():
            return d
    return None


def git_root(start):
    p = Path(start).resolve()
    for d in [p, *p.parents]:
        if (d / ".git").exists():
            return d
    return None


def ignored_reason(p):
    try:
        p = Path(p).resolve()
    except Exception:
        return "unresolvable"
    if len(p.parts) <= 1:
        return "drive root"
    home = HOME.resolve()
    if p == home:
        return "home directory"
    nm = p.name.lower()
    if nm in IGNORE_NAMES or nm.startswith("onedrive - "):
        return "container folder"
    guards = [home / ".claude", home / ".codex", home / ".gemini", POS_HOME]
    if not os.environ.get("POS_ALLOW_TEMP"):  # test hook only
        guards += [home / "AppData", Path(tempfile.gettempdir())]
    for var in ("SystemRoot", "ProgramFiles", "ProgramFiles(x86)"):
        if os.environ.get(var):
            guards.append(Path(os.environ[var]))
    for extra in cfg().get("ignore_paths", []):
        guards.append(Path(os.path.expanduser(extra)))
    for g in guards:
        try:
            g = g.resolve()
            if p == g or g in p.parents:
                return f"ignored path ({g})"
        except Exception:
            pass
    if (p / ".projectos-ignore").exists():
        return ".projectos-ignore present"
    return None


def has_marker(d):
    try:
        names = os.listdir(d)
    except Exception:
        return False
    if MARKERS.intersection(names):
        return True
    return any(n.endswith(MARKER_SUFFIXES) for n in names)


def child_projects(d, limit=80):
    out = []
    try:
        entries = [e for e in os.scandir(d) if e.is_dir(follow_symlinks=False)]
    except Exception:
        return out
    for e in entries[:limit]:
        if e.name.startswith(".") or e.name in SKIP_DIRS:
            continue
        if (Path(e.path) / DIRNAME / "project.json").exists() or has_marker(e.path):
            out.append(Path(e.path))
    return out


def classify_dir(d):
    """'project' or 'workspace' (a folder that mainly contains other projects)."""
    if has_marker(d):
        return "project"
    return "workspace" if len(child_projects(d)) >= 2 else "project"


def need_root(path=None):
    root = find_root(Path(path) if path else Path.cwd())
    if not root:
        sys.exit("pos: no ProjectOS project here (run `pos init` in the project root).")
    return root


def pdir(root):
    return Path(root) / DIRNAME


def meta(root):
    return jload(pdir(root) / "project.json", {})


# --------------------------------------------------------------------------- hub (graphify shared context)

def hub_dir():
    for c in cfg().get("hub_candidates", []):
        try:
            p = Path(os.path.expanduser(c))
            if (p / "shared" / "anchor.json").exists():
                return p
        except Exception:
            pass
    return None


def ledger(root, action, details=""):
    h = hub_dir()
    if not h:
        return
    try:
        anchor = jload(h / "shared" / "anchor.json", {})
        m = meta(root)
        entry = {
            "ts_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "anchor_id": anchor.get("anchor_id", ""),
            "project_id": m.get("id", ""),
            "project_path": str(root),
            "action": action,
            "details": one_line(details, 300),
            "user": os.environ.get("USERNAME") or os.environ.get("USER") or "",
            "host": MACHINE,
        }
        append(h / "shared" / "acts" / "ledger.ndjson", json.dumps(entry, ensure_ascii=False) + "\n")
    except Exception as e:
        log_err("ledger", e)


def register(root, extra=None):
    """Record project in local registry and in the hub (one file per project -> no sync conflicts)."""
    m = meta(root)
    rec_extra = extra or {}
    try:
        reg = jload(POS_HOME / "registry.json", {})
        reg[str(root)] = {"id": m.get("id"), "name": m.get("name"), "updated": iso(), **rec_extra}
        jsave(POS_HOME / "registry.json", reg)
    except Exception as e:
        log_err("registry", e)
    h = hub_dir()
    if not h or not m.get("id"):
        return
    try:
        p = h / "projects" / m["id"] / "projectos.json"
        rec = jload(p, {})
        paths = rec.get("paths", {})
        paths[MACHINE] = str(root)
        rec.update({"id": m["id"], "uid": m.get("uid"), "name": m.get("name"), "paths": paths,
                    "updated": iso(), **rec_extra})
        jsave(p, rec)
    except Exception as e:
        log_err("hub_register", e)


def unique_id(root, uid):
    base = slugify(Path(root).name)
    h = hub_dir()
    if h:
        rec = jload(h / "projects" / base / "projectos.json", {})
        if rec.get("uid") and rec.get("uid") != uid:
            return f"{base}-{uid[:6]}"
    return base


def maybe_graphify(root, wait=False):
    if os.name != "nt" or not cfg().get("graphify_bootstrap"):
        return "off"
    h = hub_dir()
    script = h / "shared" / "scripts" / "bootstrap-project.ps1" if h else None
    if not script or not script.exists():
        return "no-hub"
    mp = jload(h / "shared" / "project-map.json", {})
    key = str(root).lower().rstrip("\\/")
    if any(str(k).lower().rstrip("\\/") == key for k in mp):
        return "already"
    cmd = ["powershell", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File", str(script), "-ProjectPath", str(root)]
    if wait:  # bulk init: run serially so project-map.json writes don't race
        try:
            subprocess.run(cmd, stdin=subprocess.DEVNULL, capture_output=True, timeout=90)
            return "done"
        except Exception as e:
            log_err("graphify", e)
            return "error"
    try:
        subprocess.Popen(
            cmd,
            stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            creationflags=0x00000008 | 0x00000200, close_fds=True,
        )
        return "started"
    except Exception as e:
        log_err("graphify", e)
        return "error"


# --------------------------------------------------------------------------- tier-0 scan (no LLM)

SKIP_DIRS = {
    ".git", ".hg", ".svn", "node_modules", "__pycache__", ".venv", "venv", "env", ".env", ".tox", ".nox",
    ".mypy_cache", ".pytest_cache", ".ruff_cache", ".ipynb_checkpoints", "site-packages", ".idea", ".vscode",
    "graphify-out", DIRNAME, ".next", ".nuxt", ".cache", ".gradle", ".dart_tool", "Pods", ".terraform",
    "target", "dist", "build", ".eggs", "htmlcov", "renv", "packrat", ".Rproj.user", "bower_components",
    "vendor", ".claude", ".codex", ".agent", ".gemini", ".github",
}
DATA_DIRS = {
    "data", "datasets", "dataset", "raw", "processed", "out", "outputs", "output", "results", "result",
    "runs", "checkpoints", "ckpt", "ckpts", "artifacts", "figures", "figs", "plots", "images", "img",
    "logs", "log", "cache", "tmp", "temp", "wandb", "mlruns", "lightning_logs", "saved_models", "weights",
    "media",
}
LANG = {
    ".py": "Python", ".ipynb": "Notebook", ".r": "R", ".rmd": "RMarkdown", ".qmd": "Quarto", ".jl": "Julia",
    ".js": "JavaScript", ".jsx": "JavaScript", ".mjs": "JavaScript", ".cjs": "JavaScript", ".ts": "TypeScript",
    ".tsx": "TypeScript", ".c": "C", ".h": "C", ".cpp": "C++", ".cc": "C++", ".hpp": "C++", ".cu": "CUDA",
    ".java": "Java", ".kt": "Kotlin", ".go": "Go", ".rs": "Rust", ".rb": "Ruby", ".php": "PHP",
    ".swift": "Swift", ".dart": "Dart", ".m": "MATLAB", ".scala": "Scala", ".sh": "Shell", ".ps1": "PowerShell",
    ".sql": "SQL", ".tex": "LaTeX", ".html": "HTML", ".css": "CSS", ".scss": "CSS", ".vue": "Vue",
    ".svelte": "Svelte", ".nf": "Nextflow", ".smk": "Snakemake", ".lua": "Lua", ".cs": "C#",
}
DOC_EXT = {".md", ".rst", ".txt", ".docx", ".pdf", ".pptx", ".odt"}
NON_CODE_LANGS = {"LaTeX", "HTML", "CSS"}
TODO_RE = re.compile(r"\b(TODO|FIXME|XXX|HACK)\b")


def _is_data_dir(path, name):
    if name.lower() not in DATA_DIRS:
        return False
    try:
        for e in os.scandir(path):
            if e.name == "__init__.py":
                return False
            if e.is_file() and Path(e.name).suffix.lower() in LANG and Path(e.name).suffix.lower() != ".ipynb":
                return False
    except Exception:
        pass
    return True


def walk(root, budget_s=8.0, max_files=25000):
    root = Path(root)
    t0 = time.monotonic()
    files, data_dirs, truncated = [], [], False
    stack = [root]
    while stack:
        d = stack.pop()
        try:
            entries = list(os.scandir(d))
        except Exception:
            continue
        for e in entries:
            nm = e.name
            try:
                if e.is_dir(follow_symlinks=False):
                    if nm.startswith(".") or nm in SKIP_DIRS or nm.lower().endswith(".egg-info"):
                        continue
                    if _is_data_dir(e.path, nm):
                        try:
                            n = sum(1 for _ in os.scandir(e.path))
                        except Exception:
                            n = -1
                        data_dirs.append((Path(e.path).relative_to(root).as_posix(), n))
                        continue
                    stack.append(Path(e.path))
                elif e.is_file(follow_symlinks=False):
                    if nm.endswith(".tmp") or nm == "desktop.ini":
                        continue
                    st = e.stat()
                    rel = Path(e.path).relative_to(root).as_posix()
                    files.append((rel, st.st_size, Path(nm).suffix.lower()))
            except Exception:
                continue
        if len(files) >= max_files or time.monotonic() - t0 > budget_s:
            truncated = bool(stack) or len(files) >= max_files
            break
    return files, data_dirs, truncated


def group_of(rel):
    return rel.split("/", 1)[0] if "/" in rel else "."


AGENT_FILES = {"AGENTS.md", "CLAUDE.md", "GEMINI.md", ".graphify-context.json"}  # edited by tooling, not architecture


def fingerprints(files):
    by = {}
    for rel, size, _ in files:
        if rel in AGENT_FILES:
            continue
        by.setdefault(group_of(rel), []).append(f"{rel}|{size}")
    return {g: hashlib.sha1("\n".join(sorted(v)).encode()).hexdigest()[:12] for g, v in by.items()}


def _py_symbols(text):
    try:
        tree = ast.parse(text)
    except Exception:
        return None
    doc = (ast.get_docstring(tree) or "").strip().splitlines()
    items, main = [], False
    for n in tree.body:
        if isinstance(n, ast.ClassDef):
            meths = [b.name for b in n.body if isinstance(b, (ast.FunctionDef, ast.AsyncFunctionDef))
                     and not b.name.startswith("_")]
            bases = ", ".join(ast.unparse(b) for b in n.bases[:2]) if n.bases else ""
            items.append(f"class {n.name}({bases})" + (f"[{len(meths)}m]" if meths else ""))
        elif isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)) and not n.name.startswith("_"):
            args = [a.arg for a in n.args.args if a.arg not in ("self", "cls")][:4]
            items.append(f"{n.name}({', '.join(args)})")
        elif isinstance(n, ast.If):
            try:
                if "__main__" in ast.unparse(n.test):
                    main = True
            except Exception:
                pass
    return {"doc": doc[0] if doc else "", "items": items, "main": main}


def _nb_info(text):
    try:
        nb = json.loads(text)
    except Exception:
        return None
    cells = nb.get("cells", [])
    heads, imports = [], set()
    ncode = nmd = 0
    for c in cells:
        src = c.get("source", "")
        src = "".join(src) if isinstance(src, list) else str(src)
        if c.get("cell_type") == "markdown":
            nmd += 1
            for ln in src.splitlines():
                if re.match(r"^#{1,3}\s", ln) and len(heads) < 8:
                    heads.append(ln.lstrip("#").strip()[:60])
        elif c.get("cell_type") == "code":
            ncode += 1
            for m in re.finditer(r"^\s*(?:from|import)\s+([A-Za-z_][\w]*)", src, re.M):
                imports.add(m.group(1))
    return {"heads": heads, "code": ncode, "md": nmd, "imports": sorted(imports)[:12]}


def _regex_symbols(text, ext):
    if ext in (".r", ".rmd"):
        pat = r"^\s*([A-Za-z.][\w.]*)\s*(?:<-|=)\s*function"
    elif ext in (".js", ".jsx", ".mjs", ".cjs", ".ts", ".tsx"):
        pat = r"^export\s+(?:default\s+)?(?:async\s+)?(?:function|class|const|let|interface|type)\s+(\w+)"
    else:
        return []
    return [m.group(1) for m in re.finditer(pat, text, re.M)][:12]


def _deps(root):
    out = {}
    root = Path(root)
    reqs = []
    for f in list(root.glob("requirements*.txt")) + list(root.glob("requirements/*.txt")):
        for ln in read(f).splitlines():
            ln = ln.strip()
            if not ln or ln.startswith(("#", "-")):
                continue
            m = re.match(r"([A-Za-z0-9_.\-\[\]]+)", ln)
            if m and m.group(1) not in reqs:
                reqs.append(m.group(1))
    if reqs:
        out["requirements"] = reqs[:40]
    pp = root / "pyproject.toml"
    if pp.exists():
        try:
            import tomllib
            d = tomllib.loads(read(pp))
            proj = d.get("project", {})
            deps = [re.split(r"[<>=!~; \[]", x, 1)[0] for x in proj.get("dependencies", [])]
            poetry = list((d.get("tool", {}).get("poetry", {}).get("dependencies", {}) or {}).keys())
            out["pyproject"] = {
                "name": proj.get("name") or d.get("tool", {}).get("poetry", {}).get("name"),
                "description": proj.get("description", ""),
                "deps": (deps or poetry)[:40],
                "scripts": proj.get("scripts", {}) or d.get("tool", {}).get("poetry", {}).get("scripts", {}),
                "optional": list((proj.get("optional-dependencies") or {}).keys())[:10],
            }
        except Exception as e:
            out["pyproject"] = {"error": str(e)[:80]}
    for envf in ("environment.yml", "environment.yaml"):
        if (root / envf).exists():
            names = re.findall(r"^\s*-\s*([A-Za-z0-9_.\-]+)", read(root / envf), re.M)
            out["conda"] = [n for n in names if n not in ("pip",)][:40]
    pj = root / "package.json"
    if pj.exists():
        d = jload(pj, {})
        out["npm"] = {"name": d.get("name"), "deps": list((d.get("dependencies") or {}).keys())[:30],
                      "dev": list((d.get("devDependencies") or {}).keys())[:15],
                      "scripts": d.get("scripts") or {}}
    desc = root / "DESCRIPTION"
    if desc.exists():
        m = re.search(r"^Imports:(.*?)(?:^\S|\Z)", read(desc), re.M | re.S)
        if m:
            out["r"] = [x.strip().split(" ")[0] for x in m.group(1).split(",") if x.strip()][:30]
    return out


def _readme(root):
    for nm in ("README.md", "README.rst", "README.txt", "readme.md", "Readme.md", "README"):
        p = Path(root) / nm
        if p.exists():
            t = read(p)
            heads = [ln.lstrip("#").strip() for ln in t.splitlines() if re.match(r"^#{1,3}\s", ln)][:15]
            para = ""
            for block in re.split(r"\n\s*\n", t):
                b = "\n".join(ln for ln in block.splitlines() if not ln.lstrip().startswith("#")).strip()
                if b and not b.startswith(("#", "!", "[!", "<", "```", "|", "=")) and len(b) > 40:
                    para = one_line(b, 500)
                    break
            return {"file": nm, "heads": heads, "para": para}
    return None


def _git(root):
    if not (Path(root) / ".git").exists():
        return None
    info = {}
    head = read(Path(root) / ".git" / "HEAD").strip()
    info["branch"] = head.split("/")[-1] if head.startswith("ref:") else head[:10]
    try:
        r = subprocess.run(["git", "-C", str(root), "log", "--oneline", "-8", "--no-color"],
                           capture_output=True, text=True, timeout=5, encoding="utf-8", errors="replace")
        info["log"] = r.stdout.strip().splitlines()[:8]
        r = subprocess.run(["git", "-C", str(root), "status", "--porcelain"],
                           capture_output=True, text=True, timeout=5, encoding="utf-8", errors="replace")
        info["dirty"] = [ln for ln in r.stdout.splitlines() if ln.strip()]
        r = subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"],
                           capture_output=True, text=True, timeout=5)
        info["head"] = r.stdout.strip()
    except Exception:
        pass
    return info


def scan(root, write_files=True, budget_s=10.0):
    root = Path(root)
    t0 = time.monotonic()
    files, data_dirs, truncated = walk(root, budget_s=budget_s * 0.5)
    fps = fingerprints(files)
    groups = {}
    langs = {}
    symbols = []  # (group, rel, loc, desc)
    entry, tests, todos, big = [], 0, [], []
    n_read = 0
    for rel, size, ext in sorted(files):
        g = group_of(rel)
        gi = groups.setdefault(g, {"files": 0, "src": 0, "loc": 0, "langs": {}})
        gi["files"] += 1
        lang = LANG.get(ext)
        name = rel.rsplit("/", 1)[-1]
        if re.match(r"(test_.*\.py|.*_test\.py|.*\.(test|spec)\.[jt]sx?)$", name):
            tests += 1
        if not lang:
            continue
        gi["src"] += 1
        gi["langs"][lang] = gi["langs"].get(lang, 0) + 1
        loc = max(1, size // 38)
        can_read = n_read < 700 and time.monotonic() - t0 < budget_s and (
            (ext == ".ipynb" and size < 8_000_000) or (ext != ".ipynb" and size < 400_000))
        desc = ""
        if can_read and ext in (".py", ".ipynb", ".r", ".rmd", ".js", ".jsx", ".mjs", ".ts", ".tsx"):
            text = read(root / rel)
            n_read += 1
            if ext == ".ipynb":
                nb = _nb_info(text)
                if nb:
                    loc = nb["code"] * 12
                    desc = (f"nb {nb['code']}code/{nb['md']}md"
                            + (f" | {'; '.join(nb['heads'][:5])}" if nb["heads"] else "")
                            + (f" | imports {', '.join(nb['imports'][:8])}" if nb["imports"] else ""))
            else:
                loc = text.count("\n") + 1
                c = len(TODO_RE.findall(text))
                if c:
                    todos.append((c, rel))
                if ext == ".py":
                    s = _py_symbols(text)
                    if s:
                        if s["main"]:
                            entry.append(f"{rel} (__main__)")
                        desc = (f'"{clip(s["doc"], 70)}" | ' if s["doc"] else "") + "; ".join(s["items"][:10])
                        if len(s["items"]) > 10:
                            desc += f"; +{len(s['items']) - 10}"
                else:
                    syms = _regex_symbols(text, ext)
                    desc = "; ".join(syms)
        gi["loc"] += loc
        langs.setdefault(lang, [0, 0])
        langs[lang][0] += 1
        langs[lang][1] += loc
        symbols.append((g, rel, loc, desc))
        if loc > 800 and lang not in NON_CODE_LANGS:
            big.append((loc, rel))
        if "/" not in rel and not any(x.startswith(rel + " ") for x in entry) and re.match(r"(main|app|train|run|cli|serve|manage|pipeline|benchmark|bench|evaluate|eval)\w*\.(py|r|sh|jl|js|ts)$", name.lower()):
            entry.append(rel)
    deps = _deps(root)
    pp = deps.get("pyproject") or {}
    for k, v in (pp.get("scripts") or {}).items():
        entry.append(f"pyproject script `{k}` -> {v}")
    for k, v in ((deps.get("npm") or {}).get("scripts") or {}).items():
        entry.append(f"npm run {k}: {clip(str(v), 60)}")
    mk = Path(root) / "Makefile"
    if mk.exists():
        tg = [t for t in re.findall(r"^([A-Za-z0-9_.\-]+):(?!=)", read(mk), re.M) if not t.startswith(".")][:12]
        if tg:
            entry.append("make: " + ", ".join(tg))
    for nm in ("Snakefile", "main.nf", "Dockerfile", "docker-compose.yml"):
        if (Path(root) / nm).exists():
            entry.append(nm)
    entry = list(dict.fromkeys(entry))
    code_src = sum(v[0] for k, v in langs.items() if k not in NON_CODE_LANGS)
    code_loc = sum(v[1] for k, v in langs.items() if k not in NON_CODE_LANGS)
    modules = sorted(((g, i["loc"]) for g, i in groups.items() if i["src"]), key=lambda x: -x[1])
    if code_src <= 60 and code_loc <= 15000:
        size_class, n_groups = "S", 1
    elif code_src <= 300 and code_loc <= 60000:
        size_class, n_groups = "M", 3
    else:
        size_class, n_groups = "L", 5
    bins = [[] for _ in range(min(n_groups, max(1, len(modules))))]
    loads = [0] * len(bins)
    low_info = {"tests", "test", "docs", "doc", "examples", "example", "notebooks", "scripts", "audit"}
    for g, loc in sorted(modules, key=lambda x: -(x[1] * (0.3 if x[0].lower() in low_info else 1))):
        i = loads.index(min(loads))
        bins[i].append(g)
        loads[i] += loc * (0.3 if g.lower() in low_info else 1)
    pairs = sorted(((w, b) for w, b in zip(loads, bins) if b), key=lambda x: x[0])
    while len(pairs) > 1 and pairs[0][0] < 0.25 * pairs[-1][0]:  # fold tiny bins together: fewer mapper calls
        (w0, b0), (w1, b1) = pairs[0], pairs[1]
        pairs = sorted([(w0 + w1, b1 + b0)] + pairs[2:], key=lambda x: x[0])
    bins = [b for _, b in sorted(pairs, key=lambda x: -x[0])]
    if size_class == "S":
        plan = "1x pos-mapper (cheap tier: haiku / mini) on the whole project"
    elif size_class == "M":
        plan = (f"{len(bins)}x pos-mapper (cheap tier) in parallel, one per group "
                f"{' | '.join(','.join(b) for b in bins)}; then 1x pos-mapper merges parts into MAP.md")
    else:
        plan = (f"{len(bins)}x pos-mapper (cheap tier) in parallel per group "
                f"{' | '.join(','.join(b) for b in bins)}; then 1x pos-architect (mid tier: sonnet) "
                f"synthesizes MAP.md from the parts only")
    result = {
        "scanned_at": iso(), "files": len(files), "truncated": truncated, "fp": fps,
        "code_src": code_src, "code_loc": code_loc, "size_class": size_class, "plan": plan,
        "groups_plan": bins, "git": _git(root),
    }
    if write_files:
        md = _render_scan(root, result, groups, langs, symbols, data_dirs, deps, entry, tests, todos, big)
        write(pdir(root) / "scan.md", md)
        jsave(pdir(root) / ".runtime" / "scan.json", result)
    return result


def _render_scan(root, r, groups, langs, symbols, data_dirs, deps, entry, tests, todos, big):
    L = [f"# Scan - {Path(root).name}",
         f"generated {r['scanned_at']} by `pos scan` (tier-0, no LLM). files {r['files']} | code files "
         f"{r['code_src']} | ~{r['code_loc']} LOC | truncated: {'YES' if r['truncated'] else 'no'}",
         f"size class **{r['size_class']}** -> mapping plan: {r['plan']}", "", "## Layout",
         "| group | files | code | ~LOC | langs |", "|---|---|---|---|---|"]
    for g, i in sorted(groups.items(), key=lambda x: -x[1]["loc"])[:30]:
        lg = ", ".join(k for k, _ in sorted(i["langs"].items(), key=lambda x: -x[1])[:3])
        L.append(f"| {g} | {i['files']} | {i['src']} | {i['loc']} | {lg} |")
    if data_dirs:
        L.append("\nData/artifact dirs (NOT scanned - never read wholesale): " +
                 ", ".join(f"{d}/ ({n} entries)" for d, n in data_dirs[:20]))
    L.append("\n## Languages")
    L.append("; ".join(f"{k} {v[0]} files ~{v[1]} LOC" for k, v in sorted(langs.items(), key=lambda x: -x[1][1])))
    rd = _readme(root)
    if rd:
        L += [f"\n## README ({rd['file']})", "headings: " + " / ".join(rd["heads"]), "intro: " + rd["para"]]
    if deps:
        L.append("\n## Dependencies")
        for k, v in deps.items():
            if isinstance(v, dict):
                v = {kk: vv for kk, vv in v.items() if vv}
            L.append(f"- {k}: {clip(json.dumps(v, ensure_ascii=False), 700)}")
    if entry:
        L.append("\n## Entry points")
        L += [f"- {e}" for e in entry[:25]]
    L.append(f"\n## Tests\n{tests} test files" + ("" if tests else " (none detected)"))
    g = r.get("git")
    if g:
        L.append(f"\n## Git\nbranch {g.get('branch')} | uncommitted {len(g.get('dirty', []))}")
        L += [f"- {x}" for x in g.get("log", [])]
    ctx = [n for n in ("AGENTS.md", "CLAUDE.md", "GEMINI.md", "graphify-out/GRAPH_REPORT.md", "CHANGELOG.md", "docs")
           if (Path(root) / n).exists()]
    if ctx:
        L.append("\n## Existing context files\n" + ", ".join(ctx) +
                 (" (graphify report: read its summary / god-nodes first)" if "graphify-out/GRAPH_REPORT.md" in ctx else ""))
    if big:
        L.append("\n## Hotspots\nlargest code files: " + ", ".join(f"{rel} ({loc})" for loc, rel in sorted(big, reverse=True)[:8]))
    if todos:
        L.append(f"TODO/FIXME: {sum(c for c, _ in todos)} total; top: " +
                 ", ".join(f"{rel} ({c})" for c, rel in sorted(todos, reverse=True)[:5]))
    L.append("\n## Symbol index (top-level defs; notebooks: headings/imports)")
    budget = 16000 - sum(len(x) for x in L)
    cur_g = None
    shown = 0
    order = sorted(symbols, key=lambda s: (-(groups[s[0]]["loc"]), s[0], s[1]))
    for gname, rel, loc, desc in order:
        line = f"- {rel} ({loc})" + (f" - {clip(desc, 170)}" if desc else "")
        head = f"### {gname}" if gname != cur_g else None
        need = len(line) + (len(head) + 1 if head else 0) + 1
        if budget - need < 0:
            break
        if head:
            L.append(head)
            cur_g = gname
        L.append(line)
        budget -= need
        shown += 1
    if shown < len(symbols):
        L.append(f"... {len(symbols) - shown} more code files not listed (use Grep/Glob).")
    return "\n".join(L) + "\n"


def map_lock(root):
    return jload(pdir(root) / "map.lock.json", None)


def map_staleness(root, budget_s=2.5):
    lock = map_lock(root)
    if not lock:
        return {"status": "unmapped"}
    files, _, trunc = walk(root, budget_s=budget_s)
    if trunc:
        return {"status": "unknown", "reason": "walk budget exceeded", **{k: lock.get(k) for k in ("mapped_at", "mapped_by", "confidence")}}
    cur = fingerprints(files)
    old = lock.get("fp", {})
    changed = sorted(g for g in set(cur) | set(old) if cur.get(g) != old.get(g))
    return {"status": "stale" if changed else "fresh", "changed": changed,
            **{k: lock.get(k) for k in ("mapped_at", "mapped_by", "confidence")}}


# --------------------------------------------------------------------------- markdown helpers

STATUSES = ["RUNNING", "SUCCESS", "PARTIAL", "FAILED", "ABANDONED", "INCONCLUSIVE"]
HEAD_RE = re.compile(r"^###\s+([AD]-\d+)\s*\|([^|]*)\|([^|]*)\|([^|]*)\|\s*(.+?)\s*$")


def md_sections(text):
    secs, cur, buf = {}, None, []
    for line in text.splitlines():
        m = re.match(r"^##\s+(.+?)\s*$", line)
        if m:
            if cur is not None:
                secs[cur] = "\n".join(buf)
            cur, buf = m.group(1).strip(), []
        elif cur is not None:
            buf.append(line)
    if cur is not None:
        secs[cur] = "\n".join(buf)
    return secs


def clean_md(s):
    s = re.sub(r"<!--.*?-->", "", s or "", flags=re.S)
    keep = []
    for ln in s.splitlines():
        t = ln.strip()
        if not t or t in ("-", "*", "1.", "- [ ]", "- [x]"):
            continue
        keep.append(ln.rstrip())
    return "\n".join(keep)


def section_digest(text, wanted, total=1200):
    secs = md_sections(text)
    out = []
    per = max(200, total // max(1, len(wanted)))
    for w in wanted:
        for name, body in secs.items():
            if name.lower().startswith(w):
                b = clean_md(body)
                if b:
                    out.append(f"**{name}**: " + (b if "\n" not in b else "\n" + b)[:per])
                break
    return clip("\n".join(out), total)


def parse_entries(text):
    entries, cur = [], None
    for line in text.splitlines():
        m = HEAD_RE.match(line.rstrip())
        if m:
            cur = {"id": m.group(1), "date": m.group(2).strip(), "agent": m.group(3).strip(),
                   "status": m.group(4).strip(), "title": m.group(5).strip(), "fields": {}, "raw": [line]}
            entries.append(cur)
            continue
        if cur is None:
            continue
        if line.startswith("#"):
            cur = None
            continue
        cur["raw"].append(line)
        fm = re.match(r"^\s*-\s*([A-Za-z ]+?):\s*(.*)$", line)
        if fm:
            k = fm.group(1).strip().lower()
            cur["fields"][k] = (cur["fields"].get(k, "") + " " + fm.group(2)).strip()
    return entries


def next_id(entries, prefix):
    nums = [int(e["id"].split("-")[1]) for e in entries if e["id"].startswith(prefix)]
    return f"{prefix}-{(max(nums) + 1) if nums else 1:03d}"


# --------------------------------------------------------------------------- templates

T_README = """# .projectos - project memory (ProjectOS)

Shared, agent-agnostic memory for this project (Claude Code, Codex, Antigravity, humans).
Read `PROTOCOL.md` for the rules. Files:

| file | role | who edits |
|---|---|---|
| PURPOSE.md | why / done-when / constraints / non-goals | user (agents propose) |
| MAP.md | architecture map so nobody re-reads the codebase | cheap mapper model |
| STATE.md | live focus, next steps, blockers | every agent, continuously |
| ATTEMPTS.md | append-only log of approaches tried and their outcome | every agent |
| DECISIONS.md | append-only design decisions (ADR-lite) | every agent |
| checkpoints/ | immutable snapshots (manual + auto) | `pos checkpoint` / hooks |
| scan.md | tier-0 mechanical scan (no LLM) | `pos scan` |
| .runtime/ | machine-local caches (gitignored) | hooks |
"""

T_PURPOSE = """# Purpose - {{name}}
<!-- pos:template  (remove this line once the sections below are filled and confirmed by the user) -->
> Owner: the user. Agents may draft, but changes need user confirmation and a DECISIONS.md entry.

## Why this project exists
{{draft}}

## Goal (definition of done)
-

## Success criteria (measurable)
- [ ]

## Constraints (time, compute, data, tools, style)
-

## Non-goals
-

## Context (deadline, venue, stakeholders, links)
-
"""

T_STATE = """# State - {{name}}
Updated: {{date}} by {{agent}}

## Focus
<!-- what is being worked on right now, 1-3 lines -->

## Done recently
-

## In progress
-

## Next steps
1.

## Blockers / open questions
-

## Notes for the next agent
-
"""

T_ATTEMPTS = """# Attempts - {{name}}
Append-only. One entry per approach/experiment. Never delete; update status instead.
Format (keep the heading exactly, `pos attempt` parses it):

    ### A-001 | YYYY-MM-DD | agent | STATUS | short title
    - goal: what it should achieve
    - approach: what was done (config, params, commit)
    - result: outcome + metric
    - evidence: log/file/commit path
    - lesson: what to do or avoid next time
    - tags: kw1, kw2

STATUS: RUNNING | SUCCESS | PARTIAL | FAILED | ABANDONED | INCONCLUSIVE

"""

T_DECISIONS = """# Decisions - {{name}}
Append-only ADR-lite. Supersede instead of editing history.

    ### D-001 | YYYY-MM-DD | agent | ACTIVE | short title
    - decision: what was chosen
    - why: reasons / evidence
    - alternatives: what was rejected and why
    - supersedes: D-000 (optional)

STATUS: PROPOSED | ACTIVE | SUPERSEDED by D-xxx | REVERTED

"""

T_MAP = """---
status: unmapped
---
# Map - {{name}}
Not mapped yet. Tier-0 facts are in `scan.md`. Mapping: skill `pos-map` (cheap model only).
"""

ENTRY_BEGIN = "<!-- projectos:begin"
ENTRY_END = "<!-- projectos:end -->"
ENTRY_BLOCK = """<!-- projectos:begin v1 -->
## ProjectOS (read first)
This project keeps shared agent memory in `.projectos/` (purpose, map, state, attempts, decisions, checkpoints) for every agent: Claude, Codex, Antigravity, humans.
1. Start: read `.projectos/STATE.md` and `PURPOSE.md` (or `pos brief`). Hooks may have injected this already.
2. Code questions: `.projectos/MAP.md` first, raw files last.
3. Before a new approach: `pos attempt list --grep <kw>`. Never repeat a FAILED attempt without stating what differs.
4. Log outcomes (`pos attempt add`, `pos decide`). After each finished sub-task: update STATE.md, then `pos checkpoint -m "..."`.
Rules: `.projectos/PROTOCOL.md`. No `pos` on this machine? `python .projectos/kit/install.py` if `.projectos/kit/` exists (full setup), else edit the markdown files directly, same formats.
<!-- projectos:end -->"""


def fill(t, **kw):
    for k, v in kw.items():
        t = t.replace("{{" + k + "}}", str(v))
    return t


def upsert_block(path, block):
    path = Path(path)
    text = read(path, "") if path.exists() else ""
    if ENTRY_BEGIN in text and ENTRY_END in text:
        new = re.sub(re.escape(ENTRY_BEGIN) + r".*?" + re.escape(ENTRY_END), lambda _: block, text, flags=re.S)
    elif "## graphify-context" in text:
        i = text.index("## graphify-context")
        new = text[:i] + block + "\n\n" + text[i:]
    elif text.strip():
        new = text.rstrip() + "\n\n" + block + "\n"
    else:
        new = block + "\n"
    if new != text:
        write(path, new)
        return True
    return False


# --------------------------------------------------------------------------- init / upgrade

def install_entry_files(root):
    upsert_block(Path(root) / "AGENTS.md", ENTRY_BLOCK)
    upsert_block(Path(root) / "CLAUDE.md", ENTRY_BLOCK)
    rule = "---\ndescription: ProjectOS shared project memory protocol\nalwaysApply: true\n---\n" + ENTRY_BLOCK + "\n"
    rp = Path(root) / ".agent" / "rules" / "projectos.md"
    if read(rp) != rule:
        write(rp, rule)


def copy_protocol(root):
    src = ENGINE_DIR / "PROTOCOL.md"
    if src.exists():
        write(pdir(root) / "PROTOCOL.md", read(src))


def init_project(root, agent="agent", do_scan=True, graphify=True, quiet=False, graphify_wait=False):
    root = Path(root).resolve()
    pd = pdir(root)
    new = not (pd / "project.json").exists()
    (pd / "checkpoints").mkdir(parents=True, exist_ok=True)
    m = meta(root)
    if new:
        uid = uuid.uuid4().hex
        m = {"id": unique_id(root, uid), "uid": uid, "name": root.name, "created": iso(), "created_by": agent}
    m.update({"protocol_version": PROTOCOL_VERSION, "engine_version": VERSION})
    jsave(pd / "project.json", m)
    copy_protocol(root)
    write(pd / ".gitignore", ".runtime/\n*.tmp\n")
    write(pd / "README.md", T_README)
    rd = _readme(root)
    draft = (f"<!-- draft from {rd['file']}, confirm with user: --> {rd['para']}" if rd and rd["para"] else "")
    kw = {"name": root.name, "date": iso(), "agent": agent, "draft": draft}
    for fname, tpl in (("PURPOSE.md", T_PURPOSE), ("STATE.md", T_STATE), ("ATTEMPTS.md", T_ATTEMPTS),
                       ("DECISIONS.md", T_DECISIONS), ("MAP.md", T_MAP)):
        if not (pd / fname).exists():
            write(pd / fname, fill(tpl, **kw))
    install_entry_files(root)
    res = scan(root) if do_scan else None
    g = maybe_graphify(root, wait=graphify_wait) if graphify else "skipped"
    register(root, {"map": "unmapped" if not map_lock(root) else "mapped"})
    if new:
        ledger(root, "projectos-init", f"size={res and res['size_class']} graphify={g}")
    if not quiet:
        print(f"ProjectOS {'initialized' if new else 'refreshed'}: {root}")
        print(f"  id={m['id']}  graphify={g}")
        if res:
            print(f"  scan: {res['files']} files, {res['code_src']} code files, ~{res['code_loc']} LOC, class {res['size_class']}")
            print(f"  mapping plan: {res['plan']}")
    return new, res


# --------------------------------------------------------------------------- checkpoints

CK_RE = re.compile(r"^(\d{4}-\d{2}-\d{2})_(\d{6})_")


def list_checkpoints(root):
    cd = pdir(root) / "checkpoints"
    try:
        names = sorted(n for n in os.listdir(cd) if CK_RE.match(n) and n.endswith(".md"))
    except Exception:
        names = []
    return names


def ck_time(name):
    m = CK_RE.match(name)
    if not m:
        return None
    try:
        return dt.datetime.strptime(m.group(1) + m.group(2), "%Y-%m-%d%H%M%S").timestamp()
    except Exception:
        return None


def ck_front(root, name):
    t = read(pdir(root) / "checkpoints" / name)[:1500]
    fm = {}
    m = re.match(r"^---\n(.*?)\n---", t, re.S)
    if m:
        for ln in m.group(1).splitlines():
            if ":" in ln:
                k, v = ln.split(":", 1)
                fm[k.strip()] = v.strip()
    return fm


def last_checkpoint_ts(root):
    cks = list_checkpoints(root)
    return ck_time(cks[-1]) if cks else None


def tail_events(path, max_bytes=3_000_000):
    out = []
    try:
        with open(path, "rb") as f:
            f.seek(0, 2)
            size = f.tell()
            f.seek(max(0, size - max_bytes))
            data = f.read()
    except Exception:
        return out
    lines = data.split(b"\n")
    if size > max_bytes:
        lines = lines[1:]
    for ln in lines:
        if not ln.strip():
            continue
        try:
            out.append(json.loads(ln))
        except Exception:
            pass
    return out


EDIT_TOOLS = {"Edit", "Write", "MultiEdit", "NotebookEdit"}
PATCH_FILE_RE = re.compile(r"\*\*\* (?:Update|Add|Delete) File: ([^\r\n]+)")


def _strings(o):
    if isinstance(o, str):
        yield o
    elif isinstance(o, dict):
        for v in o.values():
            yield from _strings(v)
    elif isinstance(o, list):
        for v in o:
            yield from _strings(v)


def _under(root, f):
    p = Path(f)
    if not p.is_absolute():
        return True  # relative paths (e.g. Codex apply_patch) are project-relative
    try:
        p.resolve().relative_to(Path(root).resolve())
        return True
    except Exception:
        return False


def transcript_digest(agent, transcript, since=None, root=None):
    """Mechanical digest of a Claude/Codex transcript: edited files (inside root if given), prompts, last reply."""
    files, prompts, last = {}, [], ""
    if not transcript or not Path(transcript).exists():
        return {"files": files, "prompts": prompts, "last": last}
    dg = _digest(agent, transcript, since)
    if root is not None:
        dg["files"] = {f: t for f, t in dg["files"].items() if _under(root, f)}
    return dg


def _digest(agent, transcript, since):
    files, prompts, last = {}, [], ""
    for e in tail_events(transcript):
        ts = parse_ts(e.get("timestamp"))
        if since and ts and ts < since:
            continue
        if agent == "codex" or "payload" in e:
            p = e.get("payload") or {}
            pt = p.get("type")
            if pt in ("function_call", "custom_tool_call", "local_shell_call"):
                for blob in _strings(p):
                    if "*** " in blob:
                        for mm in PATCH_FILE_RE.finditer(blob):
                            files[mm.group(1).strip()] = "apply_patch"
            elif pt == "user_message" and p.get("message"):
                prompts.append(str(p["message"]))
            elif pt == "agent_message" and p.get("message"):
                last = str(p["message"])
        else:
            t = e.get("type")
            msg = e.get("message") or {}
            if t == "assistant":
                for c in msg.get("content") or []:
                    if not isinstance(c, dict):
                        continue
                    if c.get("type") == "tool_use" and c.get("name") in EDIT_TOOLS:
                        inp = c.get("input") or {}
                        fp = inp.get("file_path") or inp.get("notebook_path")
                        if fp:
                            files[fp] = c["name"]
                    elif c.get("type") == "text" and c.get("text"):
                        last = c["text"]
            elif t == "user" and not e.get("isMeta") and isinstance(msg.get("content"), str):
                s = msg["content"].strip()
                if s and not s.startswith("<"):
                    prompts.append(s)
    return {"files": files, "prompts": prompts[-3:], "last": last}


def rel_to(root, p):
    try:
        return Path(p).resolve().relative_to(Path(root).resolve()).as_posix()
    except Exception:
        return str(p)


def write_checkpoint(root, agent, summary, kind="manual", reason="milestone", next_steps=None,
                     advances=None, digest=None, usage_s=None):
    root = Path(root)
    pd = pdir(root)
    t = now()
    name = f"{t.strftime('%Y-%m-%d_%H%M%S')}_{slugify(agent, 12)}_{slugify(summary, 40)}.md"
    L = ["---", f"time: {iso(t)}", f"agent: {agent}", f"kind: {kind}", f"reason: {reason}",
         f"summary: {one_line(summary, 200)}", f"usage: {usage_s or 'n/a'}", f"machine: {MACHINE}", "---",
         f"# Checkpoint - {one_line(summary, 120)}", ""]
    if advances:
        L += [f"**Advances purpose:** {advances}", ""]
    if next_steps:
        L += ["## Next", next_steps, ""]
    if digest and (digest.get("files") or digest.get("prompts")):
        if digest.get("files"):
            L.append(f"## Files edited since previous checkpoint ({len(digest['files'])})")
            L += [f"- {rel_to(root, f)} ({tool})" for f, tool in list(digest["files"].items())[:60]]
            L.append("")
        if digest.get("prompts"):
            L.append("## Recent user requests")
            L += [f"- {one_line(p, 300)}" for p in digest["prompts"]]
            L.append("")
        if digest.get("last") and kind == "auto":
            L += ["## Last agent message (excerpt)", clip(digest["last"], 900), ""]
    g = _git(root) if kind == "manual" else None
    if g and g.get("dirty"):
        L.append(f"## Uncommitted changes ({len(g['dirty'])})")
        L += [f"- {x}" for x in g["dirty"][:25]]
        L.append("")
    state = read(pd / "STATE.md")
    L += ["## STATE.md snapshot", "", ("\n" + state).replace("\n# ", "\n#### ").replace("\n## ", "\n#### ").lstrip("\n")]
    write(pd / "checkpoints" / name, "\n".join(L) + "\n")
    rotate_checkpoints(root)
    rt = jload(pd / ".runtime" / "state.json", {})
    rt.update({"last_checkpoint": name, "last_checkpoint_ts": epoch(), "last_kind": kind})
    jsave(pd / ".runtime" / "state.json", rt)
    if kind == "manual":
        register(root, {"last_checkpoint": {"file": name, "time": iso(t), "summary": one_line(summary, 160),
                                            "agent": agent}, "focus": one_line(section_digest(state, ["focus"], 200), 200)})
        ledger(root, "projectos-checkpoint", f"{agent}: {summary}")
    return name


def rotate_checkpoints(root):
    keep = int(cfg().get("keep_checkpoints", 40))
    cks = list_checkpoints(root)
    if len(cks) <= keep:
        return
    cd = pdir(root) / "checkpoints"
    lines = []
    for n in cks[: len(cks) - keep]:
        fm = ck_front(root, n)
        lines.append(f"- {n}: [{fm.get('kind', '?')}/{fm.get('reason', '?')}] {fm.get('agent', '?')} - {fm.get('summary', '')}")
        try:
            (cd / n).unlink()
        except Exception:
            pass
    if lines:
        if not (cd / "_archive.md").exists():
            append(cd / "_archive.md", "# Archived checkpoints (condensed)\n")
        append(cd / "_archive.md", "\n".join(lines) + "\n")


# --------------------------------------------------------------------------- usage detection

def latest_codex_rollout():
    base = HOME / ".codex" / "sessions"
    try:
        days = sorted(base.glob("[0-9][0-9][0-9][0-9]/[0-9][0-9]/[0-9][0-9]"), key=lambda x: x.as_posix())[-3:]
        files = [f for d in days for f in d.glob("*.jsonl")]
        return max(files, key=lambda x: x.stat().st_mtime) if files else None
    except Exception:
        return None


def usage(agent, sid=None, transcript=None):
    u = {"context_pct": None, "five_hour_pct": None, "seven_day_pct": None}
    t = epoch()
    if agent == "claude":
        ud = RUNTIME / "usage"
        rec = jload(ud / f"claude-{slugify(sid, 80)}.json", {}) if sid else {}
        latest = jload(ud / "claude-latest.json", {})
        if not rec and not sid:
            rec = latest
        if rec and t - rec.get("ts", 0) < 900:
            u["context_pct"] = rec.get("context_pct")
        rl = rec if rec.get("five_hour_pct") is not None else latest
        if rl and t - rl.get("ts", 0) < 3600:
            if rl.get("five_hour_pct") is not None and (rl.get("five_hour_reset") or t + 1) > t:
                u["five_hour_pct"] = rl["five_hour_pct"]
            if rl.get("seven_day_pct") is not None and (rl.get("seven_day_reset") or t + 1) > t:
                u["seven_day_pct"] = rl["seven_day_pct"]
        if u["context_pct"] is None and transcript:
            # Transcript fallback (e.g. desktop app, no statusline). Only trust a % when the window size is
            # known; otherwise report tokens only - PreCompact auto-checkpoint still covers compaction.
            for e in reversed(tail_events(transcript, 600_000)):
                msg = e.get("message") or {}
                us = msg.get("usage") if isinstance(msg, dict) else None
                if e.get("type") == "assistant" and us:
                    model = str(msg.get("model", ""))
                    tok = (us.get("input_tokens") or 0) + (us.get("cache_creation_input_tokens") or 0) + (us.get("cache_read_input_tokens") or 0)
                    u["context_tokens"] = tok
                    window = (rec or {}).get("context_size")
                    if not window:
                        window = 1_000_000 if "[1m]" in model else next(
                            (v for k, v in (cfg().get("context_windows") or {}).items() if model.startswith(k)), None)
                    if window and tok <= window:
                        u["context_pct"] = round(100 * tok / window, 1)
                    break
    elif agent == "codex":
        tp = transcript or latest_codex_rollout()
        if tp:
            for e in reversed(tail_events(tp, 800_000)):
                p = e.get("payload") or {}
                if p.get("type") != "token_count":
                    continue
                info = p.get("info") or {}
                last = info.get("last_token_usage") or {}
                win = info.get("model_context_window")
                if win and last.get("input_tokens") is not None:
                    u["context_pct"] = round(100 * last["input_tokens"] / win, 1)
                rl = p.get("rate_limits") or {}
                if (rl.get("primary") or {}).get("used_percent") is not None:
                    u["five_hour_pct"] = rl["primary"]["used_percent"]
                if (rl.get("secondary") or {}).get("used_percent") is not None:
                    u["seven_day_pct"] = rl["secondary"]["used_percent"]
                break
    vals = [u[k] for k in ("context_pct", "five_hour_pct", "seven_day_pct") if isinstance(u.get(k), (int, float))]
    u["max"] = max(vals) if vals else None
    return u


def usage_str(u):
    parts = []
    for k, lab in (("context_pct", "ctx"), ("five_hour_pct", "5h"), ("seven_day_pct", "7d")):
        if u.get(k) is not None:
            parts.append(f"{lab} {u[k]:.0f}%")
    if u.get("context_pct") is None and u.get("context_tokens"):
        parts.append(f"ctx {u['context_tokens'] // 1000}k tok")
    return " | ".join(parts) if parts else "n/a"


# --------------------------------------------------------------------------- sessions (multi-agent overlap)

def sessions_touch(root, sid, agent, transcript=None, task=None, end=False, force=False):
    p = pdir(root) / ".runtime" / "sessions.json"
    data = jload(p, {})
    t = epoch()
    horizon = int(cfg()["session_active_min"]) * 60 * 4
    data = {k: v for k, v in data.items() if t - v.get("heartbeat", 0) < horizon}
    if end:
        data.pop(sid, None)
        jsave(p, data)
        return
    rec = data.get(sid) or {"agent": agent, "started": t, "machine": MACHINE}
    if not force and not task and sid in data and t - rec.get("heartbeat", 0) < 300:
        return
    rec["heartbeat"] = t
    if transcript:
        rec["transcript"] = transcript
    if task is not None:
        rec["task"] = task
    data[sid] = rec
    jsave(p, data)


def other_sessions(root, sid):
    data = jload(pdir(root) / ".runtime" / "sessions.json", {})
    t = epoch()
    act = int(cfg()["session_active_min"]) * 60
    return [(k, v) for k, v in data.items() if k != sid and t - v.get("heartbeat", 0) < act]


def own_session(root, agent):
    """Best guess of the calling agent's session (for CLI calls made from inside an agent)."""
    sid = env_session_id()
    data = jload(pdir(root) / ".runtime" / "sessions.json", {})
    if sid and sid in data:
        return sid, data[sid]
    mine = sorted(((k, v) for k, v in data.items() if v.get("agent") == agent), key=lambda kv: -kv[1].get("heartbeat", 0))
    return mine[0] if mine else (sid, {})


def sess_state(sid):
    return jload(RUNTIME / "sessions" / f"{slugify(sid or 'none', 80)}.json", {})


def sess_save(sid, st):
    jsave(RUNTIME / "sessions" / f"{slugify(sid or 'none', 80)}.json", st)


# --------------------------------------------------------------------------- brief

def map_line(root, agent, st=None):
    st = st or map_staleness(root)
    sc = jload(pdir(root) / ".runtime" / "scan.json", {})
    if st["status"] == "unmapped":
        plan = sc.get("plan", "run `pos scan` first")
        if agent == "claude":
            how = ("Dispatch now, in background: Agent(subagent_type='pos-mapper', run_in_background=true, "
                   f"prompt='Map project at {root}' + group list if plan has groups). pos-mapper runs on haiku. "
                   "Keep working on the user's request meanwhile using scan.md. Never map with Opus/main model.")
        else:
            how = ("Follow skill `pos-map`: delegate to a sub-agent on the cheapest model if you can "
                   "(e.g. a *-mini model), else do it yourself within the read budget. Never re-read the whole repo.")
        return f"MAP MISSING (size {sc.get('size_class', '?')}; plan: {plan}).\n{how}"
    info = f"mapped {st.get('mapped_at', '?')} by {st.get('mapped_by', '?')}, confidence {st.get('confidence', '?')}"
    if st["status"] == "stale":
        return (f"MAP STALE in: {', '.join(st['changed'][:12])} ({info}). Incremental remap of ONLY those groups "
                f"via pos-mapper (cheap tier) when you next need them; read MAP.md for the rest.")
    if st["status"] == "unknown":
        return f"MAP present ({info}); staleness unknown ({st.get('reason')}). Read .projectos/MAP.md before exploring code."
    return f"MAP fresh ({info}). Read .projectos/MAP.md before exploring code."


def brief(root, agent="agent", sid=None, transcript=None, fresh=False, source="startup", max_chars=None, staleness=None):
    root = Path(root)
    max_chars = max_chars or int(cfg()["brief_max_chars"])
    pd = pdir(root)
    m = meta(root)
    L = [f"# ProjectOS brief - {m.get('name', root.name)} [{m.get('id', '?')}]",
         f"root: {root} | agent: {agent} | protocol v{PROTOCOL_VERSION} | {iso()}"]
    if m.get("protocol_version", 0) < PROTOCOL_VERSION:
        L.append("(protocol outdated in this project: run `pos upgrade`)")
    if source == "compact":
        L.append("!! Context was just compacted; an auto-checkpoint captured edited files. Re-read .projectos/STATE.md; "
                 "if stale, update it and run `pos checkpoint -m ...`.")
    if fresh:
        L.append("!! Project just initialized by ProjectOS (tier-0 scan, zero LLM tokens).")
    purpose = read(pd / "PURPOSE.md")
    if "<!-- pos:draft" in purpose:
        L += ["\n## Purpose (DRAFT by an agent - confirm with the user, then delete the draft marker line)",
              section_digest(purpose, ["why", "goal", "success", "constraint"], 1100)]
    elif "<!-- pos:template" in purpose:
        L.append("\n## Purpose: NOT DEFINED\nAsk the user 2-3 short questions (why / done-when / constraints) at a natural "
                 "moment, or draft from README and get confirmation; fill PURPOSE.md, delete the template marker line.")
    else:
        L += ["\n## Purpose", section_digest(purpose, ["why", "goal", "success", "constraint"], 1100)]
    state = read(pd / "STATE.md")
    upd = re.search(r"^Updated:\s*(.+)$", state, re.M)
    L += [f"\n## State (STATE.md{', updated ' + upd.group(1).strip() if upd else ''})",
          section_digest(state, ["focus", "in progress", "next", "blocker", "notes for"], 1700) or "(empty)"]
    cks = list_checkpoints(root)
    if cks:
        manual = [c for c in cks[-12:] if ck_front(root, c).get("kind") != "auto"]
        lastm = manual[-1] if manual else None
        after_auto = len(cks) - 1 - cks.index(lastm) if lastm else len(cks)
        if lastm:
            fm = ck_front(root, lastm)
            L.append(f"\n## Last checkpoint: {lastm} ({ago(epoch() - (ck_time(lastm) or epoch()))}, {fm.get('agent', '?')}): "
                     f"{fm.get('summary', '')}")
        if after_auto > 0:
            L.append(f"{after_auto} auto-checkpoint(s) after it (latest: checkpoints/{cks[-1]}) - STATE.md may lag; "
                     "check that file's edited-files list.")
    else:
        L.append("\n## Checkpoints: none yet")
    ents = parse_entries(read(pd / "ATTEMPTS.md"))
    if ents:
        counts = {}
        for e in ents:
            k = e["status"].split()[0].upper()
            counts[k] = counts.get(k, 0) + 1
        L.append(f"\n## Attempts ({', '.join(f'{k} {v}' for k, v in counts.items())}) - check before trying anything")
        shown = set()
        for e in [e for e in ents if e["status"].upper().startswith("RUNNING")] + ents[-5:]:
            if e["id"] in shown:
                continue
            shown.add(e["id"])
            les = e["fields"].get("lesson") or e["fields"].get("result") or ""
            L.append(f"- {e['id']} {e['status']} {one_line(e['title'], 70)}" + (f" -> {one_line(les, 110)}" if les else ""))
        failed = [e for e in ents if e["status"].upper().startswith(("FAILED", "ABANDONED")) and e["id"] not in shown][-8:]
        for e in failed:
            L.append(f"- {e['id']} {e['status']} {one_line(e['title'], 60)} -> lesson: {one_line(e['fields'].get('lesson', '?'), 100)}")
    decs = [d for d in parse_entries(read(pd / "DECISIONS.md")) if d["status"].upper().startswith(("ACTIVE", "PROPOSED"))]
    if decs:
        L.append("\n## Decisions in force")
        L += [f"- {d['id']} {d['status']}: {one_line(d['title'], 80)}" for d in decs[-8:]]
    L.append("\n## Map\n" + map_line(root, agent, staleness))
    others = other_sessions(root, sid)
    if others:
        L.append("\n## Other agents active on this project (avoid overlap: different files, coordinate via STATE.md)")
        for k, v in others:
            L.append(f"- {v.get('agent')} on {v.get('machine')} (session {k[:8]}), heartbeat {ago(epoch() - v.get('heartbeat', 0))}"
                     + (f", task: {v['task']}" if v.get("task") else ""))
    u = usage(agent, sid, transcript) if agent in ("claude", "codex") else None
    if u and u.get("max") is not None:
        L.append(f"\n## Usage: {usage_str(u)} (warn {cfg()['warn_pct']}%, forced checkpoint {cfg()['force_pct']}%)")
    L.append("\n## Rules: MAP.md before code | `pos attempt list --grep kw` before new approaches | log attempts/decisions "
             "| after each sub-task: update STATE.md + `pos checkpoint -m` | cheapest capable model for mechanical work "
             "| full protocol: .projectos/PROTOCOL.md, skills: projectos, pos-map, pos-checkpoint, pos-log, pos-resume")
    return clip("\n".join(x for x in L if x is not None), max_chars)


def workspace_brief(base):
    kids = child_projects(base, limit=200)
    mapped = [k for k in kids if (k / DIRNAME / "project.json").exists()]
    L = [f"# ProjectOS: workspace folder `{base}` ({len(kids)} sub-projects, {len(mapped)} with ProjectOS).",
         "Not initializing a project here. Open a sub-project, or map all of them cheaply: "
         "`pos discover . --init` then skill `pos-map` (mode: all)."]
    for k in kids[:20]:
        L.append(f"- {k.name}: {'ProjectOS' if k in mapped else 'unmapped'}")
    return "\n".join(L)


# --------------------------------------------------------------------------- hooks

def ctx_out(event, text):
    if not text:
        return None
    return {"hookSpecificOutput": {"hookEventName": event, "additionalContext": text}}


def block_out(reason):
    return {"decision": "block", "reason": reason}


def force_thresholds():
    f = float(cfg()["force_pct"])
    return [f, min(97.0, f + 7)]


def h_session_start(pl, cwd, sid, tp, agent):
    source = pl.get("source") or "startup"
    root = find_root(cwd)
    fresh = False
    if root is None:
        if ignored_reason(cwd):
            return None
        base = git_root(cwd) or Path(cwd).resolve()
        if ignored_reason(base):
            base = Path(cwd).resolve()
        if classify_dir(base) == "workspace":
            return ctx_out("SessionStart", workspace_brief(base))
        if not cfg().get("auto_init", True):
            return ctx_out("SessionStart", f"ProjectOS: `{base}` has no .projectos/. Run `pos init` to enable project memory.")
        init_project(base, agent=agent, quiet=True)
        root, fresh = base, True
    sessions_touch(root, sid, agent, transcript=tp, force=True)
    st = sess_state(sid)
    st.setdefault("started", epoch())
    st["project"] = str(root)
    sess_save(sid, st)
    return ctx_out("SessionStart", brief(root, agent, sid, tp, fresh=fresh, source=source))


def h_prompt(pl, cwd, sid, tp, agent):
    root = find_root(cwd)
    if not root:
        return None
    sessions_touch(root, sid, agent, transcript=tp)
    u = usage(agent, sid, tp)
    lvl = u.get("max")
    if lvl is None:
        return None
    st = sess_state(sid)
    msgs = []
    for thr in force_thresholds():
        if lvl >= thr and thr not in st.get("prompt_forced", []):
            st.setdefault("prompt_forced", []).append(thr)
            msgs.append(f"ProjectOS: usage {usage_str(u)} >= {thr:.0f}%. BEFORE other work: update .projectos/STATE.md "
                        f"(focus, in progress, next steps) and run `pos checkpoint -m \"<summary>\" --reason usage` "
                        f"(skill pos-checkpoint, emergency mode, <2 min). Then continue.")
            break
    if not msgs and lvl >= float(cfg()["warn_pct"]) and not st.get("warned"):
        st["warned"] = True
        msgs.append(f"ProjectOS: usage {usage_str(u)}. Keep STATE.md current and checkpoint at the next natural break.")
    sess_save(sid, st)
    return ctx_out("UserPromptSubmit", "\n".join(msgs)) if msgs else None


def h_stop(pl, cwd, sid, tp, agent):
    if pl.get("stop_hook_active"):
        return None
    root = find_root(cwd)
    if not root:
        return None
    st = sess_state(sid)
    t = epoch()
    if "started" not in st:
        st["started"] = t
        sess_save(sid, st)
    last_ck = max(last_checkpoint_ts(root) or 0, 0)
    u = usage(agent, sid, tp)
    lvl = u.get("max")
    if lvl is not None:
        for thr in force_thresholds():
            if lvl >= thr and thr not in st.get("stop_forced", []):
                st.setdefault("stop_forced", []).append(thr)
                sess_save(sid, st)
                if t - last_ck < 600:
                    return None
                return block_out(f"ProjectOS: usage {usage_str(u)} (>= {thr:.0f}%). Write a checkpoint now so no work is lost: "
                                 "update .projectos/STATE.md (focus / in progress / next steps / notes for next agent), log any "
                                 "concluded attempt with `pos attempt add`, then run `pos checkpoint -m \"<summary>\" --reason usage`. "
                                 "Keep it short, then stop.")
    period = int(cfg()["periodic_checkpoint_min"]) * 60
    since = max(last_ck, st.get("started", t))
    if t - since > period and t - st.get("periodic_nudged", 0) > period:
        dg = transcript_digest(agent, tp, since=last_ck or None, root=root)
        if dg["files"]:
            st["periodic_nudged"] = t
            sess_save(sid, st)
            return block_out(f"ProjectOS periodic checkpoint: {len(dg['files'])} file(s) edited and no checkpoint for "
                             f"{ago(t - since).replace(' ago', '')}. Update .projectos/STATE.md, log concluded attempts/decisions, "
                             "run `pos checkpoint -m \"<summary>\"`, then stop.")
    return None


def auto_checkpoint(root, agent, sid, tp, reason):
    last_ck = last_checkpoint_ts(root)
    dg = transcript_digest(agent, tp, since=last_ck, root=root)
    if reason == "session-end" and not dg["files"]:
        return None
    u = usage(agent, sid, tp)
    summ = f"auto ({reason}): {len(dg['files'])} file(s) edited"
    return write_checkpoint(root, agent, summ, kind="auto", reason=reason, digest=dg, usage_s=usage_str(u))


def h_precompact(pl, cwd, sid, tp, agent):
    root = find_root(cwd)
    if root:
        auto_checkpoint(root, agent, sid, tp, "precompact")
    return None


def h_session_end(pl, cwd, sid, tp, agent):
    root = find_root(cwd)
    if not root:
        return None
    sessions_touch(root, sid, agent, end=True)
    auto_checkpoint(root, agent, sid, tp, "session-end")
    return None


HOOKS = {"session-start": h_session_start, "prompt": h_prompt, "stop": h_stop,
         "precompact": h_precompact, "session-end": h_session_end}


def cmd_hook(a):
    try:
        raw = sys.stdin.read() if not sys.stdin.isatty() else ""
        pl = json.loads(raw) if raw.strip() else {}
    except Exception:
        pl = {}
    cwd = Path(pl.get("cwd") or os.getcwd())
    sid = pl.get("session_id") or env_session_id() or "nosession"
    tp = pl.get("transcript_path")
    try:
        out = HOOKS[a.event](pl, cwd, sid, tp, a.agent)
    except Exception as e:
        log_err(f"hook {a.event}", e)
        out = None
    if out:
        print(json.dumps(out))
    sys.exit(0)


def cmd_statusline(a):
    try:
        data = json.loads(sys.stdin.read() or "{}")
    except Exception:
        data = {}
    sid = data.get("session_id") or "nosession"
    cw = data.get("context_window") or {}
    rl = data.get("rate_limits") or {}
    fh, sd = rl.get("five_hour") or {}, rl.get("seven_day") or {}
    rec = {"ts": epoch(), "session_id": sid, "context_pct": cw.get("used_percentage"),
           "context_size": cw.get("context_window_size"), "five_hour_pct": fh.get("used_percentage"),
           "five_hour_reset": fh.get("resets_at"), "seven_day_pct": sd.get("used_percentage"),
           "seven_day_reset": sd.get("resets_at"), "transcript": data.get("transcript_path")}
    try:
        ud = RUNTIME / "usage"
        f = ud / f"claude-{slugify(sid, 80)}.json"
        old = jload(f, {})
        if {k: v for k, v in old.items() if k != "ts"} != {k: v for k, v in rec.items() if k != "ts"} or epoch() - old.get("ts", 0) > 60:
            jsave(f, rec)
            jsave(ud / "claude-latest.json", rec)
    except Exception as e:
        log_err("statusline", e)
    cwd = Path((data.get("workspace") or {}).get("current_dir") or data.get("cwd") or os.getcwd())
    root = find_root(cwd)
    parts = []
    cav = HOME / ".claude" / ".caveman-active"
    if cav.exists():
        lvl = read(cav).strip().lower()
        parts.append("\033[38;5;172m[CAVEMAN" + ("" if lvl in ("", "full") else ":" + lvl.upper()) + "]\033[0m")
    parts.append(f"\033[36m{root.name}\033[0m" if root else f"\033[2m{cwd.name} (no pos)\033[0m")
    warn, force = float(cfg()["warn_pct"]), float(cfg()["force_pct"])
    for val, lab in ((rec["context_pct"], "ctx"), (rec["five_hour_pct"], "5h"), (rec["seven_day_pct"], "7d")):
        if val is None:
            continue
        col = "\033[31m" if val >= force else "\033[33m" if val >= warn else "\033[32m"
        parts.append(f"{col}{lab} {val:.0f}%\033[0m")
    mx = max([v for v in (rec["context_pct"], rec["five_hour_pct"], rec["seven_day_pct"]) if v is not None] or [0])
    if root and mx >= force:
        parts.append("\033[1;31mCHECKPOINT\033[0m")
    print(" | ".join(parts))


# --------------------------------------------------------------------------- commands

def cmd_init(a):
    root = Path(a.path or os.getcwd()).resolve()
    existing = find_root(root)
    if existing and existing != root and not a.here:
        print(f"Inside existing ProjectOS project {existing}; refreshing that one (use --here to nest).")
        root = existing
    init_project(root, agent=a.agent or detect_agent(), do_scan=not a.no_scan, graphify=not a.no_graphify, quiet=a.quiet)


KIT_ITEMS = ("engine/pos.py", "engine/PROTOCOL.md", "skills", "agents", "install.py", "VERSION", "README.md")


def package_source():
    """Full package dir (engine/ + skills/ + agents/ + install.py), from the installer's SOURCE note or layout."""
    for cand in (Path(read(ENGINE_DIR / "SOURCE").strip() or "."), ENGINE_DIR.parent):
        if cand.is_dir() and (cand / "install.py").exists() and (cand / "skills").is_dir():
            return cand.resolve()
    return None


def vendor_kit(root):
    src = package_source()
    if not src:
        raise SystemExit("pos: package source not found (reinstall with install.py)")
    kit = pdir(root) / "kit"
    if src == kit.resolve():
        return kit
    for item in KIT_ITEMS:
        s, d = src / item, kit / item
        if s.is_dir():
            if d.exists():
                shutil.rmtree(d)
            shutil.copytree(s, d, ignore=shutil.ignore_patterns("desktop.ini", "__pycache__", "*.pyc"))
        elif s.exists():
            d.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(s, d)
    return kit


def cmd_vendor(a):
    root = need_root()
    kit = vendor_kit(root)
    install_entry_files(root)
    ledger(root, "projectos-vendor", str(kit))
    print(f"kit vendored -> {kit}\nOn any machine: `python .projectos/kit/install.py` installs the full system from this repo.")


def cmd_upgrade(a):
    roots = []
    if a.all:
        roots = [Path(p) for p in jload(POS_HOME / "registry.json", {}).keys()]
    else:
        roots = [need_root()]
    for r in roots:
        if not (r / DIRNAME / "project.json").exists():
            continue
        copy_protocol(r)
        install_entry_files(r)
        if (pdir(r) / "kit").exists():
            vendor_kit(r)
        m = meta(r)
        m.update({"protocol_version": PROTOCOL_VERSION, "engine_version": VERSION})
        jsave(pdir(r) / "project.json", m)
        print(f"upgraded {r}")


def cmd_scan(a):
    root = need_root()
    r = scan(root)
    if not a.quiet:
        print(f"scan.md written: {r['files']} files, {r['code_src']} code files, ~{r['code_loc']} LOC, class {r['size_class']}"
              f"{' (TRUNCATED)' if r['truncated'] else ''}\nplan: {r['plan']}")


def cmd_brief(a):
    root = need_root()
    agent = a.agent or detect_agent()
    sid, _ = own_session(root, agent)
    print(brief(root, agent, sid, None, max_chars=a.max_chars))


def cmd_status(a):
    root = need_root()
    m = meta(root)
    st = map_staleness(root)
    cks = list_checkpoints(root)
    ents = parse_entries(read(pdir(root) / "ATTEMPTS.md"))
    decs = parse_entries(read(pdir(root) / "DECISIONS.md"))
    agent = a.agent or detect_agent()
    u = usage(agent) if agent in ("claude", "codex") else {}
    info = {
        "root": str(root), "id": m.get("id"), "protocol": m.get("protocol_version"), "map": st,
        "purpose": ("template" if "<!-- pos:template" in (pt := read(pdir(root) / "PURPOSE.md"))
                    else "draft" if "<!-- pos:draft" in pt else "defined"),
        "checkpoints": len(cks), "last_checkpoint": cks[-1] if cks else None,
        "last_checkpoint_age": ago(epoch() - ck_time(cks[-1])) if cks else None,
        "attempts": len(ents), "decisions": len(decs),
        "other_sessions": [{"sid": k, **v} for k, v in other_sessions(root, env_session_id())],
        "usage": u,
    }
    if a.json:
        print(json.dumps(info, indent=2, default=str))
        return
    print(f"project   {m.get('name')} [{m.get('id')}]  {root}")
    print(f"purpose   {info['purpose'].upper() if info['purpose'] != 'defined' else 'defined'}")
    print(f"map       {st['status']}" + (f" (changed: {', '.join(st.get('changed', []))})" if st.get("changed") else "")
          + (f" | {st.get('mapped_by')} @ {st.get('mapped_at')}" if st.get("mapped_at") else ""))
    print(f"ckpts     {len(cks)} | last {info['last_checkpoint'] or '-'} ({info['last_checkpoint_age'] or '-'})")
    print(f"log       {len(ents)} attempts, {len(decs)} decisions")
    for s in info["other_sessions"]:
        print(f"active    {s.get('agent')} {s['sid'][:8]} task={s.get('task', '-')}")
    if u:
        print(f"usage     {usage_str(u)}")


def cmd_map_done(a):
    root = need_root()
    files, _, _ = walk(root, budget_s=30.0)
    r = {"fp": fingerprints(files), "git": _git(root)}
    lock = map_lock(root) or {"fp": {}}
    if a.groups:
        for g in [x.strip() for x in a.groups.split(",") if x.strip()]:
            if g in r["fp"]:
                lock["fp"][g] = r["fp"][g]
            else:
                lock["fp"].pop(g, None)
    else:
        lock["fp"] = r["fp"]
    lock.update({"mapped_at": iso(), "mapped_by": a.model, "confidence": a.confidence,
                 "git_head": (r.get("git") or {}).get("head")})
    jsave(pdir(root) / "map.lock.json", lock)
    mp = pdir(root) / "MAP.md"
    t = read(mp)
    fm = f"---\nstatus: mapped\nmapped_at: {lock['mapped_at']}\nmapped_by: {a.model}\nconfidence: {a.confidence}\n---\n"
    t = re.sub(r"^---\n.*?\n---\n", "", t, count=1, flags=re.S)
    write(mp, fm + t.lstrip())
    register(root, {"map": "mapped"})
    ledger(root, "projectos-map", f"{a.model} conf={a.confidence} groups={a.groups or 'all'}")
    print(f"map recorded ({a.model}, {a.confidence}); {len(lock['fp'])} groups fingerprinted")


def cmd_checkpoint(a):
    root = need_root()
    agent = a.agent or detect_agent()
    sid, srec = own_session(root, agent)
    tp = srec.get("transcript") if srec else None
    last = last_checkpoint_ts(root)
    dg = transcript_digest(agent, tp, since=last, root=root) if tp else None
    if dg:
        dg["prompts"] = []
    u = usage(agent, sid, tp) if agent in ("claude", "codex") else {}
    name = write_checkpoint(root, agent, a.message, kind="auto" if a.auto else "manual", reason=a.reason,
                            next_steps=a.next, advances=a.advances, digest=dg, usage_s=usage_str(u) if u else None)
    st = sess_state(sid)
    if st:
        st["last_manual_ck"] = epoch()
        sess_save(sid, st)
    print(f"checkpoint written: .projectos/checkpoints/{name}")


def _append_entry(path, entry):
    t = read(path)
    if not t.endswith("\n"):
        t += "\n"
    write(path, t + "\n" + entry + "\n")


def cmd_attempt(a):
    root = need_root()
    path = pdir(root) / "ATTEMPTS.md"
    ents = parse_entries(read(path))
    agent = getattr(a, "agent", None) or detect_agent()
    if a.sub == "add":
        aid = next_id(ents, "A")
        status = a.status.upper()
        L = [f"### {aid} | {a.date or now().date()} | {agent} | {status} | {one_line(a.title, 120)}"]
        for k in ("goal", "approach", "result", "evidence", "lesson", "tags"):
            v = getattr(a, k, None)
            if v:
                L.append(f"- {k}: {v}")
        _append_entry(path, "\n".join(L))
        print(f"{aid} logged ({status})")
    elif a.sub == "update":
        text = read(path)
        m = re.search(r"^###\s+" + re.escape(a.id) + r"\s*\|.*$", text, re.M)
        if not m:
            sys.exit(f"pos: {a.id} not found")
        head = m.group(0)
        parts = [p for p in head.split("|")]
        if a.status:
            parts[3] = f" {a.status.upper()} "
        new_head = "|".join(parts)
        extra = []
        for k in ("result", "evidence", "lesson"):
            v = getattr(a, k, None)
            if v:
                extra.append(f"- {k}: {v}")
        if a.note:
            extra.append(f"- update {now().date()}: {a.note}")
        end = text.find("\n###", m.end())
        end = len(text) if end < 0 else end
        block = text[m.start():end].rstrip("\n").replace(head, new_head, 1)
        block += ("\n" + "\n".join(extra)) if extra else ""
        write(path, text[:m.start()] + block + "\n" + text[end:])
        print(f"{a.id} updated")
    else:
        out = ents
        if a.status:
            out = [e for e in out if e["status"].upper().startswith(a.status.upper())]
        if a.grep:
            kws = [k.lower() for k in re.split(r"[,\s]+", a.grep) if k]
            out = [e for e in out if any(k in "\n".join(e["raw"]).lower() for k in kws)]
        out = out[-a.last:] if a.last else out
        if not out:
            print("no matching attempts")
        for e in out:
            if a.full:
                print("\n".join(e["raw"]).rstrip() + "\n")
            else:
                les = e["fields"].get("lesson") or e["fields"].get("result") or ""
                print(f"{e['id']} {e['date']} {e['status']:<10} {one_line(e['title'], 70)}" + (f"  -> {one_line(les, 120)}" if les else ""))


def cmd_decide(a):
    root = need_root()
    path = pdir(root) / "DECISIONS.md"
    text = read(path)
    ents = parse_entries(text)
    agent = a.agent or detect_agent()
    did = next_id(ents, "D")
    status = "PROPOSED" if a.proposed else "ACTIVE"
    L = [f"### {did} | {a.date or now().date()} | {agent} | {status} | {one_line(a.title, 120)}"]
    for k in ("decision", "why", "alternatives", "supersedes"):
        v = getattr(a, k, None)
        if v:
            L.append(f"- {k}: {v}")
    if a.supersedes:
        for old in [x.strip() for x in a.supersedes.split(",")]:
            text = re.sub(r"^(###\s+" + re.escape(old) + r"\s*\|[^|]*\|[^|]*\|)\s*[^|]*\|",
                          lambda mm: f"{mm.group(1)} SUPERSEDED by {did} |", text, flags=re.M)
        write(path, text)
    _append_entry(path, "\n".join(L))
    print(f"{did} recorded ({status})")


def cmd_decisions(a):
    root = need_root()
    for d in parse_entries(read(pdir(root) / "DECISIONS.md")):
        if a.all or d["status"].upper().startswith(("ACTIVE", "PROPOSED")):
            if a.grep and a.grep.lower() not in "\n".join(d["raw"]).lower():
                continue
            print(f"{d['id']} {d['date']} {d['status']:<10} {d['title']}" + (f"  -> {one_line(d['fields'].get('decision', ''), 110)}" if d['fields'].get('decision') else ""))


def cmd_claim(a):
    root = need_root()
    agent = a.agent or detect_agent()
    sid, _ = own_session(root, agent)
    sid = sid or f"manual-{agent}-{MACHINE}"
    others = [(k, v) for k, v in other_sessions(root, sid) if v.get("task")]
    sessions_touch(root, sid, agent, task=a.task, force=True)
    print(f"claimed: {a.task} ({agent}, session {sid[:8]})")
    for k, v in others:
        print(f"  note: {v.get('agent')} ({k[:8]}) is working on: {v['task']}")


def cmd_release(a):
    root = need_root()
    agent = a.agent or detect_agent()
    sid, _ = own_session(root, agent)
    if sid:
        sessions_touch(root, sid, agent, task="", force=True)
    print("released")


def cmd_usage(a):
    agent = a.agent or detect_agent()
    tp = a.transcript
    sid = a.session or env_session_id()
    root = find_root(Path.cwd())
    if not tp and root and agent in ("claude", "codex"):
        _, rec = own_session(root, agent)
        tp = rec.get("transcript") if rec else None
    u = usage(agent, sid, tp)
    if a.json:
        print(json.dumps(u))
    else:
        lvl = u.get("max")
        verdict = ("unknown" if lvl is None else "CHECKPOINT NOW" if lvl >= cfg()["force_pct"]
                   else "checkpoint soon" if lvl >= cfg()["warn_pct"] else "ok")
        print(f"{agent}: {usage_str(u)} -> {verdict}")


def cmd_export(a):
    root = need_root()
    pd = pdir(root)
    m = meta(root)
    L = [f"# Handoff - {m.get('name')} [{m.get('id')}]", f"generated {iso()} from {root}",
         "Paste this into any agent. Source of truth stays in .projectos/.", "",
         "## Purpose", clean_md(read(pd / "PURPOSE.md").split("\n", 1)[-1]), "",
         "## State", clean_md(read(pd / "STATE.md").split("\n", 1)[-1]), "",
         "## Map", clean_md(re.sub(r"^---\n.*?\n---\n", "", read(pd / "MAP.md"), flags=re.S)), "",
         "## Decisions in force"]
    for d in parse_entries(read(pd / "DECISIONS.md")):
        if d["status"].upper().startswith(("ACTIVE", "PROPOSED")):
            L.append(f"- {d['id']} {d['title']}: {d['fields'].get('decision', '')} (why: {d['fields'].get('why', '')})")
    L += ["", "## Attempts (all, one line each)"]
    for e in parse_entries(read(pd / "ATTEMPTS.md")):
        L.append(f"- {e['id']} {e['status']} {e['title']} -> {e['fields'].get('result', '')} | lesson: {e['fields'].get('lesson', '')}")
    L += ["", "## Recent checkpoints"]
    for c in list_checkpoints(root)[-5:]:
        fm = ck_front(root, c)
        L.append(f"- {c} [{fm.get('kind')}] {fm.get('agent')}: {fm.get('summary')}")
    out = "\n".join(L) + "\n"
    out = re.sub(r"^#{2,4} [^\n]*\n(?=\s*(?:#{2,4} |\Z))", "", out, flags=re.M)  # drop empty template headings
    if a.out:
        write(Path(a.out), out)
        print(f"handoff written: {a.out}")
    else:
        print(out)


def cmd_list(a):
    recs = {}
    h = hub_dir()
    if h:
        for p in (h / "projects").glob("*/projectos.json"):
            r = jload(p, {})
            if r.get("id"):
                recs[r["id"]] = r
    for path, r in jload(POS_HOME / "registry.json", {}).items():
        if r.get("id") and r["id"] not in recs:
            recs[r["id"]] = {**r, "paths": {MACHINE: path}}
    if not recs:
        print("no ProjectOS projects registered yet")
        return
    for pid, r in sorted(recs.items(), key=lambda kv: kv[1].get("updated", ""), reverse=True):
        lc = r.get("last_checkpoint") or {}
        path = (r.get("paths") or {}).get(MACHINE) or next(iter((r.get("paths") or {}).values()), "?")
        print(f"{pid:<34} map={r.get('map', '?'):<9} last={lc.get('time', '-')[:16]:<16} {one_line(lc.get('summary', ''), 50)}\n"
              f"{'':<34} {path}")


def cmd_discover(a):
    base = Path(a.root or os.getcwd()).resolve()
    found = []

    def rec(d, depth):
        if depth > a.depth:
            return
        try:
            kids = [Path(e.path) for e in os.scandir(d) if e.is_dir(follow_symlinks=False)]
        except Exception:
            return
        for c in kids:
            n = c.name
            if n.startswith(".") or n in SKIP_DIRS or n.lower() in DATA_DIRS:
                continue
            if (c / DIRNAME / "project.json").exists():
                found.append((c, "mapped" if map_lock(c) else "initialized"))
            elif has_marker(c):
                found.append((c, "unmapped"))
            else:
                rec(c, depth + 1)

    if (base / DIRNAME / "project.json").exists() or has_marker(base):
        found.append((base, "mapped" if map_lock(base) else "initialized" if (base / DIRNAME).exists() else "unmapped"))
    else:
        rec(base, 1)
    rows = []
    for p, stt in found:
        row = {"path": str(p), "status": stt}
        if a.init and stt == "unmapped" and not ignored_reason(p):
            _, r = init_project(p, agent=detect_agent(), graphify=not a.no_graphify, quiet=True, graphify_wait=True)
            row.update(status="initialized", size_class=r["size_class"] if r else None, plan=r["plan"] if r else None)
        elif stt != "unmapped":
            sc = jload(p / DIRNAME / ".runtime" / "scan.json", {})
            row.update(size_class=sc.get("size_class"), plan=sc.get("plan"))
        rows.append(row)
    if a.json:
        print(json.dumps(rows, indent=2))
        return
    for r in rows:
        print(f"{r['status']:<12} {r.get('size_class') or '-':<2} {r['path']}")
    need = [r for r in rows if r["status"] == "initialized"]
    if need:
        print(f"\n{len(need)} project(s) need an LLM map -> skill pos-map (mode: all), cheap tier only.")


def cmd_doctor(a):
    ok = lambda b: "OK " if b else "-- "
    print(f"pos {VERSION} (protocol v{PROTOCOL_VERSION}) engine={ENGINE_DIR}")
    print(ok(sys.version_info >= (3, 9)) + f"python {sys.version.split()[0]}")
    print(ok((ENGINE_DIR / "PROTOCOL.md").exists()) + "PROTOCOL.md next to engine")
    h = hub_dir()
    print(ok(bool(h)) + f"graphify hub: {h or 'not found (optional)'}")
    cs = jload(HOME / ".claude" / "settings.json", {})
    hooks_txt = json.dumps(cs.get("hooks", {}))
    print(ok("pos.py" in hooks_txt) + "Claude Code hooks")
    print(ok("pos.py" in json.dumps(cs.get("statusLine", {}))) + "Claude statusline (rate-limit capture)")
    print(ok((HOME / ".claude" / "agents" / "pos-mapper.md").exists()) + "Claude subagents pos-mapper/pos-architect")
    print(ok("pos.py" in read(HOME / ".codex" / "hooks.json")) + "Codex hooks.json (trust them once via /hooks in Codex)")
    for nm, p in (("Claude", HOME / ".claude" / "skills" / "projectos"), ("Codex", HOME / ".codex" / "skills" / "projectos"),
                  ("Antigravity", HOME / ".gemini" / "skills" / "projectos")):
        print(ok(p.exists()) + f"{nm} skills")
    for nm, p in (("~/.claude/CLAUDE.md", HOME / ".claude" / "CLAUDE.md"), ("~/.codex/AGENTS.md", HOME / ".codex" / "AGENTS.md"),
                  ("~/.gemini/GEMINI.md", HOME / ".gemini" / "GEMINI.md")):
        print(ok("projectos:begin" in read(p)) + f"global block in {nm}")
    print(ok(shutil.which("pos") is not None) + "`pos` on PATH")
    err = RUNTIME / "errors.log"
    if err.exists():
        tail = read(err).strip().splitlines()[-3:]
        print("recent hook errors:\n  " + "\n  ".join(tail))


def main(argv=None):
    ap = argparse.ArgumentParser(prog="pos", description="ProjectOS - shared project memory for all agents")
    ap.add_argument("--version", action="version", version=f"pos {VERSION}")
    sp = ap.add_subparsers(dest="cmd", required=True)

    p = sp.add_parser("init", help="create/refresh .projectos here (tier-0 scan, no LLM)")
    p.add_argument("path", nargs="?")
    p.add_argument("--agent")
    p.add_argument("--no-scan", action="store_true")
    p.add_argument("--no-graphify", action="store_true")
    p.add_argument("--here", action="store_true", help="allow nesting inside another project")
    p.add_argument("--quiet", action="store_true")
    p.set_defaults(fn=cmd_init)

    p = sp.add_parser("upgrade", help="refresh PROTOCOL.md + entry blocks (current project or --all)")
    p.add_argument("--all", action="store_true")
    p.set_defaults(fn=cmd_upgrade)

    p = sp.add_parser("scan", help="tier-0 mechanical scan -> .projectos/scan.md")
    p.add_argument("--quiet", action="store_true")
    p.set_defaults(fn=cmd_scan)

    p = sp.add_parser("brief", help="compact context digest for starting/resuming work")
    p.add_argument("--agent")
    p.add_argument("--max-chars", type=int)
    p.set_defaults(fn=cmd_brief)

    p = sp.add_parser("status", help="project health: map freshness, checkpoints, sessions, usage")
    p.add_argument("--agent")
    p.add_argument("--json", action="store_true")
    p.set_defaults(fn=cmd_status)

    p = sp.add_parser("map-done", help="record that MAP.md was (re)written; fingerprints groups")
    p.add_argument("--model", required=True, help="model that produced the map, e.g. haiku, sonnet, gpt-5-mini")
    p.add_argument("--confidence", default="medium", choices=["high", "medium", "low"])
    p.add_argument("--groups", help="comma list of groups remapped (incremental); default all")
    p.set_defaults(fn=cmd_map_done)

    p = sp.add_parser("checkpoint", help="snapshot STATE.md + mechanical delta into checkpoints/")
    p.add_argument("-m", "--message", required=True, help="one-line summary of what was achieved")
    p.add_argument("--next", help="next steps (if not already in STATE.md)")
    p.add_argument("--advances", help="which PURPOSE success criterion this advances")
    p.add_argument("--reason", default="milestone", help="milestone|usage|periodic|handoff|end")
    p.add_argument("--agent")
    p.add_argument("--auto", action="store_true")
    p.set_defaults(fn=cmd_checkpoint)

    p = sp.add_parser("attempt", help="log/list/update attempts (ATTEMPTS.md)")
    asp = p.add_subparsers(dest="sub", required=True)
    q = asp.add_parser("add")
    q.add_argument("-t", "--title", required=True)
    q.add_argument("-s", "--status", default="RUNNING", type=str.upper, choices=STATUSES)
    for k in ("goal", "approach", "result", "evidence", "lesson", "tags", "agent", "date"):
        q.add_argument(f"--{k}")
    q.set_defaults(fn=cmd_attempt)
    q = asp.add_parser("update")
    q.add_argument("id")
    q.add_argument("-s", "--status", type=str.upper, choices=STATUSES)
    for k in ("result", "evidence", "lesson", "note"):
        q.add_argument(f"--{k}")
    q.set_defaults(fn=cmd_attempt)
    q = asp.add_parser("list")
    q.add_argument("--status")
    q.add_argument("--grep")
    q.add_argument("--last", type=int)
    q.add_argument("--full", action="store_true")
    q.set_defaults(fn=cmd_attempt)

    p = sp.add_parser("decide", help="record a decision (DECISIONS.md)")
    p.add_argument("-t", "--title", required=True)
    for k in ("decision", "why", "alternatives", "supersedes", "agent", "date"):
        p.add_argument(f"--{k}")
    p.add_argument("--proposed", action="store_true", help="needs user confirmation")
    p.set_defaults(fn=cmd_decide)

    p = sp.add_parser("decisions", help="list decisions in force")
    p.add_argument("--all", action="store_true")
    p.add_argument("--grep")
    p.set_defaults(fn=cmd_decisions)

    p = sp.add_parser("claim", help="announce what this agent works on (overlap avoidance)")
    p.add_argument("task")
    p.add_argument("--agent")
    p.set_defaults(fn=cmd_claim)
    p = sp.add_parser("release", help="clear this agent's claim")
    p.add_argument("--agent")
    p.set_defaults(fn=cmd_release)

    p = sp.add_parser("usage", help="context / rate-limit usage for claude or codex")
    p.add_argument("--agent")
    p.add_argument("--session")
    p.add_argument("--transcript")
    p.add_argument("--json", action="store_true")
    p.set_defaults(fn=cmd_usage)

    p = sp.add_parser("export", help="single-file handoff (purpose+state+map+decisions+attempts)")
    p.add_argument("--out")
    p.set_defaults(fn=cmd_export)

    p = sp.add_parser("list", help="all ProjectOS projects (hub + local registry)")
    p.set_defaults(fn=cmd_list)

    p = sp.add_parser("discover", help="find projects under a folder; --init scaffolds unmapped ones (no LLM)")
    p.add_argument("root", nargs="?")
    p.add_argument("--depth", type=int, default=3)
    p.add_argument("--init", action="store_true")
    p.add_argument("--no-graphify", action="store_true")
    p.add_argument("--json", action="store_true")
    p.set_defaults(fn=cmd_discover)

    p = sp.add_parser("vendor", help="copy the full kit into .projectos/kit/ (repo becomes self-installing)")
    p.set_defaults(fn=cmd_vendor)

    p = sp.add_parser("doctor", help="verify installation")
    p.set_defaults(fn=cmd_doctor)

    p = sp.add_parser("hook", help="(internal) agent hook entry point")
    p.add_argument("event", choices=sorted(HOOKS))
    p.add_argument("--agent", default="claude")
    p.set_defaults(fn=cmd_hook)

    p = sp.add_parser("statusline", help="(internal) Claude Code statusline; captures rate limits")
    p.set_defaults(fn=cmd_statusline)

    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
