"""Optional R bridge: find an Rscript, ask it once which packages it has, run R scripts as subprocesses.

grnbench never installs anything and never loads R into Python (no rpy2, D-032): R runs as
``Rscript --vanilla`` in a child process, by default from the dedicated conda env ``grnbench-r``
(``requirements/r-bridge.yml`` + ``requirements/r-bridge.R``, see ``docs/R_BRIDGE.md``).

* :func:`find_rscript` picks the R: ``--rscript`` > ``$GRNBENCH_RSCRIPT`` (``none`` turns R off) > the
  conda env ``$GRNBENCH_R_ENV`` (default ``grnbench-r``) > ``$R_HOME`` > ``Rscript`` on PATH. An R that
  was named, or a conda env that exists, is used or is an error: never a quiet switch to another R. A
  conda R carries the folders it needs in front of PATH (on Windows its DLLs are not found otherwise).
* :func:`probe` asks that R for its version, library folders and package versions (and the R each
  package was built for): one subprocess per R and library state, remembered in this process and on
  disk, so spawned workers do not start R again; ``fresh=True`` asks again (diagnostics).
* :func:`run_r` runs a method script through ``r/_wrapper.R``: seed and RNG kind fixed by the wrapper,
  inputs and parameters in ``params.json``, everything R prints in ``r.log``, R and package versions in
  ``session.json``; a timeout or an interrupt kills the whole process tree. :func:`run_external` does the
  same for a command-line tool.
* :func:`requirement_status` and :func:`environment_report` (``grnbench env``) say what is there and how
  to get the rest.

Not wired into the knowledge base yet: ``kb.Option.availability`` is unchanged (step 2c).
"""
from __future__ import annotations

import copy
import hashlib
import importlib.metadata
import importlib.util
import json
import numbers
import os
import re
import shutil
import signal
import subprocess
import tempfile
import time
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import TYPE_CHECKING, Any, Mapping, Sequence

if TYPE_CHECKING:
    from types import EllipsisType

R_DIR = Path(__file__).with_name("r")
WRAPPER = R_DIR / "_wrapper.R"
PROBE = R_DIR / "_probe.R"
PROBE_TIMEOUT_S = 30            # a cold R start behind an antivirus scan; the answer is cached, so this is rare
DEFAULT_R_ENV = "grnbench-r"
DOC = "docs/R_BRIDGE.md"
R_REQUIREMENTS: tuple[str, ...]     # computed from the knowledge base on first use (module __getattr__)

_WINDOWS = os.name == "nt"
_EXE = ".exe" if _WINDOWS else ""
_TAIL = 20
_LIB_VARS = ("R_LIBS_USER", "R_LIBS", "R_LIBS_SITE")
_PIP_NAMES = {"scvi": "scvi-tools", "magic": "magic-impute"}   # import name -> pip name, for install hints
_BIOC = frozenset({"DropletUtils", "scDblFinder", "celda", "miQC", "scran", "slingshot"})
_CRAN = frozenset({"SoupX", "sctransform", "jsonlite", "Matrix"})
_GITHUB = {"DoubletFinder": "chris-mcginnis-ucsf/DoubletFinder"}
_PROBE_MEMO: dict[str, dict] = {}
_CONDA_MEMO: dict[str, list[str]] = {}


def _home() -> Path | None:
    try:
        return Path.home()
    except (RuntimeError, KeyError, OSError):     # no HOME / USERPROFILE (some containers)
        return None


def _scrub(text: str) -> str:
    """``text`` with the home folder written ``~`` wherever it appears (messages people read and paste)."""
    home = _home()
    h = str(home).rstrip("\\/") if home else ""
    if len(h) < 3:
        return text
    back = h.replace("/", "\\")
    # also the repr() spelling, which doubles backslashes (OSError texts and {path!r} in messages)
    for variant in sorted({h, h.replace("\\", "/"), back, back.replace("\\", "\\\\")}, key=len, reverse=True):
        text = re.sub(re.escape(variant) + r"(?![\w.-])", "~", text, flags=re.IGNORECASE if _WINDOWS else 0)
    return text


class RBridgeError(RuntimeError):
    """R (or an external tool) could not be found, started, or did not finish its job; the message says why
    (the home folder written ``~``)."""

    def __init__(self, message: object = "") -> None:
        super().__init__(_scrub(str(message)))


@dataclass(frozen=True)
class RRuntime:
    """One usable R: the Rscript to start, its version, where it was found, its package library, the folders a
    conda R needs in front of PATH, and the conda env it belongs to (None: not a conda R)."""
    rscript: str
    version: str
    source: str                 # "--rscript" | "GRNBENCH_RSCRIPT" | "conda:<env>" | "R_HOME" | "PATH"
    library: str | None
    path_prepend: tuple[str, ...] = ()
    conda_prefix: str | None = None

    def subprocess_env(self, base: Mapping[str, str] | None = None) -> dict[str, str]:
        """The environment of every R child: ``base`` (default ``os.environ``) with ``R_HOME`` removed unless
        this R was found through it. A conda R also gets PATH prefixed by ``path_prepend``, ``R_LIBS`` and
        ``R_LIBS_SITE`` removed and ``R_LIBS_USER`` = its env library, so no other R's packages leak in; any
        other R keeps its own library settings (the layout of its installation is not known)."""
        env = dict(os.environ if base is None else base)
        if self.source != "R_HOME":
            _pop(env, "R_HOME")
        if self.path_prepend:
            old = _pop(env, "PATH")
            env["PATH"] = os.pathsep.join([*self.path_prepend, *([old] if old else [])])
        if self.conda_prefix:
            for name in _LIB_VARS:
                _pop(env, name)
            if self.library:
                env["R_LIBS_USER"] = self.library
        return env


@dataclass(frozen=True)
class ReqStatus:
    ok: bool
    reason: str                 # "" when ok; "missing: <req>" (the knowledge base's wording) or "built for R x.y"
    install_hint: str           # how to get it ("" when ok)


@dataclass(frozen=True)
class RResult:
    outputs: dict[str, Path]
    log: Path
    session: dict
    seconds: float
    runtime: RRuntime
    skipped: tuple[str, ...] = ()       # R candidates passed over before this one (why each was skipped)

    def manifest_block(self) -> dict:
        """What a run manifest records about R; the home folder is written ``~`` in every path.

        ``rng_changed`` is true when the method changed the generator kind the wrapper set or called
        ``set.seed`` itself (each call is in ``set_seed_calls``; a call written ``base::set.seed`` or made inside
        a package is not seen)."""
        s = self.session
        kind, after = s.get("rng_kind"), s.get("rng_kind_after")
        calls = list(s.get("set_seed_calls") or [])
        return {"version": self.runtime.version, "source": self.runtime.source,
                "rscript": tilde(self.runtime.rscript), "library": tilde(self.runtime.library),
                "packages": dict(s.get("packages") or {}), "seed": s.get("seed"), "rng_kind": kind,
                "rng_kind_after": after, "set_seed_calls": calls,
                "rng_changed": bool(calls) or (after is not None and after != kind),
                "seconds": round(self.seconds, 3), "log": tilde(self.log),
                "skipped": [_scrub(n) for n in self.skipped]}


@dataclass(frozen=True)
class ExternalResult:
    returncode: int
    log: Path
    seconds: float
    version: str | None


# -- requirements from the knowledge base ------------------------------------------------- #
@lru_cache(maxsize=1)
def _r_requirements() -> tuple[str, ...]:
    """The ``r:`` packages the knowledge base names, plus jsonlite (the wrapper reads params.json with it)."""
    names = {req.split(":", 1)[1] for req in _kb_uses() if req.startswith("r:")}
    return tuple(sorted(names | {"jsonlite"}))


@lru_cache(maxsize=1)
def _kb_uses() -> dict[str, tuple[str, ...]]:
    """requirement -> the ``decision=option`` pairs of the knowledge base that need it."""
    from grnbench import kb
    uses: dict[str, list[str]] = {}
    for d in kb.DECISIONS:
        for o in d.options:
            for req in o.requires:
                uses.setdefault(req, []).append(f"{d.id}={o.id}")
    return {req: tuple(pairs) for req, pairs in uses.items()}


def __getattr__(name: str):
    if name == "R_REQUIREMENTS":
        return _r_requirements()
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


# -- small helpers ------------------------------------------------------------------------ #
def tilde(path: str | os.PathLike | None) -> str | None:
    """``path`` with the home folder written ``~`` (and forward slashes), for manifests and printed reports."""
    if path is None:
        return None
    text = str(path)
    home = _home()
    if home is None:
        return text
    norm, h = text.replace("\\", "/"), str(home).replace("\\", "/").rstrip("/")
    a, b = (norm.lower(), h.lower()) if _WINDOWS else (norm, h)
    if h and a.startswith(b) and (len(a) == len(b) or a[len(b)] == "/"):
        return "~" + norm[len(h):]
    return text


def clear_memo() -> None:
    """Forget what this process remembers about R and conda (the disk cache stays): a new worker's view.
    It also forgets the R the knowledge base looked up (``kb._r_runtime``), failed lookups included."""
    _PROBE_MEMO.clear()
    _CONDA_MEMO.clear()
    from grnbench import kb
    kb._R_RUNTIME.clear()


def cache_dir() -> Path:
    """Where probe answers are kept: ``$GRNBENCH_CACHE_DIR``, else the per-user cache folder of the system
    (a folder in the temporary directory when there is no home folder)."""
    chosen = os.environ.get("GRNBENCH_CACHE_DIR", "").strip()
    if chosen:
        return Path(chosen).expanduser()
    home = _home()
    if _WINDOWS:
        local = os.environ.get("LOCALAPPDATA", "").strip()
        if local or home:
            return (Path(local) if local else home / "AppData" / "Local") / "grnbench" / "cache"
    else:
        xdg = os.environ.get("XDG_CACHE_HOME", "").strip()
        if xdg or home:
            return (Path(xdg) if xdg else home / ".cache") / "grnbench"
    return Path(tempfile.gettempdir()) / "grnbench-cache"


def _pop(env: dict[str, str], name: str) -> str | None:
    """Remove ``name`` from ``env`` in any spelling on Windows (names are case-insensitive there)."""
    keys = [k for k in env if (k.upper() == name if _WINDOWS else k == name)]
    values = [env.pop(k) for k in keys]
    return values[0] if values else None


def _norm(path: str | os.PathLike) -> str:
    return os.path.normcase(os.path.normpath(str(path)))


def _mtime_ns(path: str | os.PathLike | None) -> int | None:
    if not path:
        return None
    try:
        return os.stat(path).st_mtime_ns
    except OSError:
        return None


def _key(kind: str, *parts: Any) -> str:
    return f"{kind}-" + hashlib.sha256(json.dumps(parts, default=str).encode("utf-8")).hexdigest()[:40]


def _cache_read(key: str) -> dict | None:
    try:
        entry = json.loads((cache_dir() / f"{key}.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return entry if isinstance(entry, dict) else None


def _cache_write(key: str, entry: dict) -> None:
    tmp = None
    try:
        folder = cache_dir()
        folder.mkdir(parents=True, exist_ok=True)
        tmp = folder / f"{key}.{os.getpid()}.tmp"
        tmp.write_text(json.dumps(entry), encoding="utf-8")
        os.replace(tmp, folder / f"{key}.json")
    except OSError:     # a convenience only: without it R is asked again next time
        if tmp is not None:
            try:
                tmp.unlink()
            except OSError:
                pass


def _lines(text: str, n: int = _TAIL) -> str:
    return "\n".join(text.splitlines()[-n:])


def _log_tail(log: Path, n: int = _TAIL, max_bytes: int = 65536) -> str:
    """The last ``n`` lines of ``log``, reading only its end."""
    try:
        with open(log, "rb") as fh:
            fh.seek(0, os.SEEK_END)
            fh.seek(max(0, fh.tell() - max_bytes))
            return _lines(fh.read().decode("utf-8", errors="replace"), n)
    except OSError:
        return ""


def _exit_hint(code: int) -> str:
    if code & 0xFFFFFFFF == 0xC0000135:
        return (" (0xC0000135, a DLL was not found: a conda R on Windows needs its env folders at the front of "
                f"PATH; point at the env's Scripts\\Rscript.exe or let grnbench find the env, see {DOC})")
    return ""


def _exit_text(code: int) -> str:
    if code < 0 and not _WINDOWS:
        return f"was killed by signal {-code}"
    return f"exited with code {code}{_exit_hint(code)}"


def _skipped_lines(notes: Sequence[str]) -> str:
    return "".join(f"\n# skipped: {n}" for n in notes)


# -- subprocesses ------------------------------------------------------------------------- #
@dataclass(frozen=True)
class _Done:
    code: int | None            # None: the timeout killed it
    out: str = ""
    err: str = ""
    orphan: str = ""            # "; process may still be running (pid N)" when killing it was not confirmed


def _orphan(pid: int) -> str:
    return f"; process may still be running (pid {pid})"


def _popen_kw(env: Mapping[str, str] | None, cwd: Path | None) -> dict:
    kw: dict[str, Any] = {"env": None if env is None else dict(env), "cwd": cwd, "stdin": subprocess.DEVNULL}
    if not _WINDOWS:
        kw["start_new_session"] = True      # its own process group: a timeout kills R and what it started
    return kw


def _kill_tree(proc: subprocess.Popen) -> bool:
    """Kill ``proc`` and every process it started (Rscript starts R; a .cmd launcher starts its program).
    True when the kill is confirmed (taskkill succeeded / the process group is gone)."""
    if _WINDOWS:
        taskkill = Path(os.environ.get("SystemRoot", r"C:\Windows")) / "System32" / "taskkill.exe"
        try:
            done = subprocess.run([str(taskkill) if taskkill.is_file() else "taskkill", "/T", "/F", "/PID",
                                   str(proc.pid)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=30)
            confirmed = done.returncode == 0
        except (OSError, subprocess.SubprocessError):
            confirmed = False
    else:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
            confirmed = True
        except ProcessLookupError:
            confirmed = True                # the whole group is already gone
        except OSError:
            confirmed = False
    try:
        proc.kill()
    except OSError:
        pass
    try:
        proc.wait(timeout=30)
    except subprocess.TimeoutExpired:
        confirmed = False
    return confirmed


def _run_captured(cmd: Sequence[str], *, timeout_s: float, env: Mapping[str, str] | None = None,
                  cwd: Path | None = None) -> _Done:
    """Run and capture stdout and stderr; on a timeout the tree is killed and the partial output kept."""
    proc = subprocess.Popen(list(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                            encoding="utf-8", errors="replace", **_popen_kw(env, cwd))
    try:
        out, err = proc.communicate(timeout=timeout_s)
    except subprocess.TimeoutExpired:
        confirmed = _kill_tree(proc)
        try:
            out, err = proc.communicate(timeout=30)
        except (subprocess.TimeoutExpired, ValueError, OSError):
            out, err = "", ""
        return _Done(None, out or "", err or "", "" if confirmed else _orphan(proc.pid))
    except BaseException:           # Ctrl+C or any error while waiting: R must not outlive us
        _kill_tree(proc)
        raise
    return _Done(proc.returncode, out, err)


def _run_logged(cmd: Sequence[str], *, log: Path, timeout_s: float, env: Mapping[str, str] | None = None,
                cwd: Path | None = None) -> _Done:
    """Run with stdout and stderr in ``log``; on a timeout or an interrupt the tree is killed."""
    try:
        fh = open(log, "wb")
    except OSError as exc:
        raise RBridgeError(f"cannot write the log {log}: {exc}") from None
    with fh:
        proc = subprocess.Popen(list(cmd), stdout=fh, stderr=subprocess.STDOUT, **_popen_kw(env, cwd))
        try:
            return _Done(proc.wait(timeout=timeout_s))
        except subprocess.TimeoutExpired:
            return _Done(None, orphan="" if _kill_tree(proc) else _orphan(proc.pid))
        except BaseException:       # Ctrl+C or any error while waiting: R must not outlive us
            _kill_tree(proc)
            raise


# -- finding R ---------------------------------------------------------------------------- #
def _conda_layout(prefix: Path) -> tuple[Path, Path, tuple[str, ...]]:
    """``(rscript, library, path_prepend)`` of the R inside a conda env."""
    if _WINDOWS:
        lib = prefix / "Library"
        prepend = (prefix, lib / "mingw-w64" / "bin", lib / "usr" / "bin", lib / "bin", prefix / "Scripts",
                   prefix / "bin")
        rscript = prefix / "Scripts" / f"Rscript{_EXE}"
    else:
        prepend = (prefix / "bin",)
        rscript = prefix / "bin" / "Rscript"
    return rscript, prefix / "lib" / "R" / "library", tuple(str(p) for p in prepend)


def _conda_prefix_of(rscript: Path) -> Path | None:
    """The conda env an Rscript belongs to (a parent folder holding ``conda-meta``), if any."""
    for parent in list(rscript.parents)[:5]:
        if (parent / "conda-meta").is_dir():
            return parent
    return None


def _conda_exe() -> str | None:
    chosen = os.environ.get("CONDA_EXE", "").strip()
    if chosen and Path(chosen).is_file():
        return chosen
    for name in ("conda", "mamba", "micromamba"):
        found = shutil.which(name)
        if found:
            return found
    return None


def _conda_env_list(looked: list[str], *, fresh: bool = False) -> list[str] | None:
    """Env folders from ``conda env list --json``, asked at most once per cache (conda and its env list
    unchanged). None, with the reason in ``looked``, when there is no conda or it failed (a failure is never
    cached)."""
    exe = _conda_exe()
    if exe is None:
        looked.append("no conda, mamba or micromamba found (CONDA_EXE, PATH)")
        return None
    home = _home()
    envs_txt = home / ".conda" / "environments.txt" if home else None
    key = _key("conda", _norm(exe), _mtime_ns(exe), _mtime_ns(envs_txt))
    if not fresh:
        if key in _CONDA_MEMO:
            return _CONDA_MEMO[key]
        entry = _cache_read(key)
        if entry is not None and isinstance(entry.get("envs"), list):
            _CONDA_MEMO[key] = [str(e) for e in entry["envs"]]
            return _CONDA_MEMO[key]
    asked = f"`{tilde(exe)} env list --json`"
    try:
        done = _run_captured([exe, "env", "list", "--json"], timeout_s=60)
    except OSError as exc:
        looked.append(_scrub(f"{asked} could not start: {exc}"))
        return None
    if done.code is None:
        looked.append(f"{asked} timed out after 60 s{done.orphan}")
        return None
    if done.code != 0:
        looked.append(_scrub(f"{asked} {_exit_text(done.code)}: {_lines(done.err or done.out, 3)}"))
        return None
    try:
        envs = [str(e) for e in json.loads(done.out[done.out.index("{"):])["envs"]]
    except (ValueError, KeyError, TypeError) as exc:
        looked.append(f"{asked} gave no env list ({type(exc).__name__}: {exc})")
        return None
    _CONDA_MEMO[key] = envs
    _cache_write(key, {"envs": envs})
    return envs


def _conda_prefix(name: str, looked: list[str], *, fresh: bool = False) -> Path | None:
    """The folder of the conda env ``name``: from ``~/.conda/environments.txt`` (no subprocess), else from
    ``conda env list --json``. When it is not found, ``looked`` says where it was looked for."""
    def match(prefixes: Sequence[str]) -> Path | None:
        for p in (line.strip() for line in prefixes):
            last = Path(p).name if p else ""
            if last and (last.lower() == name.lower() if _WINDOWS else last == name) and Path(p).is_dir():
                return Path(p)
        return None

    home = _home()
    listed: list[str] = []
    if home is None:
        looked.append("no home folder, so ~/.conda/environments.txt was not read")
    else:
        try:
            listed = (home / ".conda" / "environments.txt").read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError:
            looked.append("no ~/.conda/environments.txt")
        else:
            hit = match(listed)
            if hit is not None:
                return hit
            looked.append("not in ~/.conda/environments.txt")
    envs = _conda_env_list(looked, fresh=fresh)
    hit = match(envs) if envs is not None else None
    if hit is None and envs is not None:
        looked.append("not in `conda env list --json`")
    return hit


def _usable(rscript: str, source: str, notes: list[str], *, strict: bool = False, prefix: Path | None = None,
            fresh: bool = False) -> RRuntime | None:
    """The runtime behind ``rscript`` once its probe answered; otherwise None and the reason in ``notes``
    (``strict``: raise instead, for an R that was named or a conda env that exists)."""
    def skip(why: str) -> None:
        why = _scrub(why)
        if strict:
            raise RBridgeError(why)
        notes.append(why)

    path = Path(rscript).expanduser()
    if not path.is_file():
        found = shutil.which(rscript)           # a bare command name such as "Rscript"
        if found:
            path = Path(found)
    path = Path(os.path.abspath(path))          # absolute, but symlinks kept (a conda env may be one)
    if not path.is_file():
        return skip(f"{source}: no Rscript at {tilde(path)}")
    prefix = prefix or _conda_prefix_of(path)
    library: str | None = None
    prepend: tuple[str, ...] = ()
    if prefix is not None:
        _, lib, prepend = _conda_layout(prefix)
        library = str(lib) if lib.is_dir() else None
    conda = str(prefix) if prefix is not None else None
    try:
        info = probe(RRuntime(str(path), "", source, library, prepend, conda), _r_requirements(), fresh=fresh)
    except RBridgeError as exc:
        return skip(f"{source}: {tilde(path)} does not run: {exc}")
    library = library or (info["library"][0] if info["library"] else None)
    return RRuntime(str(path), info["version"], source, library, prepend, conda)


def _find(explicit: str | None = None, *, fresh: bool = False) -> tuple[RRuntime | None, list[str]]:
    """:func:`find_rscript` plus why each passed-over candidate was skipped (for ``grnbench env`` and
    manifests)."""
    notes: list[str] = []
    if explicit is not None:
        if not str(explicit).strip():
            raise RBridgeError("--rscript is empty: give the path of an Rscript")
        return _usable(str(explicit), "--rscript", notes, strict=True, fresh=fresh), notes
    chosen = os.environ.get("GRNBENCH_RSCRIPT", "").strip()
    if chosen.lower() == "none":
        return None, ["R is turned off: GRNBENCH_RSCRIPT=none"]
    if chosen:
        return _usable(chosen, "GRNBENCH_RSCRIPT", notes, strict=True, fresh=fresh), notes
    named = os.environ.get("GRNBENCH_R_ENV", "").strip()
    name = named or DEFAULT_R_ENV
    looked: list[str] = []
    prefix = _conda_prefix(name, looked, fresh=fresh)
    if prefix is not None:      # the env exists: it is the R to use, or an error (never a quiet switch)
        rscript = str(_conda_layout(prefix)[0])
        return _usable(rscript, f"conda:{name}", notes, strict=True, prefix=prefix, fresh=fresh), notes
    if named:
        raise RBridgeError(f"GRNBENCH_R_ENV={name}: no conda env called {name} ({'; '.join(looked)}); "
                           f"create it ({DOC}) or unset GRNBENCH_R_ENV")
    notes.append(_scrub(f"conda:{name}: no conda env called {name} ({'; '.join(looked)})"))
    r_home = os.environ.get("R_HOME", "").strip()
    if r_home:
        rt = _usable(str(Path(r_home) / "bin" / f"Rscript{_EXE}"), "R_HOME", notes, fresh=fresh)
        if rt is not None:
            return rt, notes
    on_path = shutil.which("Rscript")
    if on_path is None:
        notes.append("PATH: no Rscript on PATH")
        return None, notes
    return _usable(on_path, "PATH", notes, fresh=fresh), notes


def find_rscript(explicit: str | None = None, *, fresh: bool = False) -> RRuntime | None:
    """The R to use, or None when there is none.

    Order: ``explicit`` (``grnbench env --rscript``) > ``$GRNBENCH_RSCRIPT`` (``none``, any case, turns R
    off without starting anything; empty means unset) > the conda env ``$GRNBENCH_R_ENV`` (default
    ``grnbench-r``; its folder from ``~/.conda/environments.txt``, else ``conda env list --json``) >
    ``$R_HOME/bin/Rscript`` > ``Rscript`` on PATH.

    Raises :class:`RBridgeError` rather than switching quietly to another R when: an R named by ``explicit``
    or ``GRNBENCH_RSCRIPT`` does not work; ``GRNBENCH_R_ENV`` names an env that does not exist; or the conda
    env exists (named or the default) but its Rscript is missing or does not answer. Only a missing default
    env, ``R_HOME`` and PATH candidates are skipped, with the reason kept. ``fresh=True`` ignores every cache
    (for diagnostics)."""
    return _find(explicit, fresh=fresh)[0]


# -- probing ------------------------------------------------------------------------------ #
def _valid_probe(result: Any) -> bool:
    return (isinstance(result, dict) and isinstance(result.get("version"), str)
            and isinstance(result.get("library"), list) and isinstance(result.get("packages"), dict)
            and isinstance(result.get("built", {}), dict))


def _probe_now(rt: RRuntime, packages: list[str]) -> dict:
    try:
        done = _run_captured([rt.rscript, "--vanilla", str(PROBE), *packages], timeout_s=PROBE_TIMEOUT_S,
                             env=rt.subprocess_env())
    except OSError as exc:
        raise RBridgeError(f"cannot start {tilde(rt.rscript)}: {exc}") from None
    if done.code is None:
        partial = _lines(f"{done.out}\n{done.err}", 5).strip()
        raise RBridgeError(f"R probe timed out after {PROBE_TIMEOUT_S} s ({tilde(rt.rscript)})"
                           + (f"; last output: {partial}" if partial else "") + done.orphan)
    if done.code != 0:
        tail = _lines(done.err or done.out, 5)
        raise RBridgeError(_exit_text(done.code) + (f": {tail}" if tail.strip() else ""))
    for line in reversed(done.out.splitlines()):
        if line.startswith("{"):
            try:
                result = json.loads(line)
            except ValueError:
                break
            if _valid_probe(result):
                return result
            break
    raise RBridgeError(f"unexpected answer to {PROBE.name}: {_lines(done.out + done.err, 5)}")


def probe(rt: RRuntime, packages: Sequence[str] = (), *, fresh: bool = False) -> dict:
    """R's ``version``, ``version_string``, ``library`` folders, ``packages`` (name -> installed version or
    None) and ``built`` (name -> the R version the installed package was built for) from one ``Rscript`` call.

    The answer is remembered in this process and in :func:`cache_dir`, keyed by the Rscript's path and mtime,
    the package names and (for an R outside conda) the caller's ``R_LIBS_USER``/``R_LIBS``/``R_LIBS_SITE``. It
    is reused while every folder of ``library`` keeps the modification time it had before R was asked:
    installing a package changes its library folder, so the next call asks R again. A cache that cannot be
    read or written is ignored; ``fresh=True`` ignores it on purpose (diagnostics). Raises
    :class:`RBridgeError` when R does not start, fails, or takes over ``PROBE_TIMEOUT_S``."""
    pkgs = sorted({str(p) for p in packages})
    stamp = _mtime_ns(rt.rscript)
    if stamp is None:
        raise RBridgeError(f"no Rscript at {tilde(rt.rscript)}")
    lib_env = None if rt.conda_prefix else [os.environ.get(name) for name in _LIB_VARS]
    key = _key("rprobe", _norm(rt.rscript), stamp, pkgs, lib_env)
    previous = None if fresh else (_PROBE_MEMO.get(key) or _cache_read(key))
    if previous is not None and _still_valid(previous, rt):
        _PROBE_MEMO[key] = previous
        return copy.deepcopy(previous["result"])
    known = [rt.library] if rt.library else []
    if previous is not None and isinstance(previous.get("library_mtimes"), dict):
        known += [str(p) for p in previous["library_mtimes"]]
    started_ns = time.time_ns()
    before = {_norm(p): _mtime_ns(p) for p in known}       # stamps taken BEFORE R is asked
    result = _probe_now(rt, pkgs)
    mtimes: dict[str, int | None] = {}
    for lib in (str(p) for p in result["library"]):
        if _norm(lib) in before:
            mtimes[lib] = before[_norm(lib)]
        else:
            stamp_now = _mtime_ns(lib)
            # changed after R was asked: the answer may predate it, so ask again next time
            mtimes[lib] = -1 if stamp_now is not None and stamp_now > started_ns else stamp_now
    entry = {"rscript": str(rt.rscript), "packages": pkgs, "library_mtimes": mtimes, "result": result}
    _PROBE_MEMO[key] = entry
    _cache_write(key, entry)
    return copy.deepcopy(result)


def _still_valid(entry: dict, rt: RRuntime) -> bool:
    mtimes = entry.get("library_mtimes")
    if not _valid_probe(entry.get("result")) or not isinstance(mtimes, dict):
        return False
    return all(_mtime_ns(lib) == stamp for lib, stamp in mtimes.items())


def _minor(version: str) -> str:
    return ".".join(str(version).split(".")[:2])


def _package_state(info: dict, pkg: str) -> tuple[str | None, str | None]:
    """``(installed version, problem)``; problem = "built for R x.y" when built for another R than this one."""
    version = (info.get("packages") or {}).get(pkg)
    if not version:
        return None, None
    built = (info.get("built") or {}).get(pkg)
    if built and _minor(built) != _minor(info["version"]):
        return version, f"built for R {_minor(built)}"
    return version, None


# -- requirement status ------------------------------------------------------------------- #
def _packages_for(name: str) -> tuple[str, ...]:
    """The probe's package list (shared with find_rscript, so one cached answer serves both)."""
    reqs = _r_requirements()
    return reqs if name in reqs else tuple(sorted({*reqs, name}))


def _no_r_hint(notes: Sequence[str] = ()) -> str:
    if any("GRNBENCH_RSCRIPT=none" in n for n in notes):
        return "R is turned off (GRNBENCH_RSCRIPT=none): unset it to let grnbench look for R"
    create = ("conda create -n grnbench-r -c conda-forge r-base=4.5.3 r-biocmanager r-jsonlite" if _WINDOWS
              else "conda env create -f requirements/r-bridge.yml")
    return "\n".join([f"# no R found; from a clone of the repository (details: {DOC}):", create,
                      "conda run -n grnbench-r Rscript --vanilla requirements/r-bridge.R"]) + _skipped_lines(notes)


def _r_install_hint(name: str, rt: RRuntime) -> str:
    """The R line that installs ``name`` into ``rt`` (into the env library for a conda R)."""
    lib = ", lib = .Library" if rt.conda_prefix else ""
    where = f"in R {rt.version} ({rt.source}: {tilde(rt.rscript)})"
    if name in _BIOC:
        return f'BiocManager::install("{name}"{lib})  # {where}'
    if name in _CRAN:
        return f'install.packages("{name}"{lib})  # {where}'
    if name in _GITHUB:
        return f'remotes::install_github("{_GITHUB[name]}"{lib})  # {where}'
    return f"# {name} is not on CRAN or Bioconductor as far as grnbench knows: see its own installation notes ({where})"


def _have_module(name: str) -> bool:
    try:
        return importlib.util.find_spec(name) is not None
    except (ImportError, ValueError):
        return False


def _py_version(name: str) -> str | None:
    for dist in dict.fromkeys((_PIP_NAMES.get(name, name), name)):
        try:
            return importlib.metadata.version(dist)
        except importlib.metadata.PackageNotFoundError:
            continue
    return None


def requirement_status(req: str, rt: RRuntime | None | EllipsisType = ...) -> ReqStatus:
    """Whether a knowledge-base requirement (``py:module``, ``r:package``, ``cli:tool``) is met here.

    A missing one has the reason ``missing: <req>``, the knowledge base's own text; an R package built for
    another R version has ``built for R x.y``. The detail is in ``install_hint``, which also names any R
    candidate that was skipped. ``rt`` is the R to ask (``...``: :func:`find_rscript`; None: there is no R).
    The knowledge base asks it for ``r:`` requirements when the bridge is on (``kb.r_bridge_on``)."""
    kind, _, name = req.partition(":")
    missing = f"missing: {req}"
    if kind == "py":
        if _have_module(name):
            return ReqStatus(True, "", "")
        return ReqStatus(False, missing, f"pip install {_PIP_NAMES.get(name, name)}")
    if kind == "cli":
        if shutil.which(name):
            return ReqStatus(True, "", "")
        return ReqStatus(False, missing, f"put `{name}` on PATH (install it in its own environment; see its documentation)")
    if kind != "r":
        return ReqStatus(False, missing, f"unknown requirement kind {kind!r} (expected py:, r: or cli:)")
    notes: list[str] = []
    if rt is ...:
        try:
            rt, notes = _find()
        except RBridgeError as exc:
            return ReqStatus(False, missing, str(exc))
    if rt is None:
        return ReqStatus(False, missing, _no_r_hint(notes))
    try:
        info = probe(rt, _packages_for(name))
    except RBridgeError as exc:
        return ReqStatus(False, missing, f"R does not answer: {exc}" + _skipped_lines(notes))
    version, problem = _package_state(info, name)
    if version and not problem:
        return ReqStatus(True, "", "")
    if problem:
        return ReqStatus(False, problem, f"# this R is {info['version']}: reinstall it\n"
                         + _r_install_hint(name, rt) + _skipped_lines(notes))
    return ReqStatus(False, missing, _r_install_hint(name, rt) + _skipped_lines(notes))


# -- running R and other tools ------------------------------------------------------------ #
def _check_seed(seed: Any) -> int:
    if isinstance(seed, bool) or not isinstance(seed, numbers.Integral):
        raise TypeError(f"seed must be an int, got {type(seed).__name__}")
    if not 0 <= int(seed) <= 2**31 - 1:
        raise ValueError(f"seed must be in [0, 2**31 - 1] (R's integer range), got {seed}")
    return int(seed)


def _json_default(value: Any):
    if isinstance(value, os.PathLike):
        return os.fspath(value)
    if type(value).__module__ == "numpy" and hasattr(value, "tolist"):
        return value.tolist()       # numpy integers, floats, bools and arrays
    raise TypeError(f"{type(value).__name__} is not a JSON value")


def _expected_files(work: Path, expect: Mapping[str, str]) -> dict[str, Path]:
    """name -> file inside ``work`` that R must write. Refused before anything is touched: absolute paths,
    paths that leave ``work`` (``..``, links), the bridge's own files and existing folders."""
    reserved = {_norm(work / n) for n in ("params.json", "session.json", "r.log")}
    files: dict[str, Path] = {}
    for name, rel in expect.items():
        given = Path(rel)
        if given.is_absolute() or given.anchor:
            raise RBridgeError(f"expect[{name!r}] = {rel!r}: give a path relative to the work folder")
        target = (work / given).resolve()
        if target == work:
            raise RBridgeError(f"expect[{name!r}] = {rel!r} is the work folder itself; name a file inside it")
        if not target.is_relative_to(work):
            raise RBridgeError(f"expect[{name!r}] = {rel!r} points outside the work folder {work}")
        if _norm(target) in reserved:
            raise RBridgeError(f"expect[{name!r}] = {rel!r} is one of the bridge's own files "
                               "(params.json, session.json, r.log)")
        if target.is_dir():
            raise RBridgeError(f"expect[{name!r}] = {rel!r} is a folder; an expected output is a file")
        files[name] = target
    return files


def run_r(script: str | Path, *, inputs: Mapping[str, Path] | None = None, params: Mapping | None = None,
          seed: int, workdir: Path, timeout_s: int = 3600, rt: RRuntime | None = None,
          expect: Mapping[str, str] | None = None) -> RResult:
    """Run the R ``script`` in ``workdir`` through ``r/_wrapper.R`` and return what it wrote.

    R starts in ``workdir`` with only the wrapper on its command line; everything else travels in
    ``workdir/params.json`` (UTF-8, so paths with any characters survive on Windows): ``script``, ``seed``,
    ``params``, ``inputs`` (name -> absolute path), ``expect`` (name -> file it must write, relative to
    ``workdir``) and ``workdir``. The wrapper reads it with ``jsonlite::fromJSON(simplifyVector = FALSE)``, so
    in R JSON arrays are lists and scalars are scalars. It sets ``RNGkind("Mersenne-Twister", "Inversion",
    "Rejection")`` and ``set.seed(seed)``, sources ``script`` with ``params``, ``inputs``, ``workdir`` and
    ``expect`` defined, and records R and package versions, the RNG kind before and after and the method's own
    ``set.seed`` calls in ``session.json``. Everything R prints goes to ``workdir/r.log``. Stale copies of
    ``session.json`` and of the expected files are removed first, so only this run's files count. Raises
    :class:`RBridgeError` on a timeout or an interrupt (the process tree is killed), a non-zero exit (with
    the end of the log), a missing expected file or a missing or invalid ``session.json``; ``rt=None`` uses
    :func:`find_rscript` and records the candidates it skipped in :attr:`RResult.skipped`."""
    seed = _check_seed(seed)
    notes: list[str] = []
    if rt is None:
        rt, notes = _find()
        if rt is None:
            raise RBridgeError("this step needs R and none was found.\n" + _no_r_hint(notes))
    script_path = Path(script).expanduser().resolve()
    if not script_path.is_file():
        raise RBridgeError(f"no R script at {script_path}")
    work = Path(workdir).expanduser().resolve()
    expect = {str(k): str(v) for k, v in (expect or {}).items()}
    outputs = _expected_files(work, expect)
    payload = {"script": str(script_path), "seed": seed, "params": dict(params or {}),
               "inputs": {str(k): str(Path(v).expanduser().resolve()) for k, v in (inputs or {}).items()},
               "expect": expect, "workdir": str(work)}
    try:
        text = json.dumps(payload, indent=2, default=_json_default, allow_nan=False)
    except (TypeError, ValueError) as exc:
        raise RBridgeError(f"params must be JSON values (numbers, text, lists, dicts, paths; no NaN or "
                           f"infinity): {exc}") from None
    work.mkdir(parents=True, exist_ok=True)
    session_path = work / "session.json"
    for stale in (session_path, *outputs.values()):
        if stale.is_file():
            try:
                stale.unlink()
            except OSError as exc:
                raise RBridgeError(f"cannot remove the previous {stale.name} in {work}: {exc}") from None
    (work / "params.json").write_text(text, encoding="utf-8")
    log = work / "r.log"
    started = time.monotonic()
    try:
        done = _run_logged([rt.rscript, "--vanilla", str(WRAPPER)], log=log, timeout_s=timeout_s,
                           env=rt.subprocess_env(), cwd=work)
    except OSError as exc:
        raise RBridgeError(f"cannot start {tilde(rt.rscript)}: {exc}") from None
    seconds = time.monotonic() - started
    if done.code is None:
        raise RBridgeError(f"R timed out after {timeout_s} s; log: {log}{done.orphan}")
    if done.code != 0:
        raise RBridgeError(f"R {_exit_text(done.code)}; log: {log}\n{_log_tail(log)}")
    absent = [f"{name} ({rel})" for name, rel in expect.items() if not outputs[name].is_file()]
    if absent:
        raise RBridgeError(f"R finished but did not write {', '.join(absent)}; log: {log}\n{_log_tail(log)}")
    try:
        session = json.loads(session_path.read_text(encoding="utf-8"))
        why = "not a JSON object"
    except (OSError, ValueError) as exc:
        session, why = None, f"{type(exc).__name__}: {exc}"
    if not isinstance(session, dict):
        raise RBridgeError(f"R finished but session.json is missing or invalid ({why}); log: {log}\n{_log_tail(log)}")
    return RResult(outputs, log, session, seconds, rt, tuple(notes))


def _tool_version(cmd: Sequence[str], env: Mapping[str, str] | None) -> str | None:
    try:
        done = _run_captured([str(c) for c in cmd], timeout_s=60, env=env)
    except OSError:
        return None
    lines = [line.strip() for line in f"{done.out}\n{done.err}".splitlines() if line.strip()]
    return lines[0] if done.code == 0 and lines else None


def run_external(cmd: Sequence[str], *, workdir: Path, timeout_s: int, env: Mapping[str, str] | None = None,
                 version_cmd: Sequence[str] | None = None) -> ExternalResult:
    """Run a command-line tool in ``workdir`` under the same rules as R: output in ``<tool>-<k>.log`` (the first
    free k, so every call keeps its own log), a timeout or an interrupt kills the process tree, a non-zero
    exit raises :class:`RBridgeError` with the end of the log. A ``cmd[0]`` with a folder in it is made
    absolute first (the tool runs in ``workdir``). ``env`` is the whole environment (default: this process's);
    ``version_cmd`` (e.g. ``[tool, "--version"]``) gives the first line it prints as ``version``."""
    cmd = [str(c) for c in cmd]
    if not cmd:
        raise ValueError("run_external needs a command")
    if os.sep in cmd[0] or (os.altsep and os.altsep in cmd[0]):
        cmd[0] = os.path.abspath(os.path.expanduser(cmd[0]))
    work = Path(workdir).expanduser().resolve()
    work.mkdir(parents=True, exist_ok=True)
    stem = Path(cmd[0]).stem or "external"
    k = 1
    while (work / f"{stem}-{k}.log").exists():
        k += 1
    log = work / f"{stem}-{k}.log"
    version = _tool_version(version_cmd, env) if version_cmd else None
    started = time.monotonic()
    try:
        done = _run_logged(cmd, log=log, timeout_s=timeout_s, env=env, cwd=work)
    except OSError as exc:
        raise RBridgeError(f"cannot start {cmd[0]}: {exc}") from None
    seconds = time.monotonic() - started
    if done.code is None:
        raise RBridgeError(f"{stem} timed out after {timeout_s} s; log: {log}{done.orphan}")
    if done.code != 0:
        raise RBridgeError(f"{stem} {_exit_text(done.code)}; log: {log}\n{_log_tail(log)}")
    return ExternalResult(done.code, log, seconds, version)


# -- environment report ------------------------------------------------------------------- #
def _row(tool: str, found: bool, version: str | None, source: str | None, used_for: Sequence[str],
         install: str) -> dict:
    return {"tool": tool, "found": bool(found), "version": version, "source": source,
            "used_for": list(used_for), "install": install}


def _internal(exc: BaseException) -> str:
    return _scrub(f"internal error ({type(exc).__name__}): {exc}")


def environment_report(explicit: str | None = None, *, fresh: bool = False) -> list[dict]:
    """One row per optional tool: R itself, each ``r:`` package (with jsonlite), each ``py:`` and ``cli:``
    requirement of the knowledge base. Keys: ``tool``, ``found``, ``version``, ``source``, ``used_for``
    (``decision=option`` pairs) and ``install`` (what to run when it is missing; for R also the candidates
    that were skipped). Never raises: an R problem becomes the R row's text, anything unexpected an
    "internal error" naming the exception type. ``fresh=True`` asks R and conda again instead of using the
    cache (``grnbench env`` and the interface's Check R do)."""
    return environment_status(explicit, fresh=fresh)[0]


def environment_status(explicit: str | None = None, *, fresh: bool = False) -> tuple[list[dict], bool]:
    """:func:`environment_report` plus whether looking for R failed: True when it raised (an R that was
    named does not work, ``GRNBENCH_R_ENV`` names no env, an existing env has no working R, or an internal
    error), False when R was found or is simply absent or turned off. ``grnbench env`` exits 1 on True.
    ``fresh=True`` also makes the knowledge base look R up again next time (``kb._r_runtime``)."""
    if fresh:
        from grnbench import kb
        kb._R_RUNTIME.clear()
    try:
        uses = _kb_uses()
        reqs = _r_requirements()
    except Exception as exc:  # noqa: BLE001 - a status command reports problems, it does not stop on them
        return [_row("R", False, None, None, [], _internal(exc))], True
    r_uses = sorted({p for req, pairs in uses.items() if req.startswith("r:") for p in pairs})
    rt: RRuntime | None = None
    info: dict | None = None
    notes: list[str] = []
    problem = ""
    try:
        rt, notes = _find(explicit, fresh=fresh)
        if rt is not None:
            info = probe(rt, reqs)          # the answer _find just stored: no second R start
    except RBridgeError as exc:
        problem = str(exc)
    except Exception as exc:  # noqa: BLE001
        problem = _internal(exc)
    rows = []
    if rt is not None and info is not None:
        rows.append(_row("R", True, rt.version, f"{rt.source}: {tilde(rt.rscript)}", r_uses,
                         _skipped_lines(notes).lstrip("\n")))
    else:
        rows.append(_row("R", False, None, None, r_uses, problem or _no_r_hint(notes)))
    for pkg in reqs:
        req = f"r:{pkg}"
        pairs = uses.get(req) or ("R bridge (reads params.json, writes session.json)",)
        if rt is not None and info is not None:
            version, bad = _package_state(info, pkg)
            ok = bool(version) and not bad
            install = "" if ok else ((f"# {bad}; this R is {info['version']}: reinstall it\n" if bad else "")
                                     + _r_install_hint(pkg, rt))
            rows.append(_row(req, ok, version, rt.source, pairs, install))
        else:
            rows.append(_row(req, False, None, None, pairs, "needs R first (see the R row)"))
    for req in sorted((r for r in uses if r.startswith(("py:", "cli:"))), key=lambda r: (r.startswith("cli:"), r)):
        name = req.split(":", 1)[1]
        try:
            st = requirement_status(req)
            if req.startswith("py:"):
                rows.append(_row(req, st.ok, _py_version(name) if st.ok else None, "python", uses[req], st.install_hint))
            else:
                rows.append(_row(req, st.ok, None, tilde(shutil.which(name)) if st.ok else None, uses[req],
                                 st.install_hint))
        except Exception as exc:  # noqa: BLE001
            rows.append(_row(req, False, None, None, uses[req], _internal(exc)))
    return rows, bool(problem)
