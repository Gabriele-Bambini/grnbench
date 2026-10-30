"""Repository and project-root discovery.

Finds the grnbench repository root by walking up from an explicit path or the
current working directory, checking for specific project markers.
"""
from __future__ import annotations

import tomllib
from pathlib import Path


def _has_grnbench_identity(candidate: Path) -> bool:
    """Return whether *candidate* identifies itself as the grnbench project."""
    if (candidate / "src" / "grnbench" / "__init__.py").is_file():
        return True

    pyproject = candidate / "pyproject.toml"
    if not pyproject.is_file():
        return False
    try:
        data = tomllib.loads(pyproject.read_text(encoding="utf-8"))
    except (OSError, tomllib.TOMLDecodeError):
        return False
    return data.get("project", {}).get("name") == "grnbench"


def find_project_root(start: str | Path | None = None) -> Path | None:
    """Locate the grnbench checkout root.

    Walks parent directories starting from ``start`` (or ``Path.cwd()``) and looks
    for the team state file plus a grnbench-specific identity marker: either the
    source package or ``project.name = "grnbench"`` in ``pyproject.toml``.
    Returns the resolved root Path, or None if no checkout is found. Generic
    repositories and arbitrary ``AGENTS.md`` files are intentionally not accepted.
    """
    p = Path(start).resolve() if start is not None else Path.cwd().resolve()
    if p.is_file():
        p = p.parent

    for candidate in (p, *p.parents):
        if (candidate / ".projectos" / "STATE.md").is_file() and _has_grnbench_identity(candidate):
            return candidate
    return None
