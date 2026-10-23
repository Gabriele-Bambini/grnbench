"""Markdown run report + JSON manifest generation."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import os
from datetime import datetime, timezone
from pathlib import Path

STATUS_ICON = {"pass": "PASS", "warn": "WARN", "fail": "FAIL"}


def file_fingerprint(path: str | Path, full_below_mb: int = 300) -> dict:
    """Full streamed SHA256; directories use a sorted recursive inventory.

    ``full_below_mb`` remains accepted for compatibility but no longer truncates
    hashing. Symlinks are rejected to avoid reading outside the named input.
    """
    p = Path(path)
    if p.is_symlink() or (hasattr(p, 'is_junction') and p.is_junction()):
        raise ValueError(f"Linked input is not supported: {p}")
    if p.is_dir():
        files = []
        for base, dirs, names in os.walk(p, followlinks=False):
            for name in dirs + names:
                item = Path(base) / name
                linked = item.is_symlink() or (hasattr(item, 'is_junction') and item.is_junction())
                # Python <3.12 lacks Path.is_junction; reparse directories also
                # must be rejected before os.walk can descend into them.
                attrs = getattr(item.lstat(), 'st_file_attributes', 0)
                if linked or (item.is_dir() and attrs & 0x400):
                    raise ValueError(f"Linked directory entry is not supported: {item}")
            for name in names:
                item = Path(base) / name
                fp = file_fingerprint(item)
                fp['path'] = item.relative_to(p).as_posix()
                files.append(fp)
        files.sort(key=lambda item: item['path'])
        digest = hashlib.sha256(json.dumps(files, sort_keys=True,
                                          separators=(',', ':')).encode()).hexdigest()
        return {'path': str(p), 'kind': 'directory', 'files': files,
                'size_bytes': sum(f['size_bytes'] for f in files), 'sha256': digest}
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b''):
            h.update(chunk)
    return {'path': str(p), 'size_bytes': p.stat().st_size, 'sha256': h.hexdigest()}


def sequence_identity(values) -> dict:
    """Exact ordered identifiers with an unambiguous content digest."""
    values = [str(v) for v in values]
    return {'count': len(values), 'values': values,
            'sha256': hashlib.sha256(json.dumps(values, ensure_ascii=False,
                                               separators=(',', ':')).encode()).hexdigest()}


def require_empty_output(path: str | Path, *, allow_failed_diagnostic: bool = True) -> Path:
    """Refuse output reuse rather than mix artifacts from different runs.

    A lone schema-v1 failure diagnostic is the one safe exception: it contains
    no run artifact and exists specifically to explain an attempt that stopped
    before producing one.  A retry may replace it, while every other existing
    file still makes reuse fail loudly.  ``allow_failed_diagnostic=False``
    drops that exception for outputs that never write a diagnostic (``demo``),
    which would otherwise leave the stale one next to their files.
    """
    out = Path(path)
    if out.exists():
        entries = list(out.iterdir())
        retryable = False
        if allow_failed_diagnostic and len(entries) == 1 and entries[0].name == "diagnostics.json":
            try:
                diagnostic = json.loads(entries[0].read_text(encoding="utf-8"))
                retryable = (
                    diagnostic.get("schema") == "grnbench.diagnostics"
                    and diagnostic.get("schema_version") == 1
                    and diagnostic.get("outcome") == "failed"
                )
            except (OSError, ValueError, AttributeError):
                retryable = False
        if entries and not retryable:
            raise ValueError(f"Output directory must be empty; choose a new path: {out}")
    out.mkdir(parents=True, exist_ok=True)
    return out


def thread_count(value: str) -> int:
    """argparse type for --threads: a CPU budget of at least 1."""
    n = int(value)
    if n < 1:
        raise argparse.ArgumentTypeError(f"must be at least 1, got {n}")
    return n


def _dist_version(name: str) -> str:
    """Installed version of a distribution, or 'not installed'/'unknown'.

    Uses importlib.metadata rather than ``module.__version__``: scanpy and
    anndata deprecated the attribute, and querying it emits FutureWarnings all
    over the test suite.
    """
    from importlib.metadata import PackageNotFoundError, version
    try:
        return version(name)
    except PackageNotFoundError:
        return "not installed"
    except Exception:  # noqa: BLE001 - metadata backends vary; never break a run
        return "unknown"


def environment_block() -> dict:
    import grnbench  # own package: __version__ is authoritative even from a source tree
    from importlib.metadata import distributions
    return {
        "grnbench": grnbench.__version__,
        "python": platform.python_version(),
        "scanpy": _dist_version("scanpy"),
        "anndata": _dist_version("anndata"),
        "numpy": _dist_version("numpy"),
        "pandas": _dist_version("pandas"),
        "scipy": _dist_version("scipy"),
        "arboreto": _dist_version("arboreto"),
        "platform": platform.platform(),
        "installed_distributions": dict(sorted(
            (d.metadata['Name'], d.version) for d in distributions() if d.metadata['Name'])),
    }


def build_manifest(input_fp: dict, ingest: dict, qc: dict, steps: list[dict],
                   strategy: str, recommendations: list[dict],
                   outputs: dict, extra: dict | None = None) -> dict:
    manifest = {
        "grnbench_run": True,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input": input_fp,
        "ingest": ingest,
        "qc": qc,
        "steps": steps,
        "strategy_chosen": strategy,
        "recommendations": recommendations,
        "outputs": outputs,
        "environment": environment_block(),
    }
    if extra:
        manifest.update(extra)
    return manifest


def render_report(manifest: dict, diagnostics: dict | None = None) -> str:
    """Human-readable markdown report from the run manifest."""
    qc = manifest["qc"]
    m = qc["metrics"]
    lines: list[str] = []
    add = lines.append

    add("# GRN-Bench preprocessing report")
    add("")
    add(f"- Generated: {manifest['created_utc']}")
    add(f"- Input: `{manifest['input']['path']}` "
        f"({manifest['input'].get('size_bytes', 0) / 1e6:.1f} MB)")
    add(f"- Strategy: **{manifest['strategy_chosen']}**")
    add("")
    add("## Ingest")
    add("")
    ing = manifest["ingest"]
    add(f"- Format: `{ing['fmt']}` - counts source: `{ing['counts_source']}`"
        + (" - table transposed" if ing.get("transposed") else ""))
    if ing.get("n_duplicates_merged"):
        add(f"- Duplicate genes merged by sum: **{ing['n_duplicates_merged']}**")
    if ing.get("version_suffixes_stripped"):
        add(f"- Ensembl version suffixes stripped: {ing['version_suffixes_stripped']}")
    for note in ing.get("notes", []):
        add(f"- Note: {note}")
    cp = ing.get("count_plausibility") or {}
    if "cv_cell_totals" in cp or "share_ones" in cp:
        parts = []
        if "cv_cell_totals" in cp:
            parts.append(f"CV of cell totals {cp['cv_cell_totals']:.3g}")
        if "share_ones" in cp:
            parts.append(f"1s are {cp['share_ones']:.1%} of non-zero values, "
                         f"max {cp['max_value']:g}")
        add(f"- Count plausibility: {'; '.join(parts)} - flags: "
            f"{', '.join(cp.get('flags') or []) or 'none'}")
    add("")
    add("## QC metrics")
    add("")
    add("| metric | value |")
    add("|---|---|")
    for k, v in m.items():
        add(f"| {k} | {v:.4g} |" if isinstance(v, float) else f"| {k} | {v} |")
    add("")
    add("## QC gates")
    add("")
    add("| gate | status | detail |")
    add("|---|---|---|")
    for g in qc["gates"]:
        add(f"| {g['name']} | **{STATUS_ICON[g['status']]}** | {g['detail']} |")
    add("")
    warn_list = manifest.get("warnings") or []
    if warn_list:
        add("## Warnings")
        add("")
        for w in warn_list:
            add(f"- **{w}**")
        add("")
    if manifest.get("organism_hint"):
        add(f"- Organism hint from symbols: `{manifest['organism_hint']}`")
    gid = manifest.get("gene_id_type") or {}
    if gid:
        add(f"- Gene-ID verdict: `{gid.get('verdict')}` "
            f"(ensembl frac {gid.get('ensembl_frac', 0):.1%})")
    add("## Strategy recommendation")
    add("")
    add("| rank | strategy | reason |")
    add("|---|---|---|")
    for r in manifest["recommendations"]:
        add(f"| {r['rank']} | `{r['strategy']}` | {r['reason']} |")
    add("")
    add("## Steps")
    add("")
    for s in manifest["steps"]:
        desc = s.get("detail", "")
        dims = ""
        if "n_cells_after" in s:
            dims = f" -> {s['n_cells_after']} cells x {s['n_genes_after']} genes"
        add(f"- **{s['step']}**{dims} {desc}".rstrip())
    add("")
    add("## Outputs")
    add("")
    for k, v in manifest["outputs"].items():
        add(f"- {k}: `{v}`")
    add("")
    add("## Hard rules honored")
    add("")
    add("- Raw counts verified at ingest; no imputation applied anywhere.")
    add("- Every step recorded above with parameters; input fingerprinted in manifest.json.")
    add("")
    from grnbench.diagnostics import render_diagnostics_markdown
    lines.extend(render_diagnostics_markdown(diagnostics).splitlines())
    return "\n".join(lines)


def write_json(path: Path, obj: dict) -> None:
    path.write_text(json.dumps(obj, indent=2, default=str), encoding="utf-8")
