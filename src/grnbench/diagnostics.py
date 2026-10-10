"""Portable, privacy-bounded diagnostics shared by the CLI, reports and UI.

The normal run manifests deliberately contain local paths and detailed
provenance.  ``diagnostics.json`` serves a different purpose: it is the small
artifact that can be copied between machines and rendered consistently.  It
therefore stores only outcomes, bounded messages, preprocessing choices and
cell/gene dimensions.  Raw values, barcodes, paths and machine identities are
forbidden by :func:`validate_portable`.
"""
from __future__ import annotations

import json
import getpass
import math
import numbers
import re
import socket
import warnings
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any


SCHEMA = "grnbench.diagnostics"
SCHEMA_VERSION = 1
MAX_MESSAGE_CHARS = 300

_FORBIDDEN_KEYS = {
    "barcode", "barcodes", "cli", "data", "email", "expression", "host",
    "hostname", "input", "manifest", "matrix", "out_dir", "output", "outputs",
    "path", "paths", "processed", "student_id", "user", "username", "values",
}
_QUOTED_PATH = re.compile(
    r"(?P<quote>['\"])(?:file://|[A-Za-z]:[\\/]|\\\\|~[\\/]|/)[^'\"]+(?P=quote)",
    re.IGNORECASE,
)
_PATH_TOKEN = re.compile(
    r"(?<![\w.+:/-])(?:file://[^\s'\"<>]+|[A-Za-z]:[\\/][^\s'\"<>]+|"
    r"\\\\[^\\\s'\"<>]+\\[^\s'\"<>]+|~[\\/][^\s'\"<>]+|"
    r"/[^\s'\"<>]+)",
    re.IGNORECASE,
)
_RATIO_SLASH = re.compile(r"(?<=\d\s)/(?=\d+(?:\s|$|[.,;:)]))")
_RATIO_SENTINEL = "__GRNBENCH_RATIO_SLASH__"
_IDENTITY_TOKEN = re.compile(
    r"(?<!\w)[A-Za-z0-9][\w.+-]*@[A-Za-z0-9][\w.-]*\.[A-Za-z]{2,}(?!\w)")
_OUTCOMES = {"unknown", "ok", "failed", "partial"}


class DiagnosticPrivacyError(ValueError):
    """Raised when a portable diagnostic contains forbidden information."""


def _known_identities() -> tuple[str, ...]:
    values = []
    for getter in (getpass.getuser, socket.gethostname):
        try:
            value = str(getter()).strip()
        except Exception:  # noqa: BLE001 - identity discovery is best effort
            continue
        if value:
            values.append(value)
    return tuple(values)


def _replace_literal(text: str, value: Any, replacement: str, *, boundary: bool) -> str:
    raw = str(value or "").strip()
    if not raw:
        return text
    pattern = re.escape(raw)
    if boundary:
        pattern = rf"(?<![\w.-]){pattern}(?![\w.-])"
    return re.sub(pattern, replacement, text, flags=re.IGNORECASE)


def safe_message(message: Any, limit: int = MAX_MESSAGE_CHARS,
                 *, sensitive_values: tuple[Any, ...] | list[Any] = ()) -> str:
    """Collapse, redact and bound an exception or diagnostic message."""
    text = _RATIO_SLASH.sub(_RATIO_SENTINEL, str(message))
    path_values = [str(value).strip() for value in sensitive_values
                   if value is not None and str(value).strip()]
    try:
        home = str(Path.home()).strip()
        if home:
            path_values.append(home)
    except Exception:  # noqa: BLE001 - redaction must not become a run failure
        pass

    # Replace exact run-specific values before generic path regexes.  Otherwise
    # an unquoted Windows path containing spaces is shortened at its first space
    # and the remaining private suffix can no longer match the original value.
    for value in sorted(set(path_values), key=len, reverse=True):
        text = _replace_literal(text, value, "<path>", boundary=True)
    text = re.sub(r"(['\"])<path>\1", "<path>", text)
    # A sensitive relative directory can be followed by a filename in an
    # exception (``<path>/matrix.mtx.gz``); keep the whole location redacted
    # before the generic matcher sees the leading slash as a second path.
    text = re.sub(r"<path>(?:[\\/][^\s'\"<>]+)+", "<path>", text)
    text = _QUOTED_PATH.sub("<path>", text)
    text = _PATH_TOKEN.sub("<path>", text)
    text = re.sub(r"(?:<path>){2,}", "<path>", text)
    for value in _known_identities():
        text = _replace_literal(text, value, "<identity>", boundary=True)
    text = _IDENTITY_TOKEN.sub("<identity>", text)
    text = text.replace(_RATIO_SENTINEL, "/")
    text = " ".join(text.split())
    if len(text) > limit:
        text = text[: max(0, limit - 1)].rstrip() + "…"
    return text


def _outcome(status: Any) -> str:
    raw = str(status).strip().lower().replace("-", "_")
    if raw in {"ok", "complete", "success", "succeeded"} or raw.startswith("ok "):
        return "ok"
    if raw in {"partial", "partial_failure", "running"}:
        return "partial"
    if (raw in {"failed", "failure", "error", "no edges"}
            or "fail" in raw or "error" in raw):
        return "failed"
    return "unknown"


def new_document(*, outcome: str = "unknown") -> dict:
    """Return an empty schema-v1 document."""
    return {
        "schema": SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "outcome": _outcome(outcome),
        "status": safe_message(outcome),
        "prep": {},
        "screen": {},
        "models": {},
        "usable_candidates": [],
        "selected_candidate": None,
    }


def failure(category: str, message: Any, *, sensitive_values=()) -> dict:
    return {"category": str(category),
            "message": safe_message(message, sensitive_values=sensitive_values)}


def classify_failure(message: Any, *, stage: str, category: str | None = None,
                     sensitive_values=()) -> dict:
    """Build a stable failure object from an exception or stored error string."""
    text = str(message)
    if category is None:
        low = text.lower()
        if stage == "prep" and "qc gate failure" in low:
            category = "qc_gate"
        elif stage == "prep" and ("prep refused" in low or "preperror" in low):
            category = "prep_refused"
        elif stage == "screen":
            category = "screening_error"
        elif stage == "model" and "unavailable" in low:
            category = "unavailable_dependency"
        else:
            category = f"{stage}_error"
    return failure(category, text, sensitive_values=sensitive_values)


def _dimensions(cells: Any = None, genes: Any = None) -> dict | None:
    out = {}
    for name, value in (("cells", cells), ("genes", genes)):
        if value is None or isinstance(value, bool):
            continue
        try:
            numeric = float(value)
        except (TypeError, ValueError, OverflowError):
            continue
        if math.isfinite(numeric) and numeric >= 0 and numeric.is_integer():
            out[name] = int(numeric)
    return out or None


def _portable_choices(args: dict | None) -> dict:
    """Select only behaviour-changing, path-free choices from a manifest."""
    allowed = (
        "strategy", "n_hvg", "group_key", "label_key", "stratify_key", "batch_key",
        "min_group_cells", "downcast", "subsample_cells", "orientation",
        "drop_doublets", "force", "seed", "subset_key", "subset_value", "cell_calling",
        "exclude_qc_genes", "cells_per_metacell", "thresholds",
    )
    out = {}
    for key in allowed:
        if not args or args.get(key) is None:
            continue
        value = args[key]
        out[key] = safe_message(value) if isinstance(value, str) else value
    return out


def _label(value: Any) -> str:
    """Portable dictionary key for a candidate, context or model name."""
    return safe_message(value, limit=120) or "unknown"


def context_label(*, subset_key: Any = None, subset_value: Any = None,
                  explicit: Any = None) -> str:
    """Return the context name shared by prep, screening and model diagnostics."""
    if explicit is not None:
        return str(explicit)
    if subset_key is not None and subset_value is not None:
        return f"{subset_key}={subset_value}"
    return "all_cells"


def prep_entry_from_manifest(manifest: dict) -> dict:
    """Convert a successful preprocessing manifest into a portable entry."""
    metrics = (manifest.get("qc") or {}).get("metrics") or {}
    before_cells = metrics.get("n_cells")
    before_genes = metrics.get("n_genes")
    for raw in manifest.get("steps") or []:
        if raw.get("step") == "cell_calling" and raw.get("n_barcodes") is not None:
            before_cells = raw.get("n_barcodes")
            break
    before = _dimensions(before_cells, before_genes)
    steps = []
    after = before
    for raw in manifest.get("steps") or []:
        dims = _dimensions(raw.get("n_cells_after"), raw.get("n_genes_after"))
        if dims:
            after = dims
        item = {"stage": str(raw.get("step", "unknown"))}
        if dims:
            item["after"] = dims
        if raw.get("detail"):
            item["message"] = safe_message(raw["detail"])
        steps.append(item)
    choices = _portable_choices(manifest.get("args"))
    chosen_strategy = manifest.get("strategy_chosen")
    if chosen_strategy is not None:
        if "strategy" in choices:
            choices["requested_strategy"] = choices["strategy"]
        choices["strategy"] = safe_message(chosen_strategy)
    return {
        "outcome": "ok",
        "failure": None,
        "choices": choices,
        "dimensions": {"before": before, "after": after},
        "steps": steps,
    }


def failed_entry(message: Any, *, stage: str, choices: dict | None = None,
                 category: str | None = None, sensitive_values=()) -> dict:
    portable_choices = _portable_choices(choices)
    if stage == "prep" and "strategy" in portable_choices:
        portable_choices["requested_strategy"] = portable_choices.pop("strategy")
    return {
        "outcome": "failed",
        "failure": classify_failure(message, stage=stage, category=category,
                                    sensitive_values=sensitive_values),
        "choices": portable_choices,
        "dimensions": {"before": None, "after": None},
        "steps": [],
    }


def prep_document(manifest: dict, *, candidate: str = "default",
                  context: str = "all_cells") -> dict:
    doc = new_document(outcome="ok")
    candidate, context = _label(candidate), _label(context)
    doc["prep"] = {candidate: {context: prep_entry_from_manifest(manifest)}}
    doc["usable_candidates"] = [candidate]
    doc["selected_candidate"] = candidate
    validate_portable(doc)
    return doc


def failed_prep_document(message: Any, *, choices: dict | None = None,
                         candidate: str = "default", context: str = "all_cells",
                         category: str | None = None, sensitive_values=()) -> dict:
    doc = new_document(outcome="failed")
    candidate, context = _label(candidate), _label(context)
    doc["prep"] = {candidate: {context: failed_entry(
        message, stage="prep", choices=choices, category=category,
        sensitive_values=sensitive_values)}}
    validate_portable(doc)
    return doc


def _read_json(path: Path) -> dict | None:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return value if isinstance(value, dict) else None


def run_document(run_dir: str | Path, *, outcome: str, stages: list[dict] | None = None) -> dict:
    """Build the one-shot run diagnostic from artifacts already written on disk."""
    root = Path(run_dir)
    doc = new_document(outcome=outcome)
    prep_diagnostic = _read_json(root / "prep" / "diagnostics.json")
    if prep_diagnostic and isinstance(prep_diagnostic.get("prep"), dict):
        doc["prep"] = prep_diagnostic["prep"]
        doc["usable_candidates"] = list(prep_diagnostic.get("usable_candidates") or [])
        doc["selected_candidate"] = prep_diagnostic.get("selected_candidate")
    else:
        prep = _read_json(root / "prep" / "manifest.json")
        if prep:
            doc["prep"] = {"default": {"all_cells": prep_entry_from_manifest(prep)}}
            doc["usable_candidates"] = ["default"]
            doc["selected_candidate"] = "default"
        else:
            prep_stage = next((x for x in (stages or []) if x.get("stage") == "prep"), None)
            if prep_stage and prep_stage.get("exit_code"):
                doc["prep"] = {"default": {"all_cells": failed_entry(
                    f"preprocessing exited with code {prep_stage['exit_code']}", stage="prep",
                    category="prep_failed")}}

    infer_diagnostic = _read_json(root / "infer" / "diagnostics.json")
    if infer_diagnostic and isinstance(infer_diagnostic.get("models"), dict):
        doc["models"] = infer_diagnostic["models"]
    else:
        infer = _read_json(root / "infer" / "inference_manifest.json")
        if infer:
            doc["models"] = model_document(infer)["models"]
    validate_portable(doc)
    return doc


def model_document(manifest: dict, *, context: str = "all_cells") -> dict:
    """Convert a standalone inference manifest into portable model outcomes."""
    status = str(manifest.get("status", "unknown"))
    doc = new_document(outcome=status)
    context = _label(context)
    gene_count = (manifest.get("effective_genes") or {}).get("count")
    for model in sorted(manifest.get("runs") or {}):
        model_name = _label(model)
        doc["models"].setdefault(model_name, {})[context] = {
            "outcome": "ok", "failure": None,
            "choices": {"model": model_name},
            "dimensions": {"before": _dimensions(genes=gene_count),
                           "after": _dimensions(genes=gene_count)},
            "steps": [],
        }
    skipped = manifest.get("methods_skipped") or {}
    errors = manifest.get("errors") or {}
    for model, message in sorted({**skipped, **errors}.items()):
        model_name = _label(model)
        doc["models"].setdefault(model_name, {})[context] = failed_entry(
            message, stage="model", category=(
                "unavailable_dependency" if model in skipped else "model_error"))
        doc["models"][model_name][context]["choices"] = {"model": model_name}
    validate_portable(doc)
    return doc


def screen_document(plans: dict, records: dict, rows, *, selected: str | None) -> dict:
    """Convert candidate/context screening records into the shared schema."""
    doc = new_document(outcome="ok")
    for candidate, recs in records.items():
        plan = plans[candidate]
        candidate_name = _label(candidate)
        for rec in recs:
            context = _label(rec.get("context", "all_cells"))
            manifest = _read_json(Path(rec["manifest"])) if rec.get("manifest") else None
            prep_entry = (prep_entry_from_manifest(manifest) if manifest else failed_entry(
                rec.get("error") or "preprocessing produced no manifest", stage="prep",
                choices={"strategy": plan.value("strategy")}))
            doc["prep"].setdefault(candidate_name, {})[context] = prep_entry

            row = rows[(rows["candidate"] == candidate) & (rows["context"] == rec.get("context"))]
            row_error = None if row.empty else row.iloc[0].get("error")
            if row_error is not None and str(row_error) != "nan":
                entry = failed_entry(row_error, stage="screen",
                                     choices={"strategy": plan.value("strategy")})
            elif rec.get("error"):
                entry = failed_entry(rec["error"], stage="screen",
                                     choices={"strategy": plan.value("strategy")})
            else:
                entry = {
                    "outcome": "ok", "failure": None,
                    "choices": {"strategy": plan.value("strategy")},
                    "dimensions": prep_entry["dimensions"], "steps": [],
                }
            doc["screen"].setdefault(candidate_name, {})[context] = entry

    usable = []
    for candidate, contexts in doc["screen"].items():
        if contexts and all(entry["outcome"] == "ok" for entry in contexts.values()):
            usable.append(candidate)
    doc["usable_candidates"] = sorted(usable)
    doc["selected_candidate"] = _label(selected) if selected is not None else None
    doc["outcome"] = "ok" if selected is not None else "failed"
    doc["status"] = doc["outcome"]
    validate_portable(doc)
    return doc


def _json_safe(value: Any) -> Any:
    """Convert numpy-like scalars and non-finite floats to strict JSON values."""
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, numbers.Integral):
        return int(value)
    if isinstance(value, numbers.Real):
        number = float(value)
        return number if math.isfinite(number) else None
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    item = getattr(value, "item", None)
    if callable(item):
        return _json_safe(item())
    raise TypeError(f"unsupported diagnostics value type: {type(value).__name__}")


def validate_portable(value: Any, *, _where: str = "$", _check_keys: bool = True) -> None:
    """Recursively reject fields or values forbidden in the portable artifact."""
    if _where == "$" and isinstance(value, dict) and "schema" in value:
        if value.get("schema") != SCHEMA or value.get("schema_version") != SCHEMA_VERSION:
            raise ValueError("unsupported diagnostics schema or version")
        for required in ("outcome", "status", "prep", "screen", "models",
                         "usable_candidates", "selected_candidate"):
            if required not in value:
                raise ValueError(f"diagnostics schema missing required field {required!r}")
        if value.get("outcome") not in _OUTCOMES:
            raise ValueError(f"unsupported diagnostics outcome {value.get('outcome')!r}")
    if isinstance(value, dict):
        for key, item in value.items():
            name = str(key).lower()
            dynamic_container = _where in {"$.prep", "$.screen", "$.models"}
            dynamic_key = dynamic_container or not _check_keys
            if (_check_keys and not dynamic_container and
                    (name in _FORBIDDEN_KEYS or name.endswith("_path") or name.endswith("_paths"))):
                raise DiagnosticPrivacyError(f"{_where}.{key}: forbidden field in portable diagnostics")
            validate_portable(str(key), _where=f"{_where}.<key>", _check_keys=True)
            validate_portable(item, _where=f"{_where}.{key}",
                              _check_keys=(True if not dynamic_key else not _check_keys))
        return
    if isinstance(value, (list, tuple)):
        for idx, item in enumerate(value):
            validate_portable(item, _where=f"{_where}[{idx}]")
        return
    if isinstance(value, str):
        path_probe = _RATIO_SLASH.sub(_RATIO_SENTINEL, value)
        if _QUOTED_PATH.search(path_probe) or _PATH_TOKEN.search(path_probe):
            raise DiagnosticPrivacyError(f"{_where}: embedded absolute path in portable diagnostics")
        if _IDENTITY_TOKEN.search(value):
            raise DiagnosticPrivacyError(f"{_where}: identity token in portable diagnostics")
        for identity in _known_identities():
            if _replace_literal(value, identity, "<identity>", boundary=True) != value:
                raise DiagnosticPrivacyError(f"{_where}: local identity in portable diagnostics")


def write_diagnostics(path: str | Path, document: dict) -> Path:
    """Validate and write deterministic UTF-8 JSON."""
    portable = _json_safe(document)
    validate_portable(portable)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(portable, indent=2, sort_keys=True, ensure_ascii=False,
                                 allow_nan=False) + "\n",
                      encoding="utf-8", newline="\n")
    return target


def try_write_diagnostics(path: str | Path, builder: Callable[..., dict],
                          *args, **kwargs) -> Path | None:
    """Observe a run without ever turning a diagnostic failure into a run failure."""
    try:
        return write_diagnostics(path, builder(*args, **kwargs))
    except Exception as exc:  # noqa: BLE001 - diagnostics are deliberately non-fatal
        warnings.warn(
            f"portable diagnostics unavailable: {type(exc).__name__}: {safe_message(exc)}",
            RuntimeWarning,
            stacklevel=2,
        )
        return None
