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
import html
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
_MANIFEST_CHOICE_KEYS = (
    "strategy", "n_hvg", "group_key", "label_key", "stratify_key", "batch_key",
    "min_group_cells", "downcast", "subsample_cells", "subsample_mode", "orientation",
    "drop_doublets", "force", "seed", "subset_key", "subset_value", "cell_calling",
    "exclude_qc_genes", "cells_per_metacell", "thresholds",
)
_PORTABLE_CHOICE_KEYS = frozenset((*_MANIFEST_CHOICE_KEYS, "requested_strategy", "model"))
_THRESHOLD_KEYS = frozenset({
    "min_genes_per_cell", "min_cells_per_gene", "mito_warn", "mito_fail",
    "low_genes_cell_frac_warn", "low_genes_cell_frac_fail", "library_cv_warn",
    "sparsity_warn", "sparsity_fail", "min_cells_warn", "min_genes_warn",
    "zero_cells_warn", "zero_cells_fail", "imbalance_ratio_warn",
    "min_class_cells_warn", "doublet_pct_warn", "doublet_rate_max",
    "max_genes_per_cell", "max_counts_per_cell", "max_pct_mito_cell",
})
_BARCODE_TOKEN = re.compile(r"(?<![A-Za-z0-9])[ACGTN]{16,}(?:-\d+)?(?![A-Za-z0-9])", re.IGNORECASE)
_NUMERIC_VECTOR = re.compile(
    r"(?<!\w)[\[(]\s*[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
    r"(?:\s*,\s*[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?){2,}\s*[\])]"
)
UNAVAILABLE_NOTE = (
    "Portable diagnostics unavailable: no valid diagnostics.json was supplied "
    "(the run may predate schema v1 or diagnostics generation may have failed)."
)
# Retention ledger rows (``grnbench.retention``): an optional, additive list per prep entry.
# Older documents simply have no ``retention`` field.
_RETENTION_COUNTS = ("cells_before", "cells_after", "genes_before", "genes_after", "genes_removed",
                     "tfs_before", "tfs_after", "tfs_removed", "tfs_removed_detected")
_RETENTION_NAMES = ("tfs_removed_names", "tfs_removed_detected_names")
_RETENTION_FIELDS = frozenset(("stage", "cells_removed", "aggregated",
                               *_RETENTION_COUNTS, *_RETENTION_NAMES))
_GENE_NAME_CHARS = 80
# The schema's own field names, one list per object; the reader uses the same lists.
_REQUIRED_TOP_FIELDS = ("outcome", "status", "prep", "screen", "models",
                        "usable_candidates", "selected_candidate")
_ENVELOPE_FIELDS = ("schema", "schema_version")
_ENTRY_FIELDS = ("outcome", "failure", "choices", "dimensions", "steps", "retention")
_FAILURE_FIELDS = ("category", "message")
_DIMENSION_SIDES = ("before", "after")
_DIMENSION_FIELDS = frozenset(("cells", "genes"))
_STEP_FIELDS = ("stage", "after", "message")
# Field names are code, not data: only these are exempt from the local user/host comparison of
# validate_portable. Any other key, at whatever position, is compared like a value.
_SCHEMA_FIELD_NAMES = frozenset((
    *_ENVELOPE_FIELDS, *_REQUIRED_TOP_FIELDS, *_ENTRY_FIELDS, *_FAILURE_FIELDS,
    *_DIMENSION_SIDES, *_DIMENSION_FIELDS, *_STEP_FIELDS,
    *_PORTABLE_CHOICE_KEYS, *_THRESHOLD_KEYS, *_RETENTION_FIELDS))
RETENTION_NOTE = (
    "What each cleaning step removed: cells, genes and, when a TF list was given, the TFs of that "
    "list. A removed TF counts as detected when it was non-zero in at least one cell of the "
    "context before the step; feature selection keeps every TF in the data (D-003), the gene "
    "filters before it do not.")
RETENTION_COLUMNS = ("step", "cells", "genes", "TFs", "detected TFs removed", "removed TFs")


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


_MAX_REDACTION_PASSES = 8
_WITHHELD = "<message withheld>"
# The placeholders redaction writes. A sensitive value or a local identity that is one of their words
# ('path', 'identity', 'barcode', 'values', 'message', 'withheld') must not be found inside them:
# '<path>' would become '<<path>>' on every pass until the message is withheld.
_PLACEHOLDERS = re.compile(r"(<(?:path|identity|barcode|values|message withheld)>)")


# Each placeholder as one private-use character (not a regex word character), for _replace_value.
_PLACEHOLDER_MASK = {"<path>": "", "<identity>": "", "<barcode>": "",
                     "<values>": "", "<message withheld>": ""}
_PLACEHOLDER_UNMASK = {mask: placeholder for placeholder, mask in _PLACEHOLDER_MASK.items()}


def _mask_placeholders(text: str) -> str:
    return _PLACEHOLDERS.sub(lambda match: _PLACEHOLDER_MASK[match.group(1)], text)


def _replace_value(text: str, value: Any, replacement: str) -> str:
    """``_replace_literal`` (word boundaries) outside the placeholders already in ``text``."""
    raw = str(value or "").strip()
    if "<" not in raw and ">" not in raw:
        parts = _PLACEHOLDERS.split(text)      # odd indices are placeholders, kept as they are
        return "".join(part if i % 2 else _replace_literal(part, raw, replacement, boundary=True)
                       for i, part in enumerate(parts))
    # A value with markup (a Linux path may hold '<' and '>') can itself contain a placeholder,
    # e.g. 'a<path>b'. Match it against the whole text with every placeholder masked to one
    # character, in the text and in the value alike: it is found where it occurs, never inside a
    # placeholder ('<path' or 'path>' cannot grow '<path>' into '<path>>').
    masked = _mask_placeholders(text)
    for candidate in dict.fromkeys((_mask_placeholders(raw), raw)):
        masked = _replace_literal(masked, candidate, replacement, boundary=True)
    return re.sub("[-]", lambda match: _PLACEHOLDER_UNMASK[match.group(0)], masked)


def _redact_once(text: str, limit: int, path_values: list[str], identities: tuple[str, ...]) -> str:
    """One redaction pass and the length cut (see :func:`safe_message`)."""
    text = _RATIO_SLASH.sub(_RATIO_SENTINEL, text)
    # Replace exact run-specific values before generic path regexes.  Otherwise
    # an unquoted Windows path containing spaces is shortened at its first space
    # and the remaining private suffix can no longer match the original value.
    for value in path_values:
        text = _replace_value(text, value, "<path>")
    text = re.sub(r"(['\"])<path>\1", "<path>", text)
    # A sensitive relative directory can be followed by a filename in an
    # exception (``<path>/matrix.mtx.gz``); keep the whole location redacted
    # before the generic matcher sees the leading slash as a second path.
    text = re.sub(r"<path>(?:[\\/][^\s'\"<>]+)+", "<path>", text)
    text = _QUOTED_PATH.sub("<path>", text)
    text = _PATH_TOKEN.sub("<path>", text)
    text = re.sub(r"(?:<path>){2,}", "<path>", text)
    for value in identities:
        text = _replace_value(text, value, "<identity>")
    text = _IDENTITY_TOKEN.sub("<identity>", text)
    text = _BARCODE_TOKEN.sub("<barcode>", text)
    text = _NUMERIC_VECTOR.sub("<values>", text)
    text = text.replace(_RATIO_SENTINEL, "/")
    text = " ".join(text.split())
    if len(text) > limit:
        text = text[: max(0, limit - 1)].rstrip() + "…"
    return text


def safe_message(message: Any, limit: int = MAX_MESSAGE_CHARS,
                 *, sensitive_values: tuple[Any, ...] | list[Any] = ()) -> str:
    """Collapse, redact and bound an exception or diagnostic message.

    Idempotent: ``safe_message(safe_message(x, n), n) == safe_message(x, n)``, so every value a
    writer stores passes the readers' round-trip check (``safe_message(s) == s``). One pass is not:
    the length cut can end inside a token (an A/C/G/T run followed by the '…' then reads as a
    barcode), and one replacement can expose another (a path right after a redacted identity). The
    pass therefore repeats until the text no longer changes, at most ``_MAX_REDACTION_PASSES``
    times; a text still changing then is withheld (``<message withheld>``, which passes the reader)
    rather than stored unsettled. That cap is reachable: a ' /N' is a ratio, not a path, only while
    its left neighbour is 'digit space', so in 'x /1 /2 /3' each pass undoes one link and a chain
    of n links needs n + 1 passes; 7 links settle, 8 or more are withheld.
    """
    path_values = [str(value).strip() for value in sensitive_values
                   if value is not None and str(value).strip()]
    try:
        home = str(Path.home()).strip()
        if home:
            path_values.append(home)
    except Exception:  # noqa: BLE001 - redaction must not become a run failure
        pass
    path_values = sorted(set(path_values), key=len, reverse=True)
    identities = _known_identities()
    text = str(message)
    for _ in range(_MAX_REDACTION_PASSES):
        again = _redact_once(text, limit, path_values, identities)
        if again == text:
            return text
        text = again
    return _redact_once(_WITHHELD, limit, [], ())


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
    out = {}
    for key in _MANIFEST_CHOICE_KEYS:
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
    entry = {
        "outcome": "ok",
        "failure": None,
        "choices": choices,
        "dimensions": {"before": before, "after": after},
        "steps": steps,
    }
    if isinstance(manifest.get("retention"), list):
        entry["retention"] = _portable_retention(manifest["retention"])
    return entry


def _portable_retention(rows: list) -> list[dict]:
    """Ledger rows of a prep manifest as portable rows: counts, flags and sanitised gene names."""
    from grnbench.retention import TF_NAMES_CAP

    out = []
    for raw in rows:
        if not isinstance(raw, Mapping):
            continue
        item: dict = {"stage": _label(raw.get("step", "unknown"))}
        for key in (*_RETENTION_COUNTS, "cells_removed"):
            if key in raw:
                item[key] = None if raw[key] is None else int(raw[key])
        if "aggregated" in raw:
            item["aggregated"] = bool(raw["aggregated"])
        for key in _RETENTION_NAMES:
            if key in raw:
                item[key] = [safe_message(name, limit=_GENE_NAME_CHARS)
                             for name in list(raw[key] or [])[:TF_NAMES_CAP]]
        out.append(item)
    return out


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


def validate_portable(value: Any, *, _where: str = "$", _check_keys: bool = True,
                      _field_name: bool = False) -> None:
    """Recursively reject fields or values forbidden in the portable artifact.

    The schema's own field names (``_SCHEMA_FIELD_NAMES``: ``message``, ``stage``, ...) are code,
    not data: as keys at a fixed position they are not compared with the local user and host names,
    so a machine called 'message' can still write a document. Every other key (candidate, context
    and model names, or any key at a fixed position that is not a schema field) is compared like a
    value.
    """
    if _where == "$" and isinstance(value, dict) and "schema" in value:
        if value.get("schema") != SCHEMA or value.get("schema_version") != SCHEMA_VERSION:
            raise ValueError("unsupported diagnostics schema or version")
        for required in _REQUIRED_TOP_FIELDS:
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
            validate_portable(str(key), _where=f"{_where}.<key>", _check_keys=True,
                              _field_name=not dynamic_key and str(key) in _SCHEMA_FIELD_NAMES)
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
        for identity in () if _field_name else _known_identities():
            if _replace_value(value, identity, "<identity>") != value:   # not inside a placeholder
                raise DiagnosticPrivacyError(f"{_where}: local identity in portable diagnostics")


def _validate_document_schema(document: Any) -> None:
    """Require the complete schema envelope before a document is rendered."""
    if not isinstance(document, Mapping):
        raise ValueError("diagnostics root must be a JSON object")
    for field in _ENVELOPE_FIELDS:
        if field not in document:
            raise ValueError(f"diagnostics schema missing required field {field!r}")
    if document.get("schema") != SCHEMA or document.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("unsupported diagnostics schema or version")
    for field in _REQUIRED_TOP_FIELDS:
        if field not in document:
            raise ValueError(f"diagnostics schema missing required field {field!r}")
    if document.get("outcome") not in _OUTCOMES:
        raise ValueError(f"unsupported diagnostics outcome {document.get('outcome')!r}")
    status = document.get("status")
    if not isinstance(status, str):
        raise ValueError("diagnostics status must be a string")
    if len(status) > MAX_MESSAGE_CHARS or safe_message(status) != status:
        raise DiagnosticPrivacyError("$.status must be a bounded, sanitised message")


def write_diagnostics(path: str | Path, document: dict) -> Path:
    """Validate and write deterministic UTF-8 JSON."""
    portable = _json_safe(document)
    _validate_document_schema(portable)
    validate_portable(portable)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(portable, indent=2, sort_keys=True, ensure_ascii=False,
                                 allow_nan=False) + "\n",
                      encoding="utf-8", newline="\n")
    return target


def diagnostics_file(path: str | Path) -> Path:
    """Resolve a file-or-run-directory argument to the portable artifact only."""
    source = Path(path)
    return source / "diagnostics.json" if source.is_dir() else source


def read_diagnostics(path: str | Path | None) -> dict | None:
    """Read and validate one schema-v1 portable diagnostic.

    ``None`` and a missing file mean that the run predates portable diagnostics;
    callers render :data:`UNAVAILABLE_NOTE`.  Existing malformed, unsupported or
    privacy-unsafe files raise a clear ``ValueError`` instead of being silently
    treated as historical runs.  A directory names its ``diagnostics.json``.
    """
    if path is None:
        return None
    source = diagnostics_file(path)
    if not source.exists():
        return None
    try:
        document = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"{source}: malformed diagnostics JSON: {exc}") from exc
    if not isinstance(document, dict):
        raise ValueError(f"{source}: diagnostics root must be a JSON object")
    try:
        _validate_document_schema(document)
        diagnostics_view(document)
    except (ValueError, DiagnosticPrivacyError, TypeError) as exc:
        raise ValueError(f"{source}: invalid diagnostics: {exc}") from exc
    return document


def _entry_mapping(value: Any, where: str) -> Mapping:
    if not isinstance(value, Mapping):
        raise ValueError(f"{where} must be an object")
    return value


def _entry_dimensions(value: Any, where: str) -> tuple[dict | None, dict | None]:
    dims = _entry_mapping(value, where)
    out = []
    for side in _DIMENSION_SIDES:
        item = dims.get(side)
        if item is not None:
            item = dict(_entry_mapping(item, f"{where}.{side}"))
            unknown = set(item) - _DIMENSION_FIELDS
            if unknown:
                raise ValueError(f"{where}.{side} has unsupported fields: {sorted(unknown)}")
            for key, number in item.items():
                if isinstance(number, bool) or not isinstance(number, numbers.Integral) or number < 0:
                    raise ValueError(f"{where}.{side}.{key} must be a non-negative integer")
        out.append(item)
    return out[0], out[1]


def _entry_retention(value: Any, where: str) -> list[dict]:
    """Validate the optional ledger rows of one entry as strictly as its dimensions."""
    from grnbench.retention import TF_NAMES_CAP

    if value is None:
        return []
    if not isinstance(value, list):
        raise ValueError(f"{where} must be an array")
    rows = []
    for idx, raw in enumerate(value):
        here = f"{where}[{idx}]"
        row = dict(_entry_mapping(raw, here))
        unknown = set(row) - _RETENTION_FIELDS
        if unknown:
            raise ValueError(f"{here} has unsupported fields: {sorted(unknown)}")
        stage = row.get("stage")
        if not isinstance(stage, str) or len(stage) > MAX_MESSAGE_CHARS or safe_message(stage) != stage:
            raise DiagnosticPrivacyError(f"{here}.stage must be a bounded, sanitised message")
        for key in (*_RETENTION_COUNTS, "cells_removed"):
            if key not in row or (key == "cells_removed" and row[key] is None):
                continue
            number = row[key]
            if isinstance(number, bool) or not isinstance(number, numbers.Integral) or number < 0:
                raise ValueError(f"{here}.{key} must be a non-negative integer")
        if "aggregated" in row and not isinstance(row["aggregated"], bool):
            raise ValueError(f"{here}.aggregated must be a boolean")
        for key in _RETENTION_NAMES:
            if key not in row:
                continue
            names = row[key]
            if not isinstance(names, list) or len(names) > TF_NAMES_CAP:
                raise ValueError(f"{here}.{key} must list at most {TF_NAMES_CAP} names")
            for name in names:
                if not isinstance(name, str) or safe_message(name, limit=_GENE_NAME_CHARS) != name:
                    raise DiagnosticPrivacyError(f"{here}.{key} must hold sanitised gene names")
        rows.append(row)
    return rows


def diagnostics_view(document: dict | None) -> dict:
    """Return the one normalised structure consumed by every report format."""
    if document is None:
        return {
            "available": False,
            "note": UNAVAILABLE_NOTE,
            "outcome": None,
            "status": None,
            "selected_candidate": None,
            "usable_candidates": [],
            "entries": [],
        }
    _validate_document_schema(document)
    validate_portable(document)
    entries = []
    for stage in ("prep", "screen", "models"):
        subjects = _entry_mapping(document.get(stage), f"$.{stage}")
        for subject in sorted(subjects, key=str):
            contexts = _entry_mapping(subjects[subject], f"$.{stage}.{subject}")
            for context in sorted(contexts, key=str):
                where = f"$.{stage}.{subject}.{context}"
                entry = _entry_mapping(contexts[context], where)
                outcome = entry.get("outcome")
                if outcome not in _OUTCOMES:
                    raise ValueError(f"{where}.outcome has unsupported value {outcome!r}")
                failure_value = entry.get("failure")
                if failure_value is not None:
                    failure_value = _entry_mapping(failure_value, f"{where}.failure")
                    if not isinstance(failure_value.get("category"), str) or not isinstance(
                            failure_value.get("message"), str):
                        raise ValueError(f"{where}.failure requires string category and message")
                    if len(failure_value["message"]) > MAX_MESSAGE_CHARS:
                        raise ValueError(f"{where}.failure.message exceeds {MAX_MESSAGE_CHARS} characters")
                    if safe_message(failure_value["message"]) != failure_value["message"]:
                        raise DiagnosticPrivacyError(
                            f"{where}.failure.message must be sanitised before rendering")
                choices = dict(_entry_mapping(entry.get("choices"), f"{where}.choices"))
                unknown_choices = set(choices) - _PORTABLE_CHOICE_KEYS
                if unknown_choices:
                    raise ValueError(
                        f"{where}.choices has unsupported fields: {sorted(unknown_choices)}")
                for choice, choice_value in choices.items():
                    if choice == "thresholds":
                        thresholds = _entry_mapping(choice_value, f"{where}.choices.thresholds")
                        unknown_thresholds = set(thresholds) - _THRESHOLD_KEYS
                        if unknown_thresholds:
                            raise ValueError(
                                f"{where}.choices.thresholds has unsupported fields: "
                                f"{sorted(unknown_thresholds)}")
                        choice_values = thresholds.values()
                    else:
                        choice_values = (choice_value,)
                    if any(isinstance(item, (Mapping, list, tuple)) for item in choice_values):
                        raise ValueError(f"{where}.choices.{choice} must contain scalar values")
                before, after = _entry_dimensions(entry.get("dimensions"), f"{where}.dimensions")
                steps = entry.get("steps")
                if not isinstance(steps, list):
                    raise ValueError(f"{where}.steps must be an array")
                for idx, step in enumerate(steps):
                    step = _entry_mapping(step, f"{where}.steps[{idx}]")
                    message = step.get("message")
                    if message is not None and (
                            not isinstance(message, str)
                            or len(message) > MAX_MESSAGE_CHARS
                            or safe_message(message) != message):
                        raise DiagnosticPrivacyError(
                            f"{where}.steps[{idx}].message must be a bounded, sanitised message")
                retention = _entry_retention(entry.get("retention"), f"{where}.retention")
                entries.append({
                    "stage": stage,
                    "name": str(subject),
                    "context": str(context),
                    "outcome": outcome,
                    "failure_category": failure_value.get("category") if failure_value else None,
                    "message": failure_value.get("message") if failure_value else None,
                    "choices": choices,
                    "before": before,
                    "after": after,
                    "steps": [dict(step) for step in steps],
                    "retention": retention,
                })
    usable = document.get("usable_candidates")
    if not isinstance(usable, list) or not all(isinstance(item, str) for item in usable):
        raise ValueError("$.usable_candidates must be an array of strings")
    selected = document.get("selected_candidate")
    if selected is not None and not isinstance(selected, str):
        raise ValueError("$.selected_candidate must be a string or null")
    return {
        "available": True,
        "note": None,
        "outcome": document["outcome"],
        "status": document["status"],
        "selected_candidate": selected,
        "usable_candidates": list(usable),
        "entries": entries,
    }


def format_dimensions(value: dict | None) -> str:
    """Portable human-readable dimensions shared by Markdown/HTML/LaTeX."""
    if not value:
        return "unknown"
    cells = value.get("cells")
    genes = value.get("genes")
    if cells is not None and genes is not None:
        return f"{cells:,} cells x {genes:,} genes"
    if cells is not None:
        return f"{cells:,} cells"
    if genes is not None:
        return f"{genes:,} genes"
    return "unknown"


def format_choices(value: Mapping) -> str:
    """Stable compact choices text; values remain sourced from diagnostics."""
    if not value:
        return "none"
    return "; ".join(
        f"{key}={json.dumps(item, sort_keys=True, ensure_ascii=False)}"
        for key, item in sorted(value.items())
    )


def _change(before: Any, after: Any, removed: Any, *, aggregated: bool = False) -> str:
    if before is None or after is None:
        return "-"
    text = f"{before:,} -> {after:,}"
    if aggregated:
        return text + " (aggregated)"
    return text + (f" (-{removed:,})" if removed else "")      # a zero removal reads '300 -> 300'


def _names(names: list, total: Any) -> str:
    text = ", ".join(names)
    more = (total or 0) - len(names)
    if more > 0:
        text += f" (+{more:,} more)"
    return text or "-"


def retention_display_rows(rows: list[dict]) -> list[tuple[str, ...]]:
    """One tuple per ledger row, in :data:`RETENTION_COLUMNS` order; every format renders these."""
    out = []
    for row in rows:
        cells = _change(row.get("cells_before"), row.get("cells_after"), row.get("cells_removed"),
                        aggregated=bool(row.get("aggregated")))
        genes = _change(row.get("genes_before"), row.get("genes_after"), row.get("genes_removed"))
        if "tfs_before" in row:
            tfs = _change(row.get("tfs_before"), row.get("tfs_after"), row.get("tfs_removed"))
            detected = f"{row.get('tfs_removed_detected', 0):,}"
            if row.get("tfs_removed_detected_names"):
                detected += ": " + _names(row["tfs_removed_detected_names"], row.get("tfs_removed_detected"))
            removed = _names(row.get("tfs_removed_names") or [], row.get("tfs_removed"))
        else:
            tfs = detected = removed = "no TF list"
        out.append((row["stage"], cells, genes, tfs, detected, removed))
    return out


def retention_sections(view: dict) -> list[tuple[dict, list[tuple[str, ...]]]]:
    """``(entry, display rows)`` for every entry of ``view`` that carries a ledger."""
    return [(entry, retention_display_rows(entry["retention"]))
            for entry in view["entries"] if entry.get("retention")]


def render_diagnostics_markdown(document: dict | None) -> str:
    """Render the shared diagnostic view as a Markdown report section."""
    view = diagnostics_view(document)
    lines = ["## Pipeline diagnostics", ""]
    if not view["available"]:
        return "\n".join(lines + [view["note"], ""])

    def cell(value: Any) -> str:
        return html.escape(str(value), quote=False).replace("|", "\\|").replace("\n", " ")

    def ledger_cell(value: Any) -> str:
        # '>' stays readable in '302 -> 301'; '<' and '&' cannot open markup once escaped
        return (str(value).replace("&", "&amp;").replace("<", "&lt;")
                .replace("|", "\\|").replace("\n", " "))

    lines += [
        f"- Overall outcome: **{cell(view['outcome'])}**",
        f"- Selected candidate: `{cell(view['selected_candidate'] or 'none')}`",
        "- Usable candidates: " + (", ".join(f"`{cell(x)}`" for x in view["usable_candidates"]) or "none"),
        "",
        "| stage | candidate/model | context | outcome | before | after | choices | failure |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for entry in view["entries"]:
        failure_text = "-"
        if entry["failure_category"]:
            failure_text = f"{entry['failure_category']}: {entry['message']}"
        lines.append(
            "| " + " | ".join(cell(value) for value in (
                entry["stage"], entry["name"], entry["context"], entry["outcome"],
                format_dimensions(entry["before"]), format_dimensions(entry["after"]),
                format_choices(entry["choices"]), failure_text,
            )) + " |"
        )
    lines.append("")
    sections = retention_sections(view)
    if sections:
        lines += ["### Retention ledger", "", RETENTION_NOTE, ""]
        for entry, rows in sections:
            lines += [f"{cell(entry['stage'])} `{cell(entry['name'])}`, context `{cell(entry['context'])}`:",
                      "",
                      "| " + " | ".join(RETENTION_COLUMNS) + " |",
                      "|" + "---|" * len(RETENTION_COLUMNS)]
            lines += ["| " + " | ".join(ledger_cell(value) for value in row) + " |" for row in rows]
            lines.append("")
    return "\n".join(lines)


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
