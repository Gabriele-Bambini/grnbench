from __future__ import annotations

import json
import getpass
import socket

import numpy as np
import pytest

from grnbench.diagnostics import (
    DiagnosticPrivacyError,
    MAX_MESSAGE_CHARS,
    context_label,
    failed_prep_document,
    model_document,
    new_document,
    prep_document,
    safe_message,
    try_write_diagnostics,
    validate_portable,
    write_diagnostics,
)


def _prep_manifest():
    return {
        "input": {"path": "/srv/example/counts.h5ad"},
        "qc": {"metrics": {"n_cells": 120, "n_genes": 400}},
        "args": {"strategy": "cpm_log1p", "n_hvg": 100,
                 "subset_key": "cell_type", "subset_value": "T cell",
                 "label_mapping": "/srv/example/mapping.tsv"},
        "steps": [
            {"step": "ingest", "detail": "loaded /srv/example/counts.h5ad"},
            {"step": "filter_cells", "n_cells_after": 110, "n_genes_after": 400,
             "detail": "removed 10 cells"},
            {"step": "feature_selection", "n_cells_after": 110, "n_genes_after": 103,
             "detail": "kept requested features"},
        ],
    }


def test_prep_document_is_portable_and_records_dimensions_choices_and_steps():
    doc = prep_document(_prep_manifest())
    entry = doc["prep"]["default"]["all_cells"]
    assert doc["schema_version"] == 1
    assert entry["choices"] == {
        "strategy": "cpm_log1p", "n_hvg": 100,
        "subset_key": "cell_type", "subset_value": "T cell",
    }
    assert entry["dimensions"] == {
        "before": {"cells": 120, "genes": 400},
        "after": {"cells": 110, "genes": 103},
    }
    assert entry["steps"][0]["message"] == "loaded <path>"
    validate_portable(doc)
    assert "/srv/example" not in json.dumps(doc)


def test_prep_document_distinguishes_requested_and_effective_strategy():
    manifest = _prep_manifest()
    manifest["args"]["strategy"] = "auto"
    manifest["strategy_chosen"] = "pearson_residuals"
    choices = prep_document(manifest)["prep"]["default"]["all_cells"]["choices"]
    assert choices["requested_strategy"] == "auto"
    assert choices["strategy"] == "pearson_residuals"


def test_failed_prep_labels_strategy_as_requested_not_applied():
    choices = failed_prep_document(
        "preprocessing failed", choices={"strategy": "auto"}
    )["prep"]["default"]["all_cells"]["choices"]
    assert choices == {"requested_strategy": "auto"}


def test_context_label_matches_subset_and_explicit_context_contract():
    assert context_label() == "all_cells"
    assert context_label(subset_key="cell_type", subset_value="alpha") == "cell_type=alpha"
    assert context_label(subset_key="cell_type", subset_value="alpha",
                         explicit="alpha cells") == "alpha cells"


def test_failure_messages_are_redacted_bounded_and_categorised():
    message = "QC gate failure reading /srv/example/private.h5ad " + "x" * 500
    doc = failed_prep_document(message)
    failure = doc["prep"]["default"]["all_cells"]["failure"]
    assert failure["category"] == "qc_gate"
    assert "<path>" in failure["message"] and "/srv/example" not in failure["message"]
    assert len(failure["message"]) <= MAX_MESSAGE_CHARS
    assert safe_message("owner@example.org failed") == "<identity> failed"


@pytest.mark.parametrize("message", [
    r"FileNotFoundError: [Errno 2] missing: 'Z:\restricted\nope\absent.h5ad'",
    r'failed at "D:\data set\private\data.h5ad"',
    r"failed(key=/srv/private/data.h5ad)",
    r"failed at \\server\share\data.h5ad",
    r"failed at ~\private\data.h5ad",
    r"failed at file:///restricted/private/data.h5ad",
])
def test_safe_message_redacts_quoted_embedded_unc_tilde_and_uri_paths(message):
    cleaned = safe_message(message)
    assert "data.h5ad" not in cleaned
    assert "<path>" in cleaned


def test_safe_message_redacts_real_exception_and_local_identity():
    error = FileNotFoundError(2, "No such file", r"Z:\restricted\missing.h5ad")
    cleaned = safe_message(
        f"{error}; user={getpass.getuser()} host={socket.gethostname()}")
    assert r"Z:\restricted" not in cleaned
    assert getpass.getuser().lower() not in cleaned.lower()
    assert socket.gethostname().lower() not in cleaned.lower()
    assert cleaned.count("<identity>") >= 2


def test_safe_message_does_not_treat_package_versions_or_ratios_as_identities_or_paths():
    message = "numpy@1.26 sce@assays kept 120 /400 CD4+/CD8+ https://example.org/a"
    assert safe_message(message) == message
    assert safe_message("missing /400/private/file") == "missing <path>"


@pytest.mark.parametrize("bad", [
    {"path": "relative.tsv"},
    {"nested": {"username": "person"}},
    {"message": "/srv/example/file"},
    {"message": r"Z:\\restricted\\file"},
    {"message": "owner@example.org"},
    {"/srv/example/key": "value"},
])
def test_portability_validator_rejects_paths_and_identity_fields(bad):
    with pytest.raises(DiagnosticPrivacyError):
        validate_portable(bad)


def test_writer_is_deterministic_and_refuses_unsafe_content(tmp_path):
    doc = prep_document(_prep_manifest())
    first = write_diagnostics(tmp_path / "a.json", doc).read_bytes()
    second = write_diagnostics(tmp_path / "b.json", doc).read_bytes()
    assert first == second and first.endswith(b"\n")
    doc["prep"]["default"]["all_cells"]["path"] = "/secret"
    with pytest.raises(DiagnosticPrivacyError):
        write_diagnostics(tmp_path / "bad.json", doc)


def test_writer_emits_strict_lf_json_and_normalises_numpy_and_nonfinite(tmp_path):
    doc = prep_document(_prep_manifest())
    choices = doc["prep"]["default"]["all_cells"]["choices"]
    choices["n_hvg"] = np.int64(7)
    choices["subsample_cells"] = np.nan
    raw = write_diagnostics(tmp_path / "strict.json", doc).read_bytes()
    assert b"\r\n" not in raw and b"NaN" not in raw and b"Infinity" not in raw
    parsed = json.loads(raw)
    assert parsed["prep"]["default"]["all_cells"]["choices"]["n_hvg"] == 7
    assert parsed["prep"]["default"]["all_cells"]["choices"]["subsample_cells"] is None


def test_safe_message_redacts_probable_barcodes_and_numeric_vectors():
    assert safe_message(
        "cell AAACCCAAGAAACCAT-1 had expression [0.0, 1.0, 2.0]"
    ) == "cell <barcode> had expression <values>"


def test_nonfinite_dimension_is_unknown_instead_of_crashing():
    manifest = _prep_manifest()
    manifest["qc"]["metrics"]["n_cells"] = np.nan
    manifest["qc"]["metrics"]["n_genes"] = np.int64(400)
    before = prep_document(manifest)["prep"]["default"]["all_cells"]["dimensions"]["before"]
    assert before == {"genes": 400}


def test_diagnostics_are_best_effort_and_never_fail_the_observed_run(tmp_path):
    def broken_builder():
        raise RuntimeError("diagnostic bug")

    with pytest.warns(RuntimeWarning, match="portable diagnostics unavailable"):
        assert try_write_diagnostics(tmp_path / "diagnostics.json", broken_builder) is None
    assert not (tmp_path / "diagnostics.json").exists()


def test_dynamic_candidate_context_and_model_names_are_not_schema_fields():
    doc = prep_document(_prep_manifest(), candidate="output", context="Values")
    validate_portable(doc)
    model = model_document({"status": "complete", "runs": {"matrix": {}},
                            "effective_genes": {"count": 3}})
    validate_portable(model)
    assert model["outcome"] == "ok" and model["status"] == "complete"


def test_model_document_keeps_success_failure_and_unavailability_without_paths():
    manifest = {
        "status": "partial_failure",
        "effective_genes": {"count": 321, "values": ["not copied"]},
        "runs": {"pearson_signed": {"edges_file": "/srv/example/edges.tsv"}},
        "methods_skipped": {"genie3": "unavailable: arboreto missing"},
        "errors": {"broken": "RuntimeError at /srv/example/model.py"},
    }
    doc = model_document(manifest)
    assert doc["models"]["pearson_signed"]["all_cells"]["outcome"] == "ok"
    assert doc["models"]["genie3"]["all_cells"]["failure"]["category"] == "unavailable_dependency"
    assert doc["models"]["broken"]["all_cells"]["failure"]["message"] == "RuntimeError at <path>"
    assert "/srv/example" not in json.dumps(doc)


def test_model_document_accepts_the_real_context_label():
    manifest = {
        "status": "complete",
        "effective_genes": {"count": 12},
        "runs": {"pearson_signed": {}},
        "errors": {"broken": "failed"},
    }
    doc = model_document(manifest, context="cell_type=alpha")
    assert set(doc["models"]["pearson_signed"]) == {"cell_type=alpha"}
    assert set(doc["models"]["broken"]) == {"cell_type=alpha"}


def test_safe_message_collapses_relative_subpaths_after_path_redaction():
    message = "FileNotFoundError: Did not find file rel_mtx_folder/matrix.mtx.gz"
    cleaned = safe_message(message, sensitive_values=("rel_mtx_folder",))
    assert cleaned == "FileNotFoundError: Did not find file <path>"
    doc = failed_prep_document(message, choices={"data": "rel_mtx_folder"},
                             sensitive_values=("rel_mtx_folder",))
    validate_portable(doc)


def test_safe_message_redacts_unquoted_sensitive_windows_path_with_spaces():
    path = r"D:\data\First Last\My Docs\tfs.txt"
    assert safe_message(
        "FileNotFoundError: " + path, sensitive_values=(path,)
    ) == "FileNotFoundError: <path>"
    assert safe_message(
        f"FileNotFoundError: '{path}'", sensitive_values=(path,)
    ) == "FileNotFoundError: <path>"


def test_no_edges_run_status_is_a_failure():
    assert new_document(outcome="no edges")["outcome"] == "failed"


def test_prep_document_uses_pre_calling_n_barcodes_for_dimensions_before():
    manifest = _prep_manifest()
    manifest["qc"]["metrics"]["n_cells"] = 1438
    manifest["steps"].insert(0, {
        "step": "cell_calling",
        "n_cells_after": 1438,
        "n_genes_after": 400,
        "n_barcodes": 31500,
        "n_called": 1438,
        "detail": "knee: kept 1438 of 31500 barcodes",
    })
    doc = prep_document(manifest)
    dims = doc["prep"]["default"]["all_cells"]["dimensions"]
    assert dims["before"] == {"cells": 31500, "genes": 400}
    assert dims["after"] == {"cells": 110, "genes": 103}


# --------------------------------------------------------------------------- #
# safe_message is idempotent: every writer output passes the reader (#177, B1)
# --------------------------------------------------------------------------- #

_TOKEN_PIECES = ("ACGT" * 5, "AAACCCAAGAAACCAT-1", "/srv/data/x.h5ad", "C:\\Users\\someone\\f.txt",
                 "a.b@example.org", "[0.1, 2, 3.5]", "~/x", "file:///tmp/a", " 3 / 4 ",
                 "\\\\host\\share\\f", "'/q/p'", "ACGTN" * 4)
_TOKEN_CHARS = "ACGTNacgtxyzXYZ0123456789-/\\:.@~_ ,[]()'\"\t"


def _random_texts(n: int, seed: int) -> list[str]:
    """Long strings dense in the pieces safe_message rewrites, so the cut often lands inside one."""
    import random

    rng = random.Random(seed)
    out = []
    for _ in range(n):
        parts: list[str] = []
        target = rng.randint(1, 420)
        while sum(map(len, parts)) < target:
            if rng.random() < 0.4:
                parts.append(rng.choice(_TOKEN_PIECES))
            else:
                parts.append("".join(rng.choice(_TOKEN_CHARS) for _ in range(rng.randint(1, 30))))
        out.append("".join(parts))
    return out


@pytest.mark.parametrize("text,limit", [
    ("X" * 60 + "-" + "ACGT" * 5 + "ZZZZZ", 80),                  # a gene name cut inside an A/C/G/T run
    ("x" * 280 + "-" + "ACGT" * 5 + "zzzz", MAX_MESSAGE_CHARS),   # a step message, same at the 300 bound
])
def test_a_cut_inside_a_token_gives_a_stable_message(text, limit):
    once = safe_message(text, limit=limit)
    assert len(once) <= limit
    assert safe_message(once, limit=limit) == once


@pytest.mark.parametrize("limit", [20, 80, 120, MAX_MESSAGE_CHARS])
def test_safe_message_is_idempotent_on_random_long_strings(limit):
    for text in _random_texts(1500, seed=limit):
        once = safe_message(text, limit=limit)
        assert len(once) <= limit, (text, once)
        assert safe_message(once, limit=limit) == once, (text, once)


def test_every_writer_output_passes_the_reader(tmp_path):
    from grnbench.diagnostics import read_diagnostics

    for i, text in enumerate(_random_texts(200, seed=7)):
        row = {"step": text, "cells_before": 5, "cells_after": 5, "cells_removed": 0,
               "genes_before": 3, "genes_after": 2, "genes_removed": 1, "aggregated": False,
               "tfs_before": 1, "tfs_after": 0, "tfs_removed": 1, "tfs_removed_detected": 1,
               "tfs_removed_names": [text], "tfs_removed_detected_names": [text]}
        manifest = {"qc": {"metrics": {"n_cells": 5, "n_genes": 3}},
                    "args": {"strategy": "cpm_log1p", "subset_value": text}, "strategy_chosen": text,
                    "steps": [{"step": "filter_genes", "n_cells_after": 5, "n_genes_after": 2,
                               "detail": text}],
                    "retention": [row]}
        ok = write_diagnostics(tmp_path / f"ok{i}.json", prep_document(manifest, candidate=text, context=text))
        assert read_diagnostics(ok) is not None
        failed = write_diagnostics(tmp_path / f"failed{i}.json", failed_prep_document(
            text, choices={"strategy": text}, candidate=text, context=text,
            sensitive_values=(text[:12],)))
        assert read_diagnostics(failed) is not None

def _ratio_chain(n: int) -> str:
    return "x" + "".join(f" /{k}" for k in range(1, n + 1))


def test_a_ratio_chain_settles_one_link_per_pass_and_a_long_one_is_withheld(tmp_path):
    """Each ' /N' is a ratio only while its left neighbour is 'digit space', so one pass undoes one
    link: 7 links settle inside the 8-pass cap, 8 links do not and are withheld (#177 re-review)."""
    from grnbench.diagnostics import read_diagnostics

    seven = safe_message(_ratio_chain(7))
    assert seven == "x" + " <path>" * 7
    assert safe_message(seven) == seven
    eight = safe_message(_ratio_chain(8))
    assert eight == "<message withheld>"
    assert safe_message(eight) == eight
    path = write_diagnostics(tmp_path / "d.json", failed_prep_document(_ratio_chain(8)))
    entry = read_diagnostics(path)["prep"]["default"]["all_cells"]
    assert entry["failure"]["message"] == "<message withheld>"


_PLACEHOLDER_WORDS = ("path", "identity", "barcode", "values", "message", "withheld")


@pytest.mark.parametrize("word", _PLACEHOLDER_WORDS)
@pytest.mark.parametrize("role", ["sensitive value", "local identity"])
def test_a_placeholder_word_as_a_sensitive_value_or_identity_never_nests(tmp_path, monkeypatch, word, role):
    """`--out path` (or a host called 'identity') used to turn '<path>' into '<<path>>' on every pass
    until the message was withheld (#177 re-review)."""
    from grnbench import diagnostics
    from grnbench.diagnostics import read_diagnostics

    sensitive = (word,) if role == "sensitive value" else ()
    if role == "local identity":
        monkeypatch.setattr(diagnostics, "_known_identities", lambda: (word,))
    message = ("wrote /x/y.h5ad for a.b@example.org; cell AAACCCAAGAAACCAT-1 had [1, 2, 3]; "
               f"kept 3 / 4; the {word} was {word}")
    for text in ("/x", message, "<path> <identity> <barcode> <values> <message withheld>"):
        out = safe_message(text, sensitive_values=sensitive)
        assert "<<" not in out and ">>" not in out, out
        assert safe_message(out) == out and safe_message(out, sensitive_values=sensitive) == out
    assert safe_message("/x", sensitive_values=sensitive) == "<path>"
    assert safe_message(message, sensitive_values=sensitive) != "<message withheld>"
    path = write_diagnostics(tmp_path / "d.json", failed_prep_document(message, sensitive_values=sensitive))
    stored = read_diagnostics(path)["prep"]["default"]["all_cells"]["failure"]["message"]
    assert "<<" not in stored and stored.startswith("wrote <path> for <identity>")


@pytest.mark.parametrize("host", ["message", "category", "outcome"])
def test_a_host_named_like_a_schema_field_can_still_write_a_failure_document(tmp_path, monkeypatch, host):
    from grnbench import diagnostics
    from grnbench.diagnostics import read_diagnostics

    monkeypatch.setattr(diagnostics, "_known_identities", lambda: (host,))
    path = write_diagnostics(tmp_path / "d.json", failed_prep_document("boom"))
    assert read_diagnostics(path)["prep"]["default"]["all_cells"]["failure"]["message"] == "boom"


def test_only_schema_field_names_are_exempt_from_the_local_identity_check(monkeypatch):
    """The exemption is an explicit list of schema field names, not 'any key at a fixed position':
    a key derived from data (a per-batch or per-TF dict, a future field) is still compared."""
    from grnbench import diagnostics

    monkeypatch.setattr(diagnostics, "_known_identities", lambda: ("node-17",))
    doc = failed_prep_document("boom")
    validate_portable(doc)

    def refused(change):
        broken = json.loads(json.dumps(doc))
        change(broken)
        with pytest.raises(DiagnosticPrivacyError, match="local identity"):
            validate_portable(broken)

    refused(lambda d: d["prep"]["default"]["all_cells"].update({"node-17": 1}))              # entry key
    refused(lambda d: d["prep"]["default"]["all_cells"].update({"per_batch": {"node-17": 3}}))  # data dict
    refused(lambda d: d["prep"].update({"node-17": d["prep"].pop("default")}))                # candidate

def _fixed_position_keys(value, where="$", check=True) -> set:
    """Keys at fixed (non-data) positions, walked the way validate_portable walks them."""
    keys: set = set()
    if isinstance(value, dict):
        for key, item in value.items():
            dynamic = where in {"$.prep", "$.screen", "$.models"} or not check
            if not dynamic:
                keys.add(key)
            keys |= _fixed_position_keys(item, f"{where}.{key}", True if not dynamic else not check)
    elif isinstance(value, list):
        for item in value:
            keys |= _fixed_position_keys(item, f"{where}[]")
    return keys


def test_every_field_name_the_writers_emit_is_in_the_schema_field_list():
    """The identity exemption list must cover what the writers write (and nothing is exempt by
    position alone): a writer that adds a field without listing it fails here, not on a user's host."""
    import dataclasses

    from grnbench.diagnostics import _MANIFEST_CHOICE_KEYS, _SCHEMA_FIELD_NAMES
    from grnbench.qc import Thresholds

    args = {key: "x" for key in _MANIFEST_CHOICE_KEYS}
    args["thresholds"] = dataclasses.asdict(Thresholds())
    retention = [{"step": "filter_genes", "cells_before": 5, "cells_after": 5, "cells_removed": 0,
                  "genes_before": 3, "genes_after": 2, "genes_removed": 1, "aggregated": False,
                  "tfs_before": 1, "tfs_after": 0, "tfs_removed": 1, "tfs_removed_detected": 1,
                  "tfs_removed_names": ["A"], "tfs_removed_detected_names": ["A"]}]
    manifest = {"qc": {"metrics": {"n_cells": 5, "n_genes": 3}}, "args": args, "strategy_chosen": "x",
                "steps": [{"step": "filter_genes", "n_cells_after": 5, "n_genes_after": 2, "detail": "d"}],
                "retention": retention}
    docs = [prep_document(manifest), failed_prep_document("boom", choices=args),
            model_document({"status": "partial", "runs": {"m1": {}}, "errors": {"m2": "boom"},
                            "methods_skipped": {"m3": "missing"}, "effective_genes": {"count": 3}})]
    for doc in docs:
        assert _fixed_position_keys(doc) <= _SCHEMA_FIELD_NAMES, _fixed_position_keys(doc) - _SCHEMA_FIELD_NAMES

@pytest.mark.parametrize("value,text,expected", [
    ("a<path>b", "output folder a<path>b already exists", "output folder <path> already exists"),
    ("a<path>b", "copy a<path>b/srv/data/x.h5ad now", "copy <path> now"),
    ("my <identity> dir", "open my <identity> dir/f failed", "open <path> failed"),
])
def test_a_sensitive_value_that_contains_a_placeholder_is_still_redacted(tmp_path, value, text, expected):
    """A Linux path may contain '<' and '>': a sensitive value spelled like 'a<path>b' must be found
    in the text even though it contains a placeholder (671c2c8 re-review)."""
    from grnbench.diagnostics import read_diagnostics

    out = safe_message(text, sensitive_values=(value,))
    assert out == expected
    assert safe_message(out) == out and safe_message(out, sensitive_values=(value,)) == out
    path = write_diagnostics(tmp_path / "d.json", failed_prep_document(text, sensitive_values=(value,)))
    assert read_diagnostics(path)["prep"]["default"]["all_cells"]["failure"]["message"] == expected


@pytest.mark.parametrize("value", ["<path", "path>", "<message withheld", "<identity>", "a<path>b"])
def test_a_sensitive_value_with_markup_never_matches_inside_a_placeholder(value):
    out = safe_message(f"see {value} and /x", sensitive_values=(value,))
    assert out == "see <path> and <path>"
    assert safe_message(out, sensitive_values=(value,)) == out