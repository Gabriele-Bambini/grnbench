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
