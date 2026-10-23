"""Every report format consumes the same validated portable diagnostic."""
import json

import pytest

from grnbench import texreport
from grnbench.cli_eval import _render_markdown, main as eval_main
from grnbench.diagnostics import (
    UNAVAILABLE_NOTE,
    diagnostics_view,
    new_document,
    read_diagnostics,
    write_diagnostics,
)
from grnbench.html_report import render_html


def _diagnostic():
    doc = new_document(outcome="partial")
    doc["selected_candidate"] = "recommended"
    doc["usable_candidates"] = ["recommended"]
    doc["prep"] = {
        "recommended": {
            "all_cells": {
                "outcome": "ok",
                "failure": None,
                "choices": {"strategy": "cpm_log1p", "n_hvg": 2000},
                "dimensions": {
                    "before": {"cells": 120, "genes": 5000},
                    "after": {"cells": 100, "genes": 2000},
                },
                "steps": [{"stage": "qc", "after": {"cells": 100, "genes": 5000}}],
            }
        }
    }
    doc["models"] = {
        "grnboost2": {
            "all_cells": {
                "outcome": "failed",
                "failure": {"category": "unavailable_dependency", "message": "optional package missing"},
                "choices": {"model": "grnboost2"},
                "dimensions": {"before": {"genes": 2000}, "after": None},
                "steps": [],
            }
        }
    }
    return doc


def _evaluation_report():
    return {
        "gold": "fixture", "gold_pairs": 3, "null_mode": "skipped", "n_nulls": 0,
        "methods": {}, "overlap": {},
    }


def test_reader_validates_and_all_formats_render_the_same_rows(tmp_path):
    path = write_diagnostics(tmp_path / "diagnostics.json", _diagnostic())
    document = read_diagnostics(path)
    view = diagnostics_view(document)
    assert [(row["stage"], row["name"], row["context"], row["outcome"])
            for row in view["entries"]] == [
        ("prep", "recommended", "all_cells", "ok"),
        ("models", "grnboost2", "all_cells", "failed"),
    ]

    rendered = (
        _render_markdown(_evaluation_report(), document),
        render_html(_evaluation_report(), diagnostics=document),
        texreport._diagnostics_tex(document),
    )
    for text in rendered:
        text = text.replace("\\_", "_")
        for expected in ("recommended", "grnboost2", "all_cells", "unavailable_dependency",
                         "optional package missing", "120", "5,000", "100", "2,000",
                         "cpm_log1p"):
            assert expected in text, (expected, text)


def test_historical_runs_are_explicitly_marked_unavailable_in_every_format():
    for text in (
        _render_markdown(_evaluation_report(), None),
        render_html(_evaluation_report(), diagnostics=None),
        texreport._diagnostics_tex(None),
    ):
        assert UNAVAILABLE_NOTE in text


@pytest.mark.parametrize("payload,match", [
    ("{", "malformed diagnostics JSON"),
    (json.dumps({"schema": "grnbench.diagnostics", "schema_version": 99}), "unsupported diagnostics schema"),
    (json.dumps({"prep": {}, "screen": {}, "models": {}, "usable_candidates": []}),
     "missing required field 'schema'"),
    (json.dumps([]), "root must be a JSON object"),
])
def test_existing_bad_diagnostics_fail_clearly(tmp_path, payload, match):
    path = tmp_path / "diagnostics.json"
    path.write_text(payload, encoding="utf-8")
    with pytest.raises(ValueError, match=match):
        read_diagnostics(path)


def test_missing_diagnostics_is_historical_not_malformed(tmp_path):
    assert read_diagnostics(tmp_path / "missing.json") is None


def test_eval_turns_missing_schema_into_controlled_report_failure(tmp_path, capsys):
    path = tmp_path / "diagnostics.json"
    path.write_text(json.dumps({
        "prep": {}, "screen": {}, "models": {}, "usable_candidates": [],
    }), encoding="utf-8")
    assert eval_main(["--gold", "dorotheaABC", "--diagnostics", str(path)]) == 2
    error = capsys.readouterr().err
    assert "[REPORT FAIL]" in error
    assert "missing required field 'schema'" in error


@pytest.mark.parametrize("choice,value", [
    ("barcode_list", ["AAACCCAAGAAACCAT-1"]),
    ("expr", [0.0, 1.0, 2.0]),
])
def test_reader_rejects_choice_fields_outside_writer_allowlist(tmp_path, choice, value):
    document = _diagnostic()
    document["prep"]["recommended"]["all_cells"]["choices"][choice] = value
    path = tmp_path / "diagnostics.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    with pytest.raises(ValueError, match="choices has unsupported fields"):
        read_diagnostics(path)


def test_reader_rejects_unsanitised_barcode_and_values_in_failure_message(tmp_path):
    document = _diagnostic()
    document["models"]["grnboost2"]["all_cells"]["failure"]["message"] = (
        "cell AAACCCAAGAAACCAT-1 had values [0.0, 1.0, 2.0]"
    )
    path = tmp_path / "diagnostics.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    with pytest.raises(ValueError, match="must be sanitised"):
        read_diagnostics(path)


def test_markdown_escapes_html_in_failure_messages():
    document = _diagnostic()
    document["models"]["grnboost2"]["all_cells"]["failure"]["message"] = (
        "optional <script>alert(1) missing"
    )
    rendered = _render_markdown(_evaluation_report(), document)
    assert "<script>" not in rendered
    assert "&lt;script&gt;alert(1)" in rendered
