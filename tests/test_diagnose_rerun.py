import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
import diagnose_rerun  # noqa: E402


def _row(pm, ps, candidate="baseline", seconds=1.0, extra=None):
    row = {"m_pearson_motif_auto_auprc_ratio": pm, "m_pearson_signed_auto_auprc_ratio": ps,
           "auto_candidate": candidate, "seconds": seconds, "prep_seconds": seconds,
           "auto_eval_mean": float("nan"), "n_cells": 100}
    row.update(extra or {})
    return row


def _write(folder: Path, rows: dict) -> Path:
    folder.mkdir(parents=True)
    for name, row in rows.items():
        (folder / f"{name}.json").write_text(json.dumps(row), encoding="utf-8")
    return folder


def test_identical_rows_report_no_difference(tmp_path):
    rows = {"a": _row(1.5, 1.3), "b": _row(1.2, 1.1)}
    out = diagnose_rerun.diagnose(_write(tmp_path / "c", rows), _write(tmp_path / "r", rows))
    assert out[0] == "datasets with at least one differing field: 0/2"
    assert out[1] == "numeric fields differing: 0"
    assert "2 models x 2 datasets): 0 differ" in out[-1]


def test_timing_and_nan_are_not_differences(tmp_path):
    c = _write(tmp_path / "c", {"a": _row(1.5, 1.3, seconds=1.0)})
    r = _write(tmp_path / "r", {"a": _row(1.5, 1.3, seconds=9.0)})
    assert diagnose_rerun.diagnose(c, r)[0] == "datasets with at least one differing field: 0/1"


def test_numeric_drift_planner_flip_and_primary_inputs(tmp_path):
    c = _write(tmp_path / "c", {"a": _row(2.0, 1.0), "b": _row(1.0, 1.0)})
    r = _write(tmp_path / "r", {"a": _row(2.2, 1.0, candidate="context_scope=pooled"),
                                "b": _row(1.0, 1.0, extra={"n_cells": 101})})
    out = diagnose_rerun.diagnose(c, r)
    assert out[0] == "datasets with at least one differing field: 2/2"
    assert out[1].startswith("numeric fields differing: 2;")
    assert "max 1.00e-01" in out[1]
    assert "  a  auto_candidate: baseline -> context_scope=pooled" in out
    assert "2 models x 2 datasets): 1 differ; relative difference min 1.00e-01, max 1.00e-01" in out[-2]
    assert out[-1].startswith("  1.00e-01  pearson_motif")


@pytest.mark.parametrize("bad", [float("nan"), None, "missing"])
def test_primary_nan_null_or_missing_on_one_side_is_reported(tmp_path, bad):
    rerun_row = _row(1.5, 1.3)
    if bad == "missing":
        del rerun_row["m_pearson_motif_auto_auprc_ratio"]
    else:
        rerun_row["m_pearson_motif_auto_auprc_ratio"] = bad
    c = _write(tmp_path / "c", {"a": _row(1.5, 1.3)})
    r = _write(tmp_path / "r", {"a": rerun_row})
    out = diagnose_rerun.diagnose(c, r)
    head = next(line for line in out if line.startswith("primary inputs"))
    assert "1 differ" in head and "1 not comparable" in head
    assert any(line.startswith("  not comparable  pearson_motif") for line in out)
    assert out[1] == "numeric fields differing: 0"


def test_nan_against_number_in_a_secondary_field_is_not_a_relative_difference(tmp_path):
    c = _write(tmp_path / "c", {"a": _row(1.5, 1.3, extra={"auto_eval_mean": 1.2})})
    r = _write(tmp_path / "r", {"a": _row(1.5, 1.3)})
    out = diagnose_rerun.diagnose(c, r)
    assert out[1] == "numeric fields differing: 0"
    assert "  a  auto_eval_mean: 1.2 -> nan" in out


def test_missing_rerun_dataset_exits_2(tmp_path):
    c = _write(tmp_path / "c", {"a": _row(1.5, 1.3), "b": _row(1.2, 1.1)})
    r = _write(tmp_path / "r", {"a": _row(1.5, 1.3)})
    with pytest.raises(SystemExit) as exc:
        diagnose_rerun.diagnose(c, r)
    assert exc.value.code == 2


def test_cli_prints_the_lines(tmp_path, capsys):
    c = _write(tmp_path / "c", {"a": _row(1.5, 1.3)})
    assert diagnose_rerun.main([str(c), str(c)]) == 0
    assert "0/1" in capsys.readouterr().out
