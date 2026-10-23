"""benchmarks/diff_rows_nested.py: where two rows/ folders differ, at any depth, timings apart."""
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("diff_rows_nested", ROOT / "benchmarks" / "diff_rows_nested.py")
dn = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dn)

ROW = {"task": "t1", "seconds": 10.0, "error": None,
       "candidates": [{"name": "baseline", "seconds": 3.5, "eval_mean": 1.25, "choices": {"hvg": "tf_plus_1000"}},
                      {"name": "x", "seconds": 2.0, "eval_mean": float("nan"), "error": "PrepError"}]}


def _write(folder: Path, rows: dict) -> Path:
    (folder / "rows").mkdir(parents=True)
    for name, row in rows.items():
        (folder / "rows" / f"{name}.json").write_text(json.dumps(row), encoding="utf-8")
    return folder / "rows"


def _rerun(**changes):
    row = json.loads(json.dumps(ROW))
    for path, value in changes.items():
        obj = row
        *head, last = path.split(".")
        for k in head:
            obj = obj[int(k)] if k.isdigit() else obj[k]
        obj[int(last) if last.isdigit() else last] = value
    return row


def test_nested_timings_are_counted_apart_and_do_not_differ(tmp_path, capsys):
    a = _write(tmp_path / "a", {"t1": ROW})
    b = _write(tmp_path / "b", {"t1": _rerun(**{"seconds": 11.0, "candidates.0.seconds": 9.9})})
    assert dn.main([str(a), str(b)]) == 0
    out = capsys.readouterr().out
    assert "timing values skipped: 3 (2 differ)" in out
    assert "non-timing differences: 0" in out


def test_nested_score_change_is_found(tmp_path, capsys):
    a = _write(tmp_path / "a", {"t1": ROW})
    b = _write(tmp_path / "b", {"t1": _rerun(**{"candidates.0.eval_mean": 1.2500001})})
    assert dn.main([str(a), str(b)]) == 1
    out = capsys.readouterr().out
    assert "non-timing differences: 1 (numeric 1, other 0)" in out
    assert "t1\tcandidates[0].eval_mean\t1.25 -> 1.2500001" in out
    assert dn.main([str(a), str(b), "--rtol", "1e-6"]) == 0


def test_nan_equals_nan_and_text_change_is_other(tmp_path, capsys):
    a = _write(tmp_path / "a", {"t1": ROW})
    b = _write(tmp_path / "b", {"t1": _rerun(**{"candidates.0.choices.hvg": "tf_plus_2000"})})
    assert dn.main([str(a), str(b)]) == 1
    assert "non-timing differences: 1 (numeric 0, other 1)" in capsys.readouterr().out


def test_missing_row_or_list_length_is_a_difference(tmp_path, capsys):
    a = _write(tmp_path / "a", {"t1": ROW, "t2": ROW})
    short = _rerun()
    short["candidates"] = short["candidates"][:1]
    b = _write(tmp_path / "b", {"t1": short})
    assert dn.main([str(a), str(b)]) == 1
    out = capsys.readouterr().out
    assert "rows: 2 committed, 1 rerun, 1 missing in rerun" in out
    assert "candidates[1]" in out


def test_empty_folder_is_an_error(tmp_path):
    a = _write(tmp_path / "a", {"t1": ROW})
    (tmp_path / "b" / "rows").mkdir(parents=True)
    assert dn.main([str(a), str(tmp_path / "b" / "rows")]) == 2
