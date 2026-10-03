"""paper/make_figures.py: digit-free templates, tables equal the committed reports, reruns are stable."""
import importlib.util
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("make_figures", ROOT / "paper" / "make_figures.py")
mf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mf)

pytestmark = pytest.mark.skipif(not (mf.RESULTS / mf.CLEANING_EXP / "rows").exists(), reason="results not checked out")


def test_templates_contain_no_literal_numbers():
    for t in (mf.GEN, mf.MOTIF_CAPTION, mf.SELECTION_CAPTION, mf.CLEANING_CAPTION, *mf.MOTIF_STATUS.values()):
        assert re.search(r"\d", t) is None, t


def _report_rows(exp: str, label: str) -> dict:
    """(A, B) -> (n, diff, interval, p) from the 'auto pipeline' declared-comparison table of REPORT_REAL.md."""
    text = (mf.RESULTS / exp / "REPORT_REAL.md").read_text(encoding="utf-8")
    block = text.split(f"#### Declared comparisons, {label}, auto pipeline")[1].split("####")[0].split("###")[0]
    fam = int(re.search(r"Holm over (\d+) directed", block).group(1))
    out = {}
    for line in block.splitlines():
        cells = [c.strip() for c in line.strip("|").split("|")]
        if len(cells) == 8 and cells[0] not in ("A", "---"):
            out[(cells[0], cells[1])] = (int(cells[2]), cells[3], cells[4], cells[5], fam)
    return out


@pytest.mark.parametrize("metric,label", [("auprc_ratio", "AUPRC ratio"), ("epr", "EPR")])
def test_motif_rows_equal_committed_report(metric, label):
    recs = mf.motif_lift_records(metric)
    assert len(recs) == 6  # C4: 1 pair, C5: 2, C6: 3
    for r in recs:
        rep = _report_rows(r["exp"], label)[(r["a"], r["b"])]
        got = (r["n"], f"{r['diff']:+.3f}", mf._interval(r["ci"]), f"{r['p_a_gt_b']:.3g}", r["family"])
        assert got == rep, (r["exp"], r["a"], r["b"])


def test_cleaning_rows_equal_committed_summary():
    text = (mf.RESULTS / mf.CLEANING_EXP / "SUMMARY.md").read_text(encoding="utf-8")
    sections = text.split("## Variants")
    own = {m[0]: m[1:] for m in re.findall(r"\| (base\+\S+) \| ([+-][\d.]+) \[(\d+)\] \|", sections[1])}
    shared = {m[0]: m[1:] for m in re.findall(r"\| (base\+\S+) \| ([+-][\d.]+) \[(\d+)\] \|", sections[2])}
    assert own and shared
    for cond, name, n_own, d_own, n_sh, d_sh in mf.cleaning_variants_rows():
        assert (d_own, str(n_own)) == own[name] and (d_sh, str(n_sh)) == shared[name], name


def test_selection_table_matches_json_row_count():
    import json
    data = json.loads((mf.RESULTS / mf.SELECTION_EXP / "CONFIRM.json").read_text(encoding="utf-8"))
    body = mf.selection_reference_free().split("\n\n", 1)[1].strip().splitlines()[2:]
    assert len(body) == sum(len(v) for v in data.values())


def test_running_twice_gives_no_diff():
    tables = ROOT / "paper" / "tables"
    run = lambda: subprocess.run([sys.executable, str(ROOT / "paper" / "make_figures.py")], cwd=ROOT,  # noqa: E731
                                 check=True, capture_output=True)
    run()
    first = {p.name: p.read_bytes() for p in tables.glob("*.md")}
    run()
    assert first == {p.name: p.read_bytes() for p in tables.glob("*.md")}
