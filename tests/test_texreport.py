import json
import re
from pathlib import Path

import numpy as np

from grnbench import texreport

ROOT = Path(__file__).resolve().parents[1]


def test_templates_contain_no_literal_numbers():
    """Every number in a report must come from a result file, so templates hold no digits."""
    for t in texreport.TEMPLATES:
        assert re.search(r"\d", t) is None, re.findall(r".{20}\d.{20}", t)


def _run(tmp_path, n=6):
    rng = np.random.default_rng(0)
    rows = tmp_path / "rows"
    rows.mkdir(parents=True)
    diffs = []
    for i in range(n):
        base = 1.2 + rng.normal(0, .05)
        lift = .2 + rng.normal(0, .03)
        diffs.append(lift)
        row = {"task": f"cxg__d{i}", "error": None}
        for arm in ("auto", "baseline"):
            row[f"m_pearson_signed_{arm}_auprc_ratio"] = base
            row[f"m_pearson_motif_{arm}_auprc_ratio"] = base + lift
            row[f"m_motif_prior_{arm}_auprc_ratio"] = base + lift / 2
            row[f"{arm}_doublets"] = "heuristic" if i == 0 else "scrublet"
        (rows / f"d{i}.json").write_text(json.dumps(row))
    return tmp_path, float(np.mean(diffs))


def test_report_numbers_come_from_the_rows(tmp_path):
    run, mean_diff = _run(tmp_path / "run")
    tex = texreport.build_report(str(run), exp="EXP-TEST", out=tmp_path / "out", compile_pdf=False)
    values = (tmp_path / "out" / "values.tex").read_text()
    key = texreport.macro_name("diff", "auprc_ratio", "auto", "pearson_motif", "vs", "pearson_signed")
    assert "\\newcommand{\\" + key + "}{" + f"{mean_diff:+.3f}" + "}" in values
    body = tex.read_text()
    assert "\\input{tables/models.tex}" in body
    assert "\\input{tables/comparisons_auprc_ratio_auto.tex}" in body
    assert "fell back to the heuristic in 1 of 6 datasets" in body


def test_report_matches_the_committed_real_bench_report():
    """Same statistics path as benchmarks/real_bench.py: EXP-20260925-LM-1 numbers agree."""
    run = ROOT / "benchmarks" / "results" / "engine-2026-09-25-models" / "test12"
    df = texreport.load_rows(run)
    rec = next(r for r in texreport.comparisons(df, "auprc_ratio", "auto")
               if (r["a"], r["b"]) == ("pearson_motif", "pearson_signed"))
    report = (run / "REPORT_REAL.md").read_text(encoding="utf-8")
    assert f"| pearson_motif | pearson_signed | 12 | {rec['diff']:+.3f} |" in report


def test_paper_draft_keeps_each_experiments_macros_apart(tmp_path):
    """Two confirmatory experiments in one --paper draft: no macro is defined twice, and each
    section's prose uses its own experiment's numbers."""
    root = tmp_path / "repo"
    (root / "experiments").mkdir(parents=True)
    (root / "paper").mkdir()
    ledger = ["| # | claim | experiments | files | status | reproduced by |", "|---|---|---|---|---|---|"]
    for i, (exp, status) in enumerate((("EXP-20260101-AA-1", "confirmed"),
                                       ("EXP-20260102-BB-1", "reproduced-computational"))):
        _run(root / "benchmarks" / "results" / exp, n=5 + i)
        (root / "experiments" / f"{exp}.md").write_text(
            f"# {exp}\n\n- type: confirmatory\n"
            "- primary comparison: pearson_motif > pearson_signed | auprc_ratio | auto\n")
        ledger.append(f"| C{i} | claim | {exp} | x | {status} | - |")
    (root / "paper" / "CLAIMS.md").write_text("\n".join(ledger) + "\n")
    main = texreport.build_paper(out=tmp_path / "build", root=root, compile_pdf=False)
    defined = []
    for exp in ("EXP-20260101-AA-1", "EXP-20260102-BB-1"):
        defined += re.findall(r"\\newcommand\{\\(\w+)\}", (tmp_path / "build" / exp / "values.tex").read_text())
    assert defined and len(defined) == len(set(defined))
    tex = main.read_text()
    used = set(re.findall(r"\\(gb[A-Za-z]+)\{\}", tex))
    assert used and used <= set(defined)
    tag_b = texreport.macro_name("EXP-20260102-BB-1")[2:]
    section_b = tex.split("\\section{EXP-20260102-BB-1}")[1]
    assert all(m.startswith("gb" + tag_b) for m in re.findall(r"\\(gb[A-Za-z]+)\{\}", section_b))
