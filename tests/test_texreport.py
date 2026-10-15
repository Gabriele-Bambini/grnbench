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


def test_co_primary_hypotheses_each_get_a_verdict(tmp_path):
    """A record with two primary comparisons (';'-separated, one Holm family) gets one
    pre-registered verdict paragraph per hypothesis, not only the first."""
    root = tmp_path / "repo"
    exp = "EXP-20260103-CC-1"
    run = root / "benchmarks" / "results" / exp
    _run(run)
    (root / "experiments").mkdir(parents=True)
    (root / "experiments" / f"{exp}.md").write_text(
        f"# {exp}\n\n- type: confirmatory\n"
        "- primary comparison: pearson_motif > pearson_signed | auprc_ratio | auto;"
        " pearson_motif > motif_prior | auprc_ratio | auto\n"
        f"- report input: `benchmarks/results/{exp}`\n")
    tex = texreport.build_report(exp, root=root, out=tmp_path / "out", compile_pdf=False).read_text()
    assert tex.count("Pre-registered primary hypothesis") == 2
    for b in (r"pearson\_signed", r"motif\_prior"):
        assert rf"pearson\_motif $>$ {b}" in tex


def test_record_fields_handle_commas_code_and_fences(tmp_path):
    """Keys with commas or backticks are keys, not continuation text; fenced code is not appended
    to the key above it (the Methods section of EXP-20260926-LM-1 was garbled by both)."""
    rec = tmp_path / "EXP.md"
    rec.write_text(
        "- primary endpoint: AUPRC ratio\n"
        "- secondary endpoints (reported, not used for H1/H2): EPR; the baseline arm\n"
        "- environment: conda env\n"
        "- data (login node, downloads only; `D=~/x`): fetched\n"
        "  ```\n"
        "  python fetch.py --out $D\n"
        "  ```\n"
        "- machine: HPC\n", encoding="utf-8")
    f = texreport.record_fields(rec)
    assert f["primary endpoint"] == "AUPRC ratio"
    assert f["secondary endpoints (reported, not used for h1/h2)"] == "EPR; the baseline arm"
    assert f["environment"] == "conda env"
    assert f["data (login node, downloads only; `d=~/x`)"] == "fetched"
    assert f["machine"] == "HPC"
    methods = texreport._methods(f, __import__("pandas").DataFrame())
    assert "EPR; the baseline arm" in methods and "fetch.py" not in methods


def test_report_status_strips_record_markdown(tmp_path):
    root = tmp_path / "repo"
    exp = "EXP-20260105-EE-1"
    _run(root / "benchmarks" / "results" / exp)
    (root / "experiments").mkdir(parents=True)
    (root / "experiments" / f"{exp}.md").write_text(
        f"# {exp}\n\n- type: **exploratory**.\n"
        "- status: **done**. Complete.\n"
        "- reproduced by: **pending**.\n"
        f"- report input: `benchmarks/results/{exp}`\n")

    tex = texreport.build_report(exp, root=root, out=tmp_path / "out",
                                 compile_pdf=False).read_text()

    assert "type: exploratory." in tex
    assert "Record status: done. Complete." in tex
    assert "Reproduced by: pending." in tex
    assert "**" not in tex


def test_primary_comparison_note_with_semicolon_is_not_a_second_hypothesis(tmp_path):
    """EXP-20260925-LM-2's field continues with '(...; no content change)': only fragments that
    parse as 'A > B | metric | arm' are hypotheses."""
    root = tmp_path / "repo"
    exp = "EXP-20260104-DD-1"
    _run(root / "benchmarks" / "results" / exp)
    (root / "experiments").mkdir(parents=True)
    (root / "experiments" / f"{exp}.md").write_text(
        f"# {exp}\n\n- type: confirmatory\n"
        "- primary comparison: pearson_motif > pearson_signed | auprc_ratio | auto\n"
        "  (restatement added before the run; no content change)\n"
        f"- report input: `benchmarks/results/{exp}`\n")
    tex = texreport.build_report(exp, root=root, out=tmp_path / "out", compile_pdf=False).read_text()
    assert tex.count("Pre-registered primary hypothesis") == 1
    assert "No machine-readable primary comparison" not in tex


def test_report_input_ignores_text_after_the_quoted_path(tmp_path):
    """EXP-20261002-LM-1's record ends the path with '`.'; the template adds prose after it."""
    root = tmp_path / "repo"
    exp = "EXP-20260105-DD-1"
    _run(root / "benchmarks" / "results" / exp)
    (root / "experiments").mkdir(parents=True)
    (root / "experiments" / f"{exp}.md").write_text(
        f"# {exp}\n\n- type: confirmatory\n"
        f"- report input: `benchmarks/results/{exp}`. (folder with rows/*.json)\n")
    assert texreport.resolve_input(exp, root=root) == (root / "benchmarks" / "results" / exp, exp)
    assert texreport.build_report(exp, root=root, out=tmp_path / "out", compile_pdf=False).exists()


def _fake_engines(monkeypatch, installed, works):
    """Fake TeX engines: `installed` are found on PATH, only those in `works` write the PDF."""
    import types

    monkeypatch.setattr(texreport.shutil, "which", lambda n: n if n in installed else None)

    def run(cmd, cwd, **kw):
        if cmd[0] in works:
            (Path(cwd) / "main.pdf").write_bytes(b"%PDF")
        return types.SimpleNamespace(stdout=f"! {cmd[0]} error\n", stderr="")

    monkeypatch.setattr(texreport.subprocess, "run", run)


def test_compile_tex_falls_through_to_the_next_engine(tmp_path, monkeypatch, capsys):
    """A failing first engine (latexmk without Perl) must not stop a working pdflatex."""
    _fake_engines(monkeypatch, {"latexmk", "pdflatex"}, works={"pdflatex"})
    assert texreport.compile_tex(tmp_path / "main.tex") == tmp_path / "main.pdf"
    assert "wrote no PDF" not in capsys.readouterr().out


def test_compile_tex_reports_every_failed_engine_once(tmp_path, monkeypatch, capsys):
    _fake_engines(monkeypatch, {"latexmk", "pdflatex"}, works=set())
    assert texreport.compile_tex(tmp_path / "main.tex") is None
    out = capsys.readouterr().out
    assert "latexmk ran but wrote no PDF" in out and "pdflatex ran but wrote no PDF" in out
    assert "! latexmk error" in out and "! pdflatex error" in out
    assert "no TeX engine found" not in out
