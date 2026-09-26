"""LaTeX reports generated from result files: no number is ever typed into a template.

    grnbench report <run_dir | EXP-ID> [--exp EXP-ID] [--out DIR]   # one experiment
    grnbench report --paper [--out paper/build]                      # draft from paper/CLAIMS.md

How numbers get into the text: every value is computed from the run's ``rows/*.json`` and written
as a macro into ``values.tex`` (``\\newcommand{\\gbDiffAuprcRatioAutoPearsonMotifVsPearsonSigned}{+0.229}``);
tables are written to ``tables/*.tex``; the prose templates below only reference macros and
tables (a test fails if a template contains a digit). Methods text comes from the experiment
record (``experiments/EXP-*.md``), limitations from the record's deviations and from the rows
(doublet fallbacks, model errors, missing datasets). Figures are the PNGs written by ``grnbench viz``.

Statistics are the declared-family path of ``benchmarks/real_bench.py``: ``compare_units`` per
declared pair over datasets, Holm over all directed hypotheses of a metric and arm.

Structure and mechanics follow the end-of-pipeline report of Luca's thesis pipeline
(``save_final_report`` in its ``networkv6.py``, ReportLab): the pipeline itself emits the report as
its last step (``real_bench.py report``), wrapped so a report failure never fails the run; a title
block with the generation time and run folder; a configuration summary built from dicts
(``table_from_dict``, after the thesis ``create_table_from_dict``); performance overall, per unit
(there: per fold; here: per dataset) and summarised across units; a comparative summary; then
visualizations and explainability with captions; and a plain-text fallback when the PDF engine
is missing (there: ``_save_text_report_fallback`` without ReportLab; here: ``report_fallback.md``
without a TeX engine). The output format is LaTeX instead of ReportLab so the tables and text
drop straight into the paper.
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
METRICS = {"auprc_ratio": "AUPRC ratio", "epr": "EPR"}
ARMS = ("auto", "baseline")

REPORT_TEMPLATE = r"""\documentclass{article}
\usepackage{fullpage,booktabs,graphicx,hyperref,float,array}
\input{values.tex}
\title{<<TITLE>>}
\author{grnbench report generated from <<EXP>>}
\date{}
\begin{document}
\maketitle
\section*{Status}
<<STATUS>>
\section{Configuration summary}
\input{tables/configuration.tex}
\section{Methods}
<<METHODS>>
\section{Results}
\subsection{Primary comparison}
<<PRIMARY>>
\subsection{Performance per model, across datasets}
\input{tables/models.tex}
\subsection{Performance per dataset}
\input{tables/per_dataset.tex}
\subsection{Comparative summary: declared comparisons}
<<COMPARISON_TABLES>>
\section{Visualizations and explainability}
<<FIGURES>>
\section{Limitations}
<<LIMITATIONS>>
\end{document}
"""

PAPER_TEMPLATE = r"""\documentclass{article}
\usepackage{fullpage,booktabs,graphicx,hyperref,float,array}
\title{grnbench: draft assembled from confirmatory experiments}
\author{}
\date{}
\begin{document}
\maketitle
\section*{How to read this draft}
Each section below is generated from one experiment listed in \texttt{paper/CLAIMS.md}. Numbers come
from the experiments' result files; claims marked draft are not yet reproduced by a second person.
<<SECTIONS>>
\section*{Background references}
\nocite{*}
\bibliographystyle{unsrt}
\bibliography{references}
\end{document}
"""

TEMPLATES = (REPORT_TEMPLATE, PAPER_TEMPLATE)


# --------------------------------------------------------------------------- helpers

def tex_escape(text: str) -> str:
    rep = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#", "_": r"\_",
           "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}
    out = "".join(rep.get(c, c) for c in str(text))
    return out.replace("`", "")


def macro_name(*parts: str) -> str:
    """'diff', 'auprc_ratio', 'auto', 'pearson_motif' -> 'gbDiffAuprcRatioAutoPearsonMotif'."""
    words = {"0": "Zero", "1": "One", "2": "Two", "3": "Three", "4": "Four", "5": "Five",
             "6": "Six", "7": "Seven", "8": "Eight", "9": "Nine"}
    out = "gb"
    for p in parts:
        for w in re.split(r"[^A-Za-z0-9]+", str(p)):
            w = "".join(words.get(c, c) for c in w)
            out += w[:1].upper() + w[1:]
    return out


def fmt(x: float, signed: bool = False, digits: int = 3) -> str:
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "n/a"
    return f"{x:+.{digits}f}" if signed else f"{x:.{digits}f}"


def fmt_p(p: float) -> str:
    return "n/a" if p is None or not np.isfinite(p) else (f"{p:.2g}" if p < 0.001 else f"{p:.3f}")


def load_rows(run_dir: Path) -> pd.DataFrame:
    rows = sorted((Path(run_dir) / "rows").glob("*.json"))
    if rows:
        return pd.DataFrame([json.loads(p.read_text(encoding="utf-8")) for p in rows])
    csv = next(iter(sorted(Path(run_dir).glob("all_rows_real.csv"))), None)
    if csv is None:
        raise FileNotFoundError(f"{run_dir}: no rows/*.json or all_rows_real.csv")
    return pd.read_csv(csv)


def model_matrix(df: pd.DataFrame, arm: str, metric: str) -> pd.DataFrame:
    suffix = f"_{arm}_{metric}"
    cols = {c[2:-len(suffix)]: c for c in df if c.startswith("m_") and c.endswith(suffix)}
    mat = pd.DataFrame({m: pd.to_numeric(df[c], errors="coerce") for m, c in cols.items()})
    mat.index = df["task"].astype(str)
    return mat


def declared_pairs() -> list[tuple[str, str]]:
    sys.path.insert(0, str(ROOT / "benchmarks"))
    from real_bench import DECLARED_PAIRS  # the family declared in code before a run
    return list(DECLARED_PAIRS)


def comparisons(df: pd.DataFrame, metric: str, arm: str) -> list[dict]:
    """Declared pairs present in the data: compare_units per pair, Holm over the directed family."""
    from grnbench.eval.comparison import compare_units, holm_adjust

    mat = model_matrix(df, arm, metric)
    recs, raw = [], []
    for a, b in declared_pairs():
        if a not in mat or b not in mat:
            continue
        pair = mat[[a, b]].dropna()
        if len(pair) < 3:
            continue
        res = compare_units(pair)["pairs"][0]
        recs.append({"a": a, "b": b, "n": len(pair), "diff": float(res["difference"]),
                     "ci": res["simultaneous_ci"]})
        raw += [res["p_a_gt_b"], res["p_b_gt_a"]]
    if recs:
        adj = holm_adjust(raw)
        for i, r in enumerate(recs):
            r["p_a_gt_b"], r["p_b_gt_a"], r["family"] = float(adj[2 * i]), float(adj[2 * i + 1]), len(raw)
            r["supported"] = r["a"] if r["p_a_gt_b"] < .05 else r["b"] if r["p_b_gt_a"] < .05 else "unresolved"
    return recs


def table_from_dict(data: dict, caption: str) -> str:
    """Two-column key/value table (after the thesis pipeline's ``create_table_from_dict``)."""
    rows = "\n".join(f"{tex_escape(k)} & {tex_escape(v)} \\\\" for k, v in data.items() if v not in (None, ""))
    return ("\\begin{table}[H]\\centering\\small\n"
            f"\\caption{{{tex_escape(caption)}}}\n"
            "\\begin{tabular}{p{0.3\\linewidth}p{0.62\\linewidth}}\\toprule\n"
            f"{rows}\n\\bottomrule\\end{{tabular}}\\end{{table}}\n")


# --------------------------------------------------------------------------- experiment record

def find_record(exp: str, root: Path = ROOT) -> Path | None:
    p = root / "experiments" / f"{exp}.md"
    return p if p.exists() else None


def record_fields(record: Path | None) -> dict[str, str]:
    """'- key: value' bullets of an experiment record (continuation lines appended)."""
    if record is None:
        return {}
    fields, key = {}, None
    for line in record.read_text(encoding="utf-8").splitlines():
        m = re.match(r"^- ([A-Za-z][A-Za-z ()/_-]*?):\s*(.*)$", line)
        if m:
            key = m.group(1).strip().lower()
            fields[key] = m.group(2).strip()
        elif key and line.startswith("  ") and not line.strip().startswith("```"):
            fields[key] += " " + line.strip()
        elif not line.strip() or line.startswith("#"):
            key = None
    return fields


def resolve_input(target: str, root: Path = ROOT) -> tuple[Path, str | None]:
    """A run directory, or an EXP id whose record names its `report input:` folder."""
    p = Path(target)
    if p.is_dir():
        exp = next((part for part in p.resolve().parts[::-1] if part.startswith("EXP-")), None)
        return p, exp
    rec = find_record(target, root)
    if rec is None:
        raise FileNotFoundError(f"{target}: not a folder and no experiments/{target}.md")
    f = record_fields(rec)
    folder = f.get("report input", "").strip("` ")
    if not folder:
        folder = f"benchmarks/results/{target}"
    return root / folder, target


# --------------------------------------------------------------------------- build

def _tables(df: pd.DataFrame, out: Path) -> tuple[dict, list[str], dict]:
    macros: dict[str, str] = {}
    (out / "tables").mkdir(parents=True, exist_ok=True)
    auto, base = model_matrix(df, "auto", "auprc_ratio"), model_matrix(df, "baseline", "auprc_ratio")
    n_rows = len(df)
    macros[macro_name("n", "datasets")] = str(n_rows)
    lines = [r"\begin{table}[H]\centering\small",
             r"\caption{Reference recovery per model (CollecTRI evaluation half; AUPRC ratio, "
             r"1 = random). Unit = dataset; n per model shown.}",
             r"\begin{tabular}{lrrrr}\toprule",
             r"model & n & mean (auto) & median (auto) & mean (baseline) \\ \midrule"]
    for m in auto.mean().sort_values(ascending=False).index:
        a = auto[m].dropna()
        b = base[m].dropna() if m in base else pd.Series(dtype=float)
        macros[macro_name("mean", "auprc_ratio", "auto", m)] = fmt(a.mean())
        macros[macro_name("n", "auto", m)] = str(len(a))
        lines.append(f"{tex_escape(m)} & {len(a)} & {fmt(a.mean())} & {fmt(a.median())} & {fmt(b.mean())} \\\\")
    lines += [r"\bottomrule\end{tabular}\end{table}"]
    (out / "tables" / "models.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # per unit (the thesis report's per-fold table): one row per dataset, auto arm
    cols = list(auto.mean().sort_values(ascending=False).index)
    t = [r"\begin{table}[H]\centering\scriptsize",
         r"\caption{AUPRC ratio per dataset, auto pipeline (1 = random); the unit of analysis.}",
         r"\resizebox{\linewidth}{!}{\begin{tabular}{l" + "r" * len(cols) + r"}\toprule",
         "dataset & " + " & ".join(tex_escape(c) for c in cols) + r" \\ \midrule"]
    for task, r in auto[cols].iterrows():
        name = tex_escape(str(task).replace("cxg__", "")[:48])
        t.append(name + " & " + " & ".join(fmt(v, digits=2) for v in r) + r" \\")
    t.append(r"\midrule")
    for label, fn in (("mean", "mean"), ("median", "median"), ("SD", "std")):
        t.append(label + " & " + " & ".join(fmt(getattr(auto[c].dropna(), fn)(), digits=2) for c in cols) + r" \\")
    t += [r"\bottomrule\end{tabular}}\end{table}"]
    (out / "tables" / "per_dataset.tex").write_text("\n".join(t) + "\n", encoding="utf-8")

    comp_inputs, all_comps = [], {}
    for metric, label in METRICS.items():
        for arm in ARMS:
            recs = comparisons(df, metric, arm)
            if not recs:
                continue
            all_comps[(metric, arm)] = recs
            name = f"comparisons_{metric}_{arm}"
            t = [r"\begin{table}[H]\centering\small",
                 f"\\caption{{Declared comparisons, {tex_escape(label)}, {arm} pipeline: paired t over "
                 f"datasets (compare\\_units), two-sided per direction, Holm over {recs[0]['family']} "
                 f"directed hypotheses; 95\\% interval of the mean paired difference.}}",
                 r"\begin{tabular}{llrrlrrl}\toprule",
                 r"A & B & n & mean A$-$B & 95\% interval & $p_{A>B}$ Holm & $p_{B>A}$ Holm & supported \\ \midrule"]
            for r in recs:
                lo, hi = r["ci"] if r["ci"] and r["ci"][0] is not None else (None, None)
                ci = "unbounded" if lo is None else f"[{fmt(lo, True)}, {fmt(hi, True)}]"
                t.append(f"{tex_escape(r['a'])} & {tex_escape(r['b'])} & {r['n']} & {fmt(r['diff'], True)} & "
                         f"{ci} & {fmt_p(r['p_a_gt_b'])} & {fmt_p(r['p_b_gt_a'])} & {tex_escape(r['supported'])} \\\\")
                key = (metric, arm, r["a"], "vs", r["b"])
                macros[macro_name("diff", *key)] = fmt(r["diff"], True)
                macros[macro_name("ciLow", *key)] = "n/a" if lo is None else fmt(lo, True)
                macros[macro_name("ciHigh", *key)] = "n/a" if hi is None else fmt(hi, True)
                macros[macro_name("pHolm", *key)] = fmt_p(r["p_a_gt_b"])
                macros[macro_name("pHolmRev", *key)] = fmt_p(r["p_b_gt_a"])
                macros[macro_name("n", *key)] = str(r["n"])
                macros[macro_name("family", metric, arm)] = str(r["family"])
            t += [r"\bottomrule\end{tabular}\end{table}"]
            (out / "tables" / f"{name}.tex").write_text("\n".join(t) + "\n", encoding="utf-8")
            comp_inputs.append(f"\\input{{tables/{name}.tex}}")
    return macros, comp_inputs, all_comps


def _primary(fields: dict, comps: dict) -> str:
    """Pre-registered decision rule applied to the data (confirmatory records only)."""
    spec = fields.get("primary comparison", "")
    m = re.match(r"\s*`?([\w]+)\s*>\s*([\w]+)\s*\|\s*([\w]+)\s*\|\s*([\w]+)", spec)
    kind = fields.get("type", "exploratory").split()[0].strip("*").lower()
    if not m:
        return ("No machine-readable primary comparison in the record (field "
                "\\texttt{primary comparison: A > B | metric | arm}); results are descriptive.\n")
    a, b, metric, arm = m.groups()
    rec = next((r for r in comps.get((metric, arm), []) if r["a"] == a and r["b"] == b), None)
    if rec is None:
        return f"The primary comparison {tex_escape(a)} vs {tex_escape(b)} has no data in this run.\n"
    key = (metric, arm, a, "vs", b)
    verdict = ("supported" if rec["p_a_gt_b"] < .05 else "refuted" if rec["p_b_gt_a"] < .05
               else "inconclusive (not resolved; this is not evidence of equality)")
    label = "Pre-registered primary hypothesis" if kind == "confirmatory" else "Primary comparison (exploratory)"
    return (f"\\paragraph{{{label}.}} {tex_escape(a)} $>$ {tex_escape(b)} on {tex_escape(METRICS.get(metric, metric))}, "
            f"{tex_escape(arm)} pipeline. Mean paired difference \\{macro_name('diff', *key)}{{}} "
            f"(95\\% interval \\{macro_name('ciLow', *key)}{{}} to \\{macro_name('ciHigh', *key)}{{}}), "
            f"n = \\{macro_name('n', *key)}{{}} datasets, paired t, Holm over "
            f"\\{macro_name('family', metric, arm)}{{}} directed hypotheses: "
            f"$p_{{\\mathrm{{Holm}}}}$ = \\{macro_name('pHolm', *key)}{{}}. Verdict by the rule fixed in the "
            f"record: \\textbf{{{verdict}}}."
            + ("" if kind == "confirmatory" else " Exploratory: this cannot support a headline claim.")
            + "\n")


def _methods(fields: dict, df: pd.DataFrame) -> str:
    keys = ["type", "hypothesis", "question", "primary endpoint", "secondary endpoints (reported, not used for h1)",
            "unit of analysis", "datasets", "exclusions", "exclusion rules (fixed now)", "comparison family",
            "test", "seeds", "models", "code", "environment", "data", "machine"]
    items = [f"\\item[{tex_escape(k)}] {tex_escape(fields[k])}" for k in keys if fields.get(k)]
    models = sorted({c[2:].rsplit("_auto_", 1)[0] for c in df if c.startswith("m_") and "_auto_auprc_ratio" in c})
    items.append(f"\\item[models scored] {tex_escape(', '.join(models))}")
    return "\\begin{description}\n" + "\n".join(items) + "\n\\end{description}\n"


def _limitations(fields: dict, df: pd.DataFrame, exploratory: bool) -> str:
    items = []
    if exploratory:
        items.append("Exploratory experiment: results need a confirmatory rerun before a headline claim.")
    for k in ("deviations", "deviations from section 1"):
        if fields.get(k) and not fields[k].lower().startswith("none"):
            items.append("Recorded deviation: " + tex_escape(fields[k]))
    dbl = [c for c in df if c.endswith("_doublets")]
    if dbl:
        fb = df[dbl].astype(str).apply(lambda s: s.str.contains("heuristic")).any(axis=1)
        if fb.any():
            items.append(f"Doublet calling fell back to the heuristic in {int(fb.sum())} of {len(df)} datasets.")
    else:
        items.append("The run did not record which doublet caller ran (rows predate the fallback record).")
    if "error" in df and df["error"].notna().any():
        items.append(f"{int(df['error'].notna().sum())} datasets errored: "
                     + tex_escape(", ".join(df.loc[df['error'].notna(), 'task'].astype(str))))
    errs = [c for c in df if c.startswith("m_") and c.endswith("_error")]
    for c in errs:
        n = int(df[c].notna().sum())
        if n:
            items.append(f"Model error in {n} dataset(s): {tex_escape(c[2:-6])}.")
    if not items:
        items.append("None recorded.")
    return "\\begin{itemize}\n" + "\n".join(f"\\item {i}" for i in items) + "\n\\end{itemize}\n"


def _figures(run_dir: Path, out: Path) -> str:
    viz = next((d for d in (run_dir / "viz", run_dir.parent / "viz", out / "viz") if d.is_dir()), None)
    if viz is None:
        return "No figures: run \\texttt{grnbench viz} on the run folder first.\n"
    index = (viz / "index.md").read_text(encoding="utf-8") if (viz / "index.md").exists() else ""
    blocks = []
    for png in sorted(viz.rglob("*.png")):
        rel = png.relative_to(viz).as_posix()
        line = next((l for l in index.splitlines() if rel in l or png.name in l), "")
        cap = re.sub(r"!\[[^\]]*\]\([^)]*\)|\[[^\]]*\]\([^)]*\)|[`*|#]", " ", line).strip() or png.stem
        # referenced in place (no copy): path relative to the .tex file
        src = os.path.relpath(png, out).replace(os.sep, "/")
        blocks.append("\\begin{figure}[H]\\centering\n"
                      f"\\includegraphics[width=0.85\\linewidth]{{{src}}}\n"
                      f"\\caption{{{tex_escape(cap)}}}\n\\end{{figure}}")
    return "\n".join(blocks) + "\n" if blocks else "No figures found in the viz folder.\n"


def build_report(target: str, exp: str | None = None, out: Path | None = None,
                 root: Path = ROOT, compile_pdf: bool = True) -> Path:
    run_dir, exp_found = resolve_input(target, root)
    exp = exp or exp_found or run_dir.name
    fields = record_fields(find_record(exp, root))
    df = load_rows(run_dir)
    out = Path(out or run_dir / "latex")
    out.mkdir(parents=True, exist_ok=True)
    macros, comp_inputs, comps = _tables(df, out)
    commit_file = next((f for f in (run_dir / "CODE_COMMIT", run_dir.parent / "CODE_COMMIT") if f.exists()), None)
    config = {"experiment": exp, "type": fields.get("type", "not recorded"), "owner": fields.get("owner"),
              "datasets (rows)": str(len(df)), "models": ", ".join(sorted(
                  {c[2:].rsplit("_auto_", 1)[0] for c in df if c.startswith("m_") and c.endswith("_auto_auprc_ratio")})),
              "arms": ", ".join(ARMS), "code": fields.get("code"),
              "code commit (run folder)": commit_file.read_text(encoding="utf-8").strip() if commit_file else None,
              "environment": fields.get("environment"), "machine": fields.get("machine"),
              "seeds": fields.get("seeds"), "run folder": _rel(run_dir, root),
              "generated": datetime.now().strftime("%Y-%m-%d %H:%M")}
    (out / "tables" / "configuration.tex").write_text(
        table_from_dict(config, "Configuration of the run, from the experiment record and the run folder."),
        encoding="utf-8")
    (out / "values.tex").write_text(
        "% generated by grnbench report from result files - do not edit\n"
        + "\n".join(f"\\newcommand{{\\{k}}}{{{v}}}" for k, v in sorted(macros.items())) + "\n",
        encoding="utf-8")
    kind = fields.get("type", "exploratory").split()[0].strip("*").lower()
    status = (f"Experiment {tex_escape(exp)}, type: {tex_escape(kind)}. "
              f"Record status: {tex_escape(fields.get('status', 'not recorded'))}. "
              f"Reproduced by: {tex_escape(fields.get('reproduced by', 'nobody yet'))}.")
    title = next((l[2:].strip() for l in (find_record(exp, root).read_text(encoding="utf-8").splitlines()
                                            if find_record(exp, root) else []) if l.startswith("# ")), exp)
    tex = (REPORT_TEMPLATE.replace("<<TITLE>>", tex_escape(title)).replace("<<EXP>>", tex_escape(exp))
           .replace("<<STATUS>>", status).replace("<<METHODS>>", _methods(fields, df))
           .replace("<<PRIMARY>>", _primary(fields, comps))
           .replace("<<COMPARISON_TABLES>>", "\n".join(comp_inputs))
           .replace("<<FIGURES>>", _figures(run_dir, out))
           .replace("<<LIMITATIONS>>", _limitations(fields, df, kind != "confirmatory")))
    main = out / "report.tex"
    main.write_text(tex, encoding="utf-8")
    if compile_pdf and compile_tex(main) is None:
        _text_fallback(out, exp, config, comps)
    return main


def _rel(path: Path, root: Path) -> str:
    try:
        return Path(path).resolve().relative_to(Path(root).resolve()).as_posix()
    except ValueError:
        return Path(path).name


def _text_fallback(out: Path, exp: str, config: dict, comps: dict) -> Path:
    """Plain-text report when no TeX engine is installed (the thesis pipeline did the same when
    ReportLab was missing): configuration and every declared comparison, from the same numbers."""
    lines = [f"# {exp}: report (text fallback, no TeX engine found)", ""]
    lines += [f"- {k}: {v}" for k, v in config.items() if v not in (None, "")]
    for (metric, arm), recs in comps.items():
        lines += ["", f"## {METRICS.get(metric, metric)}, {arm} pipeline (paired t, Holm over "
                      f"{recs[0]['family']} directed hypotheses)", "",
                  "| A | B | n | mean A-B | 95% interval | p(A>B) Holm | p(B>A) Holm | supported |",
                  "|---|---|---|---|---|---|---|---|"]
        for r in recs:
            lo, hi = r["ci"] if r["ci"] and r["ci"][0] is not None else (None, None)
            ci = "unbounded" if lo is None else f"[{fmt(lo, True)}, {fmt(hi, True)}]"
            lines.append(f"| {r['a']} | {r['b']} | {r['n']} | {fmt(r['diff'], True)} | {ci} | "
                         f"{fmt_p(r['p_a_gt_b'])} | {fmt_p(r['p_b_gt_a'])} | {r['supported']} |")
    path = out / "report_fallback.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_bib(path: Path) -> int:
    from grnbench.kb import REFERENCES
    entries = []
    for r in REFERENCES.values():
        field = "doi" if "doi.org/" in r.link else "url"
        val = r.link.split("doi.org/", 1)[1] if field == "doi" else r.link
        entries.append(f"@misc{{{r.key},\n  note = {{{tex_escape(r.citation)}}},\n  {field} = {{{val}}}\n}}")
    path.write_text("\n\n".join(entries) + "\n", encoding="utf-8")
    return len(entries)


def build_paper(out: Path | None = None, root: Path = ROOT, compile_pdf: bool = True) -> Path:
    """One draft from the confirmatory experiments listed in paper/CLAIMS.md."""
    out = Path(out or root / "paper" / "build")
    out.mkdir(parents=True, exist_ok=True)
    claims = (root / "paper" / "CLAIMS.md").read_text(encoding="utf-8")
    exps = []
    for line in claims.splitlines():
        cells = [c.strip() for c in line.split("|")]
        if len(cells) > 6 and cells[5] in ("confirmed", "reproduced-computational", "reproduced"):
            exps += [e for e in re.findall(r"EXP-[\w-]+", cells[3]) if e not in exps]
    sections = []
    for exp in exps:
        fields = record_fields(find_record(exp, root))
        if not fields.get("type", "").lower().startswith(("confirmatory", "**confirmatory")):
            continue
        sub = build_report(exp, out=out / exp, root=root, compile_pdf=False)
        # every experiment defines the same macro names: prefix them per experiment
        # (\gbDiff... -> \gbEXP...LMTwoDiff...) so sections do not redefine each other's numbers
        tag = macro_name(exp)[2:]
        for f in (out / exp / "values.tex", sub):
            f.write_text(re.sub(r"\\gb(?=[A-Z])", lambda _m: "\\gb" + tag, f.read_text(encoding="utf-8")),
                         encoding="utf-8")
        body = sub.read_text(encoding="utf-8")
        body = body.split("\\maketitle", 1)[1].split("\\end{document}", 1)[0]
        body = body.replace("\\input{", f"\\input{{{exp}/")
        body = re.sub(r"(\\includegraphics\[[^\]]*\]\{)", lambda m: m.group(1) + f"{exp}/", body)
        body = body.replace("\\section{", "\\subsection{").replace("\\section*{", "\\subsection*{")
        title = tex_escape(exp)
        sections.append(f"\\input{{{exp}/values.tex}}\n\\section{{{title}}}\n{body}")
    n_refs = write_bib(out / "references.bib")
    tex = PAPER_TEMPLATE.replace("<<SECTIONS>>", "\n".join(sections) if sections else
                                 "No confirmed or reproduced claim in paper/CLAIMS.md yet.\n")
    main = out / "paper.tex"
    main.write_text(tex, encoding="utf-8")
    print(f"paper draft: {len(sections)} experiment section(s), {n_refs} references")
    if compile_pdf:
        compile_tex(main)
    return main


def compile_tex(main: Path) -> Path | None:
    """latexmk, else tectonic, else pdflatex; None (with instructions) when no engine is installed."""
    main = Path(main)
    for cmd in (["latexmk", "-pdf", "-interaction=nonstopmode", "-quiet", main.name],
                ["tectonic", main.name],
                ["pdflatex", "-interaction=nonstopmode", main.name]):
        if shutil.which(cmd[0]):
            subprocess.run(cmd, cwd=main.parent, capture_output=True, text=True)
            pdf = main.with_suffix(".pdf")
            if pdf.exists():
                print(f"compiled {pdf}")
                return pdf
    print(f"no TeX engine found: wrote {main}. Build it with one of: `latexmk -pdf {main.name}` "
          f"(TeX Live / MiKTeX), `tectonic {main.name}`, or `quarto install tinytex` then latexmk, "
          f"run inside {main.parent}")
    return None
