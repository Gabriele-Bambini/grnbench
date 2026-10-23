"""Self-contained HTML comparison report for a benchmarking run.

No CDN, no JavaScript framework, no plotting dependency: charts are inline SVG
built here, so the file opens from a scratch directory on a compute node with
no network. ``matplotlib`` is never imported.

Sections, in the order a reader needs them: the verdict table, PR curves, what
the nulls say, whether methods actually agree with each other (top-k Jaccard),
subsample stability, gold reachability, the preprocessing QC gates and the full
provenance block.
"""
from __future__ import annotations

import html
import math
from datetime import datetime, timezone

# colour-blind-safe categorical palette (Okabe-Ito), stable across reports
PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7",
           "#E69F00", "#56B4E9", "#F0E442", "#000000"]

_CSS = """
:root { --fg:#1a1d21; --muted:#5b6570; --line:#dfe3e8; --bg:#ffffff;
        --head:#f5f7f9; --warn:#8a6d00; --fail:#b3261e; --pass:#0f7b34; }
* { box-sizing:border-box; }
body { margin:0; padding:2rem 1.5rem 4rem; background:var(--bg); color:var(--fg);
       font:14px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,
       "Helvetica Neue",Arial,sans-serif; }
.wrap { max-width:1080px; margin:0 auto; }
h1 { font-size:1.6rem; margin:0 0 .25rem; letter-spacing:-.01em; }
h2 { font-size:1.1rem; margin:2.4rem 0 .6rem; padding-bottom:.3rem;
     border-bottom:1px solid var(--line); }
h3 { font-size:.95rem; margin:1.4rem 0 .4rem; color:var(--muted);
     text-transform:uppercase; letter-spacing:.06em; }
p.sub { color:var(--muted); margin:.1rem 0 0; }
.meta { display:flex; flex-wrap:wrap; gap:.4rem 1.6rem; margin:1rem 0 0;
        color:var(--muted); font-size:.86rem; }
.meta b { color:var(--fg); font-weight:600; }
.scroll { overflow-x:auto; }
table { border-collapse:collapse; width:100%; font-variant-numeric:tabular-nums;
        margin:.4rem 0 0; }
th,td { padding:.42rem .6rem; text-align:right; border-bottom:1px solid var(--line);
        white-space:nowrap; }
th:first-child,td:first-child { text-align:left; }
thead th { background:var(--head); font-weight:600; font-size:.82rem;
           color:var(--muted); text-transform:uppercase; letter-spacing:.04em; }
tbody tr:hover { background:#fafbfc; }
td.best { font-weight:700; }
code { background:var(--head); padding:.08rem .3rem; border-radius:3px;
       font-size:.86em; }
.badge { display:inline-block; padding:.05rem .45rem; border-radius:999px;
         font-size:.75rem; font-weight:700; letter-spacing:.03em; }
.pass { background:#e6f4ea; color:var(--pass); }
.warn { background:#fdf3d7; color:var(--warn); }
.fail { background:#fce8e6; color:var(--fail); }
.note { background:#fdf3d7; border-left:3px solid #d9a400; padding:.6rem .8rem;
        margin:.7rem 0; border-radius:0 4px 4px 0; }
.legend { display:flex; flex-wrap:wrap; gap:.2rem 1.1rem; margin:.5rem 0 0;
          font-size:.84rem; }
.legend span { display:inline-flex; align-items:center; gap:.35rem; }
.swatch { width:11px; height:11px; border-radius:2px; display:inline-block; }
.hm td { text-align:center; color:#fff; font-size:.8rem; }
.hm td.diag { color:var(--muted); background:var(--head)!important; }
footer { margin-top:3rem; padding-top:1rem; border-top:1px solid var(--line);
         color:var(--muted); font-size:.82rem; }
details { margin:.5rem 0; }
summary { cursor:pointer; color:var(--muted); font-size:.86rem; }
pre { background:var(--head); padding:.7rem .9rem; border-radius:4px;
      overflow-x:auto; font-size:.8rem; margin:.5rem 0 0; }
"""


# --------------------------------------------------------------------------- #
# small formatting helpers
# --------------------------------------------------------------------------- #

def _e(v) -> str:
    return html.escape(str(v), quote=True)


def _is_nan(v) -> bool:
    return isinstance(v, float) and math.isnan(v)


def _num(v, spec: str = ".4f", dash: str = "&ndash;") -> str:
    if v is None or _is_nan(v):
        return dash
    try:
        return format(float(v), spec)
    except (TypeError, ValueError):
        return _e(v)


def _best(values: dict, higher_is_better: bool = True) -> set:
    """Names holding the best finite value (ties all marked)."""
    finite = {k: v for k, v in values.items()
              if v is not None and not _is_nan(v) and isinstance(v, (int, float))}
    if not finite:
        return set()
    target = max(finite.values()) if higher_is_better else min(finite.values())
    return {k for k, v in finite.items() if v == target}


# --------------------------------------------------------------------------- #
# charts (hand-rolled inline SVG)
# --------------------------------------------------------------------------- #

def _pr_curves_svg(methods: dict, width: int = 660, height: int = 330) -> str:
    """Overlaid precision-recall curves with the random-ranking baseline."""
    curves = {n: m["pr_curve"] for n, m in methods.items()
              if m.get("pr_curve") and m["pr_curve"].get("recall")}
    if not curves:
        return '<p class="sub">No PR curve: no method scored a single gold hit.</p>'

    pad_l, pad_r, pad_t, pad_b = 52, 14, 14, 40
    w = width - pad_l - pad_r
    h = height - pad_t - pad_b
    ymax = max(
        [max(c["precision"]) for c in curves.values() if c["precision"]] + [0.0])
    ymax = max(ymax * 1.15, 0.02)

    def sx(r: float) -> float:
        return pad_l + r * w

    def sy(p: float) -> float:
        return pad_t + h - min(p / ymax, 1.0) * h

    parts = [f'<svg viewBox="0 0 {width} {height}" width="100%" '
             f'role="img" aria-label="precision-recall curves" '
             f'style="max-width:{width}px;font:11px sans-serif">']
    # grid + axes
    for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
        y = pad_t + h - frac * h
        parts.append(f'<line x1="{pad_l}" y1="{y:.1f}" x2="{pad_l + w}" y2="{y:.1f}" '
                     f'stroke="#eef1f4"/>')
        parts.append(f'<text x="{pad_l - 8}" y="{y + 3.5:.1f}" text-anchor="end" '
                     f'fill="#5b6570">{frac * ymax:.3f}</text>')
        x = pad_l + frac * w
        parts.append(f'<text x="{x:.1f}" y="{pad_t + h + 16}" text-anchor="middle" '
                     f'fill="#5b6570">{frac:.2f}</text>')
    parts.append(f'<line x1="{pad_l}" y1="{pad_t + h}" x2="{pad_l + w}" '
                 f'y2="{pad_t + h}" stroke="#9aa4ae"/>')
    parts.append(f'<line x1="{pad_l}" y1="{pad_t}" x2="{pad_l}" y2="{pad_t + h}" '
                 f'stroke="#9aa4ae"/>')
    parts.append(f'<text x="{pad_l + w / 2:.0f}" y="{height - 6}" '
                 f'text-anchor="middle" fill="#5b6570">recall</text>')
    parts.append(f'<text transform="translate(13,{pad_t + h / 2:.0f}) rotate(-90)" '
                 f'text-anchor="middle" fill="#5b6570">precision</text>')

    baseline = max((c.get("baseline", 0.0) for c in curves.values()), default=0.0)
    if baseline > 0:
        yb = sy(baseline)
        parts.append(f'<line x1="{pad_l}" y1="{yb:.1f}" x2="{pad_l + w}" y2="{yb:.1f}" '
                     f'stroke="#9aa4ae" stroke-dasharray="5 4"/>')
        parts.append(f'<text x="{pad_l + w - 4}" y="{yb - 5:.1f}" text-anchor="end" '
                     f'fill="#5b6570">prevalence ({baseline:.4f})</text>')

    for i, (name, c) in enumerate(curves.items()):
        colour = PALETTE[i % len(PALETTE)]
        pts = " ".join(f"{sx(r):.1f},{sy(p):.1f}"
                       for r, p in zip(c["recall"], c["precision"]))
        parts.append(f'<polyline points="{pts}" fill="none" stroke="{colour}" '
                     f'stroke-width="2" stroke-linejoin="round"/>')
    parts.append("</svg>")

    legend = "".join(
        f'<span><i class="swatch" style="background:{PALETTE[i % len(PALETTE)]}"></i>'
        f'{_e(n)}</span>'
        for i, n in enumerate(curves)
    )
    return "".join(parts) + f'<div class="legend">{legend}</div>'


def _bars_svg(values: dict, label: str, width: int = 660,
              row_h: int = 26) -> str:
    """Horizontal bar chart for one metric across methods."""
    finite = {k: float(v) for k, v in values.items()
              if v is not None and not _is_nan(v) and isinstance(v, (int, float))}
    if not finite:
        return ""
    vmax = max(max(finite.values()), 1e-9)
    pad_l, pad_r = 150, 60
    w = width - pad_l - pad_r
    height = row_h * len(finite) + 24
    parts = [f'<svg viewBox="0 0 {width} {height}" width="100%" role="img" '
             f'aria-label="{_e(label)} by method" '
             f'style="max-width:{width}px;font:11px sans-serif">']
    for i, (name, v) in enumerate(finite.items()):
        y = 8 + i * row_h
        bw = max((v / vmax) * w, 1.0)
        colour = PALETTE[i % len(PALETTE)]
        parts.append(f'<text x="{pad_l - 8}" y="{y + 13}" text-anchor="end" '
                     f'fill="#1a1d21">{_e(name)}</text>')
        parts.append(f'<rect x="{pad_l}" y="{y + 3}" width="{bw:.1f}" height="15" '
                     f'rx="2" fill="{colour}" opacity="0.88"/>')
        parts.append(f'<text x="{pad_l + bw + 6:.1f}" y="{y + 15}" '
                     f'fill="#5b6570">{v:.3f}</text>')
    parts.append("</svg>")
    return "".join(parts)


def _heatmap(overlap: dict) -> str:
    """Top-k Jaccard matrix rendered as a coloured table."""
    names = overlap.get("methods") or []
    jac = overlap.get("jaccard") or {}
    if len(names) < 2:
        return ""
    head = "".join(f"<th>{_e(n)}</th>" for n in names)
    rows = []
    for a in names:
        cells = []
        for b in names:
            v = jac.get(a, {}).get(b)
            if a == b:
                cells.append('<td class="diag">&mdash;</td>')
                continue
            if v is None or _is_nan(v):
                cells.append('<td class="diag">&ndash;</td>')
                continue
            # single-hue ramp: pale for disagreement, saturated for overlap
            alpha = 0.10 + 0.85 * float(v)
            cells.append(f'<td style="background:rgba(0,114,178,{alpha:.2f})">'
                         f'{v:.3f}</td>')
        rows.append(f"<tr><th>{_e(a)}</th>{''.join(cells)}</tr>")
    return (f'<div class="scroll"><table class="hm"><thead><tr><th></th>{head}</tr>'
            f'</thead><tbody>{"".join(rows)}</tbody></table></div>')


# --------------------------------------------------------------------------- #
# sections
# --------------------------------------------------------------------------- #

def _tie_hint(m: dict) -> str:
    """Show the EPR interval when the cut at k landed inside tied scores.

    The headline EPR is the expectation under a random tie-break. If the method
    scored a large block of edges identically and k falls inside it, that single
    number is a summary of a range the method never resolved - so print the
    range next to it instead of letting a coin flip look like a measurement.
    """
    block = m.get("tie_block_at_k") or 0
    rng = m.get("epr_range") or []
    if block <= 1 or len(rng) != 2 or any(_is_nan(v) for v in rng):
        return ""
    if abs(rng[1] - rng[0]) < 1e-9:
        return ""
    return (f'<br><span class="sub" style="font-size:.78rem">'
            f"{rng[0]:.2f}&ndash;{rng[1]:.2f} &middot; {block:,} tied</span>")


def _verdict_table(methods: dict) -> str:
    if not methods:
        return '<p class="sub">No method produced a scorable edge list.</p>'
    best_auprc = _best({n: m.get("auprc") for n, m in methods.items()})
    best_epr = _best({n: m.get("epr") for n, m in methods.items()})
    best_ratio = _best({n: m.get("auprc_ratio") for n, m in methods.items()})
    rows = []
    for name, m in methods.items():
        null = m.get("null") or {}
        stab = m.get("stability") or {}
        p = null.get("p_value")
        p_cell = _num(p, ".3f")
        # never decorate a p-value that describes an undefined observation: a
        # method with no gold hit at all must not wear a green "significant" pill
        if null.get("observed_undefined"):
            p_cell = '<span class="badge fail">n/a</span>'
        rows.append(
            "<tr>"
            f"<td>{_e(name)}</td>"
            f"<td>{m.get('n_edges', 0):,}</td>"
            f"<td>{m.get('n_gold_in_list', 0):,}</td>"
            f"<td>{m.get('missing_predictions', 0):,}</td>"
            f"<td>{_e(m.get('support_origin', 'unspecified'))}</td>"
            f"<td class=\"{'best' if name in best_auprc else ''}\">"
            f"{_num(m.get('auprc'))}</td>"
            f"<td class=\"{'best' if name in best_ratio else ''}\">"
            f"{_num(m.get('auprc_ratio'), '.2f')}</td>"
            f"<td class=\"{'best' if name in best_epr else ''}\">"
            f"{_num(m.get('epr'), '.2f')}{_tie_hint(m)}</td>"
            f"<td>{_num(m.get('auroc'))}</td>"
            f"<td>{p_cell}</td>"
            f"<td>{_num(stab.get('jaccard_mean'), '.3f')}</td>"
            "</tr>"
        )
    return (
        '<div class="scroll"><table><thead><tr>'
        "<th>method</th><th>candidates</th><th>predicted gold hits</th><th>missing predictions</th><th>support origin</th><th>AUPRC</th>"
        "<th>AUPRC ratio</th><th>EPR</th><th>AUROC</th><th>null p</th>"
        "<th>stability</th>"
        "</tr></thead><tbody>" + "".join(rows) + "</tbody></table></div>"
        '<p class="sub">AUPRC ratio and EPR divide by positive-label prevalence. '
        "Prevalence is a conventional AP baseline, not its exact finite-sample random expectation. "
        "Null p tests the stated score-exchangeability hypothesis, unadjusted across methods; "
        "degree rewiring has no calibrated p-value. Unlabelled candidates are not certified biological negatives. "
        "Missing predictions form a tied bottom block on explicit U. Without U, scores are conditional "
        "on the shared output support. Stability uses complete boundary ties across cell subsamples.</p>"
    )


def _nulls_section(methods: dict) -> str:
    withnull = {n: m["null"] for n, m in methods.items() if m.get("null")}
    if not withnull:
        return '<p class="sub">Null models were skipped for this run.</p>'
    rows = "".join(
        "<tr>"
        f"<td>{_e(n)}</td><td>{_e(v.get('null_hypothesis', 'unspecified legacy null'))}</td><td>{_num(v.get('observed'))}</td>"
        f"<td>{_num(v.get('null_mean'))}</td><td>{_num(v.get('null_std'))}</td>"
        f"<td>{_num(v.get('p_value'), '.4f')}</td>"
        f"<td>{v.get('n_nulls_used', 0)}</td>"
        f"<td>{v.get('n_nulls_degenerate', 0)}</td>"
        "</tr>"
        for n, v in withnull.items()
    )
    degenerate = sum(v.get("n_nulls_degenerate", 0) for v in withnull.values())
    note = ""
    if degenerate:
        note = ('<div class="note"><b>Degenerate nulls counted, not hidden.</b> '
                f"{degenerate} shuffle(s) were undefined or inadmissible. "
                "Incomplete score-permutation distributions do not yield a p-value. "
                "Degree rewiring is diagnostic only.</div>")
    return (
        '<div class="scroll"><table><thead><tr><th>method</th><th>hypothesis</th><th>observed</th>'
        "<th>null mean</th><th>null sd</th><th>p</th><th>nulls used</th>"
        "<th>degenerate</th></tr></thead><tbody>" + rows + "</tbody></table></div>"
        + note
    )


def _stability_section(methods: dict) -> str:
    stab = {n: m["stability"] for n, m in methods.items() if m.get("stability")}
    if not stab:
        return ('<p class="sub">Stability not run. Pass <code>--stability-data '
                "processed.h5ad</code> to <code>grnbench-eval</code> to re-run each "
                "method on cell subsamples.</p>")
    any_cfg = next(iter(stab.values()))
    rows = "".join(
        "<tr>"
        f"<td>{_e(n)}</td><td>{_num(v.get('jaccard_mean'), '.4f')}</td>"
        f"<td>{_num(v.get('jaccard_std'), '.4f')}</td>"
        f"<td>{v.get('n_subsamples', 0)}</td>"
        f"<td>{_num(v.get('subsample_frac'), '.2f')}</td>"
        f"<td>{v.get('top_k', 0)}</td>"
        "</tr>"
        for n, v in stab.items()
    )
    return (
        f'<p class="sub">Each method re-run from scratch on '
        f'{_num(any_cfg.get("subsample_frac"), ".0%")} cell subsamples; the score '
        "is the mean pairwise Jaccard including every tied edge at the top-k boundary. "
        "This describes sensitivity to cell sampling; it does not prove biological validity. "
        "The algorithm seed is fixed by default; varying it additionally measures algorithm randomness.</p>"
        '<div class="scroll"><table><thead><tr><th>method</th><th>Jaccard mean</th>'
        "<th>sd</th><th>subsamples</th><th>frac</th><th>top-k</th>"
        "</tr></thead><tbody>" + rows + "</tbody></table></div>"
        + _bars_svg({n: v.get("jaccard_mean") for n, v in stab.items()},
                    "subsample stability")
    )


def _coverage_section(coverage: dict, gold_pairs: int) -> str:
    if not coverage:
        return ""
    rows = "".join(
        "<tr>"
        f"<td>{_e(n)}</td>"
        f"<td>{c.get('gold_reachable', 0):,} / {c.get('gold_pairs', gold_pairs):,}</td>"
        f"<td>{c.get('tf_coverage_in_data', 0):,}</td>"
        "</tr>"
        for n, c in coverage.items()
    )
    return (
        '<p class="sub">A method cannot be credited for edges its candidate space '
        "never contained. Reachability describes the chosen candidate universe; "
        "omitted predictions within explicit U still count in metric denominators.</p>"
        '<div class="scroll"><table><thead><tr><th>method</th>'
        "<th>gold reachable in candidates</th><th>regulators overlapping gold</th>"
        "</tr></thead><tbody>" + rows + "</tbody></table></div>"
    )


_ID_VERDICTS = {
    "ok": ("pass", "identifiers agree"),
    "case-mismatch": ("warn", "matched only after folding case"),
    "no-overlap": ("fail", "no shared gene identifiers &mdash; metrics vacuous"),
}


def _identifier_banner(report: dict) -> str:
    """A page-top alarm when the numbers below cannot mean anything.

    A method sharing no gene identifier with the gold still fills every column
    of the verdict table with zeros and dashes, which reads as "performed
    poorly" rather than "was never scored". Say it before the tables, not after.
    """
    ids = report.get("identifiers") or {}
    broken = [n for n, d in ids.items() if d.get("verdict") == "no-overlap"]
    if not broken:
        return ""
    return ('<div class="note" style="border-left-color:#b3261e;background:#fce8e6">'
            "<b>Every metric on this page is vacuous for "
            f"{', '.join(f'<code>{_e(n)}</code>' for n in broken)}.</b> "
            "Those edge lists and the gold standard <b>share no gene identifier "
            "at all</b>, so nothing could match. Usual causes: the wrong "
            "organism's gold, Ensembl IDs scored against gene symbols, or the "
            "wrong gold file. See <i>Gene-identifier matching</i> below.</div>")


def _collision_note(report: dict) -> str:
    """Two edge lists deriving the same method name; one was renamed to score both.

    Kept out of :func:`_identifier_section` on purpose: collisions happen on runs
    that carry no identifier diagnostic at all, and nesting them there made the
    rename silently disappear from those reports.
    """
    collisions = report.get("name_collisions") or []
    if not collisions:
        return ""
    return ('<div class="note"><b>Method names collided &mdash; one was renamed.'
            "</b><ul>"
            + "".join(f"<li><code>{_e(c.get('file', ''))}</code> also derives "
                      f"<code>{_e(c.get('derived', ''))}</code>; renamed to "
                      f"<code>{_e(c.get('renamed_to', ''))}</code></li>"
                      for c in collisions)
            + "</ul></div>")


def _identifier_section(report: dict) -> str:
    ids = report.get("identifiers") or {}
    if not ids:
        return ""
    rows = []
    for name, d in ids.items():
        verdict = d.get("verdict", "ok")
        cls, blurb = _ID_VERDICTS.get(verdict, ("warn", verdict))
        examples = ""
        if d.get("example_edge_genes") or d.get("example_gold_genes"):
            examples = (f"edges: <code>{_e(', '.join(d.get('example_edge_genes') or []))}"
                        f"</code> &middot; gold: "
                        f"<code>{_e(', '.join(d.get('example_gold_genes') or []))}</code>")
        rows.append(
            "<tr>"
            f"<td>{_e(name)}</td>"
            f'<td><span class="badge {cls}">{_e(verdict.upper())}</span></td>'
            f"<td>{d.get('matched_pairs_case_sensitive', 0):,}</td>"
            f"<td>{d.get('matched_pairs_case_insensitive', 0):,}</td>"
            f"<td>{d.get('shared_genes_case_insensitive', 0):,}</td>"
            f'<td style="text-align:left;white-space:normal">{blurb}'
            + (f"<br><span class=\"sub\">{examples}</span>" if examples else "")
            + "</td></tr>"
        )
    mode = ("case-sensitively (<code>--case-sensitive</code>)"
            if report.get("case_sensitive_symbols")
            else "case-insensitively; pass <code>--case-sensitive</code> to refuse")
    return (
        f'<p class="sub">Gene symbols were matched {mode}. Expression matrices '
        "and ground-truth networks routinely disagree on case (BEELINE ships "
        "mESC as <code>Pou5f1</code> and every mouse network as "
        "<code>POU5F1</code>), and exact matching turns that into a silent "
        "zero.</p>"
        + '<div class="scroll"><table><thead><tr><th>method</th><th>verdict</th>'
          "<th>pairs matched exactly</th><th>pairs matched folded</th>"
          "<th>shared genes</th><th>note</th>"
          "</tr></thead><tbody>" + "".join(rows) + "</tbody></table></div>"
    )


def _significance_section(report: dict) -> str:
    """Descriptive annotation removal; never render inferential tiers."""
    sig = report.get("gold_thinning_sensitivity")
    if not sig:
        return ('<p>Legacy annotation-resampling inference omitted. Recompute using '
                '--gold-thinning-sensitivity.</p>') if report.get('significance') else ''
    rows = []
    for key, pr in (sig.get("pairs") or {}).items():
        a, b = key.split("|", 1)
        ci = pr.get("diff_quantiles") or [float("nan")] * 2
        rows.append(
            f"<tr><td>{_e(a)} vs {_e(b)}</td>"
            f"<td>{_num(pr.get('fraction_a_gt_b'), '.2f')}</td>"
            f"<td>{_num(pr.get('diff_mean'), '+.3f')}</td>"
            f"<td>{_num(ci[0], '+.3f')}, {_num(ci[1], '+.3f')}</td></tr>")
    return (
        f"<p>{_e(sig.get('interpretation', 'Descriptive scenarios only.'))}</p>"
        f"<p>{_e(sig.get('resample', ''))}; {sig.get('n_draws', 0)} draws; "
        f"fixed k={sig.get('epr_k', '')}. Metric: {_e(sig.get('metric', ''))}.</p>"
        '<div class="scroll"><table><thead><tr><th>comparison</th>'
        "<th>fraction A&gt;B</th><th>mean diff</th><th>2.5% / 97.5% scenario quantiles</th>"
        "</tr></thead><tbody>" + "".join(rows) + "</tbody></table></div>"
    )


def _qc_section(prep: dict | None) -> str:
    if not prep:
        return ""
    qc = prep.get("qc") or {}
    gates = qc.get("gates") or []
    if not gates:
        return ""
    rows = "".join(
        f'<tr><td>{_e(g["name"])}</td>'
        f'<td><span class="badge {g["status"]}">{g["status"].upper()}</span></td>'
        f'<td style="text-align:left;white-space:normal">{_e(g["detail"])}</td></tr>'
        for g in gates
    )
    warn_list = prep.get("warnings") or []
    warns = ""
    if warn_list:
        warns = ('<div class="note"><b>Preprocessing warnings</b><ul>'
                 + "".join(f"<li>{_e(w)}</li>" for w in warn_list) + "</ul></div>")
    strat = prep.get("strategy_chosen")
    head = (f'<p class="sub">Normalization strategy: <code>{_e(strat)}</code>'
            + (" &middot; <b>run forced past failing gates</b>"
               if prep.get("forced_past_gates") else "") + "</p>")
    return (head + warns
            + '<div class="scroll"><table><thead><tr><th>gate</th><th>status</th>'
              "<th>detail</th></tr></thead><tbody>" + rows + "</tbody></table></div>")


def _provenance(report: dict, prep: dict | None) -> str:
    # Two environments, and they are routinely different machines: the edge
    # lists are produced on a cluster and scored on a laptop weeks later. An
    # `or` chain showed only the scoring one - unlabelled - so a report could
    # state "arboreto: not installed" directly beside the genie3 edge list that
    # arboreto had produced. Show both, and say which is which.
    envs = [("Environment &mdash; evaluation (this report)",
             report.get("environment") or {}),
            ("Environment &mdash; preprocessing (from the prep manifest)",
             (prep or {}).get("environment") or {})]
    blocks = []
    for title, env in envs:
        if not env:
            continue
        rows = "".join(f"<tr><td>{_e(k)}</td>"
                       f'<td style="text-align:left">{_e(v)}</td></tr>'
                       for k, v in env.items())
        blocks.append(f"<h3>{title}</h3><div class=\"scroll\"><table><tbody>"
                      + rows + "</tbody></table></div>")
    files = {n: m.get("file") for n, m in (report.get("methods") or {}).items()
             if m.get("file")}
    if files:
        blocks.append("<h3>Edge lists scored</h3><div class=\"scroll\"><table><tbody>"
                      + "".join(f"<tr><td>{_e(n)}</td>"
                                f'<td style="text-align:left"><code>{_e(f)}</code>'
                                "</td></tr>" for n, f in files.items())
                      + "</tbody></table></div>")
    plugins = report.get("plugins") or []
    if plugins:
        blocks.append("<h3>Plugins registered</h3><div class=\"scroll\"><table>"
                      "<thead><tr><th>name</th><th>origin</th><th>requested as</th>"
                      "</tr></thead><tbody>"
                      + "".join(f"<tr><td>{_e(p.get('name'))}</td>"
                                f"<td>{_e(p.get('origin'))}</td>"
                                f'<td style="text-align:left"><code>'
                                f"{_e(p.get('requested_as'))}</code></td></tr>"
                                for p in plugins)
                      + "</tbody></table></div>")
    if prep and prep.get("steps"):
        steps = "".join(
            f"<li><b>{_e(s.get('step'))}</b>"
            + (f" &rarr; {s['n_cells_after']:,} cells &times; {s['n_genes_after']:,} genes"
               if "n_cells_after" in s else "")
            + (f" &middot; {_e(s.get('detail'))}" if s.get("detail") else "")
            + "</li>"
            for s in prep["steps"])
        blocks.append("<details><summary>Preprocessing steps "
                      f"({len(prep['steps'])})</summary><ul>{steps}</ul></details>")
    return "".join(blocks)


# --------------------------------------------------------------------------- #
# public entry point
# --------------------------------------------------------------------------- #

def _diagnostics_section(document: dict | None) -> str:
    from grnbench.diagnostics import diagnostics_view, format_choices, format_dimensions

    view = diagnostics_view(document)
    if not view["available"]:
        return f'<p class="sub">{_e(view["note"])}</p>'
    rows = []
    for entry in view["entries"]:
        failure = "&ndash;"
        if entry["failure_category"]:
            failure = f"{_e(entry['failure_category'])}: {_e(entry['message'])}"
        rows.append(
            "<tr>"
            f"<td>{_e(entry['stage'])}</td><td>{_e(entry['name'])}</td>"
            f"<td>{_e(entry['context'])}</td><td>{_e(entry['outcome'])}</td>"
            f"<td>{_e(format_dimensions(entry['before']))}</td>"
            f"<td>{_e(format_dimensions(entry['after']))}</td>"
            f"<td>{_e(format_choices(entry['choices']))}</td><td>{failure}</td>"
            "</tr>"
        )
    selected = _e(view["selected_candidate"] or "none")
    usable = ", ".join(_e(item) for item in view["usable_candidates"]) or "none"
    return (
        f'<p class="sub">Overall outcome: <b>{_e(view["outcome"])}</b>; '
        f'selected candidate: <code>{selected}</code>; usable candidates: {usable}.</p>'
        '<div class="scroll"><table><thead><tr><th>stage</th><th>candidate/model</th>'
        '<th>context</th><th>outcome</th><th>before</th><th>after</th><th>choices</th>'
        '<th>failure</th></tr></thead><tbody>' + "".join(rows) + "</tbody></table></div>"
    )

def render_html(report: dict, prep_manifest: dict | None = None,
                diagnostics: dict | None = None) -> str:
    """Render a full evaluation report (see :mod:`grnbench.cli_eval`) as HTML."""
    methods = report.get("methods") or {}
    gold = report.get("gold", "unknown")
    generated = report.get("created_utc") or datetime.now(timezone.utc).isoformat(
        timespec="seconds")

    meta = [
        f"<span><b>Gold standard</b> {_e(gold)}</span>",
        f"<span><b>Gold pairs</b> {report.get('gold_pairs', 0):,}</span>",
        f"<span><b>Methods</b> {len(methods)}</span>",
        f"<span><b>Null model</b> {_e(report.get('null_mode', 'skipped'))}"
        + (f" &times;{report.get('n_nulls')}" if report.get("n_nulls") else "")
        + "</span>",
        f"<span><b>Generated</b> {_e(generated)}</span>",
    ]

    overlap = report.get("overlap") or {}
    overlap_html = _heatmap(overlap)
    if overlap_html:
        overlap_html = (
            f'<p class="sub">Pairwise Jaccard of the top-{overlap.get("k", 0):,} '
            "edges. Two methods with similar scores but near-zero overlap are not "
            "two confirmations of the same biology.</p>" + overlap_html)

    body = f"""<div class="wrap">
<h1>GRN-Bench comparison report</h1>
<p class="sub">preprocess &rarr; infer &rarr; evaluate &mdash; one harness, one edge schema, full provenance.</p>
<div class="meta">{''.join(meta)}</div>
{_identifier_banner(report)}
{_collision_note(report)}

<h2>Verdict</h2>
{_verdict_table(methods)}
{_bars_svg({n: m.get('auprc_ratio') for n, m in methods.items()}, 'AUPRC ratio')}

<h2>Precision&ndash;recall</h2>
{_pr_curves_svg(methods)}

<h2>Against the null</h2>
{_nulls_section(methods)}

<h2>Gold-thinning sensitivity</h2>
{_significance_section(report) or '<p class="sub">Not computed. Pass <code>--gold-thinning-sensitivity</code> with at least two methods for descriptive annotation-removal scenarios.</p>'}

<h2>Do the methods agree?</h2>
{overlap_html or '<p class="sub">Overlap needs at least two scorable methods.</p>'}

<h2>Subsample stability</h2>
{_stability_section(methods)}

<h2>Gold reachability</h2>
{_coverage_section(report.get('coverage') or {}, report.get('gold_pairs', 0))}

<h2>Gene-identifier matching</h2>
{_identifier_section(report) or '<p class="sub">No identifier diagnostic recorded (older report, or scored outside <code>grnbench-eval</code>).</p>'}

<h2>Preprocessing QC</h2>
{_qc_section(prep_manifest) or '<p class="sub">No preprocessing manifest supplied (pass <code>--prep-manifest manifest.json</code>).</p>'}

<h2>Pipeline diagnostics</h2>
{_diagnostics_section(diagnostics)}

<h2>Provenance</h2>
{_provenance(report, prep_manifest)}

<footer>Generated by GRN-Bench. Reference annotations may be incomplete and biased.
Scores and conditional null diagnostics do not establish biological truth.
Comparisons require a shared candidate universe and a reference appropriate to the scientific task.</footer>
</div>"""

    return ("<!doctype html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">"
            "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">"
            f"<title>GRN-Bench &mdash; {_e(gold)}</title><style>{_CSS}</style>"
            f"</head><body>{body}</body></html>\n")
