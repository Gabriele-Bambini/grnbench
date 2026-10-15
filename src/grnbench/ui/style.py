"""Look and feel: one stylesheet, the page header and the brand mark.

Cosmetic only. Nothing in the interface depends on this CSS to work: it styles Streamlit's stable
``data-testid`` hooks and the ``st-key-*`` classes of the widgets this app names, and every colour is either
the accent or a translucent overlay of the page's own text colour, so it holds in a light and a dark browser.
"""
from __future__ import annotations

from grnbench.ui import steps, theme

CSS = """
:root {
  --gb-accent: __ACCENT__;            /* fills: badges, the brand mark (white text on it) */
  --gb-accent-text: __ACCENT__;       /* the accent as text or outline; lighter on a dark page, see below */
  --gb-good: #15803D;
  --gb-line: rgba(128, 128, 128, 0.30);
  --gb-wash: rgba(128, 128, 128, 0.07);
}
/* Streamlit sets color-scheme on the app from the active theme, so light-dark() follows its light or dark variant. */
@supports (color: light-dark(red, blue)) {
  :root { --gb-accent-text: light-dark(__ACCENT__, __ACCENT_ON_DARK__); --gb-good: light-dark(#15803D, #4ADE80); }
}
:root { --gb-accent-soft: color-mix(in srgb, var(--gb-accent-text) 14%, transparent); }

/* ---- frame ------------------------------------------------------------------------------- */
[data-testid="stMainBlockContainer"] { max-width: 1180px; padding: 4.6rem 2.4rem 5rem; }
@media (max-width: 800px) { [data-testid="stMainBlockContainer"] { padding: 4.2rem 1rem 4rem; } }
.stAppDeployButton, [data-testid="stAppDeployButton"] { display: none; }
[data-testid="stMainBlockContainer"] h2 { margin-top: .4rem; }
[data-testid="stMainBlockContainer"] h3 { margin-top: .2rem; }
[data-testid="stCaptionContainer"] { line-height: 1.5; color: inherit; opacity: .78; }   /* a step darker than the theme's 60% grey */
[data-testid="stCaptionContainer"] p { color: inherit; }

/* ---- page header ------------------------------------------------------------------------- */
.gb-head { margin: 0 0 1.5rem; }
.gb-eyebrow { font-size: .72rem; font-weight: 700; letter-spacing: .1em; text-transform: uppercase;
  color: var(--gb-accent-text); margin: 0 0 .4rem; display: flex; gap: .55rem; align-items: center; }
.gb-opt { border: 1px solid var(--gb-line); border-radius: 999px; padding: .05rem .5rem; font-size: .66rem;
  letter-spacing: .06em; color: inherit; opacity: .75; }
.gb-title { font-size: 2rem; line-height: 1.18; font-weight: 700; letter-spacing: -.012em; margin: 0 0 .55rem; }
.gb-lede { font-size: 1.04rem; line-height: 1.6; opacity: .78; max-width: 60rem; margin: 0; }

/* ---- sidebar ----------------------------------------------------------------------------- */
[data-testid="stSidebar"] [data-testid="stSidebarContent"] { padding: 0 .9rem 1.2rem; }
[data-testid="stSidebarHeader"] { height: 2.2rem; }
.gb-brand { display: flex; align-items: center; gap: .65rem; margin: .1rem 0 1rem .15rem; }
.gb-brand svg { width: 34px; height: 34px; flex: none; }
.gb-brand b { display: block; font-size: 1.2rem; font-weight: 700; letter-spacing: -.015em; line-height: 1.15; }
.gb-brand small { display: block; font-size: .75rem; opacity: .7; line-height: 1.3; margin-top: .1rem;
  text-rendering: optimizeLegibility; }
.gb-group { font-size: .66rem; font-weight: 700; letter-spacing: .12em; text-transform: uppercase; opacity: .5;
  margin: 1.15rem 0 .25rem .55rem; }
[class*="st-key-nav-"] { margin: 0; }
[class*="st-key-nav-"] button { width: 100%; min-height: 0; padding: .4rem .6rem; border: 0; border-radius: .55rem;
  background: transparent; box-shadow: none; justify-content: flex-start; }
[class*="st-key-nav-"] button > div { width: 100%; justify-content: flex-start; }
[class*="st-key-nav-"] button p { font-size: .92rem; font-weight: 500; text-align: left; margin: 0; }
[class*="st-key-nav-"] button:hover { background: var(--gb-wash); }
[class*="st-key-nav-"] [data-testid="stIconMaterial"] { font-size: 1.15rem; margin-right: .2rem; }
[class*="st-key-nav-todo-"] button { opacity: .74; }
[class*="st-key-nav-done-"] [data-testid="stIconMaterial"] { color: var(--gb-good); }
[class*="st-key-nav-cur-"] button { background: var(--gb-accent-soft); }
[class*="st-key-nav-cur-"] button p { font-weight: 650; }
[class*="st-key-nav-cur-"] [data-testid="stIconMaterial"] { color: var(--gb-accent-text); }
.gb-side-note { font-size: .78rem; opacity: .62; line-height: 1.4; margin: .2rem .3rem 0; }

/* ---- cards and blocks -------------------------------------------------------------------- */
[data-testid="stVerticalBlockBorderWrapper"] { border-radius: .8rem; }
[class*="st-key-card"] { padding: .2rem .35rem .35rem; }
.gb-card-title { display: flex; align-items: center; gap: .6rem; font-size: 1.04rem; font-weight: 650; margin: .15rem 0 .1rem; }
.gb-num { display: inline-flex; align-items: center; justify-content: center; width: 1.55rem; height: 1.55rem;
  border-radius: 999px; background: var(--gb-accent); color: #fff; font-size: .8rem; font-weight: 700; flex: none; }
.gb-card-sub { font-size: .9rem; opacity: .7; margin: 0 0 .5rem 2.15rem; line-height: 1.45; }
.st-key-empty-state { padding: 1.8rem .6rem 1.4rem; }
.st-key-empty-state [data-testid="stVerticalBlock"] { gap: .55rem; }
.st-key-empty-state [data-testid="stMarkdownContainer"] { text-align: center; }
.st-key-empty-state [data-testid="stMarkdownContainer"] h2 { opacity: .45; margin: 0; padding: 0; }
.st-key-empty-state [data-testid="stMarkdownContainer"] h3 { margin: 0 0 .1rem; padding: 0; }
.st-key-empty-state [data-testid="stMarkdownContainer"] p { opacity: .78; max-width: 36rem; margin: 0 auto; line-height: 1.55; }

/* ---- radios that choose one of a few options look like a segmented control ---------------- */
:is([class*="st-key-"][class*="_source"], [class*="st-key-"][class*="_preset"], .st-key-explore_stage,
    .st-key-explore_method, .st-key-bench_metric) [data-testid="stRadioGroup"] { gap: .4rem; flex-wrap: wrap; }
:is([class*="st-key-"][class*="_source"], [class*="st-key-"][class*="_preset"], .st-key-explore_stage,
    .st-key-explore_method, .st-key-bench_metric) label[data-testid="stRadioOption"] {
  border: 1px solid var(--gb-line); border-radius: 999px; padding: .28rem .85rem; margin: 0; }
:is([class*="st-key-"][class*="_source"], [class*="st-key-"][class*="_preset"], .st-key-explore_stage,
    .st-key-explore_method, .st-key-bench_metric) label[data-testid="stRadioOption"] > div > div:first-child { display: none; }
:is([class*="st-key-"][class*="_source"], [class*="st-key-"][class*="_preset"], .st-key-explore_stage,
    .st-key-explore_method, .st-key-bench_metric) label[data-testid="stRadioOption"]:hover { background: var(--gb-wash); }
:is([class*="st-key-"][class*="_source"], [class*="st-key-"][class*="_preset"], .st-key-explore_stage,
    .st-key-explore_method, .st-key-bench_metric) label[data-testid="stRadioOption"][data-selected="true"] {
  background: var(--gb-accent-soft); border-color: var(--gb-accent-text); }
:is([class*="st-key-"][class*="_source"], [class*="st-key-"][class*="_preset"], .st-key-explore_stage,
    .st-key-explore_method, .st-key-bench_metric) label[data-selected="true"] p { color: var(--gb-accent-text); font-weight: 600; }

/* ---- native elements --------------------------------------------------------------------- */
[data-testid="stAlertContainer"] { border-radius: .7rem; }
[data-testid="stMetric"] { padding: .15rem .1rem; }
[data-testid="stMetricLabel"] p { font-size: .82rem; opacity: .75; }
[data-testid="stMetricValue"] { font-weight: 650; letter-spacing: -.01em; font-size: 1.75rem; }
[data-testid="stMetricValue"] > div { white-space: normal; overflow: visible; text-overflow: clip; line-height: 1.2; }
[data-testid="stExpander"] { border-radius: .7rem; }
[data-testid="stTabs"] [role="tab"] p { font-size: .96rem; font-weight: 550; }
[data-testid="stTabs"] [role="tab"][aria-selected="true"] p, [data-testid="stSliderThumbValue"] { color: var(--gb-accent-text); }
[data-testid="stFileUploaderDropzone"] { border-radius: .7rem; }
button[kind="primary"] { font-weight: 600; }
/* the theme darkens a primary button to near-black on hover; a step is enough */
[data-testid="stBaseButton-primary"]:hover:not(:disabled) { background-color: color-mix(in srgb, var(--gb-accent) 84%, #000);
  border-color: color-mix(in srgb, var(--gb-accent) 84%, #000); }
/* a row of headline numbers keeps two per line on a phone instead of one card per row */
@media (max-width: 640px) {
  [data-testid="stHorizontalBlock"]:has(> [data-testid="stColumn"] [data-testid="stMetric"]) { flex-direction: row; flex-wrap: wrap; gap: .6rem; }
  [data-testid="stHorizontalBlock"]:has(> [data-testid="stColumn"] [data-testid="stMetric"]) > [data-testid="stColumn"] {
    flex: 1 1 calc(50% - .6rem) !important; min-width: calc(50% - .6rem) !important; width: auto !important; }
  [data-testid="stMetricValue"] { font-size: 1.35rem; }
}
.gb-foot { border-top: 1px solid var(--gb-line); margin: 2.2rem 0 .9rem; }
.gb-hint { font-size: .86rem; opacity: .68; }
"""

BRAND_SVG = (
    '<svg viewBox="0 0 32 32" aria-hidden="true"><rect width="32" height="32" rx="9" fill="__ACCENT__"/>'
    '<g stroke="#fff" stroke-width="1.7" stroke-linecap="round" opacity=".92">'
    '<path d="M10 10.5 22 16M10 21.5 22 16"/></g>'
    '<g fill="#fff"><circle cx="10" cy="10.5" r="3.1"/><circle cx="10" cy="21.5" r="3.1"/>'
    '<circle cx="22" cy="16" r="3.6"/></g></svg>'
)


def stylesheet() -> str:
    return "<style>" + CSS.replace("__ACCENT_ON_DARK__", theme.ACCENT_ON_DARK).replace("__ACCENT__", theme.ACCENT) + "</style>"


def brand_html() -> str:
    return ('<div class="gb-brand">' + BRAND_SVG.replace("__ACCENT__", theme.ACCENT) +
            '<div><b>grnbench</b><small>Context-aware cleaning for<br>gene regulatory network inference</small>'
            '</div></div>')


def group_html(label: str) -> str:
    return f'<div class="gb-group">{label}</div>'


def header_html(step_key: str) -> str:
    """Eyebrow, title and lede of a screen; the stylesheet rides along so it is loaded before anything is drawn."""
    s = steps.BY_KEY[step_key]
    pos = steps.position(step_key)
    eyebrow = f"Step {pos[0]} of {pos[1]} · {s.short}" if pos else (s.group if step_key == steps.QUICK else f"{s.group} · {s.short}")
    opt = '<span class="gb-opt">Optional</span>' if s.optional else ""
    return (stylesheet() + f'<div class="gb-head"><div class="gb-eyebrow">{eyebrow}{opt}</div>'
            f'<h1 class="gb-title">{s.title}</h1><p class="gb-lede">{s.lede}</p></div>')


def card_title_html(number: int | None, title: str, sub: str = "") -> str:
    """Heading of a card: a numbered badge, a title and, below, one line saying what the card is for."""
    badge = f'<span class="gb-num">{number}</span>' if number is not None else ""
    out = f'<div class="gb-card-title">{badge}<span>{title}</span></div>'
    if sub:
        left = "2.15rem" if number is not None else "0"
        out += f'<div class="gb-card-sub" style="margin-left:{left}">{sub}</div>'
    return out

