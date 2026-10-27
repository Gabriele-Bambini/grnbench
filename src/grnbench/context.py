"""Dataset context: what the data are, read from every place that says so.

Sources, in the order they are consulted:

* the **matrix** itself (gene identifiers -> organism, Ensembl prefixes);
* the AnnData **metadata** (``uns`` text, CELLxGENE-style ``obs`` columns
  such as ``assay``/``organism``/``tissue``, batch/donor/cell-type columns);
* **README-like files** shipped next to the data (README*, *.md, *.txt, ...);
* **download pages** (GEO, Zenodo, CELLxGENE, any URL the user gives).

Extraction is rule-based and deterministic on purpose: the same inputs always
give the same context, and every fact keeps the source and the text snippet
it came from so the user can check it. Facts are weighted votes; when sources
disagree the conflict is reported instead of silently resolved.
"""
from __future__ import annotations

import html
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from urllib.parse import parse_qs, urlparse
from urllib.request import Request, urlopen

import pandas as pd

MAX_TEXT_BYTES = 1_000_000
MAX_README_FILES = 12
README_PATTERNS = ("README*", "readme*", "*.md", "*.rst", "*.txt", "*.json", "*.yaml", "*.yml")
_SKIP_NAMES = {"desktop.ini"}

# (regex, value, weight). Species binomials, taxon ids and genome builds are
# strong evidence; the bare words "human"/"mouse" also appear in passing
# ("human disease", "mouse model of"), so they vote weakly.
_ORGANISM_RULES = [
    (r"\bhomo sapiens\b|ncbitaxon:9606|\bhg19\b|\bhg38\b|\bgrch3[78]\b", "human", 3),
    (r"\bhuman\b", "human", 1),
    (r"\bmus musculus\b|ncbitaxon:10090|\bmm10\b|\bmm39\b|\bgrcm3[89]\b", "mouse", 3),
    (r"\bmouse\b|\bmice\b|\bmurine\b", "mouse", 1),
]
# (regex, platform, technology). "umi" = molecule counts from a droplet/UMI
# protocol; "full_length" = read counts from plate-based full-length protocols.
_PLATFORM_RULES = [
    (r"10x genomics|\bchromium\b|\b10x\b", "10x Chromium", "umi"),
    (r"drop-?seq", "Drop-seq", "umi"),
    (r"\bindrops?\b", "inDrop", "umi"),
    (r"\bcel-?seq2?\b", "CEL-seq", "umi"),
    (r"\bmars-?seq2?\b", "MARS-seq", "umi"),
    (r"\bseq-?well\b", "Seq-Well", "umi"),
    (r"\bsplit-?seq\b", "SPLiT-seq", "umi"),
    (r"\bsci-rna-seq3?\b", "sci-RNA-seq", "umi"),
    (r"smart-?seq ?3", "Smart-seq3", "full_length"),
    (r"smart-?seq ?2|smart-?seq\b|\bsmarter\b", "Smart-seq2", "full_length"),
    (r"fluidigm c1", "Fluidigm C1", "full_length"),
]
_TECHNOLOGY_RULES = [
    (r"\bumis?\b|unique molecular identifier", "umi"),
    (r"full[- ]length", "full_length"),
]
_VALUE_TYPE_RULES = [
    (r"raw (?:umi )?counts?|unnormali[sz]ed|umi counts|count matri(?:x|ces)|read counts", "raw_counts"),
    (r"\btpm\b|\bfpkm\b|\brpkm\b|log[- ]?normali[sz]ed|\blog2?\s*\(|normali[sz]ed (?:expression|data|counts|values)",
     "normalized"),
]
_ACCESSION_RULES = [
    r"\b(GSE\d{3,})\b", r"\b(E-MTAB-\d+)\b", r"\b(PRJNA\d+)\b", r"\b(SRP\d{4,})\b",
    r"\b(10\.5281/zenodo\.\d+)\b",
]
_FLAG_RULES = {
    "time_course": r"time[- ]?course|pseudo-?time|differentiation|time[- ]points?|trajector",
    "perturbation": r"perturb-?seq|\bcrispr\b|knock-?out|knock-?down|\bsgrnas?\b",
}
_TISSUE_RE = re.compile(r"\btissue\s*[:=]\s*([A-Za-z][^\n;,|\t]{1,60})", re.IGNORECASE)
_N_CELLS_RE = re.compile(r"\b(\d{1,3}(?:,\d{3})+|\d{3,7})\s+(?:single[- ])?(?:cells|nuclei)\b",
                         re.IGNORECASE)

_CELL_TYPE_COL_RE = re.compile(r"cell_?type|celltype|cell_ontology|annotation|cluster|leiden|louvain",
                               re.IGNORECASE)
_DONOR_COL_RE = re.compile(r"donor|patient|individual|subject", re.IGNORECASE)
# CELLxGENE-schema obs columns whose (majority) value is descriptive text.
_TEXT_OBS_COLS = ("assay", "organism", "tissue", "suspension_type", "development_stage")


@dataclass
class Fact:
    key: str
    value: str
    source: str
    snippet: str = ""
    weight: int = 1


@dataclass
class DatasetContext:
    facts: list[Fact] = field(default_factory=list)
    sources_read: list[str] = field(default_factory=list)
    sources_failed: list[dict] = field(default_factory=list)
    # The label as written for a value stored in lower case ({(key, value): label}), for a reader that needs the
    # case: the high-mito tissue matcher splits camelCase words (#161). The label kept is the one seen with the
    # highest weight (its strongest evidence), the first seen on a tie, so uns read before obs cannot hide the
    # curated obs label. Not in ``to_dict``.
    spellings: dict = field(default_factory=dict)
    # {(key, value): {label: highest weight it was seen with}}, labels in the order first seen.
    spelling_weights: dict = field(default_factory=dict)

    def add(self, key, value, source, snippet="", weight=1) -> None:
        self.facts.append(Fact(key, str(value), source, snippet.strip(), int(weight)))

    def votes(self, key: str) -> dict[str, int]:
        out: dict[str, int] = {}
        for f in self.facts:
            if f.key == key:
                out[f.value] = out.get(f.value, 0) + f.weight
        return out

    def get(self, key: str) -> str | None:
        """Highest-weighted value; ties go to the value seen first."""
        votes = self.votes(key)
        return max(votes, key=votes.get) if votes else None

    def all(self, key: str) -> list[str]:
        return list(self.votes(key))

    def evidence(self, key: str) -> list[Fact]:
        return [f for f in self.facts if f.key == key]

    def conflicts(self) -> dict[str, dict[str, int]]:
        single = {"organism", "technology", "value_type"}
        return {k: v for k in single if len(v := self.votes(k)) > 1}

    def spelling(self, key: str, value: str | None) -> str:
        """The label as written behind a value stored in lower case, else the value itself ("" for None)."""
        return self.spellings.get((key, value), value or "")

    def add_lower(self, key, label, source, snippet="", weight=1) -> None:
        """``add`` the label in lower case (as stored and shown) and remember how it was written: the spelling with
        the highest weight, the first seen on a tie."""
        label, weight = str(label), int(weight)
        self.add(key, label.lower(), source, snippet, weight)
        seen = self.spelling_weights.setdefault((key, label.lower()), {})
        seen[label] = max(seen.get(label, weight), weight)
        self.spellings[(key, label.lower())] = max(seen, key=seen.get)   # max keeps the first of equal weights

    def summary(self) -> dict:
        keys = ("organism", "technology", "platform", "value_type", "tissue",
                "cell_type_key", "batch_key", "donor_key", "time_course",
                "perturbation", "n_cells_declared")
        out = {k: self.get(k) for k in keys}
        out["accessions"] = self.all("accession")
        return out

    def to_dict(self) -> dict:
        return {"summary": self.summary(), "conflicts": self.conflicts(),
                "facts": [asdict(f) for f in self.facts],
                "sources_read": list(self.sources_read),
                "sources_failed": list(self.sources_failed)}


def _snippet(text: str, m: re.Match, pad: int = 60) -> str:
    a, b = max(0, m.start() - pad), min(len(text), m.end() + pad)
    return " ".join(text[a:b].split())


def extract_from_text(text: str, source: str, ctx: DatasetContext | None = None) -> DatasetContext:
    """Add every rule match in ``text`` to ``ctx`` (created when None)."""
    ctx = ctx if ctx is not None else DatasetContext()
    if not text:
        return ctx
    for pattern, value, weight in _ORGANISM_RULES:
        for m in re.finditer(pattern, text, re.I):
            ctx.add("organism", value, source, _snippet(text, m), weight)
    for pattern, platform, tech in _PLATFORM_RULES:
        m = re.search(pattern, text, re.I)
        if m:
            ctx.add("platform", platform, source, _snippet(text, m), 2)
            ctx.add("technology", tech, source, _snippet(text, m), 2)
    for pattern, tech in _TECHNOLOGY_RULES:
        m = re.search(pattern, text, re.I)
        if m:
            ctx.add("technology", tech, source, _snippet(text, m), 1)
    for pattern, value in _VALUE_TYPE_RULES:
        m = re.search(pattern, text, re.I)
        if m:
            ctx.add("value_type", value, source, _snippet(text, m), 1)
    seen = set()
    for pattern in _ACCESSION_RULES:
        for m in re.finditer(pattern, text, re.IGNORECASE):
            acc = m.group(1).upper() if not m.group(1).startswith("10.") else m.group(1)
            if acc not in seen:
                seen.add(acc)
                ctx.add("accession", acc, source, _snippet(text, m))
    for key, pattern in _FLAG_RULES.items():
        m = re.search(pattern, text, re.I)
        if m:
            ctx.add(key, "yes", source, _snippet(text, m))
    for m in _TISSUE_RE.finditer(text):
        ctx.add_lower("tissue", m.group(1).strip(), source, _snippet(text, m), 2)
    for m in _N_CELLS_RE.finditer(text):
        ctx.add("n_cells_declared", m.group(1).replace(",", ""), source, _snippet(text, m))
    return ctx


def species_label(text: str) -> str:
    """'human' / 'mouse' for those two, else the lower-case binomial (e.g. 'macaca mulatta')."""
    low = " ".join(str(text).lower().split())
    if "homo sapiens" in low or low in {"human", "ncbitaxon:9606"}:
        return "human"
    if "mus musculus" in low or low in {"mouse", "ncbitaxon:10090"}:
        return "mouse"
    return low


def extract_from_adata(adata, ctx: DatasetContext | None = None) -> DatasetContext:
    """Facts carried by the matrix and its AnnData metadata."""
    from grnbench.hygiene import detect_batch_column, organism_hint

    ctx = ctx if ctx is not None else DatasetContext()
    names = pd.Index(adata.var_names).astype(str)
    ens_human = names.str.startswith("ENSG").mean()
    ens_mouse = names.str.startswith("ENSMUSG").mean()
    if ens_human > 0.5:
        ctx.add("organism", "human", "matrix", f"{ens_human:.0%} of gene ids are ENSG*", 10)
    elif ens_mouse > 0.5:
        ctx.add("organism", "mouse", "matrix", f"{ens_mouse:.0%} of gene ids are ENSMUSG*", 10)
    else:
        hint = organism_hint(names)
        if hint in ("human-like", "mouse-like"):
            ctx.add("organism", hint.split("-")[0], "matrix",
                    f"gene symbol casing looks {hint}", 4)
    ctx.sources_read.append("matrix")

    uns_text = []
    for key, value in getattr(adata, "uns", {}).items():
        if isinstance(value, str) and len(value) < 20_000:
            uns_text.append(f"{key}: {value}")
    if uns_text:
        extract_from_text("\n".join(uns_text), "metadata:uns", ctx)
        ctx.sources_read.append("metadata:uns")

    obs = adata.obs
    for col in _TEXT_OBS_COLS:
        if col in obs.columns and len(obs):
            top = obs[col].astype(str).value_counts()
            value = str(top.index[0])
            if col == "organism":
                # Curated species label: strongest evidence, and it names species other than
                # human/mouse (macaque, marmoset, zebrafish...) that gene casing cannot.
                ctx.add("organism", species_label(value), f"metadata:obs['{col}']",
                        f"majority value of obs['{col}']: {value}", 12)
                continue
            extract_from_text(f"{col}: {value}", f"metadata:obs['{col}']", ctx)
            if col == "tissue":
                ctx.add_lower("tissue", value, f"metadata:obs['{col}']",
                        f"majority value of obs['{col}']", 3)
    batch, how = detect_batch_column(adata)
    if batch is not None:
        ctx.add("batch_key", batch, "metadata:obs", f"column '{batch}' ({how})", 2)
    for col in obs.columns:
        n = obs[col].nunique()
        if _DONOR_COL_RE.search(str(col)) and 2 <= n <= 500:
            ctx.add("donor_key", col, "metadata:obs", f"column '{col}' with {n} levels", 2)
        if _CELL_TYPE_COL_RE.search(str(col)) and 2 <= n <= 300:
            # Curated annotation columns outrank cluster assignments.
            low = str(col).lower()
            weight = (4 if low in {"cell_type", "celltype", "cell_ontology_class", "cell_type_label"}
                      else 3 if "type" in low or "annotation" in low or "ontology" in low else 1)
            ctx.add("cell_type_key", col, "metadata:obs", f"column '{col}' with {n} levels", weight)
    return ctx


def _looks_like_table(head: str) -> bool:
    lines = [ln for ln in head.splitlines()[:5] if ln.strip()]
    return bool(lines) and all(max(ln.count("\t"), ln.count(",")) >= 20 for ln in lines)


def discover_text_files(data_path: str | Path) -> list[Path]:
    """README-like files next to the data (and one level up for 10x folders)."""
    p = Path(data_path)
    folders = [p] if p.is_dir() else [p.parent]
    folders.append(folders[0].parent)
    found: list[Path] = []
    for folder in folders:
        if not folder.is_dir():
            continue
        for pattern in README_PATTERNS:
            for f in sorted(folder.glob(pattern)):
                if (f.is_file() and f.name.lower() not in _SKIP_NAMES and f != p
                        and f not in found and f.stat().st_size <= MAX_TEXT_BYTES):
                    found.append(f)
        if found:
            break  # the data folder's own README wins over the parent's
    return found[:MAX_README_FILES]


def read_text_file(path: str | Path) -> str | None:
    """Text content, or None for tables/binary files that only look like text."""
    raw = Path(path).read_bytes()[:MAX_TEXT_BYTES]
    if b"\x00" in raw[:4096]:
        return None
    text = raw.decode("utf-8", errors="replace")
    return None if _looks_like_table(text[:4096]) else text


def html_to_text(page: str) -> str:
    page = re.sub(r"(?is)<(script|style|noscript)\b.*?</\1>", " ", page)
    page = re.sub(r"(?i)<br\s*/?>|</(p|div|tr|li|h\d)>", "\n", page)
    page = re.sub(r"(?s)<[^>]+>", " ", page)
    return re.sub(r"[ \t\r\f\v]+", " ", html.unescape(page).replace("\xa0", " "))


def geo_text_url(url_or_accession: str) -> str | None:
    """The plain-text (SOFT, brief) GEO record for a GSE/GSM URL or accession."""
    s = url_or_accession.strip()
    acc = None
    if re.fullmatch(r"(?i)G(SE|SM)\d+", s):
        acc = s.upper()
    else:
        parsed = urlparse(s)
        if "ncbi.nlm.nih.gov" in parsed.netloc and "/geo/" in parsed.path:
            acc = (parse_qs(parsed.query).get("acc") or [None])[0]
    if not acc:
        return None
    return (f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}"
            "&targ=self&form=text&view=brief")


def fetch_url_text(url: str, timeout: float = 20.0, opener=urlopen) -> str:
    target = geo_text_url(url) or url
    req = Request(target, headers={"User-Agent": "grnbench-context/1.0"})
    with opener(req, timeout=timeout) as resp:
        raw = resp.read(MAX_TEXT_BYTES)
        ctype = resp.headers.get("Content-Type", "") if getattr(resp, "headers", None) else ""
    text = raw.decode("utf-8", errors="replace")
    return html_to_text(text) if "html" in ctype or text.lstrip()[:1] == "<" else text


def build_context(data_path: str | Path | None = None, adata=None, *,
                  text_files=None, urls=(), notes: str | None = None,
                  discover: bool = True, fetch=fetch_url_text) -> DatasetContext:
    """Collect facts from the matrix, metadata, README files, URLs and notes.

    Network or file problems never raise: they are listed in
    ``sources_failed`` with the reason, and the context is built from what
    could be read.
    """
    ctx = DatasetContext()
    if adata is not None:
        extract_from_adata(adata, ctx)
    files = list(text_files or [])
    if discover and data_path is not None:
        files += [f for f in discover_text_files(data_path) if f not in files]
    for f in files:
        try:
            text = read_text_file(f)
        except OSError as exc:
            ctx.sources_failed.append({"source": str(f), "reason": str(exc)})
            continue
        if text is None:
            continue
        extract_from_text(text, f"readme:{Path(f).name}", ctx)
        ctx.sources_read.append(f"readme:{Path(f).name}")
    for url in urls:
        url = str(url).strip()
        if not url:
            continue
        try:
            text = fetch(url)
        except Exception as exc:  # noqa: BLE001 - offline/HTTP errors are data here
            ctx.sources_failed.append({"source": url, "reason": f"{type(exc).__name__}: {exc}"})
            continue
        extract_from_text(text, f"url:{url}", ctx)
        ctx.sources_read.append(f"url:{url}")
    if notes:
        extract_from_text(notes, "user-notes", ctx)
        ctx.sources_read.append("user-notes")
    return ctx
