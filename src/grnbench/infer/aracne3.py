"""ARACNe3 adapter (Java jar) with honest availability probing.

ARACNe3 ships as a jar (github.com/califano-lab/ARACNe3) and needs ``java`` on
PATH. The wrapper never fakes availability: without java + the ``ARACNE3_JAR``
env var the method is reported unavailable and running it raises
:class:`MethodError` with the exact missing piece. The output parser targets
the published ARACNe3 network table format (tab-separated; regulator/target/
score columns, optional sign column) and is unit-tested against a fixture, so
only the subprocess call itself remains unverified until a jar is available.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd


def aracne3_availability() -> tuple[bool, str]:
    java = shutil.which("java")
    if java is None:
        return False, "java executable not found on PATH"
    jar = os.environ.get("ARACNE3_JAR", "")
    if not jar or not Path(jar).exists():
        return False, ("ARACNE3_JAR env var not set (point it at the aracne3.jar "
                       "build); download from github.com/califano-lab/ARACNe3")
    return True, ""


def run_aracne3(expr: pd.DataFrame,
                tf_names: list[str] | None = None,
                threads: int = 4,
                timeout_s: int = 3600,
                extra_args: list[str] | None = None,
                **_) -> pd.DataFrame:
    """Write the expression matrix to disk, invoke the jar, parse its network.

    Expression is written samples-as-rows / genes-as-columns (ARACNe3's
    documented input layout). A TF whitelist, when provided, is passed through
    via the jar's regulator-file argument when supported by the build; here it
    filters post hoc to keep the adapter independent of flag drift across
    releases.
    """
    from grnbench.infer.base import MethodError

    ok, why = aracne3_availability()
    if not ok:
        raise MethodError(f"aracne3 unavailable: {why}")

    with tempfile.TemporaryDirectory(prefix="grnbench_aracne3_") as td:
        inp = Path(td) / "expr.tsv"
        outp = Path(td) / "network.txt"
        expr.to_csv(inp, sep="\t", float_format="%.6f")

        cmd = [
            shutil.which("java"), "-Xmx4g", "-jar", os.environ["ARACNE3_JAR"],
            "-i", str(inp), "-o", str(outp),
            "-t", str(threads),
            *(list(extra_args) if extra_args else []),
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              timeout=int(timeout_s))  # noqa: S603 - fixed argv list
        if proc.returncode != 0:
            tail = (proc.stderr or proc.stdout or "").strip()[-800:]
            raise MethodError(f"aracne3 failed (rc={proc.returncode}): {tail}")
        if not outp.exists():
            raise MethodError(
                "aracne3 produced no network file; inspect stdout/stderr wiring "
                "for this jar release"
            )
        net = parse_aracne3_text(outp.read_text(encoding="utf-8", errors="replace"))

    if tf_names:
        allowed = {str(t).upper() for t in tf_names}
        net = net[net.iloc[:, 0].str.upper().isin(allowed)]
    return net


_HEADER_ALIASES = {
    "source": {"regulator", "tf", "source", "tf1", "gene1", "transcription_factor"},
    "target": {"target", "target_gene", "gene2"},
    "weight": {"score", "mi", "weight", "importance", "value", "mipscore"},
    "direction": {"sign", "direction", "mode", "activation", "effect"},
}


def _canonical(col: str) -> str | None:
    c = col.strip().lower()
    for canon, aliases in _HEADER_ALIASES.items():
        if c in aliases or c.startswith(canon):
            return canon
    return None


def parse_aracne3_text(text: str) -> pd.DataFrame:
    """Parse an ARACNe3-style network table into [source, target, weight(, direction)].

    Tolerant of column-name drift across jar releases: header names are matched
    case-insensitively against known aliases; sign columns accept +/-, 1/-1 and
    activation/inhibition tokens.
    """
    lines = [ln.rstrip("\n") for ln in text.splitlines()
             if ln.strip() and not ln.lstrip().startswith("#")]
    if not lines:
        raise ValueError("empty ARACNe3 output")
    header = lines[0].split("\t") if "\t" in lines[0] else lines[0].split()
    mapping = {}
    for i, col in enumerate(header):
        canon = _canonical(col)
        if canon and canon not in mapping:
            mapping[canon] = i
    rows = []
    for ln in lines[1:]:
        parts = ln.split("\t") if "\t" in ln else ln.split()
        try:
            row = {
                "source": parts[mapping.get("source", 0)],
                "target": parts[mapping.get("target", 1)],
                "weight": float(parts[mapping.get("weight", 2)]),
            }
        except (IndexError, ValueError):
            continue  # comment/footer line
        if "direction" in mapping:
            tok = parts[mapping["direction"]].strip().lower()
            row["direction"] = (-1 if tok in {"-1", "-", "inhibition", "inhibits", "repress"}
                                else 1)
        rows.append(row)
    if not rows:
        raise ValueError("no data rows parsed from ARACNe3 output")
    return pd.DataFrame(rows)
