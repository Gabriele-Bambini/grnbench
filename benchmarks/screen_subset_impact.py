"""#75: do the C4-C6 motif lifts (auto arm) move when the datasets whose auto arm is the
one-cell-type `context_scope=selected_type` candidate are dropped or given the baseline arm's values?
Uses the committed rows and real_bench.models_section (the registered comparison code)."""
import json, sys
from pathlib import Path
import pandas as pd

repo = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(repo / "benchmarks"))
import real_bench

SEL = "context_scope=selected_type"
for exp in ("EXP-20260925-LM-2", "EXP-20260926-LM-1", "EXP-20261002-LM-1"):
    rows = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((repo / "benchmarks/results" / exp / "rows").glob("*.json"))]
    df = pd.DataFrame(rows)
    hit = df["auto_candidate"].astype(str) == SEL
    print(f"\n===== {exp}: {len(df)} rows, auto = {SEL} in {int(hit.sum())}: "
          f"auto_n_contexts {df.loc[hit, 'auto_n_contexts'].tolist()} vs rules_n_contexts {df.loc[hit, 'rules_n_contexts'].tolist()}")
    swapped = df.copy()
    for c in [c for c in df if c.startswith("m_") and "_auto_" in c]:
        b = c.replace("_auto_", "_baseline_")
        if b in df:
            swapped.loc[hit, c] = df.loc[hit, b]
    for label, d in (("registered (all rows)", df), ("affected rows dropped", df[~hit]),
                     ("affected rows -> baseline arm values", swapped)):
        lines = real_bench.models_section(d.reset_index(drop=True))
        # the AUPRC-ratio block, auto arm, declared comparison lines
        start = next(i for i, l in enumerate(lines) if l.startswith("### AUPRC"))
        stop = next(i for i, l in enumerate(lines) if l.startswith("### EPR"))
        block = lines[start:stop]
        auto_i = next(i for i, l in enumerate(block) if "auto" in l.lower() and l.startswith("####"))
        base_i = next((i for i, l in enumerate(block) if "baseline" in l.lower() and l.startswith("####")), len(block))
        print(f"-- {label}")
        for l in block[auto_i:base_i]:
            if "_motif |" in l or "| motif_prior |" in l:
                print("  ", l)
