"""P9 metrics for EXP-20261016-LM-1: m1 usable units, m2 share of decreased responder rows,
m3 genes responding in at least half of the usable TFs. usage: python p9_metrics.py TRUTH_DIR"""
import csv
import json
import sys
from collections import Counter

d = sys.argv[1]
units = list(csv.DictReader(open(f"{d}/units.tsv"), delimiter="\t"))
usable = {u["tf"] for u in units if u["usable"] == "True"}
rows = [r for r in csv.DictReader(open(f"{d}/responders.tsv"), delimiter="\t") if r["tf"] in usable]
neg = sum(float(r["log2fc"]) < 0 for r in rows)
per_gene = Counter(r["gene"] for r in rows)
half = len(usable) / 2
print(json.dumps({
    "m1_usable": len(usable),
    "responder_rows": len(rows),
    "m2_decreased_share": round(neg / len(rows), 4) if rows else None,
    "m3_genes_in_half_of_tfs": sum(c >= half for c in per_gene.values()),
}))
