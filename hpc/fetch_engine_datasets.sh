#!/bin/bash
# Fetch the real datasets for the engine mass test into an ordered folder tree.
# On a cluster, RUN ON THE LOGIN NODE (compute nodes usually have no network).
# Every download is listed here with its source; nothing else is fetched.
#
#   bash hpc/fetch_engine_datasets.sh [DEST] [CXG_MAX]     # default: $HOME/grnbench_data 30
#
#   DEST/
#     INDEX.md                  every dataset: folder, source, protocol, organism, cells, size
#     references/               CollecTRI + DoRothEA (ABC), human and mouse (OmniPath via decoupler)
#     beeline/                  BEELINE v1.1, Zenodo record 19009603 (~280 MB)
#     10x_pbmc3k/filtered|raw/  10x Genomics PBMC 3k, Cell Ranger 1.1.0 (~7 MB + ~30 MB)
#     scanpy_paul15/            Paul et al. 2015 mouse myeloid progenitors (scanpy.datasets.paul15)
#     cellxgene/<protocol>/<organism>/<tissue>__<id>/   up to CXG_MAX CELLxGENE Discover datasets,
#                               raw counts, <= 100 MB each, diverse protocols/species/tissues
set -euo pipefail
DEST="${1:-$HOME/grnbench_data}"
CXG_MAX="${2:-30}"
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="${PY:-python}"
export PYTHONPATH="$REPO_DIR/src${PYTHONPATH:+:$PYTHONPATH}"
mkdir -p "$DEST"/{references,beeline,10x_pbmc3k,scanpy_paul15,cellxgene}

echo "=== references (CollecTRI, DoRothEA) ==="
cd "$DEST/references"
"$PY" - <<'EOF'
from grnbench.tfs import decoupler_network
for org in ("human", "mouse"):
    for kind in ("collectri", "dorothea"):
        out = f"{kind}_{org}.tsv"
        try:
            net = decoupler_network(kind, organism=org, levels="ABC")
            net[["source", "target"]].drop_duplicates().to_csv(out, sep="\t", index=False)
            print(f"  {out}: {len(net)} edges")
        except Exception as exc:
            print(f"  FAILED {out}: {exc}")
EOF

echo "=== BEELINE ==="
bash "$REPO_DIR/hpc/fetch_beeline.sh" "$DEST/beeline"

echo "=== PBMC 3k (10x Genomics) ==="
cd "$DEST/10x_pbmc3k"
for kind in filtered raw; do
  f="pbmc3k_${kind}_gene_bc_matrices.tar.gz"
  [ -s "$f" ] || curl -fSL --retry 5 -o "$f" "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/$f"
  mkdir -p "$kind" && tar xzf "$f" -C "$kind"
done
cat > README.md <<'EOF'
PBMC 3k from a healthy donor (Homo sapiens), 10x Genomics Chromium (v1 chemistry), raw UMI counts.
filtered/ = Cell Ranger cell calls (~2,700 cells); raw/ = every barcode (unfiltered droplets).
Source: 10x Genomics datasets, pbmc3k, Cell Ranger 1.1.0.
EOF

echo "=== paul15 (scanpy) ==="
cd "$DEST/scanpy_paul15"
[ -s paul15.h5ad ] || "$PY" - <<'EOF'
import scanpy as sc
a = sc.datasets.paul15()
a.write_h5ad("paul15.h5ad")
open("README.md", "w").write("Mouse (Mus musculus) bone marrow myeloid progenitors, MARS-seq UMI counts, "
                             "differentiation time course (Paul et al. 2015, Cell). %d cells.\n" % a.n_obs)
print(a)
EOF

echo "=== CELLxGENE (<= ${CXG_MAX} datasets, <= 100 MB each) ==="
"$PY" "$REPO_DIR/benchmarks/fetch_cellxgene.py" --out "$DEST/cellxgene" --max "$CXG_MAX" --max-mb 100

echo "=== index ==="
"$PY" "$REPO_DIR/benchmarks/data_index.py" "$DEST"
echo "done: $DEST"
