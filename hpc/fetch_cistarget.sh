#!/bin/bash
# Download the cisTarget v10 files the motif models use (about 1.3 GB) into DEST.
# Login node: downloads only. These are the files hashed in experiments/datasets.tsv
# (cistarget/*); check them afterwards with
#   python benchmarks/data_manifest.py <DATA_ROOT> --check experiments/datasets.tsv
# and point the models at them with GRNBENCH_CISTARGET_DIR=DEST.
#
#   bash hpc/fetch_cistarget.sh ~/grnbench_data/cistarget
set -euo pipefail
DEST="${1:?usage: fetch_cistarget.sh DEST}"
BASE=https://resources.aertslab.org/cistarget
FILES=(
  databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
  databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
  databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
  databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
  motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
  motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl
)
mkdir -p "$DEST"
for f in "${FILES[@]}"; do
  out="$DEST/$(basename "$f")"
  if [ -s "$out" ]; then echo "have $(basename "$f") (the sha256 check decides if it is right)"; continue; fi
  echo "=== $(basename "$f")"
  curl -fL --retry 3 -o "$out.part" "$BASE/$f" && mv "$out.part" "$out"
done
echo "done: $(ls "$DEST" | wc -l) files in $DEST"
