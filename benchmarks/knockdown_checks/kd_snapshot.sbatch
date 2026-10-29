#!/bin/bash
#SBATCH --partition=stud
#SBATCH --qos=stud
#SBATCH --job-name=kd_snapshot
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.out
# D-020 snapshot of EXP-20261016-LM-1's data root: the two Replogle 2022 screens, the registered CollecTRI
# human table and the human cisTarget v10 files, hard-linked (no copy) into one root that holds nothing else.
# Writes the snapshot manifest, checks it with --require-all, checks the registered references against the
# cumulative registry (experiments/datasets.tsv) and the screens against fetch_perturb.py's pinned Zenodo md5.
# usage, from the synced repo dir: sbatch benchmarks/knockdown_checks/kd_snapshot.sbatch
set -euo pipefail
ROOT=${KD_ROOT:-$HOME/grnbench_data_kd1}
KEEP=${KD_KEEP:-$HOME/lm_runs/kd_checks/snapshot}
PY=${PY:-$HOME/.conda/envs/grnbench/bin/python}
mkdir -p "$ROOT/perturb" "$ROOT/references" "$ROOT/cistarget" "$KEEP"
link() { [ -e "$2" ] || ln "$1" "$2"; }
link "$HOME/grnbench_data_perturb/ReplogleWeissman2022_rpe1.h5ad" "$ROOT/perturb/ReplogleWeissman2022_rpe1.h5ad"
link "$HOME/grnbench_data_perturb/ReplogleWeissman2022_K562_gwps.h5ad" "$ROOT/perturb/ReplogleWeissman2022_K562_gwps.h5ad"
link "$HOME/grnbench_data_fresh1/references/collectri_human.tsv" "$ROOT/references/collectri_human.tsv"
for f in hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
         hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
         motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl; do
  link "$HOME/grnbench_data_lm1/cistarget/$f" "$ROOT/cistarget/$f"
done
echo "host $(hostname) start $(date -Is) code $(cat GIT_COMMIT 2>/dev/null || echo unknown)"
find "$ROOT" -type f | sort
$PY benchmarks/data_manifest.py "$ROOT" > "$KEEP/EXP-20261016-LM-1.snapshot.tsv"
cat "$KEEP/EXP-20261016-LM-1.snapshot.tsv"
echo "== exact snapshot check"
$PY benchmarks/data_manifest.py "$ROOT" --check "$KEEP/EXP-20261016-LM-1.snapshot.tsv" --require-all
echo "== against the cumulative registry (the two screens are new rows, so UNLISTED is expected for them only)"
$PY benchmarks/data_manifest.py "$ROOT" --check experiments/datasets.tsv || true
echo "== screens against the pinned Zenodo md5 (fetch_perturb.py)"
( cd "$ROOT/perturb" && printf '%s  %s\n' cc7f1ec50aeb3a3e1b4a6cfa713d80fa ReplogleWeissman2022_rpe1.h5ad \
    13db594f8f1d2ccb88fec44a13e414dc ReplogleWeissman2022_K562_gwps.h5ad | md5sum -c - )
echo "end $(date -Is)"
