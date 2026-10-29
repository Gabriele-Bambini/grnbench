#!/bin/bash
#SBATCH --partition=stud
#SBATCH --qos=stud
#SBATCH --job-name=kd_tsub
#SBATCH --output=%x_%j.out
# usage: sbatch [resources] kd_truth_sub.sbatch <screen: rpe1|k562> <workers>
# P3 and P7b of EXP-20261016-LM-1 through a one-pass subset copy (code fff93ec).
# For rpe1 the truth tables are compared byte for byte with the direct run of job 703375.
set -euo pipefail
SCREEN="$1"; WORKERS="$2"
case "$SCREEN" in
  rpe1) H5=$HOME/grnbench_data_perturb/ReplogleWeissman2022_rpe1.h5ad ;;
  k562) H5=$HOME/grnbench_data_perturb/ReplogleWeissman2022_K562_gwps.h5ad ;;
esac
TFS=$HOME/grnbench_data_fresh1/references/collectri_human.tsv
CODE=$HOME/lm_runs/kd_code
KEEP=$HOME/lm_runs/kd_checks/${SCREEN}_sub
# home quota is nearly full: work on the node-local disk, copy only the small tables back
OUT=${TMPDIR:-/tmp}/kd_${SLURM_JOB_ID}
mkdir -p "$OUT" "$KEEP"
df -h "$OUT" | tail -1
trap 'rm -rf "$OUT"' EXIT
cd "$CODE"
export PYTHONPATH="$CODE/src" OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=$HOME/.conda/envs/grnbench/bin/python
echo "host $(hostname) start $(date -Is) screen $SCREEN workers $WORKERS code fff93ec"
$PY benchmarks/knockdown_bench.py subset --h5ad "$H5" --tfs "$TFS" --out "$OUT/subset.h5ad"
echo "subset done $(date -Is)"
$PY benchmarks/knockdown_bench.py truth --h5ad "$OUT/subset.h5ad" --tfs "$TFS" --out "$OUT/truth" --workers "$WORKERS"
$PY benchmarks/knockdown_bench.py units --truth "$OUT/truth" | tee "$OUT/units.json"
if [ "$SCREEN" = rpe1 ]; then
  for f in units.tsv responders.tsv; do
    cmp "$OUT/truth/$f" "$HOME/lm_runs/kd_checks/rpe1/truth/$f" && echo "IDENTICAL $f" || echo "DIFFERENT $f"
  done
else
  $PY benchmarks/knockdown_bench.py controls --h5ad "$OUT/subset.h5ad" --out "$OUT/controls" --n 3000
  $PY benchmarks/knockdown_bench.py baseline-check --h5ad "$OUT/controls/controls.h5ad" --workdir "$OUT/baseline" | tee "$OUT/baseline_check.json"
fi
mkdir -p "$KEEP/truth"
cp "$OUT/truth/units.tsv" "$OUT/truth/responders.tsv" "$OUT/truth/truth.json" "$KEEP/truth/"
cp "$OUT/units.json" "$OUT/subset.h5ad.subset.json" "$KEEP/" 2>/dev/null || cp "$OUT/units.json" "$KEEP/"
[ -f "$OUT/baseline_check.json" ] && cp "$OUT/baseline_check.json" "$KEEP/"
ls -la "$OUT" "$KEEP/truth"
echo "end $(date -Is)"
