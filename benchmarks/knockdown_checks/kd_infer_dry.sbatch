#!/bin/bash
#SBATCH --partition=stud
#SBATCH --qos=stud
#SBATCH --job-name=kd_infer
#SBATCH --output=%x_%j.out
# Dry run of the inference step of EXP-20261016-LM-1 on one screen's control cells: builds every arm,
# scores nothing (no responder is read). Checks resources, the all-genes pipeline and real fit times.
# usage: sbatch [resources] kd_infer_dry.sbatch <screen: rpe1|k562>
set -euo pipefail
SCREEN="$1"
case "$SCREEN" in
  rpe1) H5=$HOME/grnbench_data_perturb/ReplogleWeissman2022_rpe1.h5ad ;;
  k562) H5=$HOME/grnbench_data_perturb/ReplogleWeissman2022_K562_gwps.h5ad ;;
esac
TFS=$HOME/grnbench_data_fresh1/references/collectri_human.tsv
export GRNBENCH_CISTARGET_DIR=$HOME/grnbench_data_lm1/cistarget
CODE=$HOME/lm_runs/kd_code
KEEP=$HOME/lm_runs/kd_checks/infer_dry_$SCREEN
OUT=${TMPDIR:-/tmp}/kd_${SLURM_JOB_ID}
mkdir -p "$OUT" "$KEEP"
trap 'rm -rf "$OUT"' EXIT
cd "$CODE"
export PYTHONPATH="$CODE/src" OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
PY=$HOME/.conda/envs/grnbench/bin/python
echo "host $(hostname) start $(date -Is) screen $SCREEN code $(cat CODE_COMMIT)"
$PY benchmarks/knockdown_bench.py subset --h5ad "$H5" --tfs "$TFS" --out "$OUT/subset.h5ad"
$PY benchmarks/knockdown_bench.py controls --h5ad "$OUT/subset.h5ad" --out "$OUT/controls" --n 10000000
$PY benchmarks/knockdown_score.py infer --experiment EXP-20261016-LM-1 --controls "$OUT/controls/controls.h5ad" --tfs "$TFS" --out "$OUT/infer" --device cpu
cp "$OUT/infer/infer.json" "$OUT/infer/regulators.tsv" "$KEEP/"
ls -la "$OUT/infer/edges" | tee "$KEEP/edges_listing.txt"
cat "$KEEP/infer.json"
echo "end $(date -Is)"
