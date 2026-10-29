#!/bin/bash
#SBATCH --partition=stud
#SBATCH --qos=stud
#SBATCH --job-name=kd_truth
#SBATCH --output=%x_%j.out
# usage: sbatch [resources] kd_truth.sbatch <screen: rpe1|k562> <workers>
# P3 (unit counts) and P7b (baseline on all genes) of EXP-20261016-LM-1, code 6c4e0ed.
set -euo pipefail
SCREEN="$1"; WORKERS="$2"
case "$SCREEN" in
  rpe1) H5=$HOME/grnbench_data_perturb/ReplogleWeissman2022_rpe1.h5ad ;;
  k562) H5=$HOME/grnbench_data_perturb/ReplogleWeissman2022_K562_gwps.h5ad ;;
esac
TFS=$HOME/grnbench_data_fresh1/references/collectri_human.tsv
CODE=$HOME/lm_runs/kd_code
OUT=$HOME/lm_runs/kd_checks/$SCREEN
mkdir -p "$OUT"
cd "$CODE"
export PYTHONPATH="$CODE/src" OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=$HOME/.conda/envs/grnbench/bin/python
echo "host $(hostname) start $(date -Is) screen $SCREEN workers $WORKERS code 6c4e0ed"
sha256sum "$TFS"
$PY benchmarks/knockdown_bench.py truth --h5ad "$H5" --tfs "$TFS" --out "$OUT/truth" --workers "$WORKERS"
$PY benchmarks/knockdown_bench.py units --truth "$OUT/truth" | tee "$OUT/units.json"
$PY benchmarks/knockdown_bench.py controls --h5ad "$H5" --out "$OUT/controls" --n 3000
$PY benchmarks/knockdown_bench.py baseline-check --h5ad "$OUT/controls/controls.h5ad" --workdir "$OUT/baseline" | tee "$OUT/baseline_check.json"
echo "end $(date -Is)"
