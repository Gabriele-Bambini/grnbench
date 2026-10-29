#!/bin/bash
#SBATCH --partition=stud
#SBATCH --qos=stud
#SBATCH --job-name=kd_norm
#SBATCH --output=%x_%j.out
# P9 of EXP-20261016-LM-1: truth with --norm cp10k and --norm core on one screen, same code and seeds.
# usage: sbatch [resources] kd_norm.sbatch <screen: rpe1|k562> <workers>
set -euo pipefail
SCREEN="$1"; WORKERS="$2"
case "$SCREEN" in
  rpe1) H5=$HOME/grnbench_data_perturb/ReplogleWeissman2022_rpe1.h5ad ;;
  k562) H5=$HOME/grnbench_data_perturb/ReplogleWeissman2022_K562_gwps.h5ad ;;
esac
TFS=$HOME/grnbench_data_fresh1/references/collectri_human.tsv
CODE=$HOME/lm_runs/kd_code
KEEP=$HOME/lm_runs/kd_checks/p9_$SCREEN
OUT=${TMPDIR:-/tmp}/kd_${SLURM_JOB_ID}
mkdir -p "$OUT" "$KEEP"
trap 'rm -rf "$OUT"' EXIT
cd "$CODE"
export PYTHONPATH="$CODE/src" OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=$HOME/.conda/envs/grnbench/bin/python
echo "host $(hostname) start $(date -Is) screen $SCREEN code $(cat CODE_COMMIT)"
$PY benchmarks/knockdown_bench.py subset --h5ad "$H5" --tfs "$TFS" --out "$OUT/subset.h5ad"
for NORM in cp10k core; do
  $PY benchmarks/knockdown_bench.py truth --h5ad "$OUT/subset.h5ad" --tfs "$TFS" --out "$OUT/truth_$NORM" --workers "$WORKERS" --norm "$NORM"
  $PY benchmarks/knockdown_bench.py units --truth "$OUT/truth_$NORM" > "$OUT/units_$NORM.json"
  mkdir -p "$KEEP/truth_$NORM"
  cp "$OUT/truth_$NORM"/units.tsv "$OUT/truth_$NORM"/responders.tsv "$OUT/truth_$NORM"/truth.json "$KEEP/truth_$NORM/"
  cp "$OUT/units_$NORM.json" "$KEEP/"
  echo "== $NORM"; cat "$OUT/units_$NORM.json"
done
echo "end $(date -Is)"
