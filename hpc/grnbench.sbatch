#!/bin/bash
#SBATCH --job-name=grnbench
#SBATCH --partition=stud
#SBATCH --qos=stud
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.out
#
# Full grnbench loop on a compute node. CPU only — arboreto/GENIE3/GRNBoost2 are
# dask-parallel over cores, nothing here wants a GPU.
#
#   sbatch --export=ALL,DATA=/path/d.h5ad,TFS=/path/tfs.txt,GOLD=dorotheaABC \
#          hpc/grnbench.sbatch
#
# Variables (all overridable via --export):
#   DATA     input dataset                              (required)
#   OUT      run directory                              (default: $REPO/runs/$SLURM_JOB_ID)
#   TFS      TF symbol list                             (optional but recommended)
#   GOLD     gold pairs file or dorothea<levels> spec    (optional; omit to skip eval)
#   METHODS  comma-separated names or 'available'       (default: available)
#   PLUGINS  SEMICOLON-separated --plugin arguments     (optional; ';' not ':',
#            because plugin references themselves contain a colon)
#   STABILITY   1 (default) or 0 to skip subsample stability
#   STABILITY_N number of subsamples                    (default: 5)
#   NHVG     HVG count                                  (default: 3000)
#   GROUPKEY obs column for pseudobulk/recommendations  (optional)
#   ENV_NAME conda env built by hpc/setup_env.sh        (default: grnbench)
set -uo pipefail
export PYTHONUTF8=1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

REPO_DIR="${REPO_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
ENV_NAME="${ENV_NAME:-grnbench}"
CONDA_SH="${CONDA_SH:-/software/miniconda3/etc/profile.d/conda.sh}"
# shellcheck disable=SC1090
source "$CONDA_SH" && conda activate "$ENV_NAME" || {
  echo "cannot activate conda env '$ENV_NAME' — run: bash hpc/setup_env.sh" >&2
  exit 1
}

: "${DATA:?set DATA=/path/to/dataset.h5ad}"
OUT="${OUT:-$REPO_DIR/runs/${SLURM_JOB_ID:-manual}}"
METHODS="${METHODS:-available}"
NHVG="${NHVG:-3000}"

echo "=== node ==="
hostname; nproc; free -g | head -2
echo "=== versions ==="
grnbench version
grnbench methods

ARGS=(run --data "$DATA" --out "$OUT" --methods "$METHODS"
      --n-hvg "$NHVG" --threads "${SLURM_CPUS_PER_TASK:-4}")
# Stability re-runs every method on N cell subsamples: worth it, but it is an
# N-fold cost multiplier, so it must be switchable on a machine where runs are
# measured in hours. STABILITY=0 turns it off.
if [ "${STABILITY:-1}" != "0" ]; then
  ARGS+=(--stability --stability-n "${STABILITY_N:-5}")
fi
[ -n "${TFS:-}" ]      && ARGS+=(--tf-list "$TFS")
[ -n "${GOLD:-}" ]     && ARGS+=(--gold "$GOLD")
[ -n "${GROUPKEY:-}" ] && ARGS+=(--group-key "$GROUPKEY")
if [ -n "${PLUGINS:-}" ]; then
  # Separator is ';' and NOT ':' - every plugin spelling grnbench accepts
  # contains a colon ("pkg.mod:Class", "path/to/model.py:Class"), so splitting
  # on ':' would tear each reference in half.
  IFS=';' read -ra _plugs <<< "$PLUGINS"
  for p in "${_plugs[@]}"; do [ -n "$p" ] && ARGS+=(--plugin "$p"); done
fi

echo "=== grnbench ${ARGS[*]} ==="
srun --unbuffered grnbench "${ARGS[@]}"
RC=$?

echo "=== artifacts ==="
ls -la "$OUT" "$OUT"/infer 2>/dev/null
[ -f "$OUT/eval/eval_report.md" ] && cat "$OUT/eval/eval_report.md"
echo "HTML report: $OUT/eval/eval_report.html"
echo "GRNBENCH_DONE rc=$RC"
exit $RC
