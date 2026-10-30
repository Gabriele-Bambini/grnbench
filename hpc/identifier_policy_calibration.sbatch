#!/bin/bash
#SBATCH --job-name=gb-id-policy
#SBATCH --partition=stud
#SBATCH --qos=stud
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out
#
# Production-ingest calibration for D-027.  This is deliberately one sequential
# compute job: every approved file is loaded once, no work runs on the login node.
set -euo pipefail
export PYTHONUTF8=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

REPO_DIR="${REPO_DIR:-$SLURM_SUBMIT_DIR}"
ENV_NAME="${ENV_NAME:-grnbench}"
CONDA_SH="${CONDA_SH:-/software/miniconda3/etc/profile.d/conda.sh}"
EXP_ID="${EXP_ID:?set EXP_ID to the registered experiment id}"
DATA_ROOT="${DATA_ROOT:-$HOME/grnbench_data}"
OUT="${OUT:-$HOME/gb_runs/$EXP_ID}"
SNAPSHOT="${SNAPSHOT:-experiments/$EXP_ID.snapshot.tsv}"
SCRATCH_BASE="${SLURM_TMPDIR:-/tmp}"
INPUT_ROOT="$SCRATCH_BASE/grnbench-${EXP_ID}-${SLURM_JOB_ID:-manual}-inputs"

# shellcheck disable=SC1090
source "$CONDA_SH" && conda activate "$ENV_NAME" || {
  echo "cannot activate conda env '$ENV_NAME' - run: bash hpc/setup_env.sh" >&2
  exit 1
}
cd "$REPO_DIR"
export PYTHONPATH="$REPO_DIR/src${PYTHONPATH:+:$PYTHONPATH}"
mkdir -p "$OUT"
[ -f GIT_COMMIT ] && cp GIT_COMMIT "$OUT/CODE_COMMIT"
[ -f "$SNAPSHOT" ] || { echo "missing registered snapshot: $SNAPSHOT" >&2; exit 1; }
[ ! -e "$INPUT_ROOT" ] || { echo "refusing existing staging root" >&2; exit 1; }
mkdir -p "$INPUT_ROOT"

# The shared cache is intentionally broader and dynamic.  Stage only the
# preregistered files as symlinks, so --require-all verifies an exact root
# without copying expression matrices or treating unrelated cache files as data.
while IFS=$'\t' read -r rel _rest; do
  [ "$rel" = "path" ] && continue
  case "$rel" in
    ""|/*|../*|*/../*|*/..) echo "unsafe snapshot path: $rel" >&2; exit 1 ;;
  esac
  src="$DATA_ROOT/$rel"
  dst="$INPUT_ROOT/$rel"
  [ -f "$src" ] || { echo "missing snapshot file: $rel" >&2; exit 1; }
  mkdir -p "$(dirname "$dst")"
  ln -s "$src" "$dst"
done < "$SNAPSHOT"
cp "$SNAPSHOT" "$OUT/input-snapshot.tsv"

echo "code $(cat GIT_COMMIT 2>/dev/null || echo unknown)"
python --version
python benchmarks/data_manifest.py "$INPUT_ROOT" --check "$SNAPSHOT" --require-all \
  | tee "$OUT/data-manifest-check.txt"
python benchmarks/identifier_policy_calibration.py "$INPUT_ROOT" \
  --tsv "$OUT/identifier-policy.tsv" | tee "$OUT/summary.txt"
