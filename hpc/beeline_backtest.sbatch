#!/bin/bash
#SBATCH --job-name=gb-beeline
#SBATCH --partition=stud
#SBATCH --qos=stud
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.out
#
# Run pinned BEELINE-derived selection and explicit-universe grnbench metrics.
#
#   bash hpc/fetch_beeline.sh                       # login node, once
#   sbatch --export=ALL,BEELINE=$HOME/beeline hpc/beeline_backtest.sbatch
#
# Variables:
#   BEELINE   corpus root created by fetch_beeline.sh   (default $HOME/beeline)
#   DATASETS  comma-separated subset                    (default: all found)
#   METHODS   comma-separated methods                   (default: the 3 built-ins)
#   NVAR      TFs + top-N variable genes                (default 500, BEELINE's setting)
#   ENV_NAME  conda env from hpc/setup_env.sh           (default grnbench)
set -uo pipefail
export PYTHONUTF8=1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

# sbatch runs a copy of this file from the spool dir, so BASH_SOURCE is only a fallback
# for a direct `bash hpc/...` run; under Slurm the checkout is where sbatch was called
REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
if [ ! -f "$REPO_DIR/pyproject.toml" ]; then
  echo "REPO_DIR=$REPO_DIR has no pyproject.toml: submit from the grnbench checkout (bash hpc/hpc.sh submit) or set REPO_DIR" >&2
  exit 1
fi
ENV_NAME="${ENV_NAME:-grnbench}"
CONDA_SH="${CONDA_SH:-/software/miniconda3/etc/profile.d/conda.sh}"
BEELINE="${BEELINE:-$HOME/beeline}"
NVAR="${NVAR:-500}"
METHODS="${METHODS:-spearman_signed,genie3,grnboost2}"

# shellcheck disable=SC1090
source "$CONDA_SH" && conda activate "$ENV_NAME" || {
  echo "cannot activate conda env '$ENV_NAME' - run: bash hpc/setup_env.sh" >&2
  exit 1
}
# Import grnbench from this checkout, not the conda env's editable install: the synced checkout is what GIT_COMMIT names.
export PYTHONPATH="$REPO_DIR/src${PYTHONPATH:+:$PYTHONPATH}"

NETWORKS="$(find "$BEELINE" -maxdepth 3 -type d -name "Networks" -printf '%h\n' 2>/dev/null | head -1)"
if [ -z "$NETWORKS" ]; then
  echo "no Networks/ directory under $BEELINE - run hpc/fetch_beeline.sh first" >&2
  exit 1
fi
DATA_ROOT="$(find "$BEELINE" -name ExpressionData.csv -printf '%h\n' 2>/dev/null | head -1)"
if [ -z "$DATA_ROOT" ]; then
  echo "no ExpressionData.csv under $BEELINE - run hpc/fetch_beeline.sh first" >&2
  exit 1
fi

OUT="${OUT:-$REPO_DIR/runs/beeline_${SLURM_JOB_ID:-manual}}"
echo "=== node ==="; hostname; nproc; free -g | head -2
echo "=== versions ==="; grnbench version; grnbench methods
echo "=== corpus ==="; echo "networks: $NETWORKS"; echo "beeline:  $BEELINE"

ARGS=(--data "$BEELINE" --networks "$NETWORKS" --out "$OUT"
      --methods "$METHODS" --n-var "$NVAR" --data-release zenodo:19009603:v4
      --threads "${SLURM_CPUS_PER_TASK:-4}")
[ -n "${DATASETS:-}" ] && ARGS+=(--datasets "$DATASETS")
[ -n "${PLUGIN:-}" ]   && ARGS+=(--plugin "$PLUGIN")

echo "=== beeline_backtest.py ${ARGS[*]} ==="
srun --unbuffered python "$REPO_DIR/benchmarks/beeline_backtest.py" "${ARGS[@]}"
RC=$?

echo "=== results ==="
[ -f "$OUT/beeline_backtest.csv" ] && cat "$OUT/beeline_backtest.csv"
echo "BEELINE_BACKTEST_DONE rc=$RC"
exit $RC
