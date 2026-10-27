#!/bin/bash
#SBATCH --job-name=gb-aracne3-smoke
#SBATCH --partition=stud
#SBATCH --qos=stud
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out
#
# ARACNe3 smoke run (EXP-20261018-GB-1): benchmarks/aracne3_smoke.py calls the aracne3 adapter
# against the program that hpc/aracne3_build.sbatch built. 4 CPUs: two of the runs ask ARACNe3
# for --threads 4 (the script asks with aracne3_threads: the adapter does not pass on the engine's
# generic threads).
#
#   sbatch --export=ALL,OUT=$HOME/gb_runs/EXP-20261018-GB-1 hpc/aracne3_smoke.sbatch
#
# GRNBENCH_ARACNE3_BIN is set here, in the job, and nowhere else: aracne3 is a default method, so
# wherever the variable is set, --methods available runs it. OMP_NUM_THREADS is not set: ARACNe3
# takes its thread count from --threads (num_threads clauses in its source).
set -uo pipefail
export PYTHONUTF8=1 PYTHONWARNINGS=ignore

REPO_DIR="${REPO_DIR:-$SLURM_SUBMIT_DIR}"
ENV_NAME="${ENV_NAME:-grnbench}"
CONDA_SH="${CONDA_SH:-/software/miniconda3/etc/profile.d/conda.sh}"
OUT="${OUT:?set OUT to the run folder}"
DATA="${DATA:-$HOME/grnbench_data}"
ARACNE3_DIR="${ARACNE3_DIR:-$HOME/tools/ARACNe3}"

# shellcheck disable=SC1090
source "$CONDA_SH" && conda activate "$ENV_NAME" || {
  echo "cannot activate conda env '$ENV_NAME' - run: bash hpc/setup_env.sh" >&2; exit 1; }
cd "$REPO_DIR" || exit 1
# Import grnbench from this checkout, whatever the conda env's editable install points at.
export PYTHONPATH="$REPO_DIR/src${PYTHONPATH:+:$PYTHONPATH}"
export GRNBENCH_ARACNE3_BIN="$ARACNE3_DIR/build/src/app/ARACNe3_app_release"
# The adapter's temporary folders (ARACNe3's working folders) on the node, not under home.
export TMPDIR="/tmp/gb_aracne3_${SLURM_JOB_ID:-0}"
mkdir -p "$TMPDIR" "$OUT"
[ -f GIT_COMMIT ] && cp GIT_COMMIT "$OUT/CODE_COMMIT"
echo "code $(cat GIT_COMMIT 2>/dev/null || echo unknown)  python $(python --version 2>&1)  cpus ${SLURM_CPUS_PER_TASK:-?}"
python benchmarks/aracne3_smoke.py --data-root "$DATA" --aracne3-src "$ARACNE3_DIR" --out "$OUT"
status=$?
rm -rf "$TMPDIR"
echo "exit status $status"
exit $status
