#!/bin/bash
#SBATCH --job-name=gb-engine
#SBATCH --partition=stud
#SBATCH --qos=stud
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.out
#
# Mass test of the cleaning engine (benchmarks/engine_bench.py) as a SLURM array.
# Each array element runs CHUNK datasets, JOBS at a time on its cores.
#
#   python benchmarks/engine_bench.py plan --presets-reps 10 --random 300 --out $OUT/tasks.json
#   N=$(python -c "import json;print(len(json.load(open('$OUT/tasks.json'))))")
#   sbatch --array=0-$(( (N + CHUNK - 1) / CHUNK - 1 ))%20 \
#          --export=ALL,OUT=$OUT,CHUNK=8 hpc/engine_bench.sbatch
#   # real datasets (after hpc/fetch_engine_datasets.sh on the login node):
#   python benchmarks/real_bench.py plan --data-root ~/grnbench_data --out $OUT/real_tasks.json
#   sbatch --array=0-N%10 --export=ALL,OUT=$OUT,SCRIPT=benchmarks/real_bench.py,TASKS=$OUT/real_tasks.json,REFS=$HOME/grnbench_data/references,CHUNK=2 \
#          hpc/engine_bench.sbatch
#   # one job for a whole registered experiment (QOS stud: 1 running + 1 queued job), with models
#   # and the artifacts the figures need; GRNBENCH_CISTARGET_DIR for the *_motif models.
#   # --export splits on commas, so MODELS separates models with '+':
#   sbatch --time=24:00:00 --export=ALL,OUT=$OUT,SCRIPT=benchmarks/real_bench.py,TASKS=$OUT/real_tasks.json,REFS=$DATA/references,CHUNK=100000,MODELS=regdiffusion+regdiffusion_motif,KEEP_ARTIFACTS=1,GRNBENCH_CISTARGET_DIR=$DATA/cistarget \
#          hpc/engine_bench.sbatch
#   # when every element is done:
#   python benchmarks/engine_bench.py report --out $OUT
#
# Resumable: a dataset whose result file exists is skipped, so resubmitting the
# same array only runs what is missing.
set -uo pipefail
export PYTHONUTF8=1 PYTHONWARNINGS=ignore OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

REPO_DIR="${REPO_DIR:-$SLURM_SUBMIT_DIR}"
ENV_NAME="${ENV_NAME:-grnbench}"
CONDA_SH="${CONDA_SH:-/software/miniconda3/etc/profile.d/conda.sh}"
OUT="${OUT:?set OUT to the results directory holding tasks.json}"
CHUNK="${CHUNK:-8}"
JOBS="${JOBS:-${SLURM_CPUS_PER_TASK:-8}}"
WORK="${TMPDIR:-/tmp}/gb_engine_${SLURM_ARRAY_JOB_ID:-0}_${SLURM_ARRAY_TASK_ID:-0}"

# shellcheck disable=SC1090
source "$CONDA_SH" && conda activate "$ENV_NAME" || {
  echo "cannot activate conda env '$ENV_NAME' - run: bash hpc/setup_env.sh" >&2; exit 1; }
cd "$REPO_DIR" || exit 1
# Import grnbench from this checkout, whatever the conda env's editable install points at:
# CODE_COMMIT below names this checkout, so the code that runs must be this checkout's.
export PYTHONPATH="$REPO_DIR/src${PYTHONPATH:+:$PYTHONPATH}"
echo "node $(hostname)  cpus $JOBS  task ${SLURM_ARRAY_TASK_ID:-0}  chunk $CHUNK"
SCRIPT="${SCRIPT:-benchmarks/engine_bench.py}"       # or benchmarks/real_bench.py
TASKS="${TASKS:-$OUT/tasks.json}"
EXTRA=()
[ -n "${REFS:-}" ] && EXTRA+=(--refs "$REFS")         # real_bench: reference networks folder
[ -n "${MODELS:-}" ] && EXTRA+=(--models "${MODELS//+/,}")   # real_bench: 'a+b' (sbatch --export splits on ',')
[ "${KEEP_ARTIFACTS:-0}" = 1 ] && EXTRA+=(--keep-artifacts)   # real_bench: edges + processed data for viz
# The record cites the code the job ran: hpc.sh sync writes GIT_COMMIT next to the code.
mkdir -p "$OUT"
[ -f GIT_COMMIT ] && cp GIT_COMMIT "$OUT/CODE_COMMIT"
echo "code $(cat GIT_COMMIT 2>/dev/null || echo unknown)  python $(python --version 2>&1)  models ${MODELS:-default}"
python "$SCRIPT" run --tasks "$TASKS" --out "$OUT" "${EXTRA[@]}" \
       --index "${SLURM_ARRAY_TASK_ID:-0}" --chunk "$CHUNK" --jobs "$JOBS" --workroot "$WORK"
status=$?
rm -rf "$WORK"
exit $status
