#!/bin/bash
#SBATCH --job-name=gb-pyscenic
#SBATCH --partition=stud
#SBATCH --qos=stud
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out
set -euo pipefail

export PYTHONUTF8=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
REPO_DIR="${REPO_DIR:-$SLURM_SUBMIT_DIR}"
ENV_NAME="${ENV_NAME:-grnbench}"
PYSCENIC_ENV="${PYSCENIC_ENV:-grnbench-pyscenic-0.12.1}"
CONDA_SH="${CONDA_SH:-/software/miniconda3/etc/profile.d/conda.sh}"
SOURCE_RUN="${SOURCE_RUN:?set SOURCE_RUN to the completed two-dataset real_bench directory}"
OUT="${OUT:?set OUT to the EXP-20261018-LD-1 comparison directory}"
CISTARGET_DIR="${CISTARGET_DIR:?set CISTARGET_DIR to the registered cisTarget snapshot}"
UNITS="${UNITS:-$REPO_DIR/experiments/EXP-20261018-LD-1.units.tsv}"
JOBS="${SLURM_CPUS_PER_TASK:-8}"

# shellcheck disable=SC1090
source "$CONDA_SH"
PYSCENIC_PREFIX="$(conda env list --json | python -c 'import json,os,sys; n=sys.argv[1]; p=[x for x in json.load(sys.stdin)["envs"] if os.path.basename(x)==n]; print(p[0] if len(p)==1 else "")' "$PYSCENIC_ENV")"
[ -x "$PYSCENIC_PREFIX/bin/pyscenic" ] || {
  echo "isolated pySCENIC env missing; run hpc/setup_pyscenic_env.sh" >&2; exit 1; }
conda activate "$ENV_NAME"
cd "$REPO_DIR"
export PYTHONPATH="$REPO_DIR/src${PYTHONPATH:+:$PYTHONPATH}"
mkdir -p "$OUT" logs
[ -f GIT_COMMIT ] && cp GIT_COMMIT "$OUT/CODE_COMMIT"
"$PYSCENIC_PREFIX/bin/python" -m pip freeze --all | sort > "$OUT/pyscenic-0.12.1.freeze.txt"

python benchmarks/pyscenic_compare.py prepare \
  --source-run "$SOURCE_RUN" --units "$UNITS" --cistarget-dir "$CISTARGET_DIR" --out "$OUT"
python benchmarks/pyscenic_compare.py run \
  --prepared "$OUT" --cistarget-dir "$CISTARGET_DIR" \
  --pyscenic "$PYSCENIC_PREFIX/bin/pyscenic" --workers "$JOBS"
python benchmarks/pyscenic_compare.py report \
  --prepared "$OUT" --units "$UNITS" --out "$OUT/portable"

echo "PYSCENIC_COMPARISON_DONE"
