#!/bin/bash
# Build the isolated reference environment used only by EXP-20261018-LD-1.
# Safe to re-run.  Set RECREATE=1 to remove and rebuild this one named env.
set -euo pipefail

ENV_NAME="${1:-grnbench-pyscenic-0.12.1}"
REPO_DIR="${2:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CONDA_SH="${CONDA_SH:-/software/miniconda3/etc/profile.d/conda.sh}"

[ -f "$CONDA_SH" ] || { echo "conda profile not found; set CONDA_SH" >&2; exit 1; }
# shellcheck disable=SC1090
source "$CONDA_SH"

if [ "${RECREATE:-0}" = 1 ] && conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  conda env remove -y -n "$ENV_NAME"
fi
if ! conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  conda create -y -c conda-forge -n "$ENV_NAME" python=3.10 pip=23.3.2
fi
conda activate "$ENV_NAME"
python -m pip install -r "$REPO_DIR/requirements/pyscenic-0.12.1-py310.txt"

python - <<'PY'
import numpy, pandas, pyscenic
print("pyscenic", pyscenic.__version__)
print("numpy", numpy.__version__)
print("pandas", pandas.__version__)
PY
pyscenic --help >/dev/null
python -m pip freeze --all | sort
echo "PYSCENIC_SETUP_DONE env=$ENV_NAME"
