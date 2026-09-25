#!/bin/bash
# One-time environment build for grnbench on the SLURM cluster.
#
#   bash hpc/setup_env.sh [ENV_NAME] [REPO_DIR]
#
# Creates a conda env with the CPU-only stack grnbench needs and installs the
# package with the `all` extra (arboreto/GENIE3/GRNBoost2, decoupler/DoRothEA,
# pyarrow for parquet plugin exchange, matplotlib). Safe to re-run.
set -euo pipefail

ENV_NAME="${1:-grnbench}"
REPO_DIR="${2:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CONDA_SH="${CONDA_SH:-/software/miniconda3/etc/profile.d/conda.sh}"

if [ ! -f "$CONDA_SH" ]; then
  echo "conda profile not found at $CONDA_SH — set CONDA_SH=/path/to/profile.d/conda.sh" >&2
  exit 1
fi
# shellcheck disable=SC1090
source "$CONDA_SH"

if ! conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "=== creating conda env '$ENV_NAME' (python 3.11) ==="
  conda create -y -n "$ENV_NAME" python=3.11 || exit 1
fi
conda activate "$ENV_NAME" || exit 1

echo "=== installing grnbench[all] from $REPO_DIR ==="
python -m pip install --upgrade pip
python -m pip install -e "$REPO_DIR[all,dev]"

echo "=== what can actually run here ==="
grnbench version
grnbench methods
grnbench strategies

echo "=== optional: ARACNe3 ==="
if command -v java >/dev/null 2>&1; then
  echo "java: $(java -version 2>&1 | head -1)"
  echo "set ARACNE3_JAR=/path/to/aracne3.jar to enable the aracne3 adapter"
else
  echo "java not on PATH: the aracne3 adapter will report itself unavailable (by design)"
fi

echo "=== self-test ==="
python -m pytest "$REPO_DIR/tests" -q
echo SETUP_DONE
