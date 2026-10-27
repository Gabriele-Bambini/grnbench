#!/bin/bash
# One-time environment build for grnbench on the SLURM cluster.
#
#   bash hpc/setup_env.sh [ENV_NAME] [REPO_DIR]
#   WITH_TORCH=1 bash hpc/setup_env.sh      # also RegDiffusion + CPU torch (regdiffusion models)
#
# Creates a conda env (python 3.11) and installs the pinned stack from
# requirements/lock-py311.txt plus grnbench itself (editable). Safe to re-run.
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

echo "=== installing the pinned environment + grnbench from $REPO_DIR ==="
python -m pip install --upgrade pip
python -m pip install -r "$REPO_DIR/requirements/lock-py311.txt"
python -m pip install --no-deps -e "$REPO_DIR"

if [ "${WITH_TORCH:-0}" = 1 ]; then
  # RegDiffusion needs CPU torch, which is not in the lock: its pinned install line lives in the
  # lock header, so the cluster and the laptops install the same versions.
  echo "=== optional: RegDiffusion + CPU torch (install line from the lock header) ==="
  # tr: a lock synced from a Windows checkout has CRLF, and Linux bash keeps the '\r' (it would
  # end the extra index URL and pip would find nothing; Git Bash strips it, so Windows can't tell).
  torch_line="$(grep -m1 '^#[[:space:]]*python -m pip install regdiffusion' "$REPO_DIR/requirements/lock-py311.txt" \
                | tr -d '\r' | sed 's/^#[[:space:]]*//')"
  [ -n "$torch_line" ] || { echo "no RegDiffusion install line in the lock header" >&2; exit 1; }
  echo "$torch_line"
  eval "$torch_line"
  python -c "import torch, regdiffusion; print('torch', torch.__version__)"
fi

echo "=== what can actually run here ==="
grnbench version
grnbench methods
grnbench strategies

echo "=== optional: ARACNe3 ==="
if [ -n "${GRNBENCH_ARACNE3_BIN:-}" ] || command -v ARACNe3_app_release >/dev/null 2>&1 \
     || command -v ARACNe3_app >/dev/null 2>&1; then
  echo "ARACNe3 program set; 'grnbench methods' above says whether the aracne3 adapter can run"
else
  echo "ARACNe3 not built: the aracne3 adapter reports itself unavailable (by design)."
  echo "To enable it, build ARACNe3 (hpc/README.md) and set GRNBENCH_ARACNE3_BIN."
fi

echo "=== self-test ==="
cd "$REPO_DIR" && python -m pytest -q -p no:cacheprovider tests/test_engine.py tests/test_context.py tests/test_simulate.py tests/test_ui.py
echo SETUP_DONE
