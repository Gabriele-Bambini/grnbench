#!/bin/bash
#SBATCH --partition=stud
#SBATCH --qos=stud
#SBATCH --job-name=kd_p6
#SBATCH --output=%x_%j.out
# usage: sbatch [resources] kd_p6.sbatch <python> <args...>
set -euo pipefail
cd ~/lm_runs/kd_checks
echo "host $(hostname) start $(date -Is) args $*"
nvidia-smi -L 2>/dev/null || true
PY="$1"; shift
export PYTHONPATH=$HOME/lm_runs/kd_checks/src; case "$PY" in *venv-tsgnn*) export PYTHONPATH="$HOME/lm_runs/kd_checks/pylib:$PYTHONPATH";; esac
OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}" "$PY" kd_p6_timing.py "$@"
echo "end $(date -Is)"
