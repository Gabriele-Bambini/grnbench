#!/bin/bash
# Copy the working tree to the cluster (no rsync needed; works from Git Bash on
# Windows through the Windows OpenSSH client and its ssh-agent).
#
#   bash hpc/sync_to_cluster.sh [HOST] [REMOTE_DIR]     # default: bocconi ~/grnbench
set -euo pipefail
HOST="${1:-bocconi}"
REMOTE="${2:-grnbench}"
SSH="${SSH:-ssh}"
[ -x /c/Windows/System32/OpenSSH/ssh.exe ] && SSH="${SSH_WIN:-/c/Windows/System32/OpenSSH/ssh.exe}"
cd "$(dirname "${BASH_SOURCE[0]}")/.."
tar czf - --exclude='./runs' --exclude='*/__pycache__' --exclude='*.pyc' --exclude='desktop.ini' \
    --exclude='./.projectos/checkpoints' --exclude='./dist' --exclude='./OUT' --exclude='./graphify-out' \
    --exclude='./.pytest_cache' . \
  | "$SSH" "$HOST" "mkdir -p ~/$REMOTE && tar xzf - -C ~/$REMOTE && echo synced to \$(hostname):~/$REMOTE"
