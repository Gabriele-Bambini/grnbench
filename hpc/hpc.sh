#!/bin/bash
# Team helper for the Bocconi HPC cluster. Nothing user-specific lives here:
# your host alias, user and remote dir come from .env.hpc (untracked, copy
# .env.hpc.example) and your own ~/.ssh/config. See docs/HPC.md.
#
#   bash hpc/hpc.sh check                      # reach the cluster without prompting?
#   bash hpc/hpc.sh sync                       # copy the working tree to ~/$GRNBENCH_HPC_DIR
#   bash hpc/hpc.sh run '<cmd>'                # run <cmd> in the remote repo dir (login node: light work only)
#   bash hpc/hpc.sh submit <sbatch args...>    # sbatch --parsable in the remote repo dir, prints the job id
#   bash hpc/hpc.sh status [JOBID]             # your queue; with JOBID also sacct (the source of truth)
#   bash hpc/hpc.sh log JOBID [LINES]          # tail the job's log files
#   bash hpc/hpc.sh cancel JOBID               # scancel one of your jobs
#   bash hpc/hpc.sh fetch REMOTE_PATH [LOCAL_DIR]   # copy a results folder back (tar over ssh)
#
# ssh runs with BatchMode=yes: it never asks for a passphrase. If your key is
# not loaded in ssh-agent the command fails at once instead of hanging.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [ -f "$ROOT/.env.hpc" ]; then set -a; . "$ROOT/.env.hpc"; set +a; fi
HOST="${GRNBENCH_HPC_HOST:?set GRNBENCH_HPC_HOST in .env.hpc (copy .env.hpc.example, see docs/HPC.md)}"
[ -n "${GRNBENCH_HPC_USER:-}" ] && HOST="$GRNBENCH_HPC_USER@$HOST"
RDIR="${GRNBENCH_HPC_DIR:-grnbench}"

SSH="${SSH:-ssh}"
# Windows: Git Bash's own ssh cannot reach the Windows ssh-agent service, the
# Windows OpenSSH client can. SSH_AUTH_SOCK would point it at the wrong agent.
if [ -x /c/Windows/System32/OpenSSH/ssh.exe ] && [ "$SSH" = ssh ]; then
  SSH=/c/Windows/System32/OpenSSH/ssh.exe; unset SSH_AUTH_SOCK
fi
SSH_OPTS=(-o BatchMode=yes -o ConnectTimeout=25 -o ServerAliveInterval=30)

r() {
  local err rc=0
  err="$(mktemp)"
  "$SSH" "${SSH_OPTS[@]}" "$HOST" "$@" 2>"$err" || rc=$?
  if [ $rc -ne 0 ]; then
    cat "$err" >&2
    if grep -q "Not allowed at this time" "$err"; then
      echo "LOCKED OUT after failed logins. Stop. Wait 10 minutes, then fix the key before retrying." >&2
    elif grep -qi "permission denied" "$err"; then
      echo "Auth failed: load your key yourself (ssh-add ~/.ssh/<key>) and check User in ~/.ssh/config. Do not retry in a loop." >&2
    fi
  fi
  rm -f "$err"
  return $rc
}

cmd="${1:-}"; shift || true
case "$cmd" in
  check)
    r "echo connected: \$(whoami)@\$(hostname); sinfo --version; ls -d ~/$RDIR 2>/dev/null || echo 'no ~/$RDIR yet: run sync'" ;;
  sync)
    cd "$ROOT"
    # Record which code the cluster runs (experiment records need the commit).
    rev="$(git rev-parse HEAD 2>/dev/null || echo unknown)"
    [ -n "$(git status --porcelain 2>/dev/null)" ] && rev="$rev-dirty"
    tar czf - --exclude='./.git' --exclude='./.venv' --exclude='./runs' --exclude='./logs' \
        --exclude='*/__pycache__' --exclude='*.pyc' --exclude='desktop.ini' --exclude='./.env.hpc' \
        --exclude='./.projectos/checkpoints' --exclude='./.projectos/.runtime' --exclude='./dist' \
        --exclude='./OUT' --exclude='./graphify-out' --exclude='./.pytest_cache' --exclude='./data' . \
      | r "mkdir -p ~/$RDIR && tar xzf - -C ~/$RDIR && echo '$rev' > ~/$RDIR/GIT_COMMIT && echo synced $rev to \$(hostname):~/$RDIR" ;;
  run)
    r "cd ~/$RDIR && $*" ;;
  submit)
    r "cd ~/$RDIR && mkdir -p logs && sbatch --parsable $*" ;;
  status)
    r "squeue --me -o '%.12i %.20j %.9T %.10M %.10l %.5D %R'${1:+; sacct -j $1 -X -o JobID%20,JobName%20,State,Elapsed,ExitCode}" ;;
  log)
    id="${1:?log JOBID}"; n="${2:-60}"
    r "cd ~/$RDIR && for f in \$(ls logs/*_${id}*.out *_${id}.out 2>/dev/null); do echo \"== \$f\"; tail -n $n \"\$f\"; done" ;;
  cancel)
    r "scancel ${1:?cancel JOBID}" ;;
  fetch)
    src="${1:?fetch REMOTE_PATH [LOCAL_DIR]}"; dest="${2:-.}"
    mkdir -p "$dest"
    r "tar czf - -C $(dirname "$src") $(basename "$src")" | tar xzf - -C "$dest"
    echo "fetched $src -> $dest/$(basename "$src")" ;;
  *)
    sed -n '2,15p' "${BASH_SOURCE[0]}"; exit 2 ;;
esac
