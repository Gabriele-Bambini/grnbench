#!/bin/bash
#SBATCH --partition=stud
#SBATCH --qos=stud
#SBATCH --job-name=kd_run
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
# EXP-20261016-LM-1, the registered run (record section 2), both screens in one job, in the record's order:
# 1. D-020 check of the data root against the committed snapshot (any mismatch stops the job);
# 2. truth for RPE1 and K562: responder tables written and hashed (truth.json) before any network is scored;
# 3. per screen: inference on the inference-half control cells only, then scoring against the hashed truth;
# 4. the cross-screen verdicts (design RPE1, replication K562).
# Resumable: a finished step leaves its folder under $RUN and is skipped by a resubmission of the same commit
# (record exclusion (c)); a step writes to the node-local disk and moves its folder into $RUN only when done.
# $RUN/run.json (commit, snapshot hash, models, permutations, normalisation) is written before the first step;
# a resubmission that differs in any of them, or finds steps that run.json does not vouch for, stops before work.
# An inference with a failed model exits 3 and publishes nothing; the failures per arm are kept in
# $RUN/<screen>/infer_attempts.json, and only a second failure of the same arm on the same commit is published.
# An arm is written there as `started` before it runs: a job killed by Slurm inside it leaves it `started`, and the
# next run of the same commit counts that as one failed attempt (a timeout is a failure like a caught error).
# usage, from the repo dir synced by `bash hpc/hpc.sh sync` (clean tree, so GIT_COMMIT holds the commit):
#   bash hpc/hpc.sh submit hpc/knockdown_run.sbatch
# KD_ROOT, KD_RUN, KD_SNAPSHOT, KD_MODELS, KD_PERMS and PY exist for the synthetic smoke test only; the
# registered run uses the defaults below, and the log prints every value.
set -euo pipefail
ROOT=${KD_ROOT:-$HOME/grnbench_data_kd1}
RUN=${KD_RUN:-$HOME/lm_runs/EXP-20261016-LM-1}
SNAP=${KD_SNAPSHOT:-experiments/EXP-20261016-LM-1.snapshot.tsv}
PY=${PY:-$HOME/.conda/envs/grnbench/bin/python}
PERMS=${KD_PERMS:-10000}
MODELS=()
[ -n "${KD_MODELS:-}" ] && MODELS=(--models "$KD_MODELS")
CPUS=${SLURM_CPUS_PER_TASK:-8}
NORM=core
TFS=$ROOT/references/collectri_human.tsv
export GRNBENCH_CISTARGET_DIR=$ROOT/cistarget PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
TMP=${TMPDIR:-/tmp}/kd_run_${SLURM_JOB_ID:-$$}
mkdir -p "$RUN" "$TMP"
trap 'rm -rf "$TMP"' EXIT
h5ad() { case "$1" in rpe1) echo "$ROOT/perturb/ReplogleWeissman2022_rpe1.h5ad" ;;
                      k562) echo "$ROOT/perturb/ReplogleWeissman2022_K562_gwps.h5ad" ;; esac; }
other() { if [ "$1" = rpe1 ]; then echo k562; else echo rpe1; fi; }
log() { echo "$(date -Is) $*"; }
subset() { [ -f "$TMP/$1.h5ad" ] || $PY benchmarks/knockdown_bench.py subset --h5ad "$(h5ad "$1")" --tfs "$TFS" --out "$TMP/$1.h5ad"; }
publish() { rm -rf "$2.partial"; cp -r "$1" "$2.partial"; mv "$2.partial" "$2"; }   # a folder in $RUN is complete

identity() { $PY benchmarks/knockdown_score.py run-identity --run "$RUN" --snapshot "$SNAP" "${MODELS[@]}" \
               --perms "$PERMS" --norm "$NORM"; }

log "host $(hostname) code $(cat GIT_COMMIT 2>/dev/null || echo unknown) root $ROOT run $RUN snapshot $SNAP" \
    "models ${KD_MODELS:-record default} perms $PERMS cpus $CPUS"
identity
$PY benchmarks/data_manifest.py "$ROOT" --check "$SNAP" --require-all

for S in rpe1 k562; do
  if [ ! -d "$RUN/$S/truth" ]; then
    log "truth $S"
    subset "$S"
    OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 $PY benchmarks/knockdown_bench.py truth \
      --h5ad "$TMP/$S.h5ad" --tfs "$TFS" --norm "$NORM" --workers "$CPUS" --out "$TMP/truth_$S"
    $PY benchmarks/knockdown_bench.py total-shift --h5ad "$TMP/$S.h5ad" --truth "$TMP/truth_$S" \
      --out "$TMP/truth_$S/total_shift.tsv"
    mkdir -p "$RUN/$S"
    publish "$TMP/truth_$S" "$RUN/$S/truth"
  fi
  $PY benchmarks/knockdown_bench.py units --truth "$RUN/$S/truth"
done

for S in rpe1 k562; do
  if [ ! -d "$RUN/$S/infer" ]; then
    log "inference $S (control cells only)"
    subset "$S"
    $PY benchmarks/knockdown_bench.py controls --h5ad "$TMP/$S.h5ad" --out "$TMP/controls_$S" --norm "$NORM"
    OMP_NUM_THREADS=$CPUS $PY benchmarks/knockdown_score.py infer --controls "$TMP/controls_$S/controls.h5ad" \
      --experiment EXP-20261016-LM-1 --tfs "$TFS" "${MODELS[@]}" --device cpu --out "$TMP/infer_$S" \
      --attempts "$RUN/$S/infer_attempts.json"
    publish "$TMP/infer_$S" "$RUN/$S/infer"
    rm -rf "$TMP/controls_$S" "$TMP/$S.h5ad"
  fi
  if [ ! -d "$RUN/$S/score" ]; then
    log "scoring $S"
    EXTRA=(--both-with "$RUN/$(other "$S")/truth" --total-shift "$RUN/$S/truth/total_shift.tsv")
    [ "$S" = k562 ] && EXTRA+=(--essential-from "$RUN/rpe1/truth")    # check P5: K562 strata only
    $PY benchmarks/knockdown_score.py score --infer "$RUN/$S/infer" --truth "$RUN/$S/truth" --out "$TMP/score_$S" \
      --perms "$PERMS" --seed 0 --collectri "$TFS" "${EXTRA[@]}"
    publish "$TMP/score_$S" "$RUN/$S/score"
  fi
done

log "combine"
identity                                   # both score folders belong to this run's truth and inference
$PY benchmarks/knockdown_score.py combine --design "$RUN/rpe1/score" --replication "$RUN/k562/score" --out "$TMP/combined"
publish "$TMP/combined" "$RUN/combined"
log "end"
