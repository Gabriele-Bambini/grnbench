# Running grnbench on the cluster

CPU-only: arboreto (GENIE3/GRNBoost2) parallelizes over cores with dask, and
nothing in the pipeline uses a GPU. Ask for cores and memory, not `--gres=gpu`.

Defaults in these scripts match the cluster already used by the sibling projects
(`--partition=stud --qos=stud`, conda at `/software/miniconda3`). Override with
`CONDA_SH=` / `ENV_NAME=` if yours differs.

Team setup (your own account, ssh alias, `.env.hpc`, VPN, lockouts, QOS limits):
**`docs/HPC.md`**. Short version:

```bash
bash hpc/hpc.sh sync                          # copy the working tree to ~/grnbench (or GRNBENCH_HPC_DIR)
bash hpc/hpc.sh run 'bash hpc/setup_env.sh'   # once: conda env "grnbench" from requirements/lock-py311.txt
```

## Submit

```bash
sbatch --export=ALL,DATA=$HOME/data/pbmc.h5ad,TFS=$HOME/data/tfs.txt,GOLD=dorotheaABC \
       hpc/grnbench.sbatch
```

Useful extra variables:

| variable | meaning |
|---|---|
| `OUT` | run directory (default `~/grnbench/runs/$SLURM_JOB_ID`) |
| `METHODS` | `available` (default) or e.g. `genie3,grnboost2,spearman_signed` |
| `PLUGINS` | **semicolon**-separated `--plugin` arguments, e.g. `my_model.json;other.py:Cls`. Not colon-separated: every plugin reference (`pkg.mod:Class`, `file.py:Class`) already contains a colon. |
| `STABILITY` | `1` (default) or `0` to skip subsample stability. It re-runs every method N times — on a big dataset that dominates the wall clock. |
| `STABILITY_N` | subsamples per method (default 5) |
| `NHVG` | HVG count (default 3000) |
| `GROUPKEY` | obs column for pseudobulk / strategy recommendations |
| `ENV_NAME` | conda env name (default `grnbench`) |

## Read the result

```bash
cat  runs/<jobid>/eval/eval_report.md
scp  USER@CLUSTER:~/grnbench/runs/<jobid>/eval/eval_report.html .
```

The HTML report is fully self-contained — no CDN, no JavaScript — so it renders
from a scratch directory on a node with no outbound network, and survives `scp`
as a single file.

## Sizing

| dataset | methods | reasonable ask |
|---|---|---|
| ≤5k cells × 3k genes | spearman + plugins | `--cpus-per-task=8 --mem=32G --time=2:00:00` |
| ~20k cells × 3k genes | + genie3/grnboost2 | `--cpus-per-task=8 --mem=64G --time=12:00:00` (the default; 8 CPUs is the QOS cap) |
| >50k cells | + genie3/grnboost2 | pseudobulk first (`GROUPKEY=cell_type`) or `--subsample-cells` |

Tree methods scale roughly with `n_genes × n_TFs`; restricting regulators with
`TFS=` is the single biggest speedup, and it is also the biologically correct
thing to do.

## ARACNe3

```bash
export ARACNE3_JAR=$HOME/tools/aracne3.jar   # + java on PATH
```

An explicitly requested unavailable method produces a nonzero exit status;
successful method artifacts remain available with a diagnostic manifest.

## BEELINE-derived backtest

`fetch_beeline.sh` verifies archives from [Zenodo v4](https://zenodo.org/records/19009603),
the data release for BEELINE v1.1. `beeline_backtest.sbatch` records that declaration.
Selection follows [pinned official code](https://github.com/Murali-group/Beeline/blob/37464085eb8a95d6cc6a3d3a3c649d36db6052ed/utils/generateExpInputs.py):
Bonferroni-filtered p-values, significant TFs and top N significant non-TFs by variance.
Case-insensitive TF matching is a documented extension.

Evaluation uses grnbench's explicit TF-by-gene universe and average precision;
this is not a reproduction of every historical BEELINE metric convention.
All runs require an empty output directory. Per-method edges, provenance,
selected expression, effective gold and candidate universe are retained.
Missing inputs, unavailable requested methods and partial errors return nonzero.
Setup fails on installation or test errors; it never reports `SETUP_DONE` after them.
