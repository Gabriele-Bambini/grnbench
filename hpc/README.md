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

ARACNe3 is a native C++ program (CMake + OpenMP), not a jar. The adapter was written against
upstream commit `3d8791a`. Download on the login node; compile inside a job (no compute on the login node).
The cluster has GCC 11.5 with OpenMP but no cmake (no binary, no module; checked 2026-10-18), so the
download includes the official Kitware binary:

```bash
git clone https://github.com/califano-lab/ARACNe3 ~/tools/ARACNe3
git -C ~/tools/ARACNe3 checkout 3d8791a
# googletest and the whole Boost superproject (164 repositories; ARACNe3 uses only Boost.Math)
git -C ~/tools/ARACNe3 submodule update --init --recursive
V=3.31.8; cd ~/tools
curl -fsSLO https://github.com/Kitware/CMake/releases/download/v$V/cmake-$V-SHA-256.txt
curl -fsSLO https://github.com/Kitware/CMake/releases/download/v$V/cmake-$V-linux-x86_64.tar.gz
grep " cmake-$V-linux-x86_64.tar.gz$" cmake-$V-SHA-256.txt | sha256sum -c - && tar xzf cmake-$V-linux-x86_64.tar.gz
# from the repo folder: the build in a job (configure, build, check the program, print its usage line)
bash hpc/hpc.sh submit "--export=ALL,CMAKE=\$HOME/tools/cmake-3.31.8-linux-x86_64/bin/cmake hpc/aracne3_build.sbatch"
```

GitHub refuses anonymous clones after a burst. On 2026-10-18 the submodule command stopped after about
75 repositories with `could not read Username for 'https://github.com'`; 20 minutes later the rest
cloned one at a time with no refusal, and the same command then checked all of them out at their
pinned commits:

```bash
cd ~/tools/ARACNe3/deps/boost
for p in $(git submodule status | awk '$1 ~ /^-/ {print $2}'); do git submodule update --init -- "$p"; sleep 15; done
git -C ~/tools/ARACNe3 submodule update --init --recursive
```

`hpc/aracne3_build.sbatch` refuses a tree with a submodule missing or off its pin. On the cluster the
build job took about 27 minutes on 4 CPUs, most of it CMake's configure step over the Boost
superproject (record `experiments/EXP-20261018-GB-1.md`).

Set `GRNBENCH_ARACNE3_BIN` only in the job or shell that should run ARACNe3, not in your login profile:
`aracne3` is a default method, so wherever it is available `--methods available`, `grnbench run` and
`grnbench.sbatch` (`METHODS=available`) include it; `hpc/aracne3_smoke.sbatch` sets it inside its job.
That smoke run (EXP-20261018-GB-1) is the first run of the adapter against the real program: it works
as written. Each call is one ARACNe3 run in a fresh temporary folder with `--threads 1`. With more
threads the null model does not repeat from the seed (why: `src/grnbench/infer/aracne3.py`): there, two
4-thread runs with the same seed built different null models, although the network of that small subset
came out byte-identical. So `aracne3` does not take the engine's generic thread count: `grnbench infer`,
`grnbench run` and `grnbench.sbatch` pass `--threads` (4, or `SLURM_CPUS_PER_TASK` in the job) to every
method, and ARACNe3 still gets `--threads 1`. Only the aracne3-specific parameter `aracne3_threads` of
`run_method` changes that; no command-line flag sets it. The run's provenance says what ARACNe3 got:
`method_info.threads` (plus `threads_not_used` when a generic count was set aside), while `params.threads`
is what the caller asked of every method.

On the home's network file system (BeeGFS) a temporary folder ARACNe3 has just written can refuse to be
deleted (`Device or resource busy`, EXP-20261019-GB-1). The adapter then keeps the network, names the folder
left behind in the run's provenance and warns (a job of `hpc/engine_bench.sbatch` hides the warning: it sets
`PYTHONWARNINGS=ignore`), but the folder stays until you remove it from `TMPDIR`. Where the node's local disk
is big enough, set `TMPDIR=/tmp/...` inside the job, as `hpc/aracne3_smoke.sbatch` does. Slurm does not
report that size here (`TmpDisk=0` on the stud nodes, checked 2026-10-21), and in a bench job `TMPDIR` also
holds the engine's working folders.

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
