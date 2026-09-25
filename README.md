# grnbench — a plugin-first GRN benchmark with bulletproof preprocessing

> Team: read [`AGENTS.md`](AGENTS.md) before working. New member or new agent: [`docs/START_HERE_AGENT_PROMPT.md`](docs/START_HERE_AGENT_PROMPT.md).
> Runs are registered in [`experiments/`](experiments/README.md); paper material lives in [`paper/`](paper/README.md).

## New in 0.8 (in development): interface + context-aware cleaning

```
grnbench ui            # graphical interface in the browser (pip install 'grnbench[ui]')
grnbench demo --out d  # synthetic dataset with a known network to try it on
grnbench kb            # the preprocessing knowledge base: every option, availability
```

The interface walks one dataset from raw counts to a model comparison:

1. **Data** - any format, plus README files, download pages (GEO/Zenodo/...) and notes.
2. **Profile** - what the data are: protocol (from text *and* from the count values),
   organism, batches, cell types, batch x cell-type confounding, QC distributions.
3. **Pipeline** - one choice per decision point (14 stages, 15 decisions, 56 options), chosen
   by rules from a cited knowledge base (`docs/knowledge/preprocessing_kb.md`); every choice
   shows its reasons, evidence level and references and can be changed.
4. **Compare pipelines** - the recommended pipeline and alternatives are run and measured
   (retention, depth/batch leakage, edge stability, recovery of the *selection half* of a
   reference network); pros and cons; the data pick the final pipeline.
5. **Clean** - one run per cellular context, with the exact `grnbench prep` command.
6. **Models** - built-in methods and **your own model** (a `.py` file), same cleaned input.
7. **Benchmark** - scores on the *evaluation half* of the reference, never used to choose
   the pipeline.

The same engine is scriptable: `grnbench.autoprep.diagnose`, `plan_preprocessing`,
`execute_plan`, `grnbench.screen.screen_candidates`.

**v0.7.0 — the ranking is not real, and neither was the test that said so.** One command takes raw scRNA-seq counts to a
comparison report: `preprocess → infer → evaluate → report`, every step logged,
every method interchangeable. Positioning and design notes:
`docs/GRN_PLATFORM_PLAN.md` in the parent project.

```
grnbench run --data d.h5ad --out runs/today --tf-list tfs.txt \
             --gold dorotheaABC --plugin my_model.json --stability
```

```
any format ──▶ counts-verified ingest ──▶ QC gates (fail loudly) ──▶ filters
            ──▶ normalization strategy (auto-recommended) ──▶ HVG ∪ TF features
            ──▶ processed.h5ad + qc_metrics.json + manifest.json + report.md

processed.h5ad ──▶ inference adapters ──▶ <method>.edges.tsv.gz  [source,target,weight(,direction)]
                   built-in: genie3, grnboost2, spearman_signed, aracne3*
                   yours:    Python entry point OR Docker/CLI image

edges + gold ────▶ evaluation ─────────▶ eval_report.json/.md/.html
                   AUPRC, AUPRC-ratio, EPR, AUROC, PR curves,
                   score-permutation nulls, cross-method overlap,
                   subsample stability, gold reachability
(*aracne3 needs java on PATH + ARACNE3_JAR; reported as unavailable otherwise)
```

## Why this and not BEELINE, or geneRNIB

BEELINE (Pratapa et al. 2020) is the canonical GRN benchmark and its evaluation
code is sound — its early-precision handles score ties correctly and its AUPRC
delegates to scikit-learn. But its method list is fixed, its glue is Docker+R,
and it has no story for "score *my* model".

[geneRNIB](https://github.com/openproblems-bio/task_grn_inference) (Open
Problems, 2025) is the live successor: 11 datasets, 12 methods, 11 metrics, a
public leaderboard and CI that re-scores everything when a new component lands.
If what you want is a leaderboard, use geneRNIB. grnbench will not out-scale it
and does not try to.

grnbench asks the question a leaderboard cannot ask about itself:

> **How much of this ranking is real?**

A leaderboard asserts a total order over methods. That assertion goes untested,
and on real data it does not hold up. Three things this harness reports that a
ranking cannot:

- **Null-calibrated scores.** Every metric is also expressed against a
  score-permutation null on the *same* fixed candidate universe (global or
  within-regulator), with a conservative `(hits+1)/(n+1)` p-value. A method
  can lead the table and still be indistinguishable from its own null; both
  facts get printed.
- **Resolved vs unresolved comparisons — with the right units.** `grnbench
  compare` runs a paired t test with Holm correction over *independent
  experimental units* — datasets, not gold edges or cells — and refuses to
  call an unresolved pair "equivalent". On 9 methods × 7 BEELINE datasets
  (six classical, plus GENIE3, GRNBoost2 and a scikit-learn random forest) it
  resolves **0 of 216** pairs — none on EPR, none on AUPRC-ratio, on any
  reference; GENIE3 is not distinguishable from Pearson correlation. An
  earlier version of this package resampled gold edges instead and reported
  80% resolved; that was pseudo-replication and is retracted (`CHANGELOG.md`,
  0.7.0).
- **Ordering differences across ground truths, reported as differences.**
  Scored against every ground truth an organism has, the point-estimate
  winner changes between references on 7 of 7 BEELINE datasets. The harness
  says so — and also says that none of those reversals is resolved at n = 7,
  so they are not evidence of contradiction, only of a design that cannot
  rank.

Plus the parts that make those numbers trustworthy: a stable edge schema, a
two-line plugin contract, counts-verified preprocessing with hard QC gates, and
a report that tells you when a good-looking AUPRC is an artifact.

It is also checked *against* BEELINE rather than only compared to it:
`benchmarks/beeline_backtest.py` replicates the BEELINE protocol on BEELINE's
own corpus, and the pass condition is that grnbench reproduces the published
regime (AUPRC-ratio ~1.10, EPR ~1.20 against cell-type-specific ChIP-seq) — not
that it beats it. An EPR of 20 on real data means the harness is leaking, not
that biology was discovered. 83 of 84 cell-type-specific results land inside
3× that regime (125 of 126 with the tree methods), and the one that does not
is eight gold hits against two expected on a 756-edge reference — the report
says so in the same row; see `benchmarks/README.md`. The run has so far found eight defects in grnbench
itself, each listed in `CHANGELOG.md` with a regression test.

## Install

```bash
cd grnbench
python -m venv .venv && .venv\Scripts\activate     # or: uv venv .venv
pip install -e ".[dev]"                             # core + tests
pip install -e ".[all]"                             # + arboreto, decoupler, pyarrow, matplotlib
pytest                                              # 294 tests, synthetic data only
```

Extras: `infer` (arboreto/GENIE3/GRNBoost2 + dask), `dorothea` (DoRothEA gold
standards via decoupler), `plugins` (pyarrow, for parquet plugin exchange),
`plot`, `all`.

> On network drives where hardlinking fails (Google Drive, some NFS mounts),
> use `uv pip install --link-mode=copy ...` or plain `pip`.

## Quickstart

### One command

```bash
grnbench run --data data.h5ad --out runs/mydata \
  --tf-list tfs.txt --n-hvg 3000 --group-key cell_type \
  --gold gold_pairs.tsv --stability
```

Produces `runs/mydata/{prep,infer,eval}/` plus `run_manifest.json`; open
`runs/mydata/eval/eval_report.html`.

### Or stage by stage

```bash
grnbench prep  data.h5ad --out runs/mydata --tf-list tfs.txt --n-hvg 3000
grnbench infer runs/mydata/processed.h5ad --out runs/infer --tf-list tfs.txt
grnbench eval  runs/infer/*.edges.tsv.gz --gold gold_pairs.tsv --out runs/eval \
               --stability-data runs/mydata/processed.h5ad
grnbench methods       # what can run here, and why the rest cannot
```

Python API: `from grnbench import run_prep`; `from grnbench.infer import
run_method, register_method`; `from grnbench.eval import evaluate_ranking,
null_metrics, subsample_stability`; `from grnbench.html_report import render_html`.

## Bring your own model

A plugin is one object with one method. `adata` arrives already normalized
(`grnbench.infer.prepare_adata`), so `obs`, `var` and `layers["counts"]` are all
reachable.

```python
class MyGRNModel:
    name = "my-model"
    def infer(self, adata, tf_list, **params) -> "pd.DataFrame":
        # -> columns [source, target, weight]  (+ optional direction ∈ {-1,0,1})
        ...
```

Three ways to plug it in:

| how | when | command |
|---|---|---|
| entry point | you ship a package | `[project.entry-points."grnbench.methods"]`<br>`my-model = "my_pkg.models:MyGRNModel"` → auto-discovered |
| file reference | a one-off script | `--plugin path/to/model.py:MyGRNModel` |
| spec file | Docker / any CLI | `--plugin my_model.json` |

Container/CLI spec — grnbench writes `expr.parquet` (cells × genes), `tfs.csv`
and `params.json` into a work dir, runs your command, reads back `edges.parquet`:

```json
{
  "name": "my-docker-model",
  "format": "parquet",
  "command": ["docker", "run", "--rm", "-v", "{workdir}:/data",
              "my-registry/my-model:1.0",
              "--expr", "/data/{expr_name}", "--tfs", "/data/{tfs_name}",
              "--out", "/data/{out_name}", "--threads", "{threads}", "--seed", "{seed}"]
}
```

Placeholders: `{workdir} {expr} {tfs} {out} {params}` (absolute host paths) and
`{expr_name} {tfs_name} {out_name} {params_name}` (basenames as seen inside a
container mount), plus `{threads}` and `{seed}`. Set `"format": "tsv"` where
pyarrow is unavailable.

Working examples ship in [`examples/plugins/`](examples/plugins/):
`ridge_model.py` is a real scikit-learn method (not a stub) usable as
`--plugin examples/plugins/ridge.json`.

Plugins go through the same registry, edge schema, manifests and evaluation as
the built-ins — and cannot silently shadow a built-in name.

## Data-hygiene coverage (what breaks in real datasets → what we do)

| real-world problem | grnbench behavior |
|---|---|
| broken `layers["counts"]` (transposed, negative, non-integer) | validated at ingest; `IngestError` with the exact defect |
| NaN/Inf anywhere in the matrix | rejected loudly (no silent cleaning) |
| cells with zero total counts | QC gate + dropped before normalization (logged) |
| duplicate gene symbols (incl. CellRanger's own, e.g. `TBCE`) | aggregated by sum on **every** loader path — `.h5ad`, 10x `.h5`, MTX, loom, tables |
| messy labels (" Tumor ", "tumor", "TUMOR", "?", "unknown") | `--label-key --label-mapping`: whitespace/case/synonym normalization with per-change report |
| missing `--group-key` labels | rows dropped before pseudobulk/stratification, logged |
| missing `--label-key` labels | rows kept, counted in a warning, **excluded from stratified subsampling** (never a literal `"nan"` stratum) |
| **class imbalance** (max/min ratio, tiny classes) | `class_imbalance` gate on `--label-key` **and** `group_imbalance` on `--group-key` when both are given; stratified subsampling via `--stratify-key` |
| batch effects | auto-detection of batch-like columns + presence gate; `--strategy cpm_log1p_combat` corrects on request, never by default |
| doublets | heuristic gate always; `--drop-doublets` calls scrublet (per batch) and drops them, **refusing above a 20% call rate** (`--max-doublet-rate`) because that is a failed threshold, not a doublet population |
| duplicated barcodes | made unique with suffixes, originals kept in `obs['barcode_as_loaded']`, gate reports the count |
| Ensembl-only / mixed gene IDs | verdict gate: warns that TF matching will fail |
| mouse vs human symbols | organism hint recorded; mismatch warning vs TF list |
| TF list barely matching the data (<30%) | explicit warning in report + manifest |
| huge dense matrices | auto-sparsify at >=80% zeros (CSR) |
| wasteful float64 | downcast to float32 (logged step) |
| pseudobulk micro-groups | kept+flagged by default, or dropped loudly with `--min-group-cells`; a strategy that empties the matrix is refused, not written |

## Hard rules (non-negotiable by design)

1. **Counts or nothing.** Ingest rejects anything that does not look like raw
   integer counts (normalized/log inputs raise `IngestError`).
2. **No imputation.** The strategy registry cannot contain smoothing/imputation
   methods — a test enforces this.
3. **HVG never drops TFs.** Kept features = `HVGs ∪ TFs detected in data`.
4. **QC gates fail loudly.** Extreme mito fraction / low complexity abort the
   run (exit code 2) unless `--force`, which is recorded in the manifest.
5. **Full provenance.** Steps, params, thresholds, environment versions and an
   input fingerprint land in `manifest.json` / `inference_manifest.json` /
   `run_manifest.json`.
6. **Unavailability is data, not a crash.** Methods with missing external deps
   (arboreto, java/jar, a Docker image, a broken plugin import) stay listed with
   their reason and are skipped loudly.

## Inference adapters (canonical edge schema)

| method | engine | notes |
|---|---|---|
| `spearman_signed` | numpy/scipy | deterministic baseline; emits sign as direction |
| `genie3` | arboreto RF | dashboard-less dask client injected |
| `grnboost2` | arboreto GBM early-stop | same |
| `aracne3` | Java jar | parser unit-tested vs published format; runner probes java+ARACNE3_JAR |
| *yours* | Python or container | see **Bring your own model** |

Implementation note: this arboreto release crashes with modern dask when its
wrappers run with default settings (empty meta delayed-list); we drive
`arboreto.core.create_graph(include_meta=True)` directly and discard the meta.
On HPC nodes with Java, export `ARACNE3_JAR=/path/aracne3.jar` to enable ARACNe3.

## Evaluation

- **Metrics** — AP/AUPRC, AUPRC-ratio vs prevalence, precision@k, EPR
  (k = #gold hits by default), tie-safe rank-based AUROC, PR curves.
- **Nulls** — score permutation on the fixed candidate universe (`random`),
  within each regulator (`regulator`), or an uncalibrated degree-rewiring
  diagnostic (`degree`, never a p-value); conservative `(hits+1)/(n+1)`
  p-values; degenerate nulls counted, never hidden.
- **Comparison** — `grnbench compare scores.csv --unit-description ...
  --assume-independent-units`: paired t + Holm over declared independent
  units, Bonferroni simultaneous intervals; `unresolved_is_equivalent: False`.
- **Gold-thinning sensitivity** — `grnbench-eval --gold-thinning-sensitivity`:
  descriptive annotation-removal scenarios at fixed `k`. Quantiles, not
  confidence intervals; no winner.
- **Cross-method overlap** — pairwise top-k Jaccard. Two methods with equal
  AUPRC and no shared edges are not two confirmations of the same biology.
- **Stability** — `--stability-data processed.h5ad` re-runs each method on cell
  subsamples and reports the mean pairwise Jaccard of top-k edges.
- **Gold reachability** — how much of the gold each method's candidate space
  could even contain, so a filter-imposed ceiling is not read as a bad algorithm.
- **Gold standards** — any file with regulator/target columns, or DoRothEA via
  the optional `decoupler` extra (`--gold dorotheaABC`).

## Preprocessing strategies

| name | recipe | notes |
|---|---|---|
| `cpm_log1p` | CPM/1e4 → log1p | safe default |
| `pearson_residuals` | analytic Pearson residuals | variance-stable across depths |
| `pseudobulk` | sum per group (+ optional cpm_log1p) | needs `--group-key`; tiny groups kept but flagged |
| `raw_counts` | pass-through | for methods normalizing internally |
| `cpm_log1p_combat` | cpm_log1p → ComBat on `--batch-key` | never auto-recommended; removes *all* batch-aligned covariance, including biology confounded with batch; values are residuals, can be negative |
| `scran_pooled` | — | planned (R interop) |
| `scvi_latent` | — | planned |

## Input formats

`.h5ad` (counts in X or `layers["counts"]`), 10x `.h5`, 10x MTX directory,
`.loom`, `.csv`/`.tsv`/`.txt` tables with automatic orientation detection.
Duplicate gene names are aggregated by sum (never silently dropped or renamed);
Ensembl version suffixes are stripped.

## HPC

`hpc/grnbench.sbatch` is a ready SLURM job: creates the venv, installs the
`all` extra, runs the full loop and drops the HTML report in the run directory.
See [`hpc/README.md`](hpc/README.md).

## Roadmap

- **v0.6** — SERGIO/BoolODE synthetic ground truth, `scran_pooled` (rpy2),
  results store (sqlite+parquet) with a dashboard over it, PIDC and CellOracle
  adapters.
