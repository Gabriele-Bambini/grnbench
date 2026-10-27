# Map - grnbench (where things are)

Keep this short and current: update it in any PR that adds, moves or removes a module, script or folder.
Updated: 2026-10-23.

## Run it
- Tests: `python -m pytest -q -p no:cacheprovider` (floor in AGENTS.md section 4). Golden values frozen
  from main for A/B tests, with the script that wrote them: `tests/golden/`.
- Board (session start): `python -m grnbench.board`; `--digest` weekly digest, `--claims` claims status table.
- GUI: `grnbench ui` (Streamlit, `.claude/launch.json` config `grnbench-ui`). Demo data: `grnbench demo --out d`.
- One dataset end to end: `grnbench run --data d.h5ad --out runs/x --gold dorotheaABC`.
- Every command also runs as `python -m grnbench <command>` (`src/grnbench/__main__.py`).
- Optional tools (R and its packages, Python and command-line extras): `grnbench env`. R env:
  `requirements/r-bridge.yml` (Linux) + `requirements/r-bridge.R`; guide `docs/R_BRIDGE.md`.
- Benchmarks: `benchmarks/engine_bench.py` (shared/context-specific planted networks,
  `--stress-grid`, immutable task/run manifests, terminal-status completeness and breaking points),
  `benchmarks/real_bench.py` (real data vs CollecTRI halves; `--models`, `--keep-artifacts`;
  reference hashes in `RUN_CONFIG.json`; `report` also writes figures + LaTeX),
  `benchmarks/count_plausibility.py` (would ingest's non-raw-counts warning fire on a data root; #42),
  `benchmarks/compare_confirm.py` (a candidate-study rerun's registered primaries vs `CONFIRM.json`; #95),
  `benchmarks/diff_rows_nested.py` (where two `rows/` folders differ at any depth, timings apart; diagnosis, #95).
- Dataset pools: `benchmarks/fetch_cellxgene.py --one-per-collection --max-cells --exclude-collections` (one dataset per study),
  `benchmarks/split_pool.py` (seeded confirmatory / exploratory split; D-057, `experiments/EXP-20261022-GB.pool.tsv`),
  `benchmarks/used_datasets.py` (datasets and studies used by any record or result on any branch),
  `benchmarks/subgroup_breakdown.py` (descriptive per-subgroup mean difference: organism, assay family, cell vs nucleus).
- Motif-prior lift decomposition: `benchmarks/motif_sensitivity.py` (`plan`, `run`, `report`) from a
  `--keep-artifacts` run (EXP-20261003-GB-1); tests `tests/test_motif_sensitivity.py`.
- ARACNe3 on the HPC: build `hpc/aracne3_build.sbatch`, adapter smoke run `hpc/aracne3_smoke.sbatch` ->
  `benchmarks/aracne3_smoke.py` (EXP-20261018-GB-1; tests `tests/test_aracne3_smoke.py`). Never on a laptop.
- C4's EPR comparison re-read from its committed rows, with planning power: `benchmarks/epr_power.py`
  (EXP-20261020-GB-1); tests `tests/test_epr_power.py`.
- Second-person rerun outputs: `benchmarks/results/EXP-20261003-LM-1-repro-GB/` (#95; record section 4b).
  `benchmarks/diagnose_rerun.py` (where a rerun's rows differ from the committed rows; exploratory, #24).
- Second-person rerun outputs: `benchmarks/results/EXP-20261002-LM-1-repro-GB/` (#94; record section 4b).
- Second-person rerun outputs: `benchmarks/results/EXP-20260926-LM-1-repro-GB/` (#93; record section 4b).
- Second-person rerun outputs: `benchmarks/results/EXP-20260925-LM-2-repro-GB/` (#24; record section 4b).
- Figures: `grnbench viz <run_dir>`. LaTeX: `grnbench report <run_dir|EXP-ID>`, `grnbench report --paper`.
- HPC: `hpc/hpc.sh`, VPN `hpc/vpn.ps1` (Windows), guide `docs/HPC.md`, agent procedure `.claude/skills/grnbench-hpc/SKILL.md`.
  Step-by-step launch runbooks per experiment issue: `docs/runbooks/` (#24 executed 2026-10-08, #52 launched 2026-10-10).

## Code (`src/grnbench/`)
| area | modules |
|---|---|
| ingest + preprocessing | `ingest` (formats -> counts AnnData), `pipeline` (`run_prep`, every step in `manifest.json`), `retention` (ledger: cells, genes and named TFs each step removed; warns when a gene filter drops a detected TF), `qc` (loud gates), `hygiene` (doublets, metadata), `features` (HVG with TF union), `strategies` (pseudobulk, metacells, ...) |
| cleaning engine | `kb` (stages, options, evidence, references), `executors` (`PLAN_STEPS`: which options grnbench can run, which are recorded only; refuses imputation, D-002), `context` (what the data are), `profile`, `selector` (rules), `autoprep` (plan, execute), `screen` (candidate pipelines on the data, reference split) |
| R bridge (optional) | `rbridge` (find an Rscript, cached probe, `run_r` / `run_external`, `requirement_status`, `environment_report` / `environment_status`), `r/_wrapper.R` (seed, RNG kind, `session.json`), `r/_probe.R` (versions as one JSON line); `kb` availability asks it for `r:` packages when `GRNBENCH_R_BRIDGE=1` (`kb.r_bridge_on`) |
| inference | `infer/` (`base` contract + registry, `correlation`, `classical`, `tree_methods` GENIE3/GRNBoost2, `aracne3`, `regdiffusion_adapter`, `cistarget` motif prior / pruning, `association` bicor/WGCNA/propr/ARACNe-style/MRNET, `regression` TIGRESS/BBSR/pcNet/extra-trees/boosting, `community` rank aggregation), `plugins` (bring your own model); what is built in and what is not: `docs/models_catalog.md` |
| evaluation + statistics | `eval/` (`metrics`, `gold`, `universe`, `nulls`, `stability`, `comparison` = `compare_units` paired t + Holm, `power` = power and needed n of that t, for planning) |
| evaluation + statistics | `eval/` (`metrics`, `gold`, `universe`, `nulls`, `stability`, `comparison` = `compare_units` paired t + Holm) |
| reporting | `diagnostics` (portable schema, strict reader and shared normalized report rows), `report` + `html_report` (per run), `viz` + `viz.mplstyle` (figures), `texreport` (LaTeX), `board` (team board) |
| seeds | `seeds` (`DEFAULT_SEED`, `derive` named child seeds, `SeedPlan`, `LEGACY` = the seeds each entry point uses today, read back from the code by `tests/test_seeds_legacy.py`; nothing calls it yet, #126) |
| reporting | `report` + `html_report` (per run), `viz` + `viz.mplstyle` (figures), `explore` (exploratory TF atlas, figures and evidence cards from run artifacts), `texreport` (LaTeX), `board` (team board) |
| interface | `ui/app.py` (Streamlit), `ui/logic.py` (tested logic), CLIs `cli_*.py` dispatched by `cli_main.py` |
| reporting | `diagnostics` (portable schema shared by prep/screen/models), `report` + `html_report` (per run), `viz` + `viz.mplstyle` (figures), `texreport` (LaTeX), `board` (team board) |
| interface | `ui/app.py` (Streamlit pages: Quick run, steps 1-7, Explore data, Explore results), `ui/steps.py` (the list of screens and their headlines), `ui/logic.py` (tested logic and plain-language summaries), `ui/catalog.py` (model families, speeds, presets), `ui/style.py` + `ui/theme.py` (stylesheet and the Streamlit options `grnbench ui` passes), CLIs `cli_*.py` dispatched by `cli_main.py` |
| simulation | `simulate` (shared + per-cell-type planted truth and degradations), `demo` |

## Records and team memory
- Direction and rules: `AGENTS.md`. Purpose: `.projectos/PURPOSE.md`. Decisions: `.projectos/DECISIONS.md`.
- Where things left off: `.projectos/STATE.md`. What is running now: GitHub issues labelled `running`.
- Runs: `experiments/` (one record per run, `datasets.tsv` provenance). Results: `benchmarks/results/<EXP or legacy folder>/`.
- `benchmarks/results/EXP-20261003-GB-1/`: motif-lift decomposition (stage 1): rows, `SENSITIVITY.md/json`,
  `RUN_INPUTS.json`, `INPUT_SHA256SUMS.*`, `perm_seeds/`, `SHA256SUMS`.
- `benchmarks/results/EXP-20261001-LD-1/synthetic/`: complete synthetic accounting, portable
  manifests, row table, generated stress-curve figures and breaking-point tables.
- `benchmarks/results/EXP-20261001-LD-1/real/`: complete 40-task accounting, portable task/data
  manifests, generated report, global visualization summary and LaTeX source; the full per-dataset
  visualizations and 144 MB of artifacts remain on HPC.
- Reproductions: section 4 of the record; `benchmarks/compare_runs.py COMMITTED RERUN --record <record>` compares a rerun with the committed results.
- Paper: `paper/` (claims ledger, checklist, CRediT log, `make_figures.py`, `build/` LaTeX draft).
- Non-run attempts (debugging, libraries tried): `.projectos/ATTEMPTS.md`. History: `docs/history/`, `audit/`.

## Data (never in git)
- Local: a folder outside the repo (e.g. `~/grnbench_data`), filled by `hpc/fetch_engine_datasets.sh`
  and `benchmarks/fetch_cellxgene.py`. Motif files: `GRNBENCH_CISTARGET_DIR`. Cluster: `~/grnbench_data`.
- Knockdown screens (Replogle 2022, scPerturb Zenodo 7041849): `python benchmarks/fetch_perturb.py --out
  <your data root>/perturb` on the login node; pinned md5, parallel byte ranges, resumable. Files land
  directly in `--out` (nothing is appended); `perturb/` is the prefix `data_manifest.py` knows.
