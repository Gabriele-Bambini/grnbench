# Map - grnbench (where things are)

Keep this short and current: update it in any PR that adds, moves or removes a module, script or folder.
Updated: 2026-10-15.

## Run it
- Tests: `python -m pytest -q -p no:cacheprovider` (floor in AGENTS.md section 4).
- Board (session start): `python -m grnbench.board`.
- GUI: `grnbench ui` (Streamlit, `.claude/launch.json` config `grnbench-ui`). Demo data: `grnbench demo --out d`.
- One dataset end to end: `grnbench run --data d.h5ad --out runs/x --gold dorotheaABC`.
- Every command also runs as `python -m grnbench <command>` (`src/grnbench/__main__.py`).
- Benchmarks: `benchmarks/engine_bench.py` (shared/context-specific planted networks,
  `--stress-grid`, immutable task/run manifests, terminal-status completeness and breaking points),
  `benchmarks/real_bench.py` (real data vs CollecTRI halves; `--models`, `--keep-artifacts`;
  reference hashes in `RUN_CONFIG.json`; `report` also writes figures + LaTeX),
  `benchmarks/count_plausibility.py` (would ingest's non-raw-counts warning fire on a data root; #42).
- Motif-prior lift decomposition: `benchmarks/motif_sensitivity.py` (`plan`, `run`, `report`) from a
  `--keep-artifacts` run (EXP-20261003-GB-1); tests `tests/test_motif_sensitivity.py`.
- Figures: `grnbench viz <run_dir>`. LaTeX: `grnbench report <run_dir|EXP-ID>`, `grnbench report --paper`.
- HPC: `hpc/hpc.sh`, VPN `hpc/vpn.ps1` (Windows), guide `docs/HPC.md`, agent procedure `.claude/skills/grnbench-hpc/SKILL.md`.

## Code (`src/grnbench/`)
| area | modules |
|---|---|
| ingest + preprocessing | `ingest` (formats -> counts AnnData), `pipeline` (`run_prep`, every step in `manifest.json`), `qc` (loud gates), `hygiene` (doublets, metadata), `features` (HVG with TF union), `strategies` (pseudobulk, metacells, ...) |
| cleaning engine | `kb` (stages, options, evidence, references), `context` (what the data are), `profile`, `selector` (rules), `autoprep` (plan, execute), `screen` (candidate pipelines on the data, reference split) |
| inference | `infer/` (`base` contract, `correlation`, `classical`, `tree_methods` GENIE3/GRNBoost2, `aracne3`, `regdiffusion_adapter`, `cistarget` motif prior / pruning), `plugins` (bring your own model) |
| evaluation + statistics | `eval/` (`metrics`, `gold`, `universe`, `nulls`, `stability`, `comparison` = `compare_units` paired t + Holm) |
| reporting | `diagnostics` (portable schema shared by prep/screen/models), `report` + `html_report` (per run), `viz` + `viz.mplstyle` (figures), `texreport` (LaTeX), `board` (team board) |
| interface | `ui/app.py` (Streamlit pages: Quick run, steps 1-7, Explore data, Explore results), `ui/steps.py` (the list of screens and their headlines), `ui/logic.py` (tested logic and plain-language summaries), `ui/catalog.py` (model families, speeds, presets), `ui/style.py` + `ui/theme.py` (stylesheet and the Streamlit options `grnbench ui` passes), CLIs `cli_*.py` dispatched by `cli_main.py` |
| simulation | `simulate` (shared + per-cell-type planted truth and degradations), `demo` |

## Records and team memory
- Direction and rules: `AGENTS.md`. Purpose: `.projectos/PURPOSE.md`. Decisions: `.projectos/DECISIONS.md`.
- Where things left off: `.projectos/STATE.md`. What is running now: GitHub issues labelled `running`.
- Runs: `experiments/` (one record per run, `datasets.tsv` provenance). Results: `benchmarks/results/<EXP or legacy folder>/`.
- `benchmarks/results/EXP-20261001-LD-1/synthetic/`: complete synthetic accounting, portable
  manifests, row table, generated stress-curve figures and breaking-point tables.
- `benchmarks/results/EXP-20261001-LD-1/real/`: complete 40-task accounting, portable task/data
  manifests, generated report, global visualization summary and LaTeX source; the full per-dataset
  visualizations and 144 MB of artifacts remain on HPC.
- Paper: `paper/` (claims ledger, checklist, CRediT log, `make_figures.py`, `build/` LaTeX draft).
- Non-run attempts (debugging, libraries tried): `.projectos/ATTEMPTS.md`. History: `docs/history/`, `audit/`.

## Data (never in git)
- Local: a folder outside the repo (e.g. `~/grnbench_data`), filled by `hpc/fetch_engine_datasets.sh`
  and `benchmarks/fetch_cellxgene.py`. Motif files: `GRNBENCH_CISTARGET_DIR`. Cluster: `~/grnbench_data`.
