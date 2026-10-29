# Changelog

## Unreleased (0.8.0-dev) — interface and context-aware cleaning

### Project
- `fetch_cellxgene.py --one-per-collection` keeps the smallest H5AD of each CELLxGENE collection (one
  dataset per study, `collection_id` in the manifest) `--max-cells` lifts the 50,000-cell cap and `--exclude-collections` skips whole
  studies, `--listing` selects from a saved listing; `benchmarks/used_datasets.py` lists every dataset
  and study any record or result on any branch used; `benchmarks/subgroup_breakdown.py` gives the
  declared descriptive breakdown of a paired difference by organism, assay family and cell vs nucleus; the
  defaults select exactly as before (pinned by a test). `benchmarks/split_pool.py` splits a pool into a
  confirmatory and an exploratory part by sha256 of a seed and the dataset id. D-057 (ACTIVE, GB's go of 2026-10-23): a pool of
  200 never-used datasets of every kind from studies never used before, 120 confirmatory and 80 exploratory; EXP-20261022-GB-1
  (confirmatory, EPR lift of the motif prior) and EXP-20261022-GB-2 (exploratory stress test) registered,
  not run.
- C5 and C6: GB signed both section 4b reruns in his own hand (comments on #99 and #100, 2026-10-23), replacing
  the lines the agent typed on his instruction; the CLAIMS cells say the commands were typed by an agent under
  his direction. C6's one differing field (`n_warnings` 4 -> 3) is confirmed to come from the dataset README
  written by the fetch route (#156). No number or status changes.
- D-022 accepted (ACTIVE, GB and LD on #79): the next experiments test networks against CRISPRi knockdowns; a single-screen experiment may take the TF as its unit, with a dependence-keeping null and a replication set. TP53 in breast cancer is parked as an exploratory example of use, run after the engine plan.
- Symmetric statistics (Pearson, Spearman, mutual information, partial correlation) now give bit-identical weights to A->B and B->A. `numpy.corrcoef` divides the covariance by the two standard-deviation vectors one after the other, so reciprocal entries could differ in the last bit (on a seeded 300 x 60 matrix, 581 of 1,770 pairs; `np.cov` itself is symmetric), and tie-aware AUPRC then split or merged 2-edge tie blocks differently on Windows and Linux: about 1.5% of the AUPRC ratio on the measured pair of #82, around its third decimal. Registered results reproduce at their recorded commit (`CODE_COMMIT`); new runs of the same data can differ from them at about that size.
- Screening and benchmark scoring subsample cells by a hash of the barcode, not by position (D-024):
  a few cells more or fewer upstream (scrublet differing across platforms) no longer reshuffle the
  3,000-cell sample or flip the screened choice. Two effects on new runs: the cap changes only above
  3,000 cells; the 150-cell halves of `edge_stability` change on every sufficiently large matrix, below
  the cap too. `engine_bench` and `real_bench` pin `cell_sampling` (`blake2b64-seed-barcode-v1`) in
  `RUN_CONFIG.json` (`real_bench` also `max_score_cells`) and refuse to resume an output without it.
- `fetch_cellxgene.py --from-manifest` (#156) now writes the same `README.md` as the listing route: it takes
  the collection name and citation from the CELLxGENE listing by dataset id, and warns on stderr, per dataset,
  when it cannot (the download goes on). The empty README it wrote before was the cause of C6's
  `n_warnings` 4 -> 3 on rerun (EXP-20261002-LM-1, section 4b); no score was affected.
- Retention ledger (#176): every cell- or gene-removing step of `run_prep` (cell calling, unlabelled
  groups, cell and gene filters, QC-gene exclusion, doublets, context subset and its gene filter,
  subsampling, strategy, feature selection) records cells and genes before/after and, with a TF list,
  the TFs of that list before/after and the ones it removed, named (first 20, sorted, full counts),
  in `manifest.json` (`retention`) and in the portable diagnostics; `report.md`, the HTML and the LaTeX
  reports show it as a generated table. A gene filter that drops a TF detected in the context (non-zero
  in at least one of its cells) adds one `TF protection:` warning naming the TFs with the number of
  context cells each was detected in. Accounting only: the kept cells and
  genes are unchanged (A/B on every simulator preset, the demo and C6). Whether the gene filters should
  keep detected TFs is an open default decision. Portable diagnostics: `safe_message` is now
  idempotent (redaction repeats until the text is stable, at most 8 passes; a text still changing,
  such as a chain of 8 or more ` /N` ratio links, is stored as `<message withheld>`), so every message
  or gene name a writer stores passes the reader; before, a length cut inside an A/C/G/T run could
  make the reader refuse a document the writer had just written, and `run_prep` fail after writing
  its outputs (#177 review). Sensitive values and local user/host names are never matched inside a
  placeholder (`--out path` no longer nests `<path>`), and the schema's own field names (an explicit
  list) are not compared with the local user/host names, every other key still is (#177 re-review).
- AGENTS.md section 5: two review tiers (tier 1 = no change to any registered number, default, claim, decision, record or `paper/` beyond a CONTRIBUTIONS row, one approval from a team member other than the author, an independent agent review from another member's side counts; tier 2 = everything else in `paper/`, result wording, records, `eval/` and statistics, output-changing defaults, decisions, with a reviewer from the "reviews" column and the 2-working-day fallback), a challenge route for any agent-reviewed or misclassified merge, any team member may merge an approved PR, and a daily merge train that runs the suite once per batch and checks each PR against the train before its squash (MAP and STATE resolved by reading, not by keeping both sides).
- Per-cell-type plans no longer apply ComBat or per-batch centring when the batch column is the
  cell-type context column: every resulting context has one batch, so the executable plan now
  keeps the requested choice as provenance, runs the selected normalization without the impossible
  correction and explains the fallback. A CLI `--n-hvg` request above the available gene count now
  says in the terminal, manifest and report that every remaining gene was kept (#164).
- Silent behaviours made loud, nothing else changed (#149 items 1, 3, 4, 5, 6; #116): Ensembl identifiers
  with no mitochondrial gene found by symbol now say that the per-cell mito filter is off and the mito gate
  always passes (plan reason and `run_prep` warning); a missing-label group such as "unknown" that becomes a
  per-cell-type context is named in the context reason (still inferred: excluding it is a decision); a
  context that ends below 200 cells after QC gets a `run_prep` warning (a planner context or a user's
  `--subset-key`/`--subset-value` alike: "fewer than 200 cells, the minimum the planner uses for a
  per-cell-type context"); `grnbench prep`/`run` without
  `--n-hvg` say that all filtered genes are used and that the planner defaults to 1,000 HVGs plus the TFs;
  an organism guessed from gene-name casing while `uns['organism']` says another is named in the
  regulators reason; and when the knee finds no cliff on an unfiltered matrix, `diagnose` records it, the
  plan warns and the cell-calling reason says so instead of "? of ? barcodes (>= nan counts)". A/B on
  every simulator preset and the demo: identical choices, values, contexts, warnings and processed
  matrices; only these notes and warnings are added.
  Review fix (LD): `--n-hvg all` is an explicit value (every filtered gene, no note; `grnbench run` too),
  a plan's recorded command always names `--n-hvg` (the number, or `all` when the plan keeps every gene)
  and `manifest['args']['note_all_genes']` records whether the note was asked for, so the replayed command
  of a plan writes the manifest steps and warnings `execute_plan` wrote, not only the same matrix.
- Tests now restore the method and plugin registries after every case, so a UI test that registers
  a user model cannot make a later plugin test fail only because of collection order (#183).
- Recorded per-context preprocessing commands now use ``--subset-value=<label>`` so labels that
  begin with ``-`` replay as data rather than being misread as command-line options (#172).
- A stratified subsample whose context contains only missing labels now fails that context with an
  actionable `PrepError` instead of dividing by zero and aborting the remaining contexts (#169).
- Preprocessing plans now tell the truth when a requested option falls back: pseudobulk names the
  ignored batch/normalization choices, per-batch MAD names its global fallback, small matrices name
  the all-gene HVG fallback, and context plans name pooled/largest-type behaviour. Random and
  stratified cell subsampling have an explicit runtime mode recorded in the manifest, and unknown
  decision keys are refused instead of ignored (#132).
- Report diagnostics (#59): Markdown, HTML and LaTeX reports now consume one validated schema-v1
  diagnostic view for preparation, screening and model outcomes, choices, failure categories and
  before/after dimensions. Historical runs are marked unavailable; malformed or unsupported
  diagnostics fail clearly instead of being guessed from other artifacts. Readers require the
  complete schema and the writer's bounded choice allowlist; one-shot reports fingerprint and
  render the final diagnostic rather than retaining the temporary partial artifact.
- `grnbench.eval.power`: power, needed number of datasets and smallest detectable effect of the paired t
  over datasets (G*Power matched-pairs model, noncentral t; one-sided = the directed test of
  `compare_units`), tested against the G*Power manual and a simulation. `benchmarks/epr_power.py` re-reads
  C4's EPR comparison from its committed rows (gate: the registered report code must reprint the committed
  lines) and writes a descriptive diagnosis and planning tables. EXP-20261020-GB-1 (exploratory, #127): the
  EPR lift of the motif prior in C4 was underpowered at n = 20; C4 is unchanged. After the review of #129
  the report also reruns the registered code with each dataset left out (the verdict changes in 0 of 20
  drops in the auto arm, 2 of 20 in the baseline arm), no longer reports power at C4's own estimate, and
  `n_for_power` never plans fewer than 3 datasets.
- Optional R bridge (#135, D-032 proposed): `grnbench.rbridge` runs R as an `Rscript --vanilla` subprocess
  from a dedicated conda env `grnbench-r` (`requirements/r-bridge.yml` on Linux, `requirements/r-bridge.R`
  on both systems, `docs/R_BRIDGE.md`), with seed and RNG kind fixed by `r/_wrapper.R` and R and package
  versions in `session.json`; the package probe runs once per cache. `grnbench env` (and the sidebar's
  Environment popover, which starts R only on its button) lists R, its packages and the other optional
  tools with the command to install what is missing. grnbench installs nothing; `GRNBENCH_RSCRIPT=none`
  turns R off. It never switches quietly to another R: an R you named, or a conda env that exists but
  does not work, is an error, and `grnbench env` exits 1 for it. The method's own `set.seed` calls and
  RNG changes are recorded; a timeout or Ctrl+C kills R with its child processes. Not wired into the
  knowledge base's availability yet (step 2c): planner output is unchanged.
- An option of the knowledge base is available only when its dependencies are installed AND grnbench can
  run it (`kb.DecisionPoint.availability`; new `executors.py`, whose `PLAN_STEPS` table both availability
  and the plan read, and which never offers imputation or scVI's denoised expression, D-002). An option
  grnbench cannot run is never chosen by the selector, and forcing one fails loudly with the reason and
  the install hint (D-031). Closes the audit finding "external options silently dropped" (#128, step 2a).
  Declared change, only where R, scvi or CellBender is installed: the selector no longer picks
  scDblFinder, EmptyDrops, DecontX & co., which the plan then ignored, and picks the best option grnbench
  can run (e.g. Scrublet). Without them (every registered run) plans, warnings and the C7
  `skipped_variants` are unchanged. The "best option ... is unavailable" warning names what actually
  runs; what it would take to run the best option is in `plan.option_notes` and in the interface. Not
  fixed here: regulator options are recorded only (#131); six other plan/run mismatches (#132).
- The knowledge base can see R packages through the R bridge (step 2c, #135, D-032): with
  `GRNBENCH_R_BRIDGE=1`, `kb` asks `rbridge.requirement_status` for `r:` requirements (one cached probe);
  off, the default, an R package counts as missing without starting R, so plans, `grnbench kb`, the
  interface and reruns of registered experiments are unchanged. The rpy2 path is gone. No R option has a
  grnbench executor yet, so with the bridge on an installed package only changes the reason to "no grnbench
  executor ... yet". `grnbench env` says whether R is on for planning.
- `test_grnbench_ui_hands_the_theme_to_streamlit` compares the app path after resolving it, so it no longer fails
  where the home directory is a symlink, as on the cluster (#124).
- `grnbench.seeds` (#126, step 1a): `DEFAULT_SEED = 42`, `derive` (named child seeds from `SeedSequence`,
  31 bits), `SeedPlan` (`legacy` / `master` / `derived`; `seed_for` = the seed a component really uses,
  `kwargs` = what its call passes, under the real parameter name; knobs per entry point follow the call site)
  and `LEGACY`, the seed each component gets today per entry point (prep, infer, eval and run CLIs, interface,
  real_bench with and without a reference, engine_bench, candidate_study). `tests/test_seeds_legacy.py` reads
  every value back from its default, constant or call-site literal, follows the call path that hands it on,
  and checks that every RNG constructor and known RNG consumer in the functions holding it uses it (model
  internals, simulator and task enumerations included); it does not see a random call through an object
  outside its consumer list, nor arithmetic on a seed inside a callee it does not trace. Not covered
  (`seeds.NOT_COVERED`, their registered values pinned): `beeline_backtest.py`, `motif_sensitivity.py`,
  `calibration.py`, `demo.py`, and Python and entry-point plugin runners. Nothing calls it yet: no seed,
  default or output changes.
- A user override keeps what the engine chose (#113): a typed value or a knowledge-base choice other than the
  selector's carries `auto_value`, `auto_reason` and `auto_source`, the executable value the engine's own plan
  would use (e.g. the MAD bound when `fixed` is chosen), its reason and source, written to `plan.json` only on
  those decisions; a plan built without choices or with the selector's own is byte-identical; setting the
  engine's value again withdraws the override; a broken trio or an unknown key in `plan.json` is refused; the
  interface reads "auto: X (why) → you: Y (chosen: option)" and 0 switches off a threshold the engine left off.
  "Apply ... and rebuild the plan" in step 3 now keeps the values the new choices set (e.g. `fixed` 200 / 20%):
  the threshold boxes, which kept their old values on the server and in the browser and wrote them back as the
  engine's, are given the new plan's values; a value typed in a box stays typed, also one sent with the click.
- The `cli` recorded in `contexts.json` carries `--seed` (a required argument of `PrepPlan.cli_args`, passed by
  `execute_plan`, which also writes `seed` on each record), so a screened context (prepped at seed 0) no longer
  replays at the CLI default 42 (audit 2026-09-27), and `--force` when `execute_plan(force=True)` ran past a failed
  QC gate (issue #117); both tested by replaying the command. The interface's "Reproduce from a terminal" shows one
  command per context with the seed *Run cleaning* uses (`DEFAULT_SEED`, now passed explicitly by the interface), and
  every command it shows is shell-quoted (`shlex.join`: labels such as `T|reg` or `NK(bright)` paste back exactly).
  `tests/test_seed_discipline.py` fails on any new random call without an explicit seed in `src/grnbench` or
  `benchmarks` (none today). No default, seed value or result row changes (issue #111).
- Tests: every registered model (17 default, 11 extra) runs once through the registry on tiny data
  (`tests/test_registry_smoke.py`). Each must give edges to every regulator of its list and to no other gene,
  with weights that are not all equal, sorted high to low; the cisTarget methods must report motif-supported
  edges, and every motif method must rank the planted motif targets well above co-expressed genes without a
  motif. Only aracne3 (java, `ARACNE3_JAR`), regdiffusion and regdiffusion_motif (regdiffusion, torch) may
  skip, as their dependencies are outside the lock file, and only while one of them is really absent
  ("not tested here: <registry reason> (waits for: <dependency>)"); a crashed availability probe, an unusable
  synthetic cisTarget folder, a method unavailable for any other reason and a registered method without a
  smoke run all fail (#114).
- EXP-20261003-LM-1 section 4b (#95): rerun on the Bocconi HPC directed by GB, commands typed by an agent session (signature pending). Both primaries
  are identical in every field (`benchmarks/compare_confirm.py`, new), but the registered row check
  `compare_runs.py --rows-only` says DIFFERENT: it compares the per-candidate timing fields inside
  `candidates`. `benchmarks/diff_rows_nested.py` (new, diagnosis) finds every non-timing value
  identical. By #95's rule C7 stays `confirmed`.
- EXP-20261002-LM-1 section 4b (#94): GB's second-person rerun on the Bocconi HPC. The three
  declared primaries MATCH at `--rtol 1e-9` and the report is identical apart from its timing columns; one row field differs
  (`n_warnings` 4 -> 3 on one dataset, not a score; most likely the dataset README, which the data
  check does not hash), so `compare_runs.py` says DIFFERENT overall. By #94's rule, and signed on GB's
  instruction (typed by the agent, record 4b), C6 is `reproduced`.
- EXP-20260926-LM-1 section 4b (#93): GB's second-person rerun on the Bocconi HPC matches the
  committed results exactly (every non-timing row field, and the report apart from its timing
  columns; both primaries MATCH at `--rtol 1e-9`). Signed on GB's instruction (typed by the agent,
  record 4b); C5 is `reproduced`.
- EXP-20261003-GB-1 run and recorded (issue #52): every dataset of C6 and C4 passed gate 1 and gate 2;
  verdicts from the generated report: C6 inconclusive, C4 mixed, combined inconclusive, so the caveat
  on C4/C6 stays as worded. Results in `benchmarks/results/EXP-20261003-GB-1/`, record section 3.
- EXP-20260925-LM-2 section 4b (#24): GB's second-person rerun on the Bocconi HPC is DIFFERENT from
  the committed rows (exact and `--rtol 1e-9`; H1 unchanged in direction, +0.183 vs +0.182). C4 stays
  `reproduced-computational`.
- `benchmarks/compare_runs.py`: judges a rerun of a registered experiment against its committed
  results from files (task ids, every row field with timings excluded, `REPORT_REAL.md` modulo line
  endings and timing columns, and the record's `primary comparison`); `verdict: MATCH` / `DIFFERENT` /
  `NOT EVALUABLE`, optional tolerance, `--json` for the record. A primary is required (`--rows-only`
  says so in the verdict); a record or `--primary` that does not parse is an error, and a primary with
  fewer than 3 usable units on a side is `NOT EVALUABLE`, never `MATCH`. For the LM-2 human
  reproduction (#24).
- Motif methods (`motif_prior`, `*_motif`, `*_cistarget`) record cisTarget coverage (genes found, regulators with a motif, database names) in the run provenance as `motif_coverage` and warn when under 50% of genes are found or no regulator has a motif (Ensembl IDs?); scores unchanged (#76 M13).
- `grnbench.explore` (EXPLORATORY, descriptive, no claim): per-(dataset, TF) atlas with detection rate and per-TF AUPRC ratio, TF detectability and motif-vs-coexpression figures, per-TF evidence cards, from saved run artifacts.
- `docs/related_work.md`: the field's research lines and narrative arcs, open questions (2026), where our
  results fit, and the figure and explainability conventions that earned trust; every claim labelled.
- `python -m grnbench.board --digest [--since DATE] [--post ISSUE]` prints a weekly digest (merged and open
  PRs, changed experiment records, claim status changes, new decisions, running/experiment issues);
  `--claims` prints the claims status table joined with experiments and generated files. Posts only with `--post`.
- `benchmarks/fetch_perturb.py`: the Replogle 2022 CRISPRi screens (K562 genome-wide, RPE1; also K562
  essential and Norman 2019) from scPerturb (Zenodo 7041849, CC BY 4.0), pinned by size and md5, into
  any data root. Fetched as parallel byte ranges: one connection from the Bocconi login node gets about
  0.25 MB/s. The join appends the parts in place after an md5 check of the stream (no second copy of
  the file, which filled a home quota at 3 GB of 8.8 on the cluster) and reports any I/O error by file name.
- Issue #39 (#40): the default `lasso` branch centers X and y before `lasso_path` (no intercept), so
  lasso scores change; its provenance records `params.path = centered`. EXP-20260917-PRE-1's lasso
  rows (and so claim C2's comparison family) came from the uncentered path. `--threads` (>= 1, default
  4) now bounds random_forest, GENIE3 and GRNBoost2 in primary inference (`grnbench infer`/`run`,
  `beeline_backtest.py`), not only in the stability re-runs: they used every core before, so wall
  time changes on machines with more than 4 cores; scores do not. A dask Client that cannot start now
  raises with its cause. The BEELINE backtest refuses duplicate dataset/method scores;
  `grnbench.sbatch` and `beeline_backtest.sbatch` take the checkout from `SLURM_SUBMIT_DIR` and refuse
  one without `pyproject.toml`.
- Follow-ups to #40 (#46): the default `community` scores change too, because `lasso` is one of its default
  members (they moved with the centered path). `extra_trees` now honours `--threads` (it used every core
  before; scores do not change). The default lasso path calls `lasso_path(alphas=20)` instead of the
  `n_alphas=20` that scikit-learn 1.9 deprecates and 1.11 removes; the scores are bit-identical.
  `grnbench.sbatch` and `beeline_backtest.sbatch` now export `PYTHONPATH="$REPO_DIR/src..."` like
  `engine_bench.sbatch`, so the job runs the checkout `GIT_COMMIT` names, not the env's install.
- `grnbench report --paper` builds again: the record's `report input` path is the first quoted span
  (EXP-20261002-LM-1 ends it with a full stop); a run with no per-model scores (a candidate study) gets a
  note instead of an empty table that stopped LaTeX; `compile_tex` reads engine output as UTF-8 and says
  when an installed engine wrote no PDF instead of "no TeX engine found".
- Audit of 2026-09-27 brought onto main (`audit/2026-09-27/`); dated caveat on C4-C6 and in
  `docs/models_that_work.md`: the motif-prior lift is partly TF-level and database-membership (#41).
- Eleven more inference methods from the literature, native and registered as *extra*
  (`grnbench methods`, `--methods all` or by name): `bicor`, `wgcna_tom`, `propr_rho`, `aracne_dpi`,
  `mrnet`, `tigress`, `bbsr`, `pcnet`, `extra_trees`, `gradient_boosting`, `community` (mean rank of
  fixed members). `MethodSpec.default` and `list_methods()["default"]` mark the set that
  `--methods available` (and `grnbench run`) runs; it is exactly the 17 methods registered before, so no
  past experiment changes meaning. They are reimplementations from the published descriptions, not the
  original packages: `docs/models_catalog.md` says where each differs and lists what is not built in and
  why (PIDC, time-series methods, multiome methods, deep models). Exploratory until an experiment record
  declares them; `clr` now shares `_clr_scores` with `bbsr`, same output.
- Interface (`grnbench ui`): **0 · Quick run** (data, models, Run: baseline cleaning, models, scores with a
  reference), a model picker with presets, family, speed and tier (default / extra / yours, the engine's own
  vocabulary) shared by the quick run and step 6, **8 · Explore data** (PCA or UMAP of the cells coloured by any annotation, before
  vs after cleaning side by side, gene and TF table, expression violins) and **9 · Explore results** (ranked
  edges with a right/wrong flag, network of the best edges, Jaccard agreement between models, precision-recall,
  per-TF recovery of reference targets). The explorers are marked exploratory (D-010); ties are handled with
  the repository's rule (tied blocks kept whole, expected hits at the cut), detection is read from the counts
  layer, TF recovery counts only reference TFs and targets that are genes of the data, one helper
  (`reference_for_scoring`) decides which reference every page scores on, and each model's edges are
  validated and its provenance written next to them. The explorers never make a whole sparse matrix dense
  (per-gene mean and variance are summed sparsely; only the sampled cells x selected genes are converted), the
  precision-recall curve counts the regulator x gene universe instead of listing it (equal to `evaluate_ranking`
  on the listed universe, tested), and the network is directed (A -> B and B -> A are two arrows; each arrow
  runs the whole edge, so plotly keeps its head wherever the two markers do not overlap, at any plot width).
  New logic in `grnbench/ui/logic.py` and `grnbench/ui/catalog.py`, tested in `tests/test_ui_explore.py`;
  nothing in the engine changes.
- Interface, found by driving the whole tool as a newcomer (issue #70, `tests/test_ui_regressions.py`): the
  step-1 path box shared its key with the session value of the same name, so Analyse from a file path always
  crashed; the cell-by-gene orientation of a table can be named when it cannot be guessed; Stop kills the
  model process (models run in a child process); a quick run that cleaned no context stops with a message
  instead of showing "Done"; a corrupt TF list or an unusable reference is one red line, not a traceback;
  with no TF list the picker marks lasso, tigress and community as slow (every gene becomes a regulator) and
  warns before running them; step 4 keeps the context scope fixed, warns when the plan was edited in step 3 and
  refuses a no-op or repeated alternative; Apply uses the choices on screen, not a hidden pending change; after
  a screening the interface names the winning pipeline and lists the cell types it left out (the engine ranks
  candidates by the mean over their own contexts; changing that rule is a decision for the owners, not this PR).
- Interface redesign (`grnbench ui`): ten screens in three menu groups (Start here / Step by step / Explore) with
  plain-language names (the numbers 0, 8 and 9 are gone; steps 1-7 keep theirs); a menu that ticks what produced
  something; on every screen a header saying where you are and what the screen is for; numbered cards; empty
  states that explain what is missing and offer the way forward; Back / Continue buttons; after a quick run the
  results come first (cells kept, genes, models run, the highest score with "one dataset is not a ranking") and the
  form folds away; a one-sentence quality verdict and a pipeline overview in plain language; the score chart marks
  the random level and the table explains each score in a tooltip; typed values survive a visit to another page.
  The look is one accent, neutral overlays that hold in light and dark, segmented radios and colour-blind-safe
  chart colours; `grnbench ui` passes it to Streamlit as command-line options (an option you give wins). Running
  "Clean" again after a quick run no longer fails with "Output directory must be empty": the next run is written
  beside the first (`clean-2`). Results say what they were scored on (full reference or held-out half) and how to read
  the AUPRC ratio; a quality warning says what it means; a plan warning about an option that cannot run is put in words
  and the option dropdowns show names, not ids. Light and dark follow the browser; the sidebar starts closed on a phone.
  The `ui` extra now asks for streamlit>=1.54 (the pages use `width="stretch"`, `icon_position` and the light/dark theme
  options; the interface tests were also run with 1.54 in front of the lock's 1.64). New `ui/steps.py`, `ui/style.py`,
  `ui/theme.py`; tested in `tests/test_ui_frame.py`.
- Tool end to end (#74), found by driving the command line and the interface on hard inputs: a matrix with
  NaN or infinity, a reference no candidate edge can reach, and contradictory `prep` options are refused
  instead of giving all-zero or all-NaN results; cell types whose labels differ only in punctuation, case or
  a '.' at the 80-character cut, or already spell another label's folder, no longer share a folder, and
  cells with no label are not a context; missing `--tf-list` and the prep warnings are said out loud;
  `grnbench kb --write-docs` to a folder says so instead of a traceback; unreadable files, `.loom` without
  loompy, repeated gene headers and stray whitespace in gene names are handled; `python -m grnbench` works,
  built-in commands answer `--help` and refuse stray arguments (`grnbench kb` too), `grnbench demo` writes
  only into a new or empty folder (not even next to a failed run's diagnostic), `grnbench ui` listens on
  localhost, `grnbench run` checks `--gold` first and skips a
  method that cannot run here; `grnbench viz` creates nothing on a wrong path and reads candidate-study
  rows. What is left is the checklist in #76; the screening ranking rule is #75.
- EXP-20261003-GB-1 registered (exploratory, issue #52): `benchmarks/motif_sensitivity.py`
  splits the motif-prior lift of Pearson (C6, C4) into target identity, TF-level reweighting,
  database membership and a TF-generic gene propensity from the kept artifacts, with calibration
  gates against the committed rows, a keyed permutation null shared by all datasets and built-in
  checks asserted on every context (stage 1 of the 2026-09-27 audit's proposal; GRNBoost2/RegDiffusion
  refused until its gap is fixed). Agent-review fixes on PR #55 (2026-10-06, script and tests): a planned
  unit without a row is reported as excluded, a fault inside `score_models` is a runtime error and a
  context without gold is NaN (not failed code checks), `run` refuses a missing cisTarget folder or
  reference and mixed commits, inputs or seed counts, the verdict line names rule 4, untested tests,
  rule 2 through T2's reverse alone and secondary pairs not compared, every gate-1 failure is listed in
  full, and tests pin the record's definitions, every code check, multi-context, mouse and degenerate
  contexts.
  Record section 2.1 (dated note, before launch) answers LM's review: provenance corrections, the
  simulations as rationale only, the cross-set order of verdicts, the motif-sharing TFs (descriptive
  `gene_tf_np` control and the shares of TFs with a motif partner, never decisive) and the Windows
  rerun procedure; the report prints the n of L, D, L_art and L_gene. Section 2.2 (dated 2026-10-10)
  adds LM's note 7 (C6 and C4 combine only with T1-T3 run on at least 15 analysed datasets each) and
  records LM's answers on the Windows rerun and on keeping #86 and #88 out of the job's commit.
- Portable per-dataset `diagnostics.json` (#56): preprocessing, candidate screening, standalone
  inference and `grnbench run` emit one schema-versioned structure with bounded failures, applied
  choices and before/after dimensions; recursive validation excludes paths, identities, barcodes
  and expression values so reports and the UI can safely share it. Runtime diagnostic writes are
  non-fatal, strict JSON is cross-platform, and a lone failed diagnostic no longer blocks the
  CLI-recommended retry while every real pre-existing artifact remains protected. Effective and
  requested preprocessing strategies are distinct, subset runs and failures use their real context
  label, paths containing spaces remain private, and raw-droplet dimensions retain the pre-calling
  barcode count.
- Ingest (#42): integer input that does not look raw (near-identical cell totals, 0/1 only) gets a loud
  warning recorded in the manifest, not a refusal until calibrated (`benchmarks/count_plausibility.py`);
  refusals say what the data look like and list count-like layers as suggestions.
- Knee cell calling ignores barcodes under 100 UMIs (DropletUtils' `lower`), as its docstring always
  said it followed: with 10 the ambient-to-empty drop won on heart_10k_v3 (92k cells for 7.7k).
- `benchmarks/fetch_10x_raw.py` and raw10x discovery in `candidate_study.py`: 8 public 10x raw and filtered
  matrices (#50) for EXP-20261004-LM-2.
- `paper/make_figures.py` also writes the main figures (`fig_motif_lift`: C4-C6 forest plot, AUPRC and EPR;
  `fig_reference_free`: C7) as PDF and 300 dpi PNG from the same committed results; `--out DIR` redirects all output.
- `paper/make_figures.py` generates the motif-lift tables (C4-C6, AUPRC and EPR), the reference-free
  selection table (C7) and the cleaning-variant table (EXP-20261004-LM-1) from committed results.
- `candidate_study.py`: every candidate also scored against the baseline on the genes both share
  (`pair_eval_*`), so gene-set-changing candidates are compared on one edge universe.
- EXP-20261003-LM-1: without a reference the fixed baseline beats the engine's former choice (claim C7);
  D-021 makes the baseline the engine's choice when no reference is given.
- `engine_bench.py` / `real_bench.py`: a worker killed by the OOM killer now ends the run (retryable,
  resubmission resumes) instead of hanging until the Slurm time limit (`pool_results`).
- Reproduction rule aligned in paper/README.md, experiments/README.md and AGENTS.md: a second person
  who directs and checks the rerun earns `reproduced`, even if an agent types the commands (#24).
- Engine robustness (found by EXP-20261002-LM-2): screening fails one candidate, never all, when a
  context cannot be screened; strategies may not emit NaN/inf (`PrepError`); Pearson residuals give a
  zero-count gene 0 instead of NaN; ComBat drops one-cell batches with a warning.
- D-020 pins reference networks per registered experiment: historical reproductions keep their
  hashes; new runs prefer versioned upstream releases or use a dated, exact snapshot manifest.
  `data_manifest.py` argument parsing is now order-independent for reliable strict checks.
- EXP-20261001-LD-1: completed 740/740 synthetic and 40/40 real task accounting, generated stress
  curves, real-data visualizations and LaTeX, retained every refusal/failure, and added portable
  task/data manifests. Real auto cleaning usually ties the baseline; the run remains exploratory.
- `real_bench.py report --out <RUN>` falls back to `tasks.json` when `real_tasks.json`
  is absent, so runs planned before D-019 (EXP-20260926-LM-1) regenerate unchanged.
- EXP-20261002-LM-1: the motif prior lift (C4, C5) replicates on 20 never-used datasets (claim C6).
- `benchmarks/candidate_study.py`: every engine candidate, the baseline and one-change variants scored
  on real datasets under controlled degradations (EXP-20261002-LM-2).
- Registered Lorenzo (`@lorenzo-de-santis`, initials `LD`) in the team metadata.
- Canonical repo is now the private GitHub repo (D-016). Rules for people and agents in
  `AGENTS.md`; onboarding prompt in `docs/START_HERE_AGENT_PROMPT.md`.
- Experiment registry `experiments/` with pre-registration (D-017), data provenance manifest
  `experiments/datasets.tsv` (`benchmarks/data_manifest.py`), `paper/` (claims ledger, reporting
  checklist, CRediT log, availability drafts, `paper/make_figures.py`), `CITATION.cff`.
- Pinned environment `requirements/lock-py311.txt`; `hpc/setup_env.sh` installs it. CI runs on
  Ubuntu / Python 3.11 with that lock, on PRs and `main` only.
- HPC: `hpc/hpc.sh` (check/sync/run/submit/status/log/cancel/fetch; per-user `.env.hpc`, never a
  username in the repo) replaces `hpc/sync_to_cluster.sh`; `docs/HPC.md`. `grnbench.sbatch` and
  `beeline_backtest.sbatch` now ask for 8 CPUs (the `stud` QOS cap), not 16.
- `.gitattributes` keeps shell scripts LF.
- Day-one team setup: `grnbench board` (runs in flight = issues labelled `running`, open PRs, STATE
  next steps); session start/end routines in AGENTS.md; short MAP.md and STATE.md.
- `grnbench viz` (`grnbench.viz`, one style file, Okabe-Ito palette): QC and embeddings before/after,
  decisions, score distributions, paired lifts, PR curves, top-k network, TF recovery and error
  analysis, motif lift by TF, calibration, stability, attribution helpers for learned models.
  `real_bench.py run --keep-artifacts` keeps edges, processed data and decisions for it.
- `grnbench report` (`grnbench.texreport`): LaTeX report per experiment and `--paper` draft with a
  bibliography from the KB; numbers only from result files (templates are digit-free, tested);
  `real_bench.py report` now also writes `viz/` and `latex/`. Structure follows the thesis
  pipeline's end-of-run report.
- Shared HPC skill `.claude/skills/grnbench-hpc/SKILL.md` (pointer in `.agents/skills/`), Antigravity
  rule moved to `.agents/rules/`; onboarding guide extended (from zero, how to work, look inside).
- Shared-work workflow: owner merge after an independent agent review when the other person is
  away (labelled as such; can be reopened); experiments are claimed with an `experiment` issue before
  the record, and the board lists claimed experiments; a recipe for stacked PRs after a squash;
  long local runs from a `git worktree`; claim status `reproduced-computational` (agent rerun from
  the record, human reproduction pending).

### Added
- `docs/runbooks/`: launch runbooks for the experiment issues #24 (LM-2 reproduction, executed
  2026-10-08) and #52 (EXP-20261003-GB-1, launched 2026-10-10, results in #92), with the lessons of the
  run; the three cluster traps are also in the `grnbench-hpc` skill (section 3).
- Context-specific synthetic validation (D-019): the simulator writes shared, union and per-cell-type
  truth; selection/evaluation splits happen before context expansion; pooled and split predictions
  are macro-scored against each context. `engine_bench.py --stress-grid` fixes a 740-task design
  across four network-sharing regimes, eight progressive corruptions and mixed worst cases, with
  generated statistical/relative/operational breaking points. Benchmark task lists and
  score-affecting run configurations are hashed; atomic rows have explicit retry/terminal statuses;
  final reports require complete accounting (`--partial` is watermarked). HPC jobs can retain the
  12 pre-declared synthetic artifact sentinels and now propagate `MAX_ATTEMPTS`.
- `data_manifest.py --check --require-all`: optional exact-snapshot validation reports registered
  paths missing from the data root; the default subset-compatible check is unchanged. HPC and
  onboarding docs now distinguish the dynamic setup cache from registered experiment roots.
- `grnbench ui`: Streamlit interface (data -> profile -> pipeline -> compare pipelines ->
  clean -> models incl. user models -> benchmark). `grnbench demo`, `grnbench kb`.
- `grnbench.kb`: preprocessing knowledge base - 14 stages, 15 decision points, 56 options,
  51 references, evidence level per option; `docs/knowledge/preprocessing_kb.md` generated
  from it (a test keeps them in sync).
- `grnbench.context`: dataset context from the matrix, AnnData metadata, README files and
  download pages (GEO text records), rule-based with snippets and conflict reporting.
- `grnbench.profile` / `grnbench.selector`: dataset profile (protocol from count values,
  unfiltered-droplet check, batch x cell-type Cramer's V, batch silhouette) and rule-based
  selection with reasons, uncertainty and one-change candidate pipelines.
- `grnbench.autoprep`: executable plans (MAD cell bounds, per-batch bounds, strategy
  resolution, one run per cellular context) with the exact CLI; `grnbench.screen`: candidate
  pipelines run and scored on the data, reference split into selection/evaluation halves.
- Strategies `metacells` (k-means on PCA within batch, counts summed, purity reported) and
  `cpm_log1p_center_batch`.
- `run_prep`/`grnbench prep`: `--subset-key/--subset-value` (applied after doublet calling on
  all cells), `--exclude-qc-genes`, `--per-batch-bounds`, `--max-genes-per-cell`,
  `--max-counts-per-cell`, `--max-pct-mito`, `--cells-per-metacell`.

### Results
- EXP-20261001-LD-1 (exploratory): 740/740 synthetic tasks completed and 40/40 real tasks reached a
  terminal state (33 completed, seven expected BEELINE refusals). On 23 paired scored real datasets,
  auto usually ties baseline and shows no consistent recovery gain. Two completed datasets with a
  usable reference (`paul15`, `pbmc3k_raw`) produced no model score and remain explicit findings.
- EXP-20260925-LM-2 (confirmatory, pre-registered): the motif prior lifts Pearson's AUPRC ratio by
  +0.182 [+0.138, +0.227] on 20 fresh CELLxGENE datasets (p_Holm 1.4e-7); EPR lift unresolved.
  Reproduced computationally (independent agent rerun from the record, identical numbers);
  human reproduction pending.
- EXP-20260926-LM-1 (confirmatory, pre-registered, first HPC experiment): on the same 20 datasets
  the motif prior lifts the AUPRC ratio for RegDiffusion by +0.184 [+0.132, +0.237] (n = 19,
  p_Holm 1.2e-6) and for GRNBoost2 by +0.130 [+0.108, +0.153] (n = 20, p_Holm 3.6e-10). EPR lift
  resolved for GRNBoost2, not for RegDiffusion. Reproduced computationally (independent agent
  rerun on the HPC, identical numbers); human reproduction pending.
- `paper/tables/T_models_real_EXP-20260925-LM-1.md` regenerated: it still had the 45-character
  truncated dataset ids fixed in #7 (numbers unchanged).

### Added (HPC experiment tooling, EXP-20260926-LM-1)
- `hpc/engine_bench.sbatch`: `MODELS` (models joined by `+`, since `sbatch --export` splits on
  commas) and `KEEP_ARTIFACTS=1` reach `real_bench.py run`; the job copies `GIT_COMMIT` to
  `$OUT/CODE_COMMIT`. `hpc/setup_env.sh`: `WITH_TORCH=1` installs RegDiffusion + CPU torch from the
  pinned line in the lock header. `hpc/fetch_cistarget.sh`: the six cisTarget v10 files.
- `fetch_cellxgene.py --from-manifest`: fetch exactly the datasets of an earlier manifest (no
  listing, no re-selection). `data_manifest.py --check`: verify a data root against
  `experiments/datasets.tsv` (exit 1 on any mismatch or unlisted file; LF/CRLF-only differences
  in text files are reported as `ok-eol`).
- `grnbench report`: several `;`-separated co-primary comparisons in a record each get a verdict.

### Added (VPN tooling for collaborators)
- `hpc/vpn.ps1` (Windows GlobalProtect): `status`, `reconnect` (restart PanGPS, wait up to 90 s,
  no elevation after a one-time grant), `setup` (prints the admin command with your own SID) and
  `grant <SID>` (idempotent `sc sdset` of one query/start/stop entry). `hpc/hpc.sh`: when the login
  node does not resolve or answer, the first error line is the VPN command for your OS. Docs:
  `docs/HPC.md` section 1b per OS (Linux `globalprotect` CLI, macOS none, both cited from Palo
  Alto docs); agents reconnect the VPN themselves and ask only for SSO, the grant, or on macOS.
  Windows ssh-agent: `StartupType Automatic` once, no per-reboot step.

### Fixed
- High-mitochondrial tissues are matched by whole words (#149 item 2; audit 2026-09-27): the fragments
  `renal`, `atri`, `ventric`, `septum`, `cardi`, `hepat` matched inside other words, so adrenal gland,
  "pediatric"/"psychiatric", brain ventricles, the lateral septum, the cardia of the stomach and the
  hepatic flexure of the colon relaxed the mito gate to 60% with a false reason. The tissue text is
  normalised first (`autoprep.tissue_words`: words split at a lower -> Upper boundary followed by a lower-case
  letter, so "iPSC" stays one word, and at letter/digit boundaries; lower case; `_` and `-` read as spaces; runs
  of white space as one). Only the matcher splits: the context stores and shows the label in lower case as on
  main, and the planner hands the matcher the label as written (`DatasetContext.spelling`). Plurals match (hearts, livers, endocardiums), so snake_case, hyphenated, camelCase, digit-suffixed and
  plural spellings (`skeletal_muscle_tissue`, `left-ventricle`, `SkeletalMuscle`, `heart1`, `hearts`) relax as
  on main. Of the 732 CELLxGENE tissue labels, 9 lose the relaxation (all of them non-cardiac, non-renal or
  non-hepatic) and none gains it. Over the labels and their spellings (6,981 strings: upper case, trailing
  space, hyphens, underscores, plural, camelCase, digit suffixes, lower-case concatenation), 103 lose it on
  purpose (those 9 in every spelling, plus adrenal, suprarenal gland, lateral septum, the lateral ventricle,
  pararenal fat, psychiatric, the hepatic flexure in any spacing or run together, liverwort, deliver, heartwood), bare
  "ventricle" loses it (heart or brain), the brown adipose tissue spellings main missed gain it (4), and no
  high-mito spelling regresses except one trade-off kept on purpose: a label written as one lower-case word
  (`skeletalmuscletissue`, `heartleftventricle`) is not split, because finding a word inside a run of
  letters is what relaxed liverwort, deliver and heartwood. Of the 53 multi-word high-mito labels written that
  way, main relaxed 52 and this rule 5 (the `hepat`/`nephr` words), so 47 lose it. The tissues of every
  committed row and registered dataset table match as before. The label as written that the matcher reads is the
  one with the strongest evidence (highest weight; the first seen on a tie), no longer the first seen: uns
  "tissue: skeletalmuscle" (weight 2) hid the curated `obs['tissue']` "SkeletalMuscle" (weight 3), read after it,
  and kept the 25% gate (review of #161).
- `aracne3` now drives the real ARACNe3 program (issue #105). The adapter ran
  `java -jar $ARACNE3_JAR -i … -o … -t …`, a command that neither ARACNe3 (a native C++ program that
  needs `-e`, `-r` and `-o`) nor its Java predecessor ARACNe-AP (where `-t` is the TF file) accepts, and
  it wrote cells as rows. It now runs `ARACNe3_app_release -e expr.tsv -r regulators.txt -o out
  --seed 42 --threads 1 --runid grnbench` in a fresh temporary working folder, found through
  `GRNBENCH_ARACNE3_BIN` or PATH (`ARACNE3_JAR` is no longer read), with genes as rows, and reads
  `consolidated-net_grnbench.tsv` strictly (weight = MI, direction = sign of the Spearman correlation);
  `extra_args` accepts only flags ARACNe3 reads. Written against upstream commit `3d8791a` and tested
  with a stand-in for the program (`tests/test_aracne3.py`) and then run against the real program on the
  HPC, where it works as written (EXP-20261018-GB-1). The method never ran before, so no result changes.
  A relative `GRNBENCH_ARACNE3_BIN` is made absolute when it is found: the program starts inside its own
  temporary folder, where the relative path passed the availability check but failed the run (LM review).
  `grnbench.seeds` now records aracne3's seed (library row `model.aracne3` = 42, the program's `--seed`,
  traced from the command line) and lists `benchmarks/aracne3_smoke.py` as not covered.
- ARACNe3 on the HPC (EXP-20261018-GB-1, issue #109): `hpc/aracne3_build.sbatch` builds the pinned
  upstream program inside a job (the cluster has no cmake; a pinned Kitware binary in `~/tools`) and
  refuses an incomplete submodule download; `hpc/aracne3_smoke.sbatch` runs `benchmarks/aracne3_smoke.py`,
  which calls the `aracne3` adapter against the real program on a fixed PBMC3k subset and records the exit
  status, messages, edges and whether same-seed runs repeat with 1 and 4 threads (`--resummarize` rewrites
  a finished run's summary from its files). Result: the adapter works as written (exit 0, ARACNe3 read
  one sample per cell, 478 edges parsed); same-seed runs repeat byte for byte with 1 thread; with 4
  threads the null model differed between runs while the network did not. `hpc/README.md` has the
  measured build steps, including the paced retry when GitHub refuses anonymous submodule clones.
- `aracne3` keeps a finished network when its temporary folder cannot be deleted: on the cluster's BeeGFS
  home one deletion failed with `EBUSY` after a successful ARACNe3 run, and the network was lost
  (EXP-20261019-GB-1, issue #125). The network is now returned, `attrs["temporary_folder_left"]` (the run's
  provenance) names the folder left behind and a `RuntimeWarning` says so; a run that failed still
  reports its own error. In a bench job neither is visible (`hpc/engine_bench.sbatch` sets
  `PYTHONWARNINGS=ignore` and `real_bench.py` keeps no provenance), so the folder stays silently in `TMPDIR`;
  `hpc/README.md` says when to put `TMPDIR` on the node's local disk instead.
- `aracne3` no longer takes the engine's generic thread count (independent review of #106): `grnbench infer`
  and `grnbench run` pass `--threads` (4 by default; `hpc/grnbench.sbatch` the job's CPU count) to every
  method, and with more than one thread ARACNe3's null model does not repeat from the seed
  (EXP-20261018-GB-1). ARACNe3 now gets `--threads 1` unless the aracne3-specific `aracne3_threads` asks for
  more; `method_info.threads` in the run's provenance is the count it got (`threads_not_used` names a generic
  count set aside). A failed run's error writes the temporary folder as `<tmp>` and the program's path as its
  file name (the error reaches `inference_manifest.json`), and a keyword the adapter does not take is refused
  instead of dropped. `benchmarks/aracne3_smoke.py` passes its `--threads` entries as `aracne3_threads`, so
  its 4-thread runs still get 4 threads; EXP-20261018-GB-1 has a dated note (its result ran from `6927b9f`).
- Candidate screening now carries the diagnosis' raw-droplet flag into its rebuilt profile. Raw
  10x inputs therefore keep knee cell calling instead of reclassifying called cells as already
  filtered and then applying QC to hundreds of thousands of empty barcodes (issue #30).
- Ingest aggregates gene identifiers that collide under grnbench's case-insensitive matching
  contract, retaining the first spelling and recording case-fold merges separately. This prevents
  one duplicate symbol from blocking every TF-resolved pipeline while preserving strict rejection
  for ambiguous matrices passed directly to inference (issue #30).
- Real-benchmark rows retain a portable structured summary of candidate/context failures and report
  scoreless rows with screening failures explicitly instead of hiding the cause behind `completed`.
- `grnbench report`: record keys containing commas or backticks were read as continuation text and
  fenced code was appended to the key above, garbling the Methods section; a note after the
  primary comparison containing `;` was read as a second hypothesis. LaTeX of EXP-20260925-LM-2
  and EXP-20260926-LM-1 regenerated from unchanged rows (numbers identical).
- `hpc/engine_bench.sbatch` imports grnbench from the checkout it runs (`PYTHONPATH`), so
  `CODE_COMMIT` always names the code that ran, not the conda env's editable install.
- RegDiffusion adapter: genes that are constant in the scored cells (possible after the 3,000-cell
  subsample) are dropped before training instead of making RegDiffusion refuse the whole context
  (`1 genes have 0 variance`; 1/20 datasets in EXP-20260926-LM-1, issue #13). They get no edges.
- `grnbench report --paper`: every experiment section defined the same macro names, so a
  second confirmed experiment would redefine the first one's numbers. Macros are now prefixed per
  experiment (test with two experiments).
- `real_bench._model_matrix` truncated dataset ids to 45 characters; two datasets with a common
  prefix collided and `compare_units` refused the table. Full ids now.
- Doublet calling silently fell back to a crude heuristic because scanpy's scrublet needs
  scikit-image, which was not installed (issue #2). scikit-image is now a dependency (pyproject +
  lock). A fallback now raises a `RuntimeWarning`, is recorded in the manifest step
  (`fallback`, `fallback_reason`) and in real_bench rows (`<arm>_doublets`) and reports.
- Per-batch scrublet refused batches of 1-30 cells and the whole dataset fell back to the
  heuristic (4/20 datasets in EXP-20260925-LM-2, issue #6). Batches under 100 cells
  (`hygiene.SCRUBLET_MIN_CELLS`) are now pooled and scored together; a pool still under 100 cells
  is kept unscored; missing batch labels are one level. The manifest step records
  `scrublet_groups`. On those four datasets every cell is now scored by scrublet.
- `test_list_methods_reports_unavailability` assumed java on PATH; aracne3 has two
  documented prerequisites.
- Texts that said something false, no behaviour change (issue #112; no rule, option id, score or
  default changed): the `heuristic` evidence level now says "not yet confirmed by a confirmatory
  experiment on the grnbench bench", like the module docstring; the `knee` cell-calling option says it
  keeps barcodes at or above the inflection after the knee (what `hygiene.knee_cell_threshold` does),
  not those above the knee, and when it stops with an error or finds the wrong cliff; that function's
  docstring describes the A-017 chord knee, what `None` means and how it differs from each DropletUtils
  knee method instead of claiming `barcodeRanks`' landmarks; ImageDoubler is cited as 2025 (Crossref),
  key `imagedoubler2024` kept; the interface docstring says the `grnbench prep` commands reproduce a
  cleaning and `plan.json` is a record of the choices that no command or screen reads back.
  `KB_VERSION` 2026.10.2, `docs/knowledge/preprocessing_kb.md` regenerated.

## 0.7.0 — the ranking is not real, and neither was the test that said so (2026-09-15)

### Retracted

The 0.6.0 "significance" layer treated **gold edges** as independent samples.
With ~4,500 gold edges per dataset that made 80% of pairwise method comparisons
look "resolved", and 47% of the resolved orderings appeared to contradict each
other across ground truths. Both numbers were pseudo-replication: the only
independent replicates in a BEELINE-style design are the **datasets** (n = 7),
and a paired test that knows this resolves almost nothing (see below). The
0.6.0 entry is kept for the record and marked retracted. Neither number should
be quoted.

### Removed

- `eval.significance.paired_bootstrap`, `significant_tiers`,
  `contradiction_report` and the `grnbench-eval --significance*` flags.
  Removed rather than deprecated: a deprecated function that still prints
  "significant" is a trap, not a courtesy.

### Added

- **`eval.comparison.compare_units`** and `grnbench compare`: paired Student t
  on declared independent units, Holm over both directions of every pair,
  Bonferroni simultaneous intervals. Refuses fewer than three units, never
  reads "unresolved" as "equivalent", and requires the caller to state what a
  unit is. This is the reference procedure, not a new estimator.
- `benchmarks/beeline_backtest.py` now ends with **`unit_comparison`** —
  `compare_units` with datasets as units, one family per ground truth × metric
  (`epr`, `auprc_ratio`) — and **`regime_check`**, the specific-chip scores
  against the published best (AUPRC-ratio 1.10 / EPR 1.20, factor 3,
  descriptive). Both land in `beeline_backtest.json` and `summary.md`;
  `--summarize existing.json` recomputes them on a finished run.
- `eval.universe`: an explicit, fixed candidate universe. Missing predictions
  get one bottom tie and stay in every denominator, so a method that emits
  fewer edges is not rewarded for it.
- `eval.significance.gold_thinning_sensitivity`: what remains of the old
  layer, renamed for what it is — annotation-removal *scenarios* with fixed
  `k`, reported as quantiles, not confidence intervals, and never a winner.
  `reference_ordering_differences` reports order reversals between references
  as differences: ChIP binding, STRING association and lof/gof response are
  different relations, and disagreeing about them is not a contradiction.
- `eval.metrics.metric_from_labels(metric, weights, y, k)`.

### Changed

- **`null_metrics` is ~30× faster** on real edge lists. A permutation of the
  scores over fixed pairs is the same draw as a permutation of the labels over
  the fixed, already-ranked score vector — same distribution, same tie blocks —
  so the `random` and `regulator` modes no longer re-validate 370k identifier
  strings on every draw (6.5 s → 0.2 s per draw on hESC). Pinned by a test that
  compares both routes on a heavily tied ranking. Draws for a given seed differ
  from 0.6.0; `degree` mode is unchanged because it does move the pairs.
- `tfs.resolve_tf_names` raises `TFListError` when a TF list matches nothing,
  instead of silently promoting every gene to regulator.
- `infer.classical` equal-frequency binning assigns a bin per distinct value
  (first rank), so tied expression values are no longer split by cell order.

### What the corrected inference says

9 methods (`pearson_signed`, `spearman_signed`, `mutual_info`, `clr`, `ppcor`,
`lasso`, `random_forest`, `grnboost2`, `genie3`) × 7 BEELINE datasets, pinned
BEELINE selection (significant TFs + top-500 significant genes), datasets as
units, α = 0.05, Holm:

| metric | resolved pairs (3 references × 36) |
|---|---|
| EPR | **0 / 108** |
| AUPRC-ratio | **0 / 108** |

GENIE3 and GRNBoost2 are not distinguishable from Pearson correlation on any
reference. Meanwhile the point-estimate winner changes with the reference on
7 of 7 datasets. Every ordering the point estimates suggest — including every
"winner flips between ground truths" — is within what seven paired
observations cannot distinguish. That is the finding: not that the methods
are equal, but that this design cannot rank them, and a leaderboard that does
is reporting noise. Archived in `benchmarks/results/beeline-2026-09-17-9methods/`.

The only resolved comparisons ever seen (4 / 45 on AUPRC-ratio, `pearson` /
`spearman` > `clr` / `mutual_info`) came from an earlier run that kept every
TF instead of the Bonferroni-significant ones, which BEELINE's own script does
not do. Fixing the selection removed them: the gene-selection rule alone moved
the only "significant" result the design had.

### Fixed

- `benchmarks/beeline_backtest.py` kept every TF present in the expression
  matrix; BEELINE's `generateExpInputs.py` keeps only TFs passing the
  Bonferroni p-value filter. hESC's candidate universe was 2.9 M pairs instead
  of 373 k. The selection is pinned to the official script's SHA-256 and tested
  against its output on a fixture (`test_subset_agrees_with_pinned_official_fixture`).
- `benchmarks/beeline_backtest.py` `regime_check` prints, for every EPR, the
  hit count behind it against the count expected at random. mDC has 756
  reachable ChIP edges in 263 k candidates: `clr`'s EPR of 4.55 there is 8 hits
  against 1.8 expected, and `pearson`'s 0.00 is 0 hits. Both are now labelled
  small-count instead of being read as a leak or a loss.
- **`run_prep` could write a 0-cell `processed.h5ad` and report success.**
  `--strategy pseudobulk --min-group-cells N` above every group's size (or
  `--force` past a QC gate on a gene-poor fixture) filtered every cell out
  *after* the gate check, and nothing looked at the strategy's output shape
  before writing it and exiting 0. Preprocessing now refuses an empty result
  the same way ingest refuses to guess an orientation: loudly, as a
  `PrepError`, before any output is written.
- **A strategy's `warnings.warn()` (e.g. pseudobulk dropping tiny groups) went
  only to stderr.** Nothing that only reads `manifest.json`/`report.md` — the
  documented run artifacts — could see it. Captured into `manifest["warnings"]`
  like every other pipeline warning.
- **`grnbench-prep --min-genes-per-cell 0` / `--min-cells-per-gene 0` silently
  reapplied the default (200 / 3).** `args.x or Thresholds.x` reads an explicit
  `0` as falsy. This CLI entry point had zero direct test coverage (every
  other preprocessing test calls `run_prep` itself, bypassing argument
  marshalling) — added `tests/test_cli_prep.py`.
- **`--group-key` and `--label-key` together audited only the group key.**
  `qc_group = group_key or label_key`: a 500/90/6 label split came back
  "balanced" because the four donors were. Now `class_imbalance` covers the
  labels and `group_imbalance` the group key, each naming its column.
- **Missing `--label-key` labels became a literal `"nan"` stratum** in
  stratified subsampling, with a proportional quota like any real class. They
  are now excluded from stratification, counted in the step and warned about;
  the cells themselves stay in the matrix (only `--group-key` drops rows, and
  the two help texts now say which does what).
- **Duplicated cell barcodes were reported and then carried through.** They
  are made unique with numeric suffixes; the original barcode is kept in
  `obs['barcode_as_loaded']` and the gate still warns (pooled libraries).

### Added — preprocessing

- **`--drop-doublets`**: scanpy's built-in scrublet (per batch when a batch
  column exists), scores kept in `obs['doublet_score']`, removal logged as a
  step. Falls back to the joint-extreme heuristic *with a manifest warning*
  when scrublet cannot run. Never on by default. **Refuses when the call rate
  exceeds `--max-doublet-rate` (0.20)**: on a homogeneous test set scrublet
  called 497 of 596 cells (83%) doublets, which is a threshold landing inside
  the singlet mode, not a doublet population — dropping it would have deleted
  the dataset and reported success.
- **`--strategy cpm_log1p_combat`**: ComBat on the detected/declared batch
  column after `cpm_log1p`. Never auto-recommended: it removes every
  batch-aligned covariance, including biology that was run one batch per
  condition. Verified on a synthetic two-batch set with gene-wise shifts: the
  mean per-gene batch gap drops below 20% of the uncorrected one. Raw counts
  survive in `layers['counts']`.
- `grnbench-prep` prints `[PREP REFUSED]` (no `--force` hint) for refusals
  that `--force` cannot fix — an empty result, a missing batch column, an
  implausible doublet call rate — and `[QC GATE FAILURE]` only for gates.

## 0.6.0 — how much of the ranking is real? (2026-09-11)

### Added

- **`eval.significance` — the layer a leaderboard cannot provide about itself.**
  *Retracted in 0.7.0: gold edges are not independent units; see above.*
  `grnbench-eval --significance` resamples the gold standard (subsampling without
  replacement, since the metrics are defined over a *set* of gold pairs and a
  with-replacement bootstrap would be silently collapsed), re-scores every method
  on the *same* draw so the comparison is paired, and reports:
  - `P(A > B)` and a 95% interval on the paired difference, per pair;
  - `fraction_resolved` — of all the pairwise comparisons a ranking implicitly
    claims, how many the data actually supports;
  - **tiers**: methods the gold standard cannot separate are grouped, instead of
    being printed in an order the reader will believe.

  Resampling the gold rather than the predictions is deliberate: confidence that
  A beats B is limited above all by which true edges happen to be in the
  ChIP-seq network, and those networks are partial, cell-line-specific and biased
  toward well-studied TFs.

  Cost: each method's gold-hit identity and tie-block structure are computed once
  and a draw is then a boolean mask over them, so 200 draws over a million-edge
  list take seconds rather than the tens of minutes a naive re-scoring needs.

- Six inference methods in `infer/classical.py`, numpy/scipy/scikit-learn only:
  `pearson_signed`, `mutual_info`, `clr`, `ppcor`, `lasso`, `random_forest`.
  The registry goes from 4 methods (of which 2 needed arboreto **and** a live
  dask cluster, and 1 needed a Java jar) to 10, of which 9 run in a plain
  `pip install grnbench[infer]` environment. A benchmark whose only tree-based
  method is arboreto-backed cannot run at all on a login node or in a
  reviewer's fresh venv.
  - `random_forest` computes what GENIE3 computes — per-target feature
    importance over the regulator set — with scikit-learn alone. On a 200-gene
    problem it took 24s against arboreto GENIE3's 171s, and the two agree on
    the top-20 edges (Jaccard > 0.4, pinned by a test that skips when arboreto
    is absent).
  - `clr` is Faith et al., *PLoS Biology* 5(1):e8 (2007); `ppcor` uses
    Ledoit-Wolf shrinkage because scRNA-seq is always in the p > n regime where
    R's `ppcor` inverts a singular covariance.
- `tests/test_classical.py` — 29 tests. Beyond the contract checks (canonical
  schema, descending sort, no self-loops, regulator whitelist honoured,
  determinism), each method is tested against **the thing it exists to do**:
  `mutual_info` must rank a U-shaped target that Pearson misses; `clr` must
  demote a broadcast hub that raw MI promotes; `ppcor` must break an indirect
  A->B->C edge that correlation keeps. A method that returns well-formed
  garbage passes a smoke test and fails the user.

### Changed

- The HTML report gained an **"Is the ranking real?"** section, and the README
  now positions grnbench against geneRNIB honestly: if you want a leaderboard,
  use geneRNIB.

### Audit of published harnesses

The tie-order defect fixed in 0.5.1 was checked against the field rather than
assumed to be general. Reading the actual sources:

- **BEELINE** (`BLEval/EarlyPrecision.py`) is **tie-safe**: it expands the
  selection to every edge tied at the cutoff and divides by the count actually
  selected. `BLEval/AUPRC.py` delegates to `precision_recall_curve` + `auc`.
  Its trapezoidal integration differs from average precision by at most ~8% in
  the GRN regime, and in the *conservative* direction — not a defect.
- **BEAR-GRN** (`R/reproduce_early_metrics.R`) does `arrange(desc(score))` then
  `head(n_edges)`. `arrange` is a stable sort, so where the cutoff falls inside
  a tied block the result depends on the input file's row order. Real, but
  measured at ~1.05x under uniformly distributed gold — a caveat, not a scandal.

The conclusion is the honest one: the published harnesses are largely fine, and
grnbench was the one with the severe version of the bug.

## 0.5.1 — what the BEELINE replication found (2026-09-07)

Running `benchmarks/beeline_backtest.py` against the real BEELINE corpus
(Zenodo 3378975, md5-verified) rather than self-authored synthetic data. Every
entry below is a defect the real data exposed and synthetic data had hidden.

### Fixed

- **Gene symbols were compared with exact string equality**, so an expression
  matrix and a ground truth that disagree on case shared *nothing*. BEELINE
  ships mESC expression in MGI case (`Pou5f1`) and every mouse network in upper
  case (`POU5F1`): 64,557 genuinely reachable gold edges scored as **0 hits**,
  every metric returned `NaN`, and not one line was printed. A wrong number
  wearing the costume of an honest null result — the exact failure the package's
  "unavailability is data, not a crash" rule exists to prevent. Six of the seven
  datasets were unaffected only because BEELINE happened to ship them already
  upper-cased, which is why no synthetic test had ever caught it.
  - `eval.metrics.label_edges` / `evaluate_ranking` / `eval.gold.coverage_report`
    now fold case by default; `grnbench-eval --case-sensitive` opts out.
  - New `eval.metrics.identifier_overlap()` classifies every edge-list/gold pair
    as `ok`, `case-mismatch` or `no-overlap`, and `grnbench-eval` prints
    `[GOLD CASE]` / `[GOLD MISMATCH]` with example symbols from **both** sides.
    A gold standard in a different namespace (Ensembl IDs vs symbols) or for the
    wrong organism now says so instead of silently reporting `NaN`.
- `benchmarks/beeline_backtest.py` applied the published cell-type-specific
  ChIP-seq ceiling (AUPRC-ratio 1.10 / EPR 1.20) to **every** ground truth, and
  so flagged six correct mHSC results as "above the published regime". The
  ceiling is now keyed by ground truth, only `specific-chip` is a falsification
  target, and the generic networks are reported with the reason their numbers
  legitimately run higher — including that **STRING contains co-expression
  channels**, which makes scoring a correlation method against it partly
  circular.

### Fixed — reporting defects found auditing the fix above

Chasing the identifier bug through the report surfaced four more ways a number
could come out wrong and *look* fine. Each was reproduced by running the code
before it was touched, and each now has a regression test.

- **A method that matched no gold edge at all was reported as statistically
  significant.** `null_metrics` computed `hits = sum(null_scores >= observed)`,
  and `>= NaN` is elementwise `False` in NumPy, so an undefined observation
  collected zero hits and got the *smallest* p-value the estimator can emit —
  `1/(n+1) = 0.022` — which `eval_report.html` then painted as a green **pass**
  badge beside a row of empty metric cells. The report simultaneously said "this
  method scored nothing" and "this method beat the null at p = 0.022".
  `null_metrics` now reports `observed_undefined: true` and a `NaN` p-value, and
  the HTML renders `n/a` in red.
- **`eval_report.html` never rendered the identifier diagnostic.** The verdict,
  the `--case-sensitive` state and the method-name collisions all existed in
  `eval_report.json` and the Markdown, but the HTML — the file `grnbench run`
  advertises and the one people actually forward — dropped them. A vacuous run
  now opens with a red banner *above* the verdict table, plus a
  **Gene-identifier matching** section; a renamed method says so in its own note.
- **`top_k_overlap` compared edges case-sensitively**, so two byte-identical
  networks differing only in `Pou5f1` vs `POU5F1` scored Jaccard **0.000** under
  a heading asserting that methods sharing no edges are not confirming each
  other. It now folds case like gold matching, with `case_sensitive=True` to opt
  out.
- **The HTML provenance block could not display the preprocessing environment.**
  `report["environment"]` is always set at scoring time, so the `or` fallback to
  the prep manifest was unreachable: a report could state `arboreto: not
  installed` directly beside the `genie3` edge list arboreto had produced on a
  cluster weeks earlier. Both environments now render, each labelled.

### Fixed — ranking metrics were measuring their own sort order

- **AUPRC, precision@k and EPR depended on the arbitrary order of edges that
  scored *identically*.** Both walked the stored row order, so for any method
  emitting coarse or repeated scores (thresholded output, integer scores, a
  large block of zeros) the headline number was partly noise. Six orderings of
  the **same** scores gave AUPRC 0.319–0.478 and EPR 3.0–5.0, while the
  tie-safe AUROC stayed at 0.9167 throughout — proof the ranking was identical
  and only the tie-breaking moved.
  - `average_precision` now reads precision/recall at distinct-score boundaries
    (the standard tie-safe definition; identical to the old result when scores
    are unique), and `pr_curve` draws the same points instead of the sawtooth
    the within-block order happened to trace.
  - `precision_at_k` returns the expectation under random tie-breaking when the
    cut at `k` lands inside a tied block.
  - `evaluate_ranking` additionally reports `tie_block_at_k` and `epr_range`,
    the worst and best EPR reachable by reordering tied edges alone. A wide
    interval means the headline EPR is a coin flip, and now says so.
  - Re-running the BEELINE backtest, this moves AUPRC-ratio by at most 0.027 and
    EPR not at all for `spearman_signed` and `ridge`, which emit continuous
    scores — the defect was real but those two methods never triggered it.

### Fixed — ingest and preprocessing

- **A genes×cells table was silently loaded transposed.** Orientation was
  decided by gene-likeness alone, but a 10x barcode (`TGGCCAAAATGTGGTG-1`)
  matches the gene-symbol pattern perfectly, so on the standard GEO/BEELINE
  layout both axes scored ~1.0, the 0.1 margin never cleared, and a
  4000×600 matrix was ingested as 4000 cells × 600 genes — recording
  `transposed: false` as though that were a determination. Every QC metric was
  then computed on the wrong axis: median genes/cell 331 instead of 2204, and
  `n_mito_genes` **0 instead of 2**, which silently disables the mitochondrial
  gate. Detection now keys on barcode shape and on the presence of mito/ribo
  markers, and **refuses to guess** when neither decides;
  `--orientation genes_x_cells|cells_x_genes` states the layout outright.
- **Stratified subsampling deleted whole strata.** Every zero-quota class was
  bumped to one cell, overshooting the budget, and the `picked[:n_keep]` slice
  that repaired it cut in alphabetical class order. With 200 singleton classes
  and one class holding half the cells, the majority class kept **zero** cells
  while the manifest still read `stratified by 'ct'`. Allocation is now purely
  largest-remainder (that class now gets its proportional 10 of 20), and any
  stratum that cannot be represented is counted in the manifest and raised as a
  run warning.
- **`filter_cells` removal counts did not add up.** Zero-count cells are
  necessarily also low-gene cells, but both were reported as separate figures,
  so "removed 25 cells with < 200 genes, 10 cells with zero counts" described 25
  actual removals. The step now reports the union plus `n_cells_removed`.

### Added

- `eval.metrics.GoldIndex` — the gold standard normalised once and reused. Null
  models re-score the same set `n_nulls` times per method; on a ~10^6-pair
  ChIP-seq network, re-normalising per draw was the entire runtime.
- The backtest records an `id_match` verdict per row and a `ceiling_verdict`
  (`in-published-regime` / `above-published-regime` / `no-published-ceiling`).

## 0.5.0 — plugin release (2026-09-05)

Closes the v0.4 roadmap: custom-model plugin contract, HTML comparison reports,
stability integrated into evaluation.

### Added

- **Plugin contract** (`grnbench.plugins`). Bring your own GRN method three ways,
  all landing in the same registry, edge schema, manifests and evaluation as the
  built-ins:
  - a `grnbench.methods` **entry point** (auto-discovered, no flags);
  - `--plugin path/to/model.py:MyModel` for a one-off script (no packaging);
  - `--plugin spec.json` describing a **Docker/CLI** command that reads
    `expr.parquet` + `tfs.csv` + `params.json` and writes `edges.parquet`
    (`"format": "tsv"` for pyarrow-less environments).
  A plugin cannot silently shadow a built-in name; a plugin whose import,
  binary or image is missing stays listed with its reason instead of vanishing.
- **`grnbench run`** — preprocess + infer + evaluate in one command, with a
  `run_manifest.json` recording every stage, its exit code and why a stage was
  skipped. Also a `grnbench` umbrella command (`run|prep|infer|eval|methods|
  strategies|version`).
- **HTML comparison report** (`grnbench.html_report`) — self-contained, no CDN,
  no JavaScript, no plotting dependency: PR curves and bar charts are inline SVG,
  so it opens on a compute node with no network. Sections: verdict table, PR
  curves, nulls (degenerate draws surfaced, not hidden), cross-method top-k
  Jaccard heatmap, stability, gold reachability, preprocessing QC gates,
  provenance.
- **Stability in `grnbench-eval`** — `--stability-data processed.h5ad` re-runs
  each scored method on cell subsamples; `--stability-{n,frac,top-k,methods,tf-list}`.
- **Cross-method overlap** — `eval.metrics.top_k_overlap`, reported in Markdown
  and as an HTML heatmap.
- `eval.metrics.pr_curve`; `infer.prepare_adata`; `infer.register_method`;
  `MethodSpec.probe` (availability re-evaluated per listing, so a jar appearing
  on disk is picked up without a restart), `MethodSpec.takes`, `MethodSpec.origin`.
- `examples/plugins/` — `ridge_model.py`, a real scikit-learn method usable as a
  plugin, plus a Docker spec template.
- `hpc/` — SLURM job scripts and notes.
- `benchmarks/beeline_backtest.py` — replicates the BEELINE protocol (Pratapa
  et al., Nature Methods 2020) on BEELINE's own data and ground truths, and
  checks the numbers against the published regime. Because the literature
  reports a *ceiling* of AUPRC-ratio 1.10 / EPR 1.20 on real data, this is a
  falsification test: a result far above that means the harness is leaking, not
  that the method is good. Ships with `hpc/fetch_beeline.sh` and
  `hpc/beeline_backtest.sbatch`.
- 48 new tests (115 total).

### Fixed

- **10x `.h5` loader silently split duplicated genes.** The `.h5` path called
  `var_names_make_unique()` before `_clean_gene_names`, renaming CellRanger's own
  duplicated symbols (e.g. `TBCE`) to `TBCE-1` and leaving the aggregation branch
  dead. The same biological data therefore produced a different matrix than the
  MTX path, the manifest recorded `n_duplicates_merged: 0`, and the renamed copy
  was invisible to TF matching — so a duplicated TF could be dropped by HVG
  selection, violating the "HVG never drops TFs" rule. Counts are now aggregated
  by sum on every loader path. Regression-tested against a real CellRanger-v3 file.
- **`pip install '.[infer]'` was unsatisfiable**: the pin was `arboreto>=1.1`, but
  arboreto's latest PyPI release is 0.1.6. Now `>=0.1.6`, with `dask`/`distributed`
  pinned alongside.
- `environment_block()` no longer reads the deprecated `scanpy.__version__` /
  `anndata.__version__` (FutureWarning on every run); it uses
  `importlib.metadata` and also records the arboreto version.
- `register_plugins` is idempotent for the same `--plugin` argument, so
  `grnbench run` passing it to both the inference and the evaluation stage is no
  longer mistaken for a name collision.
- Null models no longer compute a PR curve per draw (`evaluate_ranking(...,
  with_curve=False)`).
- `run_method` rejects a method returning a non-DataFrame or fewer than three
  columns with a clear contract error instead of an opaque `IndexError`.
- ARACNe3 availability is now a `probe` on its spec rather than a hardcoded
  special case in `list_methods()`.

Found by replicating the BEELINE benchmark on real published data
(`benchmarks/beeline_backtest.py`):

- **Every BEELINE ground-truth network was unreadable.** BEELINE ships its
  ChIP-seq / STRING / lof-gof networks as `Gene1,Gene2`, and that convention
  propagated to most published GRN ground truths. The alias table in
  `eval/gold.py` accepted `gene2` as a target but not `gene1` as a source, so
  all ten shipped networks raised "need source/regulator + target columns" —
  grnbench could not read the reference benchmark of its own field. `gene1`,
  `g1`, `g2` and the symbol variants are now accepted, and the error message
  lists what *is* accepted. Regression-tested in the real format.
- `load_gold` read `net.csv.gz` as TSV, because the separator was chosen from
  the `.gz` suffix rather than the name underneath it, silently yielding one
  wide column. The suffix is now resolved with `.gz` stripped.

Found by the pre-release audit:

- **`grnbench-eval` silently dropped an edge list when two inputs derived the
  same method name.** The name came from `path.name.split(".")[0]`, so scoring
  the same method across two runs (`runs/a/genie3.edges.tsv.gz` and
  `runs/b/genie3.edges.tsv.gz`) fed two files in and reported one method, with
  no warning — a provenance lie. Colliding names are now qualified with the
  parent directory (`runB/genie3`), announced on stderr, and recorded in the
  report under `name_collisions`.
- `canonical_edge_frame` resolved columns by iteration order, so a plugin
  emitting BOTH `source` and `regulator` (or `weight` and `score`) could have
  its regulator column chosen by luck. An exact canonical name now always wins
  over an alias.
- `hpc/grnbench.sbatch` split the `PLUGINS` variable on `:` — but every plugin
  reference grnbench accepts contains a colon (`pkg.mod:Class`,
  `file.py:Class`), so each reference was torn in half. The separator is now
  `;`, and the documented example in `hpc/README.md` was itself broken.
- `hpc/grnbench.sbatch` hardcoded `--stability`, an N-fold runtime multiplier
  with no way to switch it off on exactly the machine where runs cost hours.
  Now `STABILITY=0` / `STABILITY_N=...`.
- `pyyaml` was imported for YAML plugin specs but declared in no extra: a
  documented capability with no supported way to install it. Added to
  `[plugins]`, `[dev]` and `[all]`, and covered by a test.
- **A misspelled method name was reported as an unavailable method.**
  `--methods geni3` printed `[SKIP] geni3: unknown` — the same channel that says
  "arboreto not installed" — so a typo looked like a missing dependency. Unknown
  names now fail the run with `[UNKNOWN METHOD]`, a close-match suggestion
  (`Did you mean: genie3?`) and the list of registered methods. A known method
  that simply cannot run here still skips loudly, and a run where *nothing* can
  run now exits 2 instead of reporting success over zero methods.
- **`grnbench-infer` crashed writing its manifest when no method ran.** The
  output directory was only created as a side effect of `save_edges`, so a run
  with zero successful methods died with `FileNotFoundError:
  inference_manifest.json` instead of failing cleanly. The directory is now
  created up front.

### Changed

- `list_methods()` reports an `origin` (`builtin` / `entry_point` / `python` /
  `command`) and triggers entry-point discovery.
- `grnbench-eval` writes `eval_report.html` by default (`--no-html` to skip) and
  accepts `--prep-manifest` to embed the preprocessing QC gates in it.
- `inference_manifest.json` records registered plugins and the reason each
  skipped method was skipped.

## 0.4.0 — data-hygiene release

Preprocessing (`grnbench-prep`), inference adapters (`grnbench-infer`),
evaluation (`grnbench-eval`). Counts-verified ingest, loud QC gates, strategy
registry with auto-recommendation, HVG ∪ TF feature selection, canonical edge
schema, AUPRC/EPR/AUROC with degree-preserving nulls, full provenance manifests.
