# Sensitivity analysis design: how much of the C4/C5 motif-prior lift is target-level sequence evidence?

> Implemented since (2026-10-09): see `experiments/EXP-20261003-GB-1.md` and `benchmarks/motif_sensitivity.py` (PR #55). This proposal is kept unchanged below as history.

## 0. Summary

- **Only partly doable from the kept artifacts.** `--keep-artifacts` keeps only the top 20,000 edges per model, arm and context (`benchmarks/real_bench.py:131-137`). The rank product (`cistarget.py:316-325`) and AUPRC with `universe=` (the rest tied bottom, `eval/universe.py:72-82`) both use the whole TF×gene list, so the registered numbers can't be rebuilt from those files.
  - **C4 (Pearson):** needs no refit. Pearson is closed form, recomputed exactly from the kept `processed.h5ad`.
  - **C5 (RegDiffusion, GRNBoost2):** needs a deterministic refit of the base model only (seeds 0 and 42, as registered). A gate checks that the refit's top 20,000 equals the kept file.
- **The analysis itself is cheap** once the full base edges exist: about 0.25 s per control variant on a 172,350-pair universe, measured locally. One 8-CPU `stud` job of about 1.5–2.5 h covers both experiments. The refits take most of that time.
- **The script was tested on synthetic run folders** made by `real_bench.score_models` itself (so they
  have the exact kept-artifact layout), and separated a planted-motif truth from a no-motif truth as
  intended. **KNOWN GAP, found re-running that validation for this report, not yet fixed:** the
  Pearson-only calibration gate reproduces cleanly (`|diff| = 0.0`), but the GRNBoost2-based gate does
  not — the recombined `grnboost2_motif` pooled AUPRC ratio differed from the mock's own row by 0.19-0.31
  on a fresh rerun. The script's internal self-consistency check (its local rank-product reimplementation
  vs. `cistarget.motif_rank_product`, given the *same* base edges) held, so the reimplementation itself is
  not the problem; the likely cause is that `with_motif_prior` fits GRNBoost2 a second, independent time
  when building the `*_motif` combination, and that fit did not reproduce the standalone `grnboost2` fit
  bit-for-bit under this environment's thread/dask defaults (both nominally seeded 42). This needs to be
  understood and fixed — and re-validated, ideally against the real, already-reproduced C5 artifacts
  (`.projectos/STATE.md` records an exact independent HPC reproduction of the real job, so the underlying
  pipeline is not in doubt) — **before running this on the cluster.** It has not been run on real data or
  the cluster at all (unreachable from here).

---

## 1. Can it run from the kept artifacts only? Which files?

**Kept files per dataset** (code: `real_bench.py:183-203`, `293-309`). Layout: `OUT/artifacts/<task>/` with:

| file | writer | used for |
|---|---|---|
| `reference_selection.tsv`, `reference_evaluation.tsv` (source, target) | `run_real_task` 296-299 | gold = evaluation half. The union of the two halves = all reference pairs (`split_reference` partitions per TF, `screen.py:36-47`), so the regulator list `tfs = sorted(sources of the union)` is recoverable without the raw data. |
| `decisions.json` | 300-309 | `chosen_plans[arm].choices` → `autoprep.resolve_strategy(choices)[0]` gives the strategy. `assume_normalized = strategy != "raw_counts"`, exactly as `score_models` line 168. |
| `<arm>/ctx<i>/processed.h5ad` (arm ∈ {baseline, auto}) | 186 | The AnnData actually scored, written after the ≤3,000-cell subsample (174-176), so no second subsample. It gives the genes, hence the universe (179-180), the `expr` for the motif matrix, and the base-model input. |
| `<arm>/ctx<i>/context.json` | 187-189 | context, n_obs, n_vars |
| `<arm>/ctx<i>/<model>.tsv.gz` | 202-203 via `_write_edges_artifact` | **Top 20,000 by weight only** (`ARTIFACT_TOP_EDGES`), columns source, target, weight. |

- `ctx<i>` is numbered by record index, so it can skip numbers. Contexts are averaged in that order.
- `--keep-artifacts` switches off the baseline==auto dedup (line 164), so both arms have their own files.
- The models kept depend on the run:
  - LM-1: `regdiffusion`, `regdiffusion_motif`, `grnboost2`, `grnboost2_motif`.
  - LM-2: `pearson_signed`, `pearson_motif`, `motif_prior`.
- For `cxg__10x-3-v2__homo-sapiens__entorhinal-cortex__06b91002`, the `regdiffusion*` files cannot exist. The model raised before `_write_edges_artifact` was called, and the committed row has NaN plus `m_regdiffusion_*_error` for both arms (checked). That dataset stays excluded, so n = 19, as registered.

Also needed:
- The run's `rows/<task>.json` (registered values).
- The cisTarget files: 4 `.feather` and 2 `.tbl`, sha256 in `experiments/datasets.tsv`.
- `references/collectri_{human,mouse}.tsv`, used for "labelable".

**Where the files are:**
- **C5:** `~/gb_runs/EXP-20260926-LM-1/` (rows, artifacts ≈ 190 MB per record section 3). An identical rerun is at `~/gb_runs/EXP-20260926-LM-1-repro/`. Data is in `~/grnbench_data_lm1/{cistarget,references}`.
- **C4:** the LM-2 artifacts are **not on the cluster**. They are "outside git in the run folder on Luca's laptop" (LM-2 record section 3), so they must be copied to the cluster first (commands in section 4).
  - LM-1's `processed.h5ad` cannot stand in for them. The 4 scrublet datasets changed with #11.
  - The committed rows also show a different auto candidate for `10x-5-v1__homo-sapiens__lung__7b3368a5`: `context_scope=pooled` in LM-2, `baseline` in LM-1.

**Truncation and what it implies:**
- The universe is |sources| × (n_vars − 1). With `MAX_TFS = 150` and `hvg = tf_plus_1000`, that is roughly 150 × 1,100 ≈ 165k pairs, so the kept 20,000 is about 12%. This is an estimate; confirm it from `context.json` n_vars.
- Pearson and RegDiffusion output the full regulator×gene block (`classical.py:124-131`, `regdiffusion_adapter.py`), so their kept files are always truncated.
- GRNBoost2 returns only links with importance > 0 (arboreto `core.py:209`, `links_df[links_df.importance > 0]`). A kept `grnboost2.tsv.gz` with fewer than 20,000 rows is therefore complete and needs no refit.
- Scoring the kept files directly measures a different quantity. On a mock with 26% of the universe kept, GRNBoost2 scored 2.6169 vs 2.8084 registered, and 3.2932 vs 3.5139 for the motif version. The script flags this as `kept_TRUNCATED` and fails the gate.
- Kept files must be read with `float_precision="round_trip"`. The default pandas parser was off by up to 9.7e-17, which broke the exact equality checks; with round_trip the difference was 0.0.

**So there are two tiers:**
- **Tier B (the one that answers the question).** Pearson is recomputed. RegDiffusion and GRNBoost2 are refit with `run_method(base, adata, tf_list=tfs, assume_normalized=assume)`, as `score_models` does.
  - The run is deterministic on this cluster: the LM-1 rerun matched all 1,260 numeric fields with max diff 0.
  - At HEAD, `src/` differs from the run commit 84f28e1 only in `regdiffusion_adapter.py` (zero-variance filter, 1bdace2; `git diff --stat 1d07221 HEAD` shows only that file among the scoring code) and `texreport.py`. The filter does nothing on the 19 datasets that ran, because they had no constant gene.
  - Gates: the refit's top 20,000 must equal the kept base file exactly; the recomputed combined top 20,000 must equal the kept `X_motif` file; and the dataset mean must equal the committed row to 1e-12.
- **Tier A (`--no-refit`).** Kept files only. It is descriptive, fails the gate by construction, and is not used for any verdict.

C4 caveat: LM-2 ran on Windows. Pearson on the cluster may differ from it in the last bit. Pre-declared fallback: if any LM-2 unit fails the gate, only those units are rerun on the machine that ran LM-2, with the same commit, and this is recorded as a deviation.

---

## 2. Metrics, controls, tests

These are computed per dataset × arm (auto primary, baseline secondary) × context × declared pair (`pearson_motif`/`pearson_signed`, `regdiffusion_motif`/`regdiffusion`, `grnboost2_motif`/`grnboost2`), then averaged over the finite contexts with `np.mean`, as `score_models` does.

**Motif matrix.** M = `cistarget.motif_ranks(tfs_in_E, expr.columns)`, the registered matrix. From the same files the script derives:
- `in_db[g]`: the gene matches a rankings column in at least one database.
- `K[t]`: the number of the TF's annotated motifs present, summed over the two databases.
- W: the worst value = min over databases of `n_total`. It is the tie value given to motif-less TFs and to genes missing from the databases, including unmapped Ensembl IDs (`ingest.py:311-321` keeps the ID when no symbol exists). Ensembl-ID genes can never be gold positives.

Each variant M′ goes through `rank_product(E, M′)`. This is a line-for-line copy of `motif_rank_product`, and the script asserts it gives identical order and weights for M′ = M.

**Metrics**, each for base, real and every control:
- (i) `pooled`: the registered metric, `evaluate_ranking(edges, GoldIndex(eval), with_curve=False, universe=U)["auprc_ratio"]`, plus EPR. The script's fast version is asserted equal to this within 1e-12.
- (ii) `macro`: the mean over TFs of the within-TF AUPRC ratio. Each TF's row keeps the ranking, with missing edges tied bottom. Only TFs with ≥ 5 evaluation positives count (`MIN_POS_TF`).
- (iii) `pooled_T` / `macro_T`: the universe restricted to U′ = targets in the database and labelable, where labelable means the target appears as a target of any TF in `collectri_<org>.tsv`. `pooled_R` / `macro_R` restrict further to U″ = U′ × TFs with K ≥ 1.

**Controls (M′):**
- (iv) Permutation controls, seeds 0–99:
  - `perm_all`: each TF's row permuted across all genes. Keeps the TF-level row distribution; destroys both target identity and the database mask.
  - `perm_db`: each row permuted only among in-database genes. Keeps the TF-level distribution and the database mask exactly; destroys only target identity. **This is the primary null.**
- (v) TF-level only (constant rows): `tf_mean` (the row's mean over database genes), `tf_kexp` ((W+1)/(K+1)), `tf_presence` (1 if K > 0, else W).
- (vi) Database membership:
  - `db_only`: 1 if the gene is in the database, else W, for every TF.
  - `worst_tie`: 1 if M < W, else W — exactly the tie set the real prior creates.
  - `tfmean_db`: the TF mean on database genes, W elsewhere. This is everything except target identity, with no noise added.

**Built-in checks (observed exactly on both mocks):**
- Every constant-row control gives a macro lift of +0.000 (zero variance, so "unbounded").
- On U″, the tie controls (`tf_presence`, `db_only`, `worst_tie`) give pooled_R and macro_R lifts of exactly 0.
- On U″, `tf_mean` and `tfmean_db` give a macro_R lift of exactly 0.
- So **macro_R is a metric on which neither the TF-level nor the database-membership mechanism can act.**

**Diagnostics per context:** number of motif-less TFs, number of genes not in the database, number of Ensembl-ID genes, share of the universe in the worst tie, gold edges in the worst tie, the sizes of U′/U″ and their gold counts, and Spearman(K, number of gold targets per TF).

**Tests.** Each variant × metric × arm runs through `compare_units` per pair, with `holm_adjust` over that experiment's 2 × (number of pairs) directed hypotheses, the same structure as `models_section`. The report also shows each variant's share of the real lift, with a 2,000-draw bootstrap over datasets.
- **Check (i):** the recomputed base and real values, plus `motif_prior` for the C4 family, are passed through `real_bench.models_section` unchanged. Its declared-comparison lines must equal the committed `REPORT_REAL.md`:
  - C4: +0.182 [+0.138, +0.227], p_Holm 1.36e-07
  - C5 RegDiffusion: +0.184 [+0.132, +0.237], p_Holm 1.18e-06
  - C5 GRNBoost2: +0.130 [+0.108, +0.153], p_Holm 3.63e-10
- **Declared primary tests** (auto arm; Holm over 3 pairs × 2 tests × 2 directions = 12):
  - **T1:** pooled(real) > pooled(perm_db), i.e. target identity beats a random assignment with the same TF-level and database structure.
  - **T2:** macro_R(real) > macro_R(base).
- Everything is labelled **exploratory**.

---

## 3. Resources (one job, QOS `stud`)

- **CPU:** 8 (`hpc/engine_bench.sbatch` defaults: 8 CPUs, 32G, 6 h, `OMP_NUM_THREADS=1`). One worker per dataset unit, 40 units (20 C4 + 20 C5).
- **Refit cost, from the committed LM-1 rows:**
  - RegDiffusion: 8,771 s (auto) + 7,185 s (baseline); the longest single dataset took 1,267 s.
  - GRNBoost2: 1,292 + 1,231 s.
  - Total ≈ 18.5k CPU-s. This drops to roughly 60% if the auto refit is reused when auto == baseline (13 of 20 LM-1 datasets), after checking that the two `processed.h5ad` files have the same sha256.
- **Analysis cost:** measured here on a 172,350-pair universe:
  - `evaluate_ranking` 0.83 s, `rank_product` 0.07 s, fast scorer 0.16 s per variant.
  - With 7 fixed variants and 2 × 100 permutations, about 60 s per context × pair, or about 3–6 CPU-h for 120–360 context-pairs. The context count is unknown; it is in `context.json`.
- **Wall time:** expected 1.5–2.5 h; request `--time=06:00:00`.
- **Memory:** 32G. The registered job, which also did screening, peaked at 11.7 GiB. Here the per-worker load is:
  - `processed.h5ad` (≤ 3,000 cells);
  - the feather columns for matched genes only (memory-mapped);
  - the `.tbl` annotation (about 0.1 GB on disk);
  - RegDiffusion training.
- **Environment:** set `TMPDIR=$HOME/gb_tmp`. The LM-1 attempt 670900 died on a full node `/tmp`.
- **Login node only runs:** `hpc.sh sync`/`check`, `mkdir`, receiving the LM-2 tarball and `sha256sum -c`, `motif_sensitivity.py plan` (lists about 40 small JSON files and checks their sha256 against the committed `SHA256SUMS`), `sbatch`/`squeue`/`sacct`, and `tar` for fetch. The refits, variants and the report run in Slurm; the report can also run on a laptop after fetch.
- **Outputs:** the rows with per-seed arrays are about 25 KB per context-pair. If the total exceeds about 1 MB, keep the `*_seeds` arrays outside git.

---

## 4. Skeleton, record and commands

**Record, committed and pushed before the run:** `experiments/EXP-20260927-LM-1.md`, or the next free LM number on the registration date. None exists today (`ls experiments/`).
- type: exploratory
- inputs: the kept artifacts of EXP-20260925-LM-2 and EXP-20260926-LM-1 (read only)
- units: the 20 datasets; RegDiffusion n = 19 by the registered rule (c)
- metrics, controls, seeds 0–99, `MIN_POS_TF` 5, tolerance 1e-12, T1/T2 with Holm over 12, and the section 5 rules
- exclusions: gate failures are reported and left out; nothing else
- results: `benchmarks/results/EXP-20260927-LM-1/`; large outputs in `~/gb_runs/EXP-20260927-LM-1`

It never writes into the LM-1 or LM-2 folders.

Ship `benchmarks/motif_sensitivity.py` together with a test built on a tiny version of the mock (gate = 0; constant-row controls give macro lift 0; tie controls give 0 on U″).

**Commands:**
```bash
gh issue create --label experiment --assignee @me --title "EXP-20260927-LM-1 motif-lift decomposition (exploratory)" --body "question, owner LM"
python -m grnbench.board        # no open running issue touching these run folders
bash hpc/hpc.sh check && bash hpc/hpc.sh sync          # clean commit -> GIT_COMMIT
# C4 inputs (on the machine holding the LM-2 run folder <LM2_RUN>):
(cd <LM2_RUN> && find rows artifacts -type f -print0 | xargs -0 sha256sum > INPUT_SHA256SUMS)
tar czf - -C <LM2_RUN> rows artifacts INPUT_SHA256SUMS | bash hpc/hpc.sh run "mkdir -p ~/gb_inputs/EXP-20260925-LM-2 && tar xzf - -C ~/gb_inputs/EXP-20260925-LM-2"
bash hpc/hpc.sh run "cd ~/gb_inputs/EXP-20260925-LM-2 && sha256sum -c --quiet INPUT_SHA256SUMS && echo inputs OK"
EXP=EXP-20260927-LM-1; A='source /software/miniconda3/etc/profile.d/conda.sh && conda activate grnbench'
bash hpc/hpc.sh run "$A && mkdir -p ~/gb_runs/$EXP ~/gb_tmp && PYTHONPATH=src python benchmarks/motif_sensitivity.py plan --run EXP-20260926-LM-1=\$HOME/gb_runs/EXP-20260926-LM-1 --run EXP-20260925-LM-2=\$HOME/gb_inputs/EXP-20260925-LM-2 --out \$HOME/gb_runs/$EXP/sens_tasks.json"
# engine_bench.sbatch calls "$SCRIPT run --tasks --out [--refs] --index --chunk --jobs --workroot": the script accepts exactly that
bash hpc/hpc.sh submit "--time=06:00:00 --export=ALL,OUT=\$HOME/gb_runs/$EXP,SCRIPT=benchmarks/motif_sensitivity.py,TASKS=\$HOME/gb_runs/$EXP/sens_tasks.json,REFS=\$HOME/grnbench_data_lm1/references,CHUNK=100000,GRNBENCH_CISTARGET_DIR=\$HOME/grnbench_data_lm1/cistarget,TMPDIR=\$HOME/gb_tmp hpc/engine_bench.sbatch"
gh issue edit N --add-label running && gh issue comment N --body "launched: HPC job <id>, 8 CPUs, ETA ~2 h, out ~/gb_runs/$EXP"
# wait with the skill's growing-interval loop; then:
bash hpc/hpc.sh submit "--partition=stud --qos=stud --time=00:20:00 --cpus-per-task=1 --mem=4G -o logs/sens_report_%j.out --wrap='$A && PYTHONPATH=src python benchmarks/motif_sensitivity.py report --out \$HOME/gb_runs/$EXP'"
bash hpc/hpc.sh run "cd ~/gb_runs/$EXP && ls rows | wc -l && find . -type f \( -name '*.json' -o -name '*.md' \) -print0 | xargs -0 sha256sum > SHA256SUMS"
bash hpc/hpc.sh fetch "~/gb_runs/$EXP" <local folder outside the repo>   # copy SENSITIVITY.md, rows/, sens_tasks.json, CODE_COMMIT, SHA256SUMS to benchmarks/results/$EXP/
```

**Script.** The full, mock-tested script is `motif_sensitivity.py` in this folder (see also
`mock_run.py`, the harness used for section 6 below). It is not wired into `benchmarks/` and has not
touched real data or the cluster.

---

## 5. Decision rules, fixed in advance (per pair, auto arm)

Definitions (means over datasets of paired differences; the registered pooled metric unless stated):
- L = pooled(real) − pooled(base). It must reproduce the registered lift.
- Δ_tgt = pooled(real) − pooled(perm_db). This is T1.
- M_R = macro_R(real) − macro_R(base). This is T2.
- L_art = pooled(tfmean_db) − pooled(base).
- Components for naming the cause: L_TF (`tf_mean`), L_pres (`tf_presence`), L_DB (`db_only`), L_tie (`worst_tie`).

| verdict | rule |
|---|---|
| **Target-level sequence evidence** | T1 supported (p_Holm(real > perm_db) < 0.05) **and** T2 supported **and** Δ_tgt ≥ 0.5·L |
| **Mostly TF-level / DB-membership artifact** | T2 not supported (or its reverse is significant) **and** L_art ≥ 0.5·L; or T1's reverse is significant (a random assignment beats the real one). The dominant cause is the largest of L_TF, L_DB, L_tie. |
| **Mixed** | T1 and T2 supported but Δ_tgt < 0.5·L, or T2 supported but L_art ≥ 0.5·L. Real within-TF sequence signal exists, but at least half of the pooled lift can be had without it; the claims must say both. |
| **Inconclusive** | anything else ("unresolved" is not "equal") |

- The baseline arm and EPR are secondary. They only have to agree in sign.
- Shares are descriptive; the verdict rests on T1 and T2. The no-motif mock showed why: the real prior's pooled lift was near zero, while the artifact controls were strongly positive.
- Because this is exploratory, it cannot overturn C4/C5 as registered: the registered metric is reproduced by design. What it can change is how they are worded: "motif prior raises pooled AUPRC" vs "sequence evidence".
  - A "mostly artifact" or "mixed" verdict leads to a note in `paper/CLAIMS.md`, and to a confirmatory fresh-set run with macro_R as primary. That run would need a D-xxx, because the metric changes.

---

## 6. What was checked here (synthetic only, not evidence about real data)

**As first reported by the agent that wrote the script** — **not reproduced on a fresh rerun, see the
GRNBoost2 caveat in section 0**: on a planted-motif mock (4 datasets; Pearson and GRNBoost2 built through
`real_bench.score_models`, with the cap patched to 3,000 so truncation is exercised), all 16 gate checks
were reported as `|diff| = 0.0`, refit top-3,000 and combined top-3,000 matched the kept files exactly,
the Check (i) lines were identical to `real_bench.py report --out` on the same rows, and:
- GRNBoost2: T1 = +1.2741 (p_Holm 9.36e-04), T2 = +1.0671 (p_Holm 0.0106); L_art (`tfmean_db`) = −0.481. Verdict: target-level, as planted.

Re-running an equivalent 4-dataset planted-motif mock for this report, the Pearson gate reproduced
(`|diff| = 0.0`) but the GRNBoost2 gate did not (0.19-0.31 off); the T1/T2/L_art numbers above should
therefore be treated as **not yet verified** until the GRNBoost2 gap in section 0 is fixed and this is
rerun. The Pearson-only result below is unaffected (no GRNBoost2 involved):

**No-motif mock** (6 datasets; random motifs; motif count K tied to TF degree, Spearman 1.0; 40 Ensembl-ID genes and 70 genes not in the database per context; 5 motif-less TFs; Pearson only):
- `motif_prior` alone scored 1.2233, with no motif information at all.
- pearson_motif − pearson_signed: +0.029 [−0.013, +0.071], unresolved.
- Controls: `tf_mean` +0.313, `tfmean_db` +0.278, `worst_tie` +0.145, `tf_presence` +0.107, `db_only` +0.039 (all p_Holm < 1e-5).
- T1 = −0.1411 (reverse p_Holm 6.19e-04); T2 = −0.3895 (reverse p_Holm 1.75e-04). Verdict: artifact, as built.

**Mechanism demos** (consistent with section 3 of `AUDIT.md`, one folder up):
- `../repro/db_membership_confound_toy.py`: `motif_prior` with random motifs scored 1.065 with no
  Ensembl-ID genes and 1.211 when they made up 15% of genes.
- `../repro/tf_level_motif_confound.py`: random motifs lifted `pearson_motif`'s pooled AUPRC ratio by
  +0.174 when the number of annotated motifs per TF is correlated with CollecTRI degree, and by −0.093
  when it is independent of degree.
