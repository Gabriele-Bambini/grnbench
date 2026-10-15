# grnbench: independent audit, 2026-09-27

Commit audited: `13106b1` (tip of `main`, includes PR #14). Scope: whole repository (engine, inference,
evaluation/statistics, benchmark harness, reporting, HPC scripts, docs/process), with emphasis on the
two committed confirmatory claims (C4, C5) and everything that could affect a future confirmatory run.

**Verdict: GO, with wording corrections to C4/C5 and a short list of fixes before the next confirmatory
run.** No finding survived adversarial verification at `blocker` or `high` severity. The registered
headline numbers of C4 and C5 are correct and reproduce exactly from the committed rows. The single
finding that matters for how the paper should read C4/C5 is mechanistic, not numerical: the motif-prior
lift is partly explained by two model-design effects unrelated to sequence evidence at the target
(section 4). This does not reverse either claim; it changes what they may say caused it.

## 0. Method

Two rounds of an independent, read-only audit, each with adversarial verification, run against a fresh
venv built from `requirements/lock-py311.txt` (`.venv/`, no edits to the working tree):

- **Round 1** — 8 area finders in parallel (statistics/eval, experiments & paper claims recomputed from
  result files, benchmark-harness design/leakage, engine, inference adapters/plugins, reporting/viz/UI,
  HPC scripts/security/secrets, docs/process) plus a completeness critic. Every finding above `low`
  severity, and every `low` finding as a batch, was re-verified by a second, independently-instructed
  agent that tried to refute it with its own reproduction. 54 agents, 1,614 tool calls.
- **Round 2** — the critic's own findings were each checked by two more agents under different lenses
  (facts vs. statistics, mechanism vs. real-data evidence); a dedicated pass reconstructed the
  pre-registration history of both confirmatory PRs from GitHub's server-side Activity API (branch/push
  timestamps are not client-settable, unlike git author dates) rather than trusting local git log; and a
  fourth pass designed and locally validated (on synthetic mocks, not real data — the HPC and the
  original datasets are unreachable from this environment) a sensitivity analysis to size the mechanistic
  effect on the real, already-collected run artifacts. 7 agents, 430 tool calls.

Total: 61 subagents, ~7.9M tokens, no repository writes during the audit (`git status` stayed clean
throughout; the only local artifacts are under this folder).

All quoted numbers below were reproduced by at least one agent from the actually committed files or by
running the actual `grnbench`/`real_bench` code, not transcribed from a report. Reproduction scripts for
the most important findings are in `repro/`; a drafted, mock-tested (not real-data-tested) sensitivity
analysis is in `proposal/`.

## 1. Test suite (fresh venv from the lock)

```
python -m venv .venv && . .venv/bin/activate
python -m pip install -r requirements/lock-py311.txt && python -m pip install --no-deps -e .
python -m pytest -q -p no:cacheprovider
```

**`395 passed, 1 skipped, 439 warnings in 162.46s`** (full log: `pytest.log`). `pip check`: no broken
requirements. The one skip is `test_motif_methods.py::...RegDiffusion...` (`torch`/`regdiffusion` not in
the lock, as documented). This runtime is Linux/Python 3.11.15; AGENTS.md's floor was recorded on
Windows/3.11.9 and predates several merged PRs, so the useful comparison is failing test IDs, not the
count (`370` vs `395`, difference explained by tests added since). The floor line and
`docs/START_HERE_AGENT_PROMPT.md`'s "expect exactly the floor" wording are stale (see 7.6); this PR
does not update either.

## 2. Status of the 2026-09-13 audit (`audit/2026-09-13/`)

All 13 blocking/high findings of that audit are fixed at HEAD, each re-verified with an adapted fixture
against the current code (`scratchpad` scripts referenced by the eval-stats agent, not shipped here —
see `findings.json` for the itemised check): EPR k frozen before gold removal, degree-preserving null
without loops/duplicates, AP with a declared universe keeps missed positives in the denominator, MI/CLR
correct on constant genes and invariant to row order, subsampling relabelled as a descriptive sensitivity
(not a bootstrap CI), tier logic replaced by pairwise `compare_units`, the Ridge plugin excludes the
target from its own regressors, BEELINE gene selection matches the pinned reference, TF case-folding,
and directory fingerprints are now full recursive inventories. One of the fixes (the top-k "boundary tie"
policy) is the direct cause of a new medium finding below (7.5).

## 3. Headline finding: the motif-prior lift is partly a model-design artefact, not only sequence evidence

This is the one finding that should change how the paper describes C4 and C5. It does **not** change
whether they are true: every registered number reproduces exactly (section 5), and both effects below
stay well short of erasing either lift.

### 3.1 Mechanism 1 — TF-level reweighting

`motif_ranks` (`src/grnbench/infer/cistarget.py:283-290`) scores every `(TF, gene)` pair as the **best
rank over all motifs annotated to that TF** (`ranks[rr].min(0)`, minimum over databases). The expected
minimum of *k* random ranks falls as *k* grows, and a TF with no annotated motif is tied at the worst
rank for every gene. `motif_rank_product` (lines 316-325) never renormalises within a TF before pooling.
Because `split_reference` (`src/grnbench/screen.py:36-47`) halves each TF's edges, the evaluation-half
degree of a TF tracks its full CollecTRI degree — so **any TF-level ordering that tracks reference
degree raises the pooled AUPRC ratio with zero information about which target a TF regulates.**

Reproduction (`repro/tf_level_motif_confound.py`, uses the real `motif_rank_product`/`motif_prior`/
`evaluate_ranking`, a fake feather DB whose motif rankings are pure random permutations — zero edge-level
information): with the number of annotated motifs per TF correlated with CollecTRI degree, `pearson_motif`
gains **+0.174** pooled AUPRC ratio (10/10 replicates positive) and `motif_prior` alone reaches 1.335 —
about the size of the registered C4 lift (+0.182) — from a prior that carries no target information at
all. With the correlation removed, the same construction gives **−0.093**. The real committed
`viz/*_false_positives.tsv` of EXP-20260925-LM-2 show the same signature on real data: `motif_prior` has
*lower* within-TF hits@20 than `pearson_signed` on the audited datasets but a *higher* pooled AUPRC ratio,
which a purely target-level mechanism would not produce.

### 3.2 Mechanism 2 — genes absent from the cisTarget rankings can never be positives

Only genes matched (case-insensitively) in the rankings feather get a real motif rank; everything else —
including genes `ingest.py:311-320` keeps as bare Ensembl IDs when no usable symbol exists — is tied at
the worst rank for every TF. `real_bench.py:81-85` builds the gold by upper-cased symbol match, so a gene
with no symbol can **never** be a CollecTRI positive. Demoting a target that can never be a true positive
cannot lower AP and typically raises it.

This was checked directly on real data, not only a toy: 12 of the 13 human LM-2 datasets were
independently refetched from the CELLxGENE census (a host reachable from this environment) and scored
with a CollecTRI snapshot from `decoupler`; the reconstruction matched the committed baseline-arm Pearson
AUPRC ratio to within 0.001 and GRNBoost2 exactly where it was also refit. On those 12 datasets,
**membership-only demotion (a control that removes just the "absent = worst rank" effect, keeping
everything else) accounts on average for about 21% of Pearson's observed lift (+0.042 of +0.202) and
about 15% of GRNBoost2's (+0.0195 of +0.131) — roughly 14% of C4's pooled lift (+0.025 of +0.182) —
and does not exceed half the lift except on 2 of the 12 datasets.** Real-data reproduction:
`repro/db_membership_confound_realdata.py`; isolated-mechanism toy (used to first size the effect
in the abstract, not as a real-data estimate): `repro/db_membership_confound_toy.py`.

### 3.3 What this means, and what it does not

- **C4 and C5 as literally worded ("adding the motif prior raises the AUPRC ratio on the CollecTRI
  evaluation half") are supported.** No committed number is wrong.
- **The mechanistic reading in `docs/models_that_work.md` and the EXP-20260925-LM-2 verdict ("sequence:
  whether the TF's binding motif sits near the target", "expression adds to sequence") is not supported
  by the current design.** Neither declared secondary comparison (`pearson_motif` vs `motif_prior`) can
  separate a target-level effect from the two mechanisms above, because `motif_prior` carries the same
  TF-level and DB-membership components.
- There is real, committed evidence of a genuine within-TF, target-level component too (hits@20 within a
  TF rise for `pearson_motif`/`grnboost2_motif` over the base model on every audited dataset), so the
  correct description is "mixed", not "artefact only" — but the design as registered cannot say how much
  of each.
- **Recommended immediate action:** add one caveat sentence to C4/C5 in `paper/CLAIMS.md` and to
  `docs/models_that_work.md` — the pooled AUPRC-ratio lift is not yet shown to be target-level sequence
  evidence; a TF-level and a DB-membership component are known to move the same metric in the same
  direction, at a magnitude comparable to the registered lift.
- **Recommended follow-up:** a small, cheap **exploratory** sensitivity experiment that reuses the
  artifacts already produced by C4/C5 (no refetching, and no refitting except the two learned base models)
  to decompose the lift with a within-TF macro metric and a target-identity-permutation control. Designed,
  and validated end-to-end on synthetic mocks, in `proposal/`; see section 6. It needs real data and HPC
  access this environment does not have (section 8).

## 4. Pre-registration & claim-integrity findings on C4/C5

**The pre-registration itself is sound — this is the good news the round-1 finders could not establish.**
A `git log` of this clone cannot show it, because both PRs were squash-merged and their branches deleted
(the repo's own rule, AGENTS.md §5): none of the branch commits the records cite (`67ca880`, `84f28e1`,
`8e091f6`, `eee1171`, …) exists as an object in a plain clone of `main`. But GitHub keeps them reachable
under `refs/pull/N/head`, and GitHub's server-side Activity API records branch-creation and push events
with timestamps the client cannot set (unlike git author/committer dates). Reconstructing both PRs from
that API and from `refs/pull/7` / `refs/pull/14`:

- **EXP-20260925-LM-2** (PR #7): registration commit `8e091f6` pushed 2026-09-25T16:25:29Z (server time).
  Section 1 (hypothesis, endpoint, arm, family, exclusions, dataset list, seeds) is **byte-identical**
  from registration through the `main` squash, except one editorial restatement of the same H1 pushed
  34 min later, still before any data was scored (first results push: 17:55:18Z).
- **EXP-20260926-LM-1** (PR #14): registration commit `eee1171` pushed 2026-09-26T19:32:39Z. Section 1 is
  **byte-identical** (same md5) at all six points checked between registration and `main`. The one
  substantive-looking addition — the "ok-eol" clause of exclusion rule (a) — was already in the
  registration commit itself, more than 5 hours before the cluster data check that needed it.
  Without it, the check would have excluded all 20 datasets (Windows-CRLF-hashed `truth.json` vs. a
  Linux fetch).

**Two minor issues remain**, both about claim wording/traceability rather than the result:

1. **`EXP-20260925-LM-2.md` misstates its own registered sample**: "human 16 / mouse 4, 8 protocols";
   the dataset list committed in the same registration commit is **13 human / 7 mouse, 7 assays**
   (`experiments/EXP-20260925-LM-2.datasets.tsv`, `tasks.json` agree). This wrong text was carried
   verbatim into the generated `benchmarks/results/EXP-20260925-LM-2/latex/report.tex` Methods section
   on a later regeneration. *Fix: correct the record with a dated erratum; regenerate the LaTeX.*
2. **No document tells a reproducer how to fetch the pinned `CODE_COMMIT`** now that the branches are
   gone (`git fetch origin pull/14/head && git checkout 84f28e1`, or a per-experiment tag pushed before
   deleting the branch). `main` has since changed RegDiffusion's zero-variance handling (#16) in exactly
   the code path that excluded one C5 dataset; a reproducer who checks out `main` instead of the pinned
   commit would silently score a different `n`. *`.projectos/STATE.md` step 3 ("reproduce
   EXP-20260926-LM-1... Issue #13 is fixed (#16)") should say explicitly: reproduce from the pinned
   commit, not from `main`.*

**A related, separate concern about the unit of analysis** (found by the round-1 completeness critic,
confirmed by two more agents): of the 20 "fresh" datasets, **3 pairs come from the same source study**
(Alzheimer's EC/SFG astrocytes, Patch-seq M1 morphological/electrophysiological, Tabula Muris Senis
kidney/diaphragm — confirmed by matching titles/gene spaces in the committed manifest), and one
confirmatory dataset (`341a0702`, Patch-seq, 636 cells) is very likely a near-duplicate of an exploratory
dataset that helped motivate H1 (`8f98c236`, 628 cells; same gene space and reference-edge count; scores
across 14 fields differ by ≤0.024, the closest pair by far in either set). `benchmarks/fetch_cellxgene.py`
excludes only by dataset-id prefix and never records the CELLxGENE *collection*, so this kind of overlap
cannot be caught mechanically. **The registered conclusions are robust to it**: collapsing each same-study
pair to one unit gives H1 +0.184 (n=17, p_Holm 1.12e-06), C5-GRNBoost2 +0.129 (n=17, 1.01e-08), C5-RegDiffusion
+0.192 (n=17, 2.42e-06) — all essentially the registered values — and dropping either member of every pair,
in all combinations tried, never lowers a result below p_Holm ≈ 2e-05. *Fix: keep `collection_id` in the
fetch manifest and in `experiments/*.datasets.tsv`; exclude by collection, not by id, on any fresh-set
rerun; report `n` as independent studies (or re-run the paired test at the study level) as a declared
sensitivity check; and word C4 as "20 previously unscored dataset ids" rather than "20 fresh… datasets"
until that is done.*

## 5. Numbers, reports and reproducibility — all checked, all correct

- `sha256sum -c` passes for every file in `SHA256SUMS` of both `benchmarks/results/EXP-20260925-LM-2/`
  and `.../EXP-20260926-LM-1/` (73 and 78 files).
- Every headline number of C4 and C5 (mean paired difference, 95% CI, n, p, p_Holm, per-pair and
  reversed) was recomputed two ways — via the repo's own `real_bench.py report` and via an independent
  `scipy`/hand-written Holm implementation — from `all_rows_real.csv`/`rows/*.json`, following exactly
  the pre-registered endpoint/arm/aggregation/family/exclusions. Both match each other and match
  `REPORT_REAL.md`, the LaTeX tables, `paper/CLAIMS.md`, `CHANGELOG.md` and `experiments/README.md`.
- Regenerating `REPORT_REAL.md` and the LaTeX tables into a scratch copy reproduces them byte-for-byte
  (line-ending differences only, which `.gitattributes` already normalises for the `.md`). One exception:
  `latex/report.tex` for both experiments is **stale text** — it was generated before later record edits
  and still prints "Reproduced by: nobody yet" although both records record an agent reproduction; the
  *numbers* in it are unaffected. *Fix: regenerate `latex/` for both experiments (`grnbench report`).*
- The test that describes the procedure as "two-sided" (`paper/CLAIMS.md`, `texreport.py:277`, both
  records) is imprecise: the code runs a **one-sided test in each direction**, Holm-adjusted over both
  directions, which is what actually produced every quoted `p_Holm`. This is more conservative than a
  genuine two-sided test would be, so no number is inflated, but a reviewer who reruns "two-sided t,
  Holm over the pairs" as literally written gets a different (still correct-favoring) p. *Fix: reword to
  "paired t, one-sided in each direction, Holm over the 2×pairs directed hypotheses" everywhere it
  appears.*

## 6. Proposed follow-up: sensitivity analysis for section 3, ready to review

`proposal/motif_sensitivity.py` (+ `proposal/mock_run.py`) is a full draft that:

- reuses only the artifacts `real_bench.py --keep-artifacts` already wrote for C4/C5 (no new download
  except copying C4's artifacts, which live outside git on the machine that ran them, to the cluster);
  refits only the two learned base models (RegDiffusion, GRNBoost2 — Pearson is closed-form) with the
  registered seeds, gated on exactly reproducing the kept top-20,000-edge files and the committed
  per-dataset rows before trusting anything derived from the refit;
- computes, per dataset/context/declared pair: the registered pooled metric (must match `REPORT_REAL.md`
  exactly — this is the calibration check), a within-TF macro metric restricted to targets that are
  labelable and DB-present (on which neither mechanism above can act, by construction), and named
  controls that isolate each mechanism (a target-identity permutation within each TF's in-database
  columns; TF-level-only constant-row variants; DB-membership-only variants);
  Holm-corrects two pre-declared tests (target-identity beats permutation; within-TF lift is real) over
  the resulting family, exactly via `compare_units`, and is explicitly `exploratory` throughout;
- was exercised on two synthetic mocks (`SENS_MOCK=1`) by the agent that wrote it: a "planted-motif" mock
  (real target-level signal) correctly called "target-level", a "no-motif" mock (pure TF-level +
  DB-membership confound, zero target information) correctly called "artefact". **Re-running that
  validation myself before writing this report, the Pearson-based calibration gate reproduces cleanly
  (`|diff| = 0.0` against the mock's own `real_bench`-generated rows), but the GRNBoost2-based gate does
  not** — on a fresh 4-dataset planted-motif mock, the recombined `grnboost2_motif` pooled AUPRC ratio
  differs from the row it should match by 0.19-0.31 (`pearson_motif`'s equivalent check passes). The
  script's own internal self-consistency assertion (its local rank-product reimplementation must equal
  `cistarget.motif_rank_product` given the *same* base edges) held throughout, so the gap is not in that
  reimplementation; the most likely cause is that `with_motif_prior` re-invokes the GRNBoost2 estimator
  independently when building the `*_motif` combination, and that second fit did not reproduce the
  standalone `grnboost2` fit bit-for-bit in this environment's thread/dask configuration (both nominally
  seed=42). This is *not* evidence that the real C5 run is unreliable — an independent agent already
  reran the actual HPC job end-to-end and got identical numbers (`.projectos/STATE.md`) — but it means
  **this draft's GRNBoost2 path is not yet trustworthy and needs to be fixed and re-validated (ideally
  against the real, already-reproduced C5 artifacts) before anyone runs it on the cluster.** The Pearson
  path (which is what C4 needs) has no such gap. **Nothing here has touched real data or the cluster** —
  that needs HPC access this environment does not have.

Estimated cost: one `stud`-QOS job, 8 CPUs, ~1.5-2.5h (dominated by the RegDiffusion/GRNBoost2 refits,
whose per-dataset cost is already known from the LM-1 job). Full design, exact commands, and the decision
rule fixed in advance are in `proposal/PLAN.md`.

This is a proposal, not a registered experiment: per AGENTS.md §7, whoever runs it opens the board issue,
commits `experiments/EXP-<date>-<initials>-n.md` from the plan before launching, and reviews the script
first (it has not been reviewed by a person, and touches `src/grnbench/eval`-adjacent logic).

## 7. Other confirmed findings (medium severity)

Everything here reproduces from a concrete command or a tiny script (see `findings.json` for the full
evidence text, verifier reasoning and suggested fix of every item, including the 61 low-severity ones
omitted here for length). None of these affects a **committed** number; they are risks for the *next*
run or a real, if narrower, bug today.

| # | file:line | finding | why it matters |
|---|---|---|---|
| 7.1 | `benchmarks/real_bench.py:203-216` | A model that errors in *some* (not all) contexts of a multi-context dataset still gets a finite dataset score, averaged over the surviving contexts only — contradicting the registered exclusion rule (c) ("a dataset where a model errors is dropped from the pairs involving that model"). Independently found by 4 of 8 area auditors. Latent in the committed runs (the only model errors were whole-dataset), so C4/C5 are unaffected. | The next multi-context confirmatory run (any per-cell-type context that fails for one model only) would silently keep a mismatched-context pair instead of dropping it as pre-registered. |
| 7.2 | `src/grnbench/qc.py:100` (+`pipeline.py:46`, `hygiene.py:328`, `profile.py:70`) | Detection-based QC (genes/cell, cells/gene, sparsity, doublet heuristic, technology guess) counts **stored sparse entries**, not nonzero values. An AnnData with explicit stored zeros (common after some scanpy operations) sails through every QC gate with the wrong metrics; a 300×2000 counts matrix that should fail `min_genes` on 299/300 cells instead keeps all 300 when zeros are stored explicitly. | Silent, wrong QC on any input with stored zeros; no test exercises this case. |
| 7.3 | `src/grnbench/autoprep.py:47-48,349` | The high-mitochondrial-tissue relaxation (A-020) matches by bare substrings (`"renal"`, `"ventric"`, `"septum"`, `"atri"`), so CELLxGENE UBERON labels like "adrenal gland", "subventricular zone" or "lateral septum" match by accident and silently get the 60% mito gate meant for liver/heart/kidney/muscle. | A brain or adrenal dataset with genuinely high mito fraction (e.g. a QC problem) passes instead of failing loudly, against D-004's stated default. |
| 7.4 | `src/grnbench/autoprep.py:450-451` (+`profile.py:186-189`) | `NaN`/`"unknown"` cell-type labels become their own per-cell-type contexts (and are counted as real cell types by the selector), although the engine already has a shared `hygiene.missing_label_mask` used everywhere else for exactly this purpose. | Context planning treats unlabelled cells as a stratum; real cell types can be pushed out of the "large enough" set by inflated counts. |
| 7.5 | `src/grnbench/autoprep.py:396` (+`kb.py:116-125`) | The selector can choose a KB option marked `implementation="external"` (e.g. `scDblFinder`, `emptydrops`) whenever its declared R/py dependency happens to be importable, but `plan_preprocessing` only knows how to execute a handful of built-ins — an unrunnable winner is silently dropped from the plan (e.g. **no doublet removal at all**) rather than falling back. | In an environment with optional R packages installed, the engine can silently run a weaker pipeline than the one it reports choosing. |
| 7.6 | `src/grnbench/selector.py:189-201` | The context-scope alternative `Selection.candidates(always=...)` is forced to screen when `per_cell_type` wins is `selected_type` — which `autoprep.py` builds as the *same first context* `per_cell_type` already runs, with identical `run_prep` arguments. The two-way "recommended vs. always-varied" screen is a no-op whenever per-cell-type wins (962 of the screened runs in `benchmarks/results/`, vs. 12 that varied `pooled`). | The screening/comparison the product's UI and papers describe (e.g. claim C3, "auto ties baseline") never actually exercises a materially different context choice in the overwhelming majority of runs. |
| 7.7 | `src/grnbench/infer/classical.py:235` | The default (`alpha=None`) `lasso` path fits `lasso_path` on **uncentered** log-normalized data; scikit-learn's `enet_path` fits through the origin in that mode, so the first regulators selected are whichever have the largest raw mean, not the largest correlation. On a toy with a true regulator (corr 0.90) and an independent high-mean gene (corr −0.01), the independent gene ranks first (0.55 vs 0.0); centering fixes it (0.79 vs 0.0). On 5 planted networks, AUPRC ratio is 5.4 as shipped vs. 8.0 centered (Pearson: 8.8). | `lasso` is one of the built-in methods in `EXP-20260917-PRE-1` (the BEELINE backtest, `benchmarks/results/beeline-2026-09-17-9methods/`); its rows likely understate what a correctly-centered Lasso would score. |
| 7.8 | `src/grnbench/infer/base.py:222-228` | `run_method` maps a Python plugin's output columns **positionally** (`raw.columns[0..2]`) and reads direction only from a column literally named `direction`; the alias-matching README promises (`regulator`/`tf`/`score`/`sign`, case-insensitive, any order) is implemented only for command-line plugins. A plugin that follows the documented Python contract but orders columns `[target, source, score]` gets every edge silently reversed; a `sign` column is silently dropped. | Anyone writing a Python plugin per the README (the "bring your own model" comparison feature, D-012/D-013) can get silently reversed edges with no error. |
| 7.9 | `src/grnbench/infer/base.py:69-74,211-212` | `prepare_adata` returns the caller's AnnData object **unchanged** (same object, not a copy) whenever the data is already normalized; `run_method` passes that object straight to any `takes="adata"` plugin — i.e. every Python plugin. A plugin that mutates its input in place (a completely ordinary scanpy pattern, e.g. `sc.pp.log1p(adata)`) changes the results of every method scored after it in the same run. | Silent cross-contamination between methods scored in the same `grnbench run`/`real_bench` invocation whenever a plugin mutates in place; order-dependent, hard to notice. |
| 7.10 | `src/grnbench/infer/regdiffusion_adapter.py:31` | The #13/#16 zero-variance guard drops genes with zero **raw** variance, but RegDiffusion 0.2.3 itself min-max-scales *per cell* before checking variance, so a gene that is the top-expressed gene in every cell (e.g. a marker like `HBB`/`ALB`) varies in raw data but is exactly 1.0 after RegDiffusion's own scaling — the guard doesn't catch it, and RegDiffusion still raises. A worse case (an all-zero cell over the retained genes) produces `NaN` throughout and the guard's check silently passes. | The fix that was supposed to close issue #13 doesn't fully close it; a future dataset can still hit the same class of failure the fix was meant to prevent. |
| 7.11 | `src/grnbench/viz.py:305-320,469` | `pr_curves`/`calibration` label a line "prevalence (random)" that is actually the mean, across models, of each model's *own* prevalence within its own top-20,000-edge dump — not the shared candidate-universe prevalence `evaluate_ranking` actually uses for the scored AUPRC ratio/EPR. In the committed `EXP-20260926-LM-1` figure for the SFG dataset, a model scored *above* random (AUPRC ratio 1.05) visibly plots *below* the dashed "random" line. | The committed exploration/explainability figures (mandatory per AGENTS.md §7) can visually contradict the very numbers they're meant to illustrate. |
| 7.12 | `src/grnbench/viz.py:72-82,597` | The "after" panel of the mandatory QC-before/after figure reads `processed.h5ad` — HVG/TF-subset genes, ≤3,000-cell-subsampled, log-normalized — and reports its column sums as if they were counts, and its gene/cell counts as if cleaning had produced them. On identical cells with no filtering, "after" total counts read 2× "before" (log-normalized values summed, not counts); the committed `EXP-20260926-LM-1` SFG figure shows total counts, genes and %mito all moving in ways cleaning never caused. | Same category as 7.11: a mandatory explainability figure shows an artefact of how it was built, not what the pipeline did. |
| 7.13 | `src/grnbench/texreport.py:145-176` (+`real_bench.py:375-386`) | Both statistics paths read the Holm family from the **current** module constant `real_bench.DECLARED_PAIRS` at report-generation time, not from anything frozen at registration or stored in the run folder. Appending one pair to that constant and regenerating `EXP-20260925-LM-2`'s report from its unchanged rows changes `p_Holm` for the registered comparison (1.36e-07 → 2.04e-07). No test pins the family. | AGENTS.md §7 requires regenerating reports "from the unchanged rows" after an analysis-code fix — which is exactly the operation that can silently change a confirmatory `p_Holm` today. |
| 7.14 | `src/grnbench/ui/app.py:207-221,373-377` | In the interactive app, editing any pipeline parameter after adopting a screened plan resets the session's reference-split state but keeps the adopted plan (which was screened *on* that split); the benchmark page then falls back to the **full** reference — including the selection half the plan was chosen on — while its own caption says "the pipeline was not chosen on it". | A user benchmarking their own model through the GUI (the product's core comparison feature) can get a benchmark run on data the shown pipeline was implicitly selected against, mislabelled as clean. |
| 7.15 | `benchmarks/real_bench.py:330-339,544` (+`engine_bench.py`) | Both HPC runners always `return 0` from `main()`, and skip any task whose `rows/<id>.json` already exists — including one that recorded an error. `hpc/engine_bench.sbatch` passes that exit code straight to Slurm. **This already happened**: job 670900 of EXP-20260926-LM-1 `COMPLETED` in 110s with all 20 rows as `OSError: No space left on device`; a person had to notice and move the error rows aside by hand before resubmitting. | A future job that errors on every task reports Slurm success; a routine resubmission (the documented recovery path) would silently skip the error rows forever instead of retrying them. |
| 7.16 | `hpc/hpc.sh:71-80` (+`engine_bench.sbatch:35,59`) | `sync` untars the working tree **in place** over the live checkout a currently-running job's spawned worker processes import from (`get_context("spawn").Pool(..., maxtasksperchild=1)` starts a fresh interpreter per task), and `CODE_COMMIT` is written once at job start — so a `sync` run while a job is in flight can make later tasks of that job silently run different code than the one recorded as `CODE_COMMIT`. Rows record no commit of their own. | Given QOS `stud` allows only 1 running + 1 queued job, this is exactly the kind of race a person preparing the *next* experiment while the current one is still running could trigger. |
| 7.17 | `hpc/beeline_backtest.sbatch:26` (+`grnbench.sbatch:34`) | These two `.sbatch` scripts locate the repo via `dirname "${BASH_SOURCE[0]}"`, which under Slurm is the **spool copy** of the script, not the submission directory (`engine_bench.sbatch` already avoids this via `$SLURM_SUBMIT_DIR`). No test exercises either script from a copied location. | Both scripts would fail to find `benchmarks/beeline_backtest.py`/write to the intended `OUT` the first time someone runs them as documented (`sbatch hpc/<script>.sbatch` from the repo). |
| 7.18 | `experiments/TEMPLATE.md:16` (+`README.md:13`, `AGENTS.md:138`, `docs/START_HERE_AGENT_PROMPT.md`) | Four docs tell the next experiment author to run the declared-family Holm test "via `grnbench compare`" — but that CLI Holm-corrects *every* pair of *every* method in the input and refuses any incomplete cell; it cannot express "the 4 declared directed hypotheses" and, run as documented on the committed LM-1 rows, exits 2 ("missing unit/method scores; incomplete pairs cannot be silently dropped"). The path that actually produced C4/C5 is `real_bench.models_section`/`texreport.comparisons`, which none of the four docs name. | The documented procedure for the statistical test at the center of every confirmatory claim does not run on real confirmatory data; the next person to register an experiment from the template would hit this immediately. |

## 8. Secrets and repository hygiene (AGENTS.md §11)

A full-history scan (`git log -p --all`, `git grep` over every blob) found **no** committed keys, tokens,
passwords, `.env.hpc`, HPC usernames or student IDs. Two pre-existing, low-impact items, reported per
§11 rather than fixed (fixing either needs history rewrite, which conflicts with the no-force-push rule,
so this needs Luca's decision, not a silent new commit):

1. **A personal email address is in the author/committer metadata** (not file content) of the first four
   commits (`4de687d`, `782da2d`, `4f2e6a6`, `f279c10`) — before the repo switched to the GitHub noreply
   address. Only matters if the repo is ever made public (D-018 anticipates this).
2. **A masked machine hostname is in `audit/2026-09-13/pytest.xml`**'s JUnit `hostname` attribute, present
   since the "sanitized" snapshot commit `4de687d` (D-016 says that snapshot was sanitized of host
   names — this one slipped through). Low impact (a vendor-default mini-PC name), but contradicts D-016.

No other absolute personal path, host name or credential was found anywhere in history. (The HPC login
and compute node names in `docs/HPC.md`/`docs/START_HERE_AGENT_PROMPT.md` are the shared cluster's own
documented names, not personal machines, and are intentionally documented there.)

## 9. Low-severity findings

61 further low-severity findings were confirmed (documentation contradicting code, stale generated files,
minor plugin/CLI robustness gaps, `--help` not honoured by several scripts, a second unreferenced lock
file, `data_manifest.py --check` not catching a missing data root, LaTeX escaping gaps, etc.). Each has
the same evidence/reproduction/fix format as above; full text, file:line, and the independent verifier's
reasoning for every one is in `findings.json` (`severity: "low"`, `verdict: "confirmed"`). None affects a
committed number or a headline claim. Two worth flagging for the next PR that touches process docs:

- `AGENTS.md:79` — two early commits changed AGENTS.md's own review/CI rules directly on `main`, before
  the branch/PR workflow the same file establishes existed (3-4 minutes after the bootstrap commit, before
  any PR). Not repeated since.
- `AGENTS.md §9` — the cisTarget motif database (1.3GB, backing both C4 and C5) is fetched by
  `hpc/fetch_cistarget.sh`, which is not `hpc/fetch_engine_datasets.sh`, the only script §9 names as
  "approved"; no licence is recorded for it anywhere `paper/README.md`'s checklist says to look.

## 10. Recommendations

**Before the next confirmatory run:**
- Fix 7.1 (partial-context scoring) — it will bite the first multi-context dataset with a partial model
  failure, and is cheap to fix (`score_models` already computes what it needs; it just needs to `NaN` the
  whole dataset score for that model instead of averaging survivors).
- Fix 7.13 (freeze the Holm family) and 7.15 (runner exit code) before relying on either path again for a
  new confirmatory experiment — both are "the next run could silently misreport itself" bugs, not
  cosmetic ones.
- Add the wording caveat to C4/C5 (section 3.3) and correct the LM-2 demographic line (section 4, item 1)
  — both are same-PR, no-code-risk text fixes.

**Needs a team decision, not a silent fix:**
- Whether/how to run the section-6 sensitivity analysis as a registered exploratory experiment, and what
  it should change about C4/C5's wording or about the primary metric of a future fresh-set confirmatory
  rerun (a metric change needs a `D-xxx`, per AGENTS.md §5).
- How to keep a pre-registration commit retrievable after squash-merge (a per-experiment tag pushed
  before deleting the branch is the cheapest fix) — currently correct but fragile (section 4).
- The two secrets-adjacent items in section 8.
- The `EXP-20260925-LM-2` dataset-independence caveat (section 4) and how it should be worded in the
  paper, and whether a fresh-set rerun should exclude by collection.

## Contents of this folder

- `AUDIT.md` — this report.
- `pytest.log` — full suite output (section 1).
- `findings.json` — all 98 round-1 findings + 3 critic findings, each with its adversarial verifier's
  independent verdict, adjusted severity, and reasoning; the round-2 re-verification of the 5 critic
  findings; and the pre-registration forensic reconstruction (section 4).
- `repro/` — runnable, independent reproductions for the findings in sections 3-5 (each was re-run against
  this PR's commit before writing this report):
  - `tf_level_motif_confound.py`, `db_membership_confound_toy.py`, `db_membership_confound_realdata.py`
    (needs network access to CELLxGENE census + `decoupler`, not available in every environment)
  - `score_models_partial_context.py`, `lasso_no_intercept.py`
- `proposal/` — the drafted sensitivity analysis of section 6: `motif_sensitivity.py`, `mock_run.py`,
  `PLAN.md`. Its Pearson path is mock-validated; its GRNBoost2 path has a known, unresolved gate failure
  (section 6). Neither has touched real data or the cluster.
