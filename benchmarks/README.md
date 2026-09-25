# Benchmarks — replicating published results

Synthetic data you wrote yourself validates nothing. A generator that builds
targets as a linear function of TF activity will always flatter correlation
methods: grnbench's own synthetic demo has `spearman_signed` beating GENIE3,
which is the *opposite* of what the literature reports on real data. That is a
property of the generator, not a discovery.

So the real regression test is the field's reference benchmark.

## `beeline_backtest.py` — BEELINE replication

Replicates the protocol of **Pratapa et al., "Benchmarking algorithms for gene
regulatory network inference from single-cell transcriptomic data", *Nature
Methods* 17, 147–154 (2020)** using BEELINE's own expression data, its own
ground-truth networks and its own gene-selection rule (TFs + top-N most variable
genes, N = 500 or 1000).

### The falsifiable prediction

On real scRNA-seq scored against **cell-type-specific ChIP-seq**, the published
literature reports that the best result across all methods and all datasets is
roughly:

| quantity | published ceiling |
|---|---|
| AUPRC ratio vs random | **1.10** |
| Early precision ratio (EPR) | **1.20** |

That is *barely above random*. So this is a falsification test, not a
confirmation test:

> If grnbench reports an EPR of 20 on real data, the harness is leaking.
> It has not discovered biology.

The script flags any result more than 3× the published ceiling and, for EPR,
prints the hit count behind it against the count expected at random — because
on a reference with a few hundred reachable edges an EPR of 4 is eight hits
against two. A clean run is one where every number is boringly close to 1.

**The ceiling applies to cell-type-specific ChIP-seq and to nothing else.** That
is the falsification target; the other ground truths are reported beside it but
are not tested against it, because there is no comparable published figure for
them and their numbers are legitimately much higher:

| ground truth | why its numbers run high |
|---|---|
| non-specific ChIP-seq | a union over many cell types — far denser and strongly hub-biased |
| STRING | aggregates evidence channels that **include co-expression**, so scoring a correlation method against it is partly circular: a high number measures agreement with co-expression, not regulation |
| lof/gof | perturbation-derived, small and high-precision; no comparable ceiling |

An earlier version of this harness compared *every* ground truth against the
1.10/1.20 figure and flagged six mHSC results as "above the published regime".
All six were scored against the generic networks, where that number does not
apply. They were false alarms — and a false alarm is not harmless, it teaches
the reader to discount the alarm that matters.

### Running it

```bash
# 1. fetch the corpus (Zenodo record 3378975, ~270 MB) - login node
bash hpc/fetch_beeline.sh ~/beeline

# 2. run the replication
python benchmarks/beeline_backtest.py \
    --data ~/beeline --networks ~/beeline/BEELINE-Networks \
    --out runs/beeline --methods spearman_signed,genie3,grnboost2

# 3. re-run only the regime check + dataset-unit comparison on a finished run
python benchmarks/beeline_backtest.py --data x --networks x --out x \
    --summarize runs/beeline/beeline_backtest.json     # rewrites summary.md
```

On SLURM: `sbatch --export=ALL,BEELINE=$HOME/beeline hpc/beeline_backtest.sbatch`.
GENIE3 across all seven datasets is hours of CPU — that is what the cluster is for.

### Datasets and ground truths

Seven experimental datasets (`hESC`, `hHEP`, `mDC`, `mESC`, `mHSC-E`, `mHSC-GM`,
`mHSC-L`). Each is scored against **every** ground truth BEELINE ships for its
organism — cell-type-specific ChIP-seq, non-specific ChIP-seq, STRING, and
loss/gain-of-function for mESC. Inference runs once per method; the same edge
list is then scored against each network, so the extra ground truths are
essentially free while the expensive part is already paid for.

### Ranking stability — and what "stability" can actually mean at n = 7

Scoring several ground truths is not thoroughness for its own sake. Recent work
(*Quantifying Ranking Instability Across Evaluation Protocol Axes*, arXiv
2603.03493) shows that GRN method rankings **flip** depending on the metric, the
ground truth, the preprocessing and the edge-ranking rule — so a benchmark that
scores one ground truth and names a winner is reporting an artifact of that
choice. In a controlled check, one unchanged edge list from one unchanged method
scored EPR 0.98 against non-specific ChIP-seq and 32.9 against the network that
generated the data — a 33× swing driven purely by the choice of ground truth.

But a flip in point estimates is only interesting if the point estimates were
resolved to begin with. The harness therefore ends with **`unit_comparison`**:
a paired t test with Holm correction in which the independent units are the
**datasets** — the only thing in this design that is replicated — one family per
ground truth × metric, and `unresolved_is_equivalent: False`. An earlier version
resampled *gold edges* instead, which is pseudo-replication: with ~4,500 edges
per reference it declared 80% of pairs resolved and 47% of them "contradictory"
across references. Both numbers are retracted (`CHANGELOG.md`, 0.7.0).

### Result of the run (7 datasets × 9 methods × 3–4 ground truths)

Six classical methods (`pearson_signed`, `spearman_signed`, `mutual_info`,
`clr`, `ppcor`, `lasso`) and the three tree ensembles (`random_forest`,
`grnboost2`, `genie3` — the latter two through arboreto, as in BEELINE); the
pinned BEELINE selection (Bonferroni p < 0.01, significant TFs + top-500
significant non-TFs by variance — verified against BEELINE's own
`generateExpInputs.py` by a test); 20 score-permutation nulls on the fixed
candidate universe. Two runs, merged after checking that every dataset's gene
set, TF set and gold sizes were byte-identical. Full outputs in
`results/beeline-2026-09-17-9methods/`.

**AUPRC-ratio vs cell-type-specific ChIP-seq** (the reference the published
regime refers to):

| dataset | pearson | spearman | MI | clr | ppcor | lasso | RF | GRNBoost2 | GENIE3 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| hESC | 1.46 | 1.45 | 1.37 | 1.25 | 1.05 | 1.11 | 1.48 | 1.39 | 1.43 |
| hHEP | 1.03 | 1.01 | 1.04 | 1.06 | 0.96 | 1.00 | 1.00 | 1.11 | 1.02 |
| mDC | 1.17 | 1.08 | 1.08 | 1.15 | 1.28 | 1.04 | 1.10 | 1.32 | 1.13 |
| mESC | 1.56 | 1.46 | 1.54 | 1.32 | 1.06 | 1.10 | 1.38 | 1.33 | 1.47 |
| mHSC-E | 1.04 | 1.07 | 1.03 | 0.91 | 1.01 | 1.03 | 1.11 | 1.01 | 1.11 |
| mHSC-GM | 1.06 | 1.03 | 1.02 | 0.92 | 1.05 | 1.12 | 1.40 | 1.11 | 1.42 |
| mHSC-L | 1.14 | 1.17 | 1.31 | 1.12 | 1.03 | 1.27 | 1.47 | 1.35 | 1.50 |

**EPR vs cell-type-specific ChIP-seq:**

| dataset | pearson | spearman | MI | clr | ppcor | lasso | RF | GRNBoost2 | GENIE3 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| hESC | 2.73 | 2.66 | 2.70 | 1.98 | 1.28 | 2.46 | 3.40 | 3.00 | 2.70 |
| hHEP | 1.07 | 0.97 | 1.31 | 1.29 | 0.82 | 0.88 | 1.31 | 1.36 | 1.26 |
| mDC | 0.00 | 1.71 | 0.00 | 4.55 | 2.85 | 1.14 | 1.14 | 2.85 | 1.71 |
| mESC | 2.22 | 2.06 | 2.23 | 1.91 | 1.10 | 1.56 | 2.22 | 2.39 | 2.47 |
| mHSC-E | 0.82 | 0.81 | 0.74 | 0.71 | 1.01 | 1.20 | 0.65 | 0.72 | 0.58 |
| mHSC-GM | 1.01 | 0.92 | 0.78 | 0.74 | 1.09 | 1.64 | 1.21 | 1.03 | 1.20 |
| mHSC-L | 1.30 | 1.36 | 1.45 | 1.17 | 1.06 | 1.83 | 1.56 | 1.54 | 1.63 |

Boringly close to 1, which is the pass condition: **grnbench does not
manufacture signal on real data.** AUPRC-ratios sit in 0.91–1.56, well inside
3× the published best. One EPR of 126 is flagged — mDC / clr at 4.55 — and the
regime check says why in the same row: mDC has 756 reachable ChIP edges among
263,220 candidates, so EPR there is **8 hits against 1.8 expected at random**,
a small-count ratio (the same row shows `pearson` at 0 hits). The flag is
correct and so is the explanation; neither is biology.

GENIE3 costs 12–24 minutes per dataset here against 2–3 seconds for Pearson
(`random_forest`, the scikit-learn re-implementation, 1.5–3 minutes;
GRNBoost2 under a minute). What that buys, at n = 7, is the next section.

Against the *generic* ground truths the numbers legitimately run much higher
(mHSC-L vs STRING: ratio 47.7, EPR 59.0 for `clr`). They are reported as
`no-published-ceiling`: STRING carries co-expression channels, so scoring a
correlation method against it is partly circular, and the pinned selection
leaves mHSC-L with 60 significant TFs and only **137** STRING pairs in 33,540
candidates — a denominator that turns a few dozen hits into a ratio of 50.

**The corrected inference, datasets as units (α = 0.05, Holm), 9 methods =
36 pairs per family:**

| reference | metric | resolved pairs |
|---|---|---|
| specific-chip | EPR | 0 / 36 |
| specific-chip | AUPRC-ratio | 0 / 36 |
| non-specific ChIP | EPR, AUPRC-ratio | 0 / 36, 0 / 36 |
| STRING | EPR, AUPRC-ratio | 0 / 36, 0 / 36 |
| lof/gof (mESC only) | — | 1 unit, no inference |

**0 of 216.** GENIE3 and GRNBoost2 — the methods BEELINE's own leaderboard
puts near the top — are not distinguishable from Pearson correlation on any
reference or metric. Meanwhile the point-estimate winner changes with the
reference on **7 of 7** datasets for AUPRC-ratio and 6 of 7 for EPR, and the
specific-chip "winner" is a tree method on five datasets and a correlation on
two. Read together: the reversals are not contradictions between references —
they are what unresolved estimates look like. Seven paired observations cannot
rank nine methods whose AUPRC-ratios differ by hundredths, and a leaderboard
that does so is printing noise. Within this design the only thing GENIE3's
extra four hundred CPU-minutes buys is a different, equally unresolved,
ordering.

A note on protocol sensitivity, because it bit this harness: an earlier run
kept *every* TF present in the expression matrix rather than only the
Bonferroni-significant ones. That over-inclusion is what BEELINE's script does
not do; with it, hESC had 2.9 M candidates instead of 373 k, and four pairs
came out resolved on AUPRC-ratio (`pearson`/`spearman` > `clr`/`mutual_info`).
Under the pinned selection those four vanish. The gene-selection rule alone
moved the only "significant" result the design had.

### What this replication has already found

Eight defects, none of which any synthetic fixture had ever caught. Each is
listed in the changelog (0.5.1, 0.7.0) with a regression test.

- **The significance layer was pseudo-replication.** Resampling gold edges and
  re-comparing methods on each draw treated ~4,500 correlated annotations as
  independent samples, so 80% of method pairs came out "resolved" and 47% of
  those orderings "contradicted" each other across references. With the
  datasets as units — the only replicated thing in the design — 0 of 216 pairs
  resolve. Removed, not deprecated; replaced by `eval.comparison.compare_units`.
- **The gene selection kept every TF, not the significant ones.** BEELINE's
  `generateExpInputs.py` filters by Bonferroni p-value *first* and keeps only
  TFs that pass. Keeping all TFs inflated hESC's candidate universe from 373 k to
  2.9 M pairs and produced the only four resolved comparisons the earlier run
  had. The selection is now pinned to the official script's SHA-256 and tested
  against its output.

- **Every BEELINE ground-truth file was unreadable** by `grnbench-eval`. They use
  the `Gene1,Gene2` convention; the alias table accepted `gene2` as a target but
  not `gene1` as a source, so all ten networks raised on load.
- **Gene symbols were matched with exact string equality.** BEELINE ships mESC
  expression in MGI case (`Pou5f1`) and every mouse network in upper case
  (`POU5F1`): 64,557 genuinely reachable gold edges scored as **0 hits**, every
  metric came back `NaN`, and nothing was printed. The worst failure mode in the
  package's whole design — a wrong number wearing the costume of an honest null
  result. Six of the seven datasets were unaffected only because BEELINE happened
  to ship them already upper-cased, which is exactly why synthetic data missed
  it. Matching now folds case (`--case-sensitive` opts out) and reports the fold;
  a `no-overlap` verdict is printed loudly instead of returning `NaN`.
- **The harness's own falsification threshold was wrong.** It applied the
  cell-type-specific ChIP-seq ceiling to *every* ground truth and flagged six
  correct mHSC results as suspicious. A benchmark that cries wolf on its own
  correct output teaches the reader to discount the alarm that matters.
- **A method matching no gold edge at all was reported as significant**, at
  `p = 0.022` with a green "significant" badge, because `null_scores >= NaN` is
  elementwise `False` in NumPy and so collected zero hits — the smallest p-value
  the estimator can emit.
- **AUPRC and EPR were reading their own sort order.** Both walked the stored row
  order, so among edges scored *identically* the arbitrary ordering decided the
  result: six permutations of the same scores gave EPR 3.0 to 5.0 while the
  tie-safe AUROC never moved. Harmless for continuous-score methods such as
  `spearman_signed` — silently wrong for any thresholded or integer-scored
  method.
- **A genes×cells CSV — the layout BEELINE itself ships — was loaded
  transposed**, because 10x barcodes match the gene-symbol pattern and both axes
  therefore scored as equally "gene-like". QC then ran on the wrong axis and
  reported `n_mito_genes: 0`, disabling the mitochondrial gate without a word.

And one finding that is not a defect:

- **BEELINE's expression inputs are already normalized**, so grnbench's
  "counts or nothing" ingest rule rejects them. The rule is right — but the
  field's standard benchmark corpus cannot enter through `grnbench prep`, so the
  documented entry point for pre-normalized data is `grnbench infer
  --assume-normalized`, which is what this harness uses. The backtest records the
  verdict per dataset rather than silently working around it.

## Do the published harnesses share these defects?

Finding six defects in your own package raises an obvious question: is the field
computing these numbers the same wrong way? Mostly **no**. The honest answer,
from reading the actual sources rather than the papers:

| harness | tie handling | gene-name case | verdict |
|---|---|---|---|
| BEELINE `BLEval/EarlyPrecision.py` | expands the cut to include every tied edge, divides by the count actually selected | n/a (closed world) | **correct** |
| BEELINE `BLEval/AUPRC.py` | `sklearn.precision_recall_curve` + `auc` | n/a | **correct** |
| [BEAR-GRN](https://github.com/UzunLab/BEAR-GRN) `R/reproduce_early_metrics.R` | `arrange(desc(score))` then `head(n_edges)` — dplyr's sort is stable, so the cut inside a tied block is decided by the row order of the input file | `toupper()` on both ends ✓ | **order-dependent, mild** |
| grnbench before 0.5.1 | sliced the first k rows | exact string match | **the severe case** |

Two things worth stating plainly because they cut against the interesting story:

- **A hypothesis that did not survive.** BEELINE integrates its PR curve with the
  trapezoid rule, and scikit-learn's own `average_precision_score` docs warn that
  this "can be too optimistic". Measured across prevalences from 0.05% to 20%,
  the trapezoid came out at 0.92–1.00x average precision — slightly
  *conservative*, never inflated. The warning does not bite in this regime, and
  the claim is withdrawn.
- **BEAR-GRN's defect is real but small.** With gold spread uniformly through the
  tied block, six orderings of identical scores moved early precision by 1.05x
  (0.098–0.103), not the 3–5x seen in grnbench. The magnitude depends on how
  concentrated the ground truth is in file order; the mechanism is there, the
  consequence is usually a few percent.

So the defects catalogued above are, with one mild exception, grnbench's own —
not the field's. That is a less exciting result than the alternative and it is
the one the evidence supports.

## Sources

- Pratapa et al. 2020, *Nature Methods* — <https://www.nature.com/articles/s41592-019-0690-6>
- BEELINE code — <https://github.com/Murali-group/Beeline>
- BEELINE inputs, Zenodo record 3378975 — <https://zenodo.org/records/3378975>
- The "best AUPRC ratio 1.10 / best EPR 1.20" figure for cell-type-specific
  ChIP-seq is quoted from *Improved gene regulatory network inference from
  single cell data with dropout augmentation*, PLOS Comput Biol (2025) —
  <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1013603>
