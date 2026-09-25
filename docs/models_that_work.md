# Models that work on real data (2026-09-25)

Question: which GRN inference methods recover real regulatory edges clearly above random
on real scRNA-seq, scored against curated TF-target collections (CollecTRI/DoRothEA) or
perturbation ground truth, and which of them use prior knowledge without leaking the
evaluation reference? Then: which 2-3 fit grnbench (open licence, pip, CPU minutes on a
few thousand cells x ~2,000 genes, the `run_method` contract)?

Measured answer on our data: `benchmarks/results/engine-2026-09-25-models/REPORT.md`.

Measured on 12 real CELLxGENE test datasets (CollecTRI evaluation half, dataset = unit,
paired t + Holm): adding the cisTarget motif prior by rank product lifts AUPRC ratio for
every model it was attached to: Pearson 1.32 -> 1.55 (+0.23 [+0.16, +0.30], p = 1e-4),
RegDiffusion 1.43 -> 1.65 (+0.22, p = 9e-4), GRNBoost2 1.10 -> 1.24 (+0.14, p = 1e-5).
GRNBoost2 alone was at random on EPR (1.02). The effect is real but still modest:
the best model reaches 1.65 x random.

## 1. What the literature says (numbers, sources)

Labels: **[read]** = checked in the primary source during this session;
**[secondary]** = from a research pass over the paper/repo, not re-read line by line.

| source | setting | finding |
|---|---|---|
| BEELINE, Pratapa et al. 2020, *Nat Methods* 17:147, doi:10.1038/s41592-019-0690-6 [secondary] | 12 methods, 5 real scRNA-seq sets (mESC, hESC, hHep, mDC, mHSC) vs cell-type-specific ChIP-seq, non-specific ChIP, STRING, LOF/GOF | EPR of most methods on cell-type-specific ChIP-seq "close to 1"; best on real data GENIE3, GRNBoost2, PIDC; PPCOR only slightly above random; SINCERITIES at or below random. Our v2 AUPRC ratio 1.2-1.3 is the field's normal range, not a pipeline bug. |
| geneRNIB, Nourisa et al. 2025, bioRxiv 10.1101/2025.02.25.640181 [read: abstract] | 12 methods x 11 real datasets, perturbation-based metrics | "simple models with fewer assumptions often outperformed more complex pipelines"; "gene expression-based correlation algorithms yielded better results than more advanced approaches incorporating prior datasets or pre-trained on large datasets". Ranking [secondary]: Pearson, GRNBoost2, SCENIC+ on top; scPRINT, CellOracle, FigR, GRaNIE, PPCOR, Portia near the bottom. |
| CausalBench, Chevalley et al. 2025, *Commun Biol*, doi:10.1038/s42003-025-07764-y [secondary] | K562/RPE1 Perturb-seq, mean Wasserstein vs false omission rate | Guanlab best; the trivial Mean Difference baseline nearly tied (K562 W=0.442/FOR=0.134 vs 0.386/0.142); interventional methods did not beat observational ones. |
| PEREGGRN, Kernfeld et al., bioRxiv 10.1101/2023.07.28.551039 [secondary] | expression forecasting after perturbation | "uncommon for expression forecasting methods to out-perform simple baselines". |
| SCENIC, Aibar et al. 2017, *Nat Methods* 14:1083; pySCENIC protocol, Van de Sande et al. 2020, *Nat Protoc* 15:2247 [secondary] | GRNBoost2/GENIE3 co-expression, then cisTarget motif enrichment keeps only targets whose promoters carry the TF's motif | The motif step is the standard way to separate direct from indirect targets. Size of its precision gain against ChIP/CollecTRI: not verified here. pySCENIC 0.12.1 (last release Nov 2022) has open install breakages with current dask/numpy (issues #561, #590). |
| RegDiffusion, Zhu & Slonim 2024, *J Comput Biol* 31:1087, doi:10.1089/cmb.2024.0607 [read: licence via GitHub API = Apache-2.0] | diffusion model over the expression matrix; learned adjacency = network | Reported better than GENIE3/GRNBoost2/DeepSEM on BEELINE-style benchmarks (exact numbers not verified here); 15k genes < 5 min on GPU. |
| CollecTRI, Müller-Dott et al. 2023, *NAR* 51:10934, doi:10.1093/nar/gkad841 [secondary] | curated signed TF-target collection (12 resources) | It is our evaluation reference, so any method using CollecTRI/DoRothEA (decoupler TF activity, networks seeded with them) is circular here. |
| LINGER, Yuan & Duren 2024, *Nat Biotechnol*; SCENIC+, Bravo González-Blas et al. 2023, *Nat Methods* [secondary] | multiome + external bulk/ATAC | Need chromatin accessibility (and for LINGER an atlas of external bulk data): out of scope for scRNA-only input. |
| scPRINT, Kalfon et al. 2025, *Nat Commun* [secondary] | foundation model, GRN from attention | GPU recommended, GPL-3.0; near the bottom in geneRNIB. |
| CellOracle, Kamimoto et al. 2023, *Nature* [secondary] | motif-scanned promoter base GRN + regularized regression | Ships promoter base GRNs for 10 species (usable without ATAC); near the bottom in geneRNIB; its benchmark is disputed (Nourisa et al. 2024, bioRxiv 10.1101/2024.10.16.618746). |
| DeepSEM, Shu et al. 2021, *Nat Comput Sci* [secondary] | VAE structural equation model | Not on PyPI; pinned to Python 3.7 / torch 1.2. |

Summary: across four independent benchmarks (2020-2025) no method clears random by a
large margin on real data from expression alone, and complex methods do not beat Pearson
or GRNBoost2. The one kind of information co-expression cannot supply is **sequence**:
whether the TF's binding motif sits near the target. That is the SCENIC idea, and it is
independent of curated TF-target collections.

## 2. Priors: allowed vs circular

| prior | where it comes from | against CollecTRI |
|---|---|---|
| TF list (regulator set) | TF catalogues | allowed; in `real_bench.py` the regulators are the 150 CollecTRI TFs with most targets in the data, identical for every method (pre-existing bench design) |
| cisTarget motif rankings + motif-to-TF annotation (Aerts lab, v10 clustered motifs) | genome sequence scored with motif collections; motif-TF links from motif databases, motif similarity and orthology | allowed: no curated TF-target edge enters |
| CellOracle promoter base GRN | motif scan of promoters | allowed (same kind as cisTarget), not used: redundant and heavy install |
| DoRothEA / CollecTRI / TRRUST / decoupler TF activity | curated TF-target edges | **circular**, excluded |

Residual overlap to keep in mind: some CollecTRI edges were curated from papers that used
motif evidence, so motif priors and CollecTRI are not fully independent evidence. That is
biology, not leakage: no edge of the evaluation half is read by any method.

## 3. Candidates chosen

1. **Motif prior + rank product** (`motif_prior`, `pearson_motif`, `grnboost2_motif`,
   `regdiffusion_motif`), in `grnbench.infer.cistarget`. Every TF x gene pair gets the best
   genome-wide rank of the gene among the TF's motifs (500 bp and 10 kb windows, the two
   standard databases); a co-expression ranking is combined with it by rank product
   (geometric mean of ranks, Breitling et al. 2004), no tuned parameter. `motif_prior`
   alone is the sequence-only control. numpy + pyarrow only; resources ~1.3 GB for human
   and mouse (four rankings files 237-313 MB, two annotation tables 99-113 MB) from
   https://resources.aertslab.org/cistarget/. Licence of the database files: not stated on
   the resource pages (public, widely redistributed with pySCENIC); verify before bundling.
2. **pySCENIC step 2 reimplemented** (`grnboost2_cistarget`, also `pearson_cistarget`):
   modules, NES >= 3 and leading edge as in pySCENIC defaults, in ~100 lines of numpy
   instead of the unmaintained `pyscenic`/`ctxcore` stack; output keeps the full ranking with
   motif-supported edges first so the evaluation universe does not change.
3. **RegDiffusion** (`regdiffusion`, pip `regdiffusion`, Apache-2.0, torch CPU): the most
   recent expression-only model that is pip-installable and CPU-feasible.
4. Kept as reference models: `pearson_signed`, `spearman_signed`, `clr` (v2 bench) and
   `grnboost2` (arboreto, already mounted; top-2 in geneRNIB).

Development check (PBMC3k only, 2,638 cells, 150 CollecTRI TFs x 1,650 genes, whole
CollecTRI as gold; PBMC3k is therefore reported apart from the test datasets):

| method | AUPRC ratio | EPR | seconds (loaded laptop) |
|---|---|---|---|
| pearson_signed | 1.494 | 3.461 | 1.5 |
| motif_prior (sequence only) | 1.609 | 2.145 | 6.0 |
| pearson_motif | 1.831 | 3.461 | 1.5 |
| pearson_cistarget | 1.424 | 2.521 | 128.8 |

Command: `GRNBENCH_CISTARGET_DIR=~\grnbench-data\cistarget PYTHONPATH=src
.venv/Scripts/python <scratchpad>/dev_methods.py pearson_signed,motif_prior,pearson_motif,pearson_cistarget`.
How the combination was picked: on PBMC3k I tried motif-only, the rank product, and
thresholded variants (motif rank < 500/1000/2000/5000 as a lexicographic tier or a x1.5
bonus). The rank product had the highest AUPRC ratio (1.833 vs 1.61-1.77) and has no
parameter, so it was adopted. Because PBMC3k informed that choice it is a development
dataset: it is reported separately and excluded from the test comparison.
