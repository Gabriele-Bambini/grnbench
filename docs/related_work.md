# Related work: where single-cell GRN inference stands (2026-10-07)

A map of the field for choosing experiments and framing the paper. It complements
`docs/models_that_work.md` (which methods to run and why) and does not repeat its benchmark numbers.

Labels: **[abstract]** = checked against the abstract or a search summary only; **[unverified]** =
from background knowledge, DOI not re-fetched. Nothing here was read in full text. Check every citation
in the primary source before it enters `paper/`.

## 1. Research lines and their narrative arcs

| line | papers | arc |
|---|---|---|
| Aerts lab | SCENIC (Aibar 2017, 10.1038/nmeth.4463); pySCENIC protocol (Van de Sande 2020, 10.1038/s41596-020-0336-2); SCENIC+ (Bravo González-Blas 2023, 10.1038/s41592-023-01938-4); enhancer design (Taskiran 2024, 10.1038/s41586-023-06936-2); CREsted (Kempynck 2026, 10.1038/s41592-026-03057-2) [abstract] | Edges are noisy, motif-pruned regulons describe cell states; then the motif itself becomes the object, tested by building DNA. The cisTarget prior we use is theirs. |
| Saez-Rodriguez lab | DoRothEA (Garcia-Alonso 2019, 10.1101/gr.240663.118); decoupleR (Badia-i-Mompel 2022, 10.1093/bioadv/vbac016); CollecTRI (Müller-Dott 2023, 10.1093/nar/gkad841); GRETA (2024, 10.1101/2024.12.20.629764) [abstract] | Curated priors plus activity scoring beat inferred edges; simple linear scoring suffices; inferred multimodal GRNs are heterogeneous and modest. Our results agree, but we score against their CollecTRI. |
| Morris lab | CellOracle (Kamimoto 2023, 10.1038/s41586-022-05688-9) | The network as a simulator: judged by whether in-silico knockouts predict phenotypes (validated in zebrafish). |
| Chromatin and multiome | FigR (Kartha 2022), scGLUE (Cao and Gao 2022), Dictys (Wang 2023), Pando (Fleck 2023), LINGER (Yuan and Duren 2024, 10.1038/s41587-024-02182-7) [abstract] | Chromatin gives mechanism and direction. LINGER reports existing methods at 17-29% above random (our numbers: `docs/models_that_work.md`) and gains with external data plus a motif prior. |
| Foundation models | Geneformer (Theodoris 2023), scGPT (Cui 2024), scFoundation (Hao 2024), scPRINT (Kalfon 2025) | Claims of emergent GRNs, then refuted by baselines: Ahlmann-Eltze 2025 (10.1038/s41592-025-02772-6), Csendes 2025, Kedzierska 2025, Systema 2025, arXiv 2602.17532 (attention is co-expression) [abstract]. |
| Benchmarks | DREAM5 2012; BEELINE 2020; McCalla 2023 (10.1093/g3journal/jkad004); CausalBench; PEREGGRN (10.1101/2023.07.28.551039); geneRNIB (10.1101/2025.02.25.640181) | "Which method is best", then "nothing beats baselines and gold standards are the problem", then "score against perturbation response". |
| Causal and interventional | DCD-FG 2022; Replogle 2022 (10.1016/j.cell.2022.05.013); inspre 2025; TRADE (Nadig 2025) [abstract]; Ota 2025 [abstract] | Genome-scale knockdowns exist, but interventional causal discovery does not beat observational methods; knockdown effects spread far downstream (about 45 genes per typical knockdown, over 500 for essential genes). |

## 2. Open questions (2026)

1. Is the gain from a prior real or circular? Motif and ChIP evidence sit inside both priors and gold
   standards. Our C4-C6 raise it; EXP-20261003-GB-1 (#52, results in #92, exploratory) found a TF-generic gene propensity reaches about all of the Pearson lift; knockdown truth would settle it.
2. Direct or indirect: knockdown responses mix direct targets with cascades (TRADE 2025; Systema).
3. What is ground truth: curated sets favour well-studied TFs; binding is not regulation (McCalla 2023).
4. Is a TF's mRNA informative about its activity? We found no abstract that plots TF detectability against recovery.
5. Context: do networks transfer across cell types and disease states (Dictys; GRETA)?
6. Why interventional data does not help causal discovery (CausalBench; inspre).
7. Do foundation models carry regulatory signal beyond co-expression (Ahlmann-Eltze 2025)?
8. Self-assessment bias in benchmarks (Norel 2011, 10.1038/msb.2011.70; Weber 2019,
   10.1186/s13059-019-1738-8); we found pre-registration in none of the abstracts we searched.

## 3. Where our results fit

This section is a proposal pending D-022 (proposed on branch `lm/d022-direction`) and the team's decision in
issue #73; it is not a project claim. Our evidence matches the field's 2025-26 position (expression-only near
random; adding the cisTarget motif prior raises the AUPRC ratio on the CollecTRI evaluation half, claims C4-C6 in
`paper/CLAIMS.md`, with the C5 caveat that its datasets were already used for C4, resolved by C6 in
`experiments/EXP-20261002-LM-1.md`; honest evaluation). A contribution we did not find in the abstracts we searched: decompose the motif prior's gain against a
truth that shares no evidence with it (knockdowns), keep pre-registration, and, as an exploratory example of use run after the engine plan,
test whether it matters in patients (TP53 status in TCGA breast tumours, within PAM50 subtype). Paper template that worked for
SCENIC, LINGER and CollecTRI: benchmark exposes a failure, a simple fix, perturbation validation, disease
application.

## 4. Figures and explainability that earned trust

- Prediction next to intervention: CellOracle predicted vs observed knockout density; SCENIC+ F1 against
  ChIP-seq and perturbation data (PMC10482700, Fig. 3).
- Per-TF case studies with independent evidence tracks (SCENIC+ Fig. 3h-n).
- Rank heatmaps with positive and negative controls and a Pearson baseline (geneRNIB Fig. 3).
- Stability: GRNBoost2 is stochastic; the pySCENIC FAQ keeps targets seen in more than 80% of runs; SCING
  shows 50/50 split reproducibility.
- Avoid: whole-network node-link drawings (Gehlenborg and Wong, 10.1038/nmeth.1862); attention maps without a
  trivial baseline.
- Not found in the abstracts we searched (so possibly new if we do it): a TF detectability plot; calibrated per-edge evidence cards.
