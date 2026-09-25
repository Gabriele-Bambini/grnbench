# grnbench preprocessing knowledge base (v2026.09.1)

Generated from `src/grnbench/kb.py` - edit the code, not this file.

## Evidence levels

- **grn**: benchmarked for GRN inference
- **coexpr**: benchmarked for gene-gene association / co-expression
- **scrna**: benchmarked for general scRNA-seq tasks (clustering, DE, integration)
- **consensus**: best-practice guideline / expert consensus
- **heuristic**: grnbench heuristic, validated on the grnbench bench

## Stages

1. **Load and verify** - Get a counts-backed matrix in the right orientation, with provenance.
   - pitfall: normalized/TPM matrices passed off as counts
   - pitfall: genes x cells tables read transposed
   - pitfall: duplicated barcodes or gene symbols
   - pitfall: processed data re-processed
2. **Gene identifiers and organism** - Make gene names match the TF list and the gold standards.
   - pitfall: Ensembl IDs vs symbols: TF lists silently match nothing
   - pitfall: version suffixes (ENSG...12)
   - pitfall: mouse vs human casing
   - pitfall: duplicated symbols after mapping
3. **Empty droplets and ambient RNA** - Keep real cells; remove the cell-free RNA soup every droplet contains.
   - pitfall: raw (unfiltered) droplet matrices loaded as if they were cells
   - pitfall: ambient marker genes appearing in every cell type
4. **Cell quality control** - Remove damaged, empty-like and outlier cells without removing biology.
   - pitfall: fixed thresholds that do not fit the tissue or protocol
   - pitfall: high-mito cell types (muscle, kidney, hepatocytes) removed wholesale
   - pitfall: one threshold across batches of different depth
5. **Doublets** - Remove droplets holding two cells.
   - pitfall: doublets form fake intermediate states and fake co-expression of two programs
   - pitfall: calling on data without cluster structure over-calls
6. **Gene filtering** - Drop genes that carry no usable signal.
   - pitfall: filtering away lowly expressed TFs
   - pitfall: mitochondrial/ribosomal genes becoming QC-driven hubs
7. **Cellular context** - Decide which cells one network is inferred from.
   - pitfall: pooled cell types: composition creates correlations that are not regulation
   - pitfall: too few cells per context for the chosen method
8. **Normalization and transformation** - Remove depth differences without inventing or destroying gene-gene association.
   - pitfall: scaling/z-scoring changes co-expression
   - pitfall: residual methods designed for UMI counts used on reads
9. **Batch and confounders** - Keep technical and nuisance structure from becoming edges.
   - pitfall: embedding-only integration used as if it corrected genes
   - pitfall: batch effects persist in covariance after mean correction
   - pitfall: regressing out biology correlated with the confounder
10. **Sparsity and dropout** - Get enough signal per gene pair without smoothing in fake correlation.
   - pitfall: imputation/smoothing inflating correlations
   - pitfall: impure metacells mixing states
11. **Features and regulators** - Choose candidate regulators and target genes (the candidate universe).
   - pitfall: HVG selection dropping TFs
   - pitfall: TF list for the wrong organism
   - pitfall: universe differs between methods
12. **Cell budget** - Keep inference tractable without biasing the sample.
   - pitfall: random subsampling losing rare states
13. **Pseudotime** - Order cells in time for methods that need it.
   - pitfall: pseudotime computed on a branching process as if linear
14. **Model-specific input** - Give each model the input it was designed for.
   - pitfall: one input for every model: tree methods, discretizing MI methods and count models want different forms

## 2. Gene identifiers and organism

### How are gene identifiers made to match TF lists and gold standards?

*Impact: high. Applies: always. Default: `keep`.*

TF lists and gold standards use symbols; a matrix of Ensembl IDs matches none of them and every regulator silently disappears.

- **Keep as loaded** (`keep`; evidence: consensus; grnbench) - Use identifiers as they are.
  - use when: identifiers are already symbols of the right organism
  - avoid when: Ensembl IDs or mixed identifiers
  - refs: [heumos2023](https://doi.org/10.1038/s41576-023-00586-w)
- **Map Ensembl IDs to symbols** (`ensembl_to_symbol`; evidence: consensus; planned) - Strip version suffixes, map with a gene annotation service, sum counts of duplicated symbols.
  - use when: Ensembl or mixed identifiers
  - pros: TF lists and gold standards match
  - cons: needs network access (mygene) or an annotation table; some IDs map to nothing or to several symbols
  - GRN: Summing counts of IDs that map to one symbol keeps counts integer.
  - refs: [heumos2023](https://doi.org/10.1038/s41576-023-00586-w)
- **Map to another organism's orthologs** (`orthologs`; evidence: consensus; planned) - Translate mouse symbols to human (or back) to use a TF list / gold of the other species.
  - use when: no TF list or gold standard exists for the data's organism
  - cons: one-to-many orthologs; regulatory edges are not always conserved

## 3. Empty droplets and ambient RNA

### Which barcodes are cells?

*Impact: high. Applies: raw droplet matrix. Default: `as_provided`.*

An unfiltered droplet matrix is mostly empty droplets; treating them as cells floods every statistic with ambient RNA.

- **Use the provided (filtered) matrix** (`as_provided`; evidence: consensus; grnbench) - Trust the cell calling already done (e.g. Cell Ranger filtered_feature_bc_matrix).
  - use when: the matrix is already filtered (typical GEO/CELLxGENE downloads)
  - refs: [amezquita2020](https://doi.org/10.1038/s41592-019-0654-x)
- **Knee / inflection point** (`knee`; evidence: consensus; grnbench) - Rank barcodes by total counts; keep those above the knee of the log-log curve.
  - use when: raw droplet matrix and no R available
  - pros: fast, dependency-free
  - cons: misses small or low-RNA cells
  - refs: [emptydrops2019](https://doi.org/10.1186/s13059-019-1662-y)
- **EmptyDrops** (`emptydrops`; evidence: scrna; external) - Test each barcode's profile against the ambient profile.
  - use when: raw droplet matrix
  - pros: keeps low-RNA cells that differ from ambient
  - cons: needs R (DropletUtils)
  - refs: [emptydrops2019](https://doi.org/10.1186/s13059-019-1662-y)

### Is cell-free (ambient) RNA removed?

*Impact: medium. Applies: always. Default: `none`.*

Ambient RNA adds the same background profile to every cell: markers of abundant cell types appear everywhere, blurring cell-type-specific expression.

- **No correction** (`none`; evidence: scrna; grnbench) - Leave counts as they are.
  - use when: plate-based protocols; no raw matrix available; low contamination
  - pros: no risk of over-correction
  - cons: background inflates low-level expression everywhere
  - refs: [ambient2026](https://doi.org/10.3389/fbinf.2026.1844838)
- **SoupX** (`soupx`; evidence: scrna; external) - Estimate the soup from empty droplets and subtract it per cell.
  - use when: droplet data with the raw matrix; CPU only
  - pros: least over-correction when contamination is low
  - cons: needs raw matrix and clusters; R
  - GRN: No GRN-specific evidence; expected to reduce spurious marker co-occurrence.
  - refs: [soupx2020](https://doi.org/10.1093/gigascience/giaa151); [ambient2026](https://doi.org/10.3389/fbinf.2026.1844838)
- **DecontX** (`decontx`; evidence: scrna; external) - Bayesian mixture of native and contaminating counts.
  - use when: droplet data; works without the raw matrix
  - cons: can over-correct when contamination is absent
  - refs: [decontx2020](https://doi.org/10.1186/s13059-020-1950-6); [ambient2026](https://doi.org/10.3389/fbinf.2026.1844838)
- **CellBender** (`cellbender`; evidence: scrna; external) - Deep generative model of ambient RNA and barcode swapping on the raw matrix.
  - use when: droplet data, raw matrix and a GPU available
  - pros: most precise background estimates (Janssen 2023)
  - cons: needs GPU and long runtime; can over-correct when contamination is absent
  - refs: [cellbender2023](https://doi.org/10.1038/s41592-023-01943-7); [janssen2023](https://doi.org/10.1186/s13059-023-02978-x); [ambient2026](https://doi.org/10.3389/fbinf.2026.1844838)

## 4. Cell quality control

### How are low-quality cells defined?

*Impact: high. Applies: always. Default: `mad`.*

Thresholds decide which cells exist for inference: too loose keeps damaged cells whose mitochondrial/stress program becomes a hub; too strict removes whole cell types.

- **Fixed thresholds** (`fixed`; evidence: consensus; grnbench) - Classic cut-offs (e.g. >=200 genes, <=20% mitochondrial).
  - use when: very small datasets where distributions cannot be estimated
  - cons: ignores depth, protocol and tissue; high-mito tissues lose whole cell types
  - refs: [luecken2019](https://doi.org/10.15252/msb.20188746); [osorio2021](https://doi.org/10.1093/bioinformatics/btaa751)
- **MAD outliers (global)** (`mad`; evidence: scrna; grnbench) - Outliers at k median absolute deviations on log scale (genes, counts) and on % mito.
  - use when: one batch or batches of similar depth
  - pros: adapts to each dataset; standard in OSCA / scater / sc-best-practices
  - cons: assumes most cells are good
  - refs: [scater2017](https://doi.org/10.1093/bioinformatics/btw777); [pipecomp2020](https://doi.org/10.1186/s13059-020-02136-7); [heumos2023](https://doi.org/10.1038/s41576-023-00586-w)
- **MAD outliers per batch** (`mad_per_batch`; evidence: scrna; grnbench) - Same rule computed inside each batch/sample.
  - use when: several batches or samples with different depth
  - pros: a shallow batch is not removed wholesale
  - cons: needs >= ~100 cells per batch
  - refs: [amezquita2020](https://doi.org/10.1038/s41592-019-0654-x); [heumos2023](https://doi.org/10.1038/s41576-023-00586-w)
- **miQC mixture model** (`miqc`; evidence: scrna; planned) - Joint model of % mito and gene count separates damaged from intact cells.
  - use when: tissues with naturally high mitochondrial content
  - pros: probabilistic, tissue-adaptive
  - cons: assumes a two-population structure
  - refs: [miqc2021](https://doi.org/10.1371/journal.pcbi.1009290)

## 5. Doublets

### Are doublets detected and removed?

*Impact: medium. Applies: always. Default: `none`.*

A doublet expresses two programs at once, creating co-expression between genes of different cell types - exactly the signal a GRN method reads as regulation.

- **No removal** (`none`; evidence: consensus; grnbench) - Keep every droplet.
  - use when: plate-based protocols; already doublet-filtered data; few cells
  - refs: [heumos2023](https://doi.org/10.1038/s41576-023-00586-w)
- **Scrublet** (`scrublet`; evidence: scrna; grnbench) - Simulate doublets from observed cells; kNN score; run per sample.
  - use when: droplet data with cluster structure
  - pros: Python, fast
  - cons: threshold fails without cluster structure (grnbench refuses >20% calls)
  - refs: [scrublet2019](https://doi.org/10.1016/j.cels.2018.11.005); [xi2021](https://doi.org/10.1016/j.cels.2020.11.008)
- **scDblFinder** (`scdblfinder`; evidence: scrna; external) - Gradient-boosted classifier on simulated doublets.
  - use when: droplet data; R available
  - pros: top accuracy in its benchmark; fewest false positives in image-verified data
  - cons: R; misses part of the doublets on some data
  - refs: [scdblfinder2021](https://doi.org/10.12688/f1000research.73600.2); [imagedoubler2024](https://doi.org/10.1038/s41467-024-55434-0)
- **DoubletFinder** (`doubletfinder`; evidence: scrna; external) - Artificial nearest-neighbour doublet proportion.
  - pros: best accuracy in Xi & Li 2021
  - cons: R/Seurat; parameter sweep needed
  - refs: [doubletfinder2019](https://doi.org/10.1016/j.cels.2019.03.003); [xi2021](https://doi.org/10.1016/j.cels.2020.11.008)
- **Solo** (`solo`; evidence: scrna; external) - Semi-supervised deep classifier on scVI latent space.
  - cons: GPU-friendly deep model; more false positives on image-verified data
  - refs: [solo2020](https://doi.org/10.1016/j.cels.2020.05.010); [imagedoubler2024](https://doi.org/10.1038/s41467-024-55434-0)

## 6. Gene filtering

### Which genes are removed before inference?

*Impact: low. Applies: always. Default: `min_cells`.*

Genes seen in a handful of cells give unstable correlations; mitochondrial and ribosomal genes track cell quality and depth, not regulation.

- **Minimum cells expressing** (`min_cells`; evidence: consensus; grnbench) - Drop genes detected in fewer than max(3, 0.1% of cells).
  - refs: [luecken2019](https://doi.org/10.15252/msb.20188746)
- **Minimum cells + drop MT-/ribosomal genes as targets** (`min_cells_drop_qc_genes`; evidence: heuristic; grnbench) - Also remove mitochondrial and ribosomal protein genes from the features.
  - pros: prevents QC-driven hubs
  - cons: loses genuine regulation of those genes (rarely the goal)
  - GRN: MT-/RP genes co-vary with cell quality and depth across all cells.

## 7. Cellular context

### Which cells does one network come from?

*Impact: high. Applies: always. Default: `pooled`.*

In mixed tissue most gene-gene co-expression reflects cell-type composition, not regulation inside a cell type; cell-type mixing is one of the stressors that break GRN inference.

- **All cells pooled** (`pooled`; evidence: coexpr; grnbench) - One network from every cell.
  - use when: a single cell type or one continuous trajectory; few cells per type
  - cons: composition-driven correlations become edges
  - refs: [farahbod2020](https://doi.org/10.1101/gr.256735.119); [retana2026](https://arxiv.org/abs/2605.04930)
- **One network per cell type** (`per_cell_type`; evidence: coexpr; grnbench) - Split by the cell-type column and infer each context separately.
  - use when: annotated, distinct cell types with enough cells each (>= ~200)
  - pros: edges reflect regulation within a state
  - cons: fewer cells per network; needs annotation
  - refs: [farahbod2020](https://doi.org/10.1101/gr.256735.119); [retana2026](https://arxiv.org/abs/2605.04930)
- **One selected cell type** (`selected_type`; evidence: coexpr; grnbench) - Restrict the analysis to the cell type of interest.
  - use when: the question is about one population
  - refs: [farahbod2020](https://doi.org/10.1101/gr.256735.119)

## 8. Normalization and transformation

### How are depth differences removed?

*Impact: high. Applies: always. Default: `cpm_log1p`.*

Normalization changes gene-gene association directly; scaling can lower co-expression recovery and residual methods assume UMI count noise.

- **Shifted log (CPM + log1p)** (`cpm_log1p`; evidence: scrna; grnbench) - Scale to 10k counts per cell, log(x+1).
  - use when: default for UMI and full-length data
  - pros: as good as more complex transforms in a large comparison; the input BEELINE used for all methods
  - cons: depth effects can remain in lowly expressed genes; any scaling step can lower co-expression recovery
  - GRN: Keep it unscaled (no z-scoring) before correlation-based methods.
  - refs: [ahlmann2023](https://doi.org/10.1038/s41592-023-01814-1); [beeline2020](https://doi.org/10.1038/s41592-019-0690-6); [sccobench2025](https://doi.org/10.1101/2025.05.26.656221); [skinnider2019](https://doi.org/10.1038/s41592-019-0372-4)
- **Analytic Pearson residuals** (`pearson_residuals`; evidence: scrna; grnbench) - Residuals of a Poisson/NB null with fixed overdispersion.
  - use when: UMI data with very variable depth
  - pros: variance-stable across depth
  - cons: designed for UMI counts, not read counts; negative values
  - refs: [lause2021](https://doi.org/10.1186/s13059-021-02451-7); [ahlmann2023](https://doi.org/10.1038/s41592-023-01814-1)
- **sctransform (regularized NB)** (`sctransform`; evidence: scrna; external) - Regularized negative-binomial regression on depth.
  - cons: R; weaker on data with heavy dropout
  - refs: [hafemeister2019](https://doi.org/10.1186/s13059-019-1874-1)
- **scran pooled size factors + log** (`scran`; evidence: scrna; external) - Size factors pooled across similar cells, then log.
  - pros: robust size factors with many zeros
  - cons: R
  - refs: [lun2016](https://doi.org/10.1186/s13059-016-0947-7); [pipecomp2020](https://doi.org/10.1186/s13059-020-02136-7)
- **Raw counts** (`raw_counts`; evidence: consensus; grnbench) - No transformation, for models that model counts themselves.
  - use when: the downstream model declares count input
  - refs: [scenic2017](https://doi.org/10.1038/nmeth.4463)

## 9. Batch and confounders

### How are batches handled?

*Impact: high. Applies: batch column with >= 2 levels. Default: `none`.*

Batch effects create shared shifts in many genes - perfect fake co-expression - and they persist in the covariance even after mean correction.

- **No correction** (`none`; evidence: consensus; grnbench) - Leave batch structure in.
  - use when: a single batch, or batch fully confounded with biology (correction would remove it)
  - refs: [heumos2023](https://doi.org/10.1038/s41576-023-00586-w)
- **Center genes per batch (regress out batch)** (`center_per_batch`; evidence: coexpr; grnbench) - Subtract each gene's batch mean on the log scale.
  - use when: batches share cell-type composition
  - pros: removes location shifts; transparent
  - cons: does not remove variance/covariance differences
  - GRN: Removes the between-batch component of correlation.
  - refs: [cobra2024](https://doi.org/10.1093/bioinformatics/btae531)
- **ComBat** (`combat`; evidence: coexpr; grnbench) - Empirical-Bayes location/scale adjustment per gene.
  - use when: batches share cell-type composition
  - cons: covariance-level batch effects persist; removes biology aligned with batch
  - refs: [combat2007](https://doi.org/10.1093/biostatistics/kxj037); [cobra2024](https://doi.org/10.1093/bioinformatics/btae531); [tran2020](https://doi.org/10.1186/s13059-019-1850-9)
- **Infer per batch, combine edges** (`per_batch_consensus`; evidence: heuristic; planned) - Run inference inside each batch and aggregate the edge ranks.
  - use when: several batches with enough cells each
  - pros: no correction assumptions; batch-reproducible edges rank higher
  - cons: fewer cells per run; more compute
- **scVI denoised expression** (`scvi_denoised`; evidence: scrna; external) - Batch-conditioned generative model; decoded expression at a reference batch.
  - cons: model-smoothed expression can inflate correlations like imputation; GPU recommended
  - refs: [scvi2018](https://doi.org/10.1038/s41592-018-0229-2); [scib2022](https://doi.org/10.1038/s41592-021-01336-8)
- **Harmony / embedding-only integration** (`harmony`; evidence: scrna; not-applicable) - Corrects a low-dimensional embedding, not gene expression.
  - GRN: Gives no corrected gene matrix, so it cannot feed gene-level GRN inference.
  - refs: [harmony2019](https://doi.org/10.1038/s41592-019-0619-0); [scib2022](https://doi.org/10.1038/s41592-021-01336-8)

### Are nuisance signals regressed out?

*Impact: medium. Applies: always. Default: `none`.*

Regression removes correlation with the confounder - including real regulation of cell-cycle or stress genes.

- **No regression** (`none`; evidence: consensus; grnbench) - Keep cell cycle, depth and stress signal.
  - use when: default
  - pros: does not remove real regulation (cell-cycle TFs are real regulators)
  - refs: [luecken2019](https://doi.org/10.15252/msb.20188746); [heumos2023](https://doi.org/10.1038/s41576-023-00586-w)
- **Regress cell-cycle scores** (`cell_cycle`; evidence: consensus; planned) - Score S/G2M (Tirosh genes) and regress out.
  - use when: cycling is a nuisance for the question asked
  - cons: removes cell-cycle regulation from the network
  - refs: [tirosh2016](https://doi.org/10.1126/science.aad0501); [luecken2019](https://doi.org/10.15252/msb.20188746)

## 10. Sparsity and dropout

### How is sparsity handled?

*Impact: high. Applies: always. Default: `single_cells`.*

Dropout and latent confounders erase the advantage of better methods; imputation does not help network inference and can create false correlations.

- **Single cells** (`single_cells`; evidence: grn; grnbench) - Infer from single-cell profiles as they are.
  - use when: moderate sparsity and enough cells
  - GRN: Dropout robustness belongs in the model (e.g. dropout augmentation), not in imputed data.
  - refs: [beeline2020](https://doi.org/10.1038/s41592-019-0690-6); [mccalla2023](https://doi.org/10.1093/g3journal/jkad004); [dazzle2025](https://doi.org/10.1371/journal.pcbi.1013603)
- **Metacells** (`metacells`; evidence: coexpr; grnbench) - Aggregate counts of small groups of transcriptionally similar cells (kNN partitions).
  - use when: very sparse data with many cells (>= ~2000)
  - pros: less dropout; better gene-gene correlation
  - cons: fewer samples; impure metacells mix states
  - GRN: Check purity against the cell-type column; keep >= ~100 metacells.
  - refs: [bilous2024](https://doi.org/10.1038/s44320-024-00045-6); [seacells2023](https://doi.org/10.1038/s41587-023-01716-9); [mcrigor2025](https://doi.org/10.1038/s41467-025-63626-5)
- **Pseudobulk per group** (`pseudobulk`; evidence: consensus; grnbench) - Sum counts per annotated group (sample x cell type).
  - use when: many groups (tens) available
  - cons: needs many groups to have enough samples
  - refs: [heumos2023](https://doi.org/10.1038/s41576-023-00586-w)
- **Imputation / smoothing (MAGIC, kNN)** (`imputation`; evidence: grn; external) - Share information between neighbouring cells to fill zeros.
  - avoid when: network inference in general
  - cons: did not improve network inference and could be detrimental; creates false correlations
  - refs: [mccalla2023](https://doi.org/10.1093/g3journal/jkad004); [andrews2018](https://doi.org/10.12688/f1000research.16613.2); [hou2020](https://doi.org/10.1186/s13059-020-02132-x); [magic2018](https://doi.org/10.1016/j.cell.2018.05.061)

## 11. Features and regulators

### Which target genes enter the network?

*Impact: medium. Applies: always. Default: `tf_plus_1000`.*

Defines the candidate universe; runtime of most methods grows with genes x TFs.

- **TFs + top 1000 HVGs** (`tf_plus_1000`; evidence: grn; grnbench) - Highly variable genes union every TF detected in the data.
  - pros: BEELINE protocol; tractable
  - refs: [beeline2020](https://doi.org/10.1038/s41592-019-0690-6); [yip2019](https://doi.org/10.1093/bib/bby011); [stuart2019](https://doi.org/10.1016/j.cell.2019.05.031)
- **TFs + top 500 HVGs** (`tf_plus_500`; evidence: grn; grnbench) - Smaller BEELINE setting.
  - use when: slow methods or few cells
  - refs: [beeline2020](https://doi.org/10.1038/s41592-019-0690-6)
- **TFs + top 2000 HVGs** (`tf_plus_2000`; evidence: heuristic; grnbench) - Broader target space.
  - use when: fast methods and >= 2000 cells
  - refs: [beeline2020](https://doi.org/10.1038/s41592-019-0690-6)
- **All expressed genes** (`all_genes`; evidence: consensus; grnbench) - No selection.
  - use when: small gene space (< ~1500 genes)
  - cons: runtime and noise grow with genes

### Where does the list of candidate regulators come from?

*Impact: high. Applies: always. Default: `user_file`.*

Only listed TFs can be sources; a list for the wrong organism or naming matches nothing.

- **User TF list** (`user_file`; evidence: consensus; grnbench) - A curated list supplied with the data (e.g. the human TF census).
  - refs: [lambert2018](https://doi.org/10.1016/j.cell.2018.01.029)
- **DoRothEA TFs** (`dorothea`; evidence: consensus; grnbench) - TFs with curated regulons (confidence levels A-E).
  - pros: regulons double as a silver standard
  - cons: ~ few hundred TFs; circular if also used as gold
  - refs: [dorothea2019](https://doi.org/10.1101/gr.240663.118)
- **CollecTRI TFs** (`collectri`; evidence: consensus; planned) - Larger curated TF-target collection.
  - cons: circular if also used as gold
  - refs: [collectri2023](https://doi.org/10.1093/nar/gkad841)
- **Every gene** (`all_genes`; evidence: consensus; grnbench) - No regulator list.
  - cons: universe grows quadratically; edges not restricted to TFs

## 12. Cell budget

### Are cells subsampled?

*Impact: low. Applies: always. Default: `none`.*

Tree ensembles scale with cells x genes x TFs; subsampling must not drop rare states.

- **All cells** (`none`; evidence: consensus; grnbench) - Use every cell.
- **Stratified subsample** (`stratified`; evidence: heuristic; grnbench) - Subsample keeping each cell type's share.
  - use when: > 20k cells with annotation
- **Random subsample** (`random`; evidence: heuristic; grnbench) - Uniform subsample.
  - cons: can lose rare states

## 13. Pseudotime

### Is pseudotime computed?

*Impact: medium. Applies: dynamic methods selected or time-course data. Default: `none`.*

Dynamic methods (SCODE, SINCERITIES, SINGE, ...) require cells ordered in time.

- **No pseudotime** (`none`; evidence: consensus; grnbench) - Static methods only.
- **Use provided pseudotime** (`provided`; evidence: consensus; grnbench) - A column already in the metadata.
  - refs: [beeline2020](https://doi.org/10.1038/s41592-019-0690-6)
- **Diffusion pseudotime** (`dpt`; evidence: scrna; planned) - Scanpy DPT from a root cell.
  - cons: needs a root; assumes connected trajectory
  - refs: [dpt2016](https://doi.org/10.1038/nmeth.3971)
- **Slingshot** (`slingshot`; evidence: grn; external) - Lineage-aware pseudotime (used by BEELINE).
  - cons: R
  - refs: [slingshot2018](https://doi.org/10.1186/s12864-018-4772-0); [beeline2020](https://doi.org/10.1038/s41592-019-0690-6)

## References

- **luecken2019**: Luecken & Theis 2019, Mol Syst Biol - current best practices in scRNA-seq analysis. https://doi.org/10.15252/msb.20188746
- **heumos2023**: Heumos et al. 2023, Nat Rev Genet - best practices for single-cell analysis across modalities. https://doi.org/10.1038/s41576-023-00586-w
- **amezquita2020**: Amezquita et al. 2020, Nat Methods - Orchestrating single-cell analysis with Bioconductor (OSCA). https://doi.org/10.1038/s41592-019-0654-x
- **beeline2020**: Pratapa et al. 2020, Nat Methods - BEELINE: benchmarking GRN inference from scRNA-seq. https://doi.org/10.1038/s41592-019-0690-6
- **mccalla2023**: McCalla et al. 2023, G3 - strengths and weaknesses of network inference from scRNA-seq (imputation does not help). https://doi.org/10.1093/g3journal/jkad004
- **retana2026**: Fernandez-de-Retana et al. 2026, arXiv - when does GRN inference break (dropout, confounders, cell-type mixing). https://arxiv.org/abs/2605.04930
- **dazzle2025**: DAZZLE 2025, PLoS Comput Biol - GRN inference with dropout augmentation. https://doi.org/10.1371/journal.pcbi.1013603
- **farahbod2020**: Farahbod & Pavlidis 2020, Genome Res - cellular composition drives coexpression. https://doi.org/10.1101/gr.256735.119
- **skinnider2019**: Skinnider et al. 2019, Nat Methods - measures of association for single-cell transcriptomics. https://doi.org/10.1038/s41592-019-0372-4
- **sccobench2025**: scCoBench 2025, bioRxiv - co-expression benchmark with promoter-reporter lines. https://doi.org/10.1101/2025.05.26.656221
- **ahlmann2023**: Ahlmann-Eltze & Huber 2023, Nat Methods - comparison of transformations for scRNA-seq. https://doi.org/10.1038/s41592-023-01814-1
- **lause2021**: Lause et al. 2021, Genome Biol - analytic Pearson residuals. https://doi.org/10.1186/s13059-021-02451-7
- **hafemeister2019**: Hafemeister & Satija 2019, Genome Biol - sctransform. https://doi.org/10.1186/s13059-019-1874-1
- **lun2016**: Lun et al. 2016, Genome Biol - scran pooled size factors. https://doi.org/10.1186/s13059-016-0947-7
- **emptydrops2019**: Lun et al. 2019, Genome Biol - EmptyDrops. https://doi.org/10.1186/s13059-019-1662-y
- **janssen2023**: Janssen et al. 2023, Genome Biol - background noise and its removal. https://doi.org/10.1186/s13059-023-02978-x
- **ambient2026**: Benchmarking computational decontamination of ambient RNA, 2026, Front Bioinform. https://doi.org/10.3389/fbinf.2026.1844838
- **soupx2020**: Young & Behjati 2020, GigaScience - SoupX. https://doi.org/10.1093/gigascience/giaa151
- **cellbender2023**: Fleming et al. 2023, Nat Methods - CellBender. https://doi.org/10.1038/s41592-023-01943-7
- **decontx2020**: Yang et al. 2020, Genome Biol - DecontX. https://doi.org/10.1186/s13059-020-1950-6
- **scater2017**: McCarthy et al. 2017, Bioinformatics - scater (MAD outlier QC). https://doi.org/10.1093/bioinformatics/btw777
- **osorio2021**: Osorio & Cai 2021, Bioinformatics - mitochondrial thresholds by species/tissue. https://doi.org/10.1093/bioinformatics/btaa751
- **miqc2021**: Hippen et al. 2021, PLoS Comput Biol - miQC. https://doi.org/10.1371/journal.pcbi.1009290
- **pipecomp2020**: Germain et al. 2020, Genome Biol - pipeComp: filtering and normalization choices. https://doi.org/10.1186/s13059-020-02136-7
- **xi2021**: Xi & Li 2021, Cell Syst - benchmarking doublet detection. https://doi.org/10.1016/j.cels.2020.11.008
- **scdblfinder2021**: Germain et al. 2021, F1000Research - scDblFinder. https://doi.org/10.12688/f1000research.73600.2
- **scrublet2019**: Wolock et al. 2019, Cell Syst - Scrublet. https://doi.org/10.1016/j.cels.2018.11.005
- **doubletfinder2019**: McGinnis et al. 2019, Cell Syst - DoubletFinder. https://doi.org/10.1016/j.cels.2019.03.003
- **solo2020**: Bernstein et al. 2020, Cell Syst - Solo. https://doi.org/10.1016/j.cels.2020.05.010
- **imagedoubler2024**: ImageDoubler 2024, Nat Commun - image-based doublet ground truth. https://doi.org/10.1038/s41467-024-55434-0
- **seacells2023**: Persad et al. 2023, Nat Biotechnol - SEACells. https://doi.org/10.1038/s41587-023-01716-9
- **bilous2024**: Bilous et al. 2024, Mol Syst Biol - building and analyzing metacells. https://doi.org/10.1038/s44320-024-00045-6
- **mcrigor2025**: mcRigor 2025, Nat Commun - rigor of metacell partitioning. https://doi.org/10.1038/s41467-025-63626-5
- **andrews2018**: Andrews & Hemberg 2018, F1000Research - false signals induced by imputation. https://doi.org/10.12688/f1000research.16613.2
- **hou2020**: Hou et al. 2020, Genome Biol - imputation benchmark. https://doi.org/10.1186/s13059-020-02132-x
- **magic2018**: van Dijk et al. 2018, Cell - MAGIC. https://doi.org/10.1016/j.cell.2018.05.061
- **cobra2024**: COBRA 2024, Bioinformatics - batch effects persist in correlation networks. https://doi.org/10.1093/bioinformatics/btae531
- **scib2022**: Luecken et al. 2022, Nat Methods - benchmarking atlas-level integration (scIB). https://doi.org/10.1038/s41592-021-01336-8
- **tran2020**: Tran et al. 2020, Genome Biol - benchmark of batch-effect correction. https://doi.org/10.1186/s13059-019-1850-9
- **combat2007**: Johnson et al. 2007, Biostatistics - ComBat. https://doi.org/10.1093/biostatistics/kxj037
- **harmony2019**: Korsunsky et al. 2019, Nat Methods - Harmony. https://doi.org/10.1038/s41592-019-0619-0
- **scvi2018**: Lopez et al. 2018, Nat Methods - scVI. https://doi.org/10.1038/s41592-018-0229-2
- **yip2019**: Yip et al. 2019, Brief Bioinform - evaluation of HVG selection. https://doi.org/10.1093/bib/bby011
- **stuart2019**: Stuart et al. 2019, Cell - Seurat v3 (VST HVG). https://doi.org/10.1016/j.cell.2019.05.031
- **lambert2018**: Lambert et al. 2018, Cell - the human transcription factors. https://doi.org/10.1016/j.cell.2018.01.029
- **dorothea2019**: Garcia-Alonso et al. 2019, Genome Res - DoRothEA. https://doi.org/10.1101/gr.240663.118
- **collectri2023**: Mueller-Dott et al. 2023, Nucleic Acids Res - CollecTRI. https://doi.org/10.1093/nar/gkad841
- **dpt2016**: Haghverdi et al. 2016, Nat Methods - diffusion pseudotime. https://doi.org/10.1038/nmeth.3971
- **slingshot2018**: Street et al. 2018, BMC Genomics - Slingshot. https://doi.org/10.1186/s12864-018-4772-0
- **tirosh2016**: Tirosh et al. 2016, Science - cell-cycle scoring. https://doi.org/10.1126/science.aad0501
- **scenic2017**: Aibar et al. 2017, Nat Methods - SCENIC. https://doi.org/10.1038/nmeth.4463
