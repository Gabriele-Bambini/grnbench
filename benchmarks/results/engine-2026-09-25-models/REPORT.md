# Modelli che funzionano - dati reali (2026-09-25, Lenovo, notte)

## Tabella principale (italiano)

Riferimento: metà di valutazione di CollecTRI (mai usata per scegliere la pipeline). Unità = dataset.
**n = 12 dataset di test** (CELLxGENE, 8 protocolli, uomo + topo). PBMC3k escluso: l'ho usato per
lo sviluppo. 1 = casuale. Pipeline "auto" (quella scelta dal motore).

| modello | AUPRC ratio (media) | EPR (media) | s/dataset | confronto dichiarato (AUPRC, Holm) |
|---|---|---|---|---|
| **regdiffusion_motif** (RegDiffusion x motivi) | **1.65** | **2.77** | ~810 | vs regdiffusion: +0.22, p = 9e-4 |
| **pearson_motif** (Pearson x motivi) | **1.55** | 2.33 | 3 | **vs Pearson: +0.23 [+0.16, +0.30], p = 1e-4** |
| motif_prior (solo sequenza, controllo) | 1.43 | 1.68 | 3 | - |
| regdiffusion | 1.43 | 2.72 | ~810 | vs Pearson: +0.11, non risolto (EPR +0.77, p = 0.003) |
| spearman_signed | 1.36 | 2.09 | 3 | - |
| pearson_signed (riferimento attuale) | 1.32 | 1.95 | 3 | - |
| clr | 1.30 | 1.91 | 19 | - |
| grnboost2_motif | 1.24 | 1.69 | 124 + 2 | vs grnboost2: +0.14, p = 1e-5 |
| grnboost2_cistarget (pySCENIC senza AUCell) | 1.17 | 1.96 | 124 + 61 | vs grnboost2: +0.08, p = 9e-4 |
| grnboost2 | 1.10 | 1.02 | 124 | peggio di Pearson: -0.22, p = 0.046 |

In breve:
- Il **prior di motivi (cisTarget, sequenza genomica)** alza il recupero per ogni modello a cui lo
  attacco: +0.23 su Pearson (12/12 dataset nella stessa direzione), +0.22 su RegDiffusion, +0.14 su
  GRNBoost2 (AUPRC ratio). Non usa CollecTRI: niente circolarità.
- Metà del guadagno viene dalla sequenza da sola (motif_prior 1.43 vs Pearson 1.32). Pearson x motivi
  vs solo motivi: +0.11 AUPRC non risolto, ma **+0.65 EPR (p = 0.01)**: i dati aggiungono precisione
  in cima alla lista.
- Resta modesto: il migliore arriva a 1.65 volte il caso (prima 1.2-1.3). In linea con la letteratura
  (BEELINE, geneRNIB: nessun metodo è molto sopra il caso sui dati reali).
- Pulizia auto vs baseline: ancora pari per tutti i modelli (differenze < 0.01).

## Detail (English)

### Setup

- Data: 12 CELLxGENE datasets (human 7, mouse 5; 10x 3' v2/v3, 10x 5' v1/v2, 10x multiome, CEL-seq,
  CEL-seq2, Smart-seq2, Patch-seq) fetched with
  `benchmarks/fetch_cellxgene.py --max 13 --max-mb 60 --organisms "Homo sapiens,Mus musculus"`,
  plus 10x PBMC3k filtered. A 13th CELLxGENE download (BD Rhapsody targeted mRNA, blood, 460 genes)
  was left out of `tasks.json` before the run: a 460-gene targeted panel is not a transcriptome.
  All data in `~\grnbench-data` (not on Drive).
- Reference: CollecTRI (decoupler 2.2.0), 150 TFs with most targets in the data, split into
  selection/evaluation halves by `screen.split_reference` (seed 0). Pipelines are chosen on the
  selection half; every number here is on the evaluation half.
- Motif resources: cisTarget v10 clustered rankings, hg38 and mm10, 500bp-up/100bp-down and
  10kbp-up/down, plus motif-to-TF tables (hgnc, mgi), from resources.aertslab.org/cistarget.
  No curated TF-target edge is used by any model.
- Arms: `baseline` (fixed standard pipeline) and `auto` (the engine's screened choice); every model is
  run on both. Contexts are averaged within a dataset; at most 3,000 cells per context.
- Code: frozen snapshot of commit 942cdea (`CODE_COMMIT`), run from
  `~\grnbench-runs\models-2026-09-25\code`, 3 workers, 03:28-05:18, 13 datasets, 0 errors:

  ```
  GRNBENCH_CISTARGET_DIR=~\grnbench-data\cistarget PYTHONPATH=code/src \
  python code/benchmarks/real_bench.py run --tasks tasks.json --refs ~/grnbench-data/references \
    --out out --jobs 3 --workroot work \
    --models pearson_signed,spearman_signed,clr,grnboost2,regdiffusion,motif_prior,pearson_motif,grnboost2_motif,regdiffusion_motif,grnboost2_cistarget
  python code/benchmarks/real_bench.py report --out <this folder>/test12   # 12 test datasets
  python code/benchmarks/real_bench.py report --out out                     # all 13 -> all13/
  ```

- Statistics: `compare_units` paired t over datasets (n = 12), Holm over the declared family (7 pairs x
  2 directions per metric and arm), declared in `real_bench.DECLARED_PAIRS` before the run. Gold edges
  are never units.

### Declared comparisons, auto pipeline, n = 12 (from `test12/REPORT_REAL.md`)

AUPRC ratio:

| A | B | mean A-B | 95% simultaneous interval | p Holm | supported |
|---|---|---|---|---|---|
| pearson_motif | pearson_signed | +0.229 | [+0.160, +0.297] | 9.9e-05 | pearson_motif |
| grnboost2_motif | grnboost2 | +0.143 | [+0.110, +0.177] | 9.6e-06 | grnboost2_motif |
| regdiffusion_motif | regdiffusion | +0.222 | [+0.136, +0.308] | 8.5e-04 | regdiffusion_motif |
| grnboost2_cistarget | grnboost2 | +0.076 | [+0.047, +0.106] | 8.5e-04 | grnboost2_cistarget |
| pearson_motif | motif_prior | +0.114 | [-0.029, +0.257] | 0.43 | unresolved |
| grnboost2 | pearson_signed | -0.221 | [-0.375, -0.067] | B>A: 0.046 | pearson_signed |
| regdiffusion | pearson_signed | +0.113 | [+0.025, +0.201] | 0.074 | unresolved |

EPR:

| A | B | mean A-B | 95% simultaneous interval | p Holm | supported |
|---|---|---|---|---|---|
| pearson_motif | pearson_signed | +0.372 | [+0.063, +0.681] | 0.10 | unresolved |
| grnboost2_motif | grnboost2 | +0.674 | [+0.511, +0.837] | 1.3e-05 | grnboost2_motif |
| regdiffusion_motif | regdiffusion | +0.043 | [-0.307, +0.392] | 1 | unresolved |
| grnboost2_cistarget | grnboost2 | +0.941 | [+0.668, +1.213] | 6.9e-05 | grnboost2_cistarget |
| pearson_motif | motif_prior | +0.648 | [+0.294, +1.002] | 0.0098 | pearson_motif |
| grnboost2 | pearson_signed | -0.934 | [-1.380, -0.489] | B>A: 0.0041 | pearson_signed |
| regdiffusion | pearson_signed | +0.772 | [+0.422, +1.122] | 0.0030 | regdiffusion |

The baseline pipeline gives the same verdicts (regdiffusion vs pearson_signed on AUPRC: +0.114,
p = 0.18, unresolved in both arms). Including PBMC3k (n = 13, `all13/REPORT_REAL.md`) changes one
verdict: regdiffusion vs pearson_signed AUPRC, auto arm, becomes resolved (p = 0.037). Quote n = 12.
"12/12 same direction" is read from the per-dataset table in `test12/REPORT_REAL.md` (auto arm).

### What this means for the product

1. Default model candidate: `pearson_motif` (3 s per dataset, numpy + pyarrow, +0.23 AUPRC ratio over
   Pearson). Best quality: `regdiffusion_motif` (1.65 / 2.77), ~13 min per dataset on this CPU (torch
   CPU, 1,000 steps). No GPU run was done.
2. The motif prior needs species data (human/mouse used here). Without `GRNBENCH_CISTARGET_DIR`, or
   for other species, the motif models are listed as unavailable with the reason.
3. pySCENIC-style pruning (`*_cistarget`) helps GRNBoost2 less than the rank product and costs ~60 s
   more per dataset; keep it as an option, not a default.

### What failed or is weak

- **GRNBoost2 is at random on EPR (1.02) and last on AUPRC (1.10)**, below Pearson on both
  (Holm-supported). This contradicts geneRNIB's ranking. Check on the dev dataset (PBMC3k, whole
  CollecTRI as gold): grnboost2 AUPRC ratio 1.227 / EPR 2.136 vs sklearn `random_forest` (an
  independent GENIE3-style implementation) 1.274 / 2.329, so the arboreto adapter does not look
  broken; tree importances are weaker than correlation on these data. Not investigated further.
- **RegDiffusion is slow on CPU**: ~810 s per dataset on average in the bench (the CPU was shared
  with other jobs). A separate probe (2,600 x 1,650, 4 torch threads, loaded machine): 30 steps in 120 s.
- **Auto vs baseline cleaning**: still ties on real data for every model (for example pearson_motif
  1.547 vs 1.538, auto > baseline on 2/12).
- The rank product was picked on PBMC3k (development) among 6 variants, so PBMC3k is excluded from
  the test comparison.
- Small scale: 12 test datasets, one laptop, human + mouse only, one reference (CollecTRI; DoRothEA
  was fetched but not scored), no GPU.
- Detection on these 12: organism 100%, technology 92% (Patch-seq returned None).

### Files

- `test12/REPORT_REAL.md`, `test12/all_rows_real.csv`, `test12/rows/*.json`: the 12 test datasets.
- `all13/`: same plus PBMC3k (development dataset).
- `tasks.json`, `CODE_COMMIT`: exact task list and code version.
- `grnbench-models.patch`: `git format-patch` of every code change since the Drive baseline.
