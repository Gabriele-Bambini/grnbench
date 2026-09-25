# Availability statements (draft)

## Code availability
grnbench is available at <public repository URL at submission> under the MIT licence (D-018 to
confirm). The exact version used for every result is the commit recorded in its experiment record
(`experiments/`); the archived release will be deposited on Zenodo with a DOI at submission.

## Data availability
All datasets are public. `experiments/datasets.tsv` lists every file used with its source URL or
accession, version and SHA-256: CELLxGENE Discover datasets (dataset IDs listed), 10x Genomics PBMC 3k
(Cell Ranger 1.1.0), BEELINE v1.1 (Zenodo 19009603), CollecTRI and DoRothEA via OmniPath/decoupler,
cisTarget v10 motif rankings (resources.aertslab.org). Results tables and per-dataset scores are in
`benchmarks/results/` and are regenerated into `paper/tables` by `python paper/make_figures.py`.

## Software versions
Python 3.11; package versions pinned in `requirements/lock-py311.txt`; optional models add
regdiffusion 0.2.3 and torch 2.14.0 (CPU). Statistics: `grnbench.eval.comparison.compare_units`
(paired t over datasets, Holm), in this repository.
