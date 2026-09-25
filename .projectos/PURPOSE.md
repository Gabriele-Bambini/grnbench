# Purpose - grnbench
> Owner: the user. Agents may draft, but changes need user confirmation and a DECISIONS.md entry.
> Direction set by the user 2026-09-23 (D-012, supersedes D-011).

## Why this project exists
Getting single-cell data ready for gene regulatory network (GRN) inference is messy, manual and different for every dataset and every method. grnbench is a **universal tool that cleans the data for GRN inference intelligently**, with a **usable interface**, and has the **GRN models already mounted on top** so the user can run them and **compare them** in one place.

## Goal (definition of done)
**A published tool** (user, 2026-09-23): released package (PyPI/conda) + software paper (target: Bioinformatics Application Note or JOSS - to confirm).
One tool where a user, without writing code:
1. loads any common scRNA-seq input (h5ad, 10x h5/MTX, loom, CSV/TSV);
2. gets an automatic diagnosis of the data and an intelligent, explained cleaning (QC, doublets, batch, gene IDs/organism, normalization, feature selection keeping TFs), with every choice shown and overridable;
3. runs the built-in GRN models (and plugins) on the cleaned data;
4. compares the models side by side (metrics vs gold standards, plots, report) and exports everything.

## Success criteria (measurable)
- [ ] Usable interface covering load -> diagnose/clean -> run models -> compare -> export
- [ ] Cleaning is data-driven (adaptive thresholds, auto-fixes where safe), explained, overridable, recorded in the manifest
- [ ] Works end-to-end on the BEELINE datasets and on a raw 10x dataset
- [ ] Models run from the interface; unavailable ones listed with the reason
- [ ] Comparison view + exportable report
- [ ] Publishable: pip/conda install, docs + tutorial, CI, case studies (BEELINE + raw 10x), paper draft
(to refine with the user)

## Constraints (time, compute, data, tools, style)
- Raw integer counts only; no imputation; TFs never dropped by HVG; QC gates fail loudly (D-001..D-004)
- Full provenance for every run (D-005); unavailable dependencies are reported, not crashes (D-006)
- Honest comparisons: statistics over independent units, never over edges/cells (D-008)
- Python, CPU; plugin-first for models (D-009)

## Non-goals
- A metascience/methods paper as the primary goal (D-011 rejected)
