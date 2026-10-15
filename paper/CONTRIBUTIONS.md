# Contribution log (CRediT)

Add a row in the PR that does the work. Roles (CRediT): Conceptualization, Data curation, Formal
analysis, Funding acquisition, Investigation, Methodology, Project administration, Resources,
Software, Supervision, Validation, Visualization, Writing - original draft, Writing - review & editing.
Agents are tools, not contributors: log the person who directed the agent.

| date | person | roles | what (PR / EXP / issue) |
|---|---|---|---|
| 2026-09-25 | LM | Software, Investigation, Formal analysis | motif-prior / cisTarget / RegDiffusion models; EXP-20260925-LM-1 |
| 2026-09-25 | LM | Project administration | GitHub repo, experiment registry, paper/ standard |
| 2026-09-25 | LM | Investigation, Formal analysis, Data curation | EXP-20260925-LM-2 confirmatory run (pre-registration, 20 datasets) |
| 2026-10-01 | LD | Conceptualization, Methodology, Software, Validation, Investigation, Formal analysis, Data curation, Visualization | context-specific simulator, controlled worst-case stress design and complete resumable synthetic + real benchmark; issues #22/#23, EXP-20261001-LD-1 |
| 2026-10-03 | GB | Software | ingest count-plausibility warning, specific refusal advice and `benchmarks/count_plausibility.py`; issue #42 |
| 2026-10-05 | GB | Methodology, Software, Formal analysis | motif-prior lift decomposition (stage 1, Pearson): registration and `benchmarks/motif_sensitivity.py`; EXP-20261003-GB-1, issue #52 |
| 2026-10-06 | GB | Software | eleven literature GRN methods (bicor, WGCNA, propr, ARACNe-style, MRNET, TIGRESS, BBSR, pcNet, extra-trees, boosting, community), extra-method registry flag, `docs/models_catalog.md`; branch `gb/more-models` |
| 2026-10-07 | GB | Software, Visualization | interface: quick run, model picker, data and results explorers (`grnbench ui`); branch `gb/ui-explorer` |
| 2026-10-09 | GB | Software, Visualization | interface redesign: screens, menu, theme and plain-language summaries (`grnbench ui`); branch `gb/ui-explorer` |
| 2026-10-09 | GB | Software, Validation | end-to-end test of the integrated tool and fixes of the input-guard, CLI, ingest and viz defects it found; issue #74 |
| 2026-10-09 | LM | Visualization, Software | main figures for C4-C7 (`paper/figures/fig_motif_lift`, `fig_reference_free`, `paper/make_figures.py`); PR #62 |
