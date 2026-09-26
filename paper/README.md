# paper/

Everything that will appear in the manuscript, and nothing typed by hand that a script can make.

| file | what |
|---|---|
| `make_figures.py` | regenerates every table in `tables/` and figure in `figures/` from committed results: `python paper/make_figures.py` |
| `CLAIMS.md` | claims ledger: each sentence-level claim -> experiment IDs -> result files -> reproduced by |
| `REPORTING_CHECKLIST.md` | the statistics and reporting standard every claim must meet (Nature reporting summary / life-sciences checklist) |
| `CONTRIBUTIONS.md` | CRediT contribution log, per PR / experiment |
| `AVAILABILITY.md` | draft code and data availability statements, software versions |

## Rules

1. A claim enters the manuscript only through `CLAIMS.md`, with experiment IDs from `experiments/`.
2. Headline claims need a **confirmatory** experiment (pre-registered, see `experiments/README.md`).
   Exploratory results may be described as exploratory, never as the finding.
3. **Second-person reproduction**: the claim's status becomes `reproduced` only after someone other
   than the owner reran it from the record (section 4 of the record filled). Only `reproduced`
   claims go into a submitted draft. An agent that did not run the experiment can rerun it from the
   record first: the status is then `reproduced-computational`, which says a second computation
   agrees and a human reproduction is still pending. It never counts as `reproduced`.
4. Tables and figures: add a function to `make_figures.py`, commit the generated output in the
   same PR. A reviewer runs `python paper/make_figures.py` and checks `git diff --exit-code paper/`.
5. Every number in the text is copied from a generated table, with the table named next to it in
   `CLAIMS.md`.
6. Licence: code is MIT (`LICENSE`), compatible with journal code-sharing policies. Open decision
   D-018: confirm MIT and the copyright holder name before submission; data we redistribute keep
   their own licences (CELLxGENE datasets are CC BY 4.0; check each one in `experiments/datasets.tsv`).
