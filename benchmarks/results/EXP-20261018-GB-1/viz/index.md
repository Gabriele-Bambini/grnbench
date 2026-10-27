# EXP-20261018-GB-1: figures

## Not produced

- Every figure of the "look inside" list (QC before/after, embeddings, pipeline decisions, score
  distributions and paired lifts, PR curves, top-k network vs a reference, TF-level errors, motif
  lift, calibration, stability, attributions). This run checks that the `aracne3` adapter drives the
  real ARACNe3 program and whether same-seed runs repeat: no network is scored against a reference,
  no cleaning is chosen and no model is compared, so those figures have nothing to show. The
  networks themselves are in `runs/<label>/edges.tsv` and `runs/<label>/aracne3_out/`.
