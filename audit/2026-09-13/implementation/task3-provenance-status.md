# Task 3 provenance and backtest changes

Implemented 2026-09-13. Source scope: src/grnbench/report.py, pipeline.py,
cli_infer.py, cli_run.py; benchmarks/beeline_backtest.py; hpc setup/fetch/backtest
scripts and README; tests/test_provenance_reliable.py plus one forced-QC test.

- Files always receive full streamed SHA256, including same-size tail changes.
  Named directories receive deterministic recursive file inventories; linked
  directories are refused before traversal. Ordered effective identifiers are
  retained and hashed. Environment includes all installed distribution versions.
- Prep records all parameters, effective genes/TFs and complete input identity.
  Infer writes full TF x gene candidate_universe.tsv.gz before method runs;
  run forwards --universe to eval. Empty candidate universes fail clearly.
- Inference catches failures per method and checkpoints diagnostic manifests,
  preserving successful edges and run_method provenance. Explicit unavailable
  methods fail; all unavailable retains code 2. --force bypasses QC but cannot
  make a zero-gene result valid for inference.
- Prep, infer, run and backtest require empty output directories. This deliberate
  stricter contract prevents mixing earlier runs; choose a new directory.
- Backtest persists selected expression, candidate U, effective gold, input
  hashes, ordered gene/TF identities, run_method provenance and edge lists.
  Missing datasets, references, methods and runtime errors return nonzero while
  retaining diagnostics. Removed empirical ceiling/pass and label-bootstrap
  significance claims. Random-ranking null is diagnostic on fixed U.
- HPC setup uses errexit and installs all+dev without a hidden fallback. A real
  Bash failure-path fixture proves pip and pytest errors cannot reach SETUP_DONE.

## Pinned primary evidence

Official code: https://github.com/Murali-group/Beeline/blob/37464085eb8a95d6cc6a3d3a3c649d36db6052ed/utils/generateExpInputs.py

Audit copy: beeline-reference-generateExpInputs.py. SHA256:
a50fbf3ff9e082f2f07f4f123f4d287ef6864be28ab695392aa54f701fc007c2

Selection: p < .01 / tested expression genes; significant TFs plus N significant
non-TFs by descending variance; alphabetical final genes. Exact fixture runs
the pinned reference and compares selected genes. Unambiguous case-insensitive
TF matching is an explicit extension. Evaluation remains grnbench's fixed
universe average precision, not all historical BEELINE metric conventions.

Data release: https://zenodo.org/records/19009603 is v4 for BEELINE v1.1.
fetch_beeline.sh checks published archive MD5s; sbatch records the release.
Custom backtest inputs default to unspecified-local-input; caller can declare
--data-release. Local input SHA256 remains the actual identity evidence.

## Verification and resume

Runtime: ~\.codex\tmp\grnbench-audit-20260913\Scripts\python.exe

Tests: python -m pytest tests/test_cli.py tests/test_provenance_reliable.py
tests/test_features_pipeline.py -q

- provenance-red.log: five expected initial regressions failed.
- provenance-green.log: initial five pass.
- provenance-replay.log: eight pass, including pinned reference and backtest.
- provenance-final.log: 21 pass with six AnnData index-conversion warnings.
- provenance-integrated.log: 37 pass with 12 AnnData index-conversion warnings.
- provenance-sklearn.log: latest 10 targeted tests pass, including independent
  AP/AUROC checks and both Bash failure paths.
- provenance-affected.log preserves the earlier two integration failures before
  correcting unavailable-only status and the scientifically unusable force case.

No large corpus downloaded. Linux cluster setup and full biological replay were
not executed. The official fixture validates selection; the persisted backtest
fixture checks AP and AUROC independently against sklearn on identical edges.
Synthetic fixtures do not establish biological validity or publication readiness.

Post-mortem: scope output freshness explicitly; silently reusing a directory
can pollute method discovery. Candidate U must be saved before inference outputs.
