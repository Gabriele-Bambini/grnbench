---
time: 2026-09-23T20:35:32+02:00
agent: claude
kind: manual
reason: milestone
summary: G0: audit fixtures mapped to regressions (tests/test_audit_fixtures.py), A-008 closed; full suite pending
usage: ctx 133k tok
machine: <host>
---
# Checkpoint - G0: audit fixtures mapped to regressions (tests/test_audit_fixtures.py), A-008 closed; full suite pending

## Files edited since previous checkpoint (2)
- tests/test_audit_fixtures.py (Write)
- .projectos/ATTEMPTS.md (Edit)

## STATE.md snapshot

#### State - grnbench
Updated: 2026-09-23 by claude (G0 audit-fixture regressions; A-008 closed)

#### Focus
Post-0.7.0: pursue the audit's pilot plan (E0-E4 / gates G0-G5) toward a metascience result - pending user confirmation of direction (D-011 PROPOSED) and PURPOSE (draft).

#### Done recently
- G0 (2026-09-23): every audit/2026-09-13 eval_findings case now mapped to a regression on its exact input; new tests/test_audit_fixtures.py (coverage guard + null_support and finite_random_ap_baseline gaps), 3/3 pass. MI fixtures already exact in tests/test_inference_reliability.py. A-008 degree null verified fixed and closed.
- 0.7.0 (2026-09-15): retracted the 0.6.0 edge-level significance layer (pseudo-replication, A-003); added `compare_units` / `grnbench compare` (paired t, Holm, Bonferroni intervals) and backtest `unit_comparison` + `regime_check`; added `eval.universe` fixed candidate universe (fixes A-009); gold thinning now fixed-k scenarios (fixes A-007); MI binning by first rank (fixes A-010); BEELINE TF selection pinned to official script SHA-256 (A-013).
- Real run: 9 methods x 7 BEELINE datasets, datasets as units, Holm: **0/108 EPR, 0/108 AUPRC-ratio resolved**; point-estimate winner changes with the reference on 7/7 datasets. Archived `benchmarks/results/beeline-2026-09-17-9methods/`.
- Score-level calibration: confirmation run (1000 reps) passes all iid scenarios; pseudoreplication fails as designed (A-011) -> G3 met at score level.
- Preprocessing fixes in 0.7.0 (0-cell output refused, strategy warnings captured in manifest, explicit 0 thresholds honored, label/group imbalance split, NaN label stratum, duplicate barcodes).

#### In progress
- Full test suite run on this machine (Python 3.12, `PYTHONPATH=src`; package not pip-installed here) - result pending, confirm next session.

#### Next steps
1. Confirm with the user: PURPOSE draft + direction D-011 (metascience study vs software-only).
2. G0: confirm full suite green here (`PYTHONPATH=src python -m pytest -q -p no:cacheprovider`); fixtures part done.
3. E1 / G1: reproduce faithfully one external published protocol that makes a concrete superiority claim, and document whether the claim survives unit-level comparison.
4. E2: controlled calibration on simulators (SERGIO / BoolODE are on the README roadmap, not integrated).
5. G4: compare against correctly applied Holm / max-T / csranks before claiming added value.

#### Blockers / open questions
- Direction and venue not confirmed by the user.
- Are cell types / organisms within a dataset independent units? Current design assumes whole datasets only (G2).

#### Notes for the next agent
- This machine (collab-pc): grnbench not installed; run tests with `PYTHONPATH=src` (python 3.12 has numpy/pandas/scipy/anndata/scanpy). Suite is slow on Google Drive.
- Read `.projectos/MAP.md` for code layout; the audit plan (Italian) is `audit/2026-09-13/AUDIT_E_PIANO_PUBBLICAZIONE.md` section 5.
- Never quote the 0.6.0 numbers (80% resolved, 47% contradictions): retracted.
- Unresolved != equivalent: the finding is that the design cannot rank, not that methods are equal.
- ProjectOS kit is vendored in `.projectos/kit/` (`python .projectos/kit/install.py` on a new machine).

