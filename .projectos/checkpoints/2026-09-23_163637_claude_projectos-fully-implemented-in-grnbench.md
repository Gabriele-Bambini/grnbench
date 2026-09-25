---
time: 2026-09-23T16:36:37+02:00
agent: claude
kind: manual
reason: milestone
summary: ProjectOS fully implemented in grnbench: kit vendored, history backfilled (13 attempts, 11 decisions), PURPOSE draft + STATE verified
usage: ctx 388k tok
machine: <host>
---
# Checkpoint - ProjectOS fully implemented in grnbench: kit vendored, history backfilled (13 attempts, 11 decisions), PURPOSE draft...

**Advances purpose:** G0 software (memory of audit findings)

## Files edited since previous checkpoint (3)
- .projectos/ATTEMPTS.md (Edit)
- .projectos/PURPOSE.md (Write)
- .projectos/STATE.md (Write)

## STATE.md snapshot

#### State - grnbench
Updated: 2026-09-23 by claude (backfill from CHANGELOG 0.7.0 + audit 2026-09-13, verified against sources)

#### Focus
Post-0.7.0: pursue the audit's pilot plan (E0-E4 / gates G0-G5) toward a metascience result - pending user confirmation of direction (D-011 PROPOSED) and PURPOSE (draft).

#### Done recently
- 0.7.0 (2026-09-15): retracted the 0.6.0 edge-level significance layer (pseudo-replication, A-003); added `compare_units` / `grnbench compare` (paired t, Holm, Bonferroni intervals) and backtest `unit_comparison` + `regime_check`; added `eval.universe` fixed candidate universe (fixes A-009); gold thinning now fixed-k scenarios (fixes A-007); MI binning by first rank (fixes A-010); BEELINE TF selection pinned to official script SHA-256 (A-013).
- Real run: 9 methods x 7 BEELINE datasets, datasets as units, Holm: **0/108 EPR, 0/108 AUPRC-ratio resolved**; point-estimate winner changes with the reference on 7/7 datasets. Archived `benchmarks/results/beeline-2026-09-17-9methods/`.
- Score-level calibration: confirmation run (1000 reps) passes all iid scenarios; pseudoreplication fails as designed (A-011) -> G3 met at score level.
- Preprocessing fixes in 0.7.0 (0-cell output refused, strategy warnings captured in manifest, explicit 0 thresholds honored, label/group imbalance split, NaN label stratum, duplicate barcodes).

#### In progress
- Nothing running (no RUNNING attempts).

#### Next steps
1. Confirm with the user: PURPOSE draft + direction D-011 (metascience study vs software-only).
2. G0: turn audit fixtures (audit/2026-09-13/eval_findings.json, reproduce_mi_findings.json) into regression tests if not already; confirm the degree null issue (A-008) - 0.7.0 says `degree` mode unchanged, so it may still be open.
3. E1 / G1: reproduce faithfully one external published protocol that makes a concrete superiority claim, and document whether the claim survives unit-level comparison.
4. E2: controlled calibration on simulators (SERGIO / BoolODE are on the README roadmap, not integrated).
5. G4: compare against correctly applied Holm / max-T / csranks before claiming added value.

#### Blockers / open questions
- Direction and venue not confirmed by the user.
- Degree-preserving null (A-008) status after 0.7.0 unverified.
- Are cell types / organisms within a dataset independent units? Current design assumes whole datasets only (G2).

#### Notes for the next agent
- Read `.projectos/MAP.md` for code layout; the audit plan (Italian) is `audit/2026-09-13/AUDIT_E_PIANO_PUBBLICAZIONE.md` section 5.
- Never quote the 0.6.0 numbers (80% resolved, 47% contradictions): retracted.
- Unresolved != equivalent: the finding is that the design cannot rank, not that methods are equal.
- ProjectOS kit is vendored in `.projectos/kit/` (`python .projectos/kit/install.py` on a new machine).

