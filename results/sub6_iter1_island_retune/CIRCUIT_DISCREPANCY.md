# HARD_STOP_OFF_CIRCUIT_HFSS_MISMATCH

RESOLVED_BY_DECISION_SESSION: series-L circuit authoritative. Corrected
HFSS_SETUP.md and the research-only circuit test expectation; previous
regression tests/fixtures unchanged. Four new checks now pass. This file
preserves the original finding and arithmetic, not a current stop.

MODEL_ONLY / NOT_HFSS_VALIDATED. 2026-11-10, S2.

Decision session's revised ENV_ONLY criterion accepted and recorded in
ENV_DIAGNOSIS.md. Subsequent pre-MoM geometry/circuit checks: 3 passed,
1 failed (geometry_checks.log). No MoM search or baseline launched.

New failure is a physical-model consistency check, not a platform hash:
the supplied BAR50-02V OFF circuit and supplied HFSS sheet reactance disagree.

At f=5.8 GHz, circuit is Z_OFF=1/(1/4000+j*omega*0.11e-12)
+j*omega*0.6e-9. Sheet impedance is terminal impedance times W/L=0.35.

- Parallel RC alone: sheet X = -86.97243773015845 ohm/sq.
- Series 0.6 nH: sheet X = +7.652919704144736 ohm/sq.
- Combined OFF: sheet X = -79.31951802601371 ohm/sq.
- Handoff HFSS OFF sheet X: -86.97 ohm/sq.

The handoff value matches the RC-only calculation, not the stated series-L
circuit. ON sheet impedance remains 1.05+j7.652919704144736 ohm/sq.
See circuit_discrepancy.json for full-precision computed values.

No src/, tests or fixtures changed this turn. The new test was prepared
earlier and now exposes this inconsistency. Do not alter the test to hide
the mismatch. The user's hard-stop rule applies; no commit/push performed.

Decision needed: confirm the requested series-L circuit is authoritative
and reconcile/update the HFSS OFF impedance, or provide a revised OFF
circuit. The decision belongs to the decision session; no assumption made.
