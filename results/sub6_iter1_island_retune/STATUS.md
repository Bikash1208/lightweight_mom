# S2 island retune: HARD_STOP_MESH_NOT_CONVERGED

Final 2026-11-10: 21 one-variable + 24 combined evaluations (40 unique),
top-three mesh96/128 and all OFF C/Rp corners completed. No qualified cell:
nominal suppression 96 -> 128 is G1 19.2455 -> 15.2031 dB,
G2 13.9615 -> 6.7164 dB, G3 13.5730 -> 1.3535 dB. ON alpha changes
176.176%, 445.889%, 1672.462%. User's >5% hard stop applies. Frequency
screening not launched; resonance/stop-band clearance NOT_ESTABLISHED.
Mesh-128 ON beta values are far-extrapolated two-node model outputs;
nominal fixed-point residuals meet 1e-5 but do not validate that extrapolation.
See REPORT.md, mesh_report.json, sheet_iteration_audit.json.

Previous-env full fast suite: 2954 passed, 129 skipped, six accepted ENV_ONLY
hash failures, one app timeout which passed isolated rerun (68.55 s).
New geometry/circuit checks: 4 passed. Solver src and original tests/fixtures
unchanged. Existing FFT assembly and exact load updates numerically audited.
HFSS OFF impedance correction documented. Prepared for commit/push github+hub.

The sections below preserve execution history.

Decision session accepted ENV_ONLY criterion and confirmed the series-L
OFF circuit. HFSS_SETUP.md and PACKAGE_9CELL_PLAN.md corrected to OFF
5.424-j79.320 ohm/sq; ON unchanged. The new research-only circuit check's
expected value now follows the authorized circuit; all four new checks pass.
Existing regression tests and fixtures unchanged. Baseline (F128) and
21-geometry mesh-96 search (F192) running in previous S2 environment with
single-threaded BLAS. No model result accepted before baseline validation.
The full previous-env fast suite is also running; six CP2 hashes fail under
the explicitly accepted ENV_ONLY criterion (initial run stopped at three).

The following sections preserve the earlier stops as history.

Latest: ENV_ONLY criterion accepted by decision session, then a new physical
circuit consistency check failed before MoM launch. See CIRCUIT_DISCREPANCY.md:
specified OFF series-L model gives Xs=-79.32, while handoff HFSS Xs=-86.97
matches RC without series L. HARD_STOP_OFF_CIRCUIT_HFSS_MISMATCH.

Follow-up diagnosis: see ENV_DIAGNOSIS.md. Polygon failure ENV_ONLY;
S11 comparisons against original entries unavailable (fixture has hashes
and norms only). Previous/new runtime S11 arrays match exactly; this does
not establish the user's <=1e-9 error threshold vs the original reference.
Retune remains stopped; no src/tests changes during diagnosis.

2026-11-10. MODEL_ONLY / NOT_HFSS_VALIDATED. No retune MoM run started.

Isolated worktree C:/lwmom_retune, branch compute/sub6-iter1-island-retune,
base 1ea2d13fe196556fef20a1ca42aed1037db0f12c. Fetched github and hub.
Original C:/lwmom_final data preserved; no use of its modified/untracked
geometry artifacts. C:/lightweight_mom untouched.

Initial fast suite command:

```
.venv313/Scripts/python.exe -m pytest -m "not slow" -q -n logical --maxfail=1
```

64 logical workers, OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=MKL_NUM_THREADS=1.
Machine WIN-D4BKN5G2BK3; Python 3.13.17; NumPy 2.5.4; pytest 9.2.0;
dependencies installed into a new isolated .venv313 with pip install
-e .[dev,app,cad] pytest-xdist threadpoolctl.

Exit 1: 4 failed, 1514 passed, 18 skipped, 73 warnings, 177.13 seconds;
suite interrupted after first failure, with already dispatched work finishing.
Full log: ../retune_initial_fast_suite.log.

Failures:

- tests/test_hologram_ports_polygon.py::test_hex4_alternative_size_input_gap_and_orientation:
  pytest.approx rejects nested vertex data structures (line 91).
- tests/test_loads_rotated_regression.py::test_preexisting_case_is_bit_identical_to_main[cp2_ideal_0100_s11]
- tests/test_loads_rotated_regression.py::test_preexisting_case_is_bit_identical_to_main[cp2_madp_1100_s11]
- tests/test_loads_rotated_regression.py::test_preexisting_case_is_bit_identical_to_main[cp2_ideal_0000_s11]

The three regression failures compare SHA256 fingerprints. Numerical error
magnitudes have NOT been computed; do not assume they are harmless platform
differences. No existing test or solver modified, no failure bypassed.

User requires hard stop on suspicious findings; CLAUDE.md requires fast
suite green before each commit. Therefore no commit or push. Prepared job,
driver and geometry/circuit checks remain uncommitted and unvalidated.
Top-three geometry, mesh convergence, sensitivity and resonance results:
NOT_AVAILABLE_FROM_CURRENT_DATA.

Read-only original diode_sheets.csv comparison against 1ea2d13:
sheet_z_thick moved immediately after sheet_centre_y; all 240 rows and
all column-named field values unchanged.
