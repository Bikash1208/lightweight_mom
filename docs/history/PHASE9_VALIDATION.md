# Physical control-variable + local-periodicity validation (Phase 9)

Two independent validation questions, resolved (or honestly left unresolved) before geometry
inversion may be considered:

```
QUESTION A: which physical complex response (S11/reflection, S21/transmission, or another
            demonstrably appropriate quantity) controls the holographic modulation?
QUESTION B: can a unit-cell response from a periodic infinite-cell simulation predict the
            response of that cell embedded in a heterogeneous metasurface?
```

## Phase 9A: control-variable validation

**No HFSS export exists anywhere in this repository** (grep-verified, the same finding every
earlier phase already made). Per Rule 1 ("if information is unavailable, explicitly report
UNRESOLVED / NOT AVAILABLE FROM CURRENT DATA"), Phase 9A builds the complete comparison
infrastructure (`hfss_reference.py`'s CSV ingestion schema, `error_metrics.py`'s complex/magnitude/
wrapped-phase statistics, `control_variable.py`'s 10-reference-cell x 7-frequency extraction from
the real Phase 7 checkpoint) but **the production run honestly reports `NOT_AVAILABLE_FROM_CURRENT_DATA`
for every HFSS-comparison metric**, never a fabricated number. Both S11 and S21 classify as
`UNVALIDATED` for that reason. `classify_control_variable` is fully implemented and unit-tested
against synthetic test data (including detecting an objective 180°-sign-mismatch pattern as
`INCONSISTENT`) so it is ready the moment a real HFSS export is supplied.

## Phase 9B: local-periodicity validation

A genuinely new, real capability, not a simulation of one: a heterogeneous N x N "supercell" is
built as one larger-period `UnitCell`-duck-typed object (`supercell_geometry.py`, multiple
sub-patches, each its own L/W) and solved through the **exact same, unmodified**
`solver.assemble_system`/`solve` MoM pipeline Phases 1-7 already use -- no coupled-cell solver was
invented, no new EM formulation, no new Green's function. The centre sub-cell's induced-current
dipole moment (`dipole_extraction.py`, reusing `basis_functions.RooftopBasis.evaluate` directly) is
compared between the isolated (periodic single-cell) and embedded (heterogeneous supercell) cases.

**Correctness proof**: a homogeneous 3x3 supercell (all 9 sub-cells identical) is mathematically
required to reproduce the isolated single-cell result exactly, since "9 identical cells repeated
periodically" and "1 cell repeated periodically" describe the same physical structure. Measured:
relative error **4e-15** (machine precision) once the isolated and embedded solves use matched
mesh/Floquet accuracy (`floquet_order = 2 x total mesh cells per side`, applied identically to
both) -- this is the load-bearing check that the whole supercell methodology is numerically sound.

### The real measured result

| scenario | relative complex error | phase error [deg] |
|---|---|---|
| 3x3 weak contrast | 0.159 | 8.46 |
| 3x3 moderate contrast | 0.060 | 3.20 |
| 3x3 symmetric (A-B-A/B-C-B/A-B-A) | 0.131 | 5.28 |
| 3x3 strong contrast | **0.997** | **52.5** |
| 5x5 moderate | 0.090 | 3.91 |
| 5x5 strong | 0.744 | 15.9 |

Strong-contrast neighbors (the sweep's own 1.0/4.4 mm extremes) change the centre cell's induced
current by an amount **comparable to the signal's own magnitude** -- a substantial, self-evident
deviation from the periodic assumption, independent of wherever an application-specific tolerance
would eventually be drawn. Classification: **`INVALIDATED`** (`classify_local_periodicity`,
Sec. 24's own criterion: not a chosen percentage, but "at least half the signal's own magnitude").

The 3x3-vs-5x5 comparison (`locality_analysis.json`) shows the strong-contrast case's deviation
*decreases* slightly from 3x3 to 5x5 (0.997 -> 0.744) -- consistent with the nearest neighbors
dominating -- but the moderate-contrast case's deviation *increases* (0.060 -> 0.090), showing
influence is not cleanly localized to the immediate neighborhood in every regime.

## Final status: `LOCAL_PERIODICITY_INVALIDATED`

Per the two-gate system (Sec. 29-30): Gate A (control-variable validation process) is `CLOSED`
(the HFSS-architecture-documented checkbox cannot be `True` when no HFSS project exists); Gate B
(local-periodicity validation process) is fully `OPEN` (every required test was actually run).
Geometry inversion requires both gates open AND a validated control variable AND
validated-or-partially-validated local periodicity -- none of that holds here, so geometry
inversion (and final STEP generation) remain **LOCKED**.

## Running it

```bash
.venv/bin/python scripts/research/run_phase9_validation.py     # results/phase9_validation/
```

## Verification

`tests/test_phase9_*.py` (68 tests): the homogeneous-supercell exact-match proof
(`test_phase9_supercell_solver.py`), the dipole-extraction cross-check against the basis's own
closed-form DC Fourier-transform integral (`test_phase9_dipole_extraction.py`), the two mandatory
Sec. 39 negative tests (unvalidated S-parameters and unvalidated local periodicity each correctly
keep geometry inversion locked, `test_phase9_gate.py`), and a direct end-to-end check that the real
run reproduces the measured strong-contrast invalidation (`test_phase9_run.py`).
