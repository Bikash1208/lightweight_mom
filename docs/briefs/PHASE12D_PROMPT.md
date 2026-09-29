# Phase 12D brief — resolve the local-periodicity step limit, and ingest the real HFSS data

Read `CLAUDE.md`, `CONVENTIONS.md`, `PHASE12C_GROUNDED_LOCALITY.md`, `results/phase12c_grounded_locality/task_b.json`
and `hfss_exports/phase12_grounded_8mm/README.md` first. Design point unchanged: grounded RT/Duroid 5880,
h = 1.575 mm, P = 8.0 mm, square PEC patch, 10–15 GHz, design frequency 12.5 GHz.

## Why

1. **Phase 12C's "local periodicity holds up to 2.7°/cell" is an artefact of the test grid, not a measured limit.**
   `task_b.json` has only four gradient cases: one at 2.66°/cell (passes) and three at 81–107°/cell (fail). The
   brief asked for ±15/30/45° steps, but the 0.4 mm L grid cannot produce small steps near resonance. The limit
   lies somewhere between 2.7° and 81°/cell and has not been measured. Also, the failing cases are not
   catastrophic in design terms: e.g. `near_resonance_1d_row` at 81°/cell gives 6.7° phase / 1.7 dB error.
2. **The metric is a proxy.** The per-cell dipole-moment deviation (Phase 9's metric) is not what the design needs.
   What matters is the reflected field: reflection phase/amplitude error per cell and, at aperture level, the
   power and phase in the intended diffraction order.
3. **Phase 12C Task F (HFSS ingestion) did not run**, but real HFSS data now exist:
   `hfss_exports/phase12_grounded_8mm/hfss_results/unit_HFSSDesign1.s36p` (G1, L = 2.0 mm). More cases will be
   added as `G<n>.s36p`.

## Tasks

### A. Fill in the step-size curve
- Use a finer L grid near resonance (e.g. 0.05–0.1 mm, with meshes aligned to it via the Phase 12B helper, and
  `floquet_order = 2 × total mesh cells` as Phase 12C found necessary) so neighbours can differ from the centre by
  **≈ 10, 20, 30, 45, 60, 80° of reflection phase**, for off-, near- and at-resonance centres, 1-D rows and 2-D
  gradients. Report the achieved steps honestly.
- For each case report, alongside Phase 9's dipole metric: the **equivalent reflection-phase and amplitude
  error** of the centre cell (define the conversion from current/dipole to reflection coefficient, validate it on
  a uniform array where it must reproduce S11 exactly).
- Plot error vs step size. Fit nothing unless the data are monotonic; report the step size at which the error
  crosses a stated tolerance (use a design-relevant one, e.g. 10° / 1 dB, and also show 5° / 0.5 dB).

### B. Aperture-level test (the design metric)
- Periodic phase-gradient supercells (anomalous reflectors) with N = 6, 8, 12 cells per 2π cycle
  (≈ 60°, 45°, 30°/cell) on the grounded 8 mm cell, cells picked from the isolated-cell lookup table.
- Solve each full supercell with MoM and compute the power in every propagating Floquet order and the phase of
  the intended order (document the method; validate on a uniform array: all power in (0,0), and on a trivial
  case with a known answer). Compare with the local-periodicity prediction from the lookup table.
- Report: fraction of reflected power in the intended order (MoM vs prediction), spurious-order levels, and the
  per-cell reflection-phase error profile. This decides Phase 13, not the dipole metric.
- Largest feasible sizes only; report runtime/memory. If N = 12 is infeasible, say so.

### C. HFSS ingestion (Phase 12C Task F, carried over)
- Touchstone N-port reader (`# GHz S MA` or `RI`, N from the extension, `! Port[k] = FloquetPortX:m` headers).
  Top port = FloquetPort1 (18 modes; modes 1/2 = TE00/TM00, co-pol for a square patch); bottom port must show
  ~0 transmission. Wire it into `phase9/hfss_reference.py` / `error_metrics.py` and the gate records.
- Provenance to record: HFSS 2025 R1; Floquet port on the air-box top at z = +5 mm, no de-embed (= MoM reference
  plane); 18 modes (up to 9.25 dB/mm at 15 GHz); PEC ground; PEC patch 35 µm thick; max ΔS 0.001;
  10–15 GHz, 401 points; exp(+jωt) both — no conjugation.
- Must reproduce the already-measured G1 result: phase difference −0.13° (10 GHz) → −0.52° (15 GHz), mean −0.29°,
  |S11| within ~1e-5, cross-pol ≤ 2e-5. Record the setup lesson: an earlier 2-mode-port / copper-ground run was
  off by up to −4.7°.
- Define and justify the tolerance up front. Evaluate `independent_hfss_validation_available` and
  `hfss_comparison_meets_tolerance` for the grounded single cells on real data; the gate should say which cases
  have HFSS data (currently G1 only) and that the resonance case (G4) is still pending. Supercell HFSS gates
  stay closed. Remaining ~0.5° trend: plausibly the 35 µm patch thickness — a hypothesis, not a finding.
- Phase 12C re-exported STEP files to `results/phase12c_hfss_export/`; the ones the user is running are in
  `hfss_exports/phase12_grounded_8mm/`. Check they're geometrically identical and say so; use the latter's IDs.

## Deliverables
- `examples/run_phase12d.py` → `results/phase12d_*` (CSV/JSON; plots: error vs step size, per-order power for
  each N, per-cell phase error, HFSS vs MoM for each available case). `PHASE12D_*.md`; README/CLAUDE.md updates.
- Statuses from fixed sets, e.g. `LOCAL_PERIODICITY_STEP_LIMIT_MEASURED` (with the number and tolerance) /
  `..._NOT_RESOLVED`; aperture test `LOOKUP_TABLE_ADEQUATE_FOR_N>=…` / `LOOKUP_TABLE_INADEQUATE`;
  `HFSS_SINGLE_CELL_VALIDATED_PARTIAL` (cases listed) etc.
- Replace Phase 12C's "≤ 2.7°/cell" recommendation with the measured one; don't edit 12C outputs in place.

## Rules
Measured numbers only; don't fabricate or extrapolate HFSS data; earlier phases untouched; no ML; property-based
tests; full suite at the end with counts and runtime.
