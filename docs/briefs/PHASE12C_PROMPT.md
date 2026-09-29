# Phase 12C brief — does local periodicity hold for the grounded 8 mm reflectarray cell?

You are working in `lightweight_mom`. Read `CLAUDE.md`, `CONVENTIONS.md`, `PHASE9_VALIDATION.md`,
`PHASE11_HETEROGENEOUS_MOM.md`, `PHASE12A_GROUND_PLANE.md` and `PHASE12B_MESH_PERIOD_STUDY.md` first.

## Context

- Phases 9–11 found local periodicity invalid and coupling nonlocal/pattern-dependent — but for the
  **ungrounded** stack, on an **8-cell/side (0.6 mm) mesh** that Phase 12B showed cannot distinguish some
  geometries (C2 ≡ C3).
- The design is now **grounded**, and Phase 12B's period study favours **P = 8.0 mm** (279.2° coverage at
  12.5 GHz, resonance near L ≈ 6.0 mm, grating onset above band even at 45°, full surface-wave margin).
- Grounded patch reflectarrays are routinely designed cell-by-cell with local periodicity when the
  geometry varies smoothly. Whether that holds for this stack is the question that decides Phase 12's
  forward model: lookup table (if it holds) vs full-aperture integral-equation solver (if it doesn't).
- Polarization: linear x/y (circular derived via the existing basis). Rotated elements are out of scope.

Working design point: RT/Duroid 5880, h = 1.575 mm, PEC ground (`backing: pec`), P = 8.0 mm, rectangular
PEC patch, 10–15 GHz, design frequency 12.5 GHz. Record it as a config; don't hard-code it.

## Tasks

### A. Grounded supercells
- Make the heterogeneous-supercell path (Phase 9's `SupercellUnitCell`, used by Phases 10–11) work with
  `backing: pec` **without editing the frozen Phase 9 class** (wrap/compose, as Phase 12A did with
  `getattr(..., "backing", "air")`). Test: a homogeneous grounded N×N supercell reproduces the single
  grounded cell to solver precision (the Phase 9 identity check, grounded).
- Mesh: use aligned meshes (Phase 12B helper). Supercells are costly — choose a per-cell resolution and
  a geometry grid that are mutually aligned, verify on a 3×3 case that the result converges against a
  finer aligned mesh, and report unknowns, runtime and memory. If a size is infeasible, say so and stop at
  the largest feasible one rather than coarsening silently.

### B. Local-periodicity test, grounded 8 mm (the Phase 9 test, redone)
Measure the centre-cell response (Phase 9's own metrics: induced-current dipole moment / complex
deviation, magnitude, phase) embedded vs in its own uniform array, for:
1. **Smooth gradients** — the realistic reflectarray case: neighbours differing from the centre by
   small steps along the phase-vs-L curve (e.g. ±15°, ±30°, ±45° of reflection phase), in 1-D (row) and
   2-D gradients, centres chosen off-resonance, near resonance and at resonance.
2. **Strong contrast** — the Phase 9 extremes, for comparison with the ungrounded result.
3. **Periodic phase-gradient supercell** — N cells spanning one full 2π phase cycle (a standard
   anomalous reflector; pick N from the design, e.g. 5–8), which is exactly periodic so no artificial
   boundary jump. Report each cell's embedded vs isolated response. Optionally, from the current's
   Floquet coefficients, report the power in the intended diffraction order vs the local-periodicity
   prediction — document the method and validate it on a trivial case before using it.
Report deviations against a measured noise floor, as a function of phase step size. The key output is
the **largest phase step per cell for which local periodicity holds within a stated tolerance** (e.g.
≤ 5–10° phase / ≤ 0.5 dB), or an explicit finding that none does.

### C. Coupling extent, grounded (short Phase 11 series)
- Growing-pattern coupling series (Phase 11 method, reusing its functions) for 2–3 grounded centre
  geometries up to the largest feasible size; include a gradient pattern, not only a checkerboard.
- Grounded slabs support a TM0 surface wave; report whether coupling decays more slowly than the
  ungrounded case and whether that correlates with the surface-wave margin. Don't fit a decay law
  unless the data are monotonic (Phase 11 rule).

### D. Spot-check Phases 9–11 on a finer aligned mesh (read-only for earlier results)
- Re-solve on an aligned mesh finer than 8/side: Phase 11's checkerboard-vs-block Case D pair, one
  Phase 11 coupling-radius series (to the largest feasible size), and Phase 9's strong-contrast 3×3.
- Re-verify Phase 10's C2 vs C3 on an aligned mesh and relabel if still coincident.
- For each, state whether the earlier conclusion holds, changes quantitatively only, or changes. Write
  results to `results/phase12c_*`; never edit earlier outputs.

### E. HFSS-ready exports
- STEP files for this set already exist in `hfss_exports/phase12_grounded_8mm/` (see its README and
  `manifest.json`): singles G1–G6 (L = 2.0/4.8/5.6/6.0/6.8/7.6 mm, with a MoM reference CSV) and 3×3 supercells
  S1 (uniform 5.6), S2 (columns 5.2/5.6/6.0), S3 (centre 5.6, neighbours 1.6). Use **exactly these layouts**
  in Tasks B/C where they fit, compute the missing MoM reference for S1–S3 (specular (0,0) S11 at the 5 mm
  plane, 10–15 GHz at 0.25 GHz; note the 24 mm period passes λ0 at 12.49 GHz), and check the existing STEP
  files against `cad_export.import_and_summarize`-style bounds rather than regenerating them.
- Otherwise, export STEP (with the ground solid) plus MoM reference results for a small, deliberate set of grounded
  8 mm cases: ~5 single cells spanning the phase curve (incl. resonance), one uniform 3×3, and one
  smooth-gradient and one strong-contrast 3×3 supercell. Use the Phase 6/9 export conventions
  (reference planes, time-convention note, reproduction recipe). Make sure the Phase 9 comparison
  tooling can ingest an HFSS CSV for these cases (`--hfss-reference`); test it on a synthetic file only.
  Do not create or imply HFSS results.

### F. Ingest the real HFSS results (now available)
Real HFSS 2025 R1 data exist in `hfss_exports/phase12_grounded_8mm/hfss_results/` (G1 done; G2–G6 to follow
as `G<n>.s36p`). Build a proper ingestion path, not a one-off script:
- Parse Touchstone N-port files (`# GHz S MA`, N from the extension, `! Port[k] = FloquetPortX:m` headers).
  Top port = FloquetPort1, 18 modes; modes 1/2 = TE00/TM00 (for a square patch both are co-pol; record which is
  x-polarised if it can be established, otherwise compare both). Bottom port FloquetPort2 must show ~0 transmission.
- HFSS setup (record it as provenance): port on the air-box top at z = +5 mm, no de-embedding (matches the MoM
  `ReferencePlanes(5.0e-3, None, h)`), 18 Floquet modes (≤ 9.25 dB/mm at 15 GHz), PEC ground, PEC patch 35 µm
  thick, max ΔS 0.001, 10–15 GHz, 401 points. exp(+jωt) in both — no conjugation.
- Feed it into Phase 9's comparison tooling (`hfss_reference.py`, `error_metrics.py`) and the gate records,
  so `independent_hfss_validation_available` / `hfss_comparison_meets_tolerance` are evaluated on real data
  for the grounded single cells. Define the tolerance up front and justify it.
- Reference result already measured for G1 (reproduce it): phase difference −0.13° (10 GHz) to −0.52° (15 GHz),
  mean −0.29°, |S11| within ~1e-5, cross-pol ≤ 2e-5. An earlier run with 2 port modes and a copper ground was
  off by up to −4.7°: document that the 2-mode port was the cause, as a setup lesson.
- The remaining ~0.5° trend is plausibly the 35 µm patch thickness (MoM is zero-thickness); state it as a
  hypothesis unless a zero-thickness HFSS run is supplied.
- Only single-cell HFSS gates can open; supercell/heterogeneous gates stay closed until S1–S3 are run.

## Deliverables
- `examples/run_phase12c.py` → `results/phase12c_grounded_locality/` (CSV/JSON, plots: deviation vs
  phase step, embedded vs isolated per cell for the phase-gradient supercell, coupling vs size grounded
  vs ungrounded, fine-vs-coarse spot-check table), plus `results/phase12c_hfss_export/`. Deterministic.
- `PHASE12C_GROUNDED_LOCALITY.md`; updates to `README.md` (status, layout, "Still open"), `CLAUDE.md`.
- Explicit statuses from fixed sets, e.g.:
  `GROUNDED_LOCAL_PERIODICITY_VALID_FOR_SMOOTH_GRADIENTS` (with the max phase step and tolerance) /
  `GROUNDED_LOCAL_PERIODICITY_INVALIDATED`; spot-check `EARLIER_CONCLUSIONS_HOLD` /
  `EARLIER_CONCLUSIONS_QUANTITATIVELY_CHANGED` / `EARLIER_CONCLUSIONS_CHANGED` (listed);
  `HFSS_EXPORT_READY`. HFSS status stays `HFSS_NOT_AVAILABLE` until real data are supplied; geometry
  inversion and STEP export of a final design stay LOCKED.
- A one-paragraph recommendation for Phase 13's forward model (lookup table with a phase-step limit vs
  full-aperture IE à la Budhu & Grbic), grounded in the measured numbers — as a recommendation, not a
  decision.

## Rules
- Measured numbers only; label synthetic data; no fabricated HFSS results.
- Don't modify earlier phases' code or results; extend by wrapping.
- Don't choose geometry or period by phase coverage alone; report tolerances you used and why.
- No ML/fitted surrogates in the locality/coupling study (Phase 11 rule).
- Keep tests property-based; run the full suite at the end and report counts and runtime.
