> **PARKED 2026-10-13 (Phase 24):** the RIS work (Phase 20) is paused, not started. Nothing below is current work until the user re-activates it.

# Phase 20 / R1 — MoM solver: lumped loads, oblique incidence, composite geometry

**Status:** active brief. Read `docs/START_HERE.md`, `docs/STATUS.md`, `CLAUDE.md`, then
`docs/briefs/PHASE20_RIS_ROADMAP.md` for context. Work on a new branch `phase20-ris` from `main`.

## Why
A 2-bit PIN-diode RIS cell (roadmap) needs three things the periodic MoM does not have: a diode (lumped RLC) across a
gap in the metal, illumination at an angle, and a patch geometry with gaps. All three are **additive extensions of
frozen core modules** — a documented exception like Phase 12A: every pre-existing config, sweep and result must
reproduce **bit-identically**, and the extension must be explicitly listed in `CLAUDE.md`'s "documented, additive
exceptions".

## What already exists (check before writing code)
- `src/excitation.py: PlaneWaveExcitation` already takes `theta`, `phi`.
- `src/green_function.py: floquet_dyads(lattice, kx0, ky0, ...)` and `src/mom_matrix.py:
  assemble_impedance_matrix(basis, green, kx0, ky0, ...)` already take the incident transverse wavenumber.
- Blockers: `src/config.py` rejects `theta_deg != 0`; `src/sparameters.py: compute_s_parameters` passes `0.0, 0.0`
  to assembly, excitation vector and `current_coefficients`; `scattering_from_fields` / `_scattering_from_fields_
  grounded` decompose the (0,0) fields with `stack_fields(0.0, 0.0, ...)`, i.e. assume normal incidence.

## Work items (one commit each, `make test-fast` green before each)

**1. Oblique incidence (grounded and ungrounded stacks).**
- Thread `(kx0, ky0)` from the excitation through assembly, excitation vector, current coefficients and the (0,0)
  field decomposition. Define the port modes at oblique incidence as **TE and TM** (with respect to the plane of
  incidence), document the basis and reference planes in `docs/theory/S_PARAMETERS.md` and
  `docs/theory/CONVENTIONS.md` (`exp(+jωt)`; record how TE/TM map to the existing x/y labels at θ = 0).
- Allow `theta_deg`, `phi_deg` in config (0 ≤ θ < 90°); warn (existing mechanism) when a higher Floquet order
  propagates at that angle; record it in sweeps.
- Tests: θ = 0 reproduces every existing result **bit-identically**; bare grounded slab and bare ungrounded slab vs
  analytic TE/TM Fresnel/Fabry–Pérot at 0°, 30°, 60° (to ~1e-10); lossless energy conservation (|S11|² summed over
  propagating orders = 1 grounded; |S11|²+|S21|² = 1 ungrounded) at oblique angles; reciprocity at oblique incidence
  (θ,φ) ↔ (θ,φ+180°) for a symmetric cell; for a square patch, TE at φ = 0 equals TE at φ = 90° rotated.
  Optional cross-check against the Phase 19 RCWA (`src/dielectric_cell`, branch `phase19-transmitarray`) for dielectric-
  only stacks — only if that branch is merged; don't merge it just for this.

**2. Lumped RLC loads.**
- A load = a small rectangular region bridging a gap, with a direction (x or y) and an impedance `Z_L(f)`
  (series R–L–C; also allow open and short). Inside the region the mesh treats it as conductor, so rooftop basis
  functions cross it; the load adds its impedance to the MoM system for the basis functions crossing the gap
  (standard lumped-element loading; derive the per-edge normalisation in `docs/theory/MOM_FORMULATION.md` — N parallel
  edges across a gap of width w each carry Z_L·N, check the length/width scaling against the rooftop normalisation).
- Config: `loads: [{x, y, length, width, direction, R, L, C}]` (SI, strict keys as elsewhere; absent = no loads =
  bit-identical old behaviour).
- Tests (physics, not "it runs"): Z_L = 0 reproduces the same geometry with the gap filled with metal (to solver
  tolerance); Z_L → open (C → 0 / |Z| → ∞) reproduces the geometry with the gap empty; a purely resistive load
  absorbs exactly 1 − |S11|² (grounded, lossless substrate) — compute absorbed power from the load current and compare;
  reactive load: lossless, |S11| = 1 grounded; mesh-refinement convergence of the phase for one loaded cell; results
  independent of how many mesh edges the load spans (within tolerance) once normalised.

**3. Composite geometry for gapped patches.**
- A geometry type built from a union of metal rectangles minus aperture rectangles (axis-aligned is enough), so a
  patch with two gaps and arms of different lengths can be described. Mesh alignment: reuse
  `discretization.feature_edges_aligned` logic and warn when any edge is off-grid. CAD export
  (`src/cad_export.py`) gets the same geometry, with each load drawn as a separate named sheet ("Load_1", …) so it can
  be assigned an RLC boundary in HFSS.
- Tests: a composite of one rectangle reproduces the `rectangle` geometry bit-identically; area/mask checks; STEP
  round-trip (one Metal body, N Load sheets, bbox).

**4. Demonstration case (not a design yet).**
- Use the roadmap's *starting values* (10 mm cell, RT/duroid 5880 1.575 mm, PEC ground) but make period, thickness and
  geometry parameters of the demo script, not constants — later phases may change them. A simple patch split by one gap with one load; sweep the load
  from short to open at 10 GHz, θ = 0, 30°, 60°, TE/TM. Report phase/|S11| vs load and angle in
  `results/phase20_r1_demo/README.md`. Write a one-cell HFSS check package (`hfss_exports/phase20_r1_loaded_cell/`:
  STEP + RLC/Floquet setup notes + MoM reference CSV) so the user can validate the loaded, oblique solver in HFSS.

**5. Close-out.** Update `docs/STATUS.md` (capability matrix: "lumped loads", "oblique incidence" rows,
model-only until the user's HFSS numbers exist; known-limitations list: remove "normal incidence only"),
`CHANGELOG.md`, `CLAUDE.md` exception list; move this brief to `docs/history/briefs/`.

## Do not
- Change any result for an existing config (bit-identical at θ = 0 and with no loads is a hard gate).
- Model vias or vertical currents (out of scope; the RIS cell's centre via is handled as a labelled approximation).
- Start R2–R6.
