# Phase 19 — 3D-printed nylon transmitarray (MTS lens), 28 GHz, dual beam ±30°

**Status:** active brief, not started. Written 2026-10-08 from the user's request. Read `docs/START_HERE.md`,
`docs/STATUS.md` and `CLAUDE.md` first; every hard rule there applies (tags, honest labels, one commit per item,
fast suite green before each commit, no silent changes to frozen solver modules).

## 1. What the user wants

A horn-fed, all-dielectric transmitting metasurface (transmitarray / MTS lens), 3D-printed in nylon, that turns one
horn's beam into **two beams at θ = ±30°** at **28 GHz**.

**The deliverable is the MTS structure as a STEP file** (plus STL for printing). The user runs all HFSS verification
themselves — this phase does not run or wait on HFSS; it produces the design, the model's own predictions to compare
against, and the STEP.

| Input | Value | Source |
|---|---|---|
| Frequency | 28 GHz (λ0 = 10.707 mm) | user |
| Material | nylon, **εr = 2.45, tanδ = 0.0185** at 28 GHz (n = 1.565) | user (measured/known) |
| Thickness | ≈ 25 mm | user |
| Unit cell | octagonal pillar, "max radius 9.8 mm" — **see §2, must be confirmed** | user |
| Feed | horn on axis, **F/D = 1.2**, linear polarisation | user (F/D) + agreed default |
| Beams | θ = ±30°, both in the x–z plane (φ = 0° and 180°), equal level | agreed default |
| Aperture outline | octagonal (assumed — the lens outline, as well as the cell) | **confirm** |
| Aperture size D | **not given** — default D = 150 mm across flats (≈ 14 λ0), parameterised | **confirm** |

### "99% accuracy" — the agreed, measurable meaning
Accuracy is measured **against HFSS, run by the user after delivery** (print tolerance and εr spread are outside the
model). The design must be built so that it can meet these targets; the user's HFSS run decides whether it does:
- Unit cell: transmission phase within **±3.6°** (1% of 360°) and |S21| within **±0.1 dB** of an HFSS Floquet-port
  solve, for every cell in the final library, at normal incidence and at the oblique angles actually used (§4).
- Full lens: each beam's peak within **±1°** of HFSS, realized gain within **±1 dB**, and the two beams' level
  difference within **±0.5 dB**.
Until the user returns HFSS results, every prediction is labelled model-only (`NOT_AVAILABLE_FROM_CURRENT_DATA` for
the HFSS comparison) — never claim the target is met without their data. Because no HFSS loop happens inside this
phase, the model itself must be as tight as possible: converged RCWA (Fourier order and layer handling shown converged
to well below 1°), the actual oblique-incidence cell response used per cell, and the actual complex transmission
(not the ideal phase) used in the pattern prediction.

## 2. Physics constraints the design must respect (computed, not assumed)

- **Cell size.** A cell bigger than λ0/(1 + sin 30°) = **7.14 mm** puts a grating lobe in visible space for a ±30° beam;
  homogenisation inside nylon also wants p < λ0/n = **6.84 mm**. So the *cell pitch* cannot be 2 × 9.8 mm = 19.6 mm.
  **Working interpretation (confirm with the user before Part B):** the cell pitch is p ≈ 5.0–5.35 mm (≈ λ0/2), each
  cell holds one octagonal nylon pillar whose circumradius varies from a printable minimum up to the cell limit, and
  9.8 mm is not a cell dimension. If the user insists on a 9.8 mm-radius octagonal cell, stop and report the grating-
  lobe consequence with numbers instead of building it.
- **Phase range.** Fixed 25 mm height with variable pillar fill: ideal span ≈ 360°·(n − 1)·h/λ0 = **475°**, so 360° is
  reachable with margin (variable *height* would need ≥ 18.9 mm). Measure the real span from the unit-cell solver;
  don't quote the ideal number as achieved.
- **Loss budget.** Solid 25 mm nylon at tanδ 0.0185 loses ≈ **1.8 dB**; Fresnel reflection per face ≈ −13 dB. The densest
  cells are the lossiest — report |S21| per cell and the aperture-averaged transmission loss.
- **Oblique incidence.** At F/D = 1.2 the edge cells see ≈ **22.6°** incidence. The existing solver is normal-incidence
  only (`docs/STATUS.md`); the new cell solver must handle oblique incidence and the design must use it.
- **Dual beam, phase-only.** A two-beam pattern needs amplitude *and* phase; a lens only sets phase. The naive
  phase arg(e^{+jk sin30 x} + e^{−jk sin30 x}) is binary (0/π) and wastes power into higher orders. Synthesize a
  continuous phase that maximises power into ±30° with equal levels (see Part C) and report the efficiency into the two
  beams vs the binary baseline.

## 3. Why the current app can't do this as-is

The solver (`src/` MoM) models a **zero-thickness PEC pattern on one thin substrate** — it cannot model a 25 mm
z-invariant dielectric pillar. Transmitarrays are parked in `docs/STATUS.md`. The surface-wave HMSA flow (monopole
`.fld` → hologram) is a different feed concept. What *is* reusable: the space-fed design logic in
`src/reflectarray/designer.py` (feed distance, required phase, illumination taper, library assignment), the vector
aperture far-field model (`src/pattern/`), the beam/SLL/XPD metrics, level matching ideas from Phase 16, and the CAD
export patterns (`src/cad_export.py`, `src/export/`).

## 4. Work items (one commit each; fast suite green before each)

**Part A — dielectric unit-cell solver (new package, do NOT modify the MoM core).**
- New package, e.g. `src/dielectric_cell/`: RCWA (Fourier modal method) for a periodic cell = air | one z-invariant layer
  (octagonal nylon pillar in air, thickness h) | air, plus an optional thin solid base plate layer. The pillar is
  extruded, so one uniform RCWA layer represents it exactly; contrast is low (εr 2.45), so convergence should be fast.
- Prefer a minimal in-repo implementation (one layer + base plate is a few hundred lines). If an external package is
  used, check its licence is acceptable for ISRO-internal use and pin it in `pyproject.toml`; no ML/fitted surrogates.
- Inputs: period, octagon circumradius, h, εr/tanδ, base-plate thickness, frequency, (θ, φ) incidence, TE/TM.
  Outputs: S11/S21 co- and cross-pol at defined reference planes (document them; mind the reference-plane lesson in
  `src/hologram/two_state_cells.py`), `exp(+jωt)` convention per `docs/theory/CONVENTIONS.md`.
- Tests (physics, not "it runs"): uniform slab (fill = 100% and 0%) vs the analytic Fabry–Pérot slab formula to
  1e-6; energy conservation |S11|² + |S21|² = 1 for tanδ = 0 and < 1 with loss; reciprocity; C4 symmetry of the
  octagon (TE/TM equal at normal incidence); Fourier-order convergence study recorded with the chosen order.

**Part B — 28 GHz nylon cell library.**
- Sweep octagon size (printable minimum feature, e.g. ≥ 0.6 mm wall/gap — make it a parameter) at the chosen pitch and
  h ≈ 25 mm, at normal incidence and at 0°/10°/20°/25° oblique, TE and TM. Pick pitch/h/base plate so the library
  spans ≥ 360° with the best worst-case |S21|. Record phase span, |S21| range and the per-cell loss.
- Also write a small, optional single-cell check set for the user (`hfss_exports/nylon_cells_28ghz/`): STEP for ~6 cells
  across the phase range + the RCWA reference CSV + Floquet-port notes, so they can spot-check cells in HFSS if they
  choose. Don't block on it.

**Part C — space-fed dual-beam lens design (`src/transmitarray/`).**
- Feed: horn modelled as cos^q(θ) with q from F/D = 1.2 (aim for ≈ −10 dB edge taper; make q a parameter), phase centre
  at F = 1.2·D on axis; per-cell incidence angle and path length.
- Required phase: compensate the feed's spherical phase + impose the dual-beam aperture phase. Synthesis: start from
  the binary baseline, then optimise a continuous phase (e.g. Gerchberg–Saxton / gradient on the aperture-integration
  pattern) for equal ±30° beams and minimum power elsewhere; report efficiency into the two beams, SLL, and level
  balance. Reuse Phase 16's level-matching measurement convention so numbers are comparable.
- Assign cells from the Part B library using each cell's *actual oblique-incidence* response, then predict the pattern
  from the assigned complex transmission (amplitude + phase, not the ideal phase) with the vector aperture model.

**Part D — export (the main deliverable).**
- **Full-lens STEP**: octagonal outline, all pillars + base plate fused into ONE valid solid, units mm, coordinates
  rounded (no 4.9999999 values), named body "Lens" so a single material assignment (nylon εr 2.45, tanδ 0.0185) covers
  it. Re-import the STEP and check: one solid, valid, bounding box = D × D × h as designed, volume matches the layout's
  own sum, pillar count matches. Also an STL for printing, a layout CSV/JSON (per-cell x, y, octagon radius), and a
  README with: horn position (phase centre at F = 1.2·D on axis), the predicted beam directions/gain/SLL/level balance
  to compare with HFSS, and suggested HFSS setup notes (radiation boundary, far-field sphere, linear co/cross).
- Package in `hfss_exports/nylon_lens_28ghz_dual30/` (STEP/STL gitignored like other packages; regenerating script recorded).

**Part E — app integration (optional, only if the user asks after receiving the STEP).**
- A "Transmitarray (dielectric, horn-fed)" design path in the app: inputs = frequency, material, thickness, D, F/D,
  beams; outputs = phase map, cell layout, predicted pattern, export. Restart the server and verify in a real browser.

**Part F — close-out.** Hand the user the STEP path first. Update `docs/STATUS.md` (capability matrix row: "model-only,
awaiting the user's HFSS"), `CHANGELOG.md`, `hfss_exports/README.md`; move this brief to `docs/history/briefs/`; tag
`phase19-complete`.

## 5. Open inputs to confirm with the user (ask before Part B, don't guess silently)
1. Cell pitch/octagon interpretation in §2 (9.8 mm cannot be the cell pitch).
2. Lens size D (default 150 mm) and whether the lens outline is octagonal.
3. Minimum printable feature size of their printer/nylon grade.
4. Horn: model or gain, polarisation, and the actual phase-centre distance.
5. Whether the lens needs mounting features (flange, holes) in the STEP.
