# lightweight_mom

A lightweight Python Method-of-Moments solver for **periodic metasurface unit
cells** (not a general-purpose full-wave solver).

```
Unit-cell geometry → Periodic MoM → Surface current → Reflected/transmitted
fields → S11/S21 (mag, phase) → Unit-cell database → Holographic synthesis
```

**Status: Phase 12D complete (foundation + Green's function + MoM + S-parameters + sweeps/database + validation export + production unit-cell library + scientific feasibility gate + holographic formulation & control-variable validation + unit-cell scattering architecture validation + universal unit-cell database architecture + physical control-variable/holographic-scattering bridge + control-variable & local-periodicity validation + environment-aware neighborhood/surrogate/inverse/beam validation + physics-based heterogeneous-MoM nonlocal-coupling study + optional PEC ground plane + mesh alignment/duplicate-geometry audit/grounded period study + grounded-supercell local-periodicity/coupling/fine-mesh spot-check + HFSS-ready exports + design-relevant step-limit re-measurement + aperture-level Floquet-order-power test + real HFSS Touchstone ingestion).**
Configuration, constants/conventions, materials, parametric geometry, Floquet bookkeeping, plane-wave excitation, the air | substrate | air spectral
dyadic Green's function ([GREEN_FUNCTION.md](../theory/GREEN_FUNCTION.md)), the rooftop MoM (`Z I = V`, [MOM_FORMULATION.md](../theory/MOM_FORMULATION.md)) and now
**S-parameter extraction** at the configurable top/bottom reference planes ([S_PARAMETERS.md](../theory/S_PARAMETERS.md)): complex S11, S21 (and S12, S22),
co- and cross-polarized, for x and y polarization at normal incidence. **Phase 5** adds frequency and geometry-parameter sweeps, automatic
mesh/Floquet validation, resumable runs, deterministic per-cell metadata and a CSV + JSON unit-cell database with linear and circular
(`e_R = (x+jy)/√2`, `e_L = (x−jy)/√2`) data ([DATABASE.md](../theory/DATABASE.md)). **Phase 6** exports 10 representative unit cells as physical STEP models
(substrate + metal solids, built directly from the same geometry the solver used) plus MoM reference results, plots and a per-case README, for
manual comparison against a full-wave solver ([VALIDATION_EXPORT.md](VALIDATION_EXPORT.md)); it changes no EM code and performs no automatic
HFSS comparison or tuning. **Phase 7** builds a production rectangular-PEC-patch unit-cell library (18×18 L/W sweep × 201 frequencies, 10–15 GHz,
on the same locked substrate/period as the brief) into `results/unit_cell_library/`: STEP geometry, full complex S-parameters, a design-frequency
phase-vs-geometry lookup table with honestly-quantified phase coverage, documented resonance/grating classification, plots and a read-only
query/phase-inversion API ([UNIT_CELL_LIBRARY.md](UNIT_CELL_LIBRARY.md)) — again no EM code changed and no HFSS agreement claimed. The library has
been built (324/324 geometries solved, 0 failed; measured — not assumed — S11xx phase coverage at 12.5 GHz is **24.8° of 360°**, largest gap
335.2°, 0 cells resonant in this 10–15 GHz band: see `results/unit_cell_library/reports/LIBRARY_REPORT.md` for the full, current numbers).
**Phase 8A** adds a scientific feasibility gate (`scripts/research/run_phase8a_feasibility.py`, [PHASE8A_FEASIBILITY.md](PHASE8A_FEASIBILITY.md)) that
checks, from the real Phase 7 database only, whether the library can support a requested holographic design before any GUI or final geometry
generation is attempted — library integrity, phase coverage/frequency-dependence/resonance, grating risk, and (design-specific) reference-field/
required-phase/matching checks, each returning one of a fixed set of explicit states, never a bare pass/fail. Run against the real library it
correctly reports **`PHASE_SOURCE_UNVALIDATED`**: this project has no implemented holographic generator, so whether S11xx is even the right
phase-control quantity cannot be checked — a genuine, expected scientific limitation, not a bug. **Phase 8B** establishes the mathematical
bridge Phase 8A found missing (`scripts/research/run_phase8b_formulation.py`, [PHASE8B_FORMULATION.md](PHASE8B_FORMULATION.md)): object-wave/multi-beam
field equations (reusing `excitation.PlaneWaveExcitation`'s convention, numerically tested), and an honest investigation of two literature-sourced
candidate holographic formulations against the real Phase 7 library. One (the classic monopole-fed, grounded holographic-impedance-surface
technique, Fong et al. 2010) is found **structurally incompatible** — this project's solver has no ground plane and no surface-wave dispersion
solver, both required by that technique's control variable. The other (a simpler reflectarray/transmitarray-style phase-difference hologram) needs
only a normal-incidence co-pol reflection/transmission phase — the *kind* of quantity Phase 7 stores — so it remains open, blocked only on the
same S11-vs-S21 architecture decision Phase 8A already flagged. Overall status: **`PHASE_SOURCE_UNVALIDATED`** (unchanged from Phase 8A, now with
a fully worked mathematical reason why). **Phase 8C** (`scripts/research/run_phase8c_validation.py`, [PHASE8C_SCATTERING.md](PHASE8C_SCATTERING.md)) puts
that question to real, full-library numeric tests: reading the full linear-basis S-matrix (not just the delivered library's "_xx" export) from the
still-present Phase 7 solve checkpoint, it measures S11xx phase coverage (24.84°) *and*, for the first time, S21xx phase coverage (43.85°, more
than S11xx — reported as data, never used to pick a winner, per the brief's explicit rule), confirms reciprocity and L↔W mirror symmetry to ~1e-15
on the real solved data, confirms negligible cross-polarization, and confirms only the fundamental Floquet mode propagates in 10–15 GHz. Neither
S11xx nor S21xx reaches `VALIDATED`; final status **`CONTROL_VARIABLE_UNVALIDATED`** — the physical interpretation (reflective vs. transmissive
MTS), not the phase coverage, is the actual open question. Holographic synthesis itself is still **not implemented**. **Phase 8D-REV**
(`scripts/research/run_phase8d_rev_universal_database.py`, [UNIVERSAL_UNIT_CELL_DATABASE.md](UNIVERSAL_UNIT_CELL_DATABASE.md)) corrects an architectural
assumption rather than a scientific one: Phase 7's 4.8×4.8 mm period was one build's locked choice, not a solver constant
(`periodic_bc.RectangularLattice` already takes any period), so nothing may keep treating it as global. `src/phase8d_rev/` adds a registry
(`UniversalDatabase`) of independent, topology/period/material/frequency/polarization/incidence-tagged `LibraryEntry` objects — read-only wrappers
around already-solved `database.Database`s — plus a physical-filter-first query engine (`search_unit_cells`, never ranking by phase coverage) that
returns either real matches or an explicit `DATABASE_COVERAGE_MISSING` with a machine-readable `LibraryBuildRequest`, never a rescaled or
interpolated substitute. The existing 324-cell/4.8 mm build is registered as exactly one `reference_validation_library` entry — proven, not
assumed, byte-identical before/after this phase's run — and used as the sole benchmark: an exact-match query recovers all 324 geometries with full
cross-polarization data intact (read from the same internal Phase 7 checkpoint Phase 8C uses); a period-mismatched query correctly returns
`DATABASE_COVERAGE_MISSING` rather than silently substituting the 4.8 mm data. No new period, topology or substrate has actually been built or
validated by this phase — the architecture supports them; only 4.8 mm/rectangular-patch/Rogers RT/Duroid 5880 is real. `CONTROL_VARIABLE_UNVALIDATED`
is carried through unchanged (this phase does not touch Phase 8C's finding, by design). **Phase 8E**
(`scripts/research/run_phase8e_control_bridge.py`, [PHASE8E_CONTROL_BRIDGE.md](PHASE8E_CONTROL_BRIDGE.md)) builds the actual reference-wave → object-wave →
required-control-variable → unit-cell-response bridge, framework-first rather than S11-or-S21-first: four candidate formulations (A: Fong et al.
2010 grounded impedance; B/C: the identical phase-difference-hologram kernel claimed against reflection or transmission; D: a new generalized 2x2
transfer-operator candidate) are tested, in a fixed order, against the real measured architecture (A and D are `ARCHITECTURE_INCOMPATIBLE` — no
ground plane; measured cross-polarization <1e-6), a real angular-spectrum reference-plane-propagation implementation (never an arbitrary phase
offset), and a genuine 2-D-complex-plane controllability test (reusing Phase 8C/8D-REV, never phase coverage alone). A real continuous-vs-discrete
beam-reconstruction test (standard aperture Fourier-optics far-field integral) shows the idealized phase-only hologram places all 3 demo beams
exactly on target (0.0° pointing error) while the real, non-unit-magnitude S11xx/S21xx library values degrade that to 4.6–8.1° once actually
matched against real geometries — a genuine, measured consequence of this library's achievable magnitude range, not asserted. Final status:
**`NO_VALIDATED_FORMULATION`** — no candidate's 9-checkbox geometry-inversion gate is fully open (B/C pass reconstruction but the S11-vs-S21
control-variable question itself remains unvalidated, Phase 8B/8C unchanged; local periodicity is always reported unvalidated, no coupled-cell
solver exists). Geometry inversion remains blocked project-wide; no STEP file is generated. **Phase 9**
(`scripts/research/run_phase9_validation.py`, [PHASE9_VALIDATION.md](PHASE9_VALIDATION.md)) resolves the two questions Phase 8E left open using real
evidence rather than declaring them permanently unresolvable. 9A builds the complete solver-vs-HFSS comparison infrastructure but — since no HFSS
export exists anywhere in this repository — honestly reports every comparison as `NOT_AVAILABLE_FROM_CURRENT_DATA` (both S11xx and S21xx classify
`UNVALIDATED`), never a fabricated number. 9B is a genuinely new, real capability: a heterogeneous N×N "supercell" (multiple different patches in
one larger-period unit cell) solved through the *exact same, unmodified* MoM pipeline Phases 1–7 already use — no coupled-cell solver was invented.
A homogeneous 3×3 supercell reproduces the isolated single-cell result to ~1e-14 (a real correctness proof), while a realistic strong-contrast 3×3
neighborhood changes the centre cell's induced current by **99.7%** of its own magnitude — a substantial, measured failure of the periodic-cell
assumption, not asserted. Final status: **`LOCAL_PERIODICITY_INVALIDATED`**; geometry inversion and final STEP generation remain locked. **Phase 10**
(`scripts/research/run_phase10_environment_model.py`, [PHASE10_ENVIRONMENT_MODEL.md](PHASE10_ENVIRONMENT_MODEL.md)) asks the question Phase 9 could not:
how much neighborhood information does the system actually need, and can a finite model replace the isolated-cell lookup? A real measured finding
across all 5 tested centre geometries: diagonal neighbours are not negligible (~25–70% of the cardinal contribution), a cardinal-only "5-cell"
model misses up to 300–490% of the centre cell's own response magnitude under strong contrast, and even the full 3×3 "9-cell" neighbourhood is
`LARGER_NEIGHBOURHOOD_REQUIRED` for 4 of 5 centres. A 640-sample environment-response dataset (real MoM solves, a geometry-based holdout centre,
symmetry-safe train/test grouping) was built to test a ridge-regularized nonlinear surrogate against two required baselines — and, honestly
reported rather than forced to look better: **neither the linear-additive nor the nonlinear surrogate beats the naive, already-invalidated
isolated-cell baseline on blind-test relative error** (0.57 vs. 1.01 vs. 1.35 mean relative error). First-stage inversion (MoM-verified, never
judged by the surrogate alone) and a beam-level far-field comparison were both implemented and run regardless, for completeness. Final status:
**`SURROGATE_INVALIDATED`**; geometry inversion and final STEP generation remain locked, independently, by both the still-missing HFSS control
variable and this new, measured surrogate failure. **Phase 11**
(`scripts/research/run_phase11_heterogeneous_mom.py`, [PHASE11_HETEROGENEOUS_MOM.md](PHASE11_HETEROGENEOUS_MOM.md)) deliberately does NOT try another
data-fitted model. It asks the physics question directly: how far does electromagnetic interaction from a centre cell actually reach, using only
the existing, unmodified `Z I = V` MoM solver as the reference (no ML anywhere in `src/phase11/`, mechanically grep-checked by its own test suite).
A fixed checkerboard pattern is grown from 1×1 through 11×11 (every smaller size an exact centred truncation of the same larger pattern, verified
parity-preserving regardless of offset); for all 5 tested centre geometries the response has NOT settled to within the solver's own measured
self-consistency floor by 9×9/11×11. Phase 10's diagonal-coupling finding is independently reproduced (ratio 0.44, literal reuse of Phase 10's own
functions). The key new finding: two spatial patterns with **identical local statistics** — a distributed checkerboard vs. a single contiguous A|B
block, same two geometries, same 50/50 split — produce substantially different centre responses (a 90% relative difference), even though nearby
pairwise perturbations are themselves close to additive (interaction error ~1e-4 relative). Final classification: **`CASE_D_PATTERN_DEPENDENT`** —
the actual spatial arrangement must be retained; neighbourhood size alone is not a sufficient description. A deterministic (non-ML, isolated-cell-
ranked, direct-MoM-verified) first-stage inversion and real heterogeneous-MoM beam-level/three-beam holographic validation were both implemented and
run to completion regardless. Geometry inversion and final STEP export remain **LOCKED** — independently, by both the still-missing HFSS control
variable and this phase's own measured pattern-dependence/non-convergence finding. **Phase 12A**
(`scripts/research/run_phase12a_ground_plane.py`, [PHASE12A_GROUND_PLANE.md](PHASE12A_GROUND_PLANE.md)) adds an optional **PEC
ground plane** at z = −h (`substrate.backing: pec`; default `air`, every existing config unchanged, verified bit-identical),
so the solver can model the grounded reflectarray/holographic-impedance-surface stack every incompatible formulation in
Phases 8B/8E needed (Fong 2010, Minatti, Pozar 1997, Budhu & Grbic 2021). With a ground, S21 ≡ 0 by construction (no port 2
exists — reading it raises rather than silently returning 0), so the S11-vs-S21 control-variable question Phases 8A-8C left
open is settled by the architecture itself for this stack. The grounded closed forms (derived by substituting the short's
Γ_L = −1 into the existing bottom-load-reflection algebra, same sigma-substitution style as the ungrounded forms) agree with
a newly-added, independently-derived PEC-terminated 3-D field-matching solver (`tests/reference_layered.py`) to ~1e-14 across
normal/oblique incidence, lossy/lossless substrates, multiple periods and thicknesses, and reproduce the lossless
`|Γ|=1`/lossy `|Γ|<1`/`h→0 ⇒ Γ→−1` textbook slab identities to the solver's own floating-point noise floor. A grounded
rectangular-PEC-patch sweep on the *same* Phase 7 substrate/period/mesh (L=W, 18 points, 1.0–4.4 mm) measures **203.6° of
360°** S11xx phase coverage at 12.5 GHz — a dramatic, physically-expected improvement over the ungrounded library's 24.8°,
since a grounded patch can approach an actual resonance the ungrounded scattering patch never reaches in this band — while a
mesh/Floquet refinement study on one representative grounded patch found a larger (~0.048) successive-difference than the
0.02 tolerance the Phase 7 mesh was chosen for, an honestly-reported, measured mesh-sensitivity finding near resonance, not
smoothed over. Final status: **`GROUNDED_STACK_PARTIAL`** (the mesh-convergence finding above; every other check —
bare-slab closed form, TM0 surface-wave-pole avoidance for both the 4.8 mm and 8 mm periods across 10–15 GHz — passed).
Geometry inversion and STEP export remain **LOCKED** — this phase does not attempt to unlock them, independently of its own
result, and HFSS status is unchanged (`HFSS_NOT_AVAILABLE`, still no HFSS export anywhere in this project). **Phase 12B**
(`scripts/research/run_phase12b.py`, [PHASE12B_MESH_PERIOD_STUDY.md](PHASE12B_MESH_PERIOD_STUDY.md)) traces Phase 12A's own
`mesh_error=0.048` finding to its root cause: a mesh cell is metal iff its centroid is metal, and Phase 7/12A's 32-cells-per-4.8mm
mesh does not divide the 0.2mm sweep step evenly, so many feature edges are staircased and — a real, measured consequence —
**155 of the Phase 7 library's 324 geometries are mesh-identical to another one** (169 distinct, 105 duplicate groups, exactly
matching an independent count of equal *solved* S11 values). A new, additive `discretization.py`/`sweep.py` capability
(`smallest_aligned_cells`, `metal_mask_key`/`duplicate_geometry_groups`, `SweepPlan.mesh_duplicate_groups`) detects and avoids
this — every pre-Phase-12B sweep, config and result is verified bit-identical. On the resulting aligned meshes, Phase 12A's
grounded refinement study now **converges monotonically and within tolerance** (successive differences 0.0058, 0.0030 vs. the
unaligned 0.048), confirming the earlier number was a staircasing artefact, and the same representative patch's phase coverage
improves from 203.6° to 237.5° at 12.5GHz. Tracing the same mesh through Phases 9–11's own (much coarser, 0.6mm) solve mesh
finds Phase 10's centers **C2 and C3 are mesh-identical** — a real, previously-unknown finding — though the held-out blind-test
center (C5) has zero leakage, so Phase 10's core conclusions stand. A grounded period study (4.8/8.0/12.0mm, aligned meshes)
measures substantially better coverage at larger periods (**279.2° at 8mm, 296.6° at 12mm**, vs. 237.5° at 4.8mm) — but also a
real trade-off: at a 45° feed angle the 12mm period's grating onset (14.6GHz) falls inside the 10–15GHz band, and its
surface-wave-pole margin shrinks from a full k0 at 10GHz to just 0.18×k0 at 15GHz. No period reaches a full 360° with a
single-layer resonant patch. Final statuses: **`ALIGNED_MESH_VALIDATED`**, **`CONCLUSIONS_AT_RISK`** (Phase 10's C2/C3
degeneracy, contained, not overturning), **`COVERAGE_INSUFFICIENT`** (best 296.6° < the 300° practical threshold). Geometry
inversion and STEP export remain **LOCKED**; HFSS status is unchanged. **Phase 12C**
(`scripts/research/run_phase12c.py`, [PHASE12C_GROUNDED_LOCALITY.md](PHASE12C_GROUNDED_LOCALITY.md)) asks, for Phase 12B's own
favoured grounded design point (P=8mm, RT/Duroid 5880, 12.5GHz, recorded as `config_phase12c_grounded_8mm.yaml`, not
hard-coded), whether local periodicity now holds. A `GroundedSupercellUnitCell` wraps Phase 9's frozen `SupercellUnitCell`
unmodified (`backing="pec"` only) and is validated to ~1e-15 against a homogeneous supercell — after finding and fixing a
genuine floquet-order bug: only `floquet_order = 2×total_mesh_cells` reaches machine precision, a merely rank-sufficient
order silently left ~1% spurious error. At that corrected convention, local periodicity is
**`GROUNDED_LOCAL_PERIODICITY_VALID_FOR_SMOOTH_GRADIENTS`** but only up to a measured **2.7°/cell** neighbour phase step —
every larger tested step (81–107°/cell, all near or at resonance) fails on phase, magnitude, or both, and an exactly
periodic 6-cell phase-gradient supercell (achieved span 260.4°) shows the same breakdown independently. Coupling is
**`GROUNDED_COUPLING_NOT_LOCAL_FOR_TESTED_SIZES`** (n≤5, gradient and checkerboard, off- and near-resonance, none
converged). Re-solving three of Phases 9–11's own benchmark cases read-only on a much finer aligned mesh finds Phase 9's
own original strong-contrast evidence for `LOCAL_PERIODICITY_INVALIDATED` (99.7% relative deviation at 8 cells/subcell)
drops to **21.8%** (below its own 50% substantiality threshold) at 48 cells/subcell — independently confirmed not a floquet
artifact — and Phase 10's C2/C3 also turn out mesh-**distinct** at a finer mesh, while Phase 11's Case D and coupling-radius
findings hold. 5 single-cell + 3 supercell grounded STEP exports (ground solid included) with real MoM references, and
one-port HFSS CSV ingestion (tested on a synthetic file only), are `HFSS_EXPORT_READY`. Geometry inversion and STEP export
remain **LOCKED**; HFSS status is unchanged. **Phase 12D**
(`scripts/research/run_phase12d.py`, [PHASE12D_LOCAL_PERIODICITY.md](PHASE12D_LOCAL_PERIODICITY.md)) re-measures Phase 12C's
step limit with the design-relevant reflection-phase/amplitude metric (not the raw dipole proxy), reusing
`sparameters.scattering_from_fields` unmodified to convert any dipole moment/current into an equivalent S11 (validated
exactly against `compute_s_parameters`). At the largest embedded-solve mesh this project can afford (24 cells/subcell,
measured: 76.8s/196.4s/409.5s at 24/28/32 cells/subcell, not finishing in 14+ minutes at 40), local periodicity is
**`LOCAL_PERIODICITY_STEP_LIMIT_MEASURED`** up to **5.8°/cell** (off-resonance only) — every near/at-resonance step
actually reachable at feasible mesh resolution (41.7–160.3°/cell, non-monotonic in step size) already fails badly
(19–173° phase error), so the true crossover remains bounded, not pinned to a single number. The real decision for
Phase 13 comes from a genuinely new capability (`phase12d/floquet_power.py`): reflected power in *every* Floquet order
of a solved supercell's own lattice, computed directly from the surface current (no dipole proxy), validated to machine
precision on a full-PEC-sheet trivial case and on a uniform-array reproduction of |S11|². Applied to real, full periodic
phase-gradient supercells (N=6, 8, 12 cells/cycle, cells from a real isolated-cell lookup table, all three sizes
cheap — 7–125s), the result is **`LOOKUP_TABLE_INADEQUATE`**: no tested N puts more than 53% of reflected power in the
intended order (29.5–47.5% measured), with 30–38% consistently lost to the *backward* (−1,0) order — because this
design's own achievable phase span (186° of 360°) is the same coverage shortfall Phases 12A–12C already measured from
a different angle. Real HFSS Floquet-port data (18 modes/port, PEC ground, HFSS 2025 R1) for one case (G1, L=2.0mm) is
ingested via a new general Touchstone reader (`phase9/hfss_reference.py`'s `read_touchstone`) and agrees with the MoM
reference to a mean phase difference of **−0.29°** and \|S11\| to ~1e-5 — **`HFSS_SINGLE_CELL_VALIDATED_PARTIAL`** (the
resonance case G4 is explicitly still pending; no supercell HFSS data exists). Geometry inversion and STEP export
remain **LOCKED**. Updated recommendation for Phase 13 (replacing Phase 12C's own): a full-aperture integral-equation
solver (à la Budhu & Grbic 2021) is the more broadly defensible default for this design point; a lookup table remains
viable only for an aperture confined to the flat, off-resonance part of the curve (steps below ~6°/cell).

## Quick start

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python -m pytest -q
.venv/bin/python examples/print_config.py
.venv/bin/python examples/run_sweep.py            # sweep + database in results/sweep_example (re-run = resume)

# Optional (Phase 6/6A/7, needs a full OpenCASCADE build — pip install -r requirements-cad.txt):
.venv/bin/python scripts/research/export_validation_cases.py   # STEP + results + plots + README per case, results/validation_cases
.venv/bin/python scripts/research/export_hfss_comparison.py     # Phase 6A: exact overlay-ready plots/CSV for case01, from its mom_results.json
.venv/bin/python examples/build_unit_cell_library.py    # Phase 7: 324-geometry x 201-frequency production library (long-running, resumable)
.venv/bin/python examples/inspect_unit_cell_library.py  # query an already-built library
.venv/bin/python scripts/research/run_phase12a_ground_plane.py  # Phase 12A: grounded (backing=pec) physics-sanity sweep, results/phase12a_ground_plane
.venv/bin/python scripts/research/run_phase12b.py               # Phase 12B: mesh alignment/duplicate audit/period study, results/phase12b_mesh_period_study
.venv/bin/python scripts/research/run_phase12c.py               # Phase 12C: grounded supercells, local-periodicity/coupling/spot-check, HFSS exports, results/phase12c_*
.venv/bin/python scripts/research/run_phase12d.py               # Phase 12D: step-limit re-measurement, aperture Floquet-power test, real HFSS ingestion, results/phase12d_*
```

Conventions (time factor `exp(+jωt)`, SI units, coordinates, kz branch) are in
[CONVENTIONS.md](../theory/CONVENTIONS.md). Read it before touching any module.

## Layout

| Path | Status | Contents |
|---|---|---|
| `src/constants.py` | done | `C0`, `MU0`, `EPS0`, `ETA0`, time-convention tag |
| `src/units.py` | done | SI ↔ mm/GHz/deg; config quantity parsing |
| `src/frequency.py` | done | frequency validation, ω, k0, λ0 |
| `src/config.py` | done | strict YAML loading → `SimulationConfig` |
| `src/materials.py` | done | `Substrate` (ε_c = ε_r(1 − j tanδ), k), `Metal` (PEC only) |
| `src/geometry.py` | done | `Geometry` (rectangle/patch/dipole/slot, rotation), `UnitCell` (**+ `backing`: `"air"`\|`"pec"`, Phase 12A**), `is_metal` |
| `src/periodic_bc.py` | done | `RectangularLattice`, Floquet wavenumbers, `normal_wavenumber` (kz branch) |
| `src/excitation.py` | done | `PlaneWaveExcitation`; `stack_tangential_field` (unpatterned-stack E at z = 0); `ReferencePlanes` (explicit S-parameter planes; **`bottom_air_thickness=None` = grounded, no port 2, Phase 12A**) |
| `src/green_function.py` | done | `LayeredGreenFunction`: TE/TM equivalent impedances, 2×2 dyad `G` (`Ẽ_t = G J̃`), Floquet-harmonic dyads, exact `kt = 0` handling; **tangential E/H at arbitrary z** (`sheet_fields`, `stack_fields`); **optional PEC ground plane at z = −h** (`backing="pec"`, default `"air"` unchanged; Phase 12A, [PHASE12A_GROUND_PLANE.md](PHASE12A_GROUND_PLANE.md)) |
| `src/basis_functions.py` | done | `MeshSettings`, `RooftopBasis`: x/y rooftops on a rectangular mesh, periodic wrap, closed-form transforms, pointwise evaluation |
| `src/mom_matrix.py` | done | `assemble_impedance_matrix` (Z), `assemble_excitation_vector` (V) |
| `src/solver.py` | done | `assemble_system` / `build_system` (geometry → mesh → V → Z), `solve` (one or several right-hand sides) |
| `src/sparameters.py` | done | `SParameters` (2×2 blocks S11/S12/S21/S22, co/cross-pol accessors, flat record; **`grounded=True` one-port result for a grounded stack — S12/S21/S22 undefined, reading them raises, Phase 12A**), `compute_s_parameters`, `stack_s_parameters` (bare slab), `at_planes` (explicit reference-plane phase), `decompose_air_fields`, `surface_wave_proximity` (Phase 12A TM0 pole diagnostic) |
| `src/sweep.py` | done | `FrequencyGrid`, `ParameterAxis`, `SweepDefinition`, policies; `plan_sweep` expands and validates every cell before any solve; `cell_id` (**includes `backing`, Phase 12A**); **+ `SweepPlan.mesh_duplicate_groups`/`n_mesh_distinct`, `plan_sweep(..., warn_duplicates=True)`, Phase 12B — additive, `run_sweep` suppresses the warning so no pre-Phase-12B sweep changes** |
| `src/discretization.py` | done | a priori mesh/Floquet checks (area error, cells/wavelength, grating onset — **+ oblique-incidence grating onset, Phase 12B**) and the refinement study (one- or two-port, Phase 12A); **+ mesh-alignment (`smallest_aligned_cells`) and duplicate-geometry detection (`metal_mask_key`/`duplicate_geometry_groups`), Phase 12B** |
| `src/runner.py` | done | `run_sweep`: solve, append-only checkpoint (`checkpoint.jsonl`), resume/restart, database assembly (one-port/grounded sweeps write JSON only, no CSV, Phase 12A) |
| `src/database.py` | done | `Database`, `Record` (**+ `port_topology`: `"two_port"`\|`"one_port"`, Phase 12A**), `CellEntry`; JSON (bit-exact) and CSV export/import with validation (CSV rejects one-port records); queries |
| `src/polarization.py` | done | circular basis, linear ↔ circular conversion, co/cross power, handedness under `exp(+jωt)` |
| `src/phase.py` | done | wrap to (−180, 180], unwrap with step diagnostics |
| `src/postprocess.py` | **empty** | docstring only |
| `src/cad_export.py` | done | STEP export of `UnitCell` solids (substrate + metal, boolean-cut slot for aperture polarity), `BRepCheck_Analyzer` validity check, STEP round-trip re-import/summarize; needs `requirements-cad.txt`; **+ `build_ground_solid`, a third `Ground` part flush below the substrate for `backing="pec"` (Phase 12A)**; **+ `export_supercell_step`/`build_supercell_solids`/`build_supercell_ground_solid`, one metal solid per sub-cell patch of a `SupercellUnitCell` plus an optional shared ground solid (Phase 12C)** |
| `src/validation_cases.py` | done | the 10 representative validation cases (`CASES`), built from unmodified `config.config_from_dict`/`sweep.FrequencyGrid` |
| `src/validation_export.py` | done | per-case folder orchestration: solve (Phase 5, unmodified) → STEP → `geometry.json`/`solver_config.json` → `mom_results.{json,csv}` → plots → `README.md`; see [VALIDATION_EXPORT.md](VALIDATION_EXPORT.md) |
| `src/hfss_comparison.py` | done | Phase 6A: reads an already-exported case's `mom_results.json` only (never re-solves) and writes 5 overlay-ready plots, `hfss_comparison_reference.csv` and `HFSS_COMPARISON.md` |
| `src/library_schema.py` | done | Phase 7: the locked 324-cell library configuration (`LibraryConfig`), CSV/record schema, documented classification thresholds |
| `src/library_builder.py` | done | Phase 7: wires the L/W × frequency sweep into `runner.run_sweep` (Phase 5, unmodified); resumable, deterministic |
| `src/phase_extraction.py` | done | Phase 7: exact-grid design-frequency extraction, the phase-vs-geometry table, phase-coverage/largest-gap and phase-collision analysis |
| `src/library_classifier.py` | done | Phase 7: documented resonance detection (peak prominence), grating-risk (reuses `discretization.py`) and phase-usability classification |
| `src/geometry_export.py` | done | Phase 7: one STEP file per solved geometry (`cad_export.py`, Phase 6, unmodified) with an independent bounds re-check |
| `src/library_plots.py` | done | Phase 7: the 8 required library plots from already-computed data |
| `src/library_validation.py` | done | Phase 7: deterministic HFSS validation subset selection and JSON/CSV internal-consistency checking; makes no HFSS-agreement claim |
| `src/unit_cell_library.py` | done | Phase 7: `unit_cell_library.{json,csv}`, `manifest.json`, `LIBRARY_REPORT.md`, and the read-only query/phase-inversion API (`UnitCellLibrary`) |
| `src/phase12a/` | done | Phase 12A: orchestration/reporting only (`sweep_study.py`, `analytic_checks.py`, `gate.py`, `plots.py`, `report.py`, `specification.py`, `run.py`); every physics change lives in the core modules above |
| `src/phase12b/` | done | Phase 12B: orchestration/reporting only (`reassessment.py` — Task B; `duplicate_audit.py` — Task C, read-only; `period_study.py` — Task D; `gate.py`, `plots.py`, `report.py`, `specification.py`, `run.py`); Task A's actual capability lives in `discretization.py`/`sweep.py` above |
| `src/phase12c/` | done | Phase 12C: `design_point.py` (loads `config_phase12c_grounded_8mm.yaml`); `grounded_supercell.py` (`GroundedSupercellUnitCell` wrapper — Task A); `locality_study.py` (Task B); `coupling_study.py` (reuses `phase11` functions unmodified — Task C); `spot_check.py` (Task D, read-only); `exports.py` (Task E); `gate.py`, `specification.py`, `plots.py`, `report.py`, `run.py`; the actual grounded-supercell/STEP-export/HFSS-ingestion capabilities live in `cad_export.py`/`phase9/hfss_reference.py` above |
| `src/phase12d/` | done | Phase 12D: `reflection_metric.py` (dipole/current → equivalent-S11 conversion, reuses `sparameters.scattering_from_fields` unmodified — Task A); `step_curve.py` (Task A orchestration, reuses `phase12c.locality_study`'s grid builders unmodified); `floquet_power.py` (per-Floquet-order reflected power from a solved current — Task B's real new capability); `aperture_test.py` (Task B orchestration); `hfss_ingest.py` (Task C, uses the new `phase9.hfss_reference.read_touchstone`); `gate.py`, `specification.py`, `plots.py`, `report.py`, `run.py` |
| `tests/` | done | 2065 tests (2061 run + 4 skipped without the optional `requirements-cad.txt`), mathematical/numerical properties only; full suite passes in 1955s (32m35s); `tests/reference_layered.py` is an independent field-matching solver (sheet and incident-wave fields at any z, plus a PEC-terminated grounded variant, Phase 12A) used only by tests |
| `examples/print_config.py` | done | loads `config.yaml`, prints in mm/GHz |
| `examples/sweep_config.yaml`, `examples/run_sweep.py` | done | example `sweep:` section; runs it (resumable) and prints a summary |
| `scripts/research/export_validation_cases.py` | done | runs `validation_export.export_all_cases` for some or all of the 10 cases, prints progress and a summary |
| `examples/build_unit_cell_library.py` | done | Phase 7: builds the full 324-geometry × 201-frequency production library into `results/unit_cell_library/` (long-running, resumable) |
| `examples/inspect_unit_cell_library.py` | done | Phase 7: read-only query/phase-inversion demo against an already-built library |
| `scripts/research/export_hfss_comparison.py` | done | runs `hfss_comparison.export_hfss_comparison` for one already-exported case (default `case01_non_resonant_patch`) |
| `scripts/research/run_phase12a_ground_plane.py` | done | Phase 12A: grounded rectangular-PEC-patch physics-sanity sweep + validation into `results/phase12a_ground_plane/` |
| `scripts/research/run_phase12b.py` | done | Phase 12B: mesh alignment, duplicate-geometry audit, grounded period study into `results/phase12b_mesh_period_study/` |
| `scripts/research/run_phase12c.py` | done | Phase 12C: grounded supercells, local-periodicity/coupling/fine-mesh spot-check, HFSS-ready exports into `results/phase12c_grounded_locality/`, `results/phase12c_hfss_export/` |
| `scripts/research/run_phase12d.py` | done | Phase 12D: design-relevant step-limit re-measurement, aperture-level Floquet-order-power test, real HFSS Touchstone ingestion into `results/phase12d_local_periodicity/` |

Implementation notes:
- Import as `from src.xxx import ...` (run pytest from the project root; `pytest.ini` sets `pythonpath = .`).
- No electrical value has a code default. `eps_r`, `tan_delta`, `thickness`,
  periods and frequency must be in the config; unknown keys are errors.
- `geometry.type: patch` is mathematically identical to `rectangle`; `dipole`
  and `slot` additionally require `length >= width`. `slot` = hole in a solid
  metal cell. Features must lie inside the cell (touching the edge is allowed).
- PyYAML reads `8e-3` as a string; the loader treats bare-number strings as SI.
- `mesh` (cells per period, Floquet order) and `reference_planes` (thicknesses of the modelled air regions; planes at
  `z_top = +top`, `z_bottom = −(h + bottom)`) are **required** config sections. Reference planes do not affect `Z` or `V`.
- Mesh rule: `2·floquet_order + 1 ≥ max(cells_x, cells_y)` (else `Z` is rank-deficient; rejected). Accuracy needs a larger order
  (≈ 2–4 × cells). Align `cells` with feature dimensions; other edges are staircased.
- Quick use: `from src.config import load_config; from src.sparameters import s_parameters_from_config;
  s = s_parameters_from_config(load_config("config.yaml")); s.component("S21_xx"), s.phase_deg("S21_xx"), s.to_record()`.
  (Lower level: `build_system(cfg)` → `Z, V`; `solve(Z, V)` → rooftop currents.)
- S-parameters: `S_ij[received, incident]`, index 0 = x, 1 = y, same global axes at both ports; the incident polarization of the config
  is ignored (both x and y are solved). Reference planes: `at_planes` re-references analytically.

## Sweeps and the database (Phase 5) — summary

Add a `sweep:` section (see `examples/sweep_config.yaml`, full reference in [DATABASE.md](../theory/DATABASE.md) Sec. 9); `run_sweep(config, out_dir)`
returns a `SweepResult` and writes `checkpoint.jsonl`, `database.json`, `database.csv`.

- Frequency grid `f_start → f_stop → step`, both end points included; a span that is not a whole number of steps is rejected.
- Sweepable: `length, width, center_x, center_y, rotation_deg` (Cartesian product, first axis slowest). Lattice, substrate, mesh,
  planes and geometry type are fixed per database. Every cell is validated **before** the first solve; `on_invalid: skip` records
  invalid cells instead of stopping.
- Validation of the discretization is *reported*, never silently adapted: a priori checks per record; opt-in refinement study
  (`floquet_error`, `mesh_error` are successive differences, not error bars).
- Resume: results are keyed by everything that determines them, so interrupted runs continue, sweeps can be extended, and stale
  results are never reused. Outputs are byte-identical for identical inputs (same platform); no timestamps.
- Phase is stored on (−180, 180]; unwrapping needs < 180° true change per sample (checked and warned at 90°).
- Circular data are derived from the linear S: `S21_RR` = R received / R incident, labels of global x–y vectors. **Handedness:** under
  `exp(+jωt)` `e_R` rotates clockwise seen from +z; it is RHCP for a wave travelling −z and LHCP for +z (right-hand rule), so a
  reflection `S11_RR` reverses handedness. See DATABASE.md Sec. 7 for the derivation and identities.

## Validation unit-cell export (Phase 6) — summary

`scripts/research/export_validation_cases.py` (needs `pip install -r requirements-cad.txt`) writes one folder per representative unit cell under
`results/validation_cases/<case_id>/`: `unit_cell.step`, `geometry.json`, `solver_config.json`, `mom_results.{json,csv}`, `plots/*.png`,
`README.md`. Full reference, coordinate/unit mapping and the case catalog rationale: [VALIDATION_EXPORT.md](VALIDATION_EXPORT.md).

- STEP solids are built directly from the same `UnitCell`/`Geometry`/`Substrate` objects the solver consumes (`src/cad_export.py`) — no
  parameter is re-entered; geometry-consistency tests independently re-import every exported file and check its recovered bounding
  box/volume/centroid against closed-form and `UnitCell`-derived expectations, and every solid's topology against OpenCASCADE's own
  `BRepCheck_Analyzer`.
- 10 cases (`src/validation_cases.py`) cover non-resonant, resonant, narrow, wide, four distinct length/width pairs, rotated (with and
  without an offset), and slot/complement (with and without rotation), chosen by measuring the existing solver's own output, not assumed.
- This phase performs **no automatic comparison against HFSS or any other tool** and applies **no tuning** of the MoM results — every
  README's HFSS reproduction recipe is manual, by design, and says so explicitly.
- Adds no EM code: `green_function.py`, `mom_matrix.py`, `solver.py`, `sparameters.py` and every number they produce are unchanged.

## Production unit-cell library (Phase 7) — summary

`examples/build_unit_cell_library.py` (needs `requirements-cad.txt`) builds a 324-geometry (18×18 rectangular-patch L/W sweep, 1.0–4.4 mm,
0.2 mm step) × 201-frequency (10–15 GHz, 25 MHz step) library into `results/unit_cell_library/`, reusing `runner.run_sweep` (Phase 5) directly.
Full reference: [UNIT_CELL_LIBRARY.md](UNIT_CELL_LIBRARY.md); build-specific numbers: `results/unit_cell_library/reports/LIBRARY_REPORT.md`.

- Mesh (32×32 cells, Floquet order 32) was chosen by directly measuring build cost and mesh-quality trade-offs, not guessed — a seemingly-ideal
  24×24 mesh (matching the 0.2 mm sweep step) was rejected after it was found to misclassify an entire mesh row for ~200/324 geometries due to
  floating-point rounding in the unmodified mesh boundary code (up to 44% area error, `tests/test_library_builder.py`).
- Phase coverage/largest-gap and resonance are *computed* from the actual solved data (never asserted): the built library measures S11xx phase
  coverage of 24.8° of 360° at 12.5 GHz (largest gap 335.2°) and 0 resonant cells in the 10–15 GHz band (the true resonance for these patches, if
  any, lies above the locked band, similar to Phase 6's dipole case needing 16.4 GHz to show its peak).
- The read-only query API (`unit_cell_library.UnitCellLibrary`) never interpolates or extrapolates: every lookup is an exact match against a
  geometry/frequency actually simulated, or an explicit `NO_VALID_CANDIDATE`/error.
- This phase performs **no automatic comparison against HFSS** and applies **no tuning** — `HFSS_VALIDATION_SUBSET.csv` is a deterministic
  recommendation only.
- Adds no EM code, same as Phase 6.

## Before Phase 13 (real HFSS data + a validated large-area/reduced-order heterogeneous formulation + GUI)

### Assumptions (physical model)
1. Rectangular lattice; one unit cell; rectangular features only (x/y-aligned meshes; others staircased).
2. Stack: air | zero-thickness PEC pattern at z = 0 | one substrate layer | air, **or** (Phase 12A, opt-in via
   `substrate.backing: pec`) air | pattern | substrate | PEC ground at z = −h. No superstrate either way.
3. Substrate homogeneous, isotropic, non-magnetic, constant ε_r and tanδ. Metal is PEC, infinitely thin. Slots = metal-complement currents.
4. **Normal incidence**, specular (0,0) order only as port mode; x/y linear (and derived circular) polarization basis, same global axes at both ports.
5. Ports = the two explicit reference planes at the outer boundaries of the modelled air regions.
6. Period < wavelength so that only (0,0) propagates (flagged per record; a warning in single solves).

### Resolved
Kernel, excitation, MoM (Phases 2–3); fields at arbitrary z, S-parameter and reference-plane conventions, validation against Fresnel slab,
full PEC sheet, energy conservation, reciprocity (Phase 4); sweeps, discretization validation, resume, database, circular basis,
phase wrapping (Phase 5); physical STEP export of representative unit cells with matching MoM reference results, plots and reproduction
instructions, for manual external comparison (Phase 6); a production rectangular-patch unit-cell library with measured phase coverage,
STEP geometry for every cell, and a read-only query/phase-inversion API (Phase 7); a scientific feasibility gate that verifies library
integrity/coverage/resonance/grating risk from the real data, refuses to validate the S11xx-as-holographic-phase assumption without an actual
formulation to check it against, and never invents a required-phase distribution in its absence (Phase 8A); object-wave/multi-beam field
equations and an honest, literature-cited investigation of two candidate holographic formulations against the real library, finding one
structurally incompatible (no ground plane, no surface-wave solver) and the other still open on the same S11-vs-S21 question (Phase 8B);
full-library S11xx *and* S21xx scattering analysis (magnitude, complex-plane, phase coverage, phase-magnitude coupling) read from the still-present
Phase 7 solve checkpoint's full linear-basis S-matrix, plus real numeric reciprocity (~1e-15) and L↔W mirror-symmetry (~1e-15) confirmation on the
solved data, concluding `CONTROL_VARIABLE_UNVALIDATED` without ever selecting a control variable by phase coverage (Phase 8C); a universal,
period-agnostic, topology-aware unit-cell database registry (`UniversalDatabase`/`LibraryEntry`/`search_unit_cells`) that wraps the existing
324-cell/4.8 mm build as one benchmark entry among possibly many, provably without modifying it, and correctly refuses a period-mismatched query
instead of ever rescaling or interpolating a substitute (Phase 8D-REV); a framework-first (never S11-or-S21-first) evaluation of four candidate
holographic formulations against the real measured architecture, a real angular-spectrum reference-plane propagation, a genuine 2-D-complex-plane
controllability test, and a continuous-vs-discrete beam-reconstruction check (standard aperture Fourier-optics far-field integral) showing the
idealized hologram reconstructs all 3 demo beams exactly while the real, achievable library responses measurably degrade that reconstruction
(Phase 8E); a real solver-vs-HFSS comparison infrastructure that honestly reports every metric `NOT_AVAILABLE_FROM_CURRENT_DATA` rather than
fabricate HFSS results that do not exist, and a genuinely new heterogeneous-supercell capability (`SupercellUnitCell`, duck-typed against the exact
same, unmodified MoM pipeline Phases 1-7 already use) proving, to ~1e-14, that a homogeneous supercell reproduces the isolated-cell result exactly,
then measuring a real 99.7%-of-signal-magnitude deviation for a realistic strong-contrast neighborhood — local periodicity `INVALIDATED` (Phase 9);
a real neighborhood-sensitivity study across 5 center geometries measuring cardinal, diagonal and distance-2 contributions independently (never
assuming rotational symmetry or that a smaller neighborhood is sufficient because it is cheaper), finding diagonal contributions are not negligible
and a 3×3 ("9-cell") neighborhood is itself `LARGER_NEIGHBOURHOOD_REQUIRED` for 4 of 5 centers; a 640-sample real-MoM environment-response dataset
with a geometry-based holdout center and symmetry-safe (never leaking a mirrored environment across train/test) grouping; three baselines
(isolated-cell, linear-additive, ridge-regularized nonlinear surrogate) honestly compared on a blind test, with neither fitted model beating the
naive isolated-cell baseline; and first-stage MoM-verified inverse geometry search plus a beam-level far-field comparison, both implemented and run
to completion regardless of the forward-model result (Phase 10); a physics-only (zero machine learning, mechanically grep-verified) coupling-radius
study growing one fixed checkerboard pattern from 1×1 through 11×11 via exact centred truncation for all 5 center geometries, none of which converges
to the solver's own self-consistency floor by 9×9/11×11; an independent reproduction of Phase 10's diagonal-coupling finding via literal reuse of its
own functions; and the discovery that two spatial patterns with identical local statistics (a distributed checkerboard vs. a contiguous A|B block)
produce a 90%-relative-difference center response despite near-neighbor pairwise interactions themselves being close to additive — `CASE_D_
PATTERN_DEPENDENT`, spatial arrangement itself must be retained, not just neighborhood size (Phase 11); and an optional PEC ground plane at
z = −h (`backing="pec"`, default `"air"` unchanged and verified bit-identical), independently derived and verified to ~1e-14 against a new
PEC-terminated 3-D field-matching solver, one-port S-parameters (S12/S21/S22 undefined, not silently zero), a one-port database record schema,
and a grounded-stack STEP export — measuring a large phase-coverage improvement (203.6° vs. 24.8° ungrounded) from a real, resonance-driven effect
(Phase 12A); and a mesh-alignment/duplicate-geometry-detection capability (`discretization.smallest_aligned_cells`,
`metal_mask_key`/`duplicate_geometry_groups`, `sweep.SweepPlan.mesh_duplicate_groups`, purely additive, every pre-Phase-12B sweep
verified bit-identical) that traced Phase 12A's own unaligned-mesh `mesh_error=0.048` to a genuine staircasing artefact — the
aligned-mesh refinement study now converges monotonically and within tolerance — and quantified real, previously-undocumented
duplicate-geometry effects in the Phase 7 library (169/324 distinct) and in Phase 10's own 5 center geometries (C2≡C3 at the
solver's own mesh, contained: no train/blind-test leakage), plus a real, measured demonstration that a larger grounded period
substantially improves achievable phase coverage (up to 296.6° at 12 mm) at the cost of a narrower oblique-incidence
grating-free/surface-wave-safe margin (Phase 12B); and, for Phase 12B's own favoured grounded, P=8 mm design point, a
`GroundedSupercellUnitCell` wrapper around Phase 9's frozen `SupercellUnitCell` (duck-typed, unmodified inner class) validated
to ~1e-15 against a homogeneous supercell once a genuine floquet-order bug was found and fixed (only `floquet_order =
2×total_mesh_cells` reaches machine precision; a merely rank-sufficient order left ~1% spurious error) — then, at that
corrected convention, `GROUNDED_LOCAL_PERIODICITY_VALID_FOR_SMOOTH_GRADIENTS` up to a measured 2.7°/cell neighbour step
(every larger tested step, all near/at resonance, fails), `GROUNDED_COUPLING_NOT_LOCAL_FOR_TESTED_SIZES` (n≤5, both gradient
and checkerboard patterns, off- and near-resonance), and — re-solving three of Phases 9–11's own read-only benchmark cases on
a much finer aligned mesh without touching their code or results — the discovery that Phase 9's own original strong-contrast
`LOCAL_PERIODICITY_INVALIDATED` evidence (99.7% relative deviation at an 8-cell/subcell mesh) drops to 21.8% (below its own
50% substantiality threshold) at a finer, aligned 48-cell/subcell mesh, independently confirmed converged and not a floquet
artifact, while Phase 11's Case D and coupling-radius findings hold at the finer mesh — plus 5 single-cell and 3 supercell
grounded STEP exports (ground solid included) with real MoM reference results, and one-port HFSS CSV ingestion support
verified on a synthetic file only (Phase 12C); and, re-measuring Phase 12C's own step limit with the design-relevant
reflection-phase/amplitude metric instead of the raw dipole proxy (a new `sparameters.scattering_from_fields`-based
conversion, validated to machine precision), `LOCAL_PERIODICITY_STEP_LIMIT_MEASURED` up to a real 5.8°/cell off-resonance
(every near/at-resonance step reachable at the largest embedded-solve mesh this project can afford, 24 cells/subcell —
finer meshes measured infeasible, up to 409.5s at 32 cells/subcell and not finishing in 14+ minutes at 40 — already
fails, 19–173° phase error, non-monotonically in step size); a genuinely new per-Floquet-order reflected-power capability
(`phase12d/floquet_power.py`, validated to machine precision on a full-PEC-sheet trivial case and a uniform-array
reproduction of \|S11\|²) applied to real periodic phase-gradient supercells (N=6/8/12 cells per cycle, all three
cheap — 7–125s), finding `LOOKUP_TABLE_INADEQUATE`: no tested N puts more than 53% of reflected power in the intended
order (29.5–47.5% measured, 30–38% consistently lost to the backward order) because this design's own achievable phase
span (186° of 360°) is the same coverage shortfall Phases 12A–12C already measured; and real HFSS Floquet-port data
(18 modes/port, PEC ground, HFSS 2025 R1) for one case (G1, L=2.0mm), ingested via a new general Touchstone reader
(`phase9.hfss_reference.read_touchstone`), agreeing with the MoM reference to a mean phase difference of −0.29° and
\|S11\| to ~1e-5 — `HFSS_SINGLE_CELL_VALIDATED_PARTIAL` (Phase 12D).

### Still open
- External validation (below, now exportable via Phase 6/7 but not yet performed); Babinet check.
- The Phase 7 library's narrow measured phase coverage (24.8° of 360° at 12.5 GHz): a holographic generator needs much broader coverage, which
  this rectangular-patch-only, ungrounded, single-band library does not yet provide — a different unit-cell topology, a ground plane, and/or a
  wider/different frequency choice are open design questions, not something this phase invents an answer to.
- Which S-parameter (S11 vs. S21) is the physically correct "holographic phase" once the generator's reflective-vs-transmissive architecture is
  decided (Phase 7's `phase_source` is currently spec-mandated to S11, not independently physically validated — `UNIT_CELL_LIBRARY.md` Sec. 4).
  Phase 8A's feasibility gate formalizes this as `PHASE_SOURCE_UNVALIDATED`; Phase 8B's `control_variable.py` investigates both S11xx and S21xx
  in full and reaches the same open conclusion (`INSUFFICIENT_EVIDENCE` for both); Phase 8C measures both across the full 324-cell library and
  confirms the conclusion still holds with real numbers — S21xx actually has *more* phase coverage (43.85° vs. 24.84°), which is reported but
  explicitly not used to pick a winner (`PHASE8C_SCATTERING.md`'s explicit rule against selecting by phase coverage). This is the one remaining
  blocker for the phase-difference-hologram candidate formulation.
- The classic monopole-fed, grounded holographic-impedance-surface formulation (Fong et al. 2010) is `LIBRARY_CONTROL_VARIABLE_INCOMPATIBLE`
  with this project's library — not an open question but a structural finding (no ground plane; no surface-wave dispersion solver). A
  reference-grid → MTS-unit-cell-centre spatial mapping method remains undefined either way (`feasibility_spatial.py`).
- Only one period (4.8 mm), one topology (rectangular patch) and one substrate are actually registered in the Phase 8D-REV `UniversalDatabase` —
  the architecture (`src/phase8d_rev/`) supports more, but a real second library (different period, topology, or substrate) has not been built or
  solved; `candidate_period.py`'s Nyquist-sampling heuristic proposes *candidate* periods only, never a `VALIDATED` one, until an actual dataset
  for it exists.
- No HFSS export (or any other independent full-wave/measured data) exists anywhere in this project. Phase 9A built the complete comparison
  infrastructure (`src/phase9/hfss_reference.py`, `error_metrics.py`, a 10-cell x 7-frequency reference extraction from the real Phase 7 checkpoint)
  and it is fully unit-tested against synthetic data, but every production classification is `UNVALIDATED` for the honest reason
  `NOT_AVAILABLE_FROM_CURRENT_DATA` — this is now the literal, primary blocker on the S11-vs-S21 question, not a numerical-agreement problem.
- The local-periodic approximation (an isolated-cell periodic solve standing in for that cell's response once embedded in a real, heterogeneous
  metasurface) is **`LOCAL_PERIODICITY_INVALIDATED`** (Phase 9B, `src/phase9/local_periodicity.py`): a realistic strong-contrast 3x3 neighborhood
  (this library's own sweep extremes, 1.0/4.4 mm) changes the centre cell's induced-current dipole moment by 99.7% of its own magnitude, measured
  via a real heterogeneous-supercell MoM solve (not a coupled-cell solver — the same, unmodified periodic MoM pipeline, applied to a larger,
  multi-patch unit cell). Per Phase 9's own stop condition, geometry inversion must not proceed on the existing isolated-cell lookup table; a
  coupled-cell or heterogeneous-response library is the required next capability (`PHASE9_VALIDATION.md`).
- No HFSS export (or any other independent full-wave/measured data) still exists anywhere in this project as of Phase 10 — this is now Gate 10A's
  own, permanent blocker, unchanged across two phases: `hfss_or_measurement_reference_exists` cannot be `True` for a project with no HFSS project.
  Geometry inversion cannot unlock until this changes, independent of anything else Phase 10 measured.
- The environment-aware surrogate (`src/phase10/surrogate.py`, a ridge-regularized nonlinear model over the 9-cell neighborhood) is
  **`SURROGATE_INVALIDATED`**: on a deliberately hard blind test (an entire held-out center geometry plus symmetry-safe grouping), it does not
  outperform the naive, already-`LOCAL_PERIODICITY_INVALIDATED` isolated-cell lookup on relative error (1.35 vs. 0.57 mean relative error). It also
  degrades sharply right at an abrupt geometry discontinuity, exactly where a real hologram changes fastest. A richer dataset (more samples, more
  centers, an actually-sufficient neighborhood size per 10A-10C's own finding) and/or a different regression architecture are the required next
  steps — not a larger STEP-generation effort on the existing surrogate.
- Even the full 3×3 ("9-cell") neighborhood representation itself is `LARGER_NEIGHBOURHOOD_REQUIRED` (not merely "9-cell is nice to have") for 4 of
  5 tested center geometries (`neighborhood_study/neighborhood_decision.json`) — a real, measured finding, not a modeling shortcut; a genuinely
  larger (5×5 or adaptive) environment representation is an open question for whatever forward model eventually replaces this phase's surrogate.
- No HFSS export still exists anywhere in this project as of Phase 11 either — the same, unchanged, now three-phases-persistent blocker
  (`physical_control_variable_validated` cannot be `True`); Phase 11's own 14-condition geometry-inversion gate is CLOSED on this basis alone,
  independent of everything else it measured.
- Phase 11's real, physics-only (no ML) coupling-radius study found the metasurface `CASE_D_PATTERN_DEPENDENT`: none of the 5 tested center
  geometries converges to the solver's own self-consistency floor by 9×9/11×11 (3 of 5 are `STILL_SUBSTANTIALLY_CHANGING`, the other 2
  `NOT_CONVERGED`), and two arrangements with identical local statistics (a checkerboard vs. a contiguous block of the same two geometries) produce
  a 90%-relative-difference center response — meaning neighborhood size alone is not a sufficient description of the coupling, and no computational
  reduction (hierarchical matrices, domain decomposition, etc.) was attempted (Sec. 19's own requirement that one first be justified by finding the
  coupling local, which it was not). Isolated-cell geometry inversion, and any inversion based only on a local neighborhood size, is not scientifically
  justified by this project's own data (`PHASE11_HETEROGENEOUS_MOM.md`).
- Phase 12A's grounded-stack physics-sanity sweep found a larger mesh/Floquet successive-difference (~0.048, refinement study on one
  representative grounded patch) than the 0.02 tolerance the reused Phase 7 mesh (32×32 cells, floquet order 32) was chosen for — grounded,
  near-resonant patches are measurably more mesh-sensitive than the ungrounded scattering patches that mesh/sweep-step combination was designed
  around. Final status `GROUNDED_STACK_PARTIAL`, not `GROUNDED_STACK_VALIDATED` (`PHASE12A_GROUND_PLANE.md`); a finer production mesh for grounded
  designs is the required next step, not this phase's to build.
- Even with a ground plane, measured S11xx phase coverage (203.6° of 360° at 12.5 GHz, the same 18-point L=W diagonal Phase 7 swept) is still short
  of full 360° coverage — a thicker/thinner substrate or a different topology is the next open design question, not something Phase 12A answers
  (its own brief's explicit instruction).
- Phase 12A's `GROUNDED_STACK_PARTIAL` mesh-convergence finding was confirmed, by Phase 12B, to be a staircasing artefact of the *unaligned* mesh,
  not genuine under-resolution: on the aligned mesh the same refinement study converges within tolerance (`ALIGNED_MESH_VALIDATED`). Coverage on
  the aligned 4.8 mm mesh improves to 237.5° (still short of 360°); larger periods do measurably better (279.2° at 8 mm, 296.6° at 12 mm,
  `COVERAGE_INSUFFICIENT` against the 300° practical threshold either way) but trade away oblique-incidence/surface-wave margin as period grows
  — no period/geometry was chosen here, by design (`PHASE12B_MESH_PERIOD_STUDY.md`).
- Phase 10's own 5 "distinct" center geometries are not all mesh-distinct at the solver's actual (8-cell/side, dx=0.6 mm) mesh: **C2=(2.0,2.8) mm
  and C3=(2.8,2.8) mm are mesh-identical** (`CONCLUSIONS_AT_RISK`, Phase 12B's read-only audit). The held-out blind-test center (C5) has zero
  leakage into training, so Phase 10's core per-center findings stand, but C2 and C3 do not currently probe two independently-verified isolated
  baselines; the minimal fix (re-verify C2 on a finer/aligned mesh, or replace it) has not been applied, and Phase 10 itself has not been re-run.
- The Phase 7 library's own 324 geometries include 105 duplicate-mesh groups (169 distinct out of 324, `UNIT_CELL_LIBRARY.md` Sec. 2A) — a
  real, quantified consequence of its 32×32 mesh choice (made for build-time reasons, Sec. 2, not reconsidered here); the delivered library itself
  is unchanged and no result in it is wrong, but fewer than 324 of its rows carry independent geometric information.
- Phase 8E's own new candidate (D: a generalized 2x2 transfer operator, for a polarization-converting hologram) is `ARCHITECTURE_INCOMPATIBLE`:
  the registered library's measured cross-polarization is negligible (<1e-6), so this specific unrotated rectangular-patch family could not realize
  it even if the formulation were otherwise validated. Final Phase 8E status: `NO_VALIDATED_FORMULATION` (`PHASE8E_CONTROL_BRIDGE.md`) — an
  explicitly acceptable outcome (Sec. 22), not a pipeline failure; geometry inversion remains blocked project-wide.
- No coupled-multi-cell/near-field-array solver exists anywhere in this project, so the local-periodic approximation (isolated-unit-cell response
  standing in for its response under the actual spatially-varying hologram field) can never be reported as validated by anything currently
  implemented (`local_periodicity.py`'s `LOCAL_PERIODICITY_UNVALIDATED`, always).
- No real reference-field data (measured or HFSS-exported) exists anywhere in this project; Phase 8E's reference/object-wave analysis runs on an
  explicitly-labelled synthetic demonstration field (`SYNTHETIC_DEMO_SPHERICAL_FEED`), never presented as real hardware data.
- Oblique incidence ports (TE/TM, `cosθ` power normalization) and grating-order ports.
- Which handedness *naming* the synthesis and any reference use (the basis vectors are fixed; the right-hand-rule table is documented).
- Reuse of `Z` across right-hand sides/cells, parallelism, dense-`Z` cost for large slot complements; adaptive mesh selection.
- Mapping S-parameters to the phase/amplitude targets of the synthesis; database interpolation/lookup.
- Guided-mode poles of the ungrounded slab (low loss); substrate dispersion; finite conductivity/thickness of metal.
- Phase 12C's own first-pass step-limit figure (≤2.7°/cell, from only 4 gradient cases and the raw dipole-moment proxy)
  is **superseded by Phase 12D's** design-relevant reflection-phase/amplitude re-measurement, 20 step cases at a finer
  (24-cell/subcell) mesh: `LOCAL_PERIODICITY_STEP_LIMIT_MEASURED` up to a real **5.8°/cell**, off-resonance only. Every
  near/at-resonance step actually reachable at the largest embedded-solve mesh this project can afford (measured
  infeasible beyond ~32 cells/subcell — 409.5 s there, not finishing in 14+ minutes at 40) already fails badly (19–173°
  phase error, non-monotonically in step size) — the true near-resonance crossover step remains bounded, not pinned to a
  single number, but no near/at-resonance step this project can afford to embed-test ever stays within tolerance
  (`PHASE12D_LOCAL_PERIODICITY.md`).
- **The aperture-level test that actually decides Phase 13** (Phase 12D Task B, a genuinely new per-Floquet-order
  reflected-power capability, `phase12d/floquet_power.py`, validated to machine precision on a full-PEC-sheet trivial
  case and a uniform-array reproduction of \|S11\|²): real periodic phase-gradient supercells (N=6, 8, 12 cells per
  assumed 2π cycle, cells from a real isolated-cell lookup table, a full supercell MoM solve each time, not the dipole
  proxy) are `LOOKUP_TABLE_INADEQUATE` — no tested N puts more than 53% of reflected power in the intended order
  (29.5–47.5% measured), with 30–38% consistently lost to the backward (−1,0) diffraction order, because this design's
  own achievable phase span (186° of 360°, fixed by its own real S11 curve) is the same coverage shortfall Phases
  12A–12C already measured from a different angle. A lookup-table forward model for Phase 13 is therefore not
  recommended for this design point except for an aperture confined entirely to the flat, off-resonance part of the
  curve (`PHASE12D_LOCAL_PERIODICITY.md`).
- Real HFSS Floquet-port data now exists for exactly one single-cell case (G1, L=2.0 mm; HFSS 2025 R1, 18 modes/port,
  PEC ground) and agrees with the MoM reference to a mean phase difference of −0.29° and \|S11\| to ~1e-5 —
  `HFSS_SINGLE_CELL_VALIDATED_PARTIAL`, ingested via a new general Touchstone reader
  (`phase9.hfss_reference.read_touchstone`). The resonance case (G4, L=6.0 mm) is explicitly still pending, and no
  supercell HFSS data exists at all (`HFSS_NOT_AVAILABLE_FOR_SUPERCELLS`) — this remains the primary blocker on any
  broader HFSS validation claim (Phase 12D).
- Coupling for the grounded 8 mm cell is `GROUNDED_COUPLING_NOT_LOCAL_FOR_TESTED_SIZES`: none of 4 tested growing-pattern
  series (off-/near-resonance × gradient/checkerboard, n≤5, the largest routinely feasible size) converges to the solver's
  own self-consistency floor — coupling extent for the grounded case is not established to be any more local than the
  ungrounded Phase 11 finding, only tested over a smaller size range (Phase 12C).
- Phase 9's own original strong-contrast evidence for `LOCAL_PERIODICITY_INVALIDATED` (99.7% relative deviation) was
  measured at an 8-cell/subcell mesh; re-measured on a finer, aligned 48-cell/subcell mesh it drops to 21.8%, below Phase
  9's own 50% substantiality threshold — Phase 12C's own `EARLIER_CONCLUSIONS_CHANGED` finding, independently confirmed not
  to be a floquet-truncation artifact. This does not reverse non-locality overall (Phase 11's Case D and coupling-radius
  findings hold at the same finer mesh, and Phase 12C's own grounded coupling study is independently `NOT_LOCAL_FOR_TESTED_
  SIZES`), but it means the *original* ungrounded evidence for it was mesh-inflated, echoing Phase 12A/12B's own unaligned-
  mesh story more dramatically; Phase 9's dataset itself has not been re-run at this finer mesh.
- Phase 10's C2/C3 mesh-identity finding (Phase 12B, at the solver's original 8-cell/side mesh) also reverses on a finer,
  aligned 24-cell mesh: **C2 and C3 are mesh-distinct there** (relative isolated-response difference 0.547) — a further,
  separate instance of the same coarse-mesh-overstates-similarity pattern (Phase 12C); Phase 10 itself has not been re-run.
- The largest supercell size actually solved to convergence anywhere in this project is n=5 (aligned 20-cell/subcell mesh,
  104.6 s); n=7 at the same, validated `floquet_order = 2×total_mesh_cells` convention was attempted and killed after 42+
  minutes without completing — a measured infeasibility, not an assumed one (Phase 12C).

### Validation still required (independent of the code under test)
1. Babinet cross-check between a patch and its complementary slot.
2. Infinite wire grid / strip grating vs published closed-form grid impedances.
3. Magnitude **and phase** agreement with an independent full-wave solver (HFSS/CST/openEMS) for a patch, a dipole and a slot, after
   converting between time conventions if the reference uses `exp(−jωt)` (conjugate S-parameters) and checking the handedness naming.
   Phase 6 exports importable STEP models plus matching MoM reference results and a manual reproduction recipe for exactly this
   (`scripts/research/export_validation_cases.py`, [VALIDATION_EXPORT.md](VALIDATION_EXPORT.md)); the comparison itself is still manual and has
   not been performed by this project.

## Phases 13–17 (production repo onward)

**Phase 13** curates the repo to production-push grade (git, `pyproject.toml`, the `docs/` layout,
logging/diagnostics/error-taxonomy/validation-labels/run-manifest, ruff, CI) and builds a local
web app (FastAPI backend, static HTML/vanilla-JS frontend) that produces reflectarray/HMSA designs
now, labelled honestly, rather than waiting on the open questions above. A live UI review then a
redesign rebuilt the app around one 3-step HMSA wizard plus a "More tools" page for the rest, and
fixed a real physics bug (Formulation B's predicted aperture field was missing the reference-wave
phasor). A second live walk-through made the surface-wave eigenmode search's analytic scaffold
real where it could be validated, fixed a two-cause pattern side-lobe bug, made STEP export opt-in
with estimated progress, and added a real DXF export. Full account: `PHASE13_WEBAPP.md` (and its
own brief `briefs/PHASE13_WEBAPP_PROMPT.md`, `briefs/PHASE13B_UI_REVIEW_FIXES.md`,
`briefs/PHASE13C_REDESIGN.md`, `briefs/PHASE13D_POLISH_AND_PHYSICS.md`).

**Phase 14** built a real surface-wave dispersion source: a direct MoM eigenmode search over the
periodic array's own Z(kx,ky) was tried and found not to resolve a loading-dependent mode at any
affordable mesh size (kept as a documented research item, `src/dispersion/mom_scaffold.py`); the
working replacement extracts a patch-sheet admittance from the already-validated normal-incidence
S11 and solves the closed-form transverse-resonance condition for beta. The app gained a "Design
method" choice (surface-wave/eigenmode vs. reflection-phase) in Step 2. Full account:
`PHASE14_EIGENMODE.md` (brief `briefs/PHASE14_EIGENMODE_DISPERSION.md`).

**Phase 15** fixed a root-cause physics gap the user found in HFSS: every monopole-fed design
radiated a vortex (doughnut) beam because the feed's radial surface current, decomposed into
circular polarization, carries an orbital-angular-momentum term the app's scalar pattern model
couldn't see. Added a per-beam polarization-compensation term in hologram synthesis
(natural/LHCP/RHCP/HP/VP/slant, plus deliberate OAM vortex requests), a genuine VECTOR far-field
model via this project's own layered Green's function (replacing the scalar model for HMSA),
multi-feed ("multipole") reference-field construction, multi-port sub-apertures (Mishra et al.,
IEEE TVT 2026), and a coded (1-bit/N-bit) design method — each verified live in a real browser.
Full account: `PHASE15_POLARIZATION_VORTEX.md` (brief
`briefs/PHASE15_POLARIZATION_VORTEX_MULTIFEED.md`).

**Phase 16** (follow-ups, no write-up of its own — see `docs/CHANGELOG.md`'s `0.16.0` entry and
`briefs/PHASE16_FOLLOWUPS.md`) fixed the multi-beam saturation that was capping achievable levels,
added default-on level matching and optional pointing correction (later re-measured on a fine
local grid after a live review found it stalling on the shared, coarser plotting grid), fixed the
2×2+1 port layout's centring and default beams, reworded the isolation-gap check as information
rather than a fault, and unified the per-beam table's own measurement with the pointing/
level-matching report so the two always agree.

**Phase 17** was a handoff/cleanup pass — no physics or solver changes. It rewrites
the onboarding docs (`docs/START_HERE.md`, `docs/STATUS.md`, this reorganized `docs/history/`),
moves research drivers and one-off scripts out of the user-facing `examples/`, and indexes every
generated-output folder (`results/README.md`, `hfss_exports/README.md`) so a fresh clone is
legible in minutes. Brief: `docs/history/briefs/PHASE17_HANDOFF_AND_CLEANUP.md` (archived to
`docs/history/briefs/` once a later phase confirms this one is fully settled).

**Phase 18** (branch `phase18-four-port-300mm`, plus `phase18-one-bit-cell-check`) scaled the 4-port 1-bit coded board
to 300.8 mm (10.6 λ0) sub-apertures at guard gaps of 0 and 28 mm (`hfss_exports/four_port_1bit_300mm_gap*`) for a
minimum-gap/isolation HFSS study, added the held-out fine-grid per-port metrics to `scripts/adhoc/coded4.py` (the x/y-cut
co-pol SLL definition every later phase cites) and hand-built the square-cell-only variants (`scripts/adhoc/square4/`,
`hfss_exports/four_port_1bit_square_*`). A single-cell check of the 10.54 GHz 1-bit pair was started and parked.

**Phase 19** (`phase19-transmitarray`, tag `phase19-partD`) added a model-only 28 GHz all-dielectric (3D-printed nylon)
transmitarray: an RCWA unit-cell solver (`src/dielectric_cell/`, checked against analytic slab formulas to ~1e-14), a cell
library, a horn-fed dual-beam ±30° lens design with local-periodicity and Monte Carlo error-bar studies (`src/transmitarray/`,
`results/transmitarray_dual30/`) and an HFSS package (`hfss_exports/nylon_lens_28ghz_dual30/`). Parts E/F (app wiring, HFSS
comparison) wait on the user's ±30.00° decision and an HFSS run. Brief: `briefs/PHASE19_DIELECTRIC_TRANSMITARRAY_28GHZ.md`.

**Phase 20** (RIS: 2-bit PIN-diode reconfigurable metasurface) was only roadmapped and then paused; the briefs are in
`briefs/_parked/`. Nothing is implemented.

**Phase 21** (`phase21-complete`) made multi-port boards buildable in the app: N×M grid and custom-rectangle layouts, coded
1-bit/N-bit multi-port synthesis, an exact-coordinate board STEP, AEDT script and per-export README; it also fixed a
float floor bug that had cost the 300 mm boards one cell column (`docs/STATUS.md` known limitations).

**Phase 22** (`phase22-complete`) replaced hand-picked coded-cell states with automatic ones: an empty cell plus
mesh-corrected (Richardson) squares interpolated from mesh-aligned MoM solves, 1/2/3-bit, with a disk cache.

**Phase 23** (`phase23-complete`) added convex polygon (hexagon) port regions and the sidelobe-cancelling 1-bit optimiser:
a cell-wise far-field model, greedy q-norm flips (single-port, all-ports with port phases, robust over frequency and phase
error) measured on a held-out grid, with rectangle layouts reproducing bit-identically
(`results/phase23/`, `hfss_exports/four_port_hex150_*`).

**Phase 24** is this curation pass: branch consolidation into `main`, briefs and untracked outputs tidied, cache/clean targets,
an export-package growth policy, docs refreshed to describe phases 18-23. No computed result changed
(`docs/briefs/PHASE24_CURATION_AND_CLEANUP.md`).
