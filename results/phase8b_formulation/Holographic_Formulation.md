# Holographic_Formulation.md

Two literature-sourced candidate holographic formulations, evaluated against this project's actual, frozen Phase 1-7 conventions and Phase 7 library -- neither is asserted as *the* project's formulation (no project document authorizes either one; `repository_scientific_audit.md`).

## Coordinate system

x: x in [-Px/2, Px/2]; y: y in [-Py/2, Py/2]; z: air (z>0, incident side) | pattern at z=0 | substrate (-h<z<0) | air (z<-h); no ground plane.

## Time convention

`exp(+j*omega*t)` (unchanged, CONVENTIONS.md).

## Reference field

Formulation A (Sec. below): a monopole-launched guided surface wave, `E_ref(rho) ~ rho^(-1/2) exp(-j k_swg rho)` (`reference_wave.monopole_reference_field`); `k_swg` has no default and cannot be computed by this project's solver. Formulation B: an arbitrary reference field (any of Sec. 9's three input types; `reference_field.py`, Phase 8A).

## Object field

`E_obj(x,y) = sum_n w_n exp(j*phi0_n) exp(-j(kx0_n x + ky0_n y))` at z=0, `kx0_n = k0 sin(theta_n) cos(phi_n)`, `ky0_n = k0 sin(theta_n) sin(phi_n)` -- the identical formula `excitation.PlaneWaveExcitation.transverse_wavenumber` implements, cross-checked by test (`src/phase8b/object_wave.py`; that class cannot be constructed directly for every (theta, phi, polarization) triple an object-wave beam needs, see its module-level comment).

## Multi-beam formulation

Field-level complex superposition (`object_wave.combine_beams`); weight = complex field-amplitude coefficient (Sec. 16); `initial_phase_deg` enters as one constant additive phase per beam (Sec. 17), no other role.

## Holographic modulation

### Formulation A: Fong et al. 2010 interference-modulation index (grounded holographic impedance surface)

Required quantity: surface reactance X(x,y) of a grounded, surface-wave-supporting structure.
Library compatibility: **INCOMPATIBLE**.
Decision trail: [Does a holographic formulation exist?: YES (Fong et al. 2010, source_registry.json fong2010)] -> [Does it define a physical control variable?: YES (surface reactance X(x,y) of a grounded structure)] -> [Can the Phase 7 library provide that quantity?: NO] -> [Result: LIBRARY GAP: no ground plane (invalidates the standard Zs<->Gamma relation) and no surface-wave dispersion solver (cannot characterize the reference wave either)]
Final status: **LIBRARY_CONTROL_VARIABLE_INCOMPATIBLE**.

### Formulation B: phase-difference hologram (reflectarray/transmitarray-style)

Required quantity: normal-incidence co-polarized reflection or transmission phase (S11xx or S21xx).
Library compatibility: **PARTIALLY_COMPATIBLE**.
Decision trail: [Does a holographic formulation exist?: YES (phase-difference/phase-conjugate hologram, standard reflectarray/transmitarray design equation)] -> [Does it define a physical control variable?: YES (normal-incidence co-pol reflection OR transmission phase -- S11xx or S21xx)] -> [Can the Phase 7 library provide that quantity?: YES (both S11xx and S21xx are stored in full, req. 5) -- but which one is correct is not yet decided] -> [Is the relationship physically validated?: NO (reflective-vs-transmissive architecture undecided; S11xx/S21xx both INSUFFICIENT_EVIDENCE)]
Final status: **PHASE_SOURCE_UNVALIDATED**.

## Physical control variable

Unresolved: two candidates (S11xx/S21xx for Formulation B; surface reactance for Formulation A, INVALID for this project). See `control_variable_analysis.md`.

## Unit-cell mapping

Not established (depends on the control-variable resolution above); no geometry interpolation/extrapolation is implemented (unchanged from Phase 8A).

## 324 geometries vs. distinct phase states

324 geometries, 169 distinct phase states (155 duplicate records). 155 of 324 geometries share a bit-identical (within 1e-06 deg) stored phase with at least one other geometry (169 distinct phase states). This is consistent with, not contradictory to, the mesh-coarseness limitation `LIBRARY_REPORT.md` Sec. 16 already documents (32x32 cells over the 4.8 mm period gives dx = dy = 0.15 mm, `library_schema.DEFAULT_MESH`): a sufficiently large L, when discretized to whole 0.15 mm mesh cells, can round the *same* number of metal cells for two nominally different W values whenever the physically relevant (x-polarized) response is dominated by L rather than W -- an ordinary staircase-quantization effect on the L axis's neighbourhood, not an independent bug. This module does not modify the mesh or solver to test this hypothesis further (out of Phase 8B's scope, Sec. 3); it is reported descriptively, as Sec. 30 requires, not asserted as certain.

