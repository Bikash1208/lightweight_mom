# control_variable_analysis.md

## 1. Candidate physical quantities

S11xx (reflection), S21xx (transmission), surface impedance Zs, reflection coefficient Gamma (Sec. 22-25).

## 2. S11xx analysis

(1) Physically: the complex ratio of the reflected to incident tangential-E wave amplitude at port 1 (top reference plane), co-polarized x-to-x (`S_PARAMETERS.md`). (2) Incident field: a normally-incident (theta=0), x-polarized plane wave, unit amplitude referenced at the pattern plane z=0 (`excitation.py`). (3) Reflected field: the x-component of the scattered wave travelling +z, referenced at z_top (CONVENTIONS.md 'Ports'). (4) Reference plane: z_top = +top_air_thickness, in air, above the pattern (not at the metal itself). (5) Polarization: linear x, co-polarized. (6) The unit cell is NOT grounded (README.md 'Assumptions' item 2, unchanged since Phase 1: air | pattern | substrate | air, no ground plane). (7) Whether the intended MTS effect is reflective or transmissive is an architecture decision this project has never made (`UNIT_CELL_LIBRARY.md` Sec. 4, `phase_extraction.py`); Phase 6's case01 measured the structure as *mostly transmissive*: |S11xx| in (-15.89, -12.79) dB vs. |S21xx| in (-0.24, -0.12) dB (recomputed directly from `results/validation_cases/case01_non_resonant_patch/mom_results.json`), the opposite of a typical grounded reflectarray element where |S11|~=0 dB (near-total reflection). (8) Whether the holographic modulation is represented by a reflection coefficient at all depends on which candidate formulation is used: Formulation B (`holographic_modulation.phase_difference_hologram`) can use *either* S11xx or S21xx phase as its control variable in principle (both are normal-incidence co-pol phases of the right *kind*); Formulation A (Fong et al. 2010) needs a **surface reactance** of a **grounded** structure, which S11xx of this ungrounded, partially-transmissive cell is not (see `investigate_surface_impedance`). (9) S11xx phase alone is not obviously sufficient without a magnitude/usability floor -- Phase 7 already applies one (`library_schema.PHASE_USABLE_MIN_MAG_DB`). (10) Only Formulation B's kernel could use reflection-coefficient phase directly; Formulation A cannot (structural reason above). (11) Phase 7 stores S11xx as `arg(S)` on (-180, 180], a fixed, consistent reference (`phase.py`, `database.CONVENTIONS`). Conclusion: S11xx is a well-defined, correctly-stored physical quantity, but whether it is *the* intended holographic control variable depends on an architecture decision (reflective-vs-transmissive; which candidate formulation) this project has not made -- neither VALIDATED nor INVALID can be honestly asserted; see `library_compatibility.py` for the per-formulation resolution this ambiguity allows.

**Status: INSUFFICIENT_EVIDENCE**

## 3. S21xx analysis

S21xx is the complex ratio of the transmitted (port 2, travelling -z) to incident (port 1, travelling -z) tangential-E wave amplitude, co-polarized x-to-x. This is NOT a grounded structure (no ground plane), so unlike a typical reflectarray/AMC unit cell, S21xx is *not* negligible here: |S21xx| in (-0.24, -0.12) dB (near-total transmission) vs. |S11xx| in (-15.89, -12.79) dB for Phase 6's case01 (same source as `investigate_s11xx`). This makes S21xx phase, not S11xx phase, the *larger*-magnitude, better-conditioned response for this specific ungrounded geometry family -- but a larger magnitude alone is not sufficient justification to prefer it (Sec. 23's explicit prohibition: 'do not select S21 merely because it offers a larger phase span/better conditioning'). Whether S21xx is the physically meaningful quantity depends on whether the intended MTS architecture is transmissive (a transmit-array-style hologram, e.g. the object wave radiates from the side opposite the feed) or reflective -- an architecture decision this project has never made.

**Status: INSUFFICIENT_EVIDENCE**

## 4. Reflection coefficient analysis

S11xx *is* mathematically the reflection coefficient Gamma_xx at port 1 by construction (`S_PARAMETERS.md` Sec. 1: b1 = S11 a1 + S12 a2; with only port 1 excited, S11 = b1/a1 = Gamma). This equivalence (Gamma_xx = S11_xx) is physically justified and requires no additional assumption -- unlike the *surface-impedance* interpretation of that same Gamma (`investigate_surface_impedance`, which is INVALID for this ungrounded structure), 'S11xx is a reflection coefficient' is simply what S11xx already, definitionally, is (`database.py` module docstring, `S_PARAMETERS.md`).

**Status: VALIDATED**

## 5. Surface impedance analysis

Formulation A (Fong et al. 2010) needs the surface reactance X(x,y) of a **grounded** modulated-impedance boundary, most commonly obtained in the literature from the closed-form normal-incidence relation Zs = eta0 (1+Gamma)/(1-Gamma) (standard transmission-line/impedance-boundary reflection-coefficient relation, `source_registry.json` `impedance_boundary_relation`), where Gamma is the reflection coefficient of a **grounded** unit cell (|Gamma|=1 for a lossless grounded boundary, so Zs is purely reactive). Phase 1-7's solver stack has **no ground plane** (`phase8b.conventions.NO_GROUND_PLANE`) and Phase 6's case01 shows the structure is mostly transmissive (|S11xx| in (-15.89, -12.79) dB, far from |Gamma|=1=0 dB). Applying Zs=eta0(1+S11xx)/(1-S11xx) to this data anyway would not fail numerically, but it would not represent a physical surface impedance: the formula's derivation assumes total reflection from a boundary condition of the form E_t = Zs (n x H_t), which this partially-transmissive, ungrounded two-port structure does not satisfy (some power genuinely leaves through the back, which the single-boundary impedance model has no way to represent -- it would appear as a large, physically meaningless 'loss' term in a structure whose actual dielectric loss is tan(delta)=0.0009, negligible). Separately, the formulation's reference wave (`reference_wave.monopole_reference_field`) is itself a *guided surface wave*, whose propagation constant k_swg this project's solver cannot compute (`phase8b.conventions.SURFACE_WAVE_DISPERSION_IMPLEMENTED = False`: no eigenmode/dispersion solver exists in `src/`). Conclusion: SURFACE_IMPEDANCE_NOT_AVAILABLE -- two independent, structural reasons (no ground plane invalidates the standard Zs<->Gamma relation; no surface-wave dispersion solver exists to characterize the reference wave itself), neither of which is a missing feature that could be worked around without changing the frozen solver.

**Status: INVALID**

## 6. Polarization compatibility

Both S11xx/S21xx are co-polarized (x-in, x-out) quantities; an x-polarized library response cannot, by itself, encode an arbitrary (e.g. dual-linear or circular) desired object-wave polarization without an additional, unvalidated assumption -- not made here.

## 7. Reference-plane implications

S11xx/S21xx are referenced at the top/bottom air reference planes (z_top = +top_air_thickness (above the pattern); z_bottom = -(h + bottom_air_thickness) (below the substrate)), not at the metal pattern itself; a holographic formulation using them must account for this fixed reference, not assume phase is referenced at z=0.

## 8. Time/phase convention

`exp(+j*omega*t)`, phase on (-180.0, 180.0) deg -- unchanged, reused throughout `src/phase8b/`.

## 9. Phase 7 library compatibility

- A: Fong et al. 2010 interference-modulation index (grounded holographic impedance surface): **INCOMPATIBLE** (final status LIBRARY_CONTROL_VARIABLE_INCOMPATIBLE).
- B: phase-difference hologram (reflectarray/transmitarray-style): **PARTIALLY_COMPATIBLE** (final status PHASE_SOURCE_UNVALIDATED).

## 10. Final scientific status

See `holographic_specification.json`'s top-level `status`.
