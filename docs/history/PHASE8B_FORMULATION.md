# Holographic formulation & control-variable validation (Phase 8B)

Establishes the mathematical bridge Phase 8A found missing: reference wave + object wave →
holographic modulation → physical unit-cell control variable → Phase 7 library quantity. Read-only
with respect to Phases 1-7 (no solver, mesh classifier, Green's function, or Phase 7 data touched);
builds no GUI and generates no final MTS geometry (per the brief's explicit scope limits).

```
CONVENTIONS.md / excitation.PlaneWaveExcitation (Phases 1-7, unmodified)
        |
        v
phase8b/conventions.py --> phase8b/object_wave.py --+--> phase8b/holographic_modulation.py
                                                      |         (2 candidate formulations)
phase8b/reference_wave.py ---------------------------+                |
                                                                       v
                                          phase8b/control_variable.py (S11xx/S21xx/Zs/Gamma investigated)
                                                                       |
                                                                       v
                                          phase8b/library_compatibility.py (decision tree, Sec. 28)
                                                                       |
                                                                       v
                                          results/phase8b_formulation/ (specification, reports, plots)
```

## The central finding

The brief's own worked example (`Psi = Re{E_ref . conj(E_obj)}`, a monopole reference wave) matches
one specific, real, citable antenna technique: **T. Fong, N. C. Colburn, J. J. Ottusch, J. L.
Visher, D. F. Sievenpiper, "Scalar and Tensor Holographic Artificial Impedance Surfaces," IEEE
Trans. Antennas Propag., vol. 58, no. 10, pp. 3212-3221, 2010** -- a monopole launches a guided
surface wave along a **grounded**, modulated-reactance surface; the physical control variable is
that surface's reactance `X(x,y)`.

This project's solver stack has **no ground plane** (an unchanged Phase 1 assumption) and **no
surface-wave/dispersion solver** of any kind. Two independent, structural reasons follow:

1. The standard `Zs = eta0(1+Gamma)/(1-Gamma)` relation this technique's unit-cell characterization
   normally uses assumes near-total reflection (`|Gamma|=1`, a grounded boundary); Phase 6's case01
   measured this project's actual patches as **mostly transmissive**
   (`|S11xx|` in [-15.89, -12.79] dB vs. `|S21xx|` in [-0.24, -0.12] dB, recomputed directly from
   `results/validation_cases/case01_non_resonant_patch/mom_results.json`), so that relation's own
   assumption does not hold for this data.
2. Even the *reference wave itself* (the guided surface wave's propagation constant `k_swg`) cannot
   be computed, because no surface-wave dispersion solver exists in `src/`.

Formulation A (Fong et al.) is therefore **`LIBRARY_CONTROL_VARIABLE_INCOMPATIBLE`** for this
project as it stands -- a structural finding, not a missing validation step.

A second, simpler, also well-established candidate -- **Formulation B**, a passive
reflectarray/transmitarray-style phase-difference hologram
(`required_phase = wrap(arg(E_obj) - arg(E_ref))`) -- needs only a normal-incidence co-polarized
reflection *or* transmission phase, exactly the *kind* of quantity Phase 7 actually stored
(S11xx/S21xx, both in full). Its status is **`PHASE_SOURCE_UNVALIDATED`**: the reflective-vs-
transmissive architecture decision (S11 vs. S21) that Phase 7 and Phase 8A already flagged as open
remains the one thing blocking it, not a structural incompatibility.

**Overall project status: `PHASE_SOURCE_UNVALIDATED`** (both formulations documented in full;
neither collapsed away -- see `holographic_specification.json`'s `holographic_modulation.formulations`).

## What is fully implemented, math-only, and numerically tested

`src/phase8b/object_wave.py`: single- and multi-beam plane-wave object fields, reusing
`excitation.PlaneWaveExcitation`'s transverse-wavenumber formula (cross-checked by test, since that
class cannot itself be constructed for every (theta, phi, polarization) triple an object-wave beam
needs -- see its module-level comment). Broadside/azimuth/phase-sign/multi-beam/polarization tests
all pass (`tests/test_phase8b_object_wave.py`, `test_phase8b_multi_beam.py`,
`test_phase8b_polarization.py`).

`src/phase8b/holographic_modulation.py`: both candidate kernels, numerically tested independent of
their (unresolved) physical validation status.

## What is investigated and honestly concluded

`src/phase8b/control_variable.py`: S11xx and S21xx are both `INSUFFICIENT_EVIDENCE` (well-defined
quantities; which one is physically correct is an open architecture decision); surface impedance is
`INVALID` (two structural reasons above); Gamma_xx=S11xx is `VALIDATED` as a definitional identity
(distinct from the surface-impedance *interpretation* of that same number).

## The 324-vs-169 finding, investigated further

`src/phase8b/distinct_states.py` confirms Phase 8A's 169-distinct-phase-states count and offers a
descriptive (not asserted, not solver-modifying) explanation: at the library's 32x32 mesh
(dx=dy=0.15 mm), the x-polarized S11xx response is dominated by L, so pairs of geometries with the
same L but different W can round to the identical discretized metal pattern -- consistent with, not
contradictory to, the mesh-coarseness limitation `LIBRARY_REPORT.md` already documents.

## Running it

```bash
.venv/bin/python scripts/research/run_phase8b_formulation.py     # results/phase8b_formulation/
```

## Verification

`tests/test_phase8b_conventions.py`, `test_object_wave.py`, `test_multi_beam.py`,
`test_polarization.py`, `test_holographic_modulation.py`, `test_reference_wave.py`,
`test_control_variable.py`, `test_library_compatibility.py`, `test_distinct_states.py`,
`test_specification.py`, `test_run.py` -- 71 tests, most against synthetic/small data, several
(skipped if absent) cross-checked against the real Phase 7 production library and the real
Phase 6 `case01` validation export.
