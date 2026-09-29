# Phase 8E -- physical control-variable and holographic-scattering bridge

## Decision tree (Sec. 29)

```
Reference field (SYNTHETIC_DEMO_SPHERICAL_FEED, no real data exists)
      |
Object wave (3-beam demo, phase8b.run.DEFAULT_DEMO_BEAMS, theta=30 deg)
      |
Candidate formulations A/B/C/D
      |
Required control variable C_req(x,y) (control_response.py)
      |
Architecture compatibility:
        A: ARCHITECTURE_INCOMPATIBLE
        B: ARCHITECTURE_COMPATIBLE_UNVALIDATED
        C: ARCHITECTURE_COMPATIBLE_UNVALIDATED
        D: ARCHITECTURE_INCOMPATIBLE
      |
Reference-plane compatibility: COMPATIBLE
      |
Complex controllability (S11xx for B, S21xx for C):
        B: mean complex error 0.7904
        C: mean complex error 0.7042
      |
Universal database query: family_exists=True
      |
Suitable family? YES
        continuous reconstruction -> discrete reconstruction -> geometry inversion permitted?
```

## Final per-formulation status

- **A**: `ARCHITECTURE_INCOMPATIBLE` (geometry_inversion_permitted=False)
  - architecture compatibility: ARCHITECTURE_INCOMPATIBLE
- **B**: `GEOMETRY_INVERSION_BLOCKED` (geometry_inversion_permitted=False)
  - architecture compatibility: ARCHITECTURE_COMPATIBLE_UNVALIDATED
  - reference-plane compatibility: COMPATIBLE
  - suitable unit-cell family registered: YES
  - continuous beam reconstruction: PASSED
  - discrete beam reconstruction: PASSED
  - physical control variable validated: NO (Phase 8B/8C: CONTROL_VARIABLE_UNVALIDATED, unchanged)
  - local periodicity acceptable: NO (local_periodicity.py: no coupled-cell solver exists to validate it)
- **C**: `GEOMETRY_INVERSION_BLOCKED` (geometry_inversion_permitted=False)
  - architecture compatibility: ARCHITECTURE_COMPATIBLE_UNVALIDATED
  - reference-plane compatibility: COMPATIBLE
  - suitable unit-cell family registered: YES
  - continuous beam reconstruction: PASSED
  - discrete beam reconstruction: PASSED
  - physical control variable validated: NO (Phase 8B/8C: CONTROL_VARIABLE_UNVALIDATED, unchanged)
  - local periodicity acceptable: NO (local_periodicity.py: no coupled-cell solver exists to validate it)
- **D**: `ARCHITECTURE_INCOMPATIBLE` (geometry_inversion_permitted=False)
  - architecture compatibility: ARCHITECTURE_INCOMPATIBLE

## Overall project status: `NO_VALIDATED_FORMULATION`

## Mandatory questions (Sec. 30)

### Q1: What exact physical quantity does the holographic formulation require?

Depends on the candidate: A needs a grounded surface reactance X(x,y); B needs a normal-incidence co-polarized reflection phase; C needs the identical phase kernel claimed against transmission; D needs a local 2x2 linear-basis transfer operator. No single answer is asserted (Sec. 5's explicit framework-not-a-choice instruction).

### Q2: What exact unit-cell electromagnetic quantity represents that physical quantity?

A: none available (no surface impedance quantity exists in this project's data, Phase 8B `investigate_surface_impedance` = INVALID). B: `S11_xx`. C: `S21_xx`. D: the full linear-basis 2x2 block (Sxx, Sxy, Syx, Syy).

### Q3: Are the reference field and unit-cell response defined at compatible reference planes?

For this run's synthetic reference field: **COMPATIBLE** (z_reference=0.0 m, pattern plane=0.0 m) -- compatible by construction (the spherical-feed model is evaluated directly at z=0). A real reference field recorded away from z=0 (e.g. Sec. 3's historical z=1.7 mm example) would instead report `BRIDGE_REFERENCE_PLANE_MISMATCH` and require the angular-spectrum propagation this run also demonstrated (`reference_plane.propagate_angular_spectrum`).

### Q4: Does the current patterned-PEC/no-ground-plane architecture satisfy the formulation assumptions?

A: **NO** (topology_audit.ground_plane = 'ABSENT': context['substrate'] = {'name': 'Rogers RT/Duroid 5880', 'eps_r': 2.2, 'tan_delta': 0.0009, 'thickness_m': 0.001575} (one dielectric layer, not a stack terminated in a conductor); context['metal'] = {'model': 'pec'} (a patterned PEC *sheet at z=0*, i.e. a finite scatterer, not a continuous ground plane); no 'ground' key exists anywhere in the context dict (checked: False); `config.py`/`materials.py`/`geometry.py`/`periodic_bc.py` contain zero occurrences of 'ground' (grep-verified) -- the config schema cannot express a ground plane at all. README.md 'Assumptions' item 2: 'No ground plane, no superstrate' (Phase 1, unchanged through Phase 7).). B: architecturally plausible but unvalidated (library-mean reflected power (co+cross) = 0.192194 (phase8c.power_balance.library_power_summary)). C: architecturally plausible but unvalidated (library-mean transmitted power (co+cross) = 0.806978 (phase8c.power_balance.library_power_summary)). D: **NO** (library-mean cross-polarized power: R_cross=8.44e-33, T_cross=1.03e-32 (threshold 1e-06) -- this unrotated rectangular-patch family does not realize a measurable polarization-converting response.).

### Q5: Is scalar S11 sufficient?

Not established. Complex-plane matching against the real 324-cell library (S11xx): mean complex error 0.7904, mean amplitude error 0.5227, mean phase error 68.93 deg -- reported as data, not collapsed into a pass/fail without a stated tolerance (none is defined by any project document).

### Q6: Is scalar S21 sufficient?

Not established. Complex-plane matching against the real 324-cell library (S21xx): mean complex error 0.7042, mean amplitude error 0.2602, mean phase error 51 deg.

### Q7: Is a matrix/vector scattering representation required?

Only if polarization conversion is actually needed by the desired beam set; for this run's demo beams (all x-polarized, x-polarized reference) the generalized operator reduces to a scalar (Txy=Tyx=Tyy=0). Independent of the specific demo, the *library's* measured cross-polarization is negligible (architecture check D: ARCHITECTURE_INCOMPATIBLE), so this unit-cell family could not realize a genuinely polarization-converting requirement even if one existed.

### Q8: What complex response region is actually required by the hologram?

See `required_response_map.csv` (the full per-aperture-point required phase for B/C) and `06_control_response_amplitude_phase.png`/`07_control_response_complex_plane.png`.

### Q9: Does the existing 4.8-mm library cover that region?

Partially, per the numbers in Q5/Q6 above (`control_variable_analysis.json`, `formulation_comparison.csv`) -- reported quantitatively, never reduced to a single coverage percentage.

### Q10: If not, what exact unit-cell library must be built?

`period_assessment.status` = `CURRENT_DATASET_LIMITATION`; see `required_library_specification.json` for the precise build request (topology/period/frequency/substrate/polarization/incidence), when one was generated.

### Q11: What period range is physically required?

Candidate (Nyquist-sampling heuristic, `CANDIDATE_HEURISTIC` only, never `VALIDATED`) periods: (13.846821, 11.077456, 8.308092, 5.538728, 2.769364) mm, from the demo beams' required transverse spatial frequency (226.882 rad/m). Registered (actually validated-by-data) periods: ((4.8, 4.8),) mm.

### Q12: Which conclusions are validated, and which remain hypotheses?

**Validated (measured, this run):** architecture is ungrounded and partially reflective/transmissive (Phase 8C, reused); cross-polarization is negligible across the whole registered library; the synthetic reference field used here is reference-plane-compatible by construction; the angular-spectrum propagation code path is demonstrated to run correctly for a nonzero reference height. **Hypotheses / unvalidated:** whether S11xx or S21xx phase is the physically correct holographic control variable (Phase 8B/8C, unchanged); whether the local-periodic approximation holds for this period/beam combination (no coupled-cell solver exists to check it); whether any candidate formulation's continuous reconstruction actually represents a physically buildable antenna (only a mathematical/Fourier-optics idealization was tested here).

## Validation status

- HFSS validation: NOT PERFORMED
- Solver / Phase 7 modified: NO
- Phase 8D-REV modified: NO
- Geometry inversion performed: NO
- Final MTS STEP generated: NO
