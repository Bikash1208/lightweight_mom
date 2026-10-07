# Physical control-variable and holographic-scattering bridge (Phase 8E)

Determines what physical complex unit-cell response the holographic formulation actually requires,
and whether that requirement can be physically connected to the real, measured unit-cell scattering
data -- *before* any geometry inversion. Builds strictly on Phase 8B (formulation math), Phase 8C
(scattering-architecture findings) and Phase 8D-REV (the universal database registry); none of
those are modified or re-run.

```
reference field (synthetic demo -- no real HFSS/measured data exists in this project)
        |
object wave (3-beam demo, phase8b.run.DEFAULT_DEMO_BEAMS)
        |
candidate formulations A/B/C/D (provenance.py)
        |
required control variable C_req(x,y) (control_response.py)
        |
architecture compatibility (architecture_compatibility.py, reuses phase8c.topology_audit)
        |
reference-plane compatibility (reference_plane.py: real angular-spectrum propagation, never an
        |                       arbitrary phase offset)
complex-plane controllability (complex_controllability.py, reuses phase8c/phase8d_rev)
        |
universal-database query (reuses phase8d_rev.query, unmodified)
        |
continuous vs. discrete beam reconstruction (beam_reconstruction.py: standard aperture Fourier-
        |                                     optics far-field integral)
geometry-inversion gate (decision_tree.py, Sec. 20's 9-checkbox gate)
```

## The four candidates

| | required quantity | architecture status (measured) |
|---|---|---|
| A (Fong et al. 2010, grounded impedance) | surface reactance of a grounded structure | `ARCHITECTURE_INCOMPATIBLE` -- no ground plane |
| B (phase-difference hologram, reflection) | `arg(S11xx)` | `ARCHITECTURE_COMPATIBLE_UNVALIDATED` |
| C (phase-difference hologram, transmission) | `arg(S21xx)` | `ARCHITECTURE_COMPATIBLE_UNVALIDATED` |
| D (generalized 2x2 transfer operator) | full linear-basis Sxx/Sxy/Syx/Syy | `ARCHITECTURE_INCOMPATIBLE` -- measured cross-pol <1e-6 |

Candidates B and C use the *identical* mathematical kernel (`phase8b.holographic_modulation.
phase_difference_hologram`); only the physical correspondence claim (reflected vs. transmitted
field) differs, and neither is validated (Phase 8B/8C's `CONTROL_VARIABLE_UNVALIDATED`, unchanged).

## The real result that came out of this bridge

Running the pipeline's demo (3 beams at theta=30 deg, phi=0/120/240, a synthetic spherical-feed
reference at the pattern plane z=0):

- **Continuous reconstruction is exact**: the idealized phase-only hologram places all three beams
  at exactly their intended directions (0.0 deg pointing error) -- the phase-difference-hologram
  *mathematics* genuinely works, independent of which unit-cell family realizes it.
- **Discrete reconstruction (using the real, non-unit-magnitude S11xx/S21xx library values)
  degrades it**: pointing error grows to 4.6-8.1 deg and the matched beam peaks drop several dB,
  because the real 324-cell library's achievable magnitude range ([0.23, 0.83] for S11xx, [0.56,
  0.97] for S21xx) never matches the idealized unit-magnitude phase-only target. Mean phase error
  against the nearest achievable complex response is ~67 deg for both S11xx and S21xx -- reported
  as data, not collapsed into a single pass/fail number (no tolerance is defined anywhere in this
  project).

## Final status: `NO_VALIDATED_FORMULATION`

No candidate's geometry-inversion gate (Sec. 20's 9 checkboxes) is fully open: A and D fail at the
architecture-compatibility step; B and C pass every reconstruction check but the *physical control
variable itself* (S11xx-vs-S21xx-as-holographic-phase) remains unvalidated (Phase 8B/8C, unchanged)
and the local-periodic approximation is always reported `LOCAL_PERIODICITY_UNVALIDATED` (no
coupled-multi-cell solver exists in this project to check it). Per Sec. 22, this is an explicitly
acceptable scientific outcome, not a failure of the pipeline -- geometry inversion remains blocked.

## Running it

```bash
.venv/bin/python scripts/research/run_phase8e_control_bridge.py     # results/phase8e_control_bridge/
```

## Verification

`tests/test_phase8e_*.py` (73 tests): pure-math unit tests for the spectral/reference-plane/
beam-reconstruction/decision-tree logic (including an exact round-trip check of the angular-
spectrum propagation and a synthetic single-point-aperture far-field sanity check), plus direct
checks against the real 324-cell library for architecture compatibility and complex-plane
controllability. `test_run_beam_reconstruction_identifies_all_three_demo_beams` is the key
end-to-end physics check: the continuous reconstruction must place a detectable peak within the
far-field grid's own angular resolution of each of the three demo beams.
