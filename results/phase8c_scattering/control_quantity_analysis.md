# control_quantity_analysis.md

## 1. Physical unit-cell architecture

PARTIALLY_REFLECTIVE_TRANSMISSIVE -- more transmissive than reflective (mean reflected power 0.1922, transmitted 0.8070). Ground plane: **ABSENT** (context['substrate'] = {'name': 'Rogers RT/Duroid 5880', 'eps_r': 2.2, 'tan_delta': 0.0009, 'thickness_m': 0.001575} (one dielectric layer, not a stack terminated in a conductor); context['metal'] = {'model': 'pec'} (a patterned PEC *sheet at z=0*, i.e. a finite scatterer, not a continuous ground plane); no 'ground' key exists anywhere in the context dict (checked: False); `config.py`/`materials.py`/`geometry.py`/`periodic_bc.py` contain zero occurrences of 'ground' (grep-verified) -- the config schema cannot express a ground plane at all. README.md 'Assumptions' item 2: 'No ground plane, no superstrate' (Phase 1, unchanged through Phase 7).)

## 2. Geometry/topology

one patterned PEC sheet at z=0 (context['metal']={'model': 'pec'})

one dielectric layer (context['substrate']={'name': 'Rogers RT/Duroid 5880', 'eps_r': 2.2, 'tan_delta': 0.0009, 'thickness_m': 0.001575})

two semi-infinite air half-spaces (z>0 and z<-h), each with an explicit reference plane (context['reference_planes']={'top_air_thickness_m': 0.005, 'bottom_air_thickness_m': 0.005})

## 3. Port definition

See `s_parameter_definition.md` (CONVENTIONS.md/S_PARAMETERS.md, reused verbatim).

## 4. S-parameter definitions

See `s_parameter_definition.md`.

## 5. Candidate control quantities

S11xx, S21xx (investigated in full below); S11yy/S21yy/cross-polarized terms investigated via the power-balance analysis (Sec. 6-8 below) using the full linear-basis matrix read from the internal Phase 7 checkpoint (`phase8c.data_source`).

## 6. Reflection analysis

Holographic control variable = reflection coefficient (S11xx)

[reference formulation defines a reflection-phase control variable (Formulation B): YES] -> [required reflection quantity = normal-incidence co-pol reflection phase: DEFINED] -> [S11xx interpretation as that quantity: INSUFFICIENT_EVIDENCE] -> [Phase 7 response available (S11 stored in full): YES]

**REFLECTION_ARCHITECTURE_UNVALIDATED**

## 7. Transmission analysis

Holographic control variable = transmission coefficient (S21xx)

[reference formulation defines a transmission-phase control variable (Formulation B): YES] -> [required transmission quantity = normal-incidence co-pol transmission phase: DEFINED] -> [S21xx interpretation as that quantity: INSUFFICIENT_EVIDENCE] -> [Phase 7 response available (S21 stored in full): YES]

**TRANSMISSION_ARCHITECTURE_UNVALIDATED**

## 8. Power balance

Library-mean (x-incidence): R_co=0.1922, R_cross=0.0000, T_co=0.8070, T_cross=0.0000, mean(R+T)=0.9992 (range [0.9990, 0.9993]).

## 9. Modal analysis

Only the fundamental (0,0) Floquet order propagates throughout the analyzed band (onset at 62.46 GHz, far above it); S11/S21 fully describe the specular scattered power.

## 10. Complex-response analysis

S11xx: 324 geometries, 169 distinct complex states, 169 distinct phase states.
S21xx: 324 geometries, 169 distinct complex states, 169 distinct phase states.

## 11. Phase coverage

S11xx: 24.84 deg of 360 deg (largest gap 335.16 deg). S21xx: 43.85 deg of 360 deg (largest gap 316.15 deg).

## 12. Magnitude/phase coupling

S11xx: magnitude in [0.2311, 0.8305] (median 0.3687, std 0.1534); phase-magnitude correlation -0.9157224804249091.
S21xx: magnitude in [0.5564, 0.9726] (median 0.9291, std 0.0897); phase-magnitude correlation 0.9849151865534151.
Phase-only approximation: S11xx PHASE_ONLY_APPROXIMATION_UNRESOLVED; S21xx PHASE_ONLY_APPROXIMATION_UNRESOLVED.

## 13. Polarization

Both S11xx/S21xx are co-polarized (x-in, x-out); an x-polarized response cannot, by itself, encode an arbitrary object-wave polarization without an additional assumption (Phase 8B `control_variable_analysis.md` Sec. 6).

## 14. Reference-plane considerations

See `reference_plane_sensitivity.py`: phase is reference-plane dependent by a documented, reused (not re-derived) formula; not corrected or offset anywhere in this analysis.

## 15. Phase 8B compatibility

Formulation A: LIBRARY_CONTROL_VARIABLE_INCOMPATIBLE (unchanged from Phase 8B).
Formulation B: corresponds to reflection phase -> REFLECTION_ARCHITECTURE_UNVALIDATED; to transmission phase -> TRANSMISSION_ARCHITECTURE_UNVALIDATED.

## 16. Library adequacy

See `compatibility_matrix.py`'s factual table below.

| Criterion | S11xx | S21xx |
|---|---|---|
| Physical interpretation | complex co-polarized (x-in, x-out) reflection coefficient at the top reference plane, normal incidence, x-polarized | complex co-polarized (x-in, x-out) transmission coefficient, port 1 to port 2, normal incidence, x-polarized |
| Architecture compatible | REFLECTION_ARCHITECTURE_UNVALIDATED | TRANSMISSION_ARCHITECTURE_UNVALIDATED |
| Phase coverage (deg) | 24.84 | 43.85 |
| Largest circular gap (deg) | 335.16 | 316.15 |
| Magnitude range | [0.2311, 0.8305] | [0.5564, 0.9726] |
| Phase-only approximation | PHASE_ONLY_APPROXIMATION_UNRESOLVED | PHASE_ONLY_APPROXIMATION_UNRESOLVED |
| Polarization | co-polarized (x-in, x-out) | co-polarized (x-in, x-out) |
| Holographic mapping (Formulation B) | PARTIAL (candidate, unresolved) | PARTIAL (candidate, unresolved) |
| Validation status | INSUFFICIENT_EVIDENCE | INSUFFICIENT_EVIDENCE |

## 17. Missing evidence

- HFSS validation data (not performed for any cell)
- A quantitative magnitude tolerance for a phase-only approximation, from either Phase 8B formulation
- A project decision on reflective-vs-transmissive intended architecture

## 18. Final scientific status

**CONTROL_VARIABLE_UNVALIDATED**

Reciprocity check: max|S21xx-S12xx|=1.144e-15, max|S21yx-S12xy|=6.678e-16 over 324 cells -- consistent with `MOM_FORMULATION.md`'s claim that Z is complex-symmetric at normal incidence (S=S^T), within solver numerical precision.

Symmetry check: S11_xx(L,W) vs S11_yy(W,L): max abs error 4.125e-15 over 171 pairs (within tolerance 1.0e-03).
