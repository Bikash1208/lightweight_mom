# repository_scientific_audit.md

Phase 8B Sec. 6: a systematic audit of this repository for every term Sec. 6 lists, distinguishing DEFINED / IMPLEMENTED / DOCUMENTED / UNDEFINED / ASSUMED. Performed by exhaustive grep across `src/`, `examples/`, `tests/`, and every top-level `.md` file before any Phase 8B code was written.

## Search terms

`holograph`, `object wave`, `object_wave`, `reference wave`, `reference_wave`, `reference field`, `reference_field`, `monopole`, `modulation`, `interference`, `surface impedance`, `surface_impedance`, `polarizability`, `S11`, `S21`, `polarization`, `phase convention`, `time convention`

## Findings

| Topic | Status | Evidence |
|---|---|---|
| Holographic generator / synthesis formulation | **UNDEFINED** | README.md:26 'Holographic synthesis is **not implemented**'; DATABASE.md:277; src/__init__.py:20. No `holograph*`/`object_wave`/`reference_wave` code anywhere in src/ before Phase 8. |
| Object wave | **UNDEFINED (before this phase)** | Zero occurrences of 'object wave' anywhere in the repository prior to Phase 8B. Now DEFINED (mathematically) by `src/phase8b/object_wave.py`, derived from `excitation.PlaneWaveExcitation`. |
| Reference wave / reference field | **PARTIALLY DEFINED (Phase 8A only)** | Zero occurrences of 'reference wave' before Phase 8A; Phase 8A's `reference_field.py` defines a generic reference-*field* CSV input format for validation purposes only, explicitly NOT an established project convention (its own module docstring says so). No monopole/guided-surface-wave reference-wave model existed before this phase. |
| Monopole feed | **UNDEFINED** | Zero occurrences of 'monopole' anywhere in the repository before this phase's own new files -- the Phase 8B brief's own text is the only place this project has mentioned it; no code, config, or prior documentation defines one. |
| Modulation / interference | **UNDEFINED** | Zero occurrences of 'modulation' or 'interference' (in the holographic sense) anywhere in the repository before this phase. |
| Surface impedance / polarizability | **UNDEFINED (as a computed quantity)** | materials.py:65 mentions surface impedance only to say finite conductivity is 'deliberately not' modelled (a PEC assumption note, not a surface-impedance-surface capability). No polarizability computation exists anywhere. |
| S11 as holographic phase source | **ASSUMED, THEN EXPLICITLY UN-ASSUMED** | `phase_extraction.resolve_phase_source()` (Phase 7) pins S11_xx as its lookup-table phase source with status='SPECIFICATION_MANDATED' (not 'PHYSICALLY_VALIDATED') and an explicit caveat that this is not an independently justified physical choice -- Phase 7 itself never treated this as an established fact. |
| S21 relevance | **DOCUMENTED** | `phase_extraction.py`/`UNIT_CELL_LIBRARY.md` Sec. 4 explicitly document that S21 is not negligible for this ungrounded structure, citing Phase 6 case01's measured magnitudes. |
| Time convention (exp(+j*omega*t)) | **DEFINED, IMPLEMENTED, DOCUMENTED** | CONVENTIONS.md 'Time convention'; `database.CONVENTIONS['time_dependence']`; used consistently throughout src/ (constants.py, excitation.py, green_function.py, ...). |
| Coordinate system / surface normal / x-y axes | **DEFINED, IMPLEMENTED, DOCUMENTED** | CONVENTIONS.md 'Coordinates and stack'; `geometry.py`. |
| Port definitions | **DEFINED, IMPLEMENTED, DOCUMENTED** | CONVENTIONS.md 'Waves, ports and S-parameters'; S_PARAMETERS.md; `excitation.ReferencePlanes`. |
| Incident-field / plane-wave convention | **DEFINED, IMPLEMENTED, DOCUMENTED** | CONVENTIONS.md 'Incident wave'; `excitation.PlaneWaveExcitation` (kx0/ky0/kz0 formula, reused -- cross-checked by test -- by `src/phase8b/object_wave.py`). |
| Polarization convention (linear x/y, derived circular) | **DEFINED, IMPLEMENTED, DOCUMENTED** | CONVENTIONS.md 'Circular polarization...'; `polarization.py`; DATABASE.md Sec. 7. |
| Phase sign / wrapping convention | **DEFINED, IMPLEMENTED, DOCUMENTED** | CONVENTIONS.md; `phase.py` (principal value on (-180, 180], delay = negative phase). |
| Oblique-incidence / surface-wave ports | **UNDEFINED / NOT IMPLEMENTED** | CONVENTIONS.md 'Not yet fixed'; README.md 'Still open': 'Oblique incidence ports... and grating-order ports' -- not yet implemented. No surface-wave dispersion solver exists. |
| Reference-grid to MTS-grid spatial mapping | **UNDEFINED** | Phase 8A `feasibility_spatial.check_mapping_method_defined()`: no such module exists in src/. |

## Conclusion

No holographic formulation, object-wave definition, reference-wave definition, monopole model, or modulation/interference equation existed anywhere in this repository before Phase 8B. Every existing solver convention needed to build one *consistently* (time convention, coordinate system, port definitions, incident-field convention, polarization convention, phase-sign convention) was, however, already fully defined, implemented and documented (Phases 1-7) -- Phase 8B reuses every one of them rather than introducing a second, parallel set (`src/phase8b/conventions.py`).
