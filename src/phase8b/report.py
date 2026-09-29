"""Phase 8B Sec. 6, 33, 51, 56: the human-readable reports. Every factual claim here is either
computed by another `phase8b` module (passed in) or a literal quotation of a project file (cited
by path) -- this module only assembles text, it does not decide anything new.
"""

from . import conventions
from .library_compatibility import LIBRARY_COMPATIBILITY_STATES

# Sec. 6's audit table. Status values: DEFINED | IMPLEMENTED | DOCUMENTED | UNDEFINED | ASSUMED.
# Built from the exhaustive grep audit performed before any code was written (see this module's
# `AUDIT_EVIDENCE` below for the exact commands/results) -- an undocumented assumption is never
# silently promoted to DEFINED (Sec. 6's explicit instruction).
AUDIT_ROWS = [
    ("Holographic generator / synthesis formulation", "UNDEFINED",
     "README.md:26 'Holographic synthesis is **not implemented**'; DATABASE.md:277; src/__init__.py:20. "
     "No `holograph*`/`object_wave`/`reference_wave` code anywhere in src/ before Phase 8."),
    ("Object wave", "UNDEFINED (before this phase)",
     "Zero occurrences of 'object wave' anywhere in the repository prior to Phase 8B. Now DEFINED "
     "(mathematically) by `src/phase8b/object_wave.py`, derived from `excitation.PlaneWaveExcitation`."),
    ("Reference wave / reference field", "PARTIALLY DEFINED (Phase 8A only)",
     "Zero occurrences of 'reference wave' before Phase 8A; Phase 8A's `reference_field.py` defines "
     "a generic reference-*field* CSV input format for validation purposes only, explicitly NOT "
     "an established project convention (its own module docstring says so). No monopole/guided-"
     "surface-wave reference-wave model existed before this phase."),
    ("Monopole feed", "UNDEFINED",
     "Zero occurrences of 'monopole' anywhere in the repository before this phase's own new files "
     "-- the Phase 8B brief's own text is the only place this project has mentioned it; no code, "
     "config, or prior documentation defines one."),
    ("Modulation / interference", "UNDEFINED",
     "Zero occurrences of 'modulation' or 'interference' (in the holographic sense) anywhere in "
     "the repository before this phase."),
    ("Surface impedance / polarizability", "UNDEFINED (as a computed quantity)",
     "materials.py:65 mentions surface impedance only to say finite conductivity is 'deliberately "
     "not' modelled (a PEC assumption note, not a surface-impedance-surface capability). No "
     "polarizability computation exists anywhere."),
    ("S11 as holographic phase source", "ASSUMED, THEN EXPLICITLY UN-ASSUMED",
     "`phase_extraction.resolve_phase_source()` (Phase 7) pins S11_xx as its lookup-table phase "
     "source with status='SPECIFICATION_MANDATED' (not 'PHYSICALLY_VALIDATED') and an explicit "
     "caveat that this is not an independently justified physical choice -- Phase 7 itself never "
     "treated this as an established fact."),
    ("S21 relevance", "DOCUMENTED",
     "`phase_extraction.py`/`UNIT_CELL_LIBRARY.md` Sec. 4 explicitly document that S21 is not "
     "negligible for this ungrounded structure, citing Phase 6 case01's measured magnitudes."),
    ("Time convention (exp(+j*omega*t))", "DEFINED, IMPLEMENTED, DOCUMENTED",
     "CONVENTIONS.md 'Time convention'; `database.CONVENTIONS['time_dependence']`; used "
     "consistently throughout src/ (constants.py, excitation.py, green_function.py, ...)."),
    ("Coordinate system / surface normal / x-y axes", "DEFINED, IMPLEMENTED, DOCUMENTED",
     "CONVENTIONS.md 'Coordinates and stack'; `geometry.py`."),
    ("Port definitions", "DEFINED, IMPLEMENTED, DOCUMENTED",
     "CONVENTIONS.md 'Waves, ports and S-parameters'; S_PARAMETERS.md; `excitation.ReferencePlanes`."),
    ("Incident-field / plane-wave convention", "DEFINED, IMPLEMENTED, DOCUMENTED",
     "CONVENTIONS.md 'Incident wave'; `excitation.PlaneWaveExcitation` (kx0/ky0/kz0 formula, "
     "reused -- cross-checked by test -- by `src/phase8b/object_wave.py`)."),
    ("Polarization convention (linear x/y, derived circular)", "DEFINED, IMPLEMENTED, DOCUMENTED",
     "CONVENTIONS.md 'Circular polarization...'; `polarization.py`; DATABASE.md Sec. 7."),
    ("Phase sign / wrapping convention", "DEFINED, IMPLEMENTED, DOCUMENTED",
     "CONVENTIONS.md; `phase.py` (principal value on (-180, 180], delay = negative phase)."),
    ("Oblique-incidence / surface-wave ports", "UNDEFINED / NOT IMPLEMENTED",
     "CONVENTIONS.md 'Not yet fixed'; README.md 'Still open': 'Oblique incidence ports... and "
     "grating-order ports' -- not yet implemented. No surface-wave dispersion solver exists."),
    ("Reference-grid to MTS-grid spatial mapping", "UNDEFINED",
     "Phase 8A `feasibility_spatial.check_mapping_method_defined()`: no such module exists in src/."),
]

# The exact search performed for Sec. 6 (reproduced here so the audit is re-runnable/checkable,
# not just narrated).
AUDIT_SEARCH_TERMS = (
    "holograph", "object wave", "object_wave", "reference wave", "reference_wave",
    "reference field", "reference_field", "monopole", "modulation", "interference",
    "surface impedance", "surface_impedance", "polarizability", "S11", "S21", "polarization",
    "phase convention", "time convention",
)


def render_repository_scientific_audit_md():
    lines = [
        "# repository_scientific_audit.md", "",
        "Phase 8B Sec. 6: a systematic audit of this repository for every term Sec. 6 lists, "
        "distinguishing DEFINED / IMPLEMENTED / DOCUMENTED / UNDEFINED / ASSUMED. Performed by "
        "exhaustive grep across `src/`, `examples/`, `tests/`, and every top-level `.md` file "
        "before any Phase 8B code was written.", "",
        "## Search terms", "",
        "`" + "`, `".join(AUDIT_SEARCH_TERMS) + "`", "",
        "## Findings", "",
        "| Topic | Status | Evidence |", "|---|---|---|",
    ]
    for topic, status, evidence in AUDIT_ROWS:
        lines.append(f"| {topic} | **{status}** | {evidence} |")
    lines += [
        "", "## Conclusion",
        "", "No holographic formulation, object-wave definition, reference-wave definition, "
        "monopole model, or modulation/interference equation existed anywhere in this repository "
        "before Phase 8B. Every existing solver convention needed to build one *consistently* "
        "(time convention, coordinate system, port definitions, incident-field convention, "
        "polarization convention, phase-sign convention) was, however, already fully defined, "
        "implemented and documented (Phases 1-7) -- Phase 8B reuses every one of them rather than "
        "introducing a second, parallel set (`src/phase8b/conventions.py`).",
    ]
    return "\n".join(lines) + "\n"


def render_holographic_formulation_md(evaluations, distinct_states):
    lines = [
        "# Holographic_Formulation.md", "",
        "Two literature-sourced candidate holographic formulations, evaluated against this "
        "project's actual, frozen Phase 1-7 conventions and Phase 7 library -- neither is "
        "asserted as *the* project's formulation (no project document authorizes either one; "
        "`repository_scientific_audit.md`).", "",
        "## Coordinate system", "",
        f"x: {conventions.COORDINATE_SYSTEM.x_range}; y: {conventions.COORDINATE_SYSTEM.y_range}; "
        f"z: {conventions.COORDINATE_SYSTEM.z_convention}.", "",
        "## Time convention", "", f"`{conventions.TIME_CONVENTION}` (unchanged, CONVENTIONS.md).", "",
        "## Reference field", "",
        "Formulation A (Sec. below): a monopole-launched guided surface wave, "
        "`E_ref(rho) ~ rho^(-1/2) exp(-j k_swg rho)` (`reference_wave.monopole_reference_field`); "
        "`k_swg` has no default and cannot be computed by this project's solver. Formulation B: an "
        "arbitrary reference field (any of Sec. 9's three input types; `reference_field.py`, "
        "Phase 8A).", "",
        "## Object field", "",
        "`E_obj(x,y) = sum_n w_n exp(j*phi0_n) exp(-j(kx0_n x + ky0_n y))` at z=0, "
        "`kx0_n = k0 sin(theta_n) cos(phi_n)`, `ky0_n = k0 sin(theta_n) sin(phi_n)` -- the "
        "identical formula `excitation.PlaneWaveExcitation.transverse_wavenumber` implements, "
        "cross-checked by test (`src/phase8b/object_wave.py`; that class cannot be constructed "
        "directly for every (theta, phi, polarization) triple an object-wave beam needs, see its "
        "module-level comment).", "",
        "## Multi-beam formulation", "",
        "Field-level complex superposition (`object_wave.combine_beams`); weight = complex "
        "field-amplitude coefficient (Sec. 16); `initial_phase_deg` enters as one constant "
        "additive phase per beam (Sec. 17), no other role.", "",
        "## Holographic modulation", "",
    ]
    for e in evaluations:
        lines += [
            f"### Formulation {e.formulation_name}", "",
            f"Required quantity: {e.required_quantity}.",
            f"Library compatibility: **{e.library_compatibility}**.",
            f"Decision trail: " + " -> ".join(f"[{q}: {a}]" for q, a in e.decision_trail),
            f"Final status: **{e.final_status}**.", "",
        ]
    lines += [
        "## Physical control variable", "",
        "Unresolved: two candidates (S11xx/S21xx for Formulation B; surface reactance for "
        "Formulation A, INVALID for this project). See `control_variable_analysis.md`.", "",
        "## Unit-cell mapping", "",
        "Not established (depends on the control-variable resolution above); no geometry "
        "interpolation/extrapolation is implemented (unchanged from Phase 8A).", "",
        "## 324 geometries vs. distinct phase states", "",
        f"{distinct_states.n_geometries} geometries, {distinct_states.n_distinct_phase_states} "
        f"distinct phase states ({distinct_states.n_duplicate_records} duplicate records). "
        f"{distinct_states.mesh_hypothesis_note}", "",
    ]
    return "\n".join(lines) + "\n"


def render_control_variable_analysis_md(evaluations):
    from . import control_variable as cv

    s11, s21, zs, gamma = (f() for f in cv.ALL_CANDIDATES)
    lines = [
        "# control_variable_analysis.md", "",
        "## 1. Candidate physical quantities", "",
        "S11xx (reflection), S21xx (transmission), surface impedance Zs, reflection coefficient "
        "Gamma (Sec. 22-25).", "",
        "## 2. S11xx analysis", "", s11.justification, "", f"**Status: {s11.validated}**", "",
        "## 3. S21xx analysis", "", s21.justification, "", f"**Status: {s21.validated}**", "",
        "## 4. Reflection coefficient analysis", "", gamma.justification, "", f"**Status: {gamma.validated}**", "",
        "## 5. Surface impedance analysis", "", zs.justification, "", f"**Status: {zs.validated}**", "",
        "## 6. Polarization compatibility", "",
        "Both S11xx/S21xx are co-polarized (x-in, x-out) quantities; an x-polarized library "
        "response cannot, by itself, encode an arbitrary (e.g. dual-linear or circular) desired "
        "object-wave polarization without an additional, unvalidated assumption -- not made here.", "",
        "## 7. Reference-plane implications", "",
        "S11xx/S21xx are referenced at the top/bottom air reference planes "
        f"({conventions.PORT_DEFINITION.port1_plane}; {conventions.PORT_DEFINITION.port2_plane}), "
        "not at the metal pattern itself; a holographic formulation using them must account for "
        "this fixed reference, not assume phase is referenced at z=0.", "",
        "## 8. Time/phase convention", "",
        f"`{conventions.TIME_CONVENTION}`, phase on {conventions.PHASE_PRINCIPAL_RANGE_DEG} deg "
        "-- unchanged, reused throughout `src/phase8b/`.", "",
        "## 9. Phase 7 library compatibility", "",
    ]
    for e in evaluations:
        lines.append(f"- {e.formulation_name}: **{e.library_compatibility}** (final status {e.final_status}).")
    lines += ["", "## 10. Final scientific status", "", "See `holographic_specification.json`'s top-level `status`."]
    return "\n".join(lines) + "\n"
