"""Phase 8B Sec. 33-36, 45: the machine-readable outputs -- `holographic_specification.json`,
`source_registry.json`, `phase_control_quantity.csv`. Every field is populated from the actual
computed/investigated objects in this package; nothing here hard-codes a conclusion (Sec. 45)."""

import json
from pathlib import Path

from . import conventions
from .library_compatibility import evaluate_all


def build_source_registry():
    """Sec. 35-36: every external claim this phase relies on, with citation, the exact
    equation/section being adopted, why it is (or is not) compatible with this project's
    conventions, and its validation status -- research-reproducibility traceability."""
    return {
        "fong2010": {
            "claim": "Holographic interference-modulation index Psi(x,y) = Re{E_ref . conj(E_obj)}, "
                     "used to modulate a grounded surface's reactance for a leaky-wave antenna fed "
                     "by a monopole-launched guided surface wave.",
            "source": "T. Fong, N. C. Colburn, J. J. Ottusch, J. L. Visher, D. F. Sievenpiper, "
                      "\"Scalar and Tensor Holographic Artificial Impedance Surfaces,\" IEEE Trans. "
                      "Antennas Propag., vol. 58, no. 10, pp. 3212-3221, 2010.",
            "equation_section": "Interference/modulation-index formulation (exact page/equation "
                                "number not independently re-verified against the paywalled full "
                                "text by this project; the general form Psi=Re{E_ref.conj(E_obj)}, "
                                "X=X_avg(1+M*Psi) is the widely-cited, standard restatement of this "
                                "paper's technique across the holographic-impedance-surface "
                                "literature -- see repository_scientific_audit.md's search log).",
            "project_convention_compatibility": "Time convention exp(+j*omega*t) is compatible "
                                                 "with the general interference-kernel form used "
                                                 "here (a real-valued correlation of two complex "
                                                 "phasors is convention-symmetric); the *reference-"
                                                 "wave* half (a monopole-launched guided surface "
                                                 "wave on a grounded structure) is structurally "
                                                 "incompatible with this project's solver (no "
                                                 "ground plane, no surface-wave dispersion solver).",
            "implementation_consequence": "src/phase8b/holographic_modulation.py::interference_modulation "
                                          "(kernel, implemented+tested); "
                                          "reactance_from_modulation (X(x,y) mapping, implemented "
                                          "but its M/X_avg calibration constants are NOT available "
                                          "from this project).",
            "validation_status": "MATHEMATICALLY DEFINED, NUMERICALLY IMPLEMENTED, NUMERICALLY "
                                 "TESTED. NOT physically validated for this project "
                                 "(control_variable.py: SURFACE_IMPEDANCE = INVALID for this "
                                 "library).",
        },
        "impedance_boundary_relation": {
            "claim": "Zs = eta0 (1+Gamma)/(1-Gamma): the surface impedance of a boundary whose "
                     "normal-incidence reflection coefficient (referenced at the boundary) is Gamma.",
            "source": "Standard transmission-line/impedance-boundary-condition relation "
                      "(load-reflection-coefficient formula Gamma=(Z_L-Z_0)/(Z_L+Z_0) applied "
                      "with Z_0=eta0, Z_L=Zs); textbook EM/microwave engineering (e.g. D. M. "
                      "Pozar, Microwave Engineering; C. A. Balanis, Advanced Engineering "
                      "Electromagnetics), used throughout the artificial-impedance-surface/AMC "
                      "literature including Sievenpiper's own earlier high-impedance-surface work.",
            "equation_section": "Zs = eta0 (1+Gamma)/(1-Gamma), Gamma referenced at the boundary "
                                "plane (not de-embedded to any other reference plane).",
            "project_convention_compatibility": "Requires |Gamma|=1 (a fully reflective, lossless "
                                                "grounded boundary) to give a purely reactive Zs. "
                                                "Phase 7's library is NOT grounded and is mostly "
                                                "transmissive (see control_variable.py's case01 "
                                                "citation), so this relation's assumption is "
                                                "violated for this project's data.",
            "implementation_consequence": "Not implemented as a conversion function in this "
                                          "package (deliberately -- Sec. 10's prohibition on "
                                          "assuming Zs<->phase equivalence without a validated "
                                          "derivation); only documented and cited here.",
            "validation_status": "Standard, well-established EM relation in general. INVALID as "
                                 "applied to this project's ungrounded, partially-transmissive "
                                 "unit-cell library (control_variable.py::investigate_surface_impedance).",
        },
        "phase_difference_hologram": {
            "claim": "required_phase(x,y) = wrap(arg(E_obj) - arg(E_ref)): the phase a passive "
                     "reflectarray/transmitarray unit cell must impose on a locally-incident "
                     "reference wave to produce a locally-desired object wave.",
            "source": "Standard phase-conjugate/phase-only hologram and reflectarray/transmitarray "
                      "antenna design equation (widely used across passive reflectarray and "
                      "coding-metasurface antenna design literature; not attributed to one single "
                      "paper here, since it is a standard, convention-independent construction).",
            "equation_section": "required_phase = wrap(arg(E_obj) - arg(E_ref))",
            "project_convention_compatibility": "Fully compatible: needs only a normal-incidence "
                                                "co-polarized reflection/transmission phase per "
                                                "cell, exactly the *kind* of quantity Phase 7 "
                                                "solved and stored (S11xx/S21xx).",
            "implementation_consequence": "src/phase8b/holographic_modulation.py::phase_difference_hologram "
                                          "(implemented+tested); can feed directly into "
                                          "`feasibility_matching.match_required_phase` (Phase 8A).",
            "validation_status": "MATHEMATICALLY DEFINED, NUMERICALLY IMPLEMENTED, NUMERICALLY "
                                 "TESTED. NOT physically validated as this project's chosen "
                                 "formulation (no project document picks it over Formulation A, "
                                 "and the S11-vs-S21 architecture question remains open).",
        },
    }


def build_holographic_specification():
    """Sec. 45's schema, extended with a `formulations` array (both candidates, transparently)
    rather than a single artificially-collapsed answer -- see `library_compatibility.overall_project_status`'s
    docstring for why."""
    evaluations, status = evaluate_all()
    return {
        "status": status,
        "time_convention": conventions.TIME_CONVENTION,
        "coordinate_system": {
            "x": conventions.COORDINATE_SYSTEM.x_range,
            "y": conventions.COORDINATE_SYSTEM.y_range,
            "z": conventions.COORDINATE_SYSTEM.z_convention,
        },
        "reference_wave": {
            "quantity": "monopole-launched guided surface wave (Formulation A) or an arbitrary "
                        "reference field (Formulation B)",
            "definition": "see Holographic_Formulation.md Sec. 'Reference field'",
            "complex_required": True,
        },
        "object_wave": {
            "type": "plane_wave (single- or multi-beam superposition)",
            "definition": "E_obj(x,y) = sum_n w_n exp(j*phi0_n) exp(-j(kx0_n x + ky0_n y)) at z=0",
            "kx": "k0 sin(theta) cos(phi)", "ky": "k0 sin(theta) sin(phi)", "kz": "k0 cos(theta) (unused at z=0)",
        },
        "multi_beam": {"definition": "field-level complex superposition, weight = complex field-amplitude coefficient"},
        "holographic_modulation": {
            "formulations": [
                {
                    "name": e.formulation_name, "required_quantity": e.required_quantity,
                    "library_compatibility": e.library_compatibility, "status": e.final_status,
                    "control_variable": {
                        "name": e.control_variable.name, "validated": e.control_variable.validated,
                    },
                }
                for e in evaluations
            ],
            "validated": False,
        },
        "control_variable": {
            "name": None,
            "physical_definition": "unresolved: see holographic_modulation.formulations for both candidates",
            "validated": False,
        },
        "library_compatibility": {"status": status},
        "hfss_validation": {"performed": False},
    }


def write_json(doc, path):
    path = Path(path)
    path.write_text(json.dumps(doc, indent=1) + "\n")
    return path


def write_phase_control_quantity_csv(evaluations, path):
    """`phase_control_quantity.csv` (Sec. 50): one row per candidate control variable investigated
    (Sec. 22-25), independent of which formulation it belongs to."""
    import csv

    from . import control_variable as cv

    rows = [f() for f in cv.ALL_CANDIDATES]
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["name", "physical_definition", "complex_or_real", "phase_relevant",
                          "magnitude_relevant", "library_quantity", "validated"])
        for r in rows:
            writer.writerow([r.name, r.physical_definition, r.complex_or_real, r.phase_relevant,
                              r.magnitude_relevant, r.library_quantity or "", r.validated])
    return path
