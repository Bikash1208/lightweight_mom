"""Phase 8A Sec. 3-5: the holographic phase-source validation gate.

Every one of Sec. 4's six checks requires comparing the library's phase-control quantity against
*the physical quantity the project's holographic generator actually requires*. This project has
no implemented holographic generator to ask -- verified by inspection, not assumed:

    README.md:26        "Holographic synthesis is **not implemented**."
    DATABASE.md:277      "Not implemented (by request): holographic synthesis, cell selection/
                          matching, interpolation of the database."
    src/__init__.py:20   "NOT implemented yet: holographic synthesis, oblique-incidence ports."

and confirmed by an exhaustive search: no `aperture`/`beam`/`reference_field`/`hologram*` module
exists anywhere in `src/`. `phase_extraction.py`'s own `resolve_phase_source()` (Phase 7) is the
closest thing this project has to an opinion on the question, and it explicitly declines to
answer it (`status="SPECIFICATION_MANDATED"`, not `"PHYSICALLY_VALIDATED"`).

Consequently Sec. 4 points 1 ("what quantity the holographic generator requires") and 3 ("whether
the two are physically equivalent") are *structurally* unanswerable today, not merely uncertain.
This module says so rather than assuming S11xx is adequate merely because Phase 7 extracted it
(brief Sec. 5's explicit prohibition).
"""

from dataclasses import dataclass

from .phase_extraction import resolve_phase_source

# Literal fact, verified as above; flipping it to True would be a lie until a holographic
# formulation actually exists in this project. `feasibility_required_phase.py` guards on the same
# constant so the two modules cannot drift out of sync with each other.
HOLOGRAPHIC_GENERATOR_IMPLEMENTED = False


@dataclass(frozen=True)
class PhaseSourceValidation:
    quantity_required: str | None  # what a holographic generator would need; None = unknown
    quantity_provided: str  # what the library's pinned phase source actually is
    physically_equivalent: bool | None  # None = cannot be assessed (quantity_required unknown)
    coordinate_convention_ok: bool | None
    polarization_convention_ok: bool | None
    sign_convention_ok: bool | None
    time_convention_ok: bool | None
    status: str  # "PHASE_SOURCE_UNVALIDATED" | "PHYSICALLY_VALIDATED"
    note: str


def validate_holographic_phase_source():
    """Sec. 4's six-point gate, run against the actual project state. Returns
    `status="PHASE_SOURCE_UNVALIDATED"` for this project today -- the honest, expected result
    (brief Sec. 48), not a bug to be worked around."""
    library_resolution = resolve_phase_source()  # phase_extraction.py, Phase 7, unmodified

    if HOLOGRAPHIC_GENERATOR_IMPLEMENTED:
        # Not reached today; left explicit so a future implementer cannot silently keep this
        # module's PHASE_SOURCE_UNVALIDATED verdict after adding a real generator.
        raise NotImplementedError(
            "a holographic generator now exists in this project but "
            "feasibility_phase_source.validate_holographic_phase_source was not updated to "
            "compare its required quantity against the library's phase source"
        )

    return PhaseSourceValidation(
        quantity_required=None,
        quantity_provided=library_resolution.component,
        physically_equivalent=None,
        coordinate_convention_ok=None,
        polarization_convention_ok=None,
        # Only one time convention is defined anywhere in this project (exp(+j*omega*t),
        # CONVENTIONS.md / constants.py TIME_CONVENTION); there is no second, independently
        # defined holographic-side convention for it to disagree with, so this one point is
        # trivially satisfied -- stated explicitly, not silently assumed satisfied in general.
        sign_convention_ok=None,
        time_convention_ok=True,
        status="PHASE_SOURCE_UNVALIDATED",
        note=(
            "This project has no implemented holographic-synthesis formulation (README.md "
            "'Holographic synthesis is not implemented'; DATABASE.md; src/__init__.py), so the "
            "physical quantity a holographic generator would require cannot be determined, and "
            "therefore cannot be checked against what the unit-cell library provides. The "
            f"library itself pins its stored phase-vs-geometry data to "
            f"{library_resolution.component} (status={library_resolution.status}): "
            f"{library_resolution.note} Both facts together mean {library_resolution.component} "
            "must not be treated as 'the' holographic phase-control quantity without further "
            "architecture decisions this project has not made."
        ),
    )
