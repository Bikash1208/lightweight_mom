"""Phase 8C Sec. 23-28: the reflection-architecture and transmission-architecture tests, and the
per-formulation Phase 8B re-evaluation -- reusing `phase8b.control_variable`/`library_compatibility`
(unmodified) rather than re-deriving the same conclusions independently, which risks silently
drifting from Phase 8B's own reasoning (Sec. 5's explicit instruction to treat Phase 8B as
authoritative).
"""

from dataclasses import dataclass

from ..phase8b import control_variable as cv
from ..phase8b import library_compatibility as lib_compat


@dataclass(frozen=True)
class ArchitectureTestResult:
    hypothesis: str
    chain: tuple            # ordered (link, status) pairs, Sec. 24-25's explicit chain
    status: str              # "...VALIDATED" | "...UNVALIDATED"


def reflection_architecture_test():
    """Sec. 24: reference formulation -> required reflection quantity -> S11 interpretation ->
    Phase 7 response."""
    formulation_b = lib_compat.evaluate_formulation_b()
    s11 = cv.investigate_s11xx()
    chain = (
        ("reference formulation defines a reflection-phase control variable (Formulation B)", "YES"),
        ("required reflection quantity = normal-incidence co-pol reflection phase", "DEFINED"),
        ("S11xx interpretation as that quantity", s11.validated),
        ("Phase 7 response available (S11 stored in full)", "YES"),
    )
    status = "REFLECTION_ARCHITECTURE_VALIDATED" if s11.validated == "VALIDATED" else "REFLECTION_ARCHITECTURE_UNVALIDATED"
    return ArchitectureTestResult(hypothesis="Holographic control variable = reflection coefficient (S11xx)", chain=chain, status=status)


def transmission_architecture_test():
    """Sec. 25: reference formulation -> required transmission quantity -> S21 interpretation ->
    Phase 7 response."""
    s21 = cv.investigate_s21xx()
    chain = (
        ("reference formulation defines a transmission-phase control variable (Formulation B)", "YES"),
        ("required transmission quantity = normal-incidence co-pol transmission phase", "DEFINED"),
        ("S21xx interpretation as that quantity", s21.validated),
        ("Phase 7 response available (S21 stored in full)", "YES"),
    )
    status = "TRANSMISSION_ARCHITECTURE_VALIDATED" if s21.validated == "VALIDATED" else "TRANSMISSION_ARCHITECTURE_UNVALIDATED"
    return ArchitectureTestResult(hypothesis="Holographic control variable = transmission coefficient (S21xx)", chain=chain, status=status)


def formulation_a_reevaluation():
    """Sec. 27: re-state (not re-derive) Phase 8B's Formulation A conclusion, with every
    requirement it lists, unmodified unless new evidence exists (none does)."""
    e = lib_compat.evaluate_formulation_a()
    return {
        "ground_plane_requirement": "grounded structure required (Fong et al. 2010)",
        "surface_wave_requirement": "guided surface-wave reference field required",
        "surface_reactance_requirement": "X(x,y) of a grounded structure",
        "reference_wave_requirement": "monopole-launched cylindrical guided surface wave",
        "object_wave_requirement": "arbitrary desired radiated beam(s)",
        "modulation_equation": "Psi = Re{E_ref . conj(E_obj)}, X = X_avg(1+M*Psi)",
        "unit_cell_response_requirement": "surface reactance of a grounded unit cell",
        "status": e.final_status,   # unchanged from Phase 8B: LIBRARY_CONTROL_VARIABLE_INCOMPATIBLE
    }


def formulation_b_reevaluation():
    """Sec. 28: does Formulation B correspond to reflection phase, transmission phase, or
    neither? Both remain open candidates (S11xx/S21xx both INSUFFICIENT_EVIDENCE, Phase 8B) --
    Formulation B is not "validated" merely because it resembles a reflectarray equation."""
    reflection = reflection_architecture_test()
    transmission = transmission_architecture_test()
    return {
        "corresponds_to_reflection_phase": reflection.status,
        "corresponds_to_transmission_phase": transmission.status,
        "note": "Formulation B (a phase-difference/phase-conjugate hologram) is mathematically "
                "compatible with either a reflection-phase or a transmission-phase control "
                "variable; this project has not determined which architecture (reflective vs. "
                "transmissive MTS) it intends, so neither correspondence is validated (Sec. 28's "
                "explicit prohibition on calling it 'validated' by resemblance alone).",
    }
