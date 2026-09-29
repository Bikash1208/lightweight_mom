"""Sec. 11-12: does the real, measured unit-cell architecture actually satisfy each candidate
formulation's assumptions? Reuses `phase8c.topology_audit`/`phase8c.power_balance` (unmodified) --
the ground-plane/architecture-classification/cross-polarization facts are not re-derived here.
"""

from dataclasses import dataclass

ARCHITECTURE_COMPATIBILITY_STATES = (
    "ARCHITECTURE_COMPATIBLE_UNVALIDATED",  # topology does not rule it out; correspondence still unvalidated
    "ARCHITECTURE_INCOMPATIBLE",  # topology structurally cannot satisfy the requirement
)


@dataclass(frozen=True)
class ArchitectureCompatibilityResult:
    formulation: str
    required_architecture: tuple
    status: str
    evidence: str


def assess_formulation_a(topology_audit):
    """Candidate A needs a grounded surface (Sec. 11); `topology_audit.ground_plane` is measured
    `ABSENT` (Phase 8C, reused)."""
    incompatible = topology_audit.ground_plane != "PRESENT"
    return ArchitectureCompatibilityResult(
        formulation="A",
        required_architecture=("grounded surface", "surface-wave excitation"),
        status="ARCHITECTURE_INCOMPATIBLE" if incompatible else "ARCHITECTURE_COMPATIBLE_UNVALIDATED",
        evidence=f"topology_audit.ground_plane = {topology_audit.ground_plane!r}: {topology_audit.ground_plane_evidence}",
    )


def assess_formulation_b(mean_power):
    """Candidate B needs a nonzero reflective response (Sec. 11: 'reflective boundary')."""
    r_total = mean_power["mean_R_co"] + mean_power["mean_R_cross"]
    compatible = r_total > 0.0
    return ArchitectureCompatibilityResult(
        formulation="B",
        required_architecture=("reflective boundary (nonzero reflection response)",),
        status="ARCHITECTURE_COMPATIBLE_UNVALIDATED" if compatible else "ARCHITECTURE_INCOMPATIBLE",
        evidence=f"library-mean reflected power (co+cross) = {r_total:.6g} (phase8c.power_balance.library_power_summary)",
    )


def assess_formulation_c(mean_power):
    """Candidate C needs a nonzero transmissive response (Sec. 11: 'transmissive boundary')."""
    t_total = mean_power["mean_T_co"] + mean_power["mean_T_cross"]
    compatible = t_total > 0.0
    return ArchitectureCompatibilityResult(
        formulation="C",
        required_architecture=("transmissive boundary (nonzero transmission response)",),
        status="ARCHITECTURE_COMPATIBLE_UNVALIDATED" if compatible else "ARCHITECTURE_INCOMPATIBLE",
        evidence=f"library-mean transmitted power (co+cross) = {t_total:.6g} (phase8c.power_balance.library_power_summary)",
    )


def assess_formulation_d(mean_power, cross_pol_threshold=1e-6):
    """Candidate D needs a polarization-converting (non-diagonal) response (Sec. 11); Phase 8C
    measured negligible cross-polarization (<1e-6 of incident power) across the whole registered
    library -- a direct, measured architecture/library mismatch, not an assumption."""
    r_cross, t_cross = mean_power["mean_R_cross"], mean_power["mean_T_cross"]
    negligible = max(r_cross, t_cross) < cross_pol_threshold
    return ArchitectureCompatibilityResult(
        formulation="D",
        required_architecture=("polarization-converting (non-diagonal) scattering response",),
        status="ARCHITECTURE_INCOMPATIBLE" if negligible else "ARCHITECTURE_COMPATIBLE_UNVALIDATED",
        evidence=f"library-mean cross-polarized power: R_cross={r_cross:.3g}, T_cross={t_cross:.3g} "
        f"(threshold {cross_pol_threshold:.0e}) -- this unrotated rectangular-patch family "
        "does not realize a measurable polarization-converting response.",
    )


def assess_all(topology_audit, mean_power):
    return {
        "A": assess_formulation_a(topology_audit),
        "B": assess_formulation_b(mean_power),
        "C": assess_formulation_c(mean_power),
        "D": assess_formulation_d(mean_power),
    }
