"""Phase 12A final status: an explicit code from a fixed set (deliverables list), never a bare
pass/fail, never silently upgraded. Geometry inversion and STEP export stay LOCKED
unconditionally -- this phase does not even attempt to unlock them (the brief's own scope).
"""

from dataclasses import dataclass

from ..phase9.hfss_reference import DEFAULT_HFSS_REFERENCE_PATH, load_hfss_reference

GROUNDED_STACK_VALIDATED = "GROUNDED_STACK_VALIDATED"
GROUNDED_STACK_PARTIAL = "GROUNDED_STACK_PARTIAL"
HFSS_VALIDATED = "HFSS_VALIDATED"
HFSS_NOT_AVAILABLE = "HFSS_NOT_AVAILABLE"
GEOMETRY_INVERSION_LOCKED = "GEOMETRY_INVERSION_LOCKED"
STEP_EXPORT_LOCKED = "STEP_EXPORT_LOCKED"

# Analytic bare-slab check: the closed form is derived algebra, so any real disagreement is a
# bug, not numerical noise -- this tolerance is generous compared to the ~1e-14 actually measured.
ANALYTIC_ERROR_TOLERANCE = 1e-9
# Same policy value as discretization.DEFAULT_REFINEMENT_TOLERANCE (reused, not re-invented) --
# a successive-difference bound on |S|, not an accuracy guarantee (discretization.py docstring).
MESH_FLOQUET_TOLERANCE = 0.02


@dataclass(frozen=True)
class GateResult:
    status: str  # GROUNDED_STACK_VALIDATED | GROUNDED_STACK_PARTIAL
    reason: str
    measured_phase_coverage_deg: float  # S11_xx circular coverage of the physics-sanity sweep
    hfss_status: str
    geometry_inversion_status: str
    step_export_status: str


def hfss_status(path=DEFAULT_HFSS_REFERENCE_PATH):
    """Reuses phase9.hfss_reference (unmodified): NOT_AVAILABLE until a real HFSS export exists
    anywhere in this project -- never fabricated, never upgraded on the strength of this phase's
    own MoM-only results."""
    return HFSS_VALIDATED if load_hfss_reference(path) is not None else HFSS_NOT_AVAILABLE


def decide(max_analytic_error, floquet_error, mesh_error, any_surface_wave_in_band, coverage_deg):
    """The Phase 12A status, from measured diagnostics only (no criterion is chosen after seeing
    the final headline number; these three are the brief's own explicit validation items)."""
    reasons = []
    if max_analytic_error > ANALYTIC_ERROR_TOLERANCE:
        reasons.append(
            f"bare grounded-slab reflection disagrees with the textbook closed form by "
            f"{max_analytic_error:.3g} (> {ANALYTIC_ERROR_TOLERANCE:g})"
        )
    if any_surface_wave_in_band:
        reasons.append(
            "a Floquet harmonic falls inside the TM0 surface-wave interval for at least one "
            "tested (period, frequency) combination"
        )
    if floquet_error > MESH_FLOQUET_TOLERANCE or mesh_error > MESH_FLOQUET_TOLERANCE:
        reasons.append(
            f"the representative grounded patch's mesh/Floquet refinement study exceeded the "
            f"{MESH_FLOQUET_TOLERANCE:g} successive-difference tolerance "
            f"(floquet_error={floquet_error:.3g}, mesh_error={mesh_error:.3g}) -- the Phase 7 "
            f"mesh choice, reused here for comparability, resolves this grounded/near-resonant "
            f"patch less precisely than the ungrounded sweep it was chosen for"
        )
    status = GROUNDED_STACK_VALIDATED if not reasons else GROUNDED_STACK_PARTIAL
    reason = (
        "bare-slab closed form, surface-wave-pole avoidance and mesh/Floquet convergence all passed within tolerance"
        if not reasons
        else "; ".join(reasons)
    )
    return GateResult(
        status=status,
        reason=reason,
        measured_phase_coverage_deg=coverage_deg,
        hfss_status=hfss_status(),
        geometry_inversion_status=GEOMETRY_INVERSION_LOCKED,
        step_export_status=STEP_EXPORT_LOCKED,
    )
