"""Phase 12B final statuses: explicit codes from fixed sets, never a bare pass/fail. Geometry
inversion and STEP export stay LOCKED unconditionally; HFSS status is unchanged (read from
phase12a.gate, not re-derived).
"""

from dataclasses import dataclass

from ..phase12a.gate import HFSS_NOT_AVAILABLE, HFSS_VALIDATED, hfss_status  # noqa: F401 (re-export)

ALIGNED_MESH_VALIDATED = "ALIGNED_MESH_VALIDATED"
ALIGNED_MESH_INVALIDATED = "ALIGNED_MESH_INVALIDATED"

NO_CONCLUSION_CHANGED = "NO_CONCLUSION_CHANGED"
CONCLUSIONS_AT_RISK = "CONCLUSIONS_AT_RISK"

COVERAGE_SUFFICIENT = "COVERAGE_SUFFICIENT"
COVERAGE_INSUFFICIENT = "COVERAGE_INSUFFICIENT"

GEOMETRY_INVERSION_LOCKED = "GEOMETRY_INVERSION_LOCKED"
STEP_EXPORT_LOCKED = "STEP_EXPORT_LOCKED"

# A full 360 deg is "sufficient" by definition; short of that, "sufficient" means close enough to
# be a practical holographic control variable on its own -- a policy threshold, reported alongside
# the actual measured number, never used to hide how far short of 360 deg the true figure is.
COVERAGE_SUFFICIENT_THRESHOLD_DEG = 300.0


@dataclass(frozen=True)
class MeshAlignmentGate:
    status: str
    reason: str
    mesh_error_aligned: float
    floquet_error_aligned: float
    tolerance: float


def decide_mesh_alignment(mesh_error_aligned, floquet_error_aligned, tolerance=0.02):
    ok = mesh_error_aligned <= tolerance and floquet_error_aligned <= tolerance
    status = ALIGNED_MESH_VALIDATED if ok else ALIGNED_MESH_INVALIDATED
    reason = (
        f"aligned-mesh refinement study converged within tolerance "
        f"(mesh_error={mesh_error_aligned:.4g}, floquet_error={floquet_error_aligned:.4g} <= {tolerance:g})"
        if ok
        else f"aligned-mesh refinement study still exceeds tolerance "
        f"(mesh_error={mesh_error_aligned:.4g}, floquet_error={floquet_error_aligned:.4g} > {tolerance:g})"
    )
    return MeshAlignmentGate(status, reason, mesh_error_aligned, floquet_error_aligned, tolerance)


@dataclass(frozen=True)
class AuditGate:
    status: str
    reason: str
    at_risk_findings: tuple


def decide_audit(phase7_duplicate_report, phase10_holdout_leakage, phase10_center_duplicates):
    at_risk = []
    if phase10_center_duplicates:
        at_risk.append(
            f"Phase 10 centers {phase10_center_duplicates} are mesh-identical at the solver's own "
            f"8-cell/side mesh (dx=0.6mm); their isolated-cell baselines coincide, though the "
            f"held-out center is unaffected (leaks={phase10_holdout_leakage['leaks']})"
        )
    status = CONCLUSIONS_AT_RISK if at_risk else NO_CONCLUSION_CHANGED
    reason = (
        "; ".join(at_risk)
        if at_risk
        else "no mesh-identical 'distinct' geometry was found among the traced Phase 9-11 usages"
    )
    return AuditGate(status, reason, tuple(at_risk))


@dataclass(frozen=True)
class PeriodStudyGate:
    status: str
    reason: str
    coverage_by_period_deg: dict  # period_mm -> coverage_deg at design frequency


def decide_period_study(coverage_by_period_deg, threshold_deg=COVERAGE_SUFFICIENT_THRESHOLD_DEG):
    best_period, best_coverage = max(coverage_by_period_deg.items(), key=lambda kv: kv[1])
    sufficient = best_coverage >= threshold_deg
    status = COVERAGE_SUFFICIENT if sufficient else COVERAGE_INSUFFICIENT
    reason = (
        f"best measured coverage is {best_coverage:.1f} deg at period={best_period}mm "
        f"({'>=' if sufficient else '<'} the {threshold_deg:g} deg threshold); "
        + (
            "a single-layer resonant patch reaching this range is usable as a holographic control variable on its own"
            if sufficient
            else "a single-layer resonant patch typically cannot reach a full 360 deg; standard remedies "
            "(thinner substrate, stacked patches, multi-resonant elements) are open options, not "
            "implemented here"
        )
    )
    return PeriodStudyGate(status, reason, dict(coverage_by_period_deg))
