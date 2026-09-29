"""Sec. 18-19, 21: what happens when a query finds nothing. The system reports
`DATABASE_COVERAGE_MISSING` and hands back a precise, machine-readable specification of the
library that *would* answer the query -- it never silently builds it, rescales an existing
library's geometry, or interpolates between unrelated periods (Sec. 21's explicit prohibitions).
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class LibraryBuildRequest:
    topology: str | None
    p_x_mm: float | None
    p_y_mm: float | None
    material_stack: dict  # {name, eps_r, tan_delta, thickness_mm} or {} if unspecified
    frequency_start_hz: float | None
    frequency_stop_hz: float | None
    frequency_points: int | None
    geometry_parameter_ranges: dict  # {} unless the query itself constrained geometry
    polarization: str | None
    incidence: str | None
    solver_settings: dict  # left mostly empty; this package does not choose mesh/Floquet order
    reason: str  # which filter stage eliminated every candidate
    status: str = "BUILD_REQUIRED"  # Sec. 18: "BUILD_REQUIRED", never "silently generated"


def generate_build_request(query_result):
    """Sec. 19: a `LibraryBuildRequest` from a `query.QueryResult` whose status is
    `DATABASE_COVERAGE_MISSING`. Raises ValueError if called on a result that actually matched
    (there is nothing to request)."""
    if query_result.status != "DATABASE_COVERAGE_MISSING":
        raise ValueError("generate_build_request: query_result was not DATABASE_COVERAGE_MISSING")
    q = query_result.query
    material_stack = {}
    if (
        q.substrate_name is not None
        or q.eps_r is not None
        or q.tan_delta is not None
        or q.substrate_thickness_mm is not None
    ):
        material_stack = {
            "name": q.substrate_name,
            "eps_r": q.eps_r,
            "tan_delta": q.tan_delta,
            "thickness_mm": q.substrate_thickness_mm,
        }
    return LibraryBuildRequest(
        topology=q.topology,
        p_x_mm=q.p_x_mm,
        p_y_mm=q.p_y_mm,
        material_stack=material_stack,
        frequency_start_hz=q.frequency_hz,
        frequency_stop_hz=q.frequency_hz,
        frequency_points=1 if q.frequency_hz is not None else None,
        geometry_parameter_ranges={},
        polarization=q.polarization,
        incidence=q.incidence,
        solver_settings={},
        reason=f"no registered library entry survived the '{query_result.eliminated_at}' filter stage "
        f"(filter trail: {query_result.filter_trail})",
    )
