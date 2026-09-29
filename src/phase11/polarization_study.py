"""Sec. 15 (Phase 11K): does spatial coupling depend on excitation polarization? Repeats the same
checkerboard coupling-radius series for x- and y-polarized excitation (the two this project's
solver actually supports, Phase 1's own documented scope -- circular polarization is a DERIVED
quantity in this project's conventions, not an independent excitation state, so it is reported as
NOT_SUPPORTED_AS_INDEPENDENT_EXCITATION rather than silently skipped or fabricated).
"""

from dataclasses import dataclass

from .coupling_radius import compute_shell_increments, run_coupling_radius_series

CIRCULAR_STATUS = "NOT_SUPPORTED_AS_INDEPENDENT_EXCITATION"


@dataclass(frozen=True)
class PolarizationPoint:
    polarization: str
    responses_by_n: dict
    max_shell_increment_relative: float


def run_polarization_study(center_lw_mm, period_m, substrate, metal, frequency_hz,
                            mesh_cells_per_subcell=8, sizes=(1, 3, 5, 9)):
    results = []
    for pol in ("x", "y"):
        points = run_coupling_radius_series(center_lw_mm, period_m, substrate, metal, frequency_hz,
                                             mesh_cells_per_subcell, sizes=sizes, polarization=pol)
        increments = compute_shell_increments(points)
        max_rel = max(inc.delta_relative for inc in increments)
        results.append(PolarizationPoint(pol, {p.n: p.response for p in points}, max_rel))
    return results
