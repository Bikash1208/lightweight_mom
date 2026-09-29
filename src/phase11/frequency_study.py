"""Sec. 14 (Phase 11J): does the coupling extent change across the operating band? A reduced
frequency grid (5 points spanning the project's 10-15 GHz band, Sec. 14's own allowance for "an
appropriately reduced frequency set if computational cost requires it") repeats one representative
coupling measurement -- the 3x3-vs-9x9 (checkerboard) shell increment from `coupling_radius.py` --
at each frequency, for the primary center.
"""

from dataclasses import dataclass

from .coupling_radius import compute_shell_increments, run_coupling_radius_series

FREQUENCIES_HZ = (10.0e9, 11.25e9, 12.5e9, 13.75e9, 15.0e9)  # 5 points, band edges + quarters + design


@dataclass(frozen=True)
class FrequencyPoint:
    frequency_hz: float
    center_name: str
    responses_by_n: dict  # n -> complex
    max_shell_increment_relative: float


def run_frequency_study(
    center_name,
    center_lw_mm,
    period_m,
    substrate,
    metal,
    mesh_cells_per_subcell=8,
    frequencies_hz=FREQUENCIES_HZ,
    sizes=(1, 3, 5, 9),
):
    results = []
    for f in frequencies_hz:
        points = run_coupling_radius_series(
            center_lw_mm, period_m, substrate, metal, f, mesh_cells_per_subcell, sizes=sizes
        )
        increments = compute_shell_increments(points)
        max_rel = max(inc.delta_relative for inc in increments)
        results.append(FrequencyPoint(f, center_name, {p.n: p.response for p in points}, max_rel))
    return results
