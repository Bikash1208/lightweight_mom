"""Sec. 10 (Phase 11F): every required controlled environment, solved directly by MoM (no ML,
no interpolation of unmeasured cases). Sec. 18 Case D (pattern-dependent behavior) is tested here:
several of these patterns share the same two geometries (A, B) and even the same 50/50 mix, but
arrange them completely differently -- if the center response differs substantially between them
despite identical local statistics, local descriptors are insufficient regardless of what the
coupling-radius series alone found.
"""

from dataclasses import dataclass

from . import patterns as pt
from .reference_solver import solve_grid

PATTERN_N = 7   # Sec. 10's patterns are 5-7 cells wide in their own ASCII diagrams; 7x7 keeps every
                # solve under ~4 s (measured) while giving each pattern a real, non-trivial extent.


@dataclass(frozen=True)
class PatternResult:
    name: str
    response: complex
    n_unknowns: int


def _solve(grid_mm, period_m, substrate, metal, frequency_hz, mesh_cells_per_subcell, n):
    grid_m = [[(l * 1e-3, w * 1e-3) for (l, w) in row] for row in grid_mm]
    sg = solve_grid(n, period_m, substrate, metal, grid_m, frequency_hz, mesh_cells_per_subcell)
    mid = n // 2
    return sg.response_at(mid, mid), sg.n_unknowns


def run_pattern_study(center_lw_mm, contrast_lw_mm, period_m, substrate, metal, frequency_hz,
                       mesh_cells_per_subcell=8, n=PATTERN_N):
    """Sec. 10's required patterns, all at the SAME size/center/contrast pair so they are directly
    comparable to each other (only the spatial ARRANGEMENT differs, isolating Sec. 18 Case D)."""
    cases = {
        "uniform": pt.uniform(n, center_lw_mm),
        "single_perturbation": pt.single_perturbation(n, center_lw_mm, contrast_lw_mm),
        "smooth_gradient": pt.smooth_gradient(n, center_lw_mm[0], contrast_lw_mm[0], center_lw_mm[1]),
        "abrupt_boundary": pt.abrupt_boundary(n, center_lw_mm, contrast_lw_mm),
        "alternating": pt.alternating(n, center_lw_mm, contrast_lw_mm),
        "checkerboard": pt.checkerboard(n, center_lw_mm, contrast_lw_mm),
        "strong_isolated_defect": pt.strong_isolated_defect(n, center_lw_mm),
    }
    results = {}
    for name, grid_mm in cases.items():
        response, n_unknowns = _solve(grid_mm, period_m, substrate, metal, frequency_hz,
                                       mesh_cells_per_subcell, n)
        results[name] = PatternResult(name, response, n_unknowns)
    return cases, results


def detect_pattern_dependence(results, substantial_relative=0.5):
    """Sec. 18 Case D: compares ``alternating`` and ``checkerboard`` (identical local 50/50
    A/B statistics, and in this project's own patterns.py literally the same arrangement -- so this
    is primarily a self-consistency check) against ``abrupt_boundary`` (same two geometries, same
    count in a 7-wide row, but a single contiguous block instead of spread out) to see whether
    arrangement alone (holding the SET of geometries present fixed) changes the center response by
    a substantial amount."""
    if "alternating" not in results or "abrupt_boundary" not in results:
        return False, {}
    base = results["checkerboard"].response
    comparisons = {}
    for name in ("alternating", "abrupt_boundary"):
        other = results[name].response
        rel = abs(other - base) / abs(base) if abs(base) > 0 else float("inf")
        comparisons[name] = rel
    detected = any(v >= substantial_relative for v in comparisons.values())
    return detected, comparisons
