"""Sec. 9 (Phase 11E): "reproduce or independently re-establish" Phase 10's own diagonal-coupling
finding (25-70% of the cardinal contribution). Literal reuse -- the exact same functions
(`phase10.neighbor_sensitivity`, `phase10.diagonal_study`), not a re-implementation -- since that
finding was already established via a real, correct MoM computation and Sec. 9 only asks that it be
confirmed still holds, not that the method be redone from scratch.
"""

from ..phase10 import diagonal_study as p10_diag
from ..phase10 import neighbor_sensitivity as p10_ns
from ..phase10.environment_geometry import strong_neighbor_low_high


def reproduce_diagonal_finding(center_lw_mm, period_m, substrate, metal, frequency_hz,
                                mesh_cells_per_subcell=8):
    low, high = strong_neighbor_low_high(center_lw_mm)
    contrast = (low, low) if abs(center_lw_mm[0] - high) < abs(center_lw_mm[0] - low) else (high, high)

    five_cell_cases = p10_ns.run_five_cell_test_matrix(center_lw_mm, contrast, period_m, substrate,
                                                        metal, frequency_hz, mesh_cells_per_subcell)
    cardinal_contribs = p10_ns.neighbor_contributions(five_cell_cases)
    s0 = five_cell_cases["5A_homogeneous"].px

    diagonal_cases = p10_diag.run_diagonal_single_tests(center_lw_mm, contrast, period_m, substrate,
                                                         metal, frequency_hz, mesh_cells_per_subcell)
    diagonal_contribs = p10_diag.diagonal_contributions(diagonal_cases, s0)

    cardinal_mean = sum(c.delta_s_abs for c in cardinal_contribs.values()) / len(cardinal_contribs)
    diagonal_mean = sum(c.delta_s_abs for c in diagonal_contribs.values()) / len(diagonal_contribs)
    ratio = diagonal_mean / cardinal_mean if cardinal_mean > 0 else float("inf")

    return {
        "cardinal_mean": cardinal_mean, "diagonal_mean": diagonal_mean, "ratio": ratio,
        "reproduces_phase10_range": 0.20 <= ratio <= 0.80,
        "cardinal_contributions": {k: v.delta_s_abs for k, v in cardinal_contribs.items()},
        "diagonal_contributions": {k: v.delta_s_abs for k, v in diagonal_contribs.items()},
    }
