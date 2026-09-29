"""Sec. 13-15: diagonal-neighbor tests, using the full 3x3 ("9-cell") environment
(`environment_geometry.build_9cell_environment`) so cardinal and diagonal positions can be
varied independently and compared directly against `neighbor_sensitivity.py`'s cardinal-only
results (same homogeneous baseline, same contrast levels, same solver settings)."""

from dataclasses import dataclass

from . import environment_geometry as eg
from .grid_solver import solve_grid

DIAGONALS = ("NE", "NW", "SE", "SW")
CARDINALS = ("N", "E", "S", "W")


@dataclass(frozen=True)
class DiagonalCase:
    name: str
    description: str
    diagonal_lw_mm: dict
    px: complex
    n_unknowns: int


def _solve(center_lw_mm, cardinal_lw_mm, diagonal_lw_mm, period_m, substrate, metal,
           frequency_hz, mesh_cells_per_subcell):
    grid_mm = eg.build_9cell_environment(center_lw_mm, cardinal_mm=cardinal_lw_mm, diagonal_mm=diagonal_lw_mm)
    grid_m = eg.grid_to_m(grid_mm)
    sg = solve_grid(3, period_m, substrate, metal, grid_m, frequency_hz, mesh_cells_per_subcell)
    return sg.dipole_at(1, 1).px, sg.n_unknowns


def run_diagonal_single_tests(center_lw_mm, contrast_lw_mm, period_m, substrate, metal,
                               frequency_hz, mesh_cells_per_subcell=8):
    """Diagonal Test 1 (Sec. 13): only one diagonal position differs from the center, repeated
    independently for NE, NW, SE, SW."""
    cases = {}
    for d in DIAGONALS:
        px, n_unk = _solve(center_lw_mm, None, {d: contrast_lw_mm}, period_m, substrate, metal,
                            frequency_hz, mesh_cells_per_subcell)
        cases[f"D1_{d}"] = DiagonalCase(f"D1_{d}", f"only {d} != C", {d: contrast_lw_mm}, px, n_unk)
    return cases


def diagonal_contributions(cases, s0):
    """Analogous to `neighbor_sensitivity.neighbor_contributions`, for the diagonal cases."""
    from .neighbor_sensitivity import NeighborContribution
    return {
        d: NeighborContribution(d, cases[f"D1_{d}"].px - s0, abs(cases[f"D1_{d}"].px - s0))
        for d in DIAGONALS
    }


@dataclass(frozen=True)
class CardinalVsDiagonal:
    contrast_label: str
    mean_abs_cardinal_delta: float
    mean_abs_diagonal_delta: float
    ratio_diagonal_to_cardinal: float


def compare_cardinal_vs_diagonal(contrast_label, cardinal_contributions, diagonal_contribs):
    """Sec. 14: mean |DeltaS| across the four cardinal single-neighbor tests vs. the four
    diagonal single-neighbor tests, at one contrast level. Does NOT assume rotational symmetry
    (each of the 4+4 values is measured independently; only the mean is compared here, the full
    per-direction values are reported separately, Sec. 10's "do not assume rotational symmetry")."""
    cardinal_mean = sum(c.delta_s_abs for c in cardinal_contributions.values()) / len(cardinal_contributions)
    diagonal_mean = sum(c.delta_s_abs for c in diagonal_contribs.values()) / len(diagonal_contribs)
    ratio = diagonal_mean / cardinal_mean if cardinal_mean > 0 else float("inf")
    return CardinalVsDiagonal(contrast_label, cardinal_mean, diagonal_mean, ratio)


@dataclass(frozen=True)
class CorrelatedEnvironmentResult:
    pattern_1_px: complex   # diagonals=A, cardinals=B  (A B A / B C B / A B A)
    pattern_2_px: complex   # diagonals=B, cardinals=A  (B A B / A C A / B A B)
    difference: complex
    difference_abs: float
    a_lw_mm: tuple
    b_lw_mm: tuple


def run_correlated_environment_test(center_lw_mm, a_lw_mm, b_lw_mm, period_m, substrate, metal,
                                     frequency_hz, mesh_cells_per_subcell=8):
    """Sec. 15: pattern 1 = A-B-A / B-C-B / A-B-A (diagonals=A, cardinals=B); pattern 2 =
    B-A-B / A-C-A / B-A-B (diagonals=B, cardinals=A). If the center response depended only on
    which geometry sits on the cardinals (i.e. diagonals were negligible), pattern 1 (cardinals=B)
    and pattern 2 (cardinals=A) would differ mainly because B != A on the cardinals, and swapping
    which set of positions gets A vs B would be symmetric information; if diagonals also matter,
    the two patterns carry genuinely different information and the difference reflects both
    swaps at once."""
    diag1 = {d: a_lw_mm for d in DIAGONALS}
    card1 = {d: b_lw_mm for d in CARDINALS}
    px1, _ = _solve(center_lw_mm, card1, diag1, period_m, substrate, metal, frequency_hz, mesh_cells_per_subcell)

    diag2 = {d: b_lw_mm for d in DIAGONALS}
    card2 = {d: a_lw_mm for d in CARDINALS}
    px2, _ = _solve(center_lw_mm, card2, diag2, period_m, substrate, metal, frequency_hz, mesh_cells_per_subcell)

    diff = px1 - px2
    return CorrelatedEnvironmentResult(px1, px2, diff, abs(diff), a_lw_mm, b_lw_mm)


def run_full_9cell_uniform(center_lw_mm, contrast_lw_mm, period_m, substrate, metal,
                            frequency_hz, mesh_cells_per_subcell=8):
    """All 8 neighbors (4 cardinal + 4 diagonal) set to the same perturbed geometry -- the
    genuine "9-cell" reference response used by `neighborhood_decision.py` to measure what a
    cardinal-only (5-cell) model misses by ignoring the diagonals."""
    uniform = {d: contrast_lw_mm for d in CARDINALS}
    diag_uniform = {d: contrast_lw_mm for d in DIAGONALS}
    px, n_unk = _solve(center_lw_mm, uniform, diag_uniform, period_m, substrate, metal,
                        frequency_hz, mesh_cells_per_subcell)
    return px, n_unk
