"""Sec. 10-11: the 5-cell (cardinal-only) test matrix -- Tests 5A-5E -- and the neighbor
contribution / (non)additivity analysis built from it. Every configuration is solved through
`grid_solver.solve_grid` on a 3x3 grid whose diagonals are held at the center geometry
(Sec. 12's 5-cell representation; see `environment_geometry.build_5cell_environment`)."""

from dataclasses import dataclass

from . import environment_geometry as eg
from .grid_solver import solve_grid

CARDINALS = ("N", "E", "S", "W")


@dataclass(frozen=True)
class FiveCellCase:
    name: str
    description: str
    cardinal_lw_mm: dict     # only the directions that differ from the center
    px: complex
    n_unknowns: int


def _solve_case(name, description, center_lw_mm, cardinal_lw_mm, period_m, substrate, metal,
                 frequency_hz, mesh_cells_per_subcell):
    grid_mm = eg.build_5cell_environment(center_lw_mm, cardinal_lw_mm)
    grid_m = eg.grid_to_m(grid_mm)
    sg = solve_grid(3, period_m, substrate, metal, grid_m, frequency_hz, mesh_cells_per_subcell)
    d = sg.dipole_at(1, 1)
    return FiveCellCase(name, description, dict(cardinal_lw_mm), d.px, sg.n_unknowns)


def run_five_cell_test_matrix(center_lw_mm, contrast_lw_mm, period_m, substrate, metal,
                               frequency_hz, mesh_cells_per_subcell=8, asymmetric_lw_mm=None):
    """``contrast_lw_mm``: the single perturbed (L, W) applied to whichever direction(s) vary at
    this contrast level (Sec. 7-9's weak/moderate/strong value for this center, already resolved
    by ``environment_geometry.weak_neighbor``/``moderate_neighbor``/``strong_neighbor_low_high``).
    ``asymmetric_lw_mm``: optional explicit {"N":..,"E":..,"S":..,"W":..} for Test 5E; defaults to
    N=small/S=large (the domain extremes) with E=W=``contrast_lw_mm`` if omitted -- either way the
    exact geometries used are recorded in the returned cases (Sec. 9's "the exact geometries must
    be recorded")."""
    cases = {}

    def solve(name, desc, card):
        return _solve_case(name, desc, center_lw_mm, card, period_m, substrate, metal,
                            frequency_hz, mesh_cells_per_subcell)

    # Test 5A: homogeneous baseline.
    cases["5A_homogeneous"] = solve("5A_homogeneous", "N=E=S=W=C", {})

    # Test 5B: one-neighbor perturbation, repeated independently for N, E, S, W.
    for d in CARDINALS:
        cases[f"5B_{d}"] = solve(f"5B_{d}", f"only {d} != C", {d: contrast_lw_mm})

    # Test 5C: opposite-neighbor perturbation.
    cases["5C_NS"] = solve("5C_NS", "N=S != C, E=W=C", {"N": contrast_lw_mm, "S": contrast_lw_mm})
    cases["5C_EW"] = solve("5C_EW", "E=W != C, N=S=C", {"E": contrast_lw_mm, "W": contrast_lw_mm})

    # Test 5D: four-neighbor (uniform) perturbation.
    cases["5D_all"] = solve("5D_all", "N=E=S=W != C", {d: contrast_lw_mm for d in CARDINALS})

    # Test 5E: asymmetric environment.
    if asymmetric_lw_mm is None:
        small_lw, large_lw = eg.strong_neighbor_low_high(center_lw_mm)
        asymmetric_lw_mm = {
            "N": (small_lw, small_lw),
            "E": contrast_lw_mm,
            "S": (large_lw, large_lw),
            "W": contrast_lw_mm,
        }
    cases["5E_asymmetric"] = solve("5E_asymmetric", "N=small, E=moderate, S=large, W=moderate",
                                    asymmetric_lw_mm)
    return cases


@dataclass(frozen=True)
class NeighborContribution:
    direction: str
    delta_s: complex
    delta_s_abs: float


def neighbor_contributions(cases):
    """Sec. 10: DeltaS_i = S(G_C, G_i) - S(G_C, G_C), one per cardinal direction."""
    s0 = cases["5A_homogeneous"].px
    return {
        d: NeighborContribution(d, cases[f"5B_{d}"].px - s0, abs(cases[f"5B_{d}"].px - s0))
        for d in CARDINALS
    }


@dataclass(frozen=True)
class NonlinearityResult:
    combo: str
    s_measured: complex
    s_additive_prediction: complex
    e_nonlinear: float
    relative_to_s0: float


def nonlinearity_tests(cases):
    """Sec. 11: is S(N, E) approx S0 + DeltaS_N + DeltaS_E? Tested for both opposite pairs
    (Test 5C) and the full four-neighbor case (Test 5D) against the sum of all four single
    contributions."""
    s0 = cases["5A_homogeneous"].px
    contributions = neighbor_contributions(cases)
    results = []
    for combo, case_name in (("N+S", "5C_NS"), ("E+W", "5C_EW")):
        d1, d2 = combo.split("+")
        predicted = s0 + contributions[d1].delta_s + contributions[d2].delta_s
        measured = cases[case_name].px
        e_nl = abs(measured - predicted)
        results.append(NonlinearityResult(combo, measured, predicted, e_nl,
                                           e_nl / abs(s0) if abs(s0) > 0 else float("inf")))

    predicted_all = s0 + sum((contributions[d].delta_s for d in CARDINALS), 0j)
    measured_all = cases["5D_all"].px
    e_nl_all = abs(measured_all - predicted_all)
    results.append(NonlinearityResult("N+E+S+W", measured_all, predicted_all, e_nl_all,
                                       e_nl_all / abs(s0) if abs(s0) > 0 else float("inf")))
    return results
