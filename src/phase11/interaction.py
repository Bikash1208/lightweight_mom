"""Sec. 12-13 (Phase 11H/11I): pairwise and higher-order interaction diagnostics -- a physics
measurement, never a "surrogate-model feature-engineering exercise" (Sec. 12's own explicit
distinction). A perturbation is one (position, geometry) pair on an otherwise-uniform grid;
``response_for_subset`` solves the grid with exactly the requested SUBSET of perturbations applied
(every other position stays at the center's own geometry) and reads out the center response.
"""

from dataclasses import dataclass
from itertools import combinations

from .patterns import uniform
from .reference_solver import solve_grid


@dataclass(frozen=True)
class Perturbation:
    name: str
    i: int
    j: int
    lw_mm: tuple


def response_for_subset(center_lw_mm, perturbations, subset_names, grid_n, period_m, substrate,
                         metal, frequency_hz, mesh_cells_per_subcell=8):
    grid_mm = uniform(grid_n, center_lw_mm)
    by_name = {p.name: p for p in perturbations}
    for name in subset_names:
        p = by_name[name]
        grid_mm[p.i][p.j] = p.lw_mm
    grid_m = [[(l * 1e-3, w * 1e-3) for (l, w) in row] for row in grid_mm]
    sg = solve_grid(grid_n, period_m, substrate, metal, grid_m, frequency_hz, mesh_cells_per_subcell)
    mid = grid_n // 2
    return sg.response_at(mid, mid)


@dataclass(frozen=True)
class PairwiseInteractionResult:
    pair_label: str
    delta_a: complex
    delta_b: complex
    delta_ab_measured: complex
    delta_ab_additive_prediction: complex
    interaction_error: float
    interaction_error_relative: float


def run_pairwise_interaction(center_lw_mm, perturbations, pair_names, grid_n, period_m, substrate,
                              metal, frequency_hz, mesh_cells_per_subcell=8):
    """Sec. 12: for each (name_a, name_b) pair in ``pair_names``, computes DeltaR(A), DeltaR(B),
    DeltaR(A+B) and the additive-model interaction_error."""

    def resp(subset):
        return response_for_subset(center_lw_mm, perturbations, subset, grid_n, period_m, substrate,
                                    metal, frequency_hz, mesh_cells_per_subcell)

    r0 = resp(())
    single_cache = {}

    def single(name):
        if name not in single_cache:
            single_cache[name] = resp((name,)) - r0
        return single_cache[name]

    results = []
    for name_a, name_b in pair_names:
        delta_a, delta_b = single(name_a), single(name_b)
        measured = resp((name_a, name_b)) - r0
        predicted = delta_a + delta_b
        err = abs(measured - predicted)
        rel = err / abs(r0) if abs(r0) > 0 else float("inf")
        results.append(PairwiseInteractionResult(f"{name_a}+{name_b}", delta_a, delta_b, measured,
                                                   predicted, err, rel))
    return r0, results


@dataclass(frozen=True)
class HigherOrderInteractionResult:
    combo_label: str
    measured: complex
    pairwise_prediction: complex
    residual_beyond_pairwise: float
    residual_relative: float


def run_higher_order_interaction(center_lw_mm, perturbations, triple_names, grid_n, period_m,
                                  substrate, metal, frequency_hz, mesh_cells_per_subcell=8):
    """Sec. 13: for a triple (A, B, C), tests whether DeltaR(A+B+C) is explained by the THREE
    pairwise interactions (inclusion-exclusion on top of the three singles) or leaves a genuine
    residual only a three-body term can explain."""

    def resp(subset):
        return response_for_subset(center_lw_mm, perturbations, subset, grid_n, period_m, substrate,
                                    metal, frequency_hz, mesh_cells_per_subcell)

    names = list(triple_names)
    if len(names) != 3:
        raise ValueError(f"run_higher_order_interaction expects exactly 3 names, got {names!r}")
    r0 = resp(())
    singles = {n: resp((n,)) - r0 for n in names}
    pairs = {}
    for a, b in combinations(names, 2):
        pairs[(a, b)] = resp((a, b)) - r0

    measured_triple = resp(tuple(names)) - r0
    # Sec. 13's "lower-order interactions": the standard inclusion-exclusion prediction for a
    # response with pairwise (but no genuine three-body) interactions --
    #   V(ABC) - V0 = sum(pairwise deltas) - sum(single deltas)
    # (each single is counted twice across the three pairs, so it must be subtracted once).
    prediction = sum(pairs.values(), 0j) - sum(singles.values(), 0j)
    residual = abs(measured_triple - prediction)
    rel = residual / abs(r0) if abs(r0) > 0 else float("inf")
    label = "+".join(names)
    return HigherOrderInteractionResult(label, measured_triple, prediction, residual, rel)
