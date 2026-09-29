"""Sec. 21-22, 28: deterministic, non-ML geometry inversion (isolated-cell ranking + direct-MoM
verification), reachable/unreachable targets."""

from src.materials import Metal, Substrate
from src.phase11 import inverse as inv_mod
from src.phase11 import patterns as pt

SUB = Substrate("Rogers RT/Duroid 5880", 2.2, 0.0009, 1.575e-3)
MET = Metal()
P = 4.8e-3
F = 12.5e9
CENTER = (2.8, 2.8)


def test_invert_target_direct_reaches_an_exactly_achievable_target():
    env = pt.uniform(3, CENTER)
    from src.phase11.reference_solver import solve_isolated

    target, _ = solve_isolated(P, CENTER[0] * 1e-3, CENTER[1] * 1e-3, SUB, MET, F, mesh_cells_per_side=8)
    result = inv_mod.invert_target_direct(
        target, env, 1, 1, P, SUB, MET, F, mesh_cells_per_subcell=8, candidate_lw_values_mm=[2.8]
    )
    assert result.status == inv_mod.FOUND
    assert result.best_candidate_lw_mm == (2.8, 2.8)
    assert result.mom_relative_error < 1e-6


def test_invert_target_direct_rejects_an_absurd_target():
    env = pt.uniform(3, CENTER)
    result = inv_mod.invert_target_direct(
        complex(1.0, 0.0), env, 1, 1, P, SUB, MET, F, mesh_cells_per_subcell=8, candidate_lw_values_mm=[1.0, 2.8, 4.4]
    )
    assert result.status == inv_mod.TARGET_UNREACHABLE


def test_grid_values_span_the_full_validated_domain():
    values = inv_mod._grid_values_mm()
    assert values[0] == 1.0 and values[-1] == 4.4
    assert len(values) == 18
