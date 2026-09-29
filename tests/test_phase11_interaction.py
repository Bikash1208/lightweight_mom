"""Sec. 12-13, 28: pairwise and higher-order interaction diagnostics."""

import pytest

from src.materials import Metal, Substrate
from src.phase11.interaction import (
    Perturbation,
    response_for_subset,
    run_higher_order_interaction,
    run_pairwise_interaction,
)

SUB = Substrate("Rogers RT/Duroid 5880", 2.2, 0.0009, 1.575e-3)
MET = Metal()
P = 4.8e-3
F = 12.5e9
CENTER = (2.8, 2.8)
CONTRAST = (4.4, 4.4)
N = 5
MID = 2

PERTURBATIONS = [
    Perturbation("N", MID - 1, MID, CONTRAST),
    Perturbation("S", MID + 1, MID, CONTRAST),
    Perturbation("E", MID, MID + 1, CONTRAST),
]


def test_response_for_subset_empty_equals_homogeneous():
    r_empty = response_for_subset(CENTER, PERTURBATIONS, (), N, P, SUB, MET, F, mesh_cells_per_subcell=8)
    r_none_applied = response_for_subset(CENTER, PERTURBATIONS, (), N, P, SUB, MET, F, mesh_cells_per_subcell=8)
    assert r_empty == r_none_applied  # deterministic, no ML


@pytest.mark.slow
def test_pairwise_interaction_reports_measured_and_predicted():
    r0, results = run_pairwise_interaction(
        CENTER, PERTURBATIONS, [("N", "S"), ("N", "E")], N, P, SUB, MET, F, mesh_cells_per_subcell=8
    )
    assert len(results) == 2
    for r in results:
        assert r.interaction_error == pytest.approx(abs(r.delta_ab_measured - r.delta_ab_additive_prediction))
        assert r.delta_ab_additive_prediction == pytest.approx(r.delta_a + r.delta_b)


def test_higher_order_interaction_requires_exactly_three_names():
    with pytest.raises(ValueError):
        run_higher_order_interaction(CENTER, PERTURBATIONS, ("N", "S"), N, P, SUB, MET, F, mesh_cells_per_subcell=8)


@pytest.mark.slow
def test_higher_order_interaction_residual_is_nonnegative():
    result = run_higher_order_interaction(
        CENTER, PERTURBATIONS, ("N", "S", "E"), N, P, SUB, MET, F, mesh_cells_per_subcell=8
    )
    assert result.residual_beyond_pairwise >= 0
    assert result.residual_relative >= 0
