"""Sec. 9, 28: literal reuse of Phase 10's own diagonal-coupling functions to reproduce/re-establish
its 25-70% finding."""

from src.materials import Metal, Substrate
from src.phase11.diagonal_study_reuse import reproduce_diagonal_finding

SUB = Substrate("Rogers RT/Duroid 5880", 2.2, 0.0009, 1.575e-3)
MET = Metal()
P = 4.8e-3
F = 12.5e9
CENTER = (2.8, 2.8)


def test_reproduce_diagonal_finding_reports_a_positive_finite_ratio():
    result = reproduce_diagonal_finding(CENTER, P, SUB, MET, F, mesh_cells_per_subcell=8)
    assert result["cardinal_mean"] > 0
    assert result["diagonal_mean"] > 0
    assert result["ratio"] > 0
    assert isinstance(result["reproduces_phase10_range"], bool)
    assert set(result["cardinal_contributions"].keys()) == {"N", "E", "S", "W"}
    assert set(result["diagonal_contributions"].keys()) == {"NE", "NW", "SE", "SW"}
