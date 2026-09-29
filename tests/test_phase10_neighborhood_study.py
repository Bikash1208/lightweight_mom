"""Sec. 5-17, 47: the full per-center orchestration and the discretization-floor flag."""

from src.materials import Metal, Substrate
from src.phase10 import neighborhood_study as study_mod

SUB = Substrate("Rogers RT/Duroid 5880", 2.2, 0.0009, 1.575e-3)
MET = Metal()
P = 4.8e-3
F = 12.5e9


def test_run_center_study_covers_all_four_contrast_levels():
    r = study_mod.run_center_study("C3", (2.8, 2.8), P, SUB, MET, F, mesh_cells_per_subcell=8)
    assert set(r.contrast_values_mm.keys()) == {"weak", "moderate", "strong_low", "strong_high"}
    assert set(r.five_cell_cases.keys()) == {"weak", "moderate", "strong_low", "strong_high"}
    assert r.decision.final_decision in ("5_CELL_SUFFICIENT", "9_CELL_REQUIRED", "LARGER_NEIGHBORHOOD_REQUIRED")


def test_run_full_neighborhood_study_covers_every_requested_center():
    centers = {"C1": (1.4, 1.4), "C3": (2.8, 2.8)}
    results = study_mod.run_full_neighborhood_study(centers, P, SUB, MET, F, mesh_cells_per_subcell=8)
    assert set(results.keys()) == {"C1", "C3"}
    for r in results.values():
        assert r.decision.final_decision in ("5_CELL_SUFFICIENT", "9_CELL_REQUIRED", "LARGER_NEIGHBORHOOD_REQUIRED")


def test_discretization_floor_flag_matches_direct_finer_mesh_check():
    """C1's weak-contrast (0.2mm) perturbation is known (measured directly, see
    PHASE10_ENVIRONMENT_MODEL.md) to fall below this project's mesh resolution at
    mesh_cells_per_subcell=8; the flag must catch this rather than silently reporting a fake zero."""
    r = study_mod.run_center_study("C1", (1.4, 1.4), P, SUB, MET, F, mesh_cells_per_subcell=8)
    assert r.discretization_floor_flags["weak"] is True
