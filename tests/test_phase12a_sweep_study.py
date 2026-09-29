"""Phase 12A: the physics-sanity sweep and its comparison with the real Phase 7 library."""

import pytest

from src.phase8c.data_source import DEFAULT_CHECKPOINT_JSON
from src.phase12a.sweep_study import (
    coverage_by_frequency,
    grounded_library_cell,
    run_grounded_patch_sweep,
    ungrounded_phase7_coverage,
    ungrounded_phase7_diagonal_table,
)

SMALL_L_MM = (1.0, 2.0, 3.0, 4.4)
DESIGN_F = 12.5e9


def test_sweep_produces_one_point_per_l_and_frequency():
    points = run_grounded_patch_sweep(l_values_mm=SMALL_L_MM, frequencies_hz=(DESIGN_F,))
    assert len(points) == len(SMALL_L_MM)
    assert all(p.frequency_hz == DESIGN_F for p in points)


def test_sweep_is_lossless_reflection_only_up_to_solver_precision():
    """The library substrate is lossy (tan_delta=0.0009), so |S11| should be < 1 but very close."""
    points = run_grounded_patch_sweep(l_values_mm=SMALL_L_MM, frequencies_hz=(DESIGN_F,))
    for p in points:
        assert 0.99 < p.s11_xx_mag <= 1.0
        assert p.s11_yx_mag < 1e-10  # unrotated square: no cross-pol


def test_coverage_by_frequency_is_a_real_circular_coverage():
    points = run_grounded_patch_sweep(l_values_mm=SMALL_L_MM, frequencies_hz=(DESIGN_F,))
    cov = coverage_by_frequency(points)[DESIGN_F]
    assert 0.0 < cov.coverage_deg <= 360.0
    assert cov.n_points == len(SMALL_L_MM)


def test_grounded_library_cell_matches_phase7_period_and_substrate():
    cell = grounded_library_cell(2.0e-3)
    assert cell.lattice.period_x == pytest.approx(4.8e-3)
    assert cell.substrate.eps_r == pytest.approx(2.2)
    assert cell.backing == "pec"


@pytest.mark.skipif(not DEFAULT_CHECKPOINT_JSON.exists(), reason="Phase 7 checkpoint not present")
def test_ungrounded_phase7_coverage_matches_the_documented_24_8_degrees():
    result = ungrounded_phase7_coverage(frequency_hz=DESIGN_F)
    assert result["status"] == "AVAILABLE"
    assert result["S11_xx"]["coverage_deg"] == pytest.approx(24.8, abs=0.1)
    assert result["S21_xx"]["coverage_deg"] == pytest.approx(43.9, abs=0.1)


@pytest.mark.skipif(not DEFAULT_CHECKPOINT_JSON.exists(), reason="Phase 7 checkpoint not present")
def test_ungrounded_phase7_diagonal_table_has_18_points():
    rows = ungrounded_phase7_diagonal_table(frequency_hz=DESIGN_F)
    assert len(rows) == 18
    assert all("L_mm" in r and "phase_deg" in r for r in rows)


def test_ungrounded_phase7_coverage_reports_not_available_for_a_missing_checkpoint(tmp_path):
    result = ungrounded_phase7_coverage(checkpoint_path=tmp_path / "does_not_exist.json")
    assert result == {"status": "NOT_AVAILABLE_FROM_CURRENT_DATA"}


def test_ungrounded_phase7_diagonal_table_empty_for_a_missing_checkpoint(tmp_path):
    rows = ungrounded_phase7_diagonal_table(checkpoint_path=tmp_path / "does_not_exist.json")
    assert rows == []
