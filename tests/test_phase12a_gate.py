"""Phase 12A: gate/status-code logic. Never a bare pass/fail; a negative result must survive
(brief's project-wide philosophy) -- this exercises both the VALIDATED and PARTIAL branches.
"""

from src.phase12a import gate


def test_all_checks_within_tolerance_gives_validated():
    result = gate.decide(
        max_analytic_error=1e-15, floquet_error=0.001, mesh_error=0.005,
        any_surface_wave_in_band=False, coverage_deg=200.0,
    )
    assert result.status == gate.GROUNDED_STACK_VALIDATED
    assert result.measured_phase_coverage_deg == 200.0


def test_analytic_error_over_tolerance_gives_partial():
    result = gate.decide(
        max_analytic_error=1e-3, floquet_error=0.001, mesh_error=0.005,
        any_surface_wave_in_band=False, coverage_deg=200.0,
    )
    assert result.status == gate.GROUNDED_STACK_PARTIAL
    assert "textbook" in result.reason


def test_surface_wave_in_band_gives_partial():
    result = gate.decide(
        max_analytic_error=1e-15, floquet_error=0.001, mesh_error=0.005,
        any_surface_wave_in_band=True, coverage_deg=200.0,
    )
    assert result.status == gate.GROUNDED_STACK_PARTIAL
    assert "surface-wave" in result.reason


def test_mesh_or_floquet_error_over_tolerance_gives_partial():
    result = gate.decide(
        max_analytic_error=1e-15, floquet_error=0.001, mesh_error=0.05,
        any_surface_wave_in_band=False, coverage_deg=200.0,
    )
    assert result.status == gate.GROUNDED_STACK_PARTIAL
    assert "mesh/Floquet" in result.reason


def test_geometry_inversion_and_step_export_always_locked():
    result = gate.decide(1e-15, 0.001, 0.005, False, 200.0)
    assert result.geometry_inversion_status == gate.GEOMETRY_INVERSION_LOCKED
    assert result.step_export_status == gate.STEP_EXPORT_LOCKED


def test_hfss_status_not_available_without_a_reference_file(tmp_path):
    assert gate.hfss_status(tmp_path / "does_not_exist.csv") == gate.HFSS_NOT_AVAILABLE
