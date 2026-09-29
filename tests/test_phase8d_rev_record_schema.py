"""Sec. 7-8: UnitCellRecord, built from a real cell/record of the Phase 7 checkpoint."""

import pytest

from src.phase8c import data_source
from src.phase8d_rev.record_schema import LINEAR_COMPONENT_NAMES, record_from_database

pytestmark = pytest.mark.skipif(not data_source.DEFAULT_CHECKPOINT_JSON.exists(), reason="Phase 7 checkpoint not present")


def _first_ok_cell_and_record(database):
    cell = next(c for c in database.cells if c.status == "ok")
    record = database.record(cell.index, 0)
    return cell, record


def test_record_from_database_preserves_the_complex_value_as_primary(phase8c_full_database):
    cell, record = _first_ok_cell_and_record(phase8c_full_database)
    r = record_from_database(phase8c_full_database, "lib1", "rectangular_patch", cell, record)
    for name in LINEAR_COMPONENT_NAMES:
        entry = r.s_parameters[name]
        z = record.component(name)
        assert entry["re"] == pytest.approx(z.real)
        assert entry["im"] == pytest.approx(z.imag)
        # magnitude/phase are DERIVED from re/im, not stored independently of them
        assert entry["mag"] == pytest.approx(abs(z))


def test_record_from_database_preserves_cross_polarization_components(phase8c_full_database):
    cell, record = _first_ok_cell_and_record(phase8c_full_database)
    r = record_from_database(phase8c_full_database, "lib1", "rectangular_patch", cell, record)
    assert "S11_xy" in r.s_parameters
    assert "S11_yx" in r.s_parameters
    assert "S21_xy" in r.s_parameters
    assert "S21_yx" in r.s_parameters


def test_record_from_database_geometry_matches_cell_parameters(phase8c_full_database):
    cell, record = _first_ok_cell_and_record(phase8c_full_database)
    r = record_from_database(phase8c_full_database, "lib1", "rectangular_patch", cell, record)
    assert r.L_mm == pytest.approx(cell.parameters["length_m"] * 1e3)
    assert r.W_mm == pytest.approx(cell.parameters["width_m"] * 1e3)
    assert r.p_x_mm == pytest.approx(4.8)
    assert r.p_y_mm == pytest.approx(4.8)


def test_record_from_database_normal_incidence_only(phase8c_full_database):
    cell, record = _first_ok_cell_and_record(phase8c_full_database)
    r = record_from_database(phase8c_full_database, "lib1", "rectangular_patch", cell, record)
    assert r.incidence_theta_deg == 0.0
    assert r.incidence_phi_deg == 0.0


def test_record_from_database_phase_coverage_metadata_never_fabricated_per_record(phase8c_full_database):
    cell, record = _first_ok_cell_and_record(phase8c_full_database)
    r = record_from_database(phase8c_full_database, "lib1", "rectangular_patch", cell, record)
    assert r.phase_coverage_metadata is None


def test_record_from_database_validation_status_is_data_driven(phase8c_full_database):
    cell, record = _first_ok_cell_and_record(phase8c_full_database)
    r = record_from_database(phase8c_full_database, "lib1", "rectangular_patch", cell, record)
    expected = "NUMERICALLY_VALIDATED" if record.diagnostics["checks_passed"] else "DATA_PRESENT"
    assert r.provenance.validation_status == expected
