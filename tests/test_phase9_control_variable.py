"""Phase 9A: reference-cell extraction and classification, against the real Phase 7 checkpoint."""

import pytest

from src.phase8c import data_source, power_balance
from src.phase9.control_variable import (
    REFERENCE_CELLS_MM,
    VALIDATION_FREQUENCIES_GHZ,
    architecture_audit,
    classify_control_variable,
    compare_against_hfss,
    extract_solver_reference_data,
)
from src.phase9.error_metrics import compute_error_stats
from src.phase9.hfss_reference import HfssRecord
from src.phase_extraction import design_frequency_index

pytestmark = pytest.mark.skipif(
    not data_source.DEFAULT_CHECKPOINT_JSON.exists(), reason="Phase 7 checkpoint not present"
)


def test_validation_frequencies_are_exact_grid_points(phase8c_full_database):
    freqs_hz = phase8c_full_database.frequencies
    for f_ghz in VALIDATION_FREQUENCIES_GHZ:
        f_hz = f_ghz * 1e9
        assert any(abs(f - f_hz) <= 1e-6 * f_hz for f in freqs_hz), f_ghz


def test_extract_solver_reference_data_covers_all_10_cells_and_7_frequencies(phase8c_full_database):
    rows = extract_solver_reference_data(phase8c_full_database)
    assert len(rows) == len(REFERENCE_CELLS_MM) * len(VALIDATION_FREQUENCIES_GHZ)
    assert {r.cell_id for r in rows} == {c[0] for c in REFERENCE_CELLS_MM}


def test_extract_solver_reference_data_raises_for_a_geometry_not_in_the_database(phase8c_full_database):
    with pytest.raises(KeyError):
        extract_solver_reference_data(phase8c_full_database, cells_mm=(("BAD", 999.0, 999.0),))


def test_architecture_audit_marks_every_hfss_field_unknown(phase8c_full_database):
    idx = design_frequency_index(phase8c_full_database, 12.5e9)
    mean_power = power_balance.library_power_summary(phase8c_full_database, idx)
    audit = architecture_audit(phase8c_full_database, mean_power)
    hfss_fields = [k for k in audit if k.startswith("hfss_")]
    assert hfss_fields
    assert all("UNKNOWN" in str(audit[k]) for k in hfss_fields)


def test_architecture_audit_documents_the_real_ground_plane_status(phase8c_full_database):
    idx = design_frequency_index(phase8c_full_database, 12.5e9)
    mean_power = power_balance.library_power_summary(phase8c_full_database, idx)
    audit = architecture_audit(phase8c_full_database, mean_power)
    assert audit["solver_ground_plane"] == "ABSENT"


def test_classify_control_variable_reports_not_available_when_stats_is_none():
    status, reason = classify_control_variable(None)
    assert status == "UNVALIDATED"
    assert "NOT_AVAILABLE_FROM_CURRENT_DATA" in reason


def test_classify_control_variable_detects_a_180_degree_sign_mismatch():
    stats = compute_error_stats([1.0 + 0j, 0.5 + 0j], [-1.0 + 0j, -0.5 + 0j])
    status, reason = classify_control_variable(stats)
    assert status == "INCONSISTENT"


def test_classify_control_variable_never_asserts_validated_without_a_sourced_tolerance():
    stats = compute_error_stats([1.0 + 0j], [1.001 + 0j])  # a tiny, "obviously good" agreement
    status, _ = classify_control_variable(stats)
    assert status == "UNVALIDATED"  # still not VALIDATED -- no tolerance is defined anywhere


def test_compare_against_hfss_raises_if_a_solver_row_has_no_matching_hfss_record(phase8c_full_database):
    rows = extract_solver_reference_data(phase8c_full_database)
    hfss_records = (HfssRecord("CELL-01", 1.0, 1.0, 12.5, 0.1 + 0.1j, 0.9 + 0.1j),)
    with pytest.raises(ValueError):
        compare_against_hfss(hfss_records, rows, "s11")
