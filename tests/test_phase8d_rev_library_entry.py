"""Sec. 3, 6, 17, 22: LibraryEntry wrapping the real Phase 7 checkpoint."""

import pytest

from src.phase8c import data_source
from src.phase8d_rev import library_entry as le

pytestmark = pytest.mark.skipif(
    not data_source.DEFAULT_CHECKPOINT_JSON.exists(), reason="Phase 7 checkpoint not present"
)


def test_benchmark_library_reports_the_real_4p8mm_period(phase8c_full_database):
    entry = le.from_database(phase8c_full_database, "lib1", "rectangular_patch", le.REFERENCE_VALIDATION_LIBRARY)
    assert entry.p_x_mm == pytest.approx(4.8)
    assert entry.p_y_mm == pytest.approx(4.8)


def test_benchmark_library_reports_324_geometries(phase8c_full_database):
    entry = le.from_database(phase8c_full_database, "lib1", "rectangular_patch", le.REFERENCE_VALIDATION_LIBRARY)
    assert entry.n_geometries == 324
    assert entry.n_ok_cells == 324


def test_from_database_rejects_unknown_role(phase8c_full_database):
    with pytest.raises(ValueError):
        le.LibraryEntry(
            library_id="x",
            topology="rectangular_patch",
            role="not_a_role",
            database=phase8c_full_database,
            p_x_mm=4.8,
            p_y_mm=4.8,
            substrate={},
            frequency_min_hz=1.0,
            frequency_max_hz=2.0,
            n_frequency_points=2,
            polarizations_supported=("x", "y"),
            incidence_supported=("normal",),
            solver_version=None,
            configuration_hash=None,
            source_type="MoM",
            validation_status="DATA_PRESENT",
        )


def test_matches_period_exact_and_mismatched(phase8c_full_database):
    entry = le.from_database(phase8c_full_database, "lib1", "rectangular_patch", le.REFERENCE_VALIDATION_LIBRARY)
    assert entry.matches_period(4.8, 4.8)
    assert not entry.matches_period(3.2, 3.2)


def test_has_frequency_exact_grid_point_only(phase8c_full_database):
    entry = le.from_database(phase8c_full_database, "lib1", "rectangular_patch", le.REFERENCE_VALIDATION_LIBRARY)
    assert entry.has_frequency(12.5e9)
    assert not entry.has_frequency(12.5001e9)


def test_record_raises_keyerror_for_a_geometry_not_in_the_library(phase8c_full_database):
    entry = le.from_database(phase8c_full_database, "lib1", "rectangular_patch", le.REFERENCE_VALIDATION_LIBRARY)
    with pytest.raises(KeyError):
        entry.record(999.0, 999.0, 12.5e9)


def test_all_ok_records_recovers_every_geometry_at_design_frequency(phase8c_full_database):
    entry = le.from_database(phase8c_full_database, "lib1", "rectangular_patch", le.REFERENCE_VALIDATION_LIBRARY)
    records = entry.all_ok_records(12.5e9)
    assert len(records) == 324


def test_benchmark_library_loader_cross_checks_and_registers_manifest_provenance():
    entry = le.benchmark_library()
    assert entry.library_id == le.BENCHMARK_LIBRARY_ID
    assert entry.role == le.REFERENCE_VALIDATION_LIBRARY
    # manifest.json exists for the real delivered library, so provenance should not be None here.
    assert entry.solver_version is not None
    assert entry.configuration_hash is not None
