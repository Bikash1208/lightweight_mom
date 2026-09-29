"""Sec. 24-25: UDB-01, UDB-04, UDB-05, UDB-11, UDB-12's integrity checks."""

import pytest

from src.phase8c import data_source
from src.phase8d_rev import library_entry as le
from src.phase8d_rev.integrity import run_integrity_checks
from src.phase8d_rev.universal_database import UniversalDatabase

pytestmark = pytest.mark.skipif(not data_source.DEFAULT_CHECKPOINT_JSON.exists(), reason="Phase 7 checkpoint not present")


def test_udb01_all_324_geometries_recovered_unchanged(phase8c_full_database):
    entry = le.from_database(phase8c_full_database, "lib1", "rectangular_patch", le.REFERENCE_VALIDATION_LIBRARY)
    assert entry.n_ok_cells == 324


def test_udb04_and_udb05_complex_and_cross_pol_survive_round_trip(phase8c_full_database, tmp_path):
    entry = le.from_database(phase8c_full_database, "lib1", "rectangular_patch", le.REFERENCE_VALIDATION_LIBRARY)
    udb = UniversalDatabase()
    udb.register(entry)
    results = run_integrity_checks(udb, tmp_path)
    by_id = {(lib, r.check_id): r for lib, r in results}
    assert by_id[("lib1", "complex_valued_response_integrity")].passed
    assert by_id[("lib1", "cross_polarization_preserved")].passed
    assert by_id[("lib1", "round_trip_serialization")].passed


def test_udb11_missing_value_handling_and_metadata_completeness(phase8c_full_database, tmp_path):
    entry = le.from_database(phase8c_full_database, "lib1", "rectangular_patch", le.REFERENCE_VALIDATION_LIBRARY)
    udb = UniversalDatabase()
    udb.register(entry)
    results = run_integrity_checks(udb, tmp_path)
    by_id = {(lib, r.check_id): r for lib, r in results}
    assert by_id[("lib1", "missing_value_handling")].passed
    assert by_id[("lib1", "metadata_completeness")].passed


def test_identity_consistency_fails_on_duplicate_library_ids_bypassing_the_registry(phase8c_full_database, tmp_path):
    """The registry itself refuses duplicates (test_phase8d_rev_universal_database.py); this
    exercises the standalone consistency check directly against a hand-built duplicate list."""
    from src.phase8d_rev.integrity import check_topology_material_period_polarization_consistency

    entry = le.from_database(phase8c_full_database, "lib1", "rectangular_patch", le.REFERENCE_VALIDATION_LIBRARY)
    result = check_topology_material_period_polarization_consistency([entry, entry])
    assert not result.passed


def test_run_integrity_checks_without_tmp_path_skips_round_trip_only(phase8c_full_database):
    entry = le.from_database(phase8c_full_database, "lib1", "rectangular_patch", le.REFERENCE_VALIDATION_LIBRARY)
    udb = UniversalDatabase()
    udb.register(entry)
    results = run_integrity_checks(udb, tmp_path=None)
    check_ids = {r.check_id for _, r in results}
    assert "round_trip_serialization" not in check_ids
    assert "unique_record_identity" in check_ids
