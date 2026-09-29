"""Sec. 6, 17, 20: UniversalDatabase registry behaviour."""

import pytest

from src.phase8c import data_source
from src.phase8d_rev import library_entry as le
from src.phase8d_rev.universal_database import UniversalDatabase, UniversalDatabaseError, with_benchmark

pytestmark = pytest.mark.skipif(
    not data_source.DEFAULT_CHECKPOINT_JSON.exists(), reason="Phase 7 checkpoint not present"
)


def test_register_rejects_duplicate_library_id(phase8c_full_database):
    entry = le.from_database(phase8c_full_database, "dup", "rectangular_patch", le.REFERENCE_VALIDATION_LIBRARY)
    udb = UniversalDatabase()
    udb.register(entry)
    with pytest.raises(UniversalDatabaseError):
        udb.register(entry)


def test_topologies_and_entries_for_topology(phase8c_full_database):
    entry = le.from_database(phase8c_full_database, "lib1", "rectangular_patch", le.REFERENCE_VALIDATION_LIBRARY)
    udb = UniversalDatabase()
    udb.register(entry)
    assert udb.topologies() == ("rectangular_patch",)
    assert udb.entries_for_topology("rectangular_patch") == (entry,)
    assert udb.entries_for_topology("nonexistent") == ()


def test_family_raises_keyerror_for_unregistered_topology():
    udb = UniversalDatabase()
    with pytest.raises(KeyError):
        udb.family("rectangular_patch")


def test_periods_for_topology_never_assumes_4p8mm(phase8c_full_database):
    entry = le.from_database(phase8c_full_database, "lib1", "rectangular_patch", le.REFERENCE_VALIDATION_LIBRARY)
    udb = UniversalDatabase()
    assert udb.periods_for_topology("rectangular_patch") == ()  # nothing registered -> nothing assumed
    udb.register(entry)
    assert udb.periods_for_topology("rectangular_patch") == ((4.8, 4.8),)


def test_with_benchmark_registers_exactly_one_reference_validation_library_entry():
    udb = with_benchmark()
    entries = udb.entries()
    assert len(entries) == 1
    assert entries[0].role == le.REFERENCE_VALIDATION_LIBRARY
    assert entries[0].n_ok_cells == 324
