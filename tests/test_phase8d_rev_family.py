"""Sec. 9: UnitCellFamily, aggregated from real LibraryEntry data (parameters computed, not
hand-listed)."""

import pytest

from src.phase8c import data_source
from src.phase8d_rev import library_entry as le
from src.phase8d_rev.family import family_from_entries

pytestmark = pytest.mark.skipif(
    not data_source.DEFAULT_CHECKPOINT_JSON.exists(), reason="Phase 7 checkpoint not present"
)


def test_family_from_entries_rejects_empty_list():
    with pytest.raises(ValueError):
        family_from_entries("rectangular_patch", [])


def test_family_from_entries_rejects_mixed_topologies(phase8c_full_database):
    entry = le.from_database(phase8c_full_database, "lib1", "rectangular_patch", le.REFERENCE_VALIDATION_LIBRARY)
    bad = le.LibraryEntry(**{**entry.__dict__, "topology": "other_topology"})
    with pytest.raises(ValueError):
        family_from_entries("rectangular_patch", [entry, bad])


def test_family_detects_length_and_width_as_the_only_varying_parameters(phase8c_full_database):
    entry = le.from_database(phase8c_full_database, "lib1", "rectangular_patch", le.REFERENCE_VALIDATION_LIBRARY)
    family = family_from_entries("rectangular_patch", [entry])
    assert set(family.parameter_schema) == {"length_m", "width_m"}
    assert set(family.fixed_parameters) == {"center_x_m", "center_y_m", "rotation_deg"}


def test_family_reports_ground_plane_absent_and_normal_incidence_only(phase8c_full_database):
    entry = le.from_database(phase8c_full_database, "lib1", "rectangular_patch", le.REFERENCE_VALIDATION_LIBRARY)
    family = family_from_entries("rectangular_patch", [entry])
    assert "ABSENT" in family.validity_constraints["ground_plane"]
    assert family.validity_constraints["normal_incidence_only"] is True


def test_family_supported_periods_reflects_the_one_registered_entry(phase8c_full_database):
    entry = le.from_database(phase8c_full_database, "lib1", "rectangular_patch", le.REFERENCE_VALIDATION_LIBRARY)
    family = family_from_entries("rectangular_patch", [entry])
    assert family.supported_periods_mm == ((4.8, 4.8),)
