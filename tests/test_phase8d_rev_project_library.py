"""Sec. 17: project-library extraction from the universal database."""

import pytest

from src.phase8c import data_source
from src.phase8d_rev import library_entry as le
from src.phase8d_rev.project_library import extract_project_library
from src.phase8d_rev.query import UnitCellQuery
from src.phase8d_rev.universal_database import UniversalDatabase

pytestmark = pytest.mark.skipif(
    not data_source.DEFAULT_CHECKPOINT_JSON.exists(), reason="Phase 7 checkpoint not present"
)


@pytest.fixture
def udb(phase8c_full_database):
    entry = le.from_database(phase8c_full_database, "lib1", "rectangular_patch", le.REFERENCE_VALIDATION_LIBRARY)
    u = UniversalDatabase()
    u.register(entry)
    return u


def test_extract_project_library_on_a_matching_query(udb):
    project = extract_project_library(udb, UnitCellQuery(topology="rectangular_patch"))
    assert project.source_library_ids == ("lib1",)
    assert len(project.entries) == 1


def test_extract_project_library_raises_on_a_missing_query(udb):
    with pytest.raises(ValueError):
        extract_project_library(udb, UnitCellQuery(topology="rectangular_patch", p_x_mm=3.2, p_y_mm=3.2))
