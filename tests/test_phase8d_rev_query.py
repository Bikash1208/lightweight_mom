"""Sec. 14, 20-21: the query engine. UDB-02 (period is a real query parameter, not global) and
UDB-03 (a nonexistent period never returns fabricated data) live here."""

import pytest

from src.phase8c import data_source
from src.phase8d_rev import library_entry as le
from src.phase8d_rev.query import UnitCellQuery, search_unit_cells
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


def test_exact_query_matches_the_registered_benchmark(udb):
    result = search_unit_cells(
        udb,
        UnitCellQuery(
            topology="rectangular_patch",
            substrate_name="Rogers RT/Duroid 5880",
            p_x_mm=4.8,
            p_y_mm=4.8,
            frequency_hz=12.5e9,
            polarization="x",
        ),
    )
    assert result.status == "MATCH"
    assert len(result.matching_entries) == 1


def test_period_mismatch_returns_database_coverage_missing_not_fabricated_data(udb):
    """UDB-02/UDB-03: changing only the query period changes the outcome, and there is never a
    substituted/rescaled match."""
    result = search_unit_cells(udb, UnitCellQuery(topology="rectangular_patch", p_x_mm=3.2, p_y_mm=3.2))
    assert result.status == "DATABASE_COVERAGE_MISSING"
    assert result.matching_entries == ()
    assert result.eliminated_at == "period"


def test_topology_mismatch_returns_database_coverage_missing(udb):
    result = search_unit_cells(udb, UnitCellQuery(topology="slot_loaded_patch"))
    assert result.status == "DATABASE_COVERAGE_MISSING"
    assert result.eliminated_at == "topology"


def test_unconstrained_query_matches_everything_registered(udb):
    result = search_unit_cells(udb, UnitCellQuery())
    assert result.status == "MATCH"
    assert len(result.matching_entries) == 1


def test_frequency_must_be_an_exact_stored_grid_point(udb):
    result = search_unit_cells(udb, UnitCellQuery(frequency_hz=12.50001e9))
    assert result.status == "DATABASE_COVERAGE_MISSING"
    assert result.eliminated_at == "frequency"


def test_filter_trail_records_every_stage_in_order(udb):
    result = search_unit_cells(udb, UnitCellQuery(topology="rectangular_patch"))
    names = [name for name, _ in result.filter_trail]
    assert names == [
        "initial",
        "topology",
        "material_stack",
        "frequency",
        "period",
        "polarization",
        "incidence",
        "response_availability",
        "validity_status",
    ]


def test_min_validation_status_filter(udb):
    # the benchmark entry's real validation status; requesting a status two levels above it
    # (HFSS_VALIDATED) must correctly eliminate it, since no HFSS validation exists anywhere.
    result = search_unit_cells(udb, UnitCellQuery(min_validation_status="HFSS_VALIDATED"))
    assert result.status == "DATABASE_COVERAGE_MISSING"
    assert result.eliminated_at == "validity_status"


def test_multiple_topologies_can_coexist_in_the_registry(phase8c_full_database):
    """UDB-07: topology is part of identity, so two entries under different topologies both
    survive registration and are independently queryable."""
    entry_a = le.from_database(phase8c_full_database, "lib_a", "rectangular_patch", le.REFERENCE_VALIDATION_LIBRARY)
    entry_b = le.from_database(phase8c_full_database, "lib_b", "slot_loaded_patch", le.REFERENCE_VALIDATION_LIBRARY)
    u = UniversalDatabase()
    u.register(entry_a)
    u.register(entry_b)
    assert set(u.topologies()) == {"rectangular_patch", "slot_loaded_patch"}
    result_a = search_unit_cells(u, UnitCellQuery(topology="rectangular_patch"))
    result_b = search_unit_cells(u, UnitCellQuery(topology="slot_loaded_patch"))
    assert result_a.matching_entries == (entry_a,)
    assert result_b.matching_entries == (entry_b,)
