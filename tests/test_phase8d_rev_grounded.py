"""Phase 12A: the stack/backing tag on ``LibraryEntry`` and its use by ``search_unit_cells`` so a
grounded and an ungrounded library are never mixed (Sec. 4's explicit requirement).
"""

import numpy as np
import pytest

from src.database import CellEntry, Database, Record
from src.phase8d_rev import library_entry as le
from src.phase8d_rev.query import UnitCellQuery, search_unit_cells
from src.phase8d_rev.universal_database import UniversalDatabase

DIAG = {
    "n_unknowns": 10, "incoming_residual": 1e-15, "reciprocity_error": 1e-16,
    "passivity_min_eig": 0.0, "cells_per_wavelength": 20.0, "specular_only": True,
    "checks_passed": True,
}
PARAMS = dict.fromkeys(("length_m", "width_m", "center_x_m", "center_y_m", "rotation_deg"), 0.0)
CONTEXT_COMMON = {
    "lattice": {"period_x_m": 4.8e-3, "period_y_m": 4.8e-3},
    "substrate": {"name": "RT5880", "eps_r": 2.2, "tan_delta": 0.0009, "thickness_m": 1.575e-3},
}


def _database(backing):
    ctx = dict(CONTEXT_COMMON, backing=backing)
    s = (
        np.array([[0.9 - 0.1j, 0.0], [0.0, 0.9 - 0.1j]])
        if backing == "pec"
        else np.eye(4, dtype=complex) * 0.5
    )
    r = Record("c0", 0, 12.5e9, s, DIAG, port_topology="one_port" if backing == "pec" else "two_port")
    cell = CellEntry(0, "c0", PARAMS, "ok")
    return Database({"context": ctx}, [cell], [r])


def _entry(backing, library_id):
    return le.from_database(_database(backing), library_id, "rectangular_patch", le.REFERENCE_VALIDATION_LIBRARY)


def test_backing_defaults_to_air_for_a_context_without_the_key():
    ctx_no_backing = dict(CONTEXT_COMMON)
    s = np.eye(4, dtype=complex) * 0.5
    r = Record("c0", 0, 12.5e9, s, DIAG)
    db = Database({"context": ctx_no_backing}, [CellEntry(0, "c0", PARAMS, "ok")], [r])
    entry = le.from_database(db, "lib", "rectangular_patch", le.REFERENCE_VALIDATION_LIBRARY)
    assert entry.backing == "air"


def test_backing_read_from_context_for_a_grounded_database():
    entry = _entry("pec", "lib_grounded")
    assert entry.backing == "pec"


def test_invalid_backing_rejected():
    with pytest.raises(ValueError):
        le.LibraryEntry(
            library_id="x", topology="t", role=le.PROJECT_LIBRARY, database=_database("air"),
            p_x_mm=4.8, p_y_mm=4.8, substrate={}, frequency_min_hz=1.0, frequency_max_hz=2.0,
            n_frequency_points=1, polarizations_supported=("x", "y"), incidence_supported=("normal",),
            solver_version=None, configuration_hash=None, source_type="MoM",
            validation_status="DATA_PRESENT", backing="copper",
        )


def test_search_with_explicit_backing_never_mixes_grounded_and_ungrounded():
    udb = UniversalDatabase()
    udb.register(_entry("air", "lib_air"))
    udb.register(_entry("pec", "lib_pec"))

    air_only = search_unit_cells(udb, UnitCellQuery(topology="rectangular_patch", backing="air"))
    assert air_only.status == "MATCH"
    assert [e.library_id for e in air_only.matching_entries] == ["lib_air"]

    pec_only = search_unit_cells(udb, UnitCellQuery(topology="rectangular_patch", backing="pec"))
    assert pec_only.status == "MATCH"
    assert [e.library_id for e in pec_only.matching_entries] == ["lib_pec"]

    unspecified = search_unit_cells(udb, UnitCellQuery(topology="rectangular_patch"))
    assert unspecified.status == "MATCH"
    assert len(unspecified.matching_entries) == 2


def test_search_for_a_backing_with_no_registered_library_reports_coverage_missing():
    udb = UniversalDatabase()
    udb.register(_entry("air", "lib_air"))
    result = search_unit_cells(udb, UnitCellQuery(topology="rectangular_patch", backing="pec"))
    assert result.status == "DATABASE_COVERAGE_MISSING"
    assert result.eliminated_at == "material_stack"
