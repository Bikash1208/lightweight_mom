"""Phase 12A: one-port (grounded) ``Record``/``Database`` support."""

import numpy as np
import pytest

from src.database import CellEntry, Database, DatabaseError, Record

DIAG = {
    "n_unknowns": 10,
    "incoming_residual": 1e-15,
    "reciprocity_error": 1e-16,
    "passivity_min_eig": 0.0,
    "cells_per_wavelength": 20.0,
    "specular_only": True,
    "checks_passed": True,
}
S2 = np.array([[0.9 - 0.1j, 0.001j], [0.001j, 0.9 - 0.1j]])
S4 = np.eye(4, dtype=complex) * 0.5


def test_default_port_topology_is_two_port():
    r = Record("c", 0, 10e9, S4, DIAG)
    assert r.port_topology == "two_port"


def test_one_port_record_accepts_2x2_matrix():
    r = Record("c", 0, 12.5e9, S2, DIAG, port_topology="one_port")
    assert r.s.shape == (2, 2)


def test_one_port_record_rejects_4x4_matrix():
    with pytest.raises(DatabaseError):
        Record("c", 0, 12.5e9, S4, DIAG, port_topology="one_port")


def test_two_port_record_rejects_2x2_matrix():
    with pytest.raises(DatabaseError):
        Record("c", 0, 12.5e9, S2, DIAG, port_topology="two_port")


def test_invalid_port_topology_rejected():
    with pytest.raises(DatabaseError):
        Record("c", 0, 12.5e9, S2, DIAG, port_topology="three_port")


def test_component_access_undefined_family_raises():
    r = Record("c", 0, 12.5e9, S2, DIAG, port_topology="one_port")
    with pytest.raises(ValueError, match="one-port"):
        r.component("S21_xx")
    assert r.component("S11_xx") == pytest.approx(S2[0, 0])


def test_json_round_trip_preserves_one_port_topology():
    r = Record("cellA", 0, 12.5e9, S2, DIAG, port_topology="one_port")
    cell = CellEntry(
        0, "cellA", dict.fromkeys(("length_m", "width_m", "center_x_m", "center_y_m", "rotation_deg"), 0.0), "ok"
    )
    db = Database({}, [cell], [r])
    doc = db.to_dict()
    assert doc["cells"][0]["records"][0]["port_topology"] == "one_port"
    db2 = Database.from_dict(doc)
    rec2 = db2.records[0]
    assert rec2.port_topology == "one_port"
    assert np.allclose(rec2.s, S2)


def test_old_style_json_without_port_topology_key_defaults_to_two_port():
    r = Record("cellB", 0, 10e9, S4, DIAG)
    cell = CellEntry(
        0, "cellB", dict.fromkeys(("length_m", "width_m", "center_x_m", "center_y_m", "rotation_deg"), 0.0), "ok"
    )
    db = Database({}, [cell], [r])
    doc = db.to_dict()
    del doc["cells"][0]["records"][0]["port_topology"]  # simulate a pre-Phase-12A file
    db2 = Database.from_dict(doc)
    assert db2.records[0].port_topology == "two_port"
    assert np.allclose(db2.records[0].s, S4)


def test_to_csv_rejects_a_database_with_one_port_records(tmp_path):
    r = Record("cellA", 0, 12.5e9, S2, DIAG, port_topology="one_port")
    cell = CellEntry(
        0, "cellA", dict.fromkeys(("length_m", "width_m", "center_x_m", "center_y_m", "rotation_deg"), 0.0), "ok"
    )
    db = Database({}, [cell], [r])
    with pytest.raises(DatabaseError, match="CSV"):
        db.to_csv(tmp_path / "db.csv")


def test_same_as_true_for_identical_two_port_records():
    r1 = Record("c", 0, 10e9, S4, DIAG, port_topology="two_port")
    r2 = Record("c", 0, 10e9, S4, DIAG, port_topology="two_port")
    assert r1.same_as(r2)


def test_same_as_false_for_a_different_component_value():
    r1 = Record("c", 0, 12.5e9, S2, DIAG, port_topology="one_port")
    r2 = Record("c", 0, 12.5e9, S2 * 1.0000001, DIAG, port_topology="one_port")
    assert not r1.same_as(r2)
