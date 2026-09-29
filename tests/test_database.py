"""Unit-cell database: JSON / CSV round trips, schema strictness, queries, metadata."""

import copy
import csv
import json

import numpy as np
import pytest

from src.database import (
    CONVENTIONS,
    CSV_COLUMNS,
    DIAGNOSTIC_TYPES,
    FORMAT_NAME,
    CellEntry,
    Database,
    DatabaseError,
    Record,
    component_names,
    component_value,
    family_block,
)
from src.phase import principal_phase_deg, wrapped_difference_deg
from src.polarization import block_to_circular
from src.runner import run_sweep
from src.sparameters import compute_s_parameters
from sweep_helpers import FREQ_8_12, sweep_config

SWEEP = {
    "frequency": {"start": "8 GHz", "stop": "10 GHz", "step": "1 GHz"},
    "parameters": {"length": {"values": ["3 mm", "9 mm", "4 mm"]}, "rotation_deg": [0, 30]},
    "on_invalid": "skip",
    "validation": {"refinement": {"enabled": True, "tolerance": 0.05}},
}


@pytest.fixture(scope="module")
def result():
    return run_sweep(sweep_config(SWEEP))


@pytest.fixture(scope="module")
def db(result):
    return result.database


@pytest.fixture()
def files(db, tmp_path):
    db.to_json(tmp_path / "db.json")
    db.to_csv(tmp_path / "db.csv")
    return tmp_path / "db.json", tmp_path / "db.csv"


# ------------------------------------------------------------------------------ content

def test_database_shape_and_skipped_cells(db):
    assert len(db.cells) == 6 and len(db.records) == 4 * 3
    assert [c.status for c in db.cells] == ["ok", "ok", "invalid", "invalid", "ok", "ok"]
    assert db.cells[2].reason and "extends outside" in db.cells[2].reason
    assert all(db.records_of(c.index) == [] for c in db.cells if c.status == "invalid")
    assert db.frequencies.tolist() == [8e9, 9e9, 10e9]


def test_metadata_is_complete_and_states_conventions(db):
    m = db.metadata
    assert m["conventions"] == CONVENTIONS and "exp(+j*omega*t)" in m["conventions"]["time_dependence"]
    assert m["context"]["lattice"] == {"period_x_m": 8e-3, "period_y_m": 8e-3}
    assert m["context"]["substrate"]["eps_r"] == 2.2 and m["context"]["mesh"]["floquet_order"] == 32
    assert m["frequency_grid"]["frequencies_hz"] == [8e9, 9e9, 10e9] and m["frequency_grid"]["count"] == 3
    assert [a["name"] for a in m["parameter_axes"]] == ["length_m", "rotation_deg"]
    assert m["validation_policy"]["refinement"]["enabled"] is True and m["on_invalid"] == "skip"
    assert m["context"]["reference_planes"] == {"top_air_thickness_m": 0.0, "bottom_air_thickness_m": 0.0}
    assert "S21_LR" in json.dumps(m["conventions"])


def test_conventions_document_the_circular_basis_and_handedness_explicitly():
    c = CONVENTIONS["circular_basis"]
    assert c["R"] == "(x + j y)/sqrt(2)" and c["L"] == "(x - j y)/sqrt(2)"
    assert "clockwise seen from +z" in c["rotation_in_time"]
    assert "RHCP when travelling -z" in c["handedness"]


def test_record_diagnostics_are_typed_and_physically_sensible(db):
    for r in db.records:
        assert set(r.diagnostics) == set(DIAGNOSTIC_TYPES)
        assert r.diagnostics["incoming_residual"] < 1e-12
        assert r.diagnostics["reciprocity_error"] < 1e-10
        assert r.diagnostics["passivity_min_eig"] > -1e-10            # lossy cell: strictly passive
        assert r.diagnostics["specular_only"] is True and r.diagnostics["checks_passed"] is True
        assert isinstance(r.diagnostics["n_unknowns"], int)


def test_component_accessors_agree_with_the_matrix(db):
    r = db.records[0]
    assert r.component("S21_yx") == r.s[3, 0] and r.component("S12_xy") == r.s[0, 3]
    assert r.component("S11_xx") == r.s[0, 0] and r.component("S22_yy") == r.s[3, 3]
    assert r.magnitude("S21_xx") == abs(r.s[2, 0])
    assert r.phase_deg("S21_xx") == pytest.approx(np.degrees(np.angle(r.s[2, 0])))
    assert r.component("S21_LR") == pytest.approx(block_to_circular(family_block(r.s, "S21"))[1, 0])
    for bad in ("S31_xx", "S21", "S21_zz", "T21_xx", "S21_xR"):
        with pytest.raises(ValueError):
            r.component(bad)
    assert len(component_names("linear")) == len(component_names("circular")) == 16


def test_record_matrices_are_immutable_and_validated(db):
    r = db.records[0]
    with pytest.raises(ValueError):
        r.s[0, 0] = 0.0
    with pytest.raises(DatabaseError):
        Record("x", 0, 1e9, np.zeros((2, 2)), r.diagnostics)
    with pytest.raises(DatabaseError):
        Record("x", 0, 1e9, np.full((4, 4), np.nan), r.diagnostics)
    for bad in ({**r.diagnostics, "extra": 1}, {k: v for k, v in r.diagnostics.items() if k != "n_unknowns"},
                {**r.diagnostics, "specular_only": 1}, {**r.diagnostics, "n_unknowns": 1.5},
                {**r.diagnostics, "reciprocity_error": float("nan")}):
        with pytest.raises(DatabaseError):
            Record("x", 0, 1e9, r.s, bad)


def test_phase_is_stored_on_the_principal_interval(db, files):
    with open(files[1], newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        for name in component_names("linear") + component_names("circular"):
            assert -180.0 < float(row[f"{name}_phase_deg"]) <= 180.0
            assert float(row[f"{name}_mag"]) >= 0.0


# ------------------------------------------------------------------------- JSON round trip

def test_json_round_trip_is_bit_exact(db, files):
    back = Database.from_json(files[0])
    assert back == db
    for a, b in zip(db.records, back.records):
        assert np.array_equal(a.s, b.s) and a.frequency == b.frequency and a.diagnostics == b.diagnostics
    assert [c.parameters for c in back.cells] == [c.parameters for c in db.cells]
    assert back.cells[0].validation == db.cells[0].validation and back.cells[2].reason == db.cells[2].reason


def test_json_export_is_idempotent_byte_for_byte(db, files, tmp_path):
    back = Database.from_json(files[0])
    back.to_json(tmp_path / "again.json")
    assert (tmp_path / "again.json").read_bytes() == files[0].read_bytes()


def test_json_contains_both_value_and_derived_fields_for_all_sixteen_entries(files):
    doc = json.loads(files[0].read_text())
    assert doc["format"] == FORMAT_NAME and doc["version"] == 1
    rec = doc["cells"][0]["records"][0]
    for fam in ("S11", "S12", "S21", "S22"):
        assert set(rec["S"][fam]) == {"xx", "yx", "xy", "yy"}
        assert set(rec["S"][fam]["xx"]) == {"re", "im", "mag", "phase_deg"}
        assert set(rec["S_circular"][fam]) == {"RR", "LR", "RL", "LL"}
        assert set(rec["S_circular"][fam]["RR"]) == {"mag", "phase_deg"}
    z = complex(rec["S"]["S21"]["yx"]["re"], rec["S"]["S21"]["yx"]["im"])
    assert rec["S"]["S21"]["yx"]["mag"] == abs(z)


def test_environment_is_provenance_only_and_excluded_from_equality(db, files):
    doc = json.loads(files[0].read_text())
    assert set(doc["environment"]) == {"python", "numpy"}
    doc["environment"] = {"python": "9.9", "numpy": "0.0"}
    assert Database.from_dict(doc) == db


def test_rebuilt_unit_cells_reproduce_the_stored_s_parameters(db, files, result):
    """The strongest round trip: metadata + parameters are enough to recompute every stored value."""
    back = Database.from_json(files[0])
    cfg = result.database.metadata["context"]
    from src.basis_functions import MeshSettings
    from src.excitation import ReferencePlanes

    mesh = MeshSettings(**cfg["mesh"])
    planes = ReferencePlanes(cfg["reference_planes"]["top_air_thickness_m"],
                             cfg["reference_planes"]["bottom_air_thickness_m"], cfg["substrate"]["thickness_m"])
    for cell in back.cells:
        if cell.status != "ok":
            continue
        uc = back.unit_cell(cell.index)
        for rec in back.records_of(cell.index):
            fresh = compute_s_parameters(uc, rec.frequency, mesh, planes, warn_grating=False).matrix()
            assert np.max(np.abs(fresh - rec.s)) < 1e-12


def test_unit_cell_of_an_invalid_cell_or_csv_import_is_refused(db, files):
    with pytest.raises(DatabaseError, match="invalid"):
        db.unit_cell(2)
    csv_db = Database.from_csv(files[1])
    with pytest.raises(DatabaseError, match="no 'context'"):
        csv_db.unit_cell(0)


# -------------------------------------------------------------------------- CSV round trip

def test_csv_layout(db, files):
    with open(files[1], newline="") as fh:
        rows = list(csv.reader(fh))
    assert tuple(rows[0]) == CSV_COLUMNS and len(rows) == 1 + len(db.records)
    assert len(CSV_COLUMNS) == 4 + 5 + 7 + 2 * 16 + 2 * 16
    assert rows[0][:9] == ["cell_id", "cell_index", "frequency_index", "frequency_hz",
                           "length_m", "width_m", "center_x_m", "center_y_m", "rotation_deg"]
    assert len(set(r[0] for r in rows[1:])) == 4                       # skipped cells are not in the CSV


def test_csv_round_trip_preserves_values_to_rounding(db, files):
    back = Database.from_csv(files[1])
    assert len(back.records) == len(db.records)
    for a, b in zip(db.records, back.records):
        assert (a.cell_id, a.frequency_index, a.frequency) == (b.cell_id, b.frequency_index, b.frequency)
        assert np.max(np.abs(a.s - b.s)) < 1e-14
        assert a.diagnostics == b.diagnostics
        for name in component_names("linear"):
            assert b.magnitude(name) == pytest.approx(a.magnitude(name), abs=1e-14, rel=1e-13)
            if a.magnitude(name) > 1e-6:
                assert abs(wrapped_difference_deg(b.phase_deg(name), a.phase_deg(name))) < 1e-9
    assert [c.parameters for c in back.cells] == [c.parameters for c in db.cells if c.status == "ok"]


def test_csv_stored_columns_round_trip_exactly_as_text(db, files, tmp_path):
    """mag / phase_deg are written with repr: reading the numbers back gives the identical float."""
    with open(files[1], newline="") as fh:
        rows = list(csv.DictReader(fh))
    first = db.records[0]
    for name in ("S21_xx", "S11_yx", "S22_yy"):
        assert float(rows[0][f"{name}_mag"]) == first.magnitude(name)
        assert float(rows[0][f"{name}_phase_deg"]) == first.phase_deg(name)
    assert float(rows[0]["frequency_hz"]) == first.frequency


def test_csv_circular_columns_are_the_conversion_of_the_linear_ones(db, files):
    with open(files[1], newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row, rec in zip(rows, db.records):
        for fam in ("S11", "S12", "S21", "S22"):
            lin = np.zeros((2, 2), dtype=complex)
            for comp, (i, j) in {"xx": (0, 0), "yx": (1, 0), "xy": (0, 1), "yy": (1, 1)}.items():
                lin[i, j] = float(row[f"{fam}_{comp}_mag"]) * np.exp(1j * np.radians(float(row[f"{fam}_{comp}_phase_deg"])))
            circ = block_to_circular(lin)
            for comp, (i, j) in {"RR": (0, 0), "LR": (1, 0), "RL": (0, 1), "LL": (1, 1)}.items():
                assert float(row[f"{fam}_{comp}_mag"]) == pytest.approx(abs(circ[i, j]), abs=1e-13)


def test_csv_import_accepts_external_metadata(db, files):
    back = Database.from_csv(files[1], metadata=db.metadata)
    assert back.unit_cell(0).geometry.length == pytest.approx(3e-3)


def test_csv_export_is_deterministic_bytes(db, tmp_path):
    db.to_csv(tmp_path / "a.csv")
    db.to_csv(tmp_path / "b.csv")
    assert (tmp_path / "a.csv").read_bytes() == (tmp_path / "b.csv").read_bytes()


# ---------------------------------------------------------------------------- import errors

def load_doc(path):
    return json.loads(path.read_text())


def write_doc(doc, tmp_path, name="bad.json"):
    p = tmp_path / name
    p.write_text(json.dumps(doc))
    return p


def test_json_import_rejects_wrong_format_version_and_garbage(files, tmp_path):
    doc = load_doc(files[0])
    for change, fragment in (({"format": "other"}, "not a"), ({"version": 2}, "unsupported database version"),
                             ({"version": None}, "unsupported")):
        with pytest.raises(DatabaseError, match=fragment):
            Database.from_json(write_doc({**doc, **change}, tmp_path))
    for key in ("metadata", "cells"):
        bad = {k: v for k, v in doc.items() if k != key}
        with pytest.raises(DatabaseError, match="missing"):
            Database.from_json(write_doc(bad, tmp_path))
    (tmp_path / "garbage.json").write_text("{not json")
    with pytest.raises(DatabaseError, match="invalid JSON"):
        Database.from_json(tmp_path / "garbage.json")
    (tmp_path / "nan.json").write_text('{"format": "x", "a": NaN}')
    with pytest.raises(DatabaseError, match="non-finite"):
        Database.from_json(tmp_path / "nan.json")
    with pytest.raises(DatabaseError, match="not a"):
        Database.from_dict([1, 2, 3])


def test_json_import_detects_edited_parameters_via_cell_id(files, tmp_path):
    doc = load_doc(files[0])
    doc["cells"][0]["parameters"]["length_m"] = 3.5e-3
    with pytest.raises(DatabaseError, match="does not match its parameters"):
        Database.from_json(write_doc(doc, tmp_path))


def test_json_import_detects_inconsistent_derived_fields(files, tmp_path):
    doc = load_doc(files[0])
    doc["cells"][0]["records"][0]["S"]["S21"]["xx"]["mag"] += 0.01
    with pytest.raises(DatabaseError, match="S21_xx.*magnitude"):
        Database.from_json(write_doc(doc, tmp_path))
    doc = load_doc(files[0])
    entry = doc["cells"][0]["records"][0]["S"]["S21"]["xx"]
    entry["phase_deg"] = float(wrapped_difference_deg(entry["phase_deg"], -20.0))
    with pytest.raises(DatabaseError, match="phase"):
        Database.from_json(write_doc(doc, tmp_path))
    doc = load_doc(files[0])
    doc["cells"][0]["records"][0]["S_circular"]["S11"]["LR"]["mag"] += 0.5
    with pytest.raises(DatabaseError, match="S11_LR"):
        Database.from_json(write_doc(doc, tmp_path))


def test_json_import_structural_errors(files, tmp_path):
    doc = load_doc(files[0])
    dup = copy.deepcopy(doc)
    dup["cells"][0]["records"].append(dup["cells"][0]["records"][0])
    with pytest.raises(DatabaseError, match="duplicate record"):
        Database.from_json(write_doc(dup, tmp_path))
    twin = copy.deepcopy(doc)
    twin["cells"][1] = copy.deepcopy(twin["cells"][0])
    with pytest.raises(DatabaseError, match="duplicate"):
        Database.from_json(write_doc(twin, tmp_path))
    bad_status = copy.deepcopy(doc)
    bad_status["cells"][2]["records"] = copy.deepcopy(doc["cells"][0]["records"])
    with pytest.raises(DatabaseError, match="not 'ok'"):
        Database.from_json(write_doc(bad_status, tmp_path))
    missing = copy.deepcopy(doc)
    del missing["cells"][0]["records"][0]["S"]["S21"]
    with pytest.raises(DatabaseError, match="malformed"):
        Database.from_json(write_doc(missing, tmp_path))
    no_diag = copy.deepcopy(doc)
    no_diag["cells"][0]["records"][0]["diagnostics"].pop("n_unknowns")
    with pytest.raises(DatabaseError, match="diagnostics"):
        Database.from_json(write_doc(no_diag, tmp_path))
    unordered = copy.deepcopy(doc)
    unordered["cells"][0], unordered["cells"][1] = unordered["cells"][1], unordered["cells"][0]
    with pytest.raises(DatabaseError, match="increasing index"):
        Database.from_json(write_doc(unordered, tmp_path))


def test_csv_import_errors(files, tmp_path):
    lines = files[1].read_text().splitlines()
    cases = {
        "header": ["cell_id,x"] + lines[1:],
        "short row": [lines[0], lines[1].rsplit(",", 1)[0]],
        "not a number": [lines[0], lines[1].replace(lines[1].split(",")[3], "abc", 1)],
        "empty": [],
    }
    for label, content in cases.items():
        p = tmp_path / f"{label.replace(' ', '_')}.csv"
        p.write_text("\n".join(content) + "\n")
        with pytest.raises(DatabaseError):
            Database.from_csv(p)


def test_csv_import_rejects_negative_magnitude_and_tampered_circular_columns(files, tmp_path):
    lines = files[1].read_text().splitlines()
    header = lines[0].split(",")
    row = lines[1].split(",")
    neg = row.copy()
    neg[header.index("S21_xx_mag")] = "-0.5"
    (tmp_path / "neg.csv").write_text("\n".join([lines[0], ",".join(neg)]) + "\n")
    with pytest.raises(DatabaseError, match="negative magnitude"):
        Database.from_csv(tmp_path / "neg.csv")
    tam = row.copy()
    tam[header.index("S21_RR_mag")] = "0.123"
    (tmp_path / "tam.csv").write_text("\n".join([lines[0], ",".join(tam)]) + "\n")
    with pytest.raises(DatabaseError, match="S21_RR"):
        Database.from_csv(tmp_path / "tam.csv")
    inconsistent = row.copy()
    inconsistent[header.index("length_m")] = "0.0042"
    (tmp_path / "inc.csv").write_text("\n".join([lines[0], lines[1], ",".join(inconsistent)]) + "\n")
    with pytest.raises(DatabaseError):
        Database.from_csv(tmp_path / "inc.csv")


# -------------------------------------------------------------------------------- queries

def test_series_and_phase_series(db):
    f, s = db.series(0, "S21_xx")
    assert f.tolist() == [8e9, 9e9, 10e9] and s.shape == (3,)
    f2, ph = db.phase_series_deg(0, "S21_xx")
    assert np.allclose(ph, principal_phase_deg(s))
    _, un = db.phase_series_deg(0, "S21_xx", unwrap=True)
    assert np.allclose(np.mod(un - ph + 180, 360) - 180, 0.0, atol=1e-9)       # same modulo 360
    with pytest.raises(KeyError):
        db.cell(99)
    with pytest.raises(KeyError):
        db.cell("nope")
    assert db.cell(0) is db.cell(db.cells[0].cell_id)
    assert db.record(0, 1).frequency == 9e9
    with pytest.raises(KeyError):
        db.record(0, 7)


def test_unwrap_requires_consecutive_frequency_indices(db):
    cell = db.cells[0]
    recs = [r for r in db.records if r.cell_id == cell.cell_id and r.frequency_index != 1]
    gappy = Database(db.metadata, [cell], recs)
    with pytest.raises(DatabaseError, match="consecutive"):
        gappy.phase_series_deg(0, "S21_xx", unwrap=True)
    assert gappy.phase_series_deg(0, "S21_xx")[0].size == 2          # principal values are fine


def test_s_array_layout_and_missing_entries(db):
    lin = db.s_array()
    assert lin.shape == (6, 3, 4, 4)
    assert np.all(np.isnan(lin[2])) and np.all(np.isnan(lin[3]))     # skipped cells
    assert np.array_equal(lin[0, 1], db.record(0, 1).s)
    circ = db.s_array("circular")
    assert np.allclose(circ[4, 2, 2:, :2], block_to_circular(family_block(db.record(4, 2).s, "S21")), atol=1e-14)
    assert np.all(np.isnan(circ[2]))
    with pytest.raises(ValueError):
        db.s_array("elliptical")


def test_parameter_array_and_failed_validation(db):
    assert np.allclose(db.parameter_array("length_m"), [3e-3, 3e-3, 9e-3, 9e-3, 4e-3, 4e-3])
    assert np.allclose(db.parameter_array("rotation_deg"), [0, 30, 0, 30, 0, 30])
    with pytest.raises(ValueError):
        db.parameter_array("period_x")
    assert db.failed_validation() == []


def test_failed_validation_reports_failed_refinement_and_failed_checks(db):
    cell = db.cells[0]
    bad_ref = CellEntry(cell.index, cell.cell_id, cell.parameters, "ok", None,
                        {**cell.validation, "refinement": {**cell.validation["refinement"], "passed": False}})
    tweaked = Database(db.metadata, [bad_ref] + list(db.cells[1:]), db.records)
    assert [c.cell_id for c in tweaked.failed_validation()] == [cell.cell_id]
    worse = [Record(r.cell_id, r.frequency_index, r.frequency, r.s, {**r.diagnostics, "checks_passed": False})
             if r.cell_id == db.cells[1].cell_id else r for r in db.records]
    assert [c.cell_id for c in Database(db.metadata, db.cells, worse).failed_validation()] == [db.cells[1].cell_id]


def test_database_constructor_validates_structure(db):
    r = db.records[0]
    with pytest.raises(DatabaseError, match="unknown cell_id"):
        Database({}, db.cells[:1], [Record("zzz", 0, 1e9, r.s, r.diagnostics)])
    with pytest.raises(DatabaseError, match="duplicate cell_id"):
        Database({}, [db.cells[0], db.cells[0]], [])
    assert Database({}, [], []).s_array().shape == (0, 0, 4, 4)


def test_component_value_helper_direct():
    m = np.arange(16, dtype=complex).reshape(4, 4)
    assert component_value(m, "S12_yx") == m[1, 2] and component_value(m, "S21_xy") == m[2, 1]
