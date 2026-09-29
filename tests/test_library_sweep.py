"""Phase 13 Part 3: `library_sweep.py` generalizes `library_builder.py`'s locked-to-Phase-7 sweep
to an arbitrary config -- in particular a GROUNDED one, which `library_schema.LibraryConfig` has
no field for at all (CLAUDE.md: "Only one library is actually built, still ungrounded")."""

import json

import pytest

from src.config import config_from_dict
from src.library_sweep import SweepSpec, extract_design_frequency, run_custom_sweep, write_app_library_json
from src.unit_cell_library import UnitCellLibrary

_SUBSTRATE = {"name": "Rogers RT/Duroid 5880", "eps_r": 2.2, "tan_delta": 0.0009, "thickness": 1.575e-3}
_MESH = {"cells_x": 8, "cells_y": 8, "floquet_order": 8}


def _ungrounded_config():
    return config_from_dict(
        {
            "period_x": 4.8e-3,
            "period_y": 4.8e-3,
            "substrate": _SUBSTRATE,
            "geometry": {"type": "rectangle", "length": 2.0e-3, "width": 2.0e-3},
            "excitation": {"frequency": 12.5e9},
            "mesh": _MESH,
            "reference_planes": {"top_air_thickness": 5.0e-3, "bottom_air_thickness": 5.0e-3},
        }
    )


def _grounded_config():
    return config_from_dict(
        {
            "period_x": 8.0e-3,
            "period_y": 8.0e-3,
            "substrate": {**_SUBSTRATE, "backing": "pec"},
            "geometry": {"type": "rectangle", "length": 2.0e-3, "width": 2.0e-3},
            "excitation": {"frequency": 12.5e9},
            "mesh": _MESH,
            "reference_planes": {"top_air_thickness": 5.0e-3},
        }
    )


_SPEC = SweepSpec(
    f_min_hz=12.0e9,
    f_max_hz=13.0e9,
    f_step_hz=1.0e9,
    f_design_hz=12.0e9,
    l_values_mm=(1.0, 2.0),
    w_values_mm=(2.0,),
)


def test_ungrounded_sweep_produces_all_four_families():
    result = run_custom_sweep(_ungrounded_config(), _SPEC, out_dir=None)
    assert len(result.database.cells) == 2
    for cell in result.database.cells:
        assert cell.status == "ok"
        records = result.database.records_of(cell.index)
        assert records[0].port_topology == "two_port"


def test_grounded_sweep_only_has_s11_and_does_not_crash_on_other_families():
    result = run_custom_sweep(_grounded_config(), _SPEC, out_dir=None)
    assert len(result.database.cells) == 2
    for cell in result.database.cells:
        records = result.database.records_of(cell.index)
        assert records[0].port_topology == "one_port"
        with pytest.raises(ValueError):
            records[0].component("S21_xx")


def test_write_app_library_json_round_trips_through_unit_cell_library_query_api(tmp_path):
    result = run_custom_sweep(_ungrounded_config(), _SPEC, out_dir=None)
    path = write_app_library_json(result.database, tmp_path / "lib.json")
    doc = json.loads(path.read_text())
    assert doc["format"] == "lightweight_mom.app_unit_cell_library"

    lib = UnitCellLibrary.from_json(path)
    entry = lib.response(1.0, 2.0, 12.0, family="S11")
    assert set(entry) == {"real", "imag", "mag", "mag_db", "phase_deg"}
    candidates = lib.find_phase_candidates(entry["phase_deg"], 12.0, tolerance_deg=1.0)
    assert any(c.l_mm == 1.0 and c.w_mm == 2.0 for c in candidates)


def test_extract_design_frequency_grounded_cell_has_only_s11():
    result = run_custom_sweep(_grounded_config(), _SPEC, out_dir=None)
    rows = extract_design_frequency(result.database, _SPEC)
    assert len(rows) == 2
    for row in rows:
        assert row.solver_status == "OK"
        assert row.s11 is not None
        assert row.s21 is None and row.s12 is None and row.s22 is None


def test_extract_design_frequency_ungrounded_cell_has_all_four_families():
    result = run_custom_sweep(_ungrounded_config(), _SPEC, out_dir=None)
    rows = extract_design_frequency(result.database, _SPEC)
    assert len(rows) == 2
    for row in rows:
        assert row.s11 is not None and row.s21 is not None


def test_write_app_library_json_grounded_cell_has_no_s21_key(tmp_path):
    result = run_custom_sweep(_grounded_config(), _SPEC, out_dir=None)
    path = write_app_library_json(result.database, tmp_path / "lib_grounded.json")
    doc = json.loads(path.read_text())
    for cell in doc["cells"]:
        assert "S11" in cell
        assert "S21" not in cell
        assert cell["backing"] == "pec"
