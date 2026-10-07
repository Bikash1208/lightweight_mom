"""Phase 13 Part 3: `library_sweep.py` generalizes `library_builder.py`'s locked-to-Phase-7 sweep
to an arbitrary config -- in particular a GROUNDED one, which `library_schema.LibraryConfig` has
no field for at all (CLAUDE.md: "Only one library is actually built, still ungrounded")."""

import json

import numpy as np
import pytest

from src.config import config_from_dict
from src.library_sweep import (
    SweepSpec,
    build_annular_slot_library_cells,
    build_square_library_cells,
    classify_cell_backing_aware,
    detect_resonance_phase_crossing,
    extract_design_frequency,
    run_custom_sweep,
    write_annular_slot_library_json,
    write_app_library_json,
    write_square_library_json,
)
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


# ------------------------------------------------------------- PHASE13B A6: grounded resonance


def test_detect_resonance_phase_crossing_finds_a_synthetic_branch_cut_jump():
    """A grounded cell's resonance is the reflection phase crossing +-180deg, essentially at
    constant |S11| -- `library_classifier.detect_resonance`'s magnitude-peak criterion cannot see
    this at all."""
    freqs = np.linspace(12.0e9, 13.0e9, 5)
    phase = np.array([-140.0, -160.0, -178.0, 175.0, 155.0])  # wraps between index 2 and 3
    result = detect_resonance_phase_crossing(freqs, phase)
    assert result.is_resonant is True
    assert result.resonance_frequency_hz == pytest.approx(0.5 * (freqs[2] + freqs[3]))


def test_detect_resonance_phase_crossing_reports_none_found_when_phase_never_wraps():
    freqs = np.linspace(12.0e9, 13.0e9, 5)
    phase = np.array([-140.0, -145.0, -150.0, -155.0, -160.0])  # smooth drift, no branch-cut jump
    result = detect_resonance_phase_crossing(freqs, phase)
    assert result.is_resonant is False


@pytest.mark.slow
def test_classify_cell_backing_aware_finds_the_grounded_design_point_resonance():
    """L=W=6.0mm, 8mm period, 12.5 GHz is this project's own recorded grounded resonance
    (PHASE12D_ADDENDUM.md/phase13/convergence.py); classify_cell_backing_aware must report it as
    resonant, unlike the ungrounded-tuned `classify_cell` (PHASE13B A6, live review). Real-measured
    sweep confirms |S11_xx| stays within 0.9957-0.9959 across the whole band while the phase falls
    smoothly through +-180deg between 12.6 and 12.7 GHz -- a magnitude-based criterion cannot see
    this design point's resonance at any mesh/frequency resolution."""
    cfg = config_from_dict(
        {
            "period_x": 8.0e-3,
            "period_y": 8.0e-3,
            "substrate": {**_SUBSTRATE, "backing": "pec"},
            "geometry": {"type": "rectangle", "length": 6.0e-3, "width": 6.0e-3},
            "excitation": {"frequency": 12.5e9},
            "mesh": {"cells_x": 40, "cells_y": 40, "floquet_order": 40},
            "reference_planes": {"top_air_thickness": 5.0e-3},
        }
    )
    spec = SweepSpec(
        f_min_hz=12.0e9,
        f_max_hz=13.0e9,
        f_step_hz=0.1e9,
        f_design_hz=12.5e9,
        l_values_mm=(6.0,),
        w_values_mm=(6.0,),
    )
    result = run_custom_sweep(cfg, spec, out_dir=None)
    rows = extract_design_frequency(result.database, spec)
    cell = result.database.cells[0]

    from src.library_classifier import classify_cell

    old_classification = classify_cell(cell, result.database, rows[0])
    new_classification = classify_cell_backing_aware(cell, result.database, rows[0], cfg.unit_cell)
    assert old_classification.resonance.is_resonant is False  # the pre-existing bug, still true
    assert new_classification.resonance.is_resonant is True


# --------------------------------------------------------------- PHASE13C: auto-built board library


def test_build_square_library_cells_solves_only_the_diagonal():
    """5 L values must give exactly 5 cells (all L=W), not the 25 a full L x W cross product
    through build_sweep_definition's two independent axes would solve."""
    base_config = {
        "period_x": 8.0e-3,
        "period_y": 8.0e-3,
        "substrate": {**_SUBSTRATE, "backing": "pec"},
        "geometry": {"type": "rectangle", "length": 2.0e-3, "width": 2.0e-3},  # placeholder
        "excitation": {"frequency": 12.5e9},  # placeholder
        "mesh": _MESH,
        "reference_planes": {"top_air_thickness": 5.0e-3},
    }
    calls = []
    cells = build_square_library_cells(
        base_config, l_values_mm=(2.0, 3.0, 4.0), frequency_hz=12.5e9, progress=calls.append
    )
    assert len(cells) == 3
    for cell, l_mm in zip(cells, (2.0, 3.0, 4.0), strict=True):
        assert cell["L_mm"] == pytest.approx(l_mm)
        assert cell["W_mm"] == pytest.approx(l_mm)
        assert cell["solver_status"] == "OK"
    assert calls[0] == {"done": 0, "total": 3}
    assert calls[-1] == {"done": 3, "total": 3}


def test_write_square_library_json_round_trips_through_unit_cell_library_query_api(tmp_path):
    base_config = {
        "period_x": 8.0e-3,
        "period_y": 8.0e-3,
        "substrate": {**_SUBSTRATE, "backing": "pec"},
        "geometry": {"type": "rectangle", "length": 2.0e-3, "width": 2.0e-3},
        "excitation": {"frequency": 12.5e9},
        "mesh": _MESH,
        "reference_planes": {"top_air_thickness": 5.0e-3},
    }
    cells = build_square_library_cells(base_config, l_values_mm=(2.0, 3.0), frequency_hz=12.5e9)
    path = write_square_library_json(cells, tmp_path / "square_lib.json")
    doc = json.loads(path.read_text())
    assert len(doc["cells"]) == 2

    lib = UnitCellLibrary.from_json(path)
    entry = lib.response(2.0, 2.0, 12.5, family="S11")
    assert set(entry) == {"real", "imag", "mag", "mag_db", "phase_deg"}


# ======================================================== PHASE15 Part E: annular-slot patch library


def _annular_base_config():
    return {
        "period_x": 6.0e-3,
        "period_y": 6.0e-3,
        "substrate": {**_SUBSTRATE, "backing": "pec"},
        "geometry": {"type": "patch_annular_slot", "a": 4.0e-3, "b": 4.0e-3, "R": 1.5e-3, "r": 1.0e-3},  # placeholder
        "excitation": {"frequency": 12.5e9},  # placeholder
        "mesh": _MESH,
        "reference_planes": {"top_air_thickness": 5.0e-3},
    }


def test_build_annular_slot_library_cells_sweeps_only_the_outer_radius():
    calls = []
    cells = build_annular_slot_library_cells(
        _annular_base_config(),
        a_mm=4.0,
        b_mm=4.0,
        r_mm=1.0,
        R_values_mm=[1.2, 1.5, 1.8],
        frequency_hz=12.5e9,
        progress=calls.append,
    )
    assert len(cells) == 3
    for cell, outer_r_mm in zip(cells, (1.2, 1.5, 1.8), strict=True):
        assert cell["a_mm"] == pytest.approx(4.0)
        assert cell["b_mm"] == pytest.approx(4.0)
        assert cell["r_mm"] == pytest.approx(1.0)
        assert cell["R_mm"] == pytest.approx(outer_r_mm)
        assert cell["solver_status"] == "OK"
        assert abs(cell["S11"][0]["mag"] - 1.0) < 0.05  # lossless grounded PEC reflector
    assert calls[0] == {"done": 0, "total": 3}
    assert calls[-1] == {"done": 3, "total": 3}


def test_annular_slot_phase_differs_from_a_plain_patch_of_the_same_outer_size():
    """A real measurement, not an assumption: the slot must actually change the reflection phase
    relative to a plain (no-slot) patch of the same outer footprint."""
    plain_config = {
        "period_x": 6.0e-3,
        "period_y": 6.0e-3,
        "substrate": {**_SUBSTRATE, "backing": "pec"},
        "geometry": {"type": "rectangle", "length": 4.0e-3, "width": 4.0e-3},
        "excitation": {"frequency": 12.5e9},
        "mesh": _MESH,
        "reference_planes": {"top_air_thickness": 5.0e-3},
    }
    plain_cells = build_square_library_cells(plain_config, l_values_mm=(4.0,), frequency_hz=12.5e9)
    annular_cells = build_annular_slot_library_cells(
        _annular_base_config(), a_mm=4.0, b_mm=4.0, r_mm=1.0, R_values_mm=[1.5], frequency_hz=12.5e9
    )
    plain_phase = plain_cells[0]["S11"][0]["phase_deg"]
    annular_phase = annular_cells[0]["S11"][0]["phase_deg"]
    assert abs(((annular_phase - plain_phase + 180.0) % 360.0) - 180.0) > 1.0


def test_write_annular_slot_library_json_round_trips(tmp_path):
    cells = build_annular_slot_library_cells(
        _annular_base_config(), a_mm=4.0, b_mm=4.0, r_mm=1.0, R_values_mm=[1.2, 1.5], frequency_hz=12.5e9
    )
    path = write_annular_slot_library_json(cells, tmp_path / "annular_lib.json")
    doc = json.loads(path.read_text())
    assert len(doc["cells"]) == 2
    assert doc["cells"][0]["R_mm"] == pytest.approx(1.2)
