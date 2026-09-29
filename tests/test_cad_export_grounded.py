"""Phase 12A: the ground-plane solid in ``cad_export.py``."""

import pytest

cadquery = pytest.importorskip("cadquery")

from src.cad_export import (
    build_ground_solid,
    build_unit_cell_solids,
    export_unit_cell_step,
    import_and_summarize,
)
from src.config import config_from_dict

BASE = {
    "period_x": 4.8e-3,
    "period_y": 4.8e-3,
    "substrate": {"name": "RT5880", "eps_r": 2.2, "tan_delta": 0.0009, "thickness": 1.575e-3},
    "geometry": {"type": "rectangle", "length": 2.0e-3, "width": 2.0e-3},
    "excitation": {"frequency": 12.5e9},
    "mesh": {"cells_x": 8, "cells_y": 8, "floquet_order": 16},
    "reference_planes": {"top_air_thickness": 5.0e-3, "bottom_air_thickness": 5.0e-3},
}


def _cfg(backing):
    import copy

    data = copy.deepcopy(BASE)
    if backing == "pec":
        data["substrate"]["backing"] = "pec"
        data["reference_planes"] = {"top_air_thickness": 5.0e-3}
    return config_from_dict(data)


def test_build_ground_solid_none_for_ungrounded():
    assert build_ground_solid(_cfg("air").unit_cell) is None


def test_build_unit_cell_solids_still_returns_a_2_tuple_when_grounded():
    substrate, metal = build_unit_cell_solids(_cfg("pec").unit_cell)
    assert substrate is not None and metal is not None


def test_ground_solid_sits_flush_below_the_substrate():
    uc = _cfg("pec").unit_cell
    ground = build_ground_solid(uc)
    from src.cad_export import SolidProperties

    props = SolidProperties.of(ground)
    h_mm = 1.575
    assert props.bbox_mm[5] == pytest.approx(-h_mm, abs=1e-9)  # zmax flush with substrate bottom
    assert props.bbox_mm[4] < props.bbox_mm[5]  # extends further down


def test_export_and_reimport_grounded_step_recovers_three_solids(tmp_path):
    uc = _cfg("pec").unit_cell
    path = tmp_path / "cell.step"
    result = export_unit_cell_step(uc, path)
    assert result.ground is not None

    summary = import_and_summarize(path)
    assert summary.ground is not None
    h_mm = 1.575
    assert summary.substrate.bbox_mm[4] == pytest.approx(-h_mm, abs=1e-6)  # substrate zmin = -h
    assert summary.ground.bbox_mm[5] == pytest.approx(-h_mm, abs=1e-6)  # ground zmax = -h
    assert summary.metal.bbox_mm[4] == pytest.approx(0.0, abs=1e-6)  # metal zmin = 0


def test_export_and_reimport_ungrounded_step_still_has_two_solids(tmp_path):
    uc = _cfg("air").unit_cell
    path = tmp_path / "cell.step"
    result = export_unit_cell_step(uc, path)
    assert result.ground is None
    summary = import_and_summarize(path)
    assert summary.ground is None
