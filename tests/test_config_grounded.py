"""Phase 12A: ``substrate.backing`` config key and the conditional ``reference_planes`` schema."""

import copy

import pytest

from src.config import ConfigError, config_from_dict

BASE = {
    "period_x": 4.8e-3,
    "period_y": 4.8e-3,
    "substrate": {"name": "RT5880", "eps_r": 2.2, "tan_delta": 0.0009, "thickness": 1.575e-3},
    "geometry": {"type": "rectangle", "length": 2.0e-3, "width": 2.0e-3},
    "excitation": {"frequency": 12.5e9},
    "mesh": {"cells_x": 8, "cells_y": 8, "floquet_order": 16},
    "reference_planes": {"top_air_thickness": 5.0e-3, "bottom_air_thickness": 5.0e-3},
}


def _grounded_base():
    data = copy.deepcopy(BASE)
    data["substrate"]["backing"] = "pec"
    data["reference_planes"] = {"top_air_thickness": 5.0e-3}
    return data


def test_absent_backing_key_gives_air_unchanged():
    cfg = config_from_dict(copy.deepcopy(BASE))
    assert cfg.unit_cell.backing == "air"
    assert cfg.reference_planes.bottom_air_thickness == 5.0e-3


def test_backing_pec_with_bottom_air_thickness_omitted():
    cfg = config_from_dict(_grounded_base())
    assert cfg.unit_cell.backing == "pec"
    assert cfg.reference_planes.bottom_air_thickness is None
    assert cfg.reference_planes.grounded is True


def test_backing_pec_rejects_bottom_air_thickness_present():
    data = _grounded_base()
    data["reference_planes"] = {"top_air_thickness": 5.0e-3, "bottom_air_thickness": 5.0e-3}
    with pytest.raises(ConfigError, match="bottom_air_thickness"):
        config_from_dict(data)


def test_backing_air_still_requires_bottom_air_thickness():
    data = copy.deepcopy(BASE)
    data["reference_planes"] = {"top_air_thickness": 5.0e-3}
    with pytest.raises(ConfigError, match="bottom_air_thickness"):
        config_from_dict(data)


def test_invalid_backing_value_rejected():
    data = copy.deepcopy(BASE)
    data["substrate"]["backing"] = "copper"
    with pytest.raises(ConfigError, match="backing"):
        config_from_dict(data)


def test_backing_is_not_a_valid_substrate_key_typo_still_caught():
    data = copy.deepcopy(BASE)
    data["substrate"]["backig"] = "pec"  # typo
    with pytest.raises(ConfigError, match="unknown key"):
        config_from_dict(data)
