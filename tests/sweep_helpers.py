"""Shared builders for the Phase 5 tests (a small, fast, mesh-aligned unit cell)."""

import copy

from src.config import config_from_dict

# 8 mm cell, 16 x 16 mesh (0.5 mm): feature edges at multiples of 1 mm lie exactly on mesh lines,
# so the staircase area error is 0 and a solve takes ~10 ms.
BASE = {
    "period_x": 8.0e-3,
    "period_y": 8.0e-3,
    "substrate": {"name": "test", "eps_r": 2.2, "tan_delta": 0.0009, "thickness": 1.575e-3},
    "geometry": {"type": "rectangle", "length": 5.0e-3, "width": 2.0e-3},
    "excitation": {"frequency": 10.0e9},
    "mesh": {"cells_x": 16, "cells_y": 16, "floquet_order": 32},
    "reference_planes": {"top_air_thickness": 0.0, "bottom_air_thickness": 0.0},
}


def config_dict(sweep=None, **overrides):
    """BASE with dotted overrides ('geometry.length': 3e-3) and an optional 'sweep' section."""
    d = copy.deepcopy(BASE)
    for path, value in overrides.items():
        keys = path.split(".")
        node = d
        for k in keys[:-1]:
            node = node[k]
        node[keys[-1]] = value
    if sweep is not None:
        d["sweep"] = copy.deepcopy(sweep)
    return d


def sweep_config(sweep=None, **overrides):
    return config_from_dict(config_dict(sweep, **overrides))


FREQ_8_12 = {"start": "8 GHz", "stop": "12 GHz", "step": "1 GHz"}
