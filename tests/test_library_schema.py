"""Phase 7 locked configuration: `LibraryConfig.default()` matches the brief's Sec. 1 exactly,
sourced from this project's own authoritative `config.yaml` where the brief allows it; validation
rejects an out-of-range sweep rather than silently truncating it."""

import math

import pytest
import yaml

from src.config import config_from_dict
from src.library_schema import (
    DEFAULT_MESH,
    L_VALUES_MM,
    W_VALUES_MM,
    LibraryConfig,
    _AUTHORITATIVE_CONFIG_PATH,
)

ROOT_CONFIG = yaml.safe_load(_AUTHORITATIVE_CONFIG_PATH.read_text())


def test_default_matches_locked_period_and_frequency_grid():
    cfg = LibraryConfig.default()
    assert cfg.period_m == pytest.approx(4.8e-3)
    assert cfg.f_min_hz == pytest.approx(10.0e9)
    assert cfg.f_max_hz == pytest.approx(15.0e9)
    assert cfg.f_step_hz == pytest.approx(25.0e6)
    assert cfg.f_design_hz == pytest.approx(12.5e9)


def test_default_substrate_matches_project_authoritative_config_yaml():
    """The brief's own locked substrate values are not invented here -- they are read from
    (and checked against) this project's root config.yaml, its authoritative source."""
    cfg = LibraryConfig.default()
    sub = ROOT_CONFIG["substrate"]
    assert cfg.substrate["name"] == sub["name"]
    assert cfg.substrate["eps_r"] == pytest.approx(sub["eps_r"])
    assert cfg.substrate["tan_delta"] == pytest.approx(sub["tan_delta"])
    assert cfg.substrate["thickness_m"] == pytest.approx(sub["thickness"])
    # And the brief's own locked value (1.575 mm) is exactly this, not a coincidence to be lost:
    assert cfg.substrate["thickness_m"] == pytest.approx(1.575e-3)


def test_frequency_grid_has_exactly_201_points_10_to_15_ghz():
    cfg = LibraryConfig.default()
    n = round((cfg.f_max_hz - cfg.f_min_hz) / cfg.f_step_hz) + 1
    assert n == 201
    freqs = [cfg.f_min_hz + i * cfg.f_step_hz for i in range(n)]
    assert freqs[0] == pytest.approx(10.0e9)
    assert freqs[-1] == pytest.approx(15.0e9)
    assert any(abs(f - 12.5e9) < 1.0 for f in freqs)   # exact grid point, not merely "close"


def test_l_w_sweep_is_18_values_1p0_to_4p4_step_0p2():
    assert len(L_VALUES_MM) == 18 and len(W_VALUES_MM) == 18
    assert L_VALUES_MM[0] == 1.0
    assert L_VALUES_MM[-1] == 4.4
    steps = [round(L_VALUES_MM[i + 1] - L_VALUES_MM[i], 10) for i in range(len(L_VALUES_MM) - 1)]
    assert all(s == 0.2 for s in steps)
    assert L_VALUES_MM == W_VALUES_MM


def test_default_has_324_geometries():
    cfg = LibraryConfig.default()
    assert cfg.n_geometries == 324
    assert len(cfg.l_values_mm) * len(cfg.w_values_mm) == 324


def test_default_mesh_is_32x32_floquet_32():
    """The mesh choice actually used (`library_schema.DEFAULT_MESH`), not a placeholder -- see
    this module's docstring for why 24x24 (which aligns with the 0.2 mm sweep step) was
    rejected: tests/test_library_builder.py::test_24x24_mesh_has_large_area_error documents the
    measured floating-point pathology directly."""
    cfg = LibraryConfig.default()
    assert cfg.mesh == DEFAULT_MESH == {"cells_x": 32, "cells_y": 32, "floquet_order": 32}
    # Rank rule (MOM_FORMULATION.md): 2*floquet_order + 1 >= max(cells_x, cells_y).
    assert 2 * cfg.mesh["floquet_order"] + 1 >= max(cfg.mesh["cells_x"], cfg.mesh["cells_y"])


def test_default_builds_a_valid_simulation_config():
    cfg = LibraryConfig.default()
    sim = cfg.build_config()
    assert sim.unit_cell.lattice.period_x == pytest.approx(4.8e-3)
    assert sim.unit_cell.substrate.eps_r == pytest.approx(2.2)
    assert sim.mesh.cells_x == 32 and sim.mesh.floquet_order == 32


def test_every_l_w_value_strictly_inside_the_period():
    cfg = LibraryConfig.default()
    for l_mm in cfg.l_values_mm:
        assert 0.0 < l_mm * 1e-3 < cfg.period_m
    for w_mm in cfg.w_values_mm:
        assert 0.0 < w_mm * 1e-3 < cfg.period_m


def test_rejects_l_w_that_does_not_fit_the_period():
    with pytest.raises(ValueError, match="period"):
        LibraryConfig(
            period_m=4.8e-3,
            substrate={"name": "x", "eps_r": 2.2, "tan_delta": 0.0, "thickness_m": 1.575e-3},
            f_min_hz=10.0e9, f_max_hz=15.0e9, f_step_hz=25.0e6, f_design_hz=12.5e9,
            l_values_mm=(1.0, 5.0), w_values_mm=(1.0, 2.0),   # 5.0 mm >= 4.8 mm period
        )


def test_rejects_design_frequency_outside_the_band():
    with pytest.raises(ValueError):
        LibraryConfig(
            period_m=4.8e-3,
            substrate={"name": "x", "eps_r": 2.2, "tan_delta": 0.0, "thickness_m": 1.575e-3},
            f_min_hz=10.0e9, f_max_hz=15.0e9, f_step_hz=25.0e6, f_design_hz=20.0e9,
            l_values_mm=(1.0,), w_values_mm=(1.0,),
        )
