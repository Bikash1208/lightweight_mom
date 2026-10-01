"""PHASE13C: server-rendered PNG plots -- just checks a real PNG file is written (styling itself
isn't tested, matching this project's convention of not asserting pixel content). Needs the
optional 'app' extra (matplotlib); skipped entirely otherwise."""

import numpy as np
import pytest

pytest.importorskip("matplotlib")

from src.app import plots

_PNG_MAGIC = b"\x89PNG\r\n\x1a\n"


def _grid(n=6, extent_m=0.02):
    coords = np.linspace(-extent_m, extent_m, n)
    return np.meshgrid(coords, coords, indexing="xy")


def _assert_real_png(path):
    assert path.exists()
    assert path.stat().st_size > 0
    assert path.read_bytes()[:8] == _PNG_MAGIC


def test_reference_phase_map_png(tmp_path):
    x_m, y_m = _grid()
    phase_deg = np.degrees(np.arctan2(y_m, x_m))
    path = plots.reference_phase_map_png(x_m, y_m, phase_deg, tmp_path / "phase.png")
    _assert_real_png(path)


def test_layout_map_png(tmp_path):
    x_m, y_m = _grid()
    value = np.hypot(x_m, y_m) * 1e3
    path = plots.layout_map_png(x_m, y_m, value, tmp_path / "layout.png", label="L (mm)")
    _assert_real_png(path)


def test_pattern_uv_map_png(tmp_path):
    theta = np.linspace(0.0, 90.0, 20)
    phi = np.linspace(0.0, 360.0, 24, endpoint=False)
    theta_grid, phi_grid = np.meshgrid(theta, phi, indexing="ij")
    power_db = -np.abs(theta_grid - 20.0)  # a fake beam near theta=20deg
    path = plots.pattern_uv_map_png(
        theta_grid, phi_grid, power_db, tmp_path / "pattern.png", beams_theta_phi_deg=[(20.0, 0.0)]
    )
    _assert_real_png(path)


def test_pattern_cut_png(tmp_path):
    theta = np.linspace(0.0, 90.0, 20)
    phi = np.linspace(0.0, 360.0, 24, endpoint=False)
    theta_grid, phi_grid = np.meshgrid(theta, phi, indexing="ij")
    power_db = -np.abs(theta_grid - 20.0)
    path = plots.pattern_cut_png(theta_grid, phi_grid, power_db, cut_phi_deg=0.0, path=tmp_path / "cut.png")
    _assert_real_png(path)


def test_vortex_phase_ring_png(tmp_path):
    alpha_deg = np.linspace(0.0, 360.0, 72, endpoint=False)
    phase_deg = np.degrees(np.radians(alpha_deg))  # a winding=+1 sawtooth
    path = plots.vortex_phase_ring_png(alpha_deg, phase_deg, winding_number=1.0, path=tmp_path / "vortex.png")
    _assert_real_png(path)


def test_patch_layout_png_renders_both_plain_and_annular_slot_cells(tmp_path):
    cells_mm = [
        (-2.0, 0.0, 2.0, 2.0),  # plain rectangle patch
        (2.0, 0.0, 2.0, 2.0, 0.7, 0.4),  # PHASE15 Part E: patch with an annular slot
    ]
    path = plots.patch_layout_png(cells_mm, period_mm=4.0, nx=2, ny=1, path=tmp_path / "mixed_layout.png")
    _assert_real_png(path)
