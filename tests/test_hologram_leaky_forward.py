"""Phase 26 A: the leakage-aware forward wrapper -- ray geometry, power accounting, and that a unit
amplitude profile leaves the existing pattern path bit-identical."""

import numpy as np
import pytest

from src.hologram import leaky_forward as lf
from src.hologram.butler_forward import board_from_cells, far_field
from src.materials import Substrate

HEXAGON = np.array([[np.cos(a), np.sin(a)] for a in np.radians(np.arange(0, 360, 60) + 30.0)]) * 0.15  # circumradius


def test_ray_distance_regular_hexagon():
    apothem = 0.15 * np.cos(np.radians(30.0))
    # vertices at 30, 90, ...: edge midpoints (apothem) at 0, 60, ...; corners (circumradius) at 30, 90
    d = lf.polygon_ray_distance(HEXAGON, (0.0, 0.0), np.radians([0.0, 30.0, 60.0, 90.0, 45.0]))
    assert d[0] == pytest.approx(apothem, rel=1e-9)
    assert d[1] == pytest.approx(0.15, rel=1e-9)
    assert d[2] == pytest.approx(apothem, rel=1e-9)
    assert d[3] == pytest.approx(0.15, rel=1e-9)
    assert d[4] == pytest.approx(apothem / np.cos(np.radians(15.0)), rel=1e-9)


def test_ray_distance_off_centre_square():
    sq = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    d = lf.polygon_ray_distance(sq, (0.25, 0.5), np.array([0.0, np.pi, np.pi / 2]))
    assert d == pytest.approx([0.75, 0.25, 0.5])


def test_converted_fraction_is_one_minus_exp_of_twice_alpha_rho():
    x = np.array([0.01, 0.02, -0.03])
    y = np.array([0.0, 0.01, 0.02])
    alpha = 3.0
    p = lf.leakage_profile(x, y, (0.0, 0.0), HEXAGON, lambda psi: np.full_like(psi, alpha), 6.4e-3, n_psi=360)
    d = lf.polygon_ray_distance(HEXAGON, (0.0, 0.0), p.psi_grid_rad)
    assert p.eta_conv == pytest.approx(np.mean(1.0 - np.exp(-2.0 * alpha * d)), rel=1e-12)
    assert p.eta_conv + p.edge_fraction == pytest.approx(1.0)
    # amplitude: rho^-1/2 exp(-alpha rho), rho floored at half a period
    rho = np.hypot(x, y)
    assert p.amplitude == pytest.approx(np.maximum(rho, 3.2e-3) ** -0.5 * np.exp(-alpha * rho))


def test_no_leakage_converts_nothing():
    p = lf.leakage_profile(np.array([0.05]), np.array([0.0]), (0.0, 0.0), HEXAGON, lambda psi: 0.0 * psi, 6.4e-3)
    assert p.eta_conv == 0.0 and p.edge_fraction == 1.0


def test_unit_amplitude_profile_is_bit_identical():
    rng = np.random.default_rng(3)
    x = rng.uniform(-0.05, 0.05, 40)
    y = rng.uniform(-0.05, 0.05, 40)
    gamma = np.exp(1j * rng.uniform(-np.pi, np.pi, 40))
    sub = Substrate("s", 2.2, 0.0009, 1.575e-3)
    board = board_from_cells(x, y, gamma, [(0.0, 0.0)], 6.4e-3, 10.54e9, sub, 1.02 * 220.6, "rhcp", "phase_only")
    scaled = lf.scaled_board(board, np.ones(40))
    th = np.array([0.0, 20.0, 45.0])
    ph = np.array([0.0, 90.0, 225.0])
    a = far_field(board, [1.0], th, ph)
    b = far_field(scaled, [1.0], th, ph)
    assert np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1])
    assert board.e_inc.shape == scaled.e_inc.shape


def test_pattern_budget_rhcp_plus_lhcp_is_total():
    rng = np.random.default_rng(5)
    theta = np.arange(0.0, 90.0 + 1e-9, 5.0)
    phi = np.arange(0.0, 360.0, 10.0)
    et = rng.normal(size=(theta.size, phi.size)) + 1j * rng.normal(size=(theta.size, phi.size))
    ep = rng.normal(size=(theta.size, phi.size)) + 1j * rng.normal(size=(theta.size, phi.size))
    b = lf.pattern_budget(theta, phi, et, ep, 30.0, 225.0)
    assert b["rhcp_fraction"] + b["lhcp_fraction"] == pytest.approx(1.0, rel=1e-9)
    assert 0.0 <= b["within_8deg_fraction"] <= b["within_20deg_fraction"] <= 1.0
