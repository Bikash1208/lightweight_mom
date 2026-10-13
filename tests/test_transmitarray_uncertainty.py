"""Phase 19 Part D-bis -- the Monte Carlo error-bar machinery (operator, neighbours, metrics) on a small lens."""

import math

import numpy as np
import pytest

from src.transmitarray import pattern as pt
from src.transmitarray.layout import LensLayout
from src.transmitarray.uncertainty import PatternOperator, neighbour_lists, percentiles, phase_steps_deg

C0 = 299_792_458.0
K0 = 2 * math.pi * 28e9 / C0


def small_lens():
    lay = LensLayout(0.040, 0.005, 0.048)
    _ix, _iy, x, y = lay.cells()
    th, _ph, r = lay.incidence(x, y)
    # a synthetic dual-beam-like field: unit amplitude, phase steering to +-30 deg (cos -> 0/pi), feed phase removed
    ph = np.where(np.cos(K0 * math.sin(math.radians(30.0)) * x) >= 0, 0.0, math.pi)
    ex = np.exp(1j * ph) * np.cos(np.radians(th)) ** 4
    return lay, x, y, th, r, ex


def test_neighbour_lists_have_axes_and_borders():
    x = [-2.5, 2.5, -2.5, 2.5]
    y = [-2.5, -2.5, 2.5, 2.5]
    nb = neighbour_lists(x, y, 5.0)
    assert all(len(n) == 2 for n in nb)  # a 2x2 block: every cell has one x- and one y-neighbour
    for n in nb:
        assert sorted(ax for _j, ax in n) == ["x", "y"]
    assert dict((j, ax) for j, ax in nb[0]) == {1: "x", 2: "y"}


def test_phase_steps_are_wrapped_absolute():
    t = np.array([np.exp(1j * np.radians(10)), np.exp(1j * np.radians(350))])
    steps = phase_steps_deg(t, [[(1, "x")], [(0, "x")]])
    assert steps[0][0] == pytest.approx(20.0) and steps[1][0] == pytest.approx(20.0)


def test_operator_matches_the_direct_far_field_at_the_beam_and_with_zero_perturbation():
    lay, x, y, th, r, ex = small_lens()
    uinc = (x / r, y / r)
    op = PatternOperator(x, y, K0, lay.pitch_m, uinc, 30.0, theta_step=1.0, phi_step=10.0)
    res = op.run(ex[:, None])
    ff = pt.far_field(x, y, ex, 0 * ex, K0, lay.pitch_m, np.arange(20, 40.01, 0.1), np.array([0.0]), u_inc=uinc)
    assert res["gain_dbi"][0, 0] == pytest.approx(ff.gain_dbi("co").max(), abs=0.05)
    # draws with eps = 0 are identical copies, and a mirror-symmetric field gives equal beams (the identity)
    two = op.run(np.repeat(ex[:, None], 3, axis=1))
    assert np.allclose(two["gain_dbi"], two["gain_dbi"][0])
    assert abs(two["level_diff_db"][0]) < 1e-6


def test_perturbation_breaks_the_level_identity_and_moves_the_beam():
    lay, x, y, th, r, ex = small_lens()
    op = PatternOperator(x, y, K0, lay.pitch_m, (x / r, y / r), 30.0, theta_step=1.0, phi_step=10.0)
    rng = np.random.default_rng(3)
    eps = 0.8 * np.exp(1j * rng.uniform(0, 2 * np.pi, (len(x), 50)))
    res = op.run(ex[:, None] * (1 + eps))
    assert (
        np.std(res["level_diff_db"]) > 0.05
    )  # the 0.00 dB of the nominal layout is an identity, not a robust property
    assert res["sll_db"].mean() > op.run(ex[:, None])["sll_db"][0]  # random errors raise the sidelobes
    assert np.all(res["pointing_error_deg"] >= 0)


def test_percentiles_keys():
    p = percentiles(np.arange(101.0))
    assert p == {"p10": 10.0, "p50": 50.0, "p90": 90.0}
