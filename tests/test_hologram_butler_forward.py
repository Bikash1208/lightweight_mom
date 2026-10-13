"""Phase 25 A: the arbitrary-excitation forward model (`src/hologram/butler_forward.py`) -- linearity in the
excitation, the single-feed degenerate case equals the existing vector far field, the designed excitation of a
real (small) hex4 design is reproduced, footprints and amplitude models behave."""

import numpy as np
import pytest

from src.green_function import LayeredGreenFunction
from src.hologram.butler_forward import (
    FeedFarField,
    board_from_cells,
    evaluate_excitation,
    far_field,
    monopole_incident_field,
)
from src.materials import Substrate
from src.pattern.vector_aperture_pattern import _vector_far_field, aperture_current_from_scalar

F = 10.54e9
SUB = Substrate("user", 2.2, 0.0009, 1.575e-3)
K_SWG = 225.0
P = 6.4e-3


def _toy_board(n=14, feeds=((0.0, 0.0), (0.02, 0.01), (-0.015, 0.03)), seed=3, **kw):
    xs = (np.arange(n) - (n - 1) / 2) * P
    xx, yy = np.meshgrid(xs, xs, indexing="xy")
    rng = np.random.default_rng(seed)
    g0, g1 = np.exp(1j * np.radians(138.1)), np.exp(1j * np.radians(-41.4))
    gamma = np.where(rng.integers(0, 2, xx.shape) == 1, g1, g0)
    return board_from_cells(xx, yy, gamma, feeds, P, F, SUB, K_SWG, **kw)


THETA = np.array([0.0, 7.0, 20.0, 33.0, 48.0, 61.0, 75.0, 20.0, 33.0])
PHI = np.array([0.0, 40.0, 80.0, 130.0, 200.0, 260.0, 300.0, 10.0, 355.0])


def test_pattern_is_linear_in_the_excitation():
    board = _toy_board()
    ff = FeedFarField(board)
    rng = np.random.default_rng(0)
    a = rng.normal(size=3) + 1j * rng.normal(size=3)
    b = rng.normal(size=3) + 1j * rng.normal(size=3)
    ea, pa = ff.fields(a, THETA, PHI)
    eb, pb = ff.fields(b, THETA, PHI)
    es, ps = ff.fields(a + 2.0j * b, THETA, PHI)
    assert np.max(np.abs(es - (ea + 2.0j * eb))) <= 1e-12 * np.max(np.abs(es))
    assert np.max(np.abs(ps - (pa + 2.0j * pb))) <= 1e-12 * np.max(np.abs(ps))
    # ... and the one-transform convenience (current summed first) agrees with the per-feed sum
    ec, pc = far_field(board, a, THETA, PHI)
    assert np.max(np.abs(ec - ea)) <= 1e-12 * np.max(np.abs(ea))
    assert np.max(np.abs(pc - pa)) <= 1e-12 * np.max(np.abs(pa))


def test_single_feed_is_the_existing_vector_far_field():
    board = _toy_board()
    ff = FeedFarField(board)
    m = 1
    green = LayeredGreenFunction(SUB, F, backing="pec")
    s = board.e_inc[m] * board.gamma
    jx, jy = aperture_current_from_scalar(board.x_m, board.y_m, s, feed_xy_m=tuple(board.feeds_xy_m[m]))
    et, ep = _vector_far_field(THETA, PHI, board.x_m, board.y_m, jx, jy, F, green, P, P)
    a = np.zeros(3)
    a[m] = 1.0
    et2, ep2 = ff.fields(a, THETA, PHI)
    assert np.max(np.abs(et2 - et)) <= 1e-12 * np.max(np.abs(et))
    assert np.max(np.abs(ep2 - ep)) <= 1e-12 * np.max(np.abs(ep))


def test_global_phase_of_the_excitation_does_not_change_the_pattern_magnitude():
    board = _toy_board()
    ff = FeedFarField(board)
    a = np.array([1.0, 0.5j, -0.3])
    e1, p1 = ff.fields(a, THETA, PHI)
    e2, p2 = ff.fields(a * np.exp(0.7j), THETA, PHI)
    assert np.allclose(np.abs(e1), np.abs(e2), rtol=1e-12)
    assert np.allclose(np.abs(p1), np.abs(p2), rtol=1e-12)


def test_footprint_zeroes_cells_a_feed_does_not_light():
    n = 8
    xs = (np.arange(n) - (n - 1) / 2) * P
    xx, yy = np.meshgrid(xs, xs, indexing="xy")
    gamma = np.ones(xx.shape, dtype=complex)
    feeds = ((-0.02, 0.0), (0.02, 0.0))
    fp = np.stack([xx.ravel() < 0, xx.ravel() >= 0])
    board = board_from_cells(xx, yy, gamma, feeds, P, F, SUB, K_SWG, footprints=fp)
    assert np.all(board.e_inc[0][~fp[0]] == 0) and np.all(board.e_inc[1][~fp[1]] == 0)
    s = board.aperture_values([1.0, 1.0])
    assert np.allclose(np.abs(s), 1.0)  # every cell lit by exactly one feed, phase-only


def test_amplitude_models():
    x = np.array([0.04, 0.16])
    y = np.zeros(2)
    ph = monopole_incident_field(x, y, (0.0, 0.0), K_SWG, "phase_only")
    assert np.allclose(np.abs(ph), 1.0)
    assert np.allclose(np.angle(ph), np.angle(np.exp(-1j * K_SWG * x)))
    phys = monopole_incident_field(x, y, (0.0, 0.0), K_SWG, "physical", period_m=P)
    assert np.isclose(np.abs(phys[0]) / np.abs(phys[1]), 2.0, rtol=1e-12)  # rho^-1/2: 4x the distance, half the field
    assert np.allclose(np.angle(phys), np.angle(ph))
    with pytest.raises(ValueError):
        monopole_incident_field(x, y, (0.0, 0.0), K_SWG, "bogus")


def test_bad_inputs_are_errors():
    board = _toy_board()
    ff = FeedFarField(board)
    with pytest.raises(ValueError):
        ff.fields([1.0, 0.0], THETA, PHI)
    with pytest.raises(ValueError):
        board.aperture_values([1.0])
    with pytest.raises(ValueError):
        board_from_cells([0.0, 1.0], [0.0, 0.0], [1.0], [(0.0, 0.0)], P, F, SUB, K_SWG)


def test_evaluate_excitation_runs_and_reports_the_beam():
    from src.hologram.sll_board import HeldOutGrid

    board = _toy_board(n=20, feeds=((0.0, 0.0),))
    ff = FeedFarField(board)
    coarse = HeldOutGrid(
        theta_step_deg=5.0, phi_step_deg=10.0, cut_step_u=0.02, patch_radius_deg=3.0, patch_step_deg=1.0
    )
    ev = evaluate_excitation(board, ff, [1.0], [(30.0, 45.0, "rhcp")], grid=coarse)
    assert len(ev["beams"]) == 1 and np.isfinite(ev["peak_directivity_dbi"])
    assert {"directivity_dbi", "pointing_error_deg", "sll_global_db", "xpd_db"} <= set(ev["beams"][0])


# ----- a real (small) hex4 design: the designed excitations are reproduced


@pytest.fixture(scope="module")
def hex_board():
    pytest.importorskip("fastapi")
    from src.app.routers.hmsa import _run_synthesize_ports_job
    from src.hologram.butler_forward import board_from_ports_result
    from src.hologram.sll_board import build_problem
    from tests.test_app_phase23_hex4 import _Job, _request

    res = _run_synthesize_ports_job(_Job("p25fwd"), _request({"preset": "hex4", "side_mm": 50.0, "gap_mm": 0.0}))
    beams = [[(30.0, 45.0 + 90.0 * k)] for k in range(4)]
    prob = build_problem(res, F, 2.2, 0.0009, 1.575, 6.4, beams, du_far=0.05, du_near=0.03)
    stored = board_from_ports_result(res, 6.4, F, 2.2, 0.0009, 1.575, incident="stored")
    analytic = board_from_ports_result(res, 6.4, F, 2.2, 0.0009, 1.575, incident="analytic")
    return res, prob, stored, analytic


def test_designed_excitation_reproduces_the_stored_aperture_field(hex_board):
    _res, _prob, stored, _analytic = hex_board
    s = stored.meta["stored_aperture_field"]
    assert np.max(np.abs(stored.aperture_values(np.ones(4)) - s)) < 1e-12


def test_designed_excitation_reproduces_the_existing_far_field(hex_board):
    from src.hologram.sll_board import _FarField

    _res, prob, stored, analytic = hex_board
    levels = {p.port_id: p.state for p in prob.opt_ports}
    alpha = {p.port_id: 0.0 for p in prob.opt_ports}
    th = np.linspace(0.0, 80.0, 25)
    ph = np.linspace(0.0, 350.0, 25)
    for active, a in (({1}, np.eye(4)[1]), ({0, 1, 2, 3}, np.ones(4))):
        et, ep = _FarField(prob, levels, alpha, F, active)(th, ph)
        scale = max(np.abs(et).max(), np.abs(ep).max())
        et2, ep2 = FeedFarField(stored).fields(a, th, ph)
        assert max(np.abs(et2 - et).max(), np.abs(ep2 - ep).max()) <= 1e-9 * scale
        # the analytic incident field (phase-only TM0 at the cell centres) is close but not identical
        et3, ep3 = FeedFarField(analytic).fields(a, th, ph)
        assert max(np.abs(et3 - et).max(), np.abs(ep3 - ep).max()) <= 0.05 * scale


def test_a_feed_alone_leaves_the_other_regions_dark(hex_board):
    _res, _prob, _stored, analytic = hex_board
    own = analytic.meta["owner_port"]
    for p in range(4):
        s = analytic.aperture_values(np.eye(4)[p])
        assert np.all(s[own != p] == 0) and np.all(np.abs(s[own == p]) > 0.9)
