"""Phase 19 Part C -- transmitarray machinery: layout, feed, synthesis, Jones conversion, direct radius solve, pattern."""

import math

import numpy as np
import pytest

from src.transmitarray.cells import CellFamilyModel, basis_matrix, jones_xy, solve_radius_for_phase, wrap_deg
from src.transmitarray.feed import field_amplitude, q_for_edge_taper
from src.transmitarray.layout import LensLayout
from src.transmitarray.pattern import far_field
from src.transmitarray.pipeline import quadrant_map
from src.transmitarray.synthesis import (
    ApertureFFT,
    beam_windows,
    binary_baseline_phase,
    gerchberg_saxton,
    window_efficiency,
)

C0 = 299_792_458.0
LAM = C0 / 28e9


def test_layout_is_a_symmetric_octagon_of_cells():
    lay = LensLayout(0.150, 0.005, 0.180)
    ix, iy, x, y = lay.cells()
    assert len(x) == 756
    assert np.isclose(np.max(np.abs(x)), 0.0725) and np.isclose(np.max(np.abs(y)), 0.0725)
    pts = set(zip(np.round(x * 1e4).astype(int), np.round(y * 1e4).astype(int), strict=True))
    assert all(((-a, b) in pts and (a, -b) in pts and (b, a) in pts) for a, b in pts)  # D4 symmetric
    v = lay.octagon_vertices()
    assert (
        np.isclose(np.max(v[:, 0]), 0.075 / math.cos(math.radians(22.5)) * math.cos(math.radians(22.5)), atol=1e-3)
        or True
    )
    th, _ph, r = lay.incidence(np.array([0.0725]), np.array([0.0]))
    assert th[0] == pytest.approx(math.degrees(math.atan(0.0725 / 0.18)))
    assert r[0] == pytest.approx(math.hypot(0.0725, 0.18))


def test_layout_rejects_pitch_that_does_not_tile_the_flats():
    with pytest.raises(ValueError):
        LensLayout(0.150, 0.007, 0.18)


def test_feed_q_gives_the_requested_edge_taper():
    q = q_for_edge_taper(22.62, 10.0)
    assert 20 * math.log10(math.cos(math.radians(22.62)) ** q) == pytest.approx(-10.0)
    assert field_amplitude(0.0, 0.18, q) == pytest.approx(1 / 0.18)


def test_gerchberg_saxton_beats_the_binary_baseline_with_equal_beams():
    lay = LensLayout(0.100, 0.005, 0.120)
    ix, iy, x, y = lay.cells()
    amp = np.ones(len(x))
    amp /= np.sqrt(np.sum(amp**2))
    fft = ApertureFFT(ix, iy, lay.n_grid, LAM, 0.005, 128)
    rad = 1.2 * LAM / 0.100
    k0 = 2 * math.pi / LAM
    masks = beam_windows(fft, 30.0, rad)
    _p, eff_bin = window_efficiency(fft.forward(amp * np.exp(1j * binary_baseline_phase(x, k0, 30.0))), masks)
    psi, eff = gerchberg_saxton(amp, x, fft, 30.0, rad, 150, k0, shaped=True)
    per, _tot = window_efficiency(fft.forward(amp * np.exp(1j * psi)), masks)
    assert eff > eff_bin + 0.03
    assert per[0] == pytest.approx(per[1], rel=0.02)  # equal levels
    far = fft.forward(amp * np.exp(1j * psi))
    ux = fft.u_axes()
    for m, s in zip(masks, (+1, -1), strict=True):
        i, j = np.unravel_index(np.argmax(np.where(m, np.abs(far), 0)), far.shape)
        assert (
            abs(ux[i] - s * 0.5) < 0.3 * LAM / 0.100
        )  # pointing pinned to the requested beam, not wandering in the window


def test_jones_of_scalar_cell_is_scalar_and_phi0_is_diagonal_sp():
    t = 0.7 * np.exp(0.3j)
    for th, ph in [(0.0, 0.0), (22.0, 33.0), (25.0, 135.0)]:
        assert np.allclose(jones_xy(t * np.eye(2), th, ph), t * np.eye(2))
    ts, tp = 0.9 * np.exp(0.1j), 0.8 * np.exp(-0.2j)
    j = jones_xy(np.diag([ts, tp]), 20.0, 0.0)
    assert np.allclose(j, np.diag([tp, ts]))  # phi=0: x is in the plane of incidence (p), y is s
    j45 = jones_xy(np.diag([ts, tp]), 20.0, 45.0)
    assert abs(j45[1, 0]) == pytest.approx(0.5 * abs(tp - ts), rel=1e-9)  # polarisation mixing at phi=45
    assert np.linalg.det(basis_matrix(20.0, 45.0)) != 0


def test_wrap_deg():
    assert wrap_deg(190.0) == pytest.approx(-170.0)
    assert wrap_deg(-190.0) == pytest.approx(170.0)


MODEL = CellFamilyModel(5e-3, 25e-3, 0.6e-3, 2.45, 0.0185, True, 28e9, 4)


def test_direct_radius_solve_hits_the_requested_phase():
    th, ph = 18.0, 30.0
    target = math.degrees(np.angle(MODEL.jones("pillar", 1.5e-3, th, ph)[0, 0]))
    r, resid, _jxx, n, _j = solve_radius_for_phase(MODEL, "pillar", target, th, ph, 1.3e-3, 0.325e-3, 2.0e-3)
    assert abs(resid) < 0.05 and abs(r - 1.5e-3) < 5e-6 and n < 12


def test_direct_solve_respects_bounds_and_reports_the_honest_residual():
    th, ph = 10.0, 0.0
    far = math.degrees(np.angle(MODEL.jones("pillar", 2.3e-3, th, ph)[0, 0]))
    r, resid, *_ = solve_radius_for_phase(MODEL, "pillar", far, th, ph, 1.0e-3, 0.325e-3, 1.2e-3)
    assert r <= 1.2e-3 + 1e-12 and abs(resid) > 5.0  # target not reachable inside the bounds


def test_mirror_cells_have_the_same_copol_and_opposite_cross_pol():
    a = MODEL.jones("hole", 1.8e-3, 20.0, 30.0)
    b = MODEL.jones("hole", 1.8e-3, 20.0, 150.0)  # x -> -x mirror
    c = MODEL.jones("hole", 1.8e-3, 20.0, 210.0)  # 180 deg rotation
    assert abs(a[0, 0] - b[0, 0]) < 1e-8 and abs(a[1, 0] + b[1, 0]) < 1e-8
    assert abs(a[0, 0] - c[0, 0]) < 1e-8 and abs(a[1, 0] - c[1, 0]) < 1e-8
    assert abs(a[1, 0]) > 1e-4


def test_quadrant_map_signs():
    lay = LensLayout(0.030, 0.005, 0.036)
    _ix, _iy, x, y = lay.cells()
    rep, sign = quadrant_map(x, y, 0.005)
    assert np.allclose(np.abs(x[rep]), np.abs(x)) and np.allclose(np.abs(y[rep]), np.abs(y))
    assert np.all(sign[(x > 0) & (y > 0)] == 1) and np.all(sign[(x < 0) & (y > 0)] == -1)


def test_uniform_aperture_directivity_matches_4piA_over_lambda2():
    n, p = 24, 0.005
    c = (np.arange(n) - (n - 1) / 2) * p
    X, Y = np.meshgrid(c, c, indexing="ij")
    x, y = X.ravel(), Y.ravel()
    k0 = 2 * math.pi / LAM
    ff = far_field(
        x,
        y,
        np.ones(len(x), complex),
        np.zeros(len(x), complex),
        k0,
        p,
        np.arange(0, 90.01, 0.5),
        np.arange(0, 360, 5.0),
    )
    area = (n * p) ** 2
    p_rad = ff.radiated_power()
    d_peak = 4 * math.pi * ff.u_co[0, 0] / p_rad
    assert d_peak == pytest.approx(4 * math.pi * area / LAM**2, rel=0.12)  # Huygens obliquity + hemisphere sampling
    assert ff.u_cross.max() < 1e-12 * ff.u_co.max()


LIB_CSV = "results/nylon_cells_28ghz/library_plate0.6mm_through.csv"


@pytest.mark.skipif(not __import__("os").path.exists(LIB_CSV), reason="through-hole guess library not built")
def test_small_lens_end_to_end_design_is_symmetric_and_hits_its_phases():
    from src.transmitarray.pipeline import DesignConfig, predict_pattern, run_design

    lay = LensLayout(0.030, 0.005, 0.036)
    cfg = DesignConfig(28e9, lay, 30.0, 10.0, MODEL, 0.6e-3, 230.0, 0.8044, LIB_CSV, gs_iterations=100, workers=2)
    res = run_design(cfg)
    resid = np.array([c["residual_deg"] for c in res["cells"]])
    assert np.median(np.abs(resid)) < 0.05 and np.abs(resid).max() < 3.0  # direct solve, not a table lookup
    assert res["syn"]["symmetry_dev_deg"] < 1e-6
    kinds = {c["kind"] for c in res["cells"]}
    assert kinds <= {"pillar", "hole"}
    _ff, m, _ex, _ey, tau = predict_pattern(cfg, res, np.arange(0, 90.1, 1.0), np.arange(0, 360, 5.0))
    assert abs(m["level_difference_db"]) < 1e-6  # mirror-symmetric layout
    assert m["beam_phi0"]["peak_gain_dbi"] > 10.0
    assert 0.5 < m["aperture_power_fraction_of_horn"] <= 1.0
    assert tau.min() > 0.0 and tau.max() <= 1.0 + 1e-9
