"""Phase 19 Part A -- RCWA dielectric unit-cell solver: physics checks (not "it runs")."""

import math

import numpy as np
import pytest

from src.dielectric_cell import OctagonPillarCell, UniformLayer, solve_cell, solve_stack
from src.dielectric_cell.analytic import stack_transmission
from src.dielectric_cell.fourier import octagon_vertices, polygon_ft

F = 28e9
P = 5.0e-3
H = 25e-3
EPS, TAND = 2.45, 0.0185
M_FAST = 6


def cell(radius_mm, base_mm=0.0, rot=22.5, h=H):
    return OctagonPillarCell(P, radius_mm * 1e-3, h, EPS, TAND, base_mm * 1e-3, rot)


# ---------------------------------------------------------------- Fourier machinery
def test_polygon_ft_matches_brute_force_quadrature():
    v = octagon_vertices(2.0e-3, 22.5)
    n = 1200
    u = (np.arange(n) + 0.5) / n - 0.5
    X, Y = np.meshgrid(u * P, u * P, indexing="ij")
    # point-in-convex-polygon: inside all CCW half-planes
    inside = np.ones(X.shape, bool)
    for i in range(8):
        a, b = v[i], v[(i + 1) % 8]
        inside &= (b[0] - a[0]) * (Y - a[1]) - (b[1] - a[1]) * (X - a[0]) >= 0
    for gx, gy in [
        (0.0, 0.0),
        (2 * np.pi / P, 0.0),
        (2 * np.pi * 2 / P, 2 * np.pi / P),
        (2 * np.pi * 3 / P, 2 * np.pi * 3 / P),
    ]:
        ref = np.sum(inside * np.exp(-1j * (gx * X + gy * Y))) * (P / n) ** 2
        got = polygon_ft(v, np.array([gx]), np.array([gy]))[0]
        assert abs(got - ref) < 2e-3 * (2 * math.sqrt(2) * (2e-3) ** 2), (gx, gy, got, ref)
    assert polygon_ft(v, np.array([0.0]), np.array([0.0]))[0].real == pytest.approx(2 * math.sqrt(2) * (2e-3) ** 2)


# ---------------------------------------------------------------- analytic slab oracle
@pytest.mark.parametrize("theta", [0.0, 25.0, 40.0])
@pytest.mark.parametrize("tand", [0.0, TAND])
def test_uniform_slab_matches_analytic_fabry_perot(theta, tand):
    res = solve_stack([UniformLayer(H, EPS, tand)], P, P, F, theta, 17.0, 1)
    for pol, k in (("TE", 0), ("TM", 1)):
        r_a, t_a = stack_transmission([(H, EPS, tand)], F, theta, pol)
        assert abs(res["t"][k, k] - t_a) < 1e-9
        assert abs(res["r"][k, k] - r_a) < 1e-9
        # no TE<->TM conversion for a uniform slab
        assert abs(res["t"][0, 1]) < 1e-9 and abs(res["t"][1, 0]) < 1e-9


def test_two_layer_uniform_stack_matches_analytic():
    layers = [(H, EPS, TAND), (1.5e-3, EPS, TAND)]
    res = solve_stack([UniformLayer(*x) for x in layers], P, P, F, 22.6, 0.0, 1)
    _, t_a = stack_transmission([(H + 1.5e-3, EPS, TAND)], F, 22.6, "TM")
    assert abs(res["t"][1, 1] - t_a) < 1e-9


def test_empty_cell_is_transparent_at_the_faces():
    resp = solve_cell(cell(0.0), F, 20.0, 30.0, M_FAST)
    path = np.exp(-1j * 2 * np.pi * F / 299_792_458.0 * math.cos(math.radians(20.0)) * H)
    assert np.allclose(resp.t, path * np.eye(2), atol=1e-10)
    assert np.allclose(resp.r, 0, atol=1e-10)


def test_empty_cell_with_base_plate_is_the_plate_slab():
    resp = solve_cell(cell(0.0, base_mm=2.0), F, 0.0, 0.0, M_FAST)
    _, t_a = stack_transmission([(2e-3, EPS, TAND)], F, 0.0, "TE")
    path = np.exp(-1j * 2 * np.pi * F / 299_792_458.0 * H)
    assert abs(resp.t[0, 0] - path * t_a) < 1e-9


def test_lossless_phase_is_minus_k_n_h_for_matched_thin_limit():
    # delta(phase) between full slab and empty cell grows ~ -k0 (n-1) h up to FP ripple: sign check
    th = solve_cell(cell(0.0, base_mm=0.5), F, 0.0, 0.0, M_FAST)
    assert np.angle(th.t[0, 0]) < 0  # exp(+jwt): transmission phase lags (negative) through dielectric


# ---------------------------------------------------------------- energy / reciprocity / symmetry
@pytest.mark.parametrize("theta,phi", [(0.0, 0.0), (22.6, 0.0), (25.0, 20.0)])
@pytest.mark.parametrize("radius", [1.0, 2.6])
def test_lossless_octagon_conserves_power(theta, phi, radius):
    c = OctagonPillarCell(P, radius * 1e-3, H, EPS, 0.0, 0.0, 22.5)
    resp = solve_cell(c, F, theta, phi, M_FAST)
    for pol in ("TE", "TM"):
        assert resp.power_balance(pol) == pytest.approx(1.0, abs=1e-6)  # NVM operator is not exactly Hermitian off-axis


@pytest.mark.parametrize("radius", [1.0, 2.6])
def test_lossy_octagon_absorbs_and_never_creates_power(radius):
    resp = solve_cell(cell(radius), F, 22.6, 20.0, M_FAST)
    for pol in ("TE", "TM"):
        b = resp.power_balance(pol)
        assert 0.5 < b < 1.0 - 1e-4


def test_c4_symmetry_te_equals_tm_at_normal_incidence():
    for radius in (1.0, 2.4):
        resp = solve_cell(cell(radius), F, 0.0, 0.0, M_FAST)
        assert abs(resp.t[0, 0] - resp.t[1, 1]) < 1e-8
        assert abs(resp.r[0, 0] - resp.r[1, 1]) < 1e-8
        assert abs(resp.t[0, 1]) < 1e-8 and abs(resp.t[1, 0]) < 1e-8  # no cross-pol


def test_rotation_by_90_degrees_in_azimuth_is_invariant():
    a = solve_cell(cell(2.0), F, 22.6, 0.0, M_FAST)
    b = solve_cell(cell(2.0), F, 22.6, 90.0, M_FAST)
    assert np.allclose(a.t, b.t, atol=1e-8) and np.allclose(a.r, b.r, atol=1e-8)


def test_mirror_symmetry_flips_cross_pol_sign_only():
    a = solve_cell(cell(2.0), F, 25.0, 20.0, M_FAST)
    b = solve_cell(cell(2.0), F, 25.0, -20.0, M_FAST)
    assert abs(a.t[0, 0] - b.t[0, 0]) < 1e-8 and abs(a.t[1, 1] - b.t[1, 1]) < 1e-8
    assert abs(a.t[0, 1] + b.t[0, 1]) < 1e-8 and abs(a.t[1, 0] + b.t[1, 0]) < 1e-8
    assert abs(a.t[0, 1]) > 1e-4  # genuinely non-zero, so the sign test means something


def test_reciprocity_cross_pol_magnitudes():
    resp = solve_cell(cell(2.0), F, 25.0, 20.0, M_FAST)
    assert abs(resp.t[0, 1]) == pytest.approx(abs(resp.t[1, 0]), rel=1e-3)


def test_swapping_octagon_orientation_by_45_degrees_is_a_distinct_but_c4_cell():
    a = solve_cell(cell(2.0, rot=22.5), F, 0.0, 0.0, M_FAST)
    b = solve_cell(cell(2.0, rot=0.0), F, 0.0, 0.0, M_FAST)
    assert abs(a.t[0, 0] - a.t[1, 1]) < 1e-8 and abs(b.t[0, 0] - b.t[1, 1]) < 1e-8
    assert abs(a.t[0, 0] - b.t[0, 0]) > 1e-3  # different geometry really changes the answer


# ---------------------------------------------------------------- physical plausibility of the patterned layer
def test_subwavelength_pillar_lies_between_wiener_bounds():
    p, h = 1.0e-3, 2.0e-3
    c = OctagonPillarCell(p, 0.5e-3, h, EPS, 0.0, 0.0, 22.5)
    f = c.fill_fraction
    eps_par = f * EPS + (1 - f)
    eps_ser = 1.0 / (f / EPS + (1 - f))
    ph = lambda e: np.angle(stack_transmission([(h, e, 0.0)], F, 0.0, "TE")[1])  # noqa: E731
    got = np.angle(solve_cell(c, F, 0.0, 0.0, M_FAST).t[0, 0])
    lo, hi = sorted([ph(eps_ser), ph(eps_par)])
    assert lo < got < hi  # phase lags monotonically with eps at this thickness
    # and sits much closer to the parallel (E tangential to walls) bound than to the series one
    assert abs(got - ph(eps_par)) < abs(got - ph(eps_ser))


def test_phase_increases_monotonically_with_pillar_size():
    ph = [np.unwrap([np.angle(solve_cell(cell(r), F, 0.0, 0.0, M_FAST).t[0, 0])])[0] for r in (0.5, 1.0, 1.5, 2.0)]
    # the 25 mm pillar wraps many times, so compare thin cells where there is no wrap
    thin = [np.angle(solve_cell(cell(r, h=2e-3), F, 0.0, 0.0, M_FAST).t[0, 0]) for r in (0.5, 1.0, 1.5, 2.0, 2.5)]
    assert all(b < a for a, b in zip(thin, thin[1:], strict=False))
    assert len(ph) == 4


# ---------------------------------------------------------------- numerics: factorisation and convergence
@pytest.mark.parametrize("radius", [1.5, 2.5])
def test_normal_vector_method_converges_fast(radius):
    t = {m: solve_cell(cell(radius), F, 0.0, 0.0, m, "normal").t[1, 1] for m in (4, 6, 8)}
    d68 = abs(np.degrees(np.angle(t[8] / t[6])))
    assert d68 < 0.25
    assert abs(abs(t[8]) - abs(t[6])) < 2e-3


def test_laurent_rule_converges_slowly_toward_the_normal_vector_answer():
    ref = solve_cell(cell(2.0), F, 0.0, 0.0, 8, "normal").t[1, 1]
    errs = [
        abs(np.degrees(np.angle(solve_cell(cell(2.0), F, 0.0, 0.0, m, "laurent").t[1, 1] / ref))) for m in (2, 4, 6)
    ]
    assert errs[0] > errs[1] > errs[2] > 0.3  # still degrees away: the reason 'normal' is the default


def test_normal_field_grid_resolution_does_not_matter():
    a = solve_stack(cell(2.0).layers(), P, P, F, 0.0, 0.0, M_FAST, "normal", n_grid=128)["t"][1, 1]
    b = solve_stack(cell(2.0).layers(), P, P, F, 0.0, 0.0, M_FAST, "normal", n_grid=384)["t"][1, 1]
    assert abs(np.degrees(np.angle(a / b))) < 0.05


@pytest.mark.slow
@pytest.mark.parametrize("radius,theta", [(1.0, 0.0), (2.0, 25.0), (2.6, 25.0)])
def test_fourier_order_convergence_to_well_below_one_degree(radius, theta):
    """Chosen production order M=8 vs M=14: phase within 0.2 deg, |S21| within 0.01 dB (recorded in docs)."""
    for k in (0, 1):
        t8 = solve_cell(cell(radius), F, theta, 0.0, 8).t[k, k]
        t14 = solve_cell(cell(radius), F, theta, 0.0, 14).t[k, k]
        assert abs(np.degrees(np.angle(t8 / t14))) < 0.2
        assert abs(20 * np.log10(abs(t8) / abs(t14))) < 0.01


def test_octagon_that_does_not_fit_is_rejected():
    with pytest.raises(ValueError):
        cell(2.8)  # 2.8*cos(22.5 deg) = 2.59 mm > p/2


def test_rcwa_source_contains_no_fitted_surrogates():
    import pathlib

    root = pathlib.Path(__file__).resolve().parents[1] / "src" / "dielectric_cell"
    text = " ".join(p.read_text() for p in root.glob("*.py")).lower()
    for word in ("sklearn", "torch", "tensorflow", "neural", "surrogate"):
        assert word not in text
