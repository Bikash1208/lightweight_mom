"""Analytic limits, TE/TM degeneracy, symmetry/reciprocity and branch checks of the dyad."""

import math

import numpy as np
import pytest

from src.constants import EPS0, ETA0, MU0
from src.green_function import LayeredGreenFunction, te_tm_unit_vectors
from src.materials import Substrate
from src.periodic_bc import RectangularLattice, normal_wavenumber

F = 10e9
H = 1.575e-3


def gf(eps_r=2.2, tand=0.0009, h=H, f=F):
    return LayeredGreenFunction(Substrate("s", eps_r, tand, h), f)


def textbook_impedances(g, kt):
    """Z_TE, Z_TM from the textbook tan(kz h) TL input impedance (independent formula path)."""
    w = g.omega
    eps = g.substrate.eps_complex
    a, b = g.normal_wavenumbers(kt, 0.0)
    t = np.tan(b * g.substrate.thickness)

    def down(Zc_air, Zc_sub):
        return Zc_sub * (Zc_air + 1j * Zc_sub * t) / (Zc_sub + 1j * Zc_air * t)

    def par(Zu, Zd):
        return Zu * Zd / (Zu + Zd)

    Zte0, Zte1 = w * MU0 / a, w * MU0 / b
    Ztm0, Ztm1 = a / (w * EPS0), b / (w * EPS0 * eps)
    return par(Zte0, down(Zte0, Zte1)), par(Ztm0, down(Ztm0, Ztm1))


REGIMES = [0.05, 0.3, 0.7, 0.95, 1.05, 1.3, 1.45, 2.5, 8.0]   # kt/k0 (k_sub/k0 = 1.483 lossless)


# ============================================================ lossless substrate

@pytest.mark.parametrize("kt_r", REGIMES)
def test_lossless_matches_textbook_tl_formula(kt_r):
    g = gf(tand=0.0)
    kt = kt_r * g.k0
    z_te, z_tm = g.equivalent_impedances(kt, 0.0)
    t_te, t_tm = textbook_impedances(g, kt)
    assert z_te == pytest.approx(t_te, rel=1e-12)
    assert z_tm == pytest.approx(t_tm, rel=1e-12)


@pytest.mark.parametrize("h", [1e-4, 1.575e-3, 0.02])
def test_air_filled_stack_is_free_space_for_any_thickness(h):
    """eps_r = 1: substrate is invisible; Z_TE = w*mu0/(2 kz), Z_TM = kz/(2 w eps0)."""
    g = gf(eps_r=1.0, tand=0.0, h=h)
    for kt_r in (0.0, 0.4, 0.9, 1.7, 6.0):
        kt = kt_r * g.k0
        kz = normal_wavenumber(g.k0, kt, 0.0)
        z_te, z_tm = g.equivalent_impedances(kt, 0.0)
        assert z_te == pytest.approx(g.omega * MU0 / (2 * kz), rel=1e-12)
        assert z_tm == pytest.approx(kz / (2 * g.omega * EPS0), rel=1e-12)


def test_free_space_sheet_radiates_minus_half_eta0_j():
    """Textbook result: a uniform current sheet in free space gives E = -(eta0/2) J."""
    g = gf(eps_r=1.0, tand=0.0)
    G = g.spectral_dyad(0.0, 0.0)
    assert G == pytest.approx(-0.5 * ETA0 * np.eye(2), rel=1e-12)
    assert G.imag == pytest.approx(0.0, abs=1e-9)


def test_zero_thickness_limit_is_free_space():
    g = gf(eps_r=4.0, tand=0.0, h=1e-12)
    ref = gf(eps_r=1.0, tand=0.0)
    for kt_r in (0.0, 0.5, 2.0):
        assert g.spectral_dyad(kt_r * g.k0, 0.0) == pytest.approx(
            ref.spectral_dyad(kt_r * g.k0, 0.0), rel=1e-9
        )


def test_lossless_real_part_only_from_radiating_harmonics():
    """Lossless: Re Z = 0 for kt > k0 (reactive only); Re Z > 0 for kt < k0 (radiation)."""
    g = gf(tand=0.0)
    for kt_r in (1.05, 1.3, 1.45, 2.0, 5.0, 40.0):
        z_te, z_tm = g.equivalent_impedances(kt_r * g.k0, 0.0)
        assert abs(z_te.real) <= 1e-12 * abs(z_te)
        assert abs(z_tm.real) <= 1e-12 * abs(z_tm)
    for kt_r in (0.0, 0.2, 0.6, 0.95):
        z_te, z_tm = g.equivalent_impedances(kt_r * g.k0, 0.0)
        assert z_te.real > 0.0 and z_tm.real > 0.0


def test_evanescent_reactance_signs_follow_exp_plus_jwt():
    """Below cutoff (free space) TE is inductive (Im Z > 0) and TM capacitive (Im Z < 0)."""
    g = gf(eps_r=1.0, tand=0.0)
    for kt_r in (1.2, 2.0, 10.0):
        z_te, z_tm = g.equivalent_impedances(kt_r * g.k0, 0.0)
        assert z_te.imag > 0.0 and z_tm.imag < 0.0


def test_large_thickness_lossless_kernel_is_finite_away_from_poles():
    g = gf(tand=0.0, h=0.05)
    z_te, z_tm = g.equivalent_impedances(0.3 * g.k0, 0.0)
    assert np.isfinite(z_te) and np.isfinite(z_tm)


# =============================================================== lossy substrate

@pytest.mark.parametrize("kt_r", REGIMES)
@pytest.mark.parametrize("tand", [0.0009, 0.02, 0.3])
def test_lossy_matches_textbook_tl_formula(kt_r, tand):
    g = gf(tand=tand)
    kt = kt_r * g.k0
    z_te, z_tm = g.equivalent_impedances(kt, 0.0)
    t_te, t_tm = textbook_impedances(g, kt)
    assert z_te == pytest.approx(t_te, rel=1e-11)
    assert z_tm == pytest.approx(t_tm, rel=1e-11)


def test_lossy_is_passive_for_every_kt():
    """-J^H G J = Z_TM|Jp|^2 + Z_TE|Jq|^2 must have Re >= 0 for a passive stack."""
    g = gf(tand=0.0009)
    kt = np.linspace(0.0, 30 * g.k0, 6001)
    z_te, z_tm = g.equivalent_impedances(kt, 0.0)
    assert np.all(z_te.real > 0.0)
    assert np.all(z_tm.real >= 0.0)      # TM is exactly 0 at kt = k0 (kz_air = 0)
    assert np.all(z_tm.real[np.abs(kt - g.k0) > 1e-6 * g.k0] > 0.0)


def test_lossy_dyad_is_positive_dissipation_operator_for_random_currents():
    g = gf(tand=0.0009)
    rng = np.random.default_rng(0)
    for _ in range(200):
        kx, ky = rng.uniform(-6 * g.k0, 6 * g.k0, 2)
        J = rng.normal(size=2) + 1j * rng.normal(size=2)
        power = -0.5 * np.real(np.conj(J) @ (g.spectral_dyad(kx, ky) @ J))
        assert power > 0.0


def test_lossy_reduces_to_lossless_as_tand_goes_to_zero():
    g0 = gf(tand=0.0)
    for kt_r in (0.3, 1.3, 3.0):      # avoid the lossless guided-mode poles
        kt = kt_r * g0.k0
        for tand in (1e-6, 1e-9, 1e-12):
            g1 = gf(tand=tand)
            assert g1.spectral_dyad(kt, 0.3 * kt) == pytest.approx(g0.spectral_dyad(kt, 0.3 * kt), rel=10 * tand + 1e-12)


def test_thick_lossy_substrate_becomes_parallel_of_two_half_spaces():
    """h -> infinity: reflections from the far side vanish, Z = Z_air || Z_substrate."""
    g = gf(tand=0.02, h=5.0)
    eps = g.substrate.eps_complex
    for kt_r in (0.0, 0.5, 1.2, 3.0):
        kt = kt_r * g.k0
        a, b = g.normal_wavenumbers(kt, 0.0)
        z_te, z_tm = g.equivalent_impedances(kt, 0.0)
        assert z_te == pytest.approx(g.omega * MU0 / (a + b), rel=1e-10)
        A, B = a, b / eps
        assert z_tm == pytest.approx((A * B / (A + B)) / (g.omega * EPS0), rel=1e-10)


def test_lossy_center_value_against_normal_incidence_transmission_line():
    """kt = 0 by hand: Z = eta0 || Zin, Zin = eta1 (eta0 + j eta1 tan(k1 h)) / (eta1 + j eta0 tan(k1 h))."""
    g = gf(tand=0.0009)
    n = g.substrate.refractive_index
    eta1 = ETA0 / n
    t = np.tan(g.k_substrate * H)
    z_in = eta1 * (ETA0 + 1j * eta1 * t) / (eta1 + 1j * ETA0 * t)
    z = ETA0 * z_in / (ETA0 + z_in)
    G = g.spectral_dyad(0.0, 0.0)
    assert G[0, 0] == pytest.approx(-z, rel=1e-10)
    assert G[1, 1] == pytest.approx(-z, rel=1e-10)


# ==================================================== kt = 0 and TE/TM degeneracy

@pytest.mark.parametrize("eps_r, tand, h", [(2.2, 0.0, H), (2.2, 0.0009, H), (10.2, 0.05, 3e-3), (1.0, 0.0, H)])
def test_te_tm_degenerate_at_normal_incidence(eps_r, tand, h):
    g = gf(eps_r, tand, h)
    z_te, z_tm = g.equivalent_impedances(0.0, 0.0)
    assert z_te == pytest.approx(z_tm, rel=1e-13)
    G = g.spectral_dyad(0.0, 0.0)
    assert G[0, 1] == 0.0 and G[1, 0] == 0.0
    assert G[0, 0] == G[1, 1]                    # exactly isotropic (identity * scalar)
    assert G[0, 0] == pytest.approx(-0.5 * (z_te + z_tm), rel=1e-13)


def test_te_tm_degeneracy_is_only_at_kt_zero():
    g = gf()
    z_te, z_tm = g.equivalent_impedances(0.2 * g.k0, 0.0)
    assert abs(z_te - z_tm) > 1e-2 * abs(z_te)
    # and the splitting closes quadratically in kt
    d1 = abs(np.subtract(*g.equivalent_impedances(1e-3 * g.k0, 0.0)))
    d2 = abs(np.subtract(*g.equivalent_impedances(1e-2 * g.k0, 0.0)))
    assert d2 / d1 == pytest.approx(100.0, rel=1e-3)


@pytest.mark.parametrize("angle", [0.0, 0.3, 1.1, math.pi / 2, 2.5, -1.0])
def test_kt_to_zero_limit_is_direction_independent_and_matches_kt_zero(angle):
    g = gf()
    G0 = g.spectral_dyad(0.0, 0.0)
    errs = []
    for kt in (1e-2 * g.k0, 1e-3 * g.k0, 1e-4 * g.k0):
        Gk = g.spectral_dyad(kt * math.cos(angle), kt * math.sin(angle))
        errs.append(np.linalg.norm(Gk - G0) / np.linalg.norm(G0))
    assert errs[0] < 1e-3 and errs[2] < 1e-7
    assert errs[1] == pytest.approx(errs[0] / 100, rel=0.05)     # O(kt^2) approach


def test_kt_zero_threshold_is_seamless():
    """No jump between the exact-zero branch and the general formula just above it."""
    g = gf()
    tiny = 1e-11 * g.k0                     # just above the isotropic-branch threshold
    G_general = g.spectral_dyad(tiny * 0.6, tiny * 0.8)
    G_zero = g.spectral_dyad(0.0, 0.0)
    assert G_general == pytest.approx(G_zero, rel=1e-12)


def test_unit_vectors_orthonormal_and_kt_zero_fallback():
    kx = np.array([0.0, 3.0, -2.0, 0.0])
    ky = np.array([0.0, 4.0, 5.0, -7.0])
    p, q = te_tm_unit_vectors(kx, ky)
    assert np.allclose(np.einsum("...i,...i->...", p, p), 1.0)
    assert np.allclose(np.einsum("...i,...i->...", q, q), 1.0)
    assert np.allclose(np.einsum("...i,...i->...", p, q), 0.0, atol=1e-15)
    # q = z x p  ->  (-py, px)
    assert np.allclose(q[..., 0], -p[..., 1]) and np.allclose(q[..., 1], p[..., 0])
    assert p[1] == pytest.approx([0.6, 0.8]) and q[1] == pytest.approx([-0.8, 0.6])
    # kt = 0 (with a floor) -> arbitrary orthonormal pair x, y
    p0, q0 = te_tm_unit_vectors(0.0, 0.0, kt_floor=1e-9)
    assert p0.tolist() == [1.0, 0.0] and q0.tolist() == [0.0, 1.0]


def test_all_outputs_finite_on_normal_incidence_floquet_grid():
    g = gf()
    lat = RectangularLattice(8e-3, 8e-3)
    kx, ky, G = g.floquet_dyads(lat, 0.0, 0.0, 12)
    assert np.all(np.isfinite(G))
    assert G[12, 12] == pytest.approx(g.spectral_dyad(0.0, 0.0))


# ============================================================ reciprocity / symmetry

def random_k(g, n=40, seed=3, scale=8.0):
    rng = np.random.default_rng(seed)
    return rng.uniform(-scale * g.k0, scale * g.k0, n), rng.uniform(-scale * g.k0, scale * g.k0, n)


@pytest.mark.parametrize("tand", [0.0, 0.0009, 0.05])
def test_dyad_is_complex_symmetric(tand):
    """Reciprocity of the layered medium: G_xy = G_yx (symmetric, not Hermitian)."""
    g = gf(tand=tand)
    kx, ky = random_k(g)
    G = g.spectral_dyad(kx, ky)
    assert np.array_equal(G[..., 0, 1], G[..., 1, 0])
    assert np.allclose(G, np.swapaxes(G, -1, -2), rtol=0, atol=0)
    assert np.any(np.abs(G[..., 0, 1]) > 0)              # not trivially diagonal
    if tand > 0:
        assert not np.allclose(G, np.conj(np.swapaxes(G, -1, -2)))   # genuinely non-Hermitian


@pytest.mark.parametrize("tand", [0.0, 0.0009])
def test_dyad_even_in_k(tand):
    """G(-k) = G(k): with symmetry this is Lorentz reciprocity G(k)^T = G(-k)."""
    g = gf(tand=tand)
    kx, ky = random_k(g)
    assert g.spectral_dyad(-kx, -ky) == pytest.approx(g.spectral_dyad(kx, ky), rel=1e-13)


def test_mirror_symmetry():
    g = gf()
    kx, ky = random_k(g)
    G, Gm = g.spectral_dyad(kx, ky), g.spectral_dyad(-kx, ky)     # x -> -x mirror
    assert Gm[..., 0, 0] == pytest.approx(G[..., 0, 0], rel=1e-13)
    assert Gm[..., 1, 1] == pytest.approx(G[..., 1, 1], rel=1e-13)
    assert Gm[..., 0, 1] == pytest.approx(-G[..., 0, 1], rel=1e-13)


def test_axis_exchange_symmetry():
    g = gf()
    kx, ky = random_k(g)
    G, Gs = g.spectral_dyad(kx, ky), g.spectral_dyad(ky, kx)
    assert Gs[..., 0, 0] == pytest.approx(G[..., 1, 1], rel=1e-13)
    assert Gs[..., 1, 1] == pytest.approx(G[..., 0, 0], rel=1e-13)
    assert Gs[..., 0, 1] == pytest.approx(G[..., 0, 1], rel=1e-13)


@pytest.mark.parametrize("theta", [0.3, 1.0, math.pi / 2, 2.6])
def test_rotational_covariance(theta):
    """Isotropic stack: G(R k) = R G(k) R^T for every in-plane rotation R."""
    g = gf()
    kx, ky = random_k(g)
    c, s = math.cos(theta), math.sin(theta)
    R = np.array([[c, -s], [s, c]])
    G = g.spectral_dyad(kx, ky)
    Gr = g.spectral_dyad(c * kx - s * ky, s * kx + c * ky)
    expected = np.einsum("ij,njk,lk->nil", R, G, R)
    assert Gr == pytest.approx(expected, rel=1e-11, abs=1e-9)


def test_floquet_grid_symmetries_at_normal_incidence():
    g = gf()
    N = 7
    lat = RectangularLattice(8e-3, 8e-3)
    _, _, G = g.floquet_dyads(lat, 0.0, 0.0, N)
    assert np.allclose(G, G[::-1, ::-1])                                   # (m,n) -> (-m,-n)
    assert np.allclose(G[..., 0, 0], np.swapaxes(G[..., 1, 1], 0, 1))      # square lattice m<->n
    assert np.allclose(G[..., 0, 1], np.swapaxes(G[..., 0, 1], 0, 1))
    assert np.allclose(G[..., 0, 1], -G[::-1, :, 0, 1])                    # m -> -m flips G_xy
    assert np.allclose(G[..., 0, 1], -G[:, ::-1, 0, 1])                    # n -> -n flips G_xy


def test_floquet_dyads_indexing_and_shapes():
    g = gf()
    lat = RectangularLattice(8e-3, 6e-3)
    kx0, ky0 = 40.0, -25.0
    kx, ky, G = g.floquet_dyads(lat, kx0, ky0, 3)
    assert kx.shape == ky.shape == (7, 7) and G.shape == (7, 7, 2, 2)
    assert kx[3, 3] == kx0 and ky[3, 3] == ky0
    assert kx[5, 1] == pytest.approx(kx0 + 2 * lat.reciprocal_x)
    assert ky[5, 1] == pytest.approx(ky0 - 2 * lat.reciprocal_y)
    assert G[5, 1] == pytest.approx(g.spectral_dyad(kx[5, 1], ky[5, 1]))


# ================================================================== kz branch

@pytest.mark.parametrize("tand", [0.0, 0.0009, 0.3])
def test_both_kz_on_decaying_branch_and_exponential_bounded(tand):
    g = gf(tand=tand)
    kt = np.linspace(0, 60 * g.k0, 4001)
    a, b = g.normal_wavenumbers(kt, 0.0)
    assert np.all(a.imag <= 0.0) and np.all(b.imag <= 0.0)
    assert np.all(a.real >= 0.0) and np.all(b.real >= 0.0)
    assert np.all(np.abs(np.exp(-2j * b * H)) <= 1.0 + 1e-15)


def test_air_cutoff_at_kt_equals_k0_is_finite_and_continuous():
    """kz_air = 0 exactly: Z_TM -> 0, Z_TE finite (no division by kz_air anywhere).

    kz_air = sqrt(k0^2 - kt^2) has a square-root branch point at kt = k0 (Rayleigh
    anomaly), so Z is continuous and *linear in kz_air*, not in kt: the deviation
    from the cutoff value must scale like |kz_air|.
    """
    g = gf()
    z_te, z_tm = g.equivalent_impedances(g.k0, 0.0)
    assert np.isfinite(z_te) and z_tm == 0.0
    assert np.all(np.isfinite(g.spectral_dyad(g.k0, 0.0)))
    slopes = []
    for eps in (1e-8, 1e-10, 1e-12, -1e-8, -1e-10):
        kt = g.k0 * (1 + eps)
        a = abs(normal_wavenumber(g.k0, kt, 0.0))
        z2_te, z2_tm = g.equivalent_impedances(kt, 0.0)
        slopes.append((abs(z2_te - z_te) / a, abs(z2_tm) / a))
    slopes = np.array(slopes)
    assert np.allclose(slopes[:, 0], slopes[0, 0], rtol=1e-2)     # TE deviation ~ |kz_air|
    assert np.allclose(slopes[:, 1], slopes[0, 1], rtol=1e-2)     # TM -> 0 ~ |kz_air|
    assert slopes[0, 0] * g.k0 < 1e2 * abs(z_te)                  # bounded slope


def test_substrate_cutoff_kz1_zero_is_finite_and_continuous():
    """Lossless kt = k_sub exactly (kz_sub = 0): removable 0/0 of the naive exponential form."""
    g = gf(tand=0.0)
    kt = g.k_substrate.real
    z_te, z_tm = g.equivalent_impedances(kt, 0.0)
    assert np.isfinite(z_te) and np.isfinite(z_tm)
    for eps in (1e-9, -1e-9):
        z2_te, z2_tm = g.equivalent_impedances(kt * (1 + eps), 0.0)
        assert z2_te == pytest.approx(z_te, rel=1e-6)
        assert z2_tm == pytest.approx(z_tm, rel=1e-6)
    # closed form of the limit: Z_TE = w mu0 (1 + j a h) / (a (2 + j a h))
    a = normal_wavenumber(g.k0, kt, 0.0)
    assert z_te == pytest.approx(g.omega * MU0 * (1 + 1j * a * H) / (a * (2 + 1j * a * H)), rel=1e-12)


def test_no_overflow_at_huge_kt_and_asymptotics():
    """Naive cos/sin or exp(+|Im|) forms overflow; this one tends to the half-space result."""
    g = gf()
    eps = g.substrate.eps_complex
    for kt in (1e4, 1e6, 1e8):
        a, b = g.normal_wavenumbers(kt, 0.0)
        z_te, z_tm = g.equivalent_impedances(kt, 0.0)
        assert np.isfinite(z_te) and np.isfinite(z_tm)
        assert z_te == pytest.approx(g.omega * MU0 / (a + b), rel=1e-9)
        A, B = a, b / eps
        assert z_tm == pytest.approx((A * B / (A + B)) / (g.omega * EPS0), rel=1e-9)
    # large-kt scaling: Z_TE ~ j w mu0 / (2 kt), Z_TM ~ -j kt / ((eps+1) w eps0)
    kt = 1e8
    z_te, z_tm = g.equivalent_impedances(kt, 0.0)
    assert z_te == pytest.approx(1j * g.omega * MU0 / (2 * kt), rel=1e-6)
    assert z_tm == pytest.approx(-1j * kt / ((eps + 1) * g.omega * EPS0), rel=1e-6)


def test_frequency_and_substrate_validation():
    sub = Substrate("s", 2.2, 0.0009, H)
    for bad in (0.0, -1e9, float("nan"), float("inf"), "10e9"):
        with pytest.raises(ValueError):
            LayeredGreenFunction(sub, bad)
    g = LayeredGreenFunction(sub, F)
    assert g.k0 == pytest.approx(2 * math.pi * F / 299792458.0, rel=1e-15)
    assert g.k_substrate.imag < 0.0
