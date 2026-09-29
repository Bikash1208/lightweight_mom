"""Tangential E/H at arbitrary z: sheet-radiated fields and unpatterned-stack fields."""

import math

import numpy as np
import pytest
from reference_layered import solve_current_sheet, solve_incident_fields

from src.constants import ETA0
from src.green_function import LayeredGreenFunction
from src.materials import Substrate

F = 10e9
H = 1.575e-3
STACKS = [(2.2, 0.0, H), (2.2, 0.0009, H), (10.2, 0.05, 3e-3)]
KT_OVER_K0 = [0.0, 0.3, 0.9, 1.3, 2.5]
ANGLES = [0.0, 0.7, 2.0]
CURRENTS = [(1.0, 0.0), (0.0, 1.0), (0.4 + 0.2j, -0.3j)]


def make(eps_r, tand, h):
    sub = Substrate("s", eps_r, tand, h)
    return sub, LayeredGreenFunction(sub, F)


def z_samples(h):
    """Points in all three regions, at both interfaces and just either side of them."""
    return [3e-2, 1e-3, 1e-12, 0.0, -0.3 * h, -0.9 * h, -h, -h - 1e-12, -h - 2e-3, -h - 4e-2]


def sz(e, h):
    """Time-average z-directed Poynting flux (1/2) Re(E x H*)_z."""
    return 0.5 * np.real(e[..., 0] * np.conj(h[..., 1]) - e[..., 1] * np.conj(h[..., 0]))


# =============================================================== sheet fields


@pytest.mark.parametrize("eps_r, tand, h", STACKS)
def test_sheet_fields_match_independent_field_matching_everywhere(eps_r, tand, h):
    sub, g = make(eps_r, tand, h)
    for kt_r in KT_OVER_K0:
        for ang in ANGLES:
            kx, ky = kt_r * g.k0 * math.cos(ang), kt_r * g.k0 * math.sin(ang)
            for j in CURRENTS:
                ref = solve_current_sheet(kx, ky, j, F, sub.eps_complex, h)["fields"]
                for z in z_samples(h):
                    for side in ("above", "below"):
                        e, hh = g.sheet_fields(kx, ky, np.array(j), z, side)
                        e_ref, h_ref = ref(z, side)
                        assert e == pytest.approx(e_ref, rel=1e-9, abs=1e-9 * max(np.abs(e_ref).max(), 1e-30))
                        assert hh == pytest.approx(h_ref, rel=1e-9, abs=1e-9 * max(np.abs(h_ref).max(), 1e-30))


@pytest.mark.parametrize("eps_r, tand, h", STACKS)
def test_sheet_E_at_z0_is_the_spectral_dyad_and_boundary_conditions_hold(eps_r, tand, h):
    """E_t(0) = G J (Phase 2); E_t continuous at z=0 and z=-h; H_t continuous at z=-h;
    z x [H(0+) - H(0-)] = J on the sheet."""
    _, g = make(eps_r, tand, h)
    rng = np.random.default_rng(5)
    for _ in range(25):
        kx, ky = rng.uniform(-3 * g.k0, 3 * g.k0, 2)
        j = rng.normal(size=2) + 1j * rng.normal(size=2)
        e_above, h_above = g.sheet_fields(kx, ky, j, 0.0, "above")
        e_below, h_below = g.sheet_fields(kx, ky, j, 0.0, "below")
        assert e_above == pytest.approx(g.spectral_dyad(kx, ky) @ j, rel=1e-11)
        assert e_above == pytest.approx(e_below, rel=1e-12)
        dh = h_above - h_below
        assert np.array([-dh[1], dh[0]]) == pytest.approx(j, rel=1e-10)  # z x dH = J
        e_in, h_in = g.sheet_fields(kx, ky, j, -h, "above")
        e_out, h_out = g.sheet_fields(kx, ky, j, -h - 1e-12, "above")
        assert e_in == pytest.approx(e_out, rel=1e-9) and h_in == pytest.approx(h_out, rel=1e-9)
        # tiny offsets around the sheet: continuous E, H jumps by J
        e_p, h_p = g.sheet_fields(kx, ky, j, 1e-13)
        e_m, h_m = g.sheet_fields(kx, ky, j, -1e-13)
        assert e_p == pytest.approx(e_m, rel=1e-8)
        assert h_p - h_m == pytest.approx(np.array([j[1], -j[0]]), rel=1e-6)


def test_normal_incidence_fields_are_isotropic_and_exactly_transverse_rotation_invariant():
    """kt = 0: E = -Z J and H = c z x J with the *same* scalars for both components."""
    _, g = make(2.2, 0.0009, H)
    j = np.array([0.7 - 0.2j, -1.1 + 0.5j])
    for z in (1e-3, 0.0, -0.4 * H, -H, -5e-3):
        e, hh = g.sheet_fields(0.0, 0.0, j, z)
        assert e[0] / j[0] == pytest.approx(e[1] / j[1], rel=1e-14)
        zxj = np.array([-j[1], j[0]])
        assert hh[0] / zxj[0] == pytest.approx(hh[1] / zxj[1], rel=1e-14)
    # far above: outgoing plane wave, H = z x E / eta0
    e, hh = g.sheet_fields(0.0, 0.0, j, 5e-2)
    assert hh == pytest.approx(np.array([-e[1], e[0]]) / ETA0, rel=1e-12)
    e, hh = g.sheet_fields(0.0, 0.0, j, -5e-2)  # far below: H = -z x E / eta0
    assert hh == pytest.approx(-np.array([-e[1], e[0]]) / ETA0, rel=1e-12)


def test_sheet_power_balance_through_the_fields_lossless():
    """Source power -1/2 Re(J* . E) = outgoing flux above + outgoing flux below (from E and H
    at z = +/- far); constant flux with z in each air region."""
    _, g = make(2.2, 0.0, H)
    for kt_r in (0.0, 0.4, 0.9):
        kx, ky = kt_r * g.k0 * 0.6, kt_r * g.k0 * 0.8
        j = np.array([0.8 + 0.1j, -0.3 + 0.6j])
        p_src = -0.5 * np.real(np.conj(j) @ g.sheet_fields(kx, ky, j, 0.0)[0])
        flux = lambda z: sz(*g.sheet_fields(kx, ky, j, z))
        up = flux(1e-3)
        down = -flux(-H - 1e-3)
        assert up > 0.0 and down > 0.0
        assert flux(0.05) == pytest.approx(up, rel=1e-10)  # no loss, no evanescent power
        assert flux(-H - 0.05) == pytest.approx(-down, rel=1e-10)
        assert up + down == pytest.approx(p_src, rel=1e-10)
        # inside the lossless substrate the net upward flux is constant too
        assert flux(-0.5 * H) == pytest.approx(flux(-0.2 * H), rel=1e-9)


def test_sheet_power_balance_lossy_includes_substrate_dissipation():
    sub, g = make(2.2, 0.05, H)
    kx, ky = 0.0, 0.0  # normal incidence: E_z = 0, so E_t carries all the dissipation
    j = np.array([1.0, 0.0])
    p_src = -0.5 * np.real(np.conj(j) @ g.sheet_fields(kx, ky, j, 0.0)[0])
    flux = lambda z: sz(*g.sheet_fields(kx, ky, j, z))
    up, down = flux(1e-3), -flux(-H - 1e-3)
    assert p_src > up + down > 0.0
    # dissipation (omega eps0/2)(-Im eps) integral |E|^2 dz  by quadrature of the fields
    z = np.linspace(-H, 0.0, 2001)
    e, _ = g.sheet_fields(kx * np.ones_like(z), ky * np.ones_like(z), j, z, "below")
    p_diss = (
        0.5 * g.omega * 8.8541878128e-12 * (-sub.eps_complex.imag) * np.trapezoid(np.sum(np.abs(e) ** 2, axis=-1), z)
    )
    assert up + down + p_diss == pytest.approx(p_src, rel=2e-6)


def test_evanescent_fields_decay_and_nothing_overflows():
    _, g = make(2.2, 0.0009, H)
    kt = 3.0 * g.k0
    j = np.array([1.0, 0.0])
    z = np.array([1e-4, 1e-3, 1e-2, 3e-2])  # exp(-kt z) stays representable
    e, hh = g.sheet_fields(kt, 0.0, j, z)
    mag = np.linalg.norm(e, axis=-1)
    assert np.all(np.diff(mag) < 0.0) and np.all(np.isfinite(hh))
    e, _ = g.sheet_fields(kt, 0.0, j, -H - z)
    assert np.all(np.diff(np.linalg.norm(e, axis=-1)) < 0.0)
    far, far_h = g.sheet_fields(kt, 0.0, j, np.array([1.0, 50.0]))  # underflows to 0, never inf/nan
    assert np.all(np.isfinite(far)) and np.all(np.isfinite(far_h)) and np.all(np.abs(far) < 1e-30)
    # extreme spectral content and distance: still finite (all exponents are <= 0)
    for kt_big in (1e4, 1e6):
        e, hh = g.sheet_fields(kt_big, 0.3 * kt_big, j, np.array([1e-9, -0.5 * H, -H - 1e-6, 10.0]))
        assert np.all(np.isfinite(e)) and np.all(np.isfinite(hh))


def test_rayleigh_point_fields_are_finite_and_linear_in_kz_air():
    """kz_air = 0 (kt = k0): E and H are finite; they are linear in kz_air (square-root branch point in kt,
    e.g. TE H_t = Y_air V with Y_air ~ kz_air -> 0), so deviations from the point scale like |kz_air|."""
    _, g = make(2.2, 0.0, H)
    j = np.array([0.6, -0.4j])
    for z in (1e-3, -0.5 * H, -H - 1e-3):
        e0, h0 = g.sheet_fields(g.k0, 0.0, j, z)
        assert np.all(np.isfinite(e0)) and np.all(np.isfinite(h0))
        slopes = []
        for eps in (1e-8, 1e-10, 1e-12, -1e-8, -1e-10):
            kt = g.k0 * (1 + eps)
            a = abs(g.normal_wavenumbers(kt, 0.0)[0])
            e1, h1 = g.sheet_fields(kt, 0.0, j, z)
            slopes.append((np.linalg.norm(e1 - e0) / a, np.linalg.norm(h1 - h0) / a))
        slopes = np.array(slopes)
        assert np.allclose(slopes[:, 0], slopes[0, 0], rtol=2e-2)
        assert np.allclose(slopes[:, 1], slopes[0, 1], rtol=2e-2)


def test_substrate_cutoff_point_is_finite_and_continuous():
    """kz_sub = 0 exactly (lossless, kt = Re k_sub): the naive 0/0 is removable and the fields are
    continuous in kt through it (E and H, all regions)."""
    _, g = make(2.2, 0.0, H)
    kt0 = g.k_substrate.real
    assert g.normal_wavenumbers(kt0, 0.0)[1] == 0.0
    j = np.array([0.6, -0.4j])
    for z in (1e-3, -0.5 * H, -H - 1e-3):
        e0, h0 = g.sheet_fields(kt0, 0.0, j, z)
        assert np.all(np.isfinite(e0)) and np.all(np.isfinite(h0))
        for eps in (1e-9, -1e-9):
            e1, h1 = g.sheet_fields(kt0 * (1 + eps), 0.0, j, z)
            assert e1 == pytest.approx(e0, rel=1e-5) and h1 == pytest.approx(h0, rel=1e-5)


def test_shapes_broadcasting_and_argument_validation():
    _, g = make(2.2, 0.0009, H)
    kx = np.linspace(0, g.k0, 4)[:, None]
    z = np.array([1e-3, -5e-4, -3e-3])[None, :]
    e, hh = g.sheet_fields(kx, 0.0, np.array([1.0, 0.0]), z)
    assert e.shape == hh.shape == (4, 3, 2)
    # per-point current arrays
    j = np.ones((4, 3, 2))
    assert g.sheet_fields(kx, 0.0, j, z)[0].shape == (4, 3, 2)
    with pytest.raises(ValueError, match="sheet_side"):
        g.sheet_fields(0.0, 0.0, np.array([1.0, 0.0]), 0.0, "left")


# ============================================================== stack fields


@pytest.mark.parametrize("eps_r, tand, h", STACKS)
@pytest.mark.parametrize("port", ["top", "bottom"])
def test_stack_fields_match_independent_field_matching_everywhere(eps_r, tand, h, port):
    sub, g = make(eps_r, tand, h)
    for th in (0.0, 0.4, 1.1):
        for ph in (0.0, 0.9, 2.5):
            k = g.k0 * math.sin(th)
            kx, ky = k * math.cos(ph), k * math.sin(ph)
            c, s, cp, sp = math.cos(th), math.sin(th), math.cos(ph), math.sin(ph)
            for tm in (False, True):
                e3 = np.array([c * cp, c * sp, s if port == "top" else -s]) if tm else np.array([-sp, cp, 0.0])
                ref = solve_incident_fields(kx, ky, e3, F, sub.eps_complex, h, port)
                for z in z_samples(h):
                    e, hh = g.stack_fields(kx, ky, e3[:2], z, port)
                    e_ref, h_ref = ref(z)
                    assert e == pytest.approx(e_ref, rel=1e-9, abs=1e-9 * max(np.abs(e_ref).max(), 1e-30))
                    assert hh == pytest.approx(h_ref, rel=1e-9, abs=1e-9 * max(np.abs(h_ref).max(), 1e-30))


@pytest.mark.parametrize("port", ["top", "bottom"])
def test_stack_fields_lossless_have_constant_net_power_flow_and_unit_amplitude(port):
    """Lossless unpatterned stack: the net z-directed flux is the same at every z (no source, no loss);
    magnitude = (1 - |Gamma|^2)/(2 eta0), directed toward the stack from the incident side."""
    _, g = make(2.2, 0.0, H)
    e_inc = np.array([1.0, 0.0])
    z = np.array([0.05, 1e-3, 0.0, -0.3 * H, -H, -H - 1e-3, -0.05])
    e, hh = g.stack_fields(0.0, 0.0, e_inc, z, port)
    flux = sz(e, hh)
    assert np.allclose(flux, flux[0], rtol=1e-10)
    sign = -1.0 if port == "top" else +1.0
    assert sign * flux[0] > 0.0
    # incident amplitude is 1 at the port's reference interface: |E_t| there is |1 + Gamma| <= 2
    z_ref = 0.0 if port == "top" else -H
    e_ref = g.stack_fields(0.0, 0.0, e_inc, z_ref, port)[0]
    r = (1 - math.sqrt(2.2)) / (1 + math.sqrt(2.2))
    u = np.exp(-2j * g.k0 * math.sqrt(2.2) * H)
    gamma = r * (1 - u) / (1 - r * r * u)
    assert e_ref[0] == pytest.approx(1 + gamma, rel=1e-12)
    assert abs(flux[0]) == pytest.approx((1 - abs(gamma) ** 2) / (2 * ETA0), rel=1e-10)


def test_bottom_port_is_the_mirror_image_of_the_top_port():
    """The unpatterned stack is symmetric under z -> -h - z: E is mirrored, H changes sign."""
    _, g = make(10.2, 0.05, 3e-3)
    h = 3e-3
    e_inc = np.array([0.3, 0.9])
    for z in (0.02, 0.0, -0.4e-3, -h, -h - 5e-3):
        et, ht = g.stack_fields(0.0, 0.0, e_inc, z, "top")
        eb, hb = g.stack_fields(0.0, 0.0, e_inc, -h - z, "bottom")
        assert eb == pytest.approx(et, rel=1e-12) and hb == pytest.approx(-ht, rel=1e-12)


def test_stack_fields_are_continuous_across_both_interfaces_and_time_convention():
    _, g = make(2.2, 0.0009, H)
    e_inc = np.array([1.0, 0.0])
    for port in ("top", "bottom"):
        for z0 in (0.0, -H):
            lo = g.stack_fields(0.0, 0.0, e_inc, z0 - 1e-12, port)
            hi = g.stack_fields(0.0, 0.0, e_inc, z0 + 1e-12, port)
            assert lo[0] == pytest.approx(hi[0], rel=1e-8) and lo[1] == pytest.approx(hi[1], rel=1e-8)
    # exp(+jwt): the incident wave (travelling -z) has phase exp(+j k0 z): E(z=1mm) = E(0) exp(+j k0 1mm)
    # for a substrate-free stack (eps = 1) where nothing is reflected
    _, g1 = make(1.0, 0.0, H)
    e0 = g1.stack_fields(0.0, 0.0, e_inc, 0.0, "top")[0]
    e1 = g1.stack_fields(0.0, 0.0, e_inc, 1e-3, "top")[0]
    assert e1 == pytest.approx(e0 * np.exp(1j * g1.k0 * 1e-3), rel=1e-12)


def test_stack_fields_require_a_propagating_incident_wave_and_valid_port():
    _, g = make(2.2, 0.0, H)
    with pytest.raises(ValueError, match="propagate"):
        g.stack_fields(1.2 * g.k0, 0.0, np.array([0.0, 1.0]), 0.0)
    with pytest.raises(ValueError, match="port"):
        g.stack_fields(0.0, 0.0, np.array([1.0, 0.0]), 0.0, "side")
