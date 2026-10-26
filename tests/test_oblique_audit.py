"""Phase 20 R1 item 1, audit: which parts of the (kx0, ky0) path are ALREADY correct at kt < k0 (no solver change).

Each test checks one stage of the chain Green's function -> stack field -> excitation vector -> impedance matrix ->
current coefficients against an independent analytic/quadrature result, at oblique transverse wavenumbers."""

import numpy as np
import pytest

from src.basis_functions import RooftopBasis
from src.constants import C0, EPS0, ETA0, MU0
from src.excitation import PlaneWaveExcitation
from src.geometry import Geometry, UnitCell
from src.green_function import LayeredGreenFunction
from src.materials import Metal, Substrate
from src.mom_matrix import assemble_excitation_vector, assemble_impedance_matrix
from src.periodic_bc import RectangularLattice

F = 10e9
SUB = Substrate("rt5880", 2.2, 0.0, 1.575e-3)
ANGLES = (0.0, 30.0, 60.0)


def _kt(theta_deg, k0):
    return k0 * np.sin(np.radians(theta_deg))


@pytest.mark.parametrize("backing", ["pec", "air"])
@pytest.mark.parametrize("theta", ANGLES)
def test_audit_green_function_impedances_equal_textbook_transmission_line_values(theta, backing):
    g = LayeredGreenFunction(SUB, F, backing=backing)
    kt = _kt(theta, g.k0)
    a, b = g.normal_wavenumbers(kt, 0.0)
    w, eps, h = g.omega, SUB.eps_r, SUB.thickness
    # air: TE Z = w mu0 / a, TM Z = a / (w eps0); substrate: TE w mu0 / b, TM b / (w eps0 eps)
    za_te, za_tm = w * MU0 / a, a / (w * EPS0)
    zs_te, zs_tm = w * MU0 / b, b / (w * EPS0 * eps)
    if backing == "pec":
        zd_te, zd_tm = 1j * zs_te * np.tan(b * h), 1j * zs_tm * np.tan(b * h)
    else:  # substrate loaded by semi-infinite air below
        zd_te = zs_te * (za_te + 1j * zs_te * np.tan(b * h)) / (zs_te + 1j * za_te * np.tan(b * h))
        zd_tm = zs_tm * (za_tm + 1j * zs_tm * np.tan(b * h)) / (zs_tm + 1j * za_tm * np.tan(b * h))
    want_te, want_tm = za_te * zd_te / (za_te + zd_te), za_tm * zd_tm / (za_tm + zd_tm)
    got_te, got_tm = g.equivalent_impedances(kt, 0.0)
    assert abs(got_te - want_te) < 1e-9 * abs(want_te) and abs(got_tm - want_tm) < 1e-9 * abs(want_tm)


@pytest.mark.parametrize("theta", ANGLES)
def test_audit_grounded_stack_reflection_equals_fresnel_airy_te_and_tm(theta):
    g = LayeredGreenFunction(SUB, F, backing="pec")
    kt = _kt(theta, g.k0)
    a, b = g.normal_wavenumbers(kt, 0.0)
    w, eps, h = g.omega, SUB.eps_r, SUB.thickness
    # incident along the plane of incidence phi = 0: TM has E_t along p = x, TE along q = y
    for e_inc, z_air, z_sub in (
        ((0.0, 1.0), w * MU0 / a, w * MU0 / b),
        ((1.0, 0.0), a / (w * EPS0), b / (w * EPS0 * eps)),
    ):
        z_in = 1j * z_sub * np.tan(b * h)
        gamma = (z_in - z_air) / (z_in + z_air)
        e, hh = g.stack_fields(kt, 0.0, np.array(e_inc, dtype=complex), 1e-3, "top")  # 1 mm above the pattern plane
        # up/down split of the mode: V = E.vec, I from H (sign conventions of _apply_modes)
        te = e_inc == (0.0, 1.0)
        v = e[1] if te else e[0]
        i = -hh[0] if te else hh[1]
        y = 1.0 / z_air
        v_up, v_dn = 0.5 * (v + i / y), 0.5 * (v - i / y)
        # 1 mm above: down-going carries exp(+j a z), up-going exp(-j a z); refer both to z = 0
        assert abs(v_up * np.exp(1j * a * 1e-3) / (v_dn * np.exp(-1j * a * 1e-3)) - gamma) < 1e-10


@pytest.mark.parametrize("theta", ANGLES)
def test_audit_ungrounded_slab_reflection_and_transmission_energy(theta):
    g = LayeredGreenFunction(SUB, F, backing="air")
    kt = _kt(theta, g.k0)
    a, b = g.normal_wavenumbers(kt, 0.0)
    for e_inc, y_air in (((0.0, 1.0), a / (g.omega * MU0)), ((1.0, 0.0), g.omega * EPS0 / a)):
        zs = [-1e-3 - SUB.thickness]
        ez, hz = g.stack_fields(kt, 0.0, np.array(e_inc, dtype=complex), zs, "top")
        te = e_inc == (0.0, 1.0)
        v = ez[0, 1] if te else ez[0, 0]
        i = -hz[0, 0] if te else hz[0, 1]
        t_wave = 0.5 * (v - i / y_air) * np.exp(-1j * a * 1e-3)  # down-going transmitted amplitude at z = -h
        # reflected amplitude from the same field above: power balance |R|^2 + |T|^2 = 1 (lossless, same medium both sides)
        e_top, h_top = g.stack_fields(kt, 0.0, np.array(e_inc, dtype=complex), [1e-3], "top")
        vt = e_top[0, 1] if te else e_top[0, 0]
        it = -h_top[0, 0] if te else h_top[0, 1]
        up = 0.5 * (vt + it / y_air) * np.exp(1j * a * 1e-3)
        dn = 0.5 * (vt - it / y_air) * np.exp(-1j * a * 1e-3)
        assert abs(dn - 1.0) < 1e-10  # unit incident tangential E at the interface
        assert abs(abs(up) ** 2 + abs(t_wave) ** 2 - 1.0) < 1e-10


def _cell(px=10e-3, size=6e-3):
    return UnitCell(
        RectangularLattice(px, px), Geometry(kind="rectangle", length=size, width=size), SUB, Metal(), backing="pec"
    )


@pytest.mark.parametrize("theta", (30.0, 60.0))
def test_audit_excitation_vector_is_the_quadrature_of_the_rooftop_against_the_incident_phase(theta):
    uc = _cell()
    basis = RooftopBasis(uc, 20, 20)
    kx0, ky0 = (
        _kt(theta, 2 * np.pi * F / C0) * np.cos(np.radians(40.0)),
        _kt(theta, 2 * np.pi * F / C0) * np.sin(np.radians(40.0)),
    )
    e_t = np.array([0.3 + 0.1j, -0.7 + 0.2j])
    v = assemble_excitation_vector(basis, e_t, kx0, ky0)
    # direct quadrature of  integral f_b(r) . E_t exp(-j k0 . r) dA  for a few interior basis functions
    n = 400
    xs = (np.arange(n) + 0.5) / n * 10e-3 - 5e-3
    gx, gy = np.meshgrid(xs, xs, indexing="ij")
    phase = np.exp(-1j * (kx0 * gx + ky0 * gy))
    for ib in (0, basis.n_x // 2, basis.n_x + 5, basis.n_unknowns - 1):
        c = np.zeros(basis.n_unknowns, complex)
        c[ib] = 1.0
        j = basis.evaluate(gx, gy, c, kx0, ky0)  # Bloch rooftop on the grid
        quad = np.sum((j[..., 0] * e_t[0] + j[..., 1] * e_t[1]) * phase) * (10e-3 / n) ** 2
        assert abs(quad - v[ib]) < 2e-3 * max(abs(v[ib]), 1e-12) + 1e-12


@pytest.mark.parametrize("theta", (30.0, 60.0))
def test_audit_current_coefficients_equal_the_floquet_integral_of_the_bloch_current(theta):
    uc = _cell()
    basis = RooftopBasis(uc, 20, 20)
    kt = _kt(theta, 2 * np.pi * F / C0)
    kx0, ky0 = kt * np.cos(0.7), kt * np.sin(0.7)
    rng = np.random.default_rng(0)
    coeffs = rng.normal(size=basis.n_unknowns) + 1j * rng.normal(size=basis.n_unknowns)
    jt = basis.current_coefficients(coeffs, kx0, ky0)
    n = 400
    xs = (np.arange(n) + 0.5) / n * 10e-3 - 5e-3
    gx, gy = np.meshgrid(xs, xs, indexing="ij")
    j = basis.evaluate(gx, gy, coeffs, kx0, ky0)
    quad = np.sum(j * np.exp(1j * (kx0 * gx + ky0 * gy))[..., None], axis=(0, 1)) * (10e-3 / n) ** 2 / 10e-3**2
    assert np.allclose(quad, jt, rtol=3e-3, atol=3e-3 * np.abs(jt).max())


@pytest.mark.parametrize("theta", (30.0, 60.0))
def test_audit_impedance_matrix_obeys_reciprocity_at_oblique_incidence(theta):
    uc = _cell()
    basis = RooftopBasis(uc, 12, 12)
    g = LayeredGreenFunction(SUB, F, backing="pec")
    kt = _kt(theta, g.k0)
    kx0, ky0 = kt * np.cos(0.5), kt * np.sin(0.5)
    z_p = assemble_impedance_matrix(basis, g, kx0, ky0, 8)
    z_m = assemble_impedance_matrix(basis, g, -kx0, -ky0, 8)
    assert np.allclose(z_p.T, z_m, rtol=1e-10, atol=1e-10 * np.abs(z_p).max())
    assert not np.allclose(z_p, z_p.T, rtol=1e-6)  # genuinely non-symmetric at kt != 0: the check above is not vacuous


def test_audit_plane_wave_excitation_accepts_oblique_x_or_y_only_in_the_plane_of_incidence_family():
    PlaneWaveExcitation(F, theta=np.radians(30.0), phi=0.0, polarization="y")  # TE
    PlaneWaveExcitation(F, theta=np.radians(30.0), phi=np.radians(90.0), polarization="x")  # TE at phi = 90
    with pytest.raises(ValueError):  # x is not perpendicular to k at phi = 0: TM is not an x/y label
        PlaneWaveExcitation(F, theta=np.radians(30.0), phi=0.0, polarization="x")
    assert ETA0 > 0
