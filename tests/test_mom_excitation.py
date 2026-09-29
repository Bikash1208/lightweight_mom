"""Excitation vector V and the unpatterned-stack field; reference planes."""

import math
from dataclasses import dataclass

import numpy as np
import pytest

from reference_layered import solve_incident_plane_wave
from src.basis_functions import MeshSettings, RooftopBasis
from src.config import config_from_dict
from src.constants import C0
from src.excitation import PlaneWaveExcitation, ReferencePlanes, stack_tangential_field
from src.geometry import Geometry, UnitCell
from src.green_function import LayeredGreenFunction
from src.materials import Metal, Substrate
from src.mom_matrix import assemble_excitation_vector
from src.periodic_bc import RectangularLattice
from src.solver import build_system

F = 10e9
PX = PY = 8e-3


@dataclass
class PolarizedStub:
    """Duck-typed excitation with an arbitrary transverse E direction (used for TM incidence,
    which PlaneWaveExcitation's x/y polarizations cannot express)."""

    frequency: float
    theta: float
    phi: float
    tm: bool

    @property
    def transverse_wavenumber(self):
        k0 = 2 * math.pi * self.frequency / C0
        return (k0 * math.sin(self.theta) * math.cos(self.phi), k0 * math.sin(self.theta) * math.sin(self.phi))

    @property
    def e_field_direction(self):
        c, s = math.cos(self.theta), math.sin(self.theta)
        cp, sp = math.cos(self.phi), math.sin(self.phi)
        if self.tm:            # E in the plane of incidence, perpendicular to k_inc
            return np.array([c * cp, c * sp, s])
        return np.array([-sp, cp, 0.0])     # TE: E perpendicular to the plane of incidence


def green(eps_r, tand, h):
    return LayeredGreenFunction(Substrate("s", eps_r, tand, h), F)


# ============================================ stack field vs independent field matching

STACKS = [(2.2, 0.0, 1.575e-3), (2.2, 0.0009, 1.575e-3), (10.2, 0.05, 3e-3), (4.4, 0.0, 0.4e-3)]


@pytest.mark.parametrize("eps_r, tand, h", STACKS)
@pytest.mark.parametrize("pol", ["x", "y"])
def test_normal_incidence_stack_field_matches_field_matching(eps_r, tand, h, pol):
    g = green(eps_r, tand, h)
    exc = PlaneWaveExcitation(F, polarization=pol)
    ref = solve_incident_plane_wave(0.0, 0.0, exc.e_field_direction, F, g.substrate.eps_complex, h)
    assert stack_tangential_field(g, exc) == pytest.approx(ref, rel=1e-12, abs=1e-14)


@pytest.mark.parametrize("eps_r, tand, h", STACKS)
@pytest.mark.parametrize("theta_deg", [15.0, 40.0, 65.0])
@pytest.mark.parametrize("tm", [False, True])
@pytest.mark.parametrize("phi_deg", [0.0, 33.0, 110.0])
def test_oblique_te_and_tm_stack_field_matches_field_matching(eps_r, tand, h, theta_deg, tm, phi_deg):
    g = green(eps_r, tand, h)
    exc = PolarizedStub(F, math.radians(theta_deg), math.radians(phi_deg), tm)
    kx, ky = exc.transverse_wavenumber
    e_dir = exc.e_field_direction
    k_inc = np.array([kx, ky, -g.k0 * math.cos(exc.theta)])
    assert abs(e_dir @ k_inc) < 1e-9 * g.k0                       # the stub really is transverse
    ref = solve_incident_plane_wave(kx, ky, e_dir, F, g.substrate.eps_complex, h)
    assert stack_tangential_field(g, exc) == pytest.approx(ref, rel=1e-11, abs=1e-13)


def test_plane_wave_excitation_te_oblique_matches_reference():
    """The real PlaneWaveExcitation class (x polarization, k in the y-z plane = TE)."""
    g = green(2.2, 0.0009, 1.575e-3)
    exc = PlaneWaveExcitation(F, theta=math.radians(35), phi=math.pi / 2, polarization="x")
    kx, ky = exc.transverse_wavenumber
    ref = solve_incident_plane_wave(kx, ky, exc.e_field_direction, F, g.substrate.eps_complex, 1.575e-3)
    assert stack_tangential_field(g, exc) == pytest.approx(ref, rel=1e-12, abs=1e-14)


# ------------------------------------------------------------- thin-film analytics

def test_half_wave_slab_is_transparent_and_quarter_wave_slab_gives_two_over_one_plus_eps():
    """Lossless normal incidence. Half-wave slab: Gamma = 0, E(0) = E_inc. Quarter-wave slab:
    Gamma = (1 - n^2)/(1 + n^2), so E(0) = 1 + Gamma = 2/(1 + eps_r) (thin-film optics, independent
    of the transmission-line machinery)."""
    eps_r = 2.2
    n = math.sqrt(eps_r)
    lam_sub = C0 / (F * n)
    exc = PlaneWaveExcitation(F)
    half = stack_tangential_field(green(eps_r, 0.0, lam_sub / 2), exc)
    quarter = stack_tangential_field(green(eps_r, 0.0, lam_sub / 4), exc)
    assert half == pytest.approx([1.0, 0.0], abs=1e-12)
    assert quarter == pytest.approx([2 / (1 + eps_r), 0.0], abs=1e-12)


def test_no_slab_limits():
    exc = PlaneWaveExcitation(F)
    assert stack_tangential_field(green(1.0, 0.0, 1.575e-3), exc) == pytest.approx([1.0, 0.0], abs=1e-14)
    assert stack_tangential_field(green(2.2, 0.0, 1e-12), exc) == pytest.approx([1.0, 0.0], abs=1e-9)


def test_thick_lossy_substrate_is_a_single_interface():
    """h -> infinity with loss: E(0) = 1 + r01, r01 = (1 - n)/(1 + n) at normal incidence."""
    g = green(2.2, 0.05, 5.0)              # far-face round trip ~ exp(-78): negligible
    n = g.substrate.refractive_index
    e = stack_tangential_field(g, PlaneWaveExcitation(F))
    assert e[0] == pytest.approx(1 + (1 - n) / (1 + n), rel=1e-12)


def test_frequency_mismatch_is_rejected():
    with pytest.raises(ValueError, match="frequency"):
        stack_tangential_field(green(2.2, 0.0, 1e-3), PlaneWaveExcitation(11e9))


# ================================================================ excitation vector

def make_basis(kind="rectangle", length=4e-3, width=2e-3):
    sub = Substrate("s", 2.2, 0.0009, 1.575e-3)
    c = UnitCell(RectangularLattice(PX, PY), Geometry(kind, length, width), sub, Metal())
    return RooftopBasis(c, 8, 8)


@pytest.mark.parametrize("pol", ["x", "y"])
def test_normal_incidence_vector_is_area_times_stack_field_on_matching_rooftops(pol):
    b = make_basis()
    g = green(2.2, 0.0009, 1.575e-3)
    exc = PlaneWaveExcitation(F, polarization=pol)
    e = stack_tangential_field(g, exc)
    v = assemble_excitation_vector(b, e, 0.0, 0.0)
    comp = 0 if pol == "x" else 1
    same = b.direction == comp
    assert v[same] == pytest.approx(b.dx * b.dy * e[comp], rel=1e-14)       # integral of f = dx*dy
    assert np.all(v[~same] == 0.0)                                           # E has no cross component
    assert v.shape == (b.n_unknowns,)


def test_oblique_vector_carries_phase_and_sinc_amplitude():
    """V_i = conj(Phi_i(k0)) E_d: a phase exp(-j k0.r_i) and the rooftop's spectral amplitude."""
    b = make_basis()
    g = green(2.2, 0.0009, 1.575e-3)
    exc = PlaneWaveExcitation(F, theta=math.radians(40), phi=math.pi / 2, polarization="x")
    kx0, ky0 = exc.transverse_wavenumber
    e = stack_tangential_field(g, exc)
    v = assemble_excitation_vector(b, e, kx0, ky0)
    xs = b.direction == 0
    amp = b.dx * b.dy * np.sinc(kx0 * b.dx / (2 * math.pi)) ** 2 * np.sinc(ky0 * b.dy / (2 * math.pi))
    phase = np.exp(-1j * (kx0 * b.center_x[xs] + ky0 * b.center_y[xs]))
    assert v[xs] == pytest.approx(amp * phase * e[0], rel=1e-13)
    assert np.all(np.abs(v[~xs]) < 1e-12 * np.abs(v[xs]).max())              # E_y = 0 for this TE wave


def test_vector_is_projection_of_the_incident_field_onto_the_basis():
    """V_i = integral f_i . E_stack e^{-j k0.r} dA, checked by quadrature for an oblique wave."""
    b = make_basis("rectangle", 4e-3, 4e-3)
    g = green(2.2, 0.0009, 1.575e-3)
    exc = PlaneWaveExcitation(F, theta=math.radians(30), phi=math.pi / 2, polarization="x")
    kx0, ky0 = exc.transverse_wavenumber
    e = stack_tangential_field(g, exc)
    v = assemble_excitation_vector(b, e, kx0, ky0)
    n = 512
    x = -PX / 2 + (np.arange(n) + 0.5) * PX / n
    X, Y = np.meshgrid(x, x, indexing="ij")
    incident = e[None, None, :] * np.exp(-1j * (kx0 * X + ky0 * Y))[..., None]
    for k in (0, b.n_x // 2, b.n_x - 1):
        unit = np.zeros(b.n_unknowns, complex)
        unit[k] = 1.0
        f = b.evaluate(X, Y, unit)
        num = np.sum(np.sum(f * incident, axis=-1)) * (PX / n) ** 2
        assert v[k] == pytest.approx(num, rel=3e-4, abs=1e-6 * abs(v[k]))


# ============================================================== reference planes

def test_reference_plane_positions_and_validation():
    rp = ReferencePlanes(10e-3, 5e-3, 1.575e-3)
    assert rp.z_top == 10e-3
    assert rp.z_bottom == pytest.approx(-(1.575e-3 + 5e-3), rel=1e-15)
    assert ReferencePlanes(0.0, 0.0, 1e-3).z_bottom == -1e-3        # planes on the interfaces
    for bad in [(-1e-3, 0, 1e-3), (0, -1e-3, 1e-3), (0, 0, 0.0), (float("nan"), 0, 1e-3), ("1", 0, 1e-3)]:
        with pytest.raises(ValueError):
            ReferencePlanes(*bad)


def config_dict(top, bottom):
    return {
        "period_x": 8e-3, "period_y": 8e-3,
        "substrate": {"eps_r": 2.2, "tan_delta": 0.0009, "thickness": 1.575e-3},
        "geometry": {"type": "rectangle", "length": 4e-3, "width": 2e-3},
        "excitation": {"frequency": 10e9},
        "mesh": {"cells_x": 8, "cells_y": 8, "floquet_order": 8},
        "reference_planes": {"top_air_thickness": top, "bottom_air_thickness": bottom},
    }


def test_reference_planes_do_not_change_Z_or_V_but_are_kept_explicit():
    """The air regions are semi-infinite in the Green's function; the planes are metadata for
    later de-embedding. Z and V must be bit-identical for different plane locations."""
    a = config_from_dict(config_dict(1e-3, 2e-3))
    b = config_from_dict(config_dict(30e-3, 0.0))
    sa, sb = build_system(a), build_system(b)
    assert np.array_equal(sa.z, sb.z) and np.array_equal(sa.v, sb.v)
    assert a.reference_planes.z_top != b.reference_planes.z_top
    assert a.reference_planes.z_bottom == pytest.approx(-(1.575e-3 + 2e-3))
    assert b.reference_planes.z_bottom == pytest.approx(-1.575e-3)


# ============================================================ bottom-port excitation (Phase 4)

@pytest.mark.parametrize("eps_r, tand, h", STACKS)
@pytest.mark.parametrize("pol", ["x", "y"])
def test_bottom_port_stack_field_at_z0_matches_field_matching(eps_r, tand, h, pol):
    """Incident from z < -h (unit amplitude at z = -h): the field on the sheet plane is the transmitted wave."""
    from reference_layered import solve_incident_fields

    g = green(eps_r, tand, h)
    exc = PlaneWaveExcitation(F, polarization=pol)
    fields = solve_incident_fields(0.0, 0.0, exc.e_field_direction, F, g.substrate.eps_complex, h, "bottom")
    assert stack_tangential_field(g, exc, "bottom") == pytest.approx(fields(0.0)[0], rel=1e-12, abs=1e-14)


def test_bottom_port_normal_incidence_equals_the_slab_transmission_coefficient():
    """E(0) = (1 - r^2) exp(-j k0 n h) / (1 - r^2 u): thin-film transmission between the two faces."""
    g = green(2.2, 0.0009, 1.575e-3)
    n = g.substrate.refractive_index
    r = (1 - n) / (1 + n)
    u = np.exp(-2j * g.k0 * n * 1.575e-3)
    expected = (1 - r * r) * np.exp(-1j * g.k0 * n * 1.575e-3) / (1 - r * r * u)
    e = stack_tangential_field(g, PlaneWaveExcitation(F), "bottom")
    assert e == pytest.approx([expected, 0.0], rel=1e-12, abs=1e-14)


def test_bottom_port_of_a_slab_free_stack_is_the_incident_wave_propagated_by_h():
    g = green(1.0, 0.0, 1.575e-3)
    e = stack_tangential_field(g, PlaneWaveExcitation(F), "bottom")
    assert e == pytest.approx([np.exp(-1j * g.k0 * 1.575e-3), 0.0], rel=1e-12)   # travelling +z: exp(-jk0 h)


def test_invalid_port_is_rejected():
    with pytest.raises(ValueError, match="port"):
        stack_tangential_field(green(2.2, 0.0, 1e-3), PlaneWaveExcitation(F), "left")
