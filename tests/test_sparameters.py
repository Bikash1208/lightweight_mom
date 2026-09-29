"""S-parameter extraction: analytic validation, energy conservation, reciprocity, plane conventions."""

import math
import warnings
from functools import cache

import numpy as np
import pytest
from reference_layered import solve_incident_fields

from src.basis_functions import MeshSettings
from src.config import config_from_dict
from src.constants import C0, ETA0
from src.excitation import ReferencePlanes
from src.geometry import Geometry, UnitCell
from src.green_function import LayeredGreenFunction
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice
from src.sparameters import (
    compute_s_parameters,
    decompose_air_fields,
    s_parameters_from_config,
    scattering_from_fields,
    stack_s_parameters,
)

F = 10e9
K0 = 2 * math.pi * F / C0
PX = PY = 8e-3
H = 1.575e-3
PLANES0 = ReferencePlanes(0.0, 0.0, H)  # planes on the substrate faces
PLANES = ReferencePlanes(5e-3, 7e-3, H)


def sub(tand=0.0, eps_r=2.2, h=H):
    return Substrate("t", eps_r, tand, h)


def green(tand=0.0, eps_r=2.2, h=H):
    return LayeredGreenFunction(sub(tand, eps_r, h), F)


def cell(kind, length, width, tand=0.0, eps_r=2.2, h=H, period=PX, **kw):
    return UnitCell(
        RectangularLattice(period, period), Geometry(kind, length, width, **kw), sub(tand, eps_r, h), Metal()
    )


MESH = MeshSettings(16, 16, 64)


@cache
def solved(kind, length, width, tand=0.0, rotation_deg=0.0, planes=PLANES0, nx=16, order=64):
    c = cell(kind, length, width, tand, rotation=math.radians(rotation_deg))
    return compute_s_parameters(c, F, MeshSettings(nx, nx, order), planes)


def airy(n, h, z_top=0.0, z_bottom=None):
    """Bare-slab S11 (=S22 up to plane shifts) and S21 (=S12): Fresnel r = (1-n)/(1+n), Airy sum.
    Amplitudes referenced to planes z_top / z_bottom (default: the two faces)."""
    z_bottom = -h if z_bottom is None else z_bottom
    r = (1 - n) / (1 + n)
    u = np.exp(-2j * K0 * n * h)
    s11 = r * (1 - u) / (1 - r * r * u) * np.exp(-2j * K0 * z_top)
    t = (1 - r * r) * np.exp(-1j * K0 * n * h) / (1 - r * r * u)
    # transmission between the faces, then free propagation along the air segments
    s21 = t * np.exp(-1j * K0 * z_top) * np.exp(+1j * K0 * (z_bottom + h))
    s22 = r * (1 - u) / (1 - r * r * u) * np.exp(+2j * K0 * (z_bottom + h))
    return s11, s21, s22


# ============================================================ decomposition of waves


def test_decompose_air_fields_recovers_up_and_down_going_plane_waves():
    e_up = np.array([0.7 - 0.2j, 0.1 + 0.4j])
    e_dn = np.array([-0.3 + 0.5j, 0.9j])
    h_up = np.array([-e_up[1], e_up[0]]) / ETA0  # +z wave:  H = z x E / eta0
    h_dn = -np.array([-e_dn[1], e_dn[0]]) / ETA0  # -z wave:  H = -z x E / eta0
    up, dn = decompose_air_fields(e_up + e_dn, h_up + h_dn)
    assert up == pytest.approx(e_up, rel=1e-13) and dn == pytest.approx(e_dn, rel=1e-13)
    up, dn = decompose_air_fields(e_up, h_up)
    assert up == pytest.approx(e_up) and np.abs(dn).max() < 1e-15
    up, dn = decompose_air_fields(e_dn, h_dn)
    assert dn == pytest.approx(e_dn) and np.abs(up).max() < 1e-15
    # batches
    e = np.stack([e_up, e_dn])
    h_ = np.stack([h_up, h_dn])
    up, dn = decompose_air_fields(e, h_)
    assert up.shape == dn.shape == (2, 2)


# ================================================ no metal: dielectric slab -> Fresnel


@pytest.mark.parametrize("tand", [0.0, 0.0009, 0.05])
@pytest.mark.parametrize("planes", [PLANES0, PLANES, ReferencePlanes(0.0, 12e-3, H), ReferencePlanes(3e-3, 0.0, H)])
def test_bare_slab_equals_fresnel_airy_closed_form(tand, planes):
    """No metal. Independent closed form: r = (1-n)/(1+n), u = exp(-2j k0 n h),
    S11 = r(1-u)/(1-r^2 u), S21 = (1-r^2) exp(-j k0 n h)/(1-r^2 u), times the free-space phase of the
    air segments between the faces and the reference planes."""
    s = stack_s_parameters(green(tand), planes)
    n = sub(tand).refractive_index
    s11, s21, s22 = airy(n, H, planes.z_top, planes.z_bottom)
    for pol in (0, 1):
        assert s.s11[pol, pol] == pytest.approx(s11, rel=1e-12, abs=1e-14)
        assert s.s21[pol, pol] == pytest.approx(s21, rel=1e-12, abs=1e-14)
        assert s.s12[pol, pol] == pytest.approx(s21, rel=1e-12, abs=1e-14)  # reciprocal
        assert s.s22[pol, pol] == pytest.approx(s22, rel=1e-12, abs=1e-14)
    for block in (s.s11, s.s12, s.s21, s.s22):
        assert block[0, 1] == 0.0 and block[1, 0] == 0.0  # no cross-pol
        assert block[0, 0] == pytest.approx(block[1, 1], rel=1e-14)  # x and y identical
    assert s.incoming_residual < 1e-12


def test_bare_slab_lossless_is_unitary_and_lossy_absorbs():
    s = stack_s_parameters(green(0.0), PLANES)
    m = s.matrix()
    assert m.conj().T @ m == pytest.approx(np.eye(4), abs=1e-12)
    assert np.allclose(s.power_out_per_incident(), 1.0, atol=1e-13)
    lossy = stack_s_parameters(green(0.05), PLANES)
    assert np.all(lossy.power_out_per_incident() < 1.0 - 1e-3)


def test_bare_slab_matches_independent_field_matching():
    """S11, S21 from the independent 3-D plane-wave solve: decompose ITS fields at the planes."""
    for tand in (0.0, 0.0009):
        sb, g = sub(tand), green(tand)
        planes = ReferencePlanes(4e-3, 6e-3, H)
        s = stack_s_parameters(g, planes)
        fields = solve_incident_fields(0.0, 0.0, np.array([1.0, 0.0, 0.0]), F, sb.eps_complex, H, "top")
        (e1, h1), (e2, h2) = fields(planes.z_top), fields(planes.z_bottom)
        up1, dn1 = decompose_air_fields(e1, h1)
        _, dn2 = decompose_air_fields(e2, h2)
        assert dn1[0] == pytest.approx(np.exp(1j * K0 * planes.z_top), rel=1e-12)  # unit at z = 0
        assert s.s11[0, 0] == pytest.approx(up1[0] / dn1[0], rel=1e-10)
        assert s.s21[0, 0] == pytest.approx(dn2[0] / dn1[0], rel=1e-10)
        fields_b = solve_incident_fields(0.0, 0.0, np.array([1.0, 0.0, 0.0]), F, sb.eps_complex, H, "bottom")
        (e1, h1), (e2, h2) = fields_b(planes.z_top), fields_b(planes.z_bottom)
        up1, _ = decompose_air_fields(e1, h1)
        up2, dn2 = decompose_air_fields(e2, h2)
        assert s.s12[0, 0] == pytest.approx(up1[0] / up2[0], rel=1e-10)
        assert s.s22[0, 0] == pytest.approx(dn2[0] / up2[0], rel=1e-10)


def test_half_wave_slab_is_transparent_and_matched_at_the_faces():
    """h = lambda_sub/2 (u = 1): S11 = 0 and S21 = exp(-j k0 * 0) times the reference phase = 1 at the faces."""
    n = math.sqrt(2.2)
    h = C0 / (F * n) / 2
    s = stack_s_parameters(green(0.0, h=h), ReferencePlanes(0.0, 0.0, h))
    assert abs(s.s11[0, 0]) < 1e-12
    assert s.s21[0, 0] == pytest.approx(-1.0, abs=1e-12)  # exp(-j k0 n h) = exp(-j pi) = -1, referenced to the faces


def test_no_dielectric_no_metal_is_free_space():
    """eps_r = 1: S11 = 0 and S21 = exp(-j k0 L), L = total distance between the planes."""
    planes = ReferencePlanes(5e-3, 7e-3, H)
    s = stack_s_parameters(green(0.0, eps_r=1.0), planes)
    assert abs(s.s11[0, 0]) < 1e-14
    assert s.s21[0, 0] == pytest.approx(np.exp(-1j * K0 * (5e-3 + H + 7e-3)), rel=1e-12)


# ============================================================ full PEC sheet


@pytest.mark.parametrize("tand", [0.0, 0.0009])
@pytest.mark.parametrize("planes", [PLANES0, PLANES])
@pytest.mark.parametrize("nx, order", [(4, 2), (8, 8), (16, 40)])
def test_full_pec_sheet_reflects_everything(tand, planes, nx, order):
    """S11 = -exp(-2j k0 z_top) * I (= -1 on the interface), S21 = S12 = 0. Uniform current is exact, so
    this holds for any mesh/order. S22 = shorted-line reflection seen through the substrate."""
    c = cell("patch", 8e-3, 8e-3, tand)
    s = compute_s_parameters(c, F, MeshSettings(nx, nx, order), planes)
    assert s.s11 == pytest.approx(-np.exp(-2j * K0 * planes.z_top) * np.eye(2), abs=1e-11)
    assert np.abs(s.s21).max() < 1e-11 and np.abs(s.s12).max() < 1e-11
    n = sub(tand).refractive_index
    z_in = 1j * (1 / n) * np.tan(K0 * n * H)  # short-circuited substrate line / eta0
    gamma = (z_in - 1) / (z_in + 1)  # reflection at the substrate's bottom face
    expected = gamma * np.exp(+2j * K0 * (planes.z_bottom + H))
    assert s.s22 == pytest.approx(expected * np.eye(2), abs=1e-11)
    if tand == 0.0:
        assert abs(gamma) == pytest.approx(1.0, abs=1e-13)
    assert s.incoming_residual < 1e-10


def test_full_pec_sheet_on_the_interface_gives_minus_one_exactly():
    s = compute_s_parameters(cell("patch", 8e-3, 8e-3, 0.0009), F, MeshSettings(8, 8, 8), PLANES0)
    assert s.component("S11_xx") == pytest.approx(-1.0, abs=1e-12)
    assert s.component("S11_yy") == pytest.approx(-1.0, abs=1e-12)
    assert abs(s.component("S11_yx")) < 1e-12
    assert s.phase_deg("S11_xx") == pytest.approx(180.0, abs=1e-6) or s.phase_deg("S11_xx") == pytest.approx(
        -180.0, abs=1e-6
    )
    assert s.magnitude_db("S21_xx") < -200.0 or s.magnitude("S21_xx") < 1e-11


# ================================================================ energy conservation

LOSSLESS_CELLS = [
    ("rectangle", 6e-3, 2e-3, 0.0),
    ("dipole", 6e-3, 1e-3, 30.0),
    ("dipole", 8e-3, 1e-3, 0.0),  # wire spanning the cell
    ("slot", 6e-3, 1e-3, 0.0),
    ("slot", 6e-3, 2e-3, 60.0),
    ("patch", 4e-3, 4e-3, 0.0),
]


@pytest.mark.parametrize("kind, length, width, rot", LOSSLESS_CELLS)
@pytest.mark.parametrize("planes", [PLANES0, PLANES])
def test_lossless_structure_conserves_energy(kind, length, width, rot, planes):
    """PEC + lossless substrate, only the (0,0) order propagates: S^H S = I over both ports and both
    polarizations, i.e. every incident wave leaves entirely as reflected + transmitted power."""
    s = solved(kind, length, width, 0.0, rot, planes)
    m = s.matrix()
    assert m.conj().T @ m == pytest.approx(np.eye(4), abs=1e-10)
    assert np.allclose(s.power_out_per_incident(), 1.0, atol=1e-10)
    assert s.incoming_residual < 1e-10


def test_lossy_structure_absorbs_and_never_creates_energy():
    s = solved("dipole", 6e-3, 1e-3, 0.05, 30.0)
    m = s.matrix()
    absorbed = np.eye(4) - m.conj().T @ m
    eig = np.linalg.eigvalsh(0.5 * (absorbed + absorbed.conj().T))
    assert eig.min() > -1e-10 and eig.max() > 1e-4
    assert np.all(s.power_out_per_incident() < 1.0)
    # a bit of loss changes little: tanδ = 0.0009 absorbs far less than 0.05
    weak = solved("dipole", 6e-3, 1e-3, 0.0009, 30.0)
    assert np.all(weak.power_out_per_incident() > s.power_out_per_incident())


# ==================================================================== reciprocity


@pytest.mark.parametrize("kind, length, width, rot", [LOSSLESS_CELLS[1], LOSSLESS_CELLS[4], LOSSLESS_CELLS[0]])
@pytest.mark.parametrize("tand", [0.0, 0.0009])
def test_reciprocal_structure_has_symmetric_S_matrix(kind, length, width, rot, tand):
    """Lorentz reciprocity: S = S^T on (port x polarization). Includes S12 = S21^T, S11 = S11^T (so
    S11_xy = S11_yx) and the same for S22, also with loss."""
    s = solved(kind, length, width, tand, rot)
    m = s.matrix()
    assert m == pytest.approx(m.T, abs=1e-10)
    assert s.s12 == pytest.approx(s.s21.T, abs=1e-10)
    assert s.s11 == pytest.approx(s.s11.T, abs=1e-10)
    assert s.s22 == pytest.approx(s.s22.T, abs=1e-10)
    if rot != 0.0:  # the test is not vacuous
        assert abs(s.s11[0, 1]) > 1e-2 and abs(s.s21[0, 1]) > 1e-2


def test_reciprocity_holds_for_the_bare_slab_and_unequal_plane_offsets():
    s = stack_s_parameters(green(0.0009), ReferencePlanes(1e-3, 9e-3, H))
    assert s.s12 == pytest.approx(s.s21.T, abs=1e-13)


# ===================================================================== symmetries


def test_mirror_symmetric_cell_has_no_cross_polarization():
    s = solved("rectangle", 6e-3, 2e-3, 0.0009)
    for block in (s.s11, s.s12, s.s21, s.s22):
        assert abs(block[0, 1]) < 1e-12 and abs(block[1, 0]) < 1e-12
    assert abs(s.s11[0, 0] - s.s11[1, 1]) > 1e-2  # anisotropic co-pol


def test_x_and_y_agree_for_a_square_symmetric_patch():
    s = solved("patch", 4e-3, 4e-3, 0.0009)
    for block in (s.s11, s.s12, s.s21, s.s22):
        assert block[0, 0] == pytest.approx(block[1, 1], rel=1e-10)
        assert abs(block[0, 1]) < 1e-12


@pytest.mark.parametrize("kind, length, width, rot", [("dipole", 6e-3, 1e-3, 30.0), ("slot", 6e-3, 2e-3, 20.0)])
def test_rotating_the_cell_by_90_degrees_rotates_the_polarization_basis(kind, length, width, rot):
    """S_B = R S_A R^T for every block (R = 90-degree rotation of the x,y basis), lossy included."""
    a = solved(kind, length, width, 0.0009, rot)
    b = solved(kind, length, width, 0.0009, rot + 90.0)
    r = np.array([[0.0, -1.0], [1.0, 0.0]])
    for name in ("s11", "s12", "s21", "s22"):
        assert getattr(b, name) == pytest.approx(r @ getattr(a, name) @ r.T, abs=1e-10)


# ============================================================ reference-plane conventions


def test_incident_amplitude_is_measured_at_the_reference_plane():
    """The amplitude the fields actually contain at the plane is exp(+j k0 z_top) (unit at z = 0, travelling -z)
    at port 1 and exp(+j k0 d_bottom) (unit at z = -h, travelling +z) at port 2."""
    planes = ReferencePlanes(5e-3, 7e-3, H)
    s = solved("dipole", 6e-3, 1e-3, 0.0, 30.0, planes)
    assert s.incident_amplitudes[0] == pytest.approx(np.exp(1j * K0 * planes.z_top) * np.ones(2), rel=1e-12)
    assert s.incident_amplitudes[1] == pytest.approx(np.exp(1j * K0 * 7e-3) * np.ones(2), rel=1e-12)
    assert s.incoming_residual < 1e-10  # nothing incoming except the specified wave


@pytest.mark.parametrize("kind, length, width, rot", [("dipole", 6e-3, 1e-3, 30.0), ("slot", 6e-3, 2e-3, 0.0)])
def test_moving_the_planes_multiplies_S_by_the_documented_phases(kind, length, width, rot):
    """S11' = S11 e^{-2jk0 D1}, S22' = S22 e^{+2jk0 D2}, S21' = S12' = S21 e^{+jk0 (D2 - D1)}; D1 = z_top' - z_top,
    D2 = z_bottom' - z_bottom. Direct recomputation at the new planes must agree, and |S| is unchanged."""
    a = solved(kind, length, width, 0.0009, rot, PLANES0)
    b = solved(kind, length, width, 0.0009, rot, PLANES)
    d1, d2 = PLANES.z_top - PLANES0.z_top, PLANES.z_bottom - PLANES0.z_bottom
    assert d1 == 5e-3 and d2 == pytest.approx(-7e-3)
    assert b.s11 == pytest.approx(a.s11 * np.exp(-2j * K0 * d1), abs=1e-11)
    assert b.s22 == pytest.approx(a.s22 * np.exp(+2j * K0 * d2), abs=1e-11)
    assert b.s21 == pytest.approx(a.s21 * np.exp(1j * K0 * (d2 - d1)), abs=1e-11)
    assert b.s12 == pytest.approx(a.s12 * np.exp(1j * K0 * (d2 - d1)), abs=1e-11)
    assert np.abs(b.matrix()) == pytest.approx(np.abs(a.matrix()), abs=1e-11)
    moved = a.at_planes(PLANES)
    for name in ("s11", "s12", "s21", "s22"):
        assert getattr(moved, name) == pytest.approx(getattr(b, name), abs=1e-11)
    assert moved.incident_amplitudes == pytest.approx(b.incident_amplitudes, abs=1e-11)
    back = moved.at_planes(PLANES0)
    assert back.matrix() == pytest.approx(a.matrix(), abs=1e-12)
    with pytest.raises(ValueError, match="thickness"):
        a.at_planes(ReferencePlanes(1e-3, 1e-3, 2 * H))


def test_exp_plus_jwt_convention_shows_as_negative_phase_slope_and_positive_group_delay():
    """Propagation delay appears as a phase that DEcreases with frequency (exp(+jwt)); moving the top plane
    outward by d subtracts 2 k0 d from arg S11 of a PEC sheet."""
    # PEC sheet: arg S11 = pi - 2 k0 z_top
    for z_top in (0.0, 1e-3, 3e-3):
        s = compute_s_parameters(cell("patch", 8e-3, 8e-3), F, MeshSettings(4, 4, 2), ReferencePlanes(z_top, 0.0, H))
        assert s.component("S11_xx") == pytest.approx(-np.exp(-2j * K0 * z_top), abs=1e-11)
    # group delay of the bare slab transmission: -d(arg S21)/d(omega) > 0 and close to (path in air + n h)/c
    freqs = np.array([9.9e9, 10.1e9])
    phases = []
    for f in freqs:
        g = LayeredGreenFunction(sub(0.0), f)
        phases.append(np.angle(stack_s_parameters(g, PLANES).s21[0, 0]))
    dphi = np.angle(np.exp(1j * (phases[1] - phases[0])))
    tau = -dphi / (2 * math.pi * (freqs[1] - freqs[0]))
    assert tau > 0.0
    assert tau == pytest.approx((5e-3 + 7e-3 + math.sqrt(2.2) * H) / C0, rel=0.1)


# =============================================================== convergence & limits


def test_S_parameters_converge_with_mesh_refinement():
    vals = [solved("dipole", 6e-3, 2e-3, 0.0009, 0.0, PLANES0, nx, 4 * nx) for nx in (8, 16, 32)]
    for name in ("S11_xx", "S21_xx"):
        v = [s.component(name) for s in vals]
        d = [abs(v[1] - v[0]), abs(v[2] - v[1])]
        assert 1.6 < d[0] / d[1] < 2.2  # first order (same as the mean current: S is linear in it)
        assert d[1] < 1.5e-2


def test_S_parameters_converge_with_floquet_order():
    vals = [solved("dipole", 6e-3, 2e-3, 0.0009, 0.0, PLANES0, 16, m).component("S11_xx") for m in (16, 32, 64, 128)]
    d = np.abs(np.diff(vals))
    assert np.all(np.diff(d) < 0.0) and np.all(d[:-1] / d[1:] > 3.0)  # ~ 1/M^2 tail
    assert d[-1] < 2e-4


def test_vanishing_metal_recovers_the_bare_slab():
    """A shrinking square patch (3, 2, 1 mm; area fraction 14% -> 1.6%) approaches the no-metal Fresnel result."""
    slab = stack_s_parameters(green(0.0009), PLANES0)
    errs = []
    for size in (3e-3, 2e-3, 1e-3):
        s = compute_s_parameters(cell("patch", size, size, 0.0009), F, MeshSettings(16, 16, 64), PLANES0)
        errs.append(abs(s.s21[0, 0] - slab.s21[0, 0]) + abs(s.s11[0, 0] - slab.s11[0, 0]))
    assert errs[0] > errs[1] > errs[2]
    assert errs[2] < 0.01


# ====================================================================== data model


def test_storage_layout_names_and_records():
    s = solved("dipole", 6e-3, 1e-3, 0.0, 30.0, PLANES)
    assert s.matrix().shape == (4, 4)
    assert s.matrix()[:2, :2] == pytest.approx(s.s11) and s.matrix()[2:, :2] == pytest.approx(s.s21)
    assert s.matrix()[:2, 2:] == pytest.approx(s.s12) and s.matrix()[2:, 2:] == pytest.approx(s.s22)
    # [received, incident]: S21_yx = y received at port 2 for x incident at port 1
    assert s.component("S21_yx") == s.s21[1, 0] and s.component("S21_xy") == s.s21[0, 1]
    assert s.component("S11_xx") == s.s11[0, 0] and s.component("S22_yy") == s.s22[1, 1]
    # co-polarized vs cross-polarized components are both stored and distinct for a rotated dipole
    assert abs(s.component("S21_xx")) > 0.1 and abs(s.component("S21_yx")) > 1e-2
    rec = s.to_record()
    assert rec["frequency_hz"] == F
    assert len(rec) == 1 + 4 * 4 * 2
    assert rec["S21_yx_mag"] == pytest.approx(abs(s.s21[1, 0]))
    assert rec["S21_yx_phase_deg"] == pytest.approx(math.degrees(np.angle(s.s21[1, 0])))
    assert -180.0 <= rec["S11_xx_phase_deg"] <= 180.0
    assert s.magnitude_db("S21_xx") == pytest.approx(20 * math.log10(abs(s.s21[0, 0])))
    for bad in ("S31_xx", "S11_zz", "S11", "s11_xx", "S11xx"):
        with pytest.raises(ValueError):
            s.component(bad)
    with pytest.raises(AttributeError):
        s.frequency = 1.0  # frozen
    assert s.planes == PLANES and s.k0 == pytest.approx(K0)


def test_mean_currents_are_stored_and_zero_without_metal():
    assert np.all(stack_s_parameters(green(), PLANES0).mean_currents == 0.0)
    s = solved("patch", 4e-3, 4e-3, 0.0)
    assert s.mean_currents.shape == (2, 2, 2) and np.abs(s.mean_currents[0, 0, 0]) > 1e-4


def test_input_validation_and_diagnostics():
    with pytest.raises(ValueError, match="thickness"):
        stack_s_parameters(green(), ReferencePlanes(1e-3, 1e-3, 2 * H))
    with pytest.raises(ValueError, match="shape"):
        scattering_from_fields(green(), PLANES0, np.zeros((2, 2)))
    # grating orders: period >= lambda0 (30 mm at 10 GHz) -> warn that S is specular-only
    big = UnitCell(RectangularLattice(40e-3, 40e-3), Geometry("patch", 20e-3, 20e-3), sub(0.0), Metal())
    with pytest.warns(UserWarning, match="higher Floquet order"):
        compute_s_parameters(big, F, MeshSettings(8, 8, 8), PLANES0)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        compute_s_parameters(cell("patch", 4e-3, 4e-3), F, MeshSettings(8, 8, 8), PLANES0)  # no warning


def test_config_route():
    d = {
        "period_x": 8e-3,
        "period_y": 8e-3,
        "substrate": {"eps_r": 2.2, "tan_delta": 0.0009, "thickness": 1.575e-3},
        "geometry": {"type": "dipole", "length": 6e-3, "width": 1e-3, "rotation_deg": 30.0},
        "excitation": {"frequency": 10e9},
        "mesh": {"cells_x": 16, "cells_y": 16, "floquet_order": 64},
        "reference_planes": {"top_air_thickness": 5e-3, "bottom_air_thickness": 7e-3},
    }
    cfg = config_from_dict(d)
    s = s_parameters_from_config(cfg)
    assert s.planes == cfg.reference_planes and s.planes.z_top == 5e-3
    direct = solved("dipole", 6e-3, 1e-3, 0.0009, 30.0, PLANES)
    assert s.matrix() == pytest.approx(direct.matrix(), abs=1e-12)
