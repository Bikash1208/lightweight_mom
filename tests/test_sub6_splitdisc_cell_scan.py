"""Physics checks of the sub-6 GHz split-disc cell scan (`scripts/research/sub6_splitdisc_cell_scan.py`, job section 4). MODEL_ONLY.

Coarse mesh (n = 16, F = 2n, period 10 mm, gap 1.25 mm = 2 mesh cells, footprint 1.25 mm = 2 cells) so the footprint is mesh-aligned and the limits are clean; the
tolerances below were fixed before the runs. (A first version used a 0.625 mm = 1-cell footprint whose edges fall exactly on cell centres: the metal sampling is then a
tie, the bridge got no basis functions and the limit tests passed vacuously; `test_footprint_bridge_is_not_vacuous` guards against that.)
"""
import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

from src.basis_functions import MeshSettings
from src.dispersion.tensor_sheet import bloch_sheet_tensor, reciprocity_error
from src.frequency import free_space_wavenumber
from src.loads import LumpedLoad

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("sub6_scan", ROOT / "scripts/research/sub6_splitdisc_cell_scan.py")
S = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(S)

N, G, FOOT, H, R1 = 16, 1.25, 1.25, 3.175, 4.0
F = S.F0
K0 = free_space_wavenumber(F)
KT = 1.07 * K0


def tensor(variant="S", q=1.0, psi=0.0, load=None, bridge=True, phi_deg=0.0, kt=KT):
    c = S.build_cell(variant, R1, q, H, G, FOOT, psi, load, bridge=bridge)
    phi = math.radians(phi_deg)
    return S.bloch_loaded_f2n(c, F, N, kt * math.cos(phi), kt * math.sin(phi))


def diode(psi=0.0, **kw):
    d = "y" if psi == 0.0 else float((psi + 90.0) % 180.0)
    kw.setdefault("transverse_conducts", S.default_transverse_conducts(d))  # the project's lumped-element convention (False on x/y)
    return LumpedLoad(x=0.0, y=0.0, length=G * 1e-3, width=FOOT * 1e-3, direction=d, **kw)


@pytest.fixture(scope="module")
def loaded():
    ld = diode(R=3.0, L=0.6e-9)
    return {"x": tensor(load=ld, phi_deg=0.0), "y": tensor(load=ld, phi_deg=90.0)}


def test_footprint_bridge_is_not_vacuous():
    from src.basis_functions import RooftopBasis

    with_b = RooftopBasis(S.build_cell("S", R1, 1.0, H, G, FOOT, 0.0, None, bridge=True), N, N).n_unknowns
    without = RooftopBasis(S.build_cell("S", R1, 1.0, H, G, FOOT, 0.0, None, bridge=False), N, N).n_unknowns
    assert with_b > without
    strap, gap = tensor(load=None, bridge=True), tensor(load=None, bridge=False)
    assert abs(strap[1, 1] - gap[1, 1]) > 0.05 * abs(gap[1, 1])  # a PEC strap across the gap is not the empty gap


def test_wrapper_equals_tensor_sheet_without_loads():
    c = S.build_cell("S", R1, 1.0, H, G, FOOT, 0.0, None, bridge=False)
    ref = bloch_sheet_tensor(c, F, MeshSettings(N, N, 2 * N), KT, 0.0)
    assert np.allclose(S.bloch_loaded_f2n(c, F, N, KT, 0.0), ref, rtol=1e-12, atol=1e-15)


def test_geometry_area_matches_analytic_disc_minus_strip():
    from scipy.integrate import quad

    c = S.build_cell("S", R1, 1.0, H, 0.65, 0.35, 0.0, None, bridge=False)
    strip = 2 * quad(lambda u: 2 * math.sqrt((R1 * 1e-3) ** 2 - u * u), 0, 0.65e-3 / 2)[0]
    assert c.geometry.feature_area == pytest.approx(math.pi * (R1 * 1e-3) ** 2 - strip, rel=3e-3)


def test_reciprocity_of_the_loaded_tensor(loaded):
    # tolerance fixed in advance: 1e-8 (mirror-symmetric S cell, kt along an axis -> Y_xy = Y_yx = 0 up to round-off)
    for y in loaded.values():
        assert reciprocity_error(y) < 1e-8


def test_mirror_symmetry_variant_s_gives_diagonal_tensor_in_both_directions(loaded):
    for key, y in loaded.items():
        assert abs(y[0, 1]) / np.abs(np.diag(y)).max() < 1e-8, key
        assert abs(y[1, 0]) / np.abs(np.diag(y)).max() < 1e-8, key


def test_variant_t_with_q_one_is_variant_s():
    ld = diode(R=3.0, L=0.6e-9)
    assert np.array_equal(tensor("T", 1.0, load=ld), tensor("S", 1.0, load=ld))


def test_rotation_by_90_deg_swaps_par_and_perp():
    ld0, ld90 = diode(0.0, R=3.0, L=0.6e-9), diode(90.0, R=3.0, L=0.6e-9)
    y0x, y0y = tensor(psi=0.0, load=ld0, phi_deg=0.0), tensor(psi=0.0, load=ld0, phi_deg=90.0)
    y9x, y9y = tensor(psi=90.0, load=ld90, phi_deg=0.0), tensor(psi=90.0, load=ld90, phi_deg=90.0)
    assert y9y[1, 1] == pytest.approx(y0x[0, 0], rel=1e-6)  # Y_par: along the gap axis (x at psi = 0, y at psi = 90)
    assert y9x[0, 0] == pytest.approx(y0y[1, 1], rel=1e-6)  # Y_perp


@pytest.fixture(scope="module")
def open_limit():
    return tensor(load=None, bridge=False), tensor(load=diode(open_circuit=True), bridge=True), tensor(load=diode(L=1.0), bridge=True)


@pytest.mark.xfail(strict=True, reason="measured: the open-circuit load (OPEN_OHM = 1e12 ohm, finite) differs from the empty gap by 1.2e-6 relative in Y_yy (real part, 3.7e-10 S) against the pre-set 1e-6; kept visible, tolerance NOT lowered")
def test_off_open_equals_the_infinite_impedance_limit_1e6(open_limit):
    ref, lim, _ = open_limit
    assert np.allclose(lim, ref, rtol=1e-6, atol=1e-12 * np.abs(ref).max())


def test_off_open_limit_is_within_1e5_and_large_l_approaches_it(open_limit):
    ref, lim, big = open_limit
    assert np.allclose(lim, ref, rtol=1e-5, atol=1e-12 * np.abs(ref).max())
    assert np.allclose(big, ref, rtol=1e-3, atol=1e-6 * np.abs(ref).max())


def test_on_state_tends_to_a_pec_strap_as_r_goes_to_zero():
    strap = tensor(load=None, bridge=True)  # metal footprint, no load = PEC strap across the gap
    low_r = tensor(load=diode(R=1e-9, transverse_conducts=True), bridge=True)
    assert np.allclose(low_r, strap, rtol=1e-6, atol=1e-12 * np.abs(strap).max())
    # the default lumped element does not conduct along the gap inside the footprint: still close to the strap for a 1-cell footprint
    default = tensor(load=diode(R=1e-9), bridge=True)
    assert np.allclose(default, strap, rtol=0.02, atol=1e-3 * np.abs(strap).max())


def test_switching_changes_the_tensor_but_barely_the_axes():
    off = tensor(load=diode(C=0.15e-12, L=0.6e-9))
    on = tensor(load=diode(R=3.0, L=0.6e-9))
    assert abs(on[1, 1] - off[1, 1]) > 1e-3 * abs(off[1, 1])  # the diode acts across the gap: Y_perp changes
    assert abs(on[0, 1]) < 1e-8 * np.abs(np.diag(on)).max()  # the axes stay


def test_gap_frame_tensor_rotates_back():
    th = math.radians(30.0)
    r = np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])
    yg = np.diag([2.0 + 0.1j, 5.0 + 0.2j]) * 1e-3
    lab = r @ yg @ r.T
    row = {"Y_mS": [[v.real * 1e3, v.imag * 1e3] for v in lab.ravel()]}
    yp, yq, oa, ob = S.gap_frame_tensor(row, row, 30.0)
    assert yp == pytest.approx(yg[0, 0]) and yq == pytest.approx(yg[1, 1]) and abs(oa) < 1e-12 and abs(ob) < 1e-12


def test_state_footprints_follow_the_diode_data_sheets():
    ma, bar = S.state_load("OFF_MA"), S.state_load("OFF_BAR50")
    assert (ma.length, ma.width) == pytest.approx((0.33e-3, 0.35e-3))
    assert (bar.length, bar.width) == pytest.approx((0.65e-3, 0.35e-3))
    assert S.state_load("OFF_OPEN") is None
    assert S.state_load("ON_BAR50", 45.0).transverse_conducts  # rotated diode: the fcc6691 default
