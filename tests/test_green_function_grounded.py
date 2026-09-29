"""Phase 12A: PEC ground plane (``backing="pec"``) in ``LayeredGreenFunction``.

Two kinds of checks, matching the rigor of the existing ungrounded tests:
(1) analytic/property checks against the closed forms and the physics they must obey;
(2) agreement with the independent 3-D plane-wave field-matching solver
    (``reference_layered.solve_current_sheet_grounded`` / ``solve_incident_fields_grounded``),
    which shares no derivation with ``src/green_function.py``.
"""

import numpy as np
import pytest
from reference_layered import solve_current_sheet_grounded, solve_incident_fields_grounded

from src.green_function import LayeredGreenFunction
from src.materials import Substrate

F = 12.5e9
H = 1.575e-3
LOSSLESS = Substrate("x", 2.2, 0.0, H)
LOSSY = Substrate("x", 2.2, 0.0009, H)


def green(sub=LOSSLESS, f=F, backing="pec"):
    return LayeredGreenFunction(sub, f, backing=backing)


# --------------------------------------------------------------------- backward compatibility


def test_default_backing_is_air():
    assert LayeredGreenFunction(LOSSLESS, F).backing == "air"


def test_default_backing_bit_identical_to_explicit_air():
    g1 = LayeredGreenFunction(LOSSLESS, F)
    g2 = LayeredGreenFunction(LOSSLESS, F, backing="air")
    kx, ky = np.array([0.0, 50.0, 300.0]), np.array([0.0, 10.0, -20.0])
    assert np.array_equal(g1.spectral_dyad(kx, ky), g2.spectral_dyad(kx, ky))
    z = np.array([0.0, 0.005, -0.001, -0.02])
    e_inc = np.array([1.0, 0.0], dtype=complex)
    e1, h1 = g1.stack_fields(0.0, 0.0, e_inc, z)
    e2, h2 = g2.stack_fields(0.0, 0.0, e_inc, z)
    assert np.array_equal(e1, e2) and np.array_equal(h1, h2)


def test_invalid_backing_rejected():
    with pytest.raises(ValueError):
        LayeredGreenFunction(LOSSLESS, F, backing="copper")


# --------------------------------------------------------------------------------- reflection


@pytest.mark.parametrize("f", [10.0e9, 11.25e9, 12.5e9, 13.75e9, 15.0e9])
def test_lossless_grounded_slab_is_totally_reflective(f):
    g = green(LOSSLESS, f)
    e, _ = g.stack_fields(0.0, 0.0, np.array([1.0, 0.0], dtype=complex), 0.0, port="top")
    gamma = e[0] - 1.0
    assert abs(gamma) == pytest.approx(1.0, abs=1e-9)


def test_lossy_grounded_slab_reflects_less_than_unity():
    g = green(LOSSY, F)
    e, _ = g.stack_fields(0.0, 0.0, np.array([1.0, 0.0], dtype=complex), 0.0, port="top")
    gamma = e[0] - 1.0
    assert abs(gamma) < 1.0


def test_vanishing_thickness_shorts_the_top_interface():
    thin = Substrate("x", 2.2, 0.0, 1e-9)
    g = green(thin, F)
    e, _ = g.stack_fields(0.0, 0.0, np.array([1.0, 0.0], dtype=complex), 0.0, port="top")
    gamma = e[0] - 1.0
    assert abs(gamma - (-1.0)) < 1e-6


# ------------------------------------------------------------------------- ground-plane physics


def test_bottom_port_does_not_exist():
    g = green(LOSSLESS, F)
    with pytest.raises(ValueError, match="bottom"):
        g.stack_fields(0.0, 0.0, np.array([1.0, 0.0], dtype=complex), -H, port="bottom")


@pytest.mark.parametrize("kt", [0.0, 100.0, 500.0, 2000.0])
def test_tangential_e_vanishes_at_the_ground(kt):
    g = green(LOSSLESS, F)
    j = np.array([1.0 + 0.2j, -0.3j], dtype=complex)
    e, _ = g.sheet_fields(kt, 0.0, j, -H, sheet_side="above")
    assert np.max(np.abs(e)) < 1e-9


@pytest.mark.parametrize("z", [-H - 1e-6, -H - 1.0, -2 * H])
def test_sheet_fields_vanish_below_the_ground(z):
    g = green(LOSSLESS, F)
    j = np.array([1.0, -0.5j], dtype=complex)
    e, h = g.sheet_fields(0.0, 0.0, j, z)
    assert np.allclose(e, 0.0) and np.allclose(h, 0.0)


@pytest.mark.parametrize("z", [-H - 1e-6, -H - 1.0])
def test_stack_fields_vanish_below_the_ground(z):
    g = green(LOSSLESS, F)
    e, h = g.stack_fields(0.0, 0.0, np.array([1.0, 0.0], dtype=complex), z, port="top")
    assert np.allclose(e, 0.0) and np.allclose(h, 0.0)


def test_grounded_and_air_backed_impedances_differ():
    g_air = LayeredGreenFunction(LOSSLESS, F, backing="air")
    g_pec = green(LOSSLESS, F)
    z_air, _ = g_air.equivalent_impedances(np.array([100.0]), np.array([0.0]))
    z_pec, _ = g_pec.equivalent_impedances(np.array([100.0]), np.array([0.0]))
    assert not np.isclose(z_air[0], z_pec[0])


# ------------------------------------------------------------------------------------ passivity


def test_passivity_holds_for_lossy_grounded_stack_across_kt():
    g = green(LOSSY, F)
    kt = np.array([0.0, 50.0, 200.0, 209.0, 300.0, 500.0, 1000.0, 5000.0])
    z_te, z_tm = g.equivalent_impedances(kt, np.zeros_like(kt))
    assert np.all(z_te.real >= -1e-9)
    assert np.all(z_tm.real >= -1e-9)


# ---------------------------------------------------------------- independent-reference agreement

SHEET_CASES = [
    (0.0, 0.0, (1.0 + 0.3j, -0.2j), F, 2.2, 0.0009, H, "normal incidence, lossy"),
    (150.0, -80.0, (0.5 - 0.1j, 1.2 + 0.4j), F, 2.2, 0.0009, H, "oblique, lossy"),
    (300.0, 0.0, (1.0, 0.0), 10.0e9, 4.4, 0.02, 0.8e-3, "evanescent-ish, different eps/h"),
    (0.0, 0.0, (1.0, 1.0j), 15.0e9, 2.2, 0.0, H, "lossless"),
]


@pytest.mark.parametrize("kx,ky,j,f,eps_r,tand,h,label", SHEET_CASES, ids=[c[-1] for c in SHEET_CASES])
def test_sheet_fields_match_independent_reference(kx, ky, j, f, eps_r, tand, h, label):
    sub = Substrate("x", eps_r, tand, h)
    g = LayeredGreenFunction(sub, f, backing="pec")
    ref = solve_current_sheet_grounded(kx, ky, list(j), f, sub.eps_complex, h)
    for z in (0.5 * h, 0.1 * h, 0.0, -0.3 * h, -0.7 * h, -h, -h - 0.2 * h, -h - 1.0):
        side = "above" if z >= 0 else "below"
        e_mine, h_mine = g.sheet_fields(kx, ky, np.array(j, dtype=complex), z, sheet_side=side)
        e_ref, h_ref = ref["fields"](z, sheet_side=side)
        assert np.max(np.abs(e_mine - e_ref)) < 1e-9
        assert np.max(np.abs(h_mine - h_ref)) < 1e-9


STACK_CASES = [
    (0.0, 0.0, F, 2.2, 0.0009, H, "normal incidence, lossy"),
    (20.0, 30.0, F, 2.2, 0.0009, H, "oblique 20/30 deg, lossy"),
    (0.0, 0.0, 15.0e9, 2.2, 0.0, H, "normal incidence, lossless"),
]


@pytest.mark.parametrize("theta_deg,phi_deg,f,eps_r,tand,h,label", STACK_CASES, ids=[c[-1] for c in STACK_CASES])
def test_stack_fields_match_independent_reference(theta_deg, phi_deg, f, eps_r, tand, h, label):
    sub = Substrate("x", eps_r, tand, h)
    k0 = 2 * np.pi * f / 299792458.0
    th, ph = np.radians(theta_deg), np.radians(phi_deg)
    kx, ky = k0 * np.sin(th) * np.cos(ph), k0 * np.sin(th) * np.sin(ph)
    kt = np.hypot(kx, ky)
    if kt == 0:
        e_inc = np.array([1.0, 0.0, 0.0])
    else:
        p = np.array([kx, ky, 0.0]) / kt
        e_inc = np.cross(np.array([0.0, 0.0, 1.0]), p)  # TE polarization
    ref = solve_incident_fields_grounded(kx, ky, e_inc, f, sub.eps_complex, h)
    g = LayeredGreenFunction(sub, f, backing="pec")
    for z in (0.5 * h, 0.0, -0.3 * h, -h, -h - 0.5):
        e_mine, h_mine = g.stack_fields(kx, ky, e_inc[:2].astype(complex), z, port="top")
        e_ref, h_ref = ref(z)
        assert np.max(np.abs(e_mine - e_ref)) < 1e-9
        assert np.max(np.abs(h_mine - h_ref)) < 1e-9


def test_independent_reference_sheet_condition_at_grounded_stack():
    """z_hat x [H(0+) - H(0-)] = J still holds with a PEC ground below (unaffected by backing)."""
    j = (0.7 - 0.1j, -0.4 + 0.2j)
    ref = solve_current_sheet_grounded(80.0, -30.0, list(j), F, LOSSY.eps_complex, H)
    _, h_above = ref["fields"](0.0, sheet_side="above")
    _, h_below = ref["fields"](0.0, sheet_side="below")
    dh = h_above - h_below
    # z_hat x dH = (-dHy, dHx) = J
    zxdh = np.array([-dh[1], dh[0]])
    assert np.max(np.abs(zxdh - np.array(j))) < 1e-9
