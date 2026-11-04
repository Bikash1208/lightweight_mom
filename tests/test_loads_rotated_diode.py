"""Phase 28 rotated-load fix: a measured-diode load on a ROTATED axis must switch the cell (it used to freeze to an open circuit, ON == OFF).
Physics checks on a split patch (two plates + bridge carrying the MA4AGP907 load) turned by an angle: ON != OFF, rotation covariance of the
sheet tensor (principal axis, principal values), reciprocity, mirror symmetry, axis-aligned results unchanged."""

import math
import warnings
from dataclasses import replace

import numpy as np
import pytest

from src.composite_geometry import AxisRect, CompositeGeometry
from src.geometry import UnitCell
from src.hologram.coded_cell_eval import sheet_tensor
from src.loads import LumpedLoad, check_load_geometry, default_transverse_conducts, load_matrix
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice
from src.polygon_geometry import PolyComposite, oriented_rect
from src.ris import diodes_measured as dm
from src.basis_functions import MeshSettings, RooftopBasis

P = 6.4e-3
SUB = Substrate("rt5880", 2.2, 0.0009, 1.575e-3)
F0 = 10.54e9
ST = {"OFF": "OFF_0V", "ON": "ON_10mA"}
A, H, G, WD = 4.0e-3, 2.4e-3, 0.4e-3, 1.2e-3


def _cell(state, ang, tc=None, with_load=True):
    t = math.radians(ang)
    u = np.array([math.cos(t), math.sin(t)])
    off = (G + (A - G) / 2) / 2
    halves = tuple(oriented_rect(s * off * u[0], s * off * u[1], (A - G) / 2, H, t) for s in (-1, 1))
    geo = PolyComposite(metal=halves + (oriented_rect(0.0, 0.0, G, WD, t),))
    direction = "x" if ang % 180 == 0 else "y" if ang % 180 == 90 else float(ang)
    loads = (dm.load(ST[state], 0.0, 0.0, G, WD, direction, None, tc),) if with_load else ()
    return UnitCell(RectangularLattice(P, P), geo, SUB, Metal(), backing="pec", loads=loads)


def _y(state, ang, **kw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return sheet_tensor(_cell(state, ang, **kw), F0, 32)


def _principal(y):
    w, v = np.linalg.eig(y)
    i = int(np.argmax(w.imag))
    return w[i], math.degrees(math.atan2(v[1, i].real, v[0, i].real)) % 180.0


def test_default_transverse_rule():
    assert [default_transverse_conducts(d) for d in ("x", "y", 0, 90, 180, 45, 30.0)] == [False, False, False, False, False, True, True]
    assert dm.load("ON_10mA", 0, 0, 1e-3, 1e-3, "x").transverse_conducts is False
    assert dm.load("ON_10mA", 0, 0, 1e-3, 1e-3, 90).transverse_conducts is False
    assert dm.load("ON_10mA", 0, 0, 1e-3, 1e-3, 45).transverse_conducts is True
    assert dm.load("ON_10mA", 0, 0, 1e-3, 1e-3, 45, None, False).transverse_conducts is False  # explicit choice still honoured


def test_axis_aligned_load_matrix_is_unchanged():
    old = LumpedLoad(0.0, 0.0, G, WD, "x", z_of_f=dm.z_of_f("ON_10mA", None), transverse_conducts=False)
    new = dm.load("ON_10mA", 0.0, 0.0, G, WD, "x")
    assert (new.direction, new.transverse_conducts, new.length, new.width) == (old.direction, old.transverse_conducts, old.length, old.width)
    b = RooftopBasis(_cell("ON", 0), 32, 32)
    w = 2 * math.pi * F0
    assert np.array_equal(load_matrix(b, [old], w), load_matrix(b, [new], w))


def test_blocked_rotated_load_is_warned():
    c = _cell("ON", 45, tc=False)
    with pytest.warns(UserWarning, match="freezes the region"):
        check_load_geometry(c.loads, c, 32, 32)


@pytest.mark.parametrize("ang", [30, 45, 60])
def test_rotated_diode_switches(ang):
    off, on = _y("OFF", ang), _y("ON", ang)
    d_on, _ = _principal(on)
    d_off, _ = _principal(off)
    assert d_on.imag - d_off.imag > 0.5e-3  # axis-aligned: ~0.75 mS here
    # the old behaviour (blocked) is the regression: ON == OFF
    assert np.allclose(_y("ON", ang, tc=False), _y("OFF", ang, tc=False), rtol=1e-6, atol=1e-9)


def _uu(y, ang):
    t = math.radians(ang)
    u = np.array([math.cos(t), math.sin(t)])
    return u @ y @ u


def test_rotation_covariance_of_the_loaded_tensor():
    """Component along the diode axis and the ON-OFF contrast are rotation covariant to the staircase error of a rotated patch."""
    ref_on, ref_off = _uu(_y("ON", 0), 0), _uu(_y("OFF", 0), 0)
    for ang in (45, 90, 135):
        on, off = _y("ON", ang), _y("OFF", ang)
        assert abs(_uu(on, ang) - ref_on) / abs(ref_on) < 0.10, (ang, _uu(on, ang), ref_on)
        assert abs((_uu(on, ang) - _uu(off, ang)) - (ref_on - ref_off)) / abs(ref_on - ref_off) < 0.10
        assert abs((_principal(on)[1] - ang + 90) % 180 - 90) < 1.0  # principal axis follows the diode axis
    y0, y90 = _y("ON", 0), _y("ON", 90)  # 0 vs 90 deg: exact axis swap
    assert np.allclose(y0, y90[::-1, ::-1], rtol=1e-9, atol=1e-12)


def test_unloaded_rotation_invariance():
    """The patch with the bridge as plain metal (no load): 90 deg is an exact axis swap, 45 deg keeps the principal values to the staircase error."""
    y0, y90, y45 = (_y("ON", a, with_load=False) for a in (0, 90, 45))
    assert np.allclose(y0, y90[::-1, ::-1], rtol=1e-9, atol=1e-12)
    w0 = sorted(np.linalg.eigvals(y0).imag)
    w45 = sorted(np.linalg.eigvals(y45).imag)
    assert np.allclose(w0, w45, rtol=0.2)


def test_reciprocity_and_mirror_symmetry():
    y = _y("ON", 30)
    assert abs(y[0, 1] - y[1, 0]) < 1e-9 * abs(y).max()
    ya, yb = _y("ON", 30), _y("ON", 150)  # mirror image about the x axis: y_xy flips sign, the diagonal is unchanged
    assert np.allclose(np.diag(ya), np.diag(yb), rtol=1e-6)
    assert np.allclose(ya[0, 1], -yb[0, 1], rtol=1e-6, atol=1e-9)


def test_near_axis_rotated_path_agrees_with_axis_aligned():
    on0 = _principal(_y("ON", 0))[0]
    near = _principal(_y("ON", 0.01, tc=True))[0]
    assert abs(near - on0) / abs(on0) < 0.05
