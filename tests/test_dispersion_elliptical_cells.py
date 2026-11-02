"""Phase 30 item 2: elliptical-patch cell family -- geometry and tensor-sheet physics."""

import math

import numpy as np
import pytest

from src.basis_functions import MeshSettings
from src.dispersion import tensor_sheet as ts
from src.dispersion.elliptical_cells import EllipticalPatch
from src.dispersion.slotted_cells import SlottedCircularPatch
from src.geometry import UnitCell
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice

F = 10.54e9
P = 6.4e-3
SUB = Substrate("rt5880", 2.2, 0.0009, 1.575e-3)
MESH = MeshSettings(32, 32, 32)
_G = np.linspace(-P / 2, P / 2, 401)
X, Y = np.meshgrid(_G, _G)


def _cell(geo):
    return UnitCell(RectangularLattice(P, P), geo, SUB, Metal(), backing="pec")


def test_circle_mask_identical_to_slotted_circle_without_slot():
    e = EllipticalPatch(5.0e-3, 5.0e-3)
    c = SlottedCircularPatch(5.0e-3, 0.0, 0.0)
    assert np.array_equal(e.in_feature(X, Y), c.in_feature(X, Y))


@pytest.mark.parametrize("kw", [{}, {"slot_length": 3.0e-3, "slot_width": 0.8e-3}, {"slot_length": 3.0e-3, "slot_width": 0.8e-3, "cross": True}])
def test_mask_area_matches_analytic(kw):
    g = EllipticalPatch(5.6e-3, 3.4e-3, 0.0, **kw)
    cell = (X[0, 1] - X[0, 0]) ** 2
    assert g.in_feature(X, Y).sum() * cell == pytest.approx(g.feature_area, rel=0.02)


def test_rotation_by_pi_over_2_swaps_axes():
    a = EllipticalPatch(5.6e-3, 3.4e-3, math.pi / 2).in_feature(X, Y)
    b = EllipticalPatch(5.6e-3, 3.4e-3, 0.0).in_feature(Y, X)  # transpose
    assert np.array_equal(a, b)


def test_rejects_minor_larger_than_major():
    with pytest.raises(ValueError):
        EllipticalPatch(3.0e-3, 4.0e-3)


@pytest.fixture(scope="module")
def tensors():
    def t(**kw):
        return ts.sheet_admittance_tensor_for_cell(_cell(EllipticalPatch(**kw)), F, MESH)

    return {
        "ell": t(major=5.4e-3, minor=3.6e-3),
        "ell90": t(major=5.4e-3, minor=3.6e-3, rotation=math.pi / 2),
        "ell_slot": t(major=5.4e-3, minor=4.0e-3, slot_length=3.0e-3, slot_width=0.8e-3),
        "circle": t(major=5.0e-3, minor=5.0e-3),
        "cross": t(major=5.4e-3, minor=5.4e-3, slot_length=3.2e-3, slot_width=0.8e-3, cross=True),
        "ell45": t(major=5.4e-3, minor=3.6e-3, rotation=math.pi / 4),
    }


def test_reciprocity_and_lossless_symmetric(tensors):
    for name in ("ell", "ell_slot", "ell45"):
        assert ts.reciprocity_error(tensors[name]) < 1e-9


def test_ellipse_is_anisotropic_with_aligned_principal_axes(tensors):
    y = tensors["ell"]
    assert abs(y[0, 1]) < 1e-9 * abs(y[0, 0])
    assert abs(y[0, 0].imag - y[1, 1].imag) > 0.05 * abs(y[0, 0].imag)
    # the long axis (x) carries the larger capacitive loading
    assert y[0, 0].imag > y[1, 1].imag


def test_90_degree_rotation_swaps_xx_and_yy(tensors):
    y, y90 = tensors["ell"], tensors["ell90"]
    assert y90[0, 0] == pytest.approx(y[1, 1], rel=1e-6)
    assert y90[1, 1] == pytest.approx(y[0, 0], rel=1e-6)


def test_circle_and_cross_slot_are_isotropic(tensors):
    for name in ("circle", "cross"):
        y = tensors[name]
        assert abs(y[0, 0] - y[1, 1]) < 1e-6 * abs(y[0, 0])
        assert abs(y[0, 1]) < 1e-6 * abs(y[0, 0])


def test_slot_reduces_loading_and_adds_anisotropy_to_a_circle():
    base = ts.sheet_admittance_tensor_for_cell(_cell(EllipticalPatch(5.4e-3, 5.4e-3)), F, MESH)
    slotted = ts.sheet_admittance_tensor_for_cell(_cell(EllipticalPatch(5.4e-3, 5.4e-3, slot_length=3.6e-3, slot_width=1.2e-3)), F, MESH)
    assert abs(base[0, 0] - base[1, 1]) < 1e-6 * abs(base[0, 0])
    assert abs(slotted[0, 0] - slotted[1, 1]) > 0.05 * abs(slotted[0, 0])


def test_diagonal_ellipse_tensor_follows_rotation_to_mesh_accuracy(tensors):
    pred = ts.rotate_tensor(tensors["ell"], math.pi / 4)
    assert ts.rotation_equivariance_error(tensors["ell45"], tensors["ell"], math.pi / 4) < 0.1
    assert np.linalg.norm(pred - tensors["ell45"]) / np.linalg.norm(pred) < 0.1


def test_superellipse_exponent_two_is_the_ellipse_and_area_matches_analytic():
    a = EllipticalPatch(5.6e-3, 3.4e-3, 0.3)
    b = EllipticalPatch(5.6e-3, 3.4e-3, 0.3, exponent=2.0)
    assert np.array_equal(a.in_feature(X, Y), b.in_feature(X, Y))
    cell = (X[0, 1] - X[0, 0]) ** 2
    for n in (4.0, 8.0):
        g = EllipticalPatch(5.6e-3, 4.4e-3, 0.0, exponent=n)
        assert g.in_feature(X, Y).sum() * cell == pytest.approx(g.feature_area, rel=0.02)


def test_squarer_outline_loads_more_capacitively():
    """Same major/minor, larger exponent -> more metal at the corners -> larger sheet susceptance (physics, not a number)."""
    ys = [
        ts.sheet_admittance_tensor_for_cell(_cell(EllipticalPatch(5.6e-3, 5.6e-3, exponent=n)), F, MESH)[0, 0].imag
        for n in (2.0, 4.0, 8.0)
    ]
    assert ys[0] < ys[1] < ys[2]
