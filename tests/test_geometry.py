import math

import numpy as np
import pytest

from src.geometry import Geometry, UnitCell
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice

SUB = Substrate("t", 2.2, 0.0009, 1.575e-3)


def cell(kind="rectangle", length=5e-3, width=2e-3, px=8e-3, py=8e-3, **kw):
    return UnitCell(RectangularLattice(px, py), Geometry(kind, length, width, **kw), SUB, Metal())


def grid_fill_fraction(c, n=1200):
    px, py = c.lattice.period_x, c.lattice.period_y
    x = (np.arange(n) + 0.5) / n * px - px / 2
    y = (np.arange(n) + 0.5) / n * py - py / 2
    xx, yy = np.meshgrid(x, y, indexing="ij")
    return c.is_metal(xx, yy).mean()


# ------------------------------------------------------------- construction


def test_kinds_and_polarity():
    assert Geometry("rectangle", 5e-3, 2e-3).polarity == "metal"
    assert Geometry("patch", 3e-3, 3e-3).polarity == "metal"
    assert Geometry("dipole", 5e-3, 0.5e-3).polarity == "metal"
    assert Geometry("slot", 5e-3, 0.5e-3).polarity == "aperture"


@pytest.mark.parametrize(
    "kw",
    [
        dict(kind="circle", length=1e-3, width=1e-3),
        dict(kind="rectangle", length=0.0, width=1e-3),
        dict(kind="rectangle", length=1e-3, width=-1e-3),
        dict(kind="rectangle", length=float("nan"), width=1e-3),
        dict(kind="rectangle", length=1e-3, width=1e-3, rotation=float("inf")),
        dict(kind="dipole", length=1e-3, width=2e-3),  # dipole must be elongated
        dict(kind="slot", length=1e-3, width=2e-3),  # slot must be elongated
        dict(kind="rectangle", length="5e-3", width=1e-3),
    ],
)
def test_geometry_rejects_invalid(kw):
    with pytest.raises(ValueError):
        Geometry(**kw)


# ------------------------------------------------------------ shape maths


def test_corners_unrotated():
    c = Geometry("rectangle", 5e-3, 2e-3).corners()
    assert set(map(tuple, np.round(c * 1e3, 9))) == {(2.5, 1.0), (-2.5, 1.0), (-2.5, -1.0), (2.5, -1.0)}


def test_corners_rotated_90_swap_extents():
    c = Geometry("rectangle", 5e-3, 2e-3, rotation=math.pi / 2).corners()
    assert np.max(np.abs(c[:, 0])) == pytest.approx(1e-3, abs=1e-15)
    assert np.max(np.abs(c[:, 1])) == pytest.approx(2.5e-3, abs=1e-15)


@pytest.mark.parametrize("rot_deg", [0, 17, 45, 90, 133, 270])
def test_rotation_preserves_edge_lengths_and_centroid(rot_deg):
    g = Geometry("rectangle", 5e-3, 2e-3, center_x=0.3e-3, center_y=-0.2e-3, rotation=math.radians(rot_deg))
    c = g.corners()
    edges = np.linalg.norm(np.roll(c, -1, axis=0) - c, axis=1)
    assert edges == pytest.approx([5e-3, 2e-3, 5e-3, 2e-3], rel=1e-12)
    assert c.mean(axis=0) == pytest.approx([0.3e-3, -0.2e-3], abs=1e-15)


def test_in_feature_basic_points():
    g = Geometry("rectangle", 5e-3, 2e-3)
    assert g.in_feature(0.0, 0.0)
    assert g.in_feature(2.4e-3, 0.9e-3)
    assert not g.in_feature(2.6e-3, 0.0)
    assert not g.in_feature(0.0, 1.1e-3)
    r = Geometry("rectangle", 5e-3, 2e-3, rotation=math.pi / 2)
    assert r.in_feature(0.0, 2.4e-3) and not r.in_feature(2.4e-3, 0.0)


def test_rotated_90_equals_swapped_dimensions():
    a = cell("rectangle", 5e-3, 2e-3, rotation=math.pi / 2)
    b = cell("rectangle", 2e-3, 5e-3)
    n = 400
    x = np.linspace(-4e-3, 4e-3, n) + 1.234e-6  # avoid exact edge points
    xx, yy = np.meshgrid(x, x, indexing="ij")
    assert np.array_equal(a.is_metal(xx, yy), b.is_metal(xx, yy))


def test_centered_unrotated_symmetry():
    c = cell("rectangle", 5e-3, 2e-3)
    rng = np.random.default_rng(3)
    x = rng.uniform(-4e-3, 4e-3, 2000)
    y = rng.uniform(-4e-3, 4e-3, 2000)
    m = c.is_metal(x, y)
    assert np.array_equal(m, c.is_metal(-x, y))
    assert np.array_equal(m, c.is_metal(x, -y))


# ---------------------------------------------------------- cell boundaries


def test_feature_inside_cell_ok():
    cell("rectangle", 5e-3, 2e-3)
    cell("rectangle", 5e-3, 2e-3, rotation=math.radians(30))  # bbox 5.33 x 4.0 mm


def test_feature_touching_cell_edge_allowed():
    c = cell("dipole", 8e-3, 1e-3)  # spans full period: continuous wire
    assert c.metal_area == pytest.approx(8e-6, rel=1e-14)
    cell("rectangle", 8e-3, 8e-3)  # fills cell exactly
    cell("rectangle", 8e-3, 1e-3, rotation=math.pi)  # rotation rounding noise
    cell("rectangle", 8e-3, 1e-3, px=1e-3, py=8e-3, rotation=math.pi / 2)  # cos(pi/2) ~ 6e-17


@pytest.mark.parametrize(
    "kw",
    [
        dict(length=8.001e-3, width=1e-3),  # too long
        dict(length=5e-3, width=2e-3, center_x=3e-3),  # shifted out
        dict(length=5e-3, width=2e-3, center_y=-3.5e-3),
        dict(length=8e-3, width=1e-3, rotation=math.radians(10)),  # 4cos10+0.5sin10 = 4.03 mm > 4 mm
        dict(length=5e-3, width=2e-3, rotation=math.pi / 2, center_x=3.5e-3),
    ],
)
def test_feature_outside_cell_rejected(kw):
    with pytest.raises(ValueError, match="outside the unit cell"):
        cell("rectangle", **kw)


def test_non_square_lattice_bounds():
    cell("rectangle", 5e-3, 2e-3, px=6e-3, py=3e-3)
    with pytest.raises(ValueError):
        cell("rectangle", 5e-3, 4e-3, px=6e-3, py=3e-3)  # width > Py
    with pytest.raises(ValueError):
        cell("rectangle", 5e-3, 2e-3, px=6e-3, py=3e-3, rotation=math.pi / 2)


def test_slot_covering_whole_cell_leaves_no_metal():
    with pytest.raises(ValueError, match="no metal"):
        cell("slot", 8e-3, 8e-3)


# ------------------------------------------------------------ area / fill


def test_metal_area_rectangle_and_slot_are_complementary():
    r = cell("rectangle", 5e-3, 2e-3)
    s = cell("slot", 5e-3, 2e-3)
    assert r.metal_area == pytest.approx(10e-6, rel=1e-14)
    assert s.metal_area == pytest.approx(64e-6 - 10e-6, rel=1e-14)
    assert r.metal_area + s.metal_area == pytest.approx(r.lattice.cell_area, rel=1e-14)
    assert r.metal_fill_fraction == pytest.approx(10 / 64, rel=1e-14)


def test_slot_is_exact_complement_of_rectangle_pointwise():
    kw = dict(length=4e-3, width=1e-3, center_x=0.5e-3, rotation=math.radians(20))
    r, s = cell("rectangle", **kw), cell("slot", **kw)
    rng = np.random.default_rng(5)
    x = rng.uniform(-4e-3, 4e-3, 4000)
    y = rng.uniform(-4e-3, 4e-3, 4000)
    assert np.array_equal(s.is_metal(x, y), ~r.is_metal(x, y))


@pytest.mark.parametrize(
    "kind, kw",
    [
        ("rectangle", dict(length=5e-3, width=2e-3)),
        ("patch", dict(length=4e-3, width=4e-3, rotation=math.radians(30))),
        ("dipole", dict(length=6e-3, width=0.7e-3, rotation=math.radians(63), center_x=0.4e-3)),
        ("slot", dict(length=5e-3, width=1e-3, rotation=math.radians(45), center_y=-0.5e-3)),
    ],
)
def test_sampled_fill_fraction_matches_analytic_area(kind, kw):
    """is_metal() (pointwise) must integrate to the closed-form area."""
    c = cell(kind, **kw)
    assert grid_fill_fraction(c) == pytest.approx(c.metal_fill_fraction, abs=3e-3)
