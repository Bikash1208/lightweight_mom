"""Phase 20 R3 item 2: polygon-composite geometry (octagon, oriented rectangles) -- masks, area, alignment, bit-identity with CompositeGeometry."""

import math
import warnings

import numpy as np
import pytest

from src.basis_functions import MeshSettings
from src.composite_geometry import AxisRect, CompositeGeometry
from src.excitation import ReferencePlanes
from src.geometry import UnitCell
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice
from src.polygon_geometry import Poly, PolyComposite, oriented_rect, regular_polygon
from src.sparameters import compute_s_parameters_oblique

SUB = Substrate("rt5880", 2.2, 0.0009, 1.575e-3)


def test_poly_contains_and_area():
    sq = Poly(((0, 0), (2, 0), (2, 2), (0, 2)))
    assert sq.area == 4 and sq.is_manhattan() and sq.bbox() == (0, 2, 0, 2)
    x, y = np.meshgrid([-0.1, 0.0, 1.0, 2.0, 2.1], [-0.1, 0.0, 1.0, 2.0, 2.1], indexing="ij")
    got = sq.contains(x, y)
    want = (x >= 0) & (x <= 2) & (y >= 0) & (y <= 2)
    assert np.array_equal(got, want)  # boundary inclusive
    ccw = Poly(((0, 0), (2, 0), (2, 2), (0, 2)))
    cw = Poly(tuple(reversed(ccw.vertices)))
    assert np.array_equal(ccw.contains(x, y), cw.contains(x, y)) and cw.area == 4
    with pytest.raises(ValueError):
        Poly(((0, 0), (1, 1)))
    with pytest.raises(ValueError):
        Poly(((0, 0), (1, 1), (2, 2)))  # zero area


def test_regular_octagon_and_oriented_rect():
    r = 4e-3
    octo = regular_polygon(0.0, 0.0, r, 8, rotation=math.pi / 8)
    assert octo.area == pytest.approx(2 * math.sqrt(2) * r * r, rel=1e-12)
    assert not octo.is_manhattan()
    a = r * math.cos(math.pi / 8)  # apothem: flat edges face 0, 45, 90, ... deg
    assert octo.contains(a - 1e-9, 0.0) and not octo.contains(a + 1e-6, 0.0)
    rect = oriented_rect(1e-3, 2e-3, 3e-3, 1e-3, math.radians(45))
    assert rect.area == pytest.approx(3e-6, rel=1e-12) and rect.contains(1e-3, 2e-3)
    c = np.array(rect.vertices).mean(axis=0)
    assert c == pytest.approx([1e-3, 2e-3])
    with pytest.raises(ValueError):
        oriented_rect(0, 0, 0.0, 1e-3, 0.3)


def test_polycomposite_interface_area_and_alignment():
    octo = regular_polygon(0.0, 0.0, 4e-3, 8, rotation=math.pi / 8)
    g = PolyComposite((octo,))
    assert g.kind == "poly_composite" and g.polarity == "metal"
    assert g.feature_area == pytest.approx(octo.area, rel=3e-3)
    x0, x1, y0, y1 = g.bbox()
    assert x1 - x0 == pytest.approx(g.length) and y1 - y0 == pytest.approx(g.width)
    assert g.corners().shape == (4, 2)
    assert g.non_manhattan_edges() and not g.edges_aligned(RectangularLattice(12e-3, 12e-3), MeshSettings(48, 48, 25))
    hole = PolyComposite((octo,), apertures=(AxisRect(0.0, 0.0, 1e-3, 1e-3),))
    assert not hole.in_feature(0.0, 0.0) and hole.in_feature(2e-3, 0.0)
    with pytest.raises(ValueError):
        PolyComposite(())
    with pytest.raises(ValueError):
        PolyComposite(("not a shape",))


def test_axis_aligned_poly_rectangle_is_bit_identical_to_the_composite_rectangle():
    rect = AxisRect(0.0, 0.0, 6e-3, 4e-3)
    poly = Poly(((-3e-3, -2e-3), (3e-3, -2e-3), (3e-3, 2e-3), (-3e-3, 2e-3)))
    lat = RectangularLattice(12e-3, 12e-3)
    mesh = MeshSettings(48, 48, 25)
    uc_a = UnitCell(lat, CompositeGeometry((rect,)), SUB, Metal(), backing="pec")
    uc_b = UnitCell(lat, PolyComposite((poly,)), SUB, Metal(), backing="pec")
    assert PolyComposite((poly,)).edges_aligned(lat, mesh) and CompositeGeometry((rect,)).edges_aligned(lat, mesh)
    planes = ReferencePlanes(0.0, None, SUB.thickness)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a = compute_s_parameters_oblique(uc_a, 10e9, mesh, planes, math.radians(20), 0.3)
        b = compute_s_parameters_oblique(uc_b, 10e9, mesh, planes, math.radians(20), 0.3)
    assert np.array_equal(a.s11, b.s11)


def test_octagon_with_fourfold_symmetry_gives_equal_xy_response_and_conserves_energy():
    lat = RectangularLattice(13.5e-3, 13.5e-3)
    uc = UnitCell(
        lat,
        PolyComposite((regular_polygon(0.0, 0.0, 4.5e-3, 8, rotation=math.pi / 8),)),
        Substrate("s", 2.2, 0.0, 1.575e-3),
        Metal(),
        backing="pec",
    )
    planes = ReferencePlanes(0.0, None, 1.575e-3)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s = compute_s_parameters_oblique(uc, 10e9, MeshSettings(54, 54, 28), planes, 0.0, 0.0)
    assert abs(s.s11[0, 0] - s.s11[1, 1]) < 1e-9  # TE (y) == TM (x): the octagon has the symmetry of the square lattice
    assert abs(s.s11[0, 1]) < 1e-9
    assert np.allclose(s.power_out_per_incident(), 1.0, atol=1e-9)  # lossless
