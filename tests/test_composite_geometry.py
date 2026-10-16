"""Phase 20 R1 item 3: composite (union minus apertures) geometry -- identity with the rectangle, area/mask, alignment, STEP."""

import math
import os
import warnings

import numpy as np
import pytest

from src.basis_functions import MeshSettings, RooftopBasis
from src.composite_geometry import AxisRect, CompositeGeometry
from src.discretization import geometry_edges_aligned
from src.excitation import ReferencePlanes
from src.geometry import Geometry, UnitCell
from src.loads import LumpedLoad
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice
from src.sparameters import compute_s_parameters, compute_s_parameters_oblique

F = 10e9
SUB = Substrate("rt5880", 2.2, 0.0009, 1.575e-3)
LAT = RectangularLattice(10e-3, 8e-3)


def _cell(geo, backing="pec", loads=()):
    return UnitCell(LAT, geo, SUB, Metal(), backing=backing, loads=loads)


def _planes(backing):
    return ReferencePlanes(2e-3, None if backing == "pec" else 2e-3, SUB.thickness)


@pytest.mark.parametrize("backing", ["pec", "air"])
def test_composite_of_one_rectangle_is_bit_identical_to_the_rectangle_geometry(backing):
    rect = Geometry(kind="rectangle", length=6e-3, width=4.5e-3, center_x=0.5e-3, center_y=-0.25e-3)
    comp = CompositeGeometry((AxisRect(0.5e-3, -0.25e-3, 6e-3, 4.5e-3),))
    ca, cb = _cell(rect, backing), _cell(comp, backing)
    mesh = MeshSettings(16, 16, 10)
    assert np.array_equal(RooftopBasis(ca, 16, 16).metal_mask, RooftopBasis(cb, 16, 16).metal_mask)
    assert comp.feature_area == pytest.approx(rect.feature_area, rel=1e-14)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a = compute_s_parameters(ca, F, mesh, _planes(backing))
        b = compute_s_parameters(cb, F, mesh, _planes(backing))
        ao = compute_s_parameters_oblique(ca, F, mesh, _planes(backing), math.radians(35.0), 0.4)
        bo = compute_s_parameters_oblique(cb, F, mesh, _planes(backing), math.radians(35.0), 0.4)
    assert np.array_equal(a.s11, b.s11) and np.array_equal(ao.s11, bo.s11)


def test_area_and_mask_of_a_union_with_overlap_and_apertures():
    comp = CompositeGeometry(
        metal=(
            AxisRect(-2e-3, 0, 4e-3, 4e-3),
            AxisRect(2e-3, 0, 4e-3, 4e-3),
            AxisRect(0, 0, 4e-3, 2e-3),
        ),  # two squares + an overlapping bridge
        apertures=(AxisRect(0, 0, 1e-3, 1e-3),),
    )
    # exact: two 4x4 squares touch at x=0 (area 32e-6) , the bridge lies inside them (no extra area), minus the 1x1 aperture
    assert comp.feature_area == pytest.approx(32e-6 - 1e-6, rel=1e-12)
    x = np.linspace(-4e-3, 4e-3, 801)
    gx, gy = np.meshgrid(x, x, indexing="ij")
    frac = comp.in_feature(gx, gy).mean() * (8e-3) ** 2
    assert frac == pytest.approx(comp.feature_area, rel=3e-3)
    assert not comp.in_feature(0.0, 0.0) and comp.in_feature(1e-3, 0.0) and not comp.in_feature(5e-3, 0.0)
    assert np.allclose(comp.corners().max(axis=0), [4e-3, 2e-3]) and np.allclose(
        comp.corners().min(axis=0), [-4e-3, -2e-3]
    )
    with pytest.raises(ValueError):
        CompositeGeometry(())
    with pytest.raises(ValueError):
        AxisRect(0, 0, -1e-3, 1e-3)


def test_gapped_patch_with_arms_of_different_lengths_and_alignment_warning():
    arms = (
        AxisRect(-2.5e-3, 0, 3e-3, 6e-3),
        AxisRect(2e-3, 0, 4e-3, 6e-3),
    )  # left arm x in [-4, -1], right arm [0, 4]: gap [-1, 0]
    comp = CompositeGeometry(arms)
    cell = _cell(comp)
    assert comp.feature_area == pytest.approx((3e-3 + 4e-3) * 6e-3)
    mesh = MeshSettings(20, 16, 10)  # dx = 0.5 mm, dy = 0.5 mm: all edges on mesh lines
    assert geometry_edges_aligned(comp, LAT, mesh) and comp.edges_aligned(LAT, mesh)
    off = CompositeGeometry((AxisRect(0, 0, 6.3e-3, 4e-3),))
    assert not geometry_edges_aligned(off, LAT, mesh)
    with pytest.warns(UserWarning, match="off the 20 x 16 mesh grid"):
        compute_s_parameters_oblique(_cell(off), F, mesh, _planes("pec"), 0.0, 0.0)
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # an aligned composite raises no warning of this kind
        warnings.filterwarnings("ignore", message=".*Floquet harmonic.*")
        compute_s_parameters_oblique(cell, F, mesh, _planes("pec"), 0.0, 0.0)


def test_loaded_gapped_patch_step_roundtrip_has_one_metal_part_and_named_load_sheets(tmp_path):
    cq = pytest.importorskip("cadquery")
    del cq
    from src.cad_export import export_unit_cell_step, import_named_parts

    comp = CompositeGeometry((AxisRect(-2.5e-3, 0, 3e-3, 6e-3), AxisRect(2e-3, 0, 4e-3, 6e-3)))
    loads = (LumpedLoad(-0.5e-3, 0.0, 1e-3, 6e-3, "x", R=10.0), LumpedLoad(-0.5e-3, 0.0, 1e-3, 6e-3, "x", C=1e-13))
    cell = UnitCell(LAT, comp, SUB, Metal(), backing="pec", loads=loads)
    res = export_unit_cell_step(cell, tmp_path / "cell.step")
    assert len(res.load_bboxes_mm) == 2
    parts = import_named_parts(tmp_path / "cell.step")
    assert {"Substrate", "Metal", "Ground", "Load_1", "Load_2"} <= set(parts)
    assert parts["Metal"]["n_solids"] == 2  # the two arms, ONE part
    bb = parts["Metal"]["bbox_mm"]
    assert bb[:4] == pytest.approx((-4.0, 4.0, -3.0, 3.0), abs=1e-9)
    assert parts["Metal"]["volume_mm3"] == pytest.approx(comp.feature_area * 1e6 * 0.035, rel=1e-6)
    for k in (1, 2):
        ld = parts[f"Load_{k}"]
        assert ld["n_solids"] == 0 and ld["n_faces"] == 1
        assert ld["bbox_mm"][:4] == pytest.approx((-1.0, 0.0, -3.0, 3.0), abs=1e-9) and ld["area_mm2"] == pytest.approx(
            6.0
        )
        assert ld["bbox_mm"][4] == pytest.approx(0.0175)  # mid-height of the 35 um slab: bridges the arms' side faces
    # cells without composite geometry or loads export exactly as before: no Load parts, one Metal solid
    plain = UnitCell(LAT, Geometry(kind="rectangle", length=6e-3, width=4e-3), SUB, Metal(), backing="pec")
    export_unit_cell_step(plain, tmp_path / "plain.step")
    names = set(import_named_parts(tmp_path / "plain.step"))
    assert names == {"Substrate", "Metal", "Ground"} and os.path.exists(tmp_path / "plain.step")
