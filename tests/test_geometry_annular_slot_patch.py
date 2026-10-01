"""Phase 15 Part E: the "patch with annular slot" unit-cell geometry (Mishra, Barad, Sarkar, IEEE
TVT 2026's own two-state 1-bit cell) -- construction/validation, metal-mask correctness, and a real
MoM solve (`@pytest.mark.slow`), the same way `test_geometry.py` validates the square patch."""

import math

import numpy as np
import pytest

from src.geometry import AnnularSlotPatch, Geometry, UnitCell
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice

SUB = Substrate("t", 2.2, 0.0009, 1.575e-3)


def _grid_fill_fraction(unit_cell, n=800):
    px, py = unit_cell.lattice.period_x, unit_cell.lattice.period_y
    x = (np.arange(n) + 0.5) / n * px - px / 2
    y = (np.arange(n) + 0.5) / n * py - py / 2
    xx, yy = np.meshgrid(x, y, indexing="ij")
    return unit_cell.is_metal(xx, yy).mean()


# ------------------------------------------------------------- construction/validation


def test_rejects_non_positive_a_or_b():
    with pytest.raises(ValueError):
        AnnularSlotPatch(a=0.0, b=2e-3, R=0.5e-3, r=0.2e-3)


def test_rejects_r_not_less_than_R():
    with pytest.raises(ValueError):
        AnnularSlotPatch(a=2e-3, b=2e-3, R=0.5e-3, r=0.5e-3)
    with pytest.raises(ValueError):
        AnnularSlotPatch(a=2e-3, b=2e-3, R=0.3e-3, r=0.5e-3)


def test_rejects_slot_bigger_than_the_patch():
    with pytest.raises(ValueError):
        AnnularSlotPatch(a=2e-3, b=2e-3, R=1.5e-3, r=0.1e-3)  # R > min(a,b)/2 = 1e-3


def test_polarity_is_metal():
    assert AnnularSlotPatch(a=2e-3, b=2e-3, R=0.7e-3, r=0.4e-3).polarity == "metal"


def test_feature_area_is_patch_minus_ring():
    a, b, outer_r, inner_r = 2e-3, 2e-3, 0.7e-3, 0.4e-3
    geo = AnnularSlotPatch(a=a, b=b, R=outer_r, r=inner_r)
    expected = a * b - math.pi * (outer_r**2 - inner_r**2)
    assert geo.feature_area == pytest.approx(expected)


def test_in_feature_matches_patch_and_not_ring():
    geo = AnnularSlotPatch(a=2e-3, b=2e-3, R=0.7e-3, r=0.4e-3)
    assert geo.in_feature(0.0, 0.0)  # centre: inside the inner disk (r=0.4mm), so metal
    assert not geo.in_feature(0.5e-3, 0.0)  # inside the ring [0.4, 0.7]mm: slot, not metal
    assert geo.in_feature(0.9e-3, 0.0)  # outside the ring but inside the patch: metal
    assert not geo.in_feature(1.5e-3, 0.0)  # outside the whole patch: not metal


def test_corners_are_the_outer_bounding_box():
    geo = AnnularSlotPatch(a=2e-3, b=3e-3, R=0.5e-3, r=0.2e-3)
    corners = geo.corners()
    assert np.allclose(np.abs(corners[:, 0]), 1e-3)
    assert np.allclose(np.abs(corners[:, 1]), 1.5e-3)


# ------------------------------------------------------------------ UnitCell integration


def _annular_cell(period=6e-3, a=4e-3, b=4e-3, outer_r=1.5e-3, inner_r=1.0e-3, backing="pec"):
    geo = AnnularSlotPatch(a=a, b=b, R=outer_r, r=inner_r)
    return UnitCell(RectangularLattice(period, period), geo, SUB, Metal(), backing=backing)


def test_unit_cell_rejects_a_slot_feature_extending_outside_the_lattice():
    with pytest.raises(ValueError):
        _annular_cell(period=3e-3, a=4e-3, b=4e-3, outer_r=1.5e-3, inner_r=1.0e-3)  # a > period


def test_unit_cell_metal_fill_fraction_is_between_plain_patch_and_zero():
    plain = UnitCell(RectangularLattice(6e-3, 6e-3), Geometry("rectangle", 4e-3, 4e-3), SUB, Metal(), backing="pec")
    annular = _annular_cell()
    assert 0.0 < annular.metal_fill_fraction < plain.metal_fill_fraction


def test_discretized_fill_fraction_matches_the_continuous_feature_area():
    uc = _annular_cell()
    measured = _grid_fill_fraction(uc)
    assert measured == pytest.approx(uc.metal_fill_fraction, abs=0.01)


def test_zero_inner_radius_is_a_plain_disk_slot_with_no_separate_island():
    geo = AnnularSlotPatch(a=2e-3, b=2e-3, R=0.6e-3, r=0.0)
    assert not geo.in_feature(0.0, 0.0)  # centre now INSIDE the (disk, r=0) "ring" -> slot


# =================================================================================== real MoM solve


@pytest.mark.slow
def test_annular_slot_cell_solves_and_differs_from_a_plain_patch_of_the_same_outer_size():
    """"validated like the square patch" (the brief): a real MoM solve on this new geometry is
    passive/reciprocal-scale (|S11| close to 1 for a lossless PEC grounded reflector) and gives a
    MEASURABLY different reflection phase than a plain patch of the same outer footprint -- the
    slot is doing something physically real, not a geometry that silently reduces to the plain
    patch."""
    from src.basis_functions import MeshSettings
    from src.excitation import ReferencePlanes
    from src.sparameters import compute_s_parameters

    period, outer_side = 6e-3, 4e-3
    mesh = MeshSettings(cells_x=16, cells_y=16, floquet_order=9)
    planes = ReferencePlanes(top_air_thickness=5e-3, bottom_air_thickness=None, substrate_thickness=SUB.thickness)

    annular = _annular_cell(period=period, a=outer_side, b=outer_side, outer_r=1.5e-3, inner_r=1.0e-3)
    plain = UnitCell(
        RectangularLattice(period, period), Geometry("rectangle", outer_side, outer_side), SUB, Metal(), backing="pec"
    )

    s_annular = compute_s_parameters(annular, 12.5e9, mesh, planes, warn_grating=False, warn_surface_wave=False)
    s_plain = compute_s_parameters(plain, 12.5e9, mesh, planes, warn_grating=False, warn_surface_wave=False)

    assert abs(s_annular.s11[0, 0]) == pytest.approx(1.0, abs=0.05)  # lossless PEC: |S11| ~ 1
    phase_annular = math.degrees(np.angle(s_annular.s11[0, 0]))
    phase_plain = math.degrees(np.angle(s_plain.s11[0, 0]))
    phase_diff = abs(((phase_annular - phase_plain + 180.0) % 360.0) - 180.0)
    assert phase_diff > 1.0, f"annular-slot phase ({phase_annular:.1f}deg) barely differs from plain patch ({phase_plain:.1f}deg)"
