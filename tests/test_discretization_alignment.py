"""Phase 12B: mesh-alignment and duplicate-geometry-detection helpers in discretization.py."""

import pytest

from src.basis_functions import MeshSettings
from src.discretization import (
    duplicate_geometry_groups,
    feature_edges_aligned,
    geometry_edges_aligned,
    is_position_mesh_aligned,
    metal_mask_key,
    smallest_aligned_cells,
)
from src.geometry import Geometry, UnitCell
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice

P = 4.8e-3
LATTICE = RectangularLattice(P, P)
SUBSTRATE = Substrate("x", 2.2, 0.0009, 1.575e-3)
METAL = Metal()


def test_mesh_line_positions_are_exactly_aligned():
    # N=48 over 4.8mm -> dx=0.1mm; -2.4 + k*0.1 for k=0..48 are exactly the mesh lines
    for k in range(49):
        pos = -P / 2 + k * (P / 48)
        assert is_position_mesh_aligned(pos, P, 48)


def test_off_grid_position_is_not_aligned():
    assert not is_position_mesh_aligned(-P / 2 + 0.5 * (P / 48), P, 48)


@pytest.mark.parametrize("n,expected", [(8, False), (16, False), (24, True), (32, False), (48, True), (96, True)])
def test_feature_edges_aligned_matches_hand_derivation(n, expected):
    assert feature_edges_aligned(P, n, 2.8e-3, 0.0) == expected


def test_smallest_aligned_cells_for_the_full_phase7_family_is_48():
    lengths_mm = [round(1.0 + 0.2 * i, 1) for i in range(18)]
    features = [(l * 1e-3, 0.0) for l in lengths_mm]
    assert smallest_aligned_cells(P, features) == 48


def test_smallest_aligned_cells_returns_none_when_unreachable_within_bound():
    # An irrational-like (in practice, very awkward) length relative to the period, over a tiny
    # search bound, should not silently find a "close enough" mesh.
    result = smallest_aligned_cells(P, [(1.0e-3 + 1e-7, 0.0)], max_cells=16)
    assert result is None


def test_smallest_aligned_cells_rejects_a_too_small_bound():
    with pytest.raises(ValueError):
        smallest_aligned_cells(P, [(1.0e-3, 0.0)], max_cells=1)


def test_geometry_edges_aligned_false_for_any_rotation():
    geo = Geometry("rectangle", 2.8e-3, 2.8e-3, rotation=0.1)
    mesh = MeshSettings(48, 48, 48)
    assert geometry_edges_aligned(geo, LATTICE, mesh) is False


def test_geometry_edges_aligned_true_on_the_aligned_mesh():
    geo = Geometry("rectangle", 2.8e-3, 2.8e-3)
    assert geometry_edges_aligned(geo, LATTICE, MeshSettings(48, 48, 48)) is True
    assert geometry_edges_aligned(geo, LATTICE, MeshSettings(32, 32, 32)) is False


def test_metal_mask_key_identifies_the_known_duplicate_pairs_at_32x32():
    mesh = MeshSettings(32, 32, 32)
    items = []
    for l_mm in (1.4, 1.6, 2.0, 2.2, 2.8):
        uc = UnitCell(LATTICE, Geometry("rectangle", l_mm * 1e-3, l_mm * 1e-3), SUBSTRATE, METAL)
        items.append((l_mm, metal_mask_key(uc, mesh)))
    groups = duplicate_geometry_groups(items)
    assert groups == ((1.4, 1.6), (2.0, 2.2))


def test_metal_mask_key_no_duplicates_on_the_aligned_mesh():
    mesh = MeshSettings(48, 48, 48)
    items = []
    for l_mm in (1.4, 1.6, 2.0, 2.2, 2.8):
        uc = UnitCell(LATTICE, Geometry("rectangle", l_mm * 1e-3, l_mm * 1e-3), SUBSTRATE, METAL)
        items.append((l_mm, metal_mask_key(uc, mesh)))
    assert duplicate_geometry_groups(items) == ()


def test_duplicate_geometry_groups_omits_singletons():
    groups = duplicate_geometry_groups([("a", b"x"), ("b", b"y"), ("c", b"x")])
    assert groups == (("a", "c"),)
