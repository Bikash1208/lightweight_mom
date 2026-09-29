"""Sec. 16: the heterogeneous supercell geometry, duck-typed against the existing solver's own
interface."""

import numpy as np
import pytest

from src.materials import Metal, Substrate
from src.phase9.supercell_geometry import SupercellUnitCell, build_supercell

SUB = Substrate("Rogers RT/Duroid 5880", 2.2, 0.0009, 1.575e-3)
MET = Metal()
P = 4.8e-3


def test_build_supercell_period_equals_n_times_subperiod():
    grid = [[(2.8e-3, 2.8e-3)] * 3 for _ in range(3)]
    sc = build_supercell(3, P, P, SUB, MET, grid)
    assert sc.lattice.period_x == pytest.approx(3 * P)
    assert sc.lattice.period_y == pytest.approx(3 * P)


def test_build_supercell_rejects_mismatched_period():
    from src.periodic_bc import RectangularLattice
    from src.geometry import Geometry

    bad_lattice = RectangularLattice(3 * P, 3 * P)
    patches = tuple(tuple(Geometry(kind="rectangle", length=2.8e-3, width=2.8e-3) for _ in range(3)) for _ in range(3))
    with pytest.raises(ValueError):
        SupercellUnitCell(bad_lattice, SUB, MET, 3, P, P * 1.1, patches)   # sub_period_y mismatched


def test_build_supercell_rejects_patch_larger_than_its_subcell():
    grid = [[(5.0e-3, 5.0e-3)] * 3 for _ in range(3)]   # subcell is 4.8mm; 5.0mm patch cannot fit
    with pytest.raises(ValueError):
        build_supercell(3, P, P, SUB, MET, grid)


def test_is_metal_identifies_the_correct_subcell_center():
    # a 3x3 supercell where only the centre subcell has a large patch, others have a tiny one.
    small, big = (0.5e-3, 0.5e-3), (4.0e-3, 4.0e-3)
    grid = [[small, small, small], [small, big, small], [small, small, small]]
    sc = build_supercell(3, P, P, SUB, MET, grid)
    # centre of the centre subcell must be metal (inside the 4mm patch)
    assert bool(sc.is_metal(np.array([0.0]), np.array([0.0]))[0])
    # centre of the corner (0,0) subcell must NOT be metal (only a tiny 0.5mm patch there)
    cx, cy = sc.sub_cell_center(0, 0)
    corner_probe_x, corner_probe_y = cx + 1.5e-3, cy + 1.5e-3   # outside the tiny patch, inside the subcell
    assert not bool(sc.is_metal(np.array([corner_probe_x]), np.array([corner_probe_y]))[0])


def test_metal_area_sums_every_subcell():
    grid = [[(1.0e-3, 1.0e-3)] * 3 for _ in range(3)]
    sc = build_supercell(3, P, P, SUB, MET, grid)
    assert sc.metal_area == pytest.approx(9 * 1.0e-3 * 1.0e-3)


def test_patches_must_be_centred_at_origin():
    from src.geometry import Geometry

    from src.periodic_bc import RectangularLattice

    lattice = RectangularLattice(3 * P, 3 * P)
    off_center = Geometry(kind="rectangle", length=1.0e-3, width=1.0e-3, center_x=1.0e-3)
    patches = tuple(
        tuple(off_center if (i, j) == (0, 0) else Geometry(kind="rectangle", length=1.0e-3, width=1.0e-3) for j in range(3))
        for i in range(3)
    )
    with pytest.raises(ValueError):
        SupercellUnitCell(lattice, SUB, MET, 3, P, P, patches)
