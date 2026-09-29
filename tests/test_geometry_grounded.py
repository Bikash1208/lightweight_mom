"""Phase 12A: ``UnitCell.backing`` (the config-facing stack/backing choice)."""

import pytest

from src.geometry import Geometry, UnitCell
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice

LATTICE = RectangularLattice(4.8e-3, 4.8e-3)
GEOMETRY = Geometry("rectangle", 2.0e-3, 2.0e-3)
SUBSTRATE = Substrate("x", 2.2, 0.0009, 1.575e-3)
METAL = Metal()


def test_default_backing_is_air():
    uc = UnitCell(LATTICE, GEOMETRY, SUBSTRATE, METAL)
    assert uc.backing == "air"


def test_explicit_air_and_pec_accepted():
    assert UnitCell(LATTICE, GEOMETRY, SUBSTRATE, METAL, backing="air").backing == "air"
    assert UnitCell(LATTICE, GEOMETRY, SUBSTRATE, METAL, backing="pec").backing == "pec"


def test_invalid_backing_rejected():
    with pytest.raises(ValueError):
        UnitCell(LATTICE, GEOMETRY, SUBSTRATE, METAL, backing="copper")


def test_positional_construction_still_works_unchanged():
    """Every existing call site constructs UnitCell(lattice, geometry, substrate, metal); the new
    field must not disturb that (Phase 12A design requirement 2's spirit, applied to UnitCell)."""
    uc = UnitCell(LATTICE, GEOMETRY, SUBSTRATE, METAL)
    assert uc.lattice is LATTICE and uc.geometry is GEOMETRY
    assert uc.substrate is SUBSTRATE and uc.metal is METAL
    assert uc.backing == "air"
