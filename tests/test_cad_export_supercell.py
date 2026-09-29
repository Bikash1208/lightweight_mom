"""Phase 12C: supercell STEP export (cad_export.build_supercell_solids/export_supercell_step),
additive to the existing single-UnitCell export path. Skips without cadquery."""

import pytest

cadquery = pytest.importorskip("cadquery")

from src import cad_export
from src.phase12c.grounded_supercell import build_grounded_supercell
from src.phase9.supercell_geometry import build_supercell
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice

SUBSTRATE = Substrate("RT5880", 2.2, 0.0009, 1.575e-3)
METAL = Metal()
PERIOD_M = 8.0e-3


def _grid(n, l_m=2.0e-3):
    return [[(l_m, l_m)] * n for _ in range(n)]


def test_build_supercell_solids_gives_one_metal_solid_per_subcell():
    sc = build_supercell(3, PERIOD_M, PERIOD_M, SUBSTRATE, METAL, _grid(3))
    substrate, metals = cad_export.build_supercell_solids(sc)
    assert len(metals) == 9
    assert cad_export.SolidProperties.of(substrate).n_solids == 1


def test_grounded_supercell_ground_solid_exists_ungrounded_does_not():
    grounded = build_grounded_supercell(3, PERIOD_M, PERIOD_M, SUBSTRATE, METAL, _grid(3))
    ungrounded = build_supercell(3, PERIOD_M, PERIOD_M, SUBSTRATE, METAL, _grid(3))
    assert cad_export.build_supercell_ground_solid(grounded) is not None
    assert cad_export.build_supercell_ground_solid(ungrounded) is None


def test_export_supercell_step_writes_a_valid_file(tmp_path):
    sc = build_grounded_supercell(3, PERIOD_M, PERIOD_M, SUBSTRATE, METAL, _grid(3))
    path = tmp_path / "supercell.step"
    result = cad_export.export_supercell_step(sc, path)
    assert path.exists()
    assert len(result.metals) == 9
    assert result.ground is not None
    imported = cadquery.importers.importStep(str(path))
    assert cad_export.solid_validity(imported)
    assert len(imported.val().Solids()) == 1 + 9 + 1   # substrate + 9 metals + ground


def test_export_supercell_step_ungrounded_has_no_ground_solid(tmp_path):
    sc = build_supercell(3, PERIOD_M, PERIOD_M, SUBSTRATE, METAL, _grid(3))
    path = tmp_path / "supercell_ungrounded.step"
    result = cad_export.export_supercell_step(sc, path)
    assert result.ground is None
    imported = cadquery.importers.importStep(str(path))
    assert len(imported.val().Solids()) == 1 + 9
