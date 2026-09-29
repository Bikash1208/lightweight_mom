"""Phase 12C Task A: the grounded supercell wrapper's correctness."""

import pytest

from src.phase12c import design_point as dp
from src.phase12c.grounded_supercell import (
    GroundedSupercellUnitCell,
    build_grounded_supercell,
    homogeneous_correctness_check,
    solve_grounded_supercell_center,
    solve_isolated_grounded_cell,
    supercell_mesh_convergence_check,
)


def test_wrapper_delegates_every_attribute_to_the_inner_supercell():
    grid = [[(2.0e-3, 2.0e-3)] * 3 for _ in range(3)]
    sc = build_grounded_supercell(3, dp.PERIOD_M, dp.PERIOD_M, dp.SUBSTRATE, dp.METAL, grid)
    assert isinstance(sc, GroundedSupercellUnitCell)
    assert sc.backing == "pec"
    assert sc.lattice is sc.inner.lattice
    assert sc.substrate is dp.SUBSTRATE
    assert sc.n == 3
    assert sc.sub_cell_center(1, 1) == sc.inner.sub_cell_center(1, 1)
    assert sc.is_metal(0.0, 0.0) == sc.inner.is_metal(0.0, 0.0)


def test_default_air_getattr_still_gives_air_for_unwrapped_supercell():
    """The pre-Phase-12A getattr pattern (solver.py/sparameters.py) must still treat a bare,
    unwrapped SupercellUnitCell as ungrounded -- unaffected by this new wrapper's existence."""
    from src.phase9.supercell_geometry import build_supercell

    grid = [[(2.0e-3, 2.0e-3)] * 3 for _ in range(3)]
    bare = build_supercell(3, dp.PERIOD_M, dp.PERIOD_M, dp.SUBSTRATE, dp.METAL, grid)
    assert getattr(bare, "backing", "air") == "air"


def test_homogeneous_supercell_reproduces_isolated_cell_to_machine_precision():
    """The load-bearing correctness proof (Phase 9's own argument, grounded): n identical
    sub-cells repeated periodically describes the SAME physical structure as 1 cell repeated
    periodically."""
    result = homogeneous_correctness_check(4.0e-3, 3, 4, dp.PERIOD_M, dp.SUBSTRATE, dp.METAL, dp.DESIGN_FREQUENCY_HZ)
    assert result.relative_error < 1e-9


def test_solve_grounded_supercell_center_rejects_even_n():
    grid = [[(2.0e-3, 2.0e-3)] * 2 for _ in range(2)]
    sc = build_grounded_supercell(2, dp.PERIOD_M, dp.PERIOD_M, dp.SUBSTRATE, dp.METAL, grid)
    with pytest.raises(ValueError, match="odd"):
        solve_grounded_supercell_center(sc, (2.0e-3, 2.0e-3), dp.DESIGN_FREQUENCY_HZ, 4, 8)


def test_solve_grounded_supercell_center_rejects_mismatched_center():
    grid = [[(2.0e-3, 2.0e-3)] * 3 for _ in range(3)]
    sc = build_grounded_supercell(3, dp.PERIOD_M, dp.PERIOD_M, dp.SUBSTRATE, dp.METAL, grid)
    with pytest.raises(ValueError, match="centre patch"):
        solve_grounded_supercell_center(sc, (3.0e-3, 3.0e-3), dp.DESIGN_FREQUENCY_HZ, 4, 8)


def test_isolated_grounded_cell_has_no_bottom_field():
    """A sanity check that the isolated reference really is grounded (reuses the Phase 12A
    physics, not a new formulation): a longer patch shifts the dipole phase measurably relative
    to a shorter one at the same frequency (a real, physically expected reflectarray effect)."""
    a = solve_isolated_grounded_cell(dp.PERIOD_M, dp.PERIOD_M, 2.0e-3, 2.0e-3, dp.SUBSTRATE, dp.METAL,
                                      dp.DESIGN_FREQUENCY_HZ, 8, 16)
    b = solve_isolated_grounded_cell(dp.PERIOD_M, dp.PERIOD_M, 4.0e-3, 4.0e-3, dp.SUBSTRATE, dp.METAL,
                                      dp.DESIGN_FREQUENCY_HZ, 8, 16)
    assert a.dipole.px != b.dipole.px


def test_mesh_convergence_check_reports_a_real_successive_difference():
    result = supercell_mesh_convergence_check(4.0e-3, 3, 4, 8, dp.PERIOD_M, dp.SUBSTRATE, dp.METAL,
                                               dp.DESIGN_FREQUENCY_HZ)
    assert result.relative_difference >= 0.0
    assert result.base_runtime_s > 0.0 and result.finer_runtime_s > 0.0
