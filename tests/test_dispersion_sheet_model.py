"""Phase 14 Part A: the homogenized sheet-admittance dispersion model (`src/dispersion/sheet_model.py`),
the real replacement for the direct MoM eigenmode search that turned out not to resolve a
loading-dependent mode (see `src/dispersion/mom_scaffold.py`'s own updated docstring)."""

import pytest

from src.basis_functions import MeshSettings
from src.dispersion.analytic_tm0 import tm0_beta
from src.dispersion.mom_scaffold import square_patch_unit_cell
from src.dispersion.sheet_model import (
    SCAFFOLD_LABELS,
    build_dispersion_table_from_sheet_model,
    find_beta_for_cell,
    find_beta_from_sheet_admittance,
    sheet_admittance_from_normal_incidence,
)
from src.errors import LwmomError
from src.geometry import Geometry, UnitCell
from src.labels import DISPERSION_FROM_SHEET_MODEL, EXPERIMENTAL
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice

# RT/duroid 5880, 1.575 mm, at 10.54 GHz -- the walk-through board: analytic TM0 beta = 224.98 rad/m.
FREQUENCY_HZ = 10.54e9
SUBSTRATE = Substrate("x", 2.2, 0.0009, 1.575e-3)
METAL = Metal()
PERIOD_M = 4.0e-3
MESH = MeshSettings(cells_x=16, cells_y=16, floquet_order=10)


def _cell(patch_side_m):
    return square_patch_unit_cell(PERIOD_M, SUBSTRATE, METAL, patch_side_m)


@pytest.mark.slow
def test_small_patch_reproduces_the_analytic_bare_slab_tm0():
    beta_analytic = tm0_beta(FREQUENCY_HZ, SUBSTRATE.eps_r, SUBSTRATE.thickness)
    mode = find_beta_for_cell(_cell(0.5e-3), FREQUENCY_HZ, MESH, period_m=PERIOD_M)
    assert mode.beta_rad_per_m == pytest.approx(beta_analytic, abs=0.5)
    assert mode.x_s_ohm > 0
    assert abs(mode.residual) < 1e-3
    assert set(mode.labels) == set(SCAFFOLD_LABELS)
    assert DISPERSION_FROM_SHEET_MODEL in mode.labels and EXPERIMENTAL in mode.labels


@pytest.mark.slow
def test_beta_increases_monotonically_with_patch_size():
    sizes_m = [0.5e-3, 1.0e-3, 1.6e-3, 2.0e-3, 2.4e-3, 2.8e-3, 3.6e-3]
    betas = [find_beta_for_cell(_cell(l), FREQUENCY_HZ, MESH, period_m=PERIOD_M).beta_rad_per_m for l in sizes_m]
    assert all(b2 >= b1 for b1, b2 in zip(betas, betas[1:]))
    assert betas[-1] > betas[0] + 1.0  # a real, non-trivial spread, not float noise


@pytest.mark.slow
def test_find_beta_for_cell_requires_grounded_backing():
    lattice = RectangularLattice(PERIOD_M, PERIOD_M)
    geometry = Geometry("rectangle", PERIOD_M * 0.5, PERIOD_M * 0.5)
    ungrounded = UnitCell(lattice, geometry, SUBSTRATE, METAL, backing="air")
    with pytest.raises(LwmomError) as excinfo:
        find_beta_for_cell(ungrounded, FREQUENCY_HZ, MESH)
    assert excinfo.value.code == "CONFIG_INVALID"


def test_find_beta_from_sheet_admittance_raises_when_no_root_exists():
    """A synthetic, deterministic exercise of the no-root path: a huge, clearly-non-physical
    Y_sheet pushes the resonance condition's imaginary part to the same sign everywhere in the
    scanned interval."""
    with pytest.raises(LwmomError) as excinfo:
        find_beta_from_sheet_admittance(1e6 + 1e6j, SUBSTRATE, FREQUENCY_HZ, n_scan=200)
    assert excinfo.value.code == "DISPERSION_SHEET_MODEL_NO_ROOT"


@pytest.mark.slow
def test_build_dispersion_table_from_sheet_model_has_beta_and_xs_columns():
    sizes_m = [1.0e-3, 2.0e-3, 3.0e-3]
    table = build_dispersion_table_from_sheet_model(lambda l: _cell(l), sizes_m, FREQUENCY_HZ, MESH, period_m=PERIOD_M)
    assert table.label == DISPERSION_FROM_SHEET_MODEL
    assert len(table.samples) == len(sizes_m)
    for sample in table.samples:
        assert sample.beta_rad_per_m is not None
        assert sample.x_s_ohm is not None
    lo, hi = table.geometry_range
    assert table.beta_at((lo + hi) / 2.0) > 0


def test_build_dispersion_table_from_sheet_model_raises_when_every_sample_is_skipped(monkeypatch):
    import src.dispersion.sheet_model as mod

    def _always_no_root(*args, **kwargs):
        raise mod.dispersion_sheet_model_no_root("synthetic: forced skip for this test")

    monkeypatch.setattr(mod, "find_beta_for_cell", _always_no_root)
    with pytest.raises(LwmomError) as excinfo:
        build_dispersion_table_from_sheet_model(lambda l: None, [1.0, 2.0], FREQUENCY_HZ, MESH)
    assert excinfo.value.code == "DISPERSION_SHEET_MODEL_NO_ROOT"


@pytest.mark.slow
def test_sheet_admittance_is_small_for_a_lightly_loaded_cell():
    """The patch layer's own contribution should vanish as the patch shrinks -- Y_sheet's real
    part in particular, since it should be near-zero (bound, low-loss regime) well before the
    patch is large enough to meaningfully perturb the normal-incidence reflection."""
    y_sheet = sheet_admittance_from_normal_incidence(_cell(0.5e-3), FREQUENCY_HZ, MESH)
    assert abs(y_sheet.real) < 1e-4
    assert abs(y_sheet) < 1e-3
