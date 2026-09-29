"""Phase 13D: the real MoM surface-wave eigenmode search (`src/dispersion/mom_scaffold.py`)."""

import numpy as np
import pytest

from src.dispersion.analytic_tm0 import tm0_beta
from src.dispersion.mom_scaffold import (
    DEFAULT_MESH,
    SCAFFOLD_LABELS,
    _locate_peak,
    build_dispersion_table_from_mom,
    find_surface_wave_mode,
    square_patch_unit_cell,
)
from src.errors import LwmomError
from src.geometry import Geometry, UnitCell
from src.labels import DISPERSION_FROM_MOM, EXPERIMENTAL
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice

# RT/duroid 5880, 1.575 mm, at 10.54 GHz -- the walk-through board (PHASE13D brief): analytic TM0
# beta = 224.98 rad/m.
FREQUENCY_HZ = 10.54e9
SUBSTRATE = Substrate("x", 2.2, 0.0009, 1.575e-3)
METAL = Metal()
PERIOD_M = 4.0e-3


def _grounded_cell(patch_side_m):
    return square_patch_unit_cell(PERIOD_M, SUBSTRATE, METAL, patch_side_m)


@pytest.mark.slow
def test_find_surface_wave_mode_matches_analytic_tm0_for_a_lightly_loaded_cell():
    """A modest patch (0.5x the period) weakly loads the bare-slab TM0 pole -- the found beta
    should sit close to the closed-form bare-slab value (the scaffold's own requirement 3)."""
    beta_analytic = tm0_beta(FREQUENCY_HZ, SUBSTRATE.eps_r, SUBSTRATE.thickness)
    cell = _grounded_cell(PERIOD_M * 0.5)

    mode = find_surface_wave_mode(cell, FREQUENCY_HZ, mesh=DEFAULT_MESH, n_scan=500)

    assert mode.beta_rad_per_m == pytest.approx(beta_analytic, abs=1.0)
    assert mode.prominence > 3.0
    assert mode.x_s_ohm is not None
    assert set(mode.labels) == set(SCAFFOLD_LABELS)
    assert DISPERSION_FROM_MOM in mode.labels and EXPERIMENTAL in mode.labels


@pytest.mark.slow
def test_find_surface_wave_mode_requires_grounded_backing():
    lattice = RectangularLattice(PERIOD_M, PERIOD_M)
    geometry = Geometry("rectangle", PERIOD_M * 0.5, PERIOD_M * 0.5)
    ungrounded = UnitCell(lattice, geometry, SUBSTRATE, METAL, backing="air")
    with pytest.raises(LwmomError) as excinfo:
        find_surface_wave_mode(ungrounded, FREQUENCY_HZ)
    assert excinfo.value.code == "CONFIG_INVALID"


def test_locate_peak_raises_on_flat_scan():
    """A synthetic, deterministic exercise of the no-clear-mode path -- not dependent on finding a
    real physical configuration that happens to be ambiguous."""
    kts = np.linspace(100.0, 200.0, 50)
    vals = np.full_like(kts, 1.0)
    with pytest.raises(LwmomError) as excinfo:
        _locate_peak(kts, vals, min_prominence=3.0)
    assert excinfo.value.code == "DISPERSION_MOM_NO_CLEAR_MODE"


def test_locate_peak_raises_when_maximum_is_at_the_scan_boundary():
    kts = np.linspace(100.0, 200.0, 50)
    vals = np.linspace(1.0, 10.0, 50)  # monotonic -- max sits at the last point, not an interior peak
    with pytest.raises(LwmomError) as excinfo:
        _locate_peak(kts, vals, min_prominence=3.0)
    assert excinfo.value.code == "DISPERSION_MOM_NO_CLEAR_MODE"


def test_locate_peak_accepts_a_clear_interior_peak():
    kts = np.linspace(100.0, 200.0, 51)
    vals = np.full_like(kts, 1.0)
    vals[25] = 50.0
    lo, hi, prominence = _locate_peak(kts, vals, min_prominence=3.0)
    assert lo == kts[24] and hi == kts[26]
    assert prominence > 3.0


@pytest.mark.slow
def test_build_dispersion_table_from_mom_has_beta_and_xs_columns():
    geometry_params = [PERIOD_M * 0.4, PERIOD_M * 0.5, PERIOD_M * 0.6]
    table = build_dispersion_table_from_mom(
        lambda gp: _grounded_cell(gp),
        geometry_params,
        FREQUENCY_HZ,
        mesh=DEFAULT_MESH,
    )
    assert table.label == DISPERSION_FROM_MOM
    assert len(table.samples) >= 1
    for sample in table.samples:
        assert sample.beta_rad_per_m is not None
        assert sample.x_s_ohm is not None
    # beta_at/geometry_for_x_s (used by Formulation A) both work off this table
    lo, hi = table.geometry_range
    mid_beta = table.beta_at((lo + hi) / 2.0) if lo != hi else table.beta_at(lo)
    assert mid_beta > 0


def test_build_dispersion_table_from_mom_raises_when_every_sample_is_skipped(monkeypatch):
    import src.dispersion.mom_scaffold as mod

    def _always_no_mode(*args, **kwargs):
        raise mod.dispersion_mom_no_clear_mode("synthetic: forced skip for this test")

    monkeypatch.setattr(mod, "find_surface_wave_mode", _always_no_mode)
    with pytest.raises(LwmomError) as excinfo:
        build_dispersion_table_from_mom(lambda gp: None, [1.0, 2.0], FREQUENCY_HZ)
    assert excinfo.value.code == "DISPERSION_MOM_NO_CLEAR_MODE"
