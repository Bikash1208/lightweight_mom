"""Phase 8C Sec. 14-15, 55: power balance and propagating-mode analysis."""

import numpy as np
import pytest

from src.constants import C0
from src.phase8c import power_balance as pb


class _FakeRecord:
    def __init__(self, s):
        self.s = np.asarray(s, dtype=complex)


def test_power_balance_matches_known_magnitudes_squared():
    s = np.zeros((4, 4), dtype=complex)
    s[0, 0] = 0.6 + 0.0j  # S11_xx
    s[1, 0] = 0.1 + 0.0j  # S11_yx
    s[2, 0] = 0.8 + 0.0j  # S21_xx
    s[3, 0] = 0.0 + 0.0j  # S21_yx
    result = pb.power_balance_for_record(_FakeRecord(s))
    assert result.R_co == pytest.approx(0.36)
    assert result.R_cross == pytest.approx(0.01)
    assert result.T_co == pytest.approx(0.64)
    assert result.T_cross == pytest.approx(0.0)
    assert result.R_plus_T == pytest.approx(0.36 + 0.01 + 0.64 + 0.0)


def test_power_balance_note_warns_about_mode_dependence():
    s = np.zeros((4, 4), dtype=complex)
    s[0, 0] = 0.5
    s[2, 0] = 0.5
    result = pb.power_balance_for_record(_FakeRecord(s))
    assert "propagating diffraction order" in result.note


# ------------------------------------------------------------------------------- Test: multi-mode


def test_propagating_mode_check_reports_only_fundamental_for_locked_period():
    from src.geometry import Geometry, UnitCell
    from src.materials import Metal, Substrate
    from src.periodic_bc import RectangularLattice

    lattice = RectangularLattice(4.8e-3, 4.8e-3)
    substrate = Substrate(name="Rogers RT/Duroid 5880", eps_r=2.2, tan_delta=0.0009, thickness=1.575e-3)
    geometry = Geometry(kind="rectangle", length=1.0e-3, width=1.0e-3)
    unit_cell = UnitCell(lattice, geometry, substrate, Metal())
    result = pb.propagating_mode_check(unit_cell, [10.0e9, 15.0e9])
    assert result["only_fundamental_mode_propagates"] is True
    assert "Only the fundamental" in result["note"]


def test_propagating_mode_check_warns_when_higher_orders_exist():
    from src.geometry import Geometry, UnitCell
    from src.materials import Metal, Substrate
    from src.periodic_bc import RectangularLattice

    lattice = RectangularLattice(30.0e-3, 30.0e-3)  # large period -> low grating onset
    substrate = Substrate(name="Rogers RT/Duroid 5880", eps_r=2.2, tan_delta=0.0009, thickness=1.575e-3)
    geometry = Geometry(kind="rectangle", length=1.0e-3, width=1.0e-3)
    unit_cell = UnitCell(lattice, geometry, substrate, Metal())
    onset_hz = C0 / 30.0e-3
    result = pb.propagating_mode_check(unit_cell, [onset_hz * 0.5, onset_hz * 1.5])
    assert result["only_fundamental_mode_propagates"] is False
    assert "NOT fully describe" in result["note"]


def test_library_power_summary_on_real_data_shows_negligible_cross_pol(phase8c_full_database):
    summary = pb.library_power_summary(phase8c_full_database, 100)
    assert summary["n_cells"] == 324
    assert summary["mean_R_cross"] < 1e-6
    assert summary["mean_T_cross"] < 1e-6
    assert 0.0 < summary["mean_R_plus_T"] <= 1.0


def test_real_library_reciprocity_and_lossless_bound(phase8c_full_database):
    summary = pb.library_power_summary(phase8c_full_database, 100)
    # tan_delta = 0.0009 (small loss): R+T must be close to, but not exceed, 1.
    assert summary["max_R_plus_T"] <= 1.0 + 1e-6
    assert summary["min_R_plus_T"] > 0.9
