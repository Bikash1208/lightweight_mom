"""Part 0 task 3 (PHASE12D_ADDENDUM.md): G4/G6 HFSS ingestion and the resonance error budget.
Every assertion is checked against numbers already recorded (by hand) in PHASE12D_ADDENDUM.md, so
this module's own, independently-built pipeline is shown to reproduce them, not merely assert a
new answer."""

import pytest

from src.phase13 import hfss_resonance_ingest as ing

pytestmark = pytest.mark.skipif(
    not (ing.G4_SHEET_TOUCHSTONE.exists() and ing.G6_TOUCHSTONE.exists()),
    reason="real HFSS export (g4.s36p / g6.s36p) not present in this checkout",
)


@pytest.mark.hfss_data
def test_ingest_g4_sheet_reproduces_addendum_dphase_near_resonance():
    """PHASE12D_ADDENDUM.md's 'Update': g4.s36p is the sheet run; at 12.5 GHz the sheet
    comparison shows about -9.8 deg (see G4_sheet_comparison.csv row, independently re-derived
    here via `compare_single_cell`, not copied)."""
    comparison = ing.ingest_g4_sheet()
    row = dict(zip(comparison.frequency_ghz.tolist(), comparison.dphase_mode1_deg.tolist()))
    assert row[12.5] == pytest.approx(-9.849, abs=0.05)
    assert comparison.max_dmag_db < 0.05


@pytest.mark.hfss_data
def test_ingest_g6_reproduces_addendum_band_edge_trend():
    """PHASE12D_ADDENDUM.md: G6 (L=7.6mm, near ITS OWN resonance below the swept band):
    -9.6 deg at 10 GHz -> -0.8 deg at 15 GHz."""
    comparison = ing.ingest_g6()
    row = dict(zip(comparison.frequency_ghz.tolist(), comparison.dphase_mode1_deg.tolist()))
    assert row[10.0] == pytest.approx(-9.6, abs=0.5)
    assert row[15.0] == pytest.approx(-0.8, abs=0.5)
    assert abs(row[10.0]) > abs(row[15.0])  # monotonic improvement moving away from resonance


@pytest.mark.slow
def test_converged_mom_curve_reproduces_addendum_g4_table():
    """PHASE12D_ADDENDUM.md's converged-MoM table at G4: -106.79/-174.65/+126.73 deg at
    11.5/12.5/13.5 GHz."""
    points = {p.frequency_ghz: p for p in ing.converged_mom_curve(ing.L_MM["G4"])}
    assert points[11.5].extrapolated_phase_deg == pytest.approx(-106.79, abs=1.5)
    assert points[12.5].extrapolated_phase_deg == pytest.approx(-174.65, abs=1.5)
    assert points[13.5].extrapolated_phase_deg == pytest.approx(126.73, abs=1.5)
    for p in points.values():
        assert p.status == "CONVERGED"


@pytest.mark.slow
def test_g4_error_budget_reproduces_addendum_numbers_at_12p5ghz():
    """PHASE12D_ADDENDUM.md's budget at resonance (12.5 GHz, total 13.5 deg):
    mesh ~ 5.9 deg, thickness ~ 3.7 deg, residual ~ 3-4 deg."""
    budgets = {b.frequency_ghz: b for b in ing.error_budget_for_g4()}
    b = budgets[12.5]
    assert b.mesh_deg == pytest.approx(5.9, abs=0.5)
    assert b.thickness_deg == pytest.approx(3.7, abs=0.5)
    assert 2.5 <= b.residual_deg <= 5.0
    assert b.total_hfss_vs_0p2mm_deg == pytest.approx(13.5, abs=0.5)
