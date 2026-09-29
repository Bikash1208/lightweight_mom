"""Resonance detection (req. 14), grating-risk checking (req. 15) and phase usability -- each
threshold tested at its documented boundary (`library_schema.py`), not merely "does it run"."""

import warnings

import numpy as np
import pytest

from src.constants import C0
from src.geometry import Geometry, UnitCell
from src.library_builder import build_library_database
from src.library_classifier import (
    ResonanceResult,
    classify_cell,
    detect_resonance,
    grating_risk,
    phase_usability,
)
from src.library_schema import PHASE_USABLE_MIN_MAG_DB, RESONANCE_MIN_PROMINENCE, LibraryConfig
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice
from src.phase_extraction import extract_design_frequency

warnings.filterwarnings("ignore", category=UserWarning)

_SUBSTRATE = {"name": "Rogers RT/Duroid 5880", "eps_r": 2.2, "tan_delta": 0.0009, "thickness_m": 1.575e-3}


# --------------------------------------------------------------------------- resonance detection


def test_detect_resonance_flat_curve_is_not_resonant():
    freqs = np.linspace(10e9, 15e9, 21)
    mag = np.full(21, 0.2)
    result = detect_resonance(freqs, mag)
    assert result.is_resonant is False
    assert result.resonance_frequency_hz is None


def test_detect_resonance_monotonic_curve_is_not_resonant():
    """No interior local maximum at all -- a curve that simply rises toward one edge (like
    Phase 6's case01) must not be flagged, since its true peak (if any) is outside the band."""
    freqs = np.linspace(10e9, 15e9, 21)
    mag = np.linspace(0.1, 0.3, 21)
    result = detect_resonance(freqs, mag)
    assert result.is_resonant is False


def test_detect_resonance_clear_interior_peak_is_resonant():
    freqs = np.linspace(10e9, 15e9, 21)
    mag = 0.1 + 0.85 * np.exp(-0.5 * ((freqs - 12.5e9) / 0.3e9) ** 2)  # peak ~0.95 at 12.5 GHz
    result = detect_resonance(freqs, mag)
    assert result.is_resonant is True
    assert result.resonance_frequency_hz == pytest.approx(12.5e9, abs=0.3e9)
    assert result.resonance_metric > RESONANCE_MIN_PROMINENCE


def test_detect_resonance_edge_peak_is_excluded():
    """A maximum at the very first or last sample (band edge) is not counted -- only strictly
    interior samples are candidates (req. 14: "identify whether it occurs inside 10-15 GHz")."""
    freqs = np.linspace(10e9, 15e9, 11)
    mag = np.linspace(0.9, 0.1, 11)  # monotonically decreasing: "peak" only at index 0 (edge)
    result = detect_resonance(freqs, mag)
    assert result.is_resonant is False


def test_detect_resonance_boundary_prominence_classification():
    """A peak whose prominence sits exactly at (just above/below) RESONANCE_MIN_PROMINENCE
    flips classification there, and nowhere else -- the threshold is meaningful, not vacuous."""
    freqs = np.linspace(10e9, 15e9, 11)
    base = 0.1
    for prominence, expect_resonant in (
        (RESONANCE_MIN_PROMINENCE + 0.01, True),
        (RESONANCE_MIN_PROMINENCE - 0.01, False),
    ):
        mag = np.full(11, base)
        mag[5] = base + prominence
        result = detect_resonance(freqs, mag)
        assert result.is_resonant is expect_resonant, f"prominence={prominence}"
        assert result.resonance_metric == pytest.approx(prominence, abs=1e-9)


def test_detect_resonance_too_few_points_is_not_resonant():
    result = detect_resonance([10e9, 11e9], [0.5, 0.9])
    assert result.is_resonant is False
    assert isinstance(result, ResonanceResult)
    assert "prominence" in result.criterion


# --------------------------------------------------------------------------------- grating risk


def _cell(length_mm, width_mm, period_mm=4.8):
    sub = Substrate(name="test", eps_r=2.2, tan_delta=0.0009, thickness=1.575e-3)
    lat = RectangularLattice(period_x=period_mm * 1e-3, period_y=period_mm * 1e-3)
    geo = Geometry(kind="rectangle", length=length_mm * 1e-3, width=width_mm * 1e-3)
    return UnitCell(lattice=lat, geometry=geo, substrate=sub, metal=Metal())


def test_grating_risk_false_for_locked_period_across_the_full_band():
    """period=4.8 mm -> onset = c0/4.8mm ~= 62.5 GHz, far above the 10-15 GHz sweep: computed,
    not assumed."""
    uc = _cell(2.0, 2.0, period_mm=4.8)
    result = grating_risk(uc, np.linspace(10e9, 15e9, 11))
    assert result["is_grating_risk"] is False
    assert result["grating_onset_hz"] == pytest.approx(C0 / 4.8e-3, rel=1e-6)


def test_grating_risk_true_when_a_frequency_exceeds_onset():
    """A large period pushed the onset below the swept band -> flagged, using the existing
    discretization.is_specular_only check (not a new period/frequency computation)."""
    uc = _cell(2.0, 2.0, period_mm=30.0)  # onset = c0/30mm = 10 GHz
    result = grating_risk(uc, np.linspace(10e9, 15e9, 11))
    assert result["is_grating_risk"] is True
    assert result["grating_onset_hz"] == pytest.approx(C0 / 30.0e-3, rel=1e-6)


# ---------------------------------------------------------------------------- phase usability


def test_phase_usability_boundary():
    assert phase_usability(PHASE_USABLE_MIN_MAG_DB + 0.01) == "PHASE_USABLE"
    assert phase_usability(PHASE_USABLE_MIN_MAG_DB - 0.01) == "PHASE_UNUSABLE"
    assert phase_usability(PHASE_USABLE_MIN_MAG_DB) == "PHASE_USABLE"  # >= is usable


def test_phase_usability_none_is_unclassified():
    assert phase_usability(None) == "UNCLASSIFIED"


# -------------------------------------------------------------------------------- classify_cell


def test_classify_cell_invalid_geometry_is_unclassified():
    from src.database import CellEntry

    cell_entry = CellEntry(0, "bad", {"length_m": 1e-3}, "invalid", "extends outside the cell")
    from src.phase_extraction import DesignFrequencyRow

    row = DesignFrequencyRow("bad", 1.0, 1.0, "INVALID", "FAILED", None, None, None, None)
    result = classify_cell(cell_entry, None, row)
    assert result.geometry_status == "INVALID"
    assert result.solver_status == "FAILED"
    assert result.phase_status == "UNCLASSIFIED"
    assert result.resonance is None


def test_classify_cell_on_a_real_small_database():
    cfg = LibraryConfig(
        period_m=4.8e-3,
        substrate=_SUBSTRATE,
        f_min_hz=10.0e9,
        f_max_hz=15.0e9,
        f_step_hz=1.0e9,
        f_design_hz=12.0e9,
        l_values_mm=(1.0, 2.4, 4.2),
        w_values_mm=(1.0, 2.4, 4.2),
    )
    database = build_library_database(cfg, out_dir=None).database
    rows = extract_design_frequency(database, cfg)
    by_lw = {(r.l_mm, r.w_mm): r for r in rows}
    for cell in database.cells:
        row = by_lw[(round(cell.parameters["length_m"] * 1e3, 9), round(cell.parameters["width_m"] * 1e3, 9))]
        result = classify_cell(cell, database, row)
        assert result.geometry_status == "VALID"
        assert result.solver_status == "OK"
        assert result.phase_status in ("PHASE_USABLE", "PHASE_UNUSABLE")
        assert result.resonance is not None
        assert result.grating is not None
        assert result.grating["is_grating_risk"] is False  # locked period, full band
