"""Additive HFSS delta table of the CP-2 cell: approval gate, exact reproduction of the HFSS points, interpolation, extrapolation label, reim reader and plane fit."""

import math
from pathlib import Path

import numpy as np
import pytest

from src.ris.hfss_calibration import (
    CalibrationNotApproved,
    DeltaTable,
    fit_reference_plane_mm,
    read_hfss_reim,
)
from src.ris.polarization import circular_reflection, to_structure_frame
from src.ris.state_library import StateLibrary

ROOT = Path(__file__).resolve().parents[1]
LIB = ROOT / "results" / "phase20_r4" / "cp2_library_v1.npz"
HF = ROOT / "hfss_results" / "phase20_r3_cp_cell"
PKG = ROOT / "hfss_exports" / "phase20_r3_cp_cell"


@pytest.fixture(scope="module")
def setup():
    lib = StateLibrary.load(LIB)
    d, _ = fit_reference_plane_mm(HF / "hfss_bare_slab_theta0_30.csv", PKG / "bare_grounded_slab_reference.csv")
    return lib, d, read_hfss_reim(HF / "hfss_cp2_states1-4_theta0_30.csv", d)


def test_plane_offset_matches_the_advisor_fit(setup):
    _, d, h = setup
    assert d == pytest.approx(1.236, abs=0.002)
    assert len(h) == 72


def test_building_needs_approval_and_provenance(setup):
    lib, _, h = setup
    with pytest.raises(CalibrationNotApproved):
        DeltaTable.build(lib, h, "src")
    with pytest.raises(ValueError):
        DeltaTable.build(lib, h, "", approved_by_user=True)


def test_delta_reproduces_hfss_at_the_hfss_points_and_labels_extrapolation(setup, tmp_path):
    lib, _, h = setup
    t = DeltaTable.build(lib, h, "hfss_results/phase20_r3_cp_cell 2026-10-23", approved_by_user=True)
    for st, f, th, ph in [(1, 10.5, 0.0, 0.0), (3, 10.0, 30.0, 0.0), (4, 9.0, 0.0, 0.0)]:
        c, label = t.apply(lib.structure_circular(f, th, ph, st)[0], f, th, st)
        want = to_structure_frame(circular_reflection(h[(st, f, th, ph)]), 0.0)
        np.testing.assert_allclose(c, want, atol=1e-9)
        assert label == "HFSS_CALIBRATED"
    assert t.at(10.5, 45.0, 1)[1] == "HFSS_CALIBRATED_EXTRAPOLATED"
    # round trip and interpolation between HFSS frequencies is the mean of the neighbours at the midpoint
    p = tmp_path / "t.npz"
    t.save(p)
    t2 = DeltaTable.load(p)
    np.testing.assert_allclose(t2.at(10.25, 15.0, 2)[0], t.at(10.25, 15.0, 2)[0])
    assert math.isclose(t.f[4], 10.0) and math.isclose(t.f[5], 10.5)
    mid = 0.25 * (t.delta[4, 0, 1] + t.delta[5, 0, 1] + t.delta[4, 1, 1] + t.delta[5, 1, 1])
    np.testing.assert_allclose(t.at(10.25, 15.0, 2)[0], mid, atol=1e-12)
