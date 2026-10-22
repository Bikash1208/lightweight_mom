"""HFSS comparison / approval-gated calibration hook (src/ris/hfss_calibration.py) on synthetic data with known answers."""

import csv
import math

import numpy as np
import pytest

from src.ris import hfss_calibration as hc
from src.ris.polarization import circular_reflection, from_structure_frame, te_tm_from_circular
from src.ris.state_library import StateLibrary


def _lib():
    f, th, phi = np.linspace(8, 12, 17), np.array([0.0, 15.0, 30.0, 45.0]), np.arange(0.0, 360, 30.0)
    FF, TT, PP = np.meshgrid(f, th, phi, indexing="ij")
    c = np.zeros((*FF.shape, 4, 2, 2), dtype=complex)
    for k in range(4):
        c[:, :, :, k, 0, 0] = 0.9 * np.exp(1j * (np.radians(90 * k - 100) + 0.3 * (FF - 10)))
        c[:, :, :, k, 1, 1] = 0.9 * np.exp(-1j * (np.radians(90 * k - 100) + 0.3 * (FF - 10)))
        c[:, :, :, k, 1, 0] = 0.15
        c[:, :, :, k, 0, 1] = 0.12
    arrays = {
        "f_ghz": f,
        "theta_deg": th,
        "phi_deg": phi,
        "state": np.arange(1, 5),
        "c_raw": c,
        "c_corr": c,
        "uncertainty": np.full(c.shape, 0.02),
    }
    from src.ris.state_library import te_tm_from_structure_circular

    arrays["s_corr"] = te_tm_from_structure_circular(c, phi)
    return StateLibrary(arrays), arrays


def _write_hfss(path, lib, points, phase_off=None, tm_flip=False):
    """Long-format CSV as the user's HFSS export would be, generated from the library (+ optional per-state r_RR phase offset applied through the circular matrix)."""
    rows = []
    for state, f, th in points:
        c = lib.structure_circular(f, th, 0.0, state)[0].copy()
        if phase_off:
            c[0, 0] *= np.exp(1j * math.radians(phase_off.get(state, 0.0)))
        s = te_tm_from_circular(from_structure_frame(c, 0.0))
        for q, (i, j) in {"S_TE_TE": (0, 0), "S_TM_TM": (1, 1), "S_TE_TM": (0, 1), "S_TM_TE": (1, 0)}.items():
            z = s[i, j] * (-1 if tm_flip and (i == 1 or j == 1) else 1)
            rows.append(
                {
                    "state": state,
                    "f_ghz": f,
                    "theta_deg": th,
                    "phi_deg": 0.0,
                    "quantity": q,
                    "hfss_mag_dB": 20 * math.log10(abs(z)),
                    "hfss_phase_deg": math.degrees(np.angle(z)),
                }
            )
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


POINTS = [(k, f, th) for k in (1, 2, 3, 4) for f in (9.0, 10.0, 11.0) for th in (0.0, 30.0)]


def test_hfss_identical_to_the_library_agrees_and_proposes_no_offset(tmp_path):
    lib, _ = _lib()
    _write_hfss(tmp_path / "h.csv", lib, POINTS)
    comps, summary = hc.compare(lib, hc.read_hfss_rows(tmp_path / "h.csv"))
    assert summary["verdict"] == hc.VERDICT_AGREE and summary["phase_diff_deg"]["max"] < 1e-6
    off = hc.proposed_offsets(comps)
    assert all(abs(v["phase_offset_deg"]) < 1e-6 and abs(v["level_offset_db"]) < 1e-6 for v in off.values())


def test_injected_phase_offset_is_found_and_flagged(tmp_path):
    lib, _ = _lib()
    _write_hfss(tmp_path / "h.csv", lib, POINTS, phase_off={2: 12.0, 4: -7.0})
    comps, summary = hc.compare(lib, hc.read_hfss_rows(tmp_path / "h.csv"))
    off = hc.proposed_offsets(comps)
    assert off[2]["phase_offset_deg"] == pytest.approx(12.0, abs=1e-6) and off[4]["phase_offset_deg"] == pytest.approx(
        -7.0, abs=1e-6
    )
    assert off[1]["phase_offset_deg"] == pytest.approx(0.0, abs=1e-6)
    assert summary["verdict"] == hc.VERDICT_DIFFER


def test_tm_sign_shift_is_undone_by_the_documented_option(tmp_path):
    lib, _ = _lib()
    _write_hfss(tmp_path / "h.csv", lib, POINTS, tm_flip=True)
    bad = hc.compare(lib, hc.read_hfss_rows(tmp_path / "h.csv"))[1]
    good = hc.compare(lib, hc.read_hfss_rows(tmp_path / "h.csv", tm_shift_180=True))[1]
    assert good["verdict"] == hc.VERDICT_AGREE and bad["verdict"] == hc.VERDICT_DIFFER


def test_calibration_needs_approval_and_provenance_and_applies_the_offsets(tmp_path):
    lib, arrays = _lib()
    off = {2: {"phase_offset_deg": 10.0, "level_offset_db": -0.5}}
    with pytest.raises(hc.CalibrationNotApproved):
        hc.apply_calibration(arrays, off)
    with pytest.raises(ValueError):
        hc.apply_calibration(arrays, off, approved_by_user=True)  # no provenance
    cal = hc.apply_calibration(arrays, off, approved_by_user=True, source="hfss_cp_cell_export_test.csv")
    assert str(cal["calibration"][0]).startswith("HFSS_CALIBRATED")
    ratio = cal["c_corr"][..., 1, 0, 0] / arrays["c_corr"][..., 1, 0, 0]
    assert np.allclose(ratio, 10 ** (-0.5 / 20) * np.exp(1j * math.radians(10.0)))
    assert np.allclose(cal["c_corr"][..., 0, 0, 0], arrays["c_corr"][..., 0, 0, 0])  # other states untouched
    assert np.allclose(cal["c_raw"], arrays["c_raw"])  # the raw table is never modified
    assert not np.allclose(cal["s_corr"], arrays["s_corr"])
    assert circular_reflection(np.eye(2)).shape == (2, 2)
