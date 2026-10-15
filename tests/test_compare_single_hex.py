"""Phase 26: the single-hexagon HFSS comparison script on synthetic exports with known answers."""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

_p = Path(__file__).resolve().parents[1] / "scripts/research/phase26_compare_single_hex.py"
_s = importlib.util.spec_from_file_location("p26_compare", _p)
cmp_ = importlib.util.module_from_spec(_s)
sys.modules["p26_compare"] = cmp_
_s.loader.exec_module(cmp_)


def _write_gain(path, peak_dbi=22.0, lhcp_db=-20.0, step=2.0):
    th = np.arange(0.0, 180.0 + 1e-9, step)
    ph = np.arange(-180.0, 180.0 + 1e-9, step)
    lines = ['"Freq [GHz]","Phi [deg]","Theta [deg]","dB(RealizedGainLHCP) []","dB(RealizedGainRHCP) []"']
    for p in ph:
        for t in th:
            d = cmp_._angle_to(t, np.mod(p, 360.0), 30.0, 225.0)
            r = peak_dbi - 12.0 * (d / 8.0) ** 2  # Gaussian-ish beam, -12 dB at 8 deg
            lines.append(f"10.54,{p:g},{t:g},{peak_dbi + lhcp_db - 12.0 * (d / 8.0) ** 2:.4f},{max(r, -30.0):.4f}")
    Path(path).write_text("\n".join(lines) + "\n")


def test_peak_direction_xpd_and_budget_of_a_synthetic_beam(tmp_path):
    f = tmp_path / "gain.csv"
    _write_gain(f)
    g = cmp_.read_gain_csv(f)
    a, _ = cmp_.analyse(g)
    assert a["peak_rhcp_gain_dbi"] == pytest.approx(22.0, abs=1e-3)
    assert a["pointing_error_deg"] <= 1.5
    assert a["xpd_at_peak_db"] == pytest.approx(20.0, abs=0.1)
    b = a["budget_of_radiated"]
    parts = b["main_beam_le_8deg_total"] + b["near_sidelobes_8_to_20deg"] + b["far_sidelobes_theta_le_70"] + b["grazing_70_to_90"] + b["back_hemisphere_gt_90"]
    assert parts == pytest.approx(1.0, abs=1e-9)
    assert b["lhcp_total_all_angles"] == pytest.approx(0.01 / 1.01, rel=0.05)


def test_compare_states_and_missing_inputs(tmp_path):
    f = tmp_path / "gain.csv"
    _write_gain(f, peak_dbi=17.0)
    args = type("A", (), {"gain": f, "s11": None, "port": 1, "params": None, "ar": None, "ez": None, "pin": "connector", "pred": str(cmp_.PRED)})()
    rep = cmp_.compare(args)
    assert rep["checks"]["peak_gain"]["state"] == "MODEL_FALSIFIED_PEAK_GAIN_BELOW_1BIT_BOARD"
    assert rep["s_parameter"] == "NOT_PROVIDED" and rep["pin_model"] == "connector"
    txt = cmp_.markdown(rep)
    assert "MODEL_FALSIFIED" in txt


def test_params_parser_takes_last_column(tmp_path):
    f = tmp_path / "p.txt"
    f.write_text("Quantity   Freq   Value   Freq   Value\nPeak Realized Gain      4.096466        97.157858\nTotal Efficiency      0.70      0.977420\n")
    p = cmp_.read_params(f)
    assert p["Peak Realized Gain"] == pytest.approx(97.157858)
    assert p["Total Efficiency"] == pytest.approx(0.97742)


def test_decay_fit_recovers_alpha(tmp_path):
    rho = np.linspace(10, 140, 60)
    e = rho**-0.5 * np.exp(-3.0 * rho * 1e-3)
    f = tmp_path / "ez.csv"
    f.write_text("Distance [mm],Mag_E [V_per_m]\n" + "\n".join(f"{r:.3f},{v:.6e}" for r, v in zip(rho, e)))
    assert cmp_.fit_decay(f)["alpha_np_per_m"] == pytest.approx(3.0, rel=1e-3)
