"""Phase 20 R2 item 4: the demo's committed outputs are internally consistent, and a quick live solve reproduces them."""

import csv
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phase20_r2"


def _demo():
    spec = importlib.util.spec_from_file_location(
        "phase20_r2_demo", ROOT / "scripts" / "research" / "phase20_r2_demo.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_committed_sweep_reproduces_with_a_live_solve_and_has_the_physical_structure():
    mod = _demo()
    with open(OUT / "demo_sweep.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 9 * 3 * 2
    row = next(r for r in rows if float(r["f_ghz"]) == 10.0 and float(r["theta_deg"]) == 30.0 and r["mode"] == "TM")
    from src.ris.diodes import smp1345_off, smp1345_on

    on, off = mod.solve([smp1345_on(), smp1345_off()], 10e9, 30.0, 40)
    assert np.degrees(np.angle(on.s11[1, 1])) == pytest.approx(float(row["on_phase_deg"]), abs=1e-6)
    assert np.degrees(np.angle(off.s11[1, 1])) == pytest.approx(float(row["off_phase_deg"]), abs=1e-6)
    te = [r for r in rows if r["mode"] == "TE" and float(r["f_ghz"]) < 12]
    assert (
        max(abs(float(r["on_off_phase_deg"])) for r in te) < 1.0
    )  # a diode-sized bridge: TE does not see the diode state
    assert max(float(r["on_loss"]) for r in rows) < 0.05 and max(float(r["off_loss"]) for r in rows) < 0.05


def test_committed_sensitivity_and_bias():
    sens = json.loads((OUT / "demo_sensitivity.json").read_text())
    assert abs(sens["on:L_s x1.2|theta0|TM"]["d_phase_deg"]) < 3.0  # ON is robust
    assert abs(sens["off:C_T x1.2|theta0|TM"]["d_phase_deg"]) > 10.0  # OFF is not
    bud = json.loads((OUT / "bias_budget.json").read_text())["rows"]
    r = next(b for b in bud if b["n_cells"] == 1000 and b["fraction_on"] == 0.5 and b["i_forward_a"] == 10e-3)
    assert r["total_current_a"] == pytest.approx(5.0) and r["dc_power_w"] == pytest.approx(5.0 * 0.89)
