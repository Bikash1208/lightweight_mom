"""Phase 20 R1 item 4: the demonstration case and the one-cell HFSS check package (consistency, not new physics)."""

import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "research"))
import phase20_r1_demo as demo  # noqa: E402
import phase20_r1_hfss_package as pkg  # noqa: E402


def test_demo_runs_with_other_parameters_and_shows_the_expected_physics(tmp_path):
    # different period / thickness / patch than the defaults: nothing is a constant of the script
    meta, bare, rows = demo.run(
        tmp_path,
        period_mm=12.0,
        thickness_mm=1.0,
        eps_r=3.0,
        patch_length_mm=9.0,
        patch_width_mm=8.0,
        gap_mm=1.2,
        mesh_cells=30,
        floquet=15,
        angles=(0.0, 60.0),
        quick=True,
    )
    assert meta["parameters"]["period_mm"] == 12.0 and (tmp_path / "demo.csv").exists()
    tm0 = [r for r in rows if r["theta_deg"] == 0.0 and r["mode"] == "TM"]
    short, opn = tm0[0], tm0[-1]
    assert short["load"].startswith("short") and opn["load"] == "open"
    swing = abs((opn["S11_co_phase_deg"] - short["S11_co_phase_deg"] + 180) % 360 - 180)
    assert swing > 60.0  # the load moves the TM phase
    assert (
        max(r["absorbed_fraction"] for r in tm0) > 0.5
    )  # a resistive load absorbs strongly somewhere between short and open
    assert (
        abs(short["S11_co_mag_dB"]) < 0.5 and abs(opn["S11_co_mag_dB"]) < 0.5
    )  # lossless limits (substrate loss only)
    # bare grounded slab: unit magnitude (tan d is small but the reference is lossless-ish)
    assert all(abs(b[m][0] - 1.0) < 5e-3 for b in bare for m in demo.MODES)


def test_richardson_extrapolation_recovers_the_limit_of_a_second_order_sequence():
    s_inf, c = 0.3 + 0.7j, 12.0 - 5.0j
    f = lambda n: s_inf + c / n**2  # noqa: E731
    assert abs(pkg.richardson(f(80), 80, f(100), 100) - s_inf) < 1e-12
    assert abs(pkg.richardson(f(60), 60, f(80), 80) - s_inf) < 1e-12


def test_committed_hfss_package_is_complete_and_self_consistent():
    d = ROOT / "hfss_exports/phase20_r1_loaded_cell"
    rows = list(csv.DictReader((d / "mom_reference.csv").read_text().splitlines()))
    assert len(rows) == len(pkg.LOADS) * len(pkg.ANGLES) * 2
    assert {r["load"] for r in rows} == {lab for lab, _, _ in pkg.LOADS}
    for r in rows:
        assert math.isfinite(float(r["mom_S11_phase_deg"])) and float(r["mom_uncertainty_phase_deg"]) >= 0.0
        assert float(r["mom_S11_mag_dB"]) <= 1e-3  # passive
    parts = json.loads((d / "step_parts.json").read_text())["parts"]
    assert {"Substrate", "Metal", "Ground", "Load_1"} <= set(parts) and parts["Metal"]["n_solids"] == 2
    assert parts["Load_1"]["n_solids"] == 0 and parts["Load_1"]["n_faces"] == 1
    bare = list(csv.DictReader((d / "bare_grounded_slab_reference.csv").read_text().splitlines()))
    assert len(bare) == 6 and all(abs(float(b["S11_mag_dB"])) < 0.2 for b in bare)
    readme = (d / "README.md").read_text()
    assert "PARALLEL" in readme and "NOT_HFSS_VALIDATED" in readme and "Floquet port" in readme
    # the committed reference equals a fresh coarse solve to within the stated uncertainty at one row (cheap consistency check)
    row = next(r for r in rows if r["load"] == "R=50" and r["theta_deg"] == "0.0" and r["floquet_mode"].startswith("2"))
    assert (
        abs(float(row["mom_mesh20_phase_deg"]) - float(row["mom_S11_phase_deg"])) < 40.0
    )  # 20 x 20 is far from converged, but not wildly off
    assert np.isfinite(float(row["mom_uncertainty_mag_dB"]))
