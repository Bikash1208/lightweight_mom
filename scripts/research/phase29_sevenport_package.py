"""Phase 29: 7-port hexagonal board package `hfss_exports/phase29_sevenport/` (MODEL_ONLY + HFSS-checked margin), same method/settings as the two-feed board: tensor model (n_psi 360), 150 mm hexagon, one fixed port per tile,
10.54 GHz, STEP files (twofeed_step.build) + AEDT script + port/excitation table + predictions (raw and calibrated with phase29_calibration). Design: GA genes `tensor_ga_n7_pass2.json`, outer tilt 9.6 deg (tensor_summary R2).
    python scripts/research/phase29_sevenport_package.py"""

import csv
import json
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts/research"))
OUT = ROOT / "hfss_exports/phase29_sevenport"
SYNC = Path(r"C:\hmsa_sync\data\phase29_sevenport")
TILT = 9.6


def main():
    import phase29_calibration as cal
    import phase29_compact7_export as ce
    import phase29_tensor_ladder as L
    import phase29_twofeed_step as ts
    from phase29_pointB_sweep import _hpbw

    from src.export import shared7_export as ex
    from src.hologram import seven_hex_cluster as c
    from src.hologram import shared_aperture_multiport as sa

    L.init()
    g = np.array(json.loads((ROOT / "results/phase29/tensor_ga_n7_pass2.json").read_text())["best_genes"])
    d = L.design_for(g, "n7", TILT, 150.0, n_psi=360)
    th, ph = np.arange(0.0, 90.0 + 1e-9, 1.0), np.arange(0.0, 360.0, 1.0)
    F = sa.port_fields(d, L._S["green"], th, ph)
    beams = [(float(np.degrees(b[0])), float(np.degrees(b[1]))) for b in d["beams"]]
    OUT.mkdir(parents=True, exist_ok=True)
    cells = [
        (float(x * 1e3), float(y * 1e3), float(a), float(b), float(np.degrees(t_)))
        for x, y, a, b, t_ in zip(d["x_m"], d["y_m"], d["slot_L_mm"], d["slot_W_mm"], d["axis_rad"], strict=True)
    ]
    feeds = [(f[0] * 1e3, f[1] * 1e3) for f in d["feeds_m"]]
    cdict = [dict(zip(("x_mm", "y_mm", "slot_L_mm", "slot_W_mm", "slot_angle_deg"), r, strict=True)) for r in cells]
    (OUT / "build_board.py").write_text(ex.aedt_build_script(cdict, feeds, name="phase29_sevenport"), encoding="utf-8")
    ports = []
    for j, f in enumerate(F, start=1):
        m = ce.metrics(*f, th, ph, beams[j - 1], _hpbw, c._unit)
        gr, sl = cal.calibrate(m["peak_rhcp_gain_dbi"], m["sll_rel_db_outside_1p5_hpbw"])
        ports.append(
            {
                "port": j,
                "tile": "centre" if j == 1 else f"ring {j - 1}",
                "feed_x_mm": round(feeds[j - 1][0], 3),
                "feed_y_mm": round(feeds[j - 1][1], 3),
                "design_beam_theta_deg": round(beams[j - 1][0], 2),
                "design_beam_phi_deg": round(beams[j - 1][1], 1),
                "pred_peak_theta_phi_deg": [round(x, 1) for x in m["peak_theta_phi_deg"]],
                "hpbw_deg": round(m["hpbw_deg"], 1),
                "gain_raw_dbi": round(m["peak_rhcp_gain_dbi"], 1),
                "sll_raw_db": round(m["sll_rel_db_outside_1p5_hpbw"], 1),
                "gain_calibrated_dbi": round(gr, 1),
                "sll_calibrated_db": round(sl, 1),
            }
        )
    with open(OUT / "ports_and_excitations.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(list(ports[0]) + [f"excite_P{k}_W_deg" for k in range(1, 8)])
        for p in ports:
            w.writerow(list(p.values()) + [("1.0 W / 0 deg" if k == p["port"] else "0 (matched)") for k in range(1, 8)])
    wt = np.ones(7, dtype=complex) / np.sqrt(7)
    mt = ce.metrics(*sa.combine(F, wt), th, ph, (0.0, 0.0), _hpbw, c._unit)
    pred = {
        "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "INTER_PORT_COUPLING_NOT_MODELLED", cal.LABEL],
        "frequency_ghz": 10.54,
        "design": f"tensor model, GA genes tensor_ga_n7_pass2.json, outer tilt {TILT} deg, n_cells {len(cells)}",
        "margin": {"gain_db": cal.GAIN_OFFSET_DB, "sll_db": cal.SLL_OFFSET_DB, "source": "two-feed board HFSS at 10.54 GHz (docs/STATUS.md)"},
        "ports": ports,
        "all_in_phase_1_over_7_W_each": {k: m for k, m in mt.items() if k in ("peak_rhcp_gain_dbi", "peak_theta_phi_deg", "sll_rel_db_outside_1p5_hpbw")},
        "beam_direction": "raw only (the two-feed margin has no pointing offset)",
    }
    (OUT / "predictions.json").write_text(json.dumps(pred, indent=1, default=float), encoding="utf-8")
    SYNC.mkdir(parents=True, exist_ok=True)
    stepf = OUT / "sevenport_board.step"
    ts.build(cells, feeds, 270.0, 310.0, stepf)
    names, ok = ts.check(stepf, len(cells), len(feeds))
    shutil.copy2(stepf, SYNC / stepf.name)
    shutil.copy2(OUT / "ports_and_excitations.csv", SYNC / "ports_and_excitations.csv")
    (OUT / "README.md").write_text(
        "# phase29_sevenport (MODEL_ONLY, HFSS-checked margin)\n\nSame method and settings as `hfss_exports/phase29_twofeed/`: 150 mm hexagon (270 x 310 mm board), 1.575 mm RT5880, 6 mm circles, pins r 0.635 mm to z = 7.575 mm, "
        "one fixed port per tile (P1 centre, P2..P7 ring), 10.54 GHz single frequency, Driven Modal, one solve then Edit Sources per port.\n\n"
        "- `sevenport_board.step` (gitignored; copy in `C:\hmsa_sync\data\phase29_sevenport\`): bodies Substrate, Ground, Metal_pattern, Pin_1..7, Port_1..7. Manual setup: the STEP section of `../phase29_twofeed/README.md` (feeds/ports 1..7 instead of 2).\n"
        "- `build_board.py`: AEDT script, same geometry.\n- `ports_and_excitations.csv`: port positions, design beam, predicted gain / peak direction / SLL raw and calibrated (gain -2.5 dB, SLL +3.5 dB, from the two-feed HFSS), and the excitation per case (port k alone 1 W).\n"
        "- `predictions.json`: same plus the all-in-phase case. Inter-port coupling not modelled.\n",
        encoding="utf-8",
    )
    print("cells", len(cells), "check_ok", ok, {k: v for k, v in names.items() if k and not k.startswith(("Pin", "Port"))})
    for p in ports:
        print(p["port"], p["design_beam_theta_deg"], p["design_beam_phi_deg"], p["pred_peak_theta_phi_deg"], p["gain_raw_dbi"], p["sll_raw_db"], p["gain_calibrated_dbi"], p["sll_calibrated_db"])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
