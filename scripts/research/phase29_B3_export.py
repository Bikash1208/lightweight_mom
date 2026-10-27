"""Phase 29 design point B3: write the HFSS package `hfss_exports/phase29_B3/` for one cone (predictions first, weight table, 7-port AEDT scripts, README).

    python scripts/research/phase29_B3_export.py --cone a|b          # needs results/phase29/pointB3_<cone>.json (phase29_pointB3.py final)

MODEL_ONLY: every number is a model prediction (inter-port coupling not modelled; gain anchored to the single-hexagon HFSS run only)."""

import argparse
import csv
import json
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts/research"))
OUT = ROOT / "hfss_exports/phase29_B3"


def main():
    import phase29_compact7_export as ce
    import phase29_pointB3 as b3
    from phase29_pointB_sweep import _hpbw

    from src.export import shared7_export as ex
    from src.hologram import seven_hex_cluster as c
    from src.hologram import shared_aperture_multiport as sa

    ap = argparse.ArgumentParser()
    ap.add_argument("--cone", choices=["a", "b"], required=True)
    cone = ap.parse_args().cone
    cfg = b3.CONES[cone]
    res = json.loads((ROOT / f"results/phase29/pointB3_{cone}.json").read_text())
    g = np.array(res["genes"])
    b3._init()
    th, ph = np.arange(0.0, 90.0 + 1e-9, 1.0), np.arange(0.0, 360.0, 1.0)
    d = b3.design(g, cone, n_psi=240)
    F = sa.port_fields(d, b3._S["green"], th, ph)
    OUT.mkdir(parents=True, exist_ok=True)
    rows = res["tuned_rows"]
    cells, feeds = res["cells"], res["feeds_mm"]
    cuts, ports = [], []
    for j, f in enumerate(F, start=1):
        beam = (0.0, 0.0) if j == 1 else (cfg["tilt"], 60.0 * (j - 2))
        m = ce.metrics(*f, th, ph, beam, _hpbw, c._unit)
        m.update({"port": j, "feed_xy_mm": feeds[j - 1], "beam_design_theta_phi_deg": list(beam)})
        ports.append(m)
        cuts += [(f"port{j}_alone", beam[1], *r) for r in ce.cut(*f, th, ph, beam[1])]
    az = np.radians([60.0 * k for k in range(6)])
    special = {
        "sum": np.ones(7),
        "delta_x": np.concatenate([[0.0], np.cos(az)]),
        "delta_y": np.concatenate([[0.0], np.sin(az)]),
    }
    sd = {}
    for name, w in special.items():
        w = w.astype(complex) / np.linalg.norm(w)
        et, ep = sa.combine(F, w)
        sd[name] = ce.metrics(et, ep, th, ph, (0.0, 0.0), _hpbw, c._unit) | {
            "port_power_W": [float(abs(x) ** 2) for x in w],
            "port_phase_deg": [float(np.degrees(np.angle(x))) if abs(x) > 1e-9 else 0.0 for x in w],
        }
        cuts += [(name, 0.0, *r) for r in ce.cut(et, ep, th, ph, 0.0)]
    examples = {}
    for name, (tt, pp) in {
        "broadside": (0.0, 0.0),
        "mid_cone_phi0": (cfg["cone"] / 2, 0.0),
        "cone_edge_phi30": (cfg["cone"], 30.0),
    }.items():
        r = min(rows, key=lambda q: c.angle_between((q["theta_deg"], q["phi_deg"]), (tt, pp)))
        w = np.array([complex(a, b) for a, b in r["weights"]])
        et, ep = sa.combine(F, w)
        m = ce.metrics(et, ep, th, ph, (r["theta_deg"], r["phi_deg"]), _hpbw, c._unit)
        m.update(
            {
                "target_theta_phi_deg": [tt, pp],
                "table_row_theta_phi_deg": [r["theta_deg"], r["phi_deg"]],
                "method": r["method"],
                "active_ports": [int(k) + 1 for k in np.nonzero(np.abs(w) > 1e-9)[0]],
                "port_power_W": [float(abs(x) ** 2) for x in w],
                "port_phase_deg": [float(np.degrees(np.angle(x))) if abs(x) > 1e-9 else 0.0 for x in w],
            }
        )
        examples[name] = m
        cuts += [(f"weighted_{name}", pp, *q) for q in ce.cut(et, ep, th, ph, pp)]
    with open(OUT / "pattern_cuts.csv", "w", newline="", encoding="utf-8") as f:
        w_ = csv.writer(f)
        w_.writerow(["case", "cut_phi_deg", "theta_signed_deg", "rhcp_gain_dbi", "lhcp_gain_dbi"])
        for cs, pc, tsg, gr, gl in cuts:
            w_.writerow([cs, f"{pc:g}", tsg, f"{gr:.2f}", f"{gl:.2f}"])
    with open(OUT / "patterns_3deg.csv", "w", newline="", encoding="utf-8") as f:  # full hemisphere, all 7 ports
        w_ = csv.writer(f)
        w_.writerow(["port", "theta_deg", "phi_deg", "rhcp_gain_dbi", "lhcp_gain_dbi"])
        for j, fj in enumerate(F, start=1):
            gr, gl = sa.gain_db(*fj), ce.lhcp_gain_db(*fj)
            for i in range(0, 91, 3):
                for k in range(0, 360, 6):
                    w_.writerow([j, i, k, f"{gr[i, k]:.2f}", f"{gl[i, k]:.2f}"])
    shutil.copy(ROOT / f"results/phase29/pointB3_{cone}_weights.csv", OUT / "weights_table.csv")
    (OUT / "build_board.py").write_text(ex.aedt_build_script(cells, feeds), encoding="utf-8")
    (OUT / "set_weights.py").write_text(ex.aedt_weights_script(rows), encoding="utf-8")
    with open(OUT / "layout.csv", "w", newline="", encoding="utf-8") as f:
        w_ = csv.writer(f)
        w_.writerow(["x_mm", "y_mm", "circle_D_mm", "slot_L_mm", "slot_W_mm", "slot_angle_deg"])
        for q in cells:
            w_.writerow(
                [
                    f"{q['x_mm']:.3f}",
                    f"{q['y_mm']:.3f}",
                    6.0,
                    f"{q['slot_L_mm']:.2f}",
                    f"{q['slot_W_mm']:.2f}",
                    f"{q['slot_angle_deg']:.3f}",
                ]
            )
    (OUT / "layout.json").write_text(
        json.dumps({"cone": cone, "feeds_mm": feeds, "n_cells": len(cells), "cells": cells}, indent=1), encoding="utf-8"
    )
    pred = {
        "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "SOLVER_SUPERCELL", "INTER_PORT_COUPLING_NOT_MODELLED"],
        "state": "PREDICTION_NOT_HFSS_VALIDATED",
        "cone": cone,
        "config": cfg,
        "frequency_ghz": 10.54,
        "n_cells": len(cells),
        "feeds_mm": feeds,
        "ring_over_lambda0": res["ring_over_lambda0"],
        "gain_reference": "realized gain, RHCP, per 1 W accepted (sum of port powers = 1 W for weighted cases); anchored to the single-hexagon HFSS run only (OFFSET_DB)",
        "per_port_alone": ports,
        "crossover": res["crossover"],
        "sum_diff": sd,
        "weighted_examples": examples,
        "scan_tuned": {k: v for k, v in res["scan_tuned"].items() if k != "method_counts"},
        "coupling_db_estimate": res["coupling_db"],
        "coupling_state": res["coupling_state"],
        "ceiling_gap_db": res["ceiling_gap_db"],
        "s11": "NOT_AVAILABLE_FROM_CURRENT_DATA (the single 150 mm hexagon: -15 dB at 10.5 GHz in HFSS; this board is not simulated)",
        "check_items": [
            "max |S_ij| between feeds (the estimate above is NOT calibrated at this spacing: > -15 dB makes the uncoupled weights unreliable)",
            "each port alone: peak within ~3 deg of the table, gain within 3 dB",
            "weighted examples: peak direction within ~4 deg of the target",
            "sum beam and delta null",
        ],
    }
    (OUT / "predictions.json").write_text(json.dumps(pred, indent=1, default=float), encoding="utf-8")
    (OUT / "README.md").write_text(readme(cone, cfg, res, pred), encoding="utf-8")
    print("wrote", OUT, len(cells), "cells")


def readme(cone, cfg, res, p):
    t = res["scan_tuned"]
    tab = "\n".join(
        f"| P{q['port']} | ({q['feed_xy_mm'][0]:.1f}, {q['feed_xy_mm'][1]:.1f}) | ({q['beam_design_theta_phi_deg'][0]:g}, {q['beam_design_theta_phi_deg'][1]:g}) | {q['peak_rhcp_gain_dbi']:.1f} | "
        f"({q['peak_theta_phi_deg'][0]:.0f}, {q['peak_theta_phi_deg'][1]:.0f}) | {q['hpbw_deg']:.0f} | {q['sll_rel_db_outside_1p5_hpbw']:.1f} | {q['axial_ratio_at_peak_db']:.1f} |"
        for q in p["per_port_alone"]
    )
    ex = "\n".join(
        f"| {k} | {v['active_ports']} | {v['peak_rhcp_gain_dbi']:.1f} | ({v['peak_theta_phi_deg'][0]:.0f}, {v['peak_theta_phi_deg'][1]:.0f}) | {v['pointing_error_deg']:.1f} | {v['hpbw_deg']:.0f} |"
        for k, v in p["weighted_examples"].items()
    )
    sd = p["sum_diff"]
    return f"""# phase29_B3: one shared 150 mm hexagon, 7 monopole ports, continuous scanning by port weights (cone {cone}: +-{cfg["cone"]:g} deg)

One fixed passive hologram on the 150 mm-circumradius hexagon (vertex at +/-y), RT/duroid 5880 1.575 mm, 10.54 GHz, 7 monopoles (centre + ring of 6 at {res["ring_mm"]:.1f} mm = {res["ring_over_lambda0"]:.2f} lambda0),
{len(res["cells"])} slotted-circle cells (D 6 mm) on a 6.4 mm lattice. Port 1: broadside beam; ports 2..7: tilt {cfg["tilt"]:g} deg, azimuth 0, 60, ..., 300 deg. Synthesis: Frey et al. (IEEE OJAP 2023, 10.1109/OJAP.2023.3261903) angle-dependent
superposition in the repo cell model + GA over feed positions / hologram offsets / defocus. **MODEL_ONLY / NOT_HFSS_VALIDATED**: all numbers are predictions; inter-port coupling is NOT modelled (its estimate is not calibrated at this spacing); gain is anchored
only to the single-hexagon HFSS run (23.03 dBi); S11 is NOT_AVAILABLE_FROM_CURRENT_DATA.

## Files
`predictions.json` (all numbers below), `patterns_3deg.csv` (RHCP/LHCP gain of ALL 7 ports, full hemisphere, 3 x 6 deg), `pattern_cuts.csv` (principal cuts: each port alone, sum, delta_x, delta_y, 3 weighted examples), `weights_table.csv` (direction -> power W and phase deg of the 7 ports, <= 3 active),
`build_board.py` (AEDT build), `set_weights.py` (sets the sources of one table row, edit THETA_DEG/PHI_DEG), `layout.csv/json`.

## Run order
1. AEDT: Tools > Run Script > `build_board.py` (several minutes, ~{len(res["cells"])} slotted circles). Check the messages; the project is saved as `phase29_B3.aedt`.
2. Validate, Analyze Setup1 (10.54 GHz single frequency; adaptive passes <= 15, max delta S 0.02). **One solve gives all 7 port patterns**: the fields of every excitation are in the one solution; do NOT re-solve per port.
3. Read the S-matrix first (S11..S77 and all |S_ij|). If a neighbouring |S_ij| is above about -15 dB the uncoupled weights are unreliable.
4. Port k alone: Results > Edit Sources: Port k = 1 W, all other ports 0 W (terminated, 50 ohm); Radiation > Far Fields (infinite sphere 3D, theta 0-180, phi -180..180, 1 deg): Realized Gain RHCP/LHCP, Axial Ratio. Compare with `patterns_3deg.csv` / `pattern_cuts.csv`.
5. Sum / delta / weighted: set the sources (below), no new solve, same far-field export.

## HFSS settings
Driven modal; lumped 50 ohm ports on the coax annulus (r_pin 0.635 .. r_out 2.05 mm) at z = 0; radiation boundary on the padded air box (15 mm sides/bottom, 30 mm top; assigned to the Air object, as in the Phase 26 board that ran);
perfect E on ground and patches; frequency 10.54 GHz, no sweep needed (`ADD_SWEEP = True` at the top of the script adds 10.0-11.1 GHz); far-field sphere 1 deg; lambda refinement on.

## How to apply the port weights (Edit Sources)
Results > Edit Sources... > type Power: for Port1..Port7 set Magnitude = `P<k>_power_W` (W, sum = 1 W) and Phase = `P<k>_phase_deg` (deg) from the row of `weights_table.csv` (nearest direction). The weights are field weights, |w|^2 = incident power, phase = angle(w), exp(+j w t); the
phase reference is the identical port plane of all 7 ports, so a common offset cancels. If beams come out mirrored in phi, negate all phases. `set_weights.py` does the same (EditSources is NOT verified against an AEDT version: wrapped, it prints the values if the call fails).
Sum beam: all seven ports 1/7 W, phase 0. delta_x / delta_y: powers proportional to cos^2 / sin^2 of the ring azimuth with the sign of cos / sin as phase 0 / 180 deg (rows in `predictions.json` `sum_diff`).

## Predictions (MODEL_ONLY, 10.54 GHz)
| Port | Feed (mm) | Design beam (theta, phi) | Peak gain (dBi) | Peak (theta, phi) | HPBW (deg) | SLL (dB) | AR (dB) |
|---|---|---|---|---|---|---|---|
{tab}

Crossover (single-port level at the midpoint of neighbouring beams, dB re own peak): centre-outer {p["crossover"]["centre_outer"]["mean_db"]:.1f}, outer-outer {p["crossover"]["outer_outer"]["mean_db"]:.1f}.
Sum beam (all ports in phase): {sd["sum"]["peak_rhcp_gain_dbi"]:.1f} dBi at ({sd["sum"]["peak_theta_phi_deg"][0]:.0f}, {sd["sum"]["peak_theta_phi_deg"][1]:.0f}), HPBW {sd["sum"]["hpbw_deg"]:.0f} deg, SLL {sd["sum"]["sll_rel_db_outside_1p5_hpbw"]:.1f} dB;
delta_x peak {sd["delta_x"]["peak_rhcp_gain_dbi"]:.1f} dBi, delta_y {sd["delta_y"]["peak_rhcp_gain_dbi"]:.1f} dBi.

| Weighted example | Active ports | Peak gain (dBi) | Peak (theta, phi) | Pointing error (deg) | HPBW (deg) |
|---|---|---|---|---|---|
{ex}

Scan over the whole cone (tuned weights, <= 3 ports): worst pointing error {t["max_pointing_error_deg"]:.1f} deg (target 2), median {t["median_pointing_error_deg"]:.1f}, gain dip {t["max_gain_dip_db"]:.1f} dB (target 3), minimum peak gain {t["min_gain_dbi"]:.1f} dBi.
Estimated coupling (raw capture-width estimate, NOT calibrated at this spacing): see `predictions.json`. Gap to the given DoF ceiling {p["ceiling_gap_db"]["ceiling_dbi_as_given"]} dBi: best port {p["ceiling_gap_db"]["gap_best_port_db"]:.1f} dB, mean {p["ceiling_gap_db"]["gap_mean_port_db"]:.1f} dB.

What would falsify the model: a port alone more than ~5 deg off its predicted peak or more than 4 dB from the table; a weighted example more than ~6 deg off its target; neighbouring-port coupling above -10 dB. With coupling the weights should be re-tuned from the measured single-port patterns (the model is linear).
"""


if __name__ == "__main__":
    main()
