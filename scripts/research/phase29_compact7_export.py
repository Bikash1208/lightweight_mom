"""Phase 29 compact 7-tile board: write the HFSS package `hfss_exports/phase29_compact7/` (predictions first, weight table CSV, AEDT scripts, README).

    python scripts/research/phase29_compact7_export.py            # reads results/phase29/compact7_tune.json (written by phase29_compact7_tune.py)

MODEL_ONLY: every number in the package is a model prediction; inter-tile coupling is not modelled; gain anchored only to the single-hexagon HFSS run."""

import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts/research"))
OUT = ROOT / "hfss_exports/phase29_compact7"
EXAMPLES = {"broadside": (0.0, 0.0), "theta15_phi0": (15.0, 0.0), "theta30_phi30": (30.0, 30.0)}


def lhcp_gain_db(et, ep):
    from src.hologram import tensor_cp_design as t
    from src.pattern.vector_aperture_pattern import copol_component

    g = np.abs(copol_component(et, ep, "lhcp")) ** 2 * 4 * np.pi * 10 ** (t.OFFSET_DB / 10)
    return 10 * np.log10(np.maximum(g, 1e-30))


def cut(et, ep, theta, phi, phi_cut):
    """RHCP / LHCP gain (dBi) along the plane phi_cut / phi_cut + 180, signed theta -90..90 (negative = the opposite half-plane)."""
    from src.hologram import shared_aperture_multiport as sa

    g_r, g_l = sa.gain_db(et, ep), lhcp_gain_db(et, ep)
    i0 = int(round(phi_cut / (phi[1] - phi[0]))) % phi.size
    i1 = int(round((phi_cut + 180.0) / (phi[1] - phi[0]))) % phi.size
    out = []
    for th in range(-90, 91):
        i, ip = abs(th), (i0 if th >= 0 else i1)
        out.append((th, float(g_r[i, ip]), float(g_l[i, ip])))
    return out


def metrics(et, ep, theta, phi, beam, hpbw_fn, unit_fn):
    from src.hologram import seven_hex_cluster as c
    from src.hologram import shared_aperture_multiport as sa
    from src.hologram import tensor_cp_design as t

    g = sa.gain_db(et, ep)
    lb = c.lobes(g, theta, phi, within_db=3.0)
    pk = (lb[0][1], lb[0][2])
    hp = hpbw_fn(g, theta, phi)
    tg, pg = np.meshgrid(theta, phi, indexing="ij")
    ang = np.degrees(np.arccos(np.clip(np.sum(unit_fn(tg, pg) * unit_fn(*pk), axis=-1), -1, 1)))
    pm = t.pol_metrics(et, ep, theta, phi, beam)
    return {
        "peak_rhcp_gain_dbi": float(lb[0][0]),
        "peak_theta_phi_deg": [float(pk[0]), float(pk[1])],
        "pointing_error_deg": c.angle_between(pk, beam),
        "hpbw_deg": hp,
        "sll_rel_db_outside_1p5_hpbw": float(g[ang > 1.5 * hp].max() - g.max()),
        "n_lobes_within_3db": len(lb),
        "axial_ratio_at_peak_db": pm["ar_peak_db"],
        "xpd_at_peak_db": pm["xpd_db"],
    }


def main():
    from phase29_compact7_tune import tile_fields

    from src.export import compact7_export as ex
    from src.hologram import seven_hex_cluster as c
    from src.hologram import shared_aperture_multiport as sa
    from src.hologram import tensor_cp_design as t
    from src.hologram.tensor_leaky_cp import SupercellTable

    tune = json.loads((ROOT / "results/phase29/compact7_tune.json").read_text())
    R = float(tune["chosen_R_mm"])
    rows = tune["sizes"][f"R{R:g}"]["rows"]
    summ, base = tune["sizes"][f"R{R:g}"]["tuned"], tune["sizes"][f"R{R:g}"]["mrc3_untuned"]
    table, green = SupercellTable.from_files(t.TABLE_GLOB), c.green_function()
    th, ph = np.arange(0.0, 90.0 + 1e-9, 1.0), np.arange(0.0, 360.0, 1.0)
    tiles, F = tile_fields(table, green, R, th, ph)
    cells = []
    for j, tl in enumerate(tiles, start=1):
        d = tl["design"]
        for i in range(tl["x_mm"].size):
            cells.append(
                {
                    "x_mm": float(tl["x_mm"][i]),
                    "y_mm": float(tl["y_mm"][i]),
                    "slot_L_mm": float(d["slot_L_mm"][i]),
                    "slot_W_mm": float(d["slot_W_mm"][i]),
                    "slot_angle_deg": float(np.degrees(d["slot_angle_rad"][i])),
                    "port": j,
                }
            )
    OUT.mkdir(parents=True, exist_ok=True)
    cen = c.tile_centres_mm(R)
    ports = []
    cuts = []
    from phase29_pointB_sweep import _hpbw as hpbw_fn

    for j, f in enumerate(F, start=1):
        beam = (0.0, 0.0) if j == 1 else (30.0, 60.0 * (j - 2))
        m = metrics(*f, th, ph, beam, hpbw_fn, c._unit)
        m.update(
            {
                "port": j,
                "feed_xy_mm": [cen[j - 1][0], cen[j - 1][1]],
                "beam_design_theta_phi_deg": list(beam),
                "n_cells": int(tiles[j - 1]["x_mm"].size),
                "eta_conv": float(tiles[j - 1]["eta_conv"]),
            }
        )
        ports.append(m)
        cuts += [(f"port{j}_alone", beam[1], *r) for r in cut(*f, th, ph, beam[1])]
    examples = {}
    for name, (tt, pp) in EXAMPLES.items():
        r = min(rows, key=lambda q: c.angle_between((q["theta_deg"], q["phi_deg"]), (tt, pp)))
        w = np.array([complex(a, b) for a, b in r["weights"]])
        et, ep = sa.combine(F, w)
        m = metrics(et, ep, th, ph, (tt, pp), hpbw_fn, c._unit)
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
        cuts += [(f"weighted_{name}", pp, *q) for q in cut(et, ep, th, ph, pp)]
    with open(OUT / "pattern_cuts.csv", "w", newline="", encoding="utf-8") as f:
        w_ = csv.writer(f)
        w_.writerow(["case", "cut_phi_deg", "theta_signed_deg", "rhcp_gain_dbi", "lhcp_gain_dbi"])
        for cs, pc, tsg, gr, gl in cuts:
            w_.writerow([cs, f"{pc:g}", tsg, f"{gr:.2f}", f"{gl:.2f}"])
    with open(OUT / "weights_table.csv", "w", newline="", encoding="utf-8") as f:
        w_ = csv.writer(f)
        w_.writerow(
            ["theta_deg", "phi_deg"]
            + [f"P{k}_{q}" for k in range(1, 8) for q in ("power_W", "phase_deg")]
            + [
                "method",
                "active_ports",
                "pred_peak_gain_dbi",
                "pred_peak_theta_deg",
                "pred_peak_phi_deg",
                "pred_pointing_error_deg",
            ]
        )
        for r in rows:
            w = np.array([complex(a, b) for a, b in r["weights"]])
            cells_ = []
            for x in w:
                cells_ += [f"{abs(x) ** 2:.5f}", f"{np.degrees(np.angle(x)) if abs(x) > 1e-9 else 0.0:.2f}"]
            w_.writerow(
                [f"{r['theta_deg']:g}", f"{r['phi_deg']:g}"]
                + cells_
                + [
                    r["method"],
                    int(np.sum(np.abs(w) > 1e-9)),
                    f"{r['peak_gain_dbi']:.2f}",
                    f"{r['peak_theta_phi_deg'][0]:.1f}",
                    f"{r['peak_theta_phi_deg'][1]:.1f}",
                    f"{r['pointing_error_deg']:.2f}",
                ]
            )
    (OUT / "build_board.py").write_text(ex.aedt_build_script(cells, R), encoding="utf-8")
    (OUT / "set_weights.py").write_text(ex.aedt_weights_script(rows), encoding="utf-8")
    ex.write_layout(cells, R, OUT)
    pred = {
        "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "SOLVER_SUPERCELL", "INTER_TILE_COUPLING_NOT_MODELLED"],
        "state": "PREDICTION_NOT_HFSS_VALIDATED",
        "frequency_ghz": 10.54,
        "tile_radius_mm": R,
        "n_cells": len(cells),
        "board_mm": ex.board_size_mm(R),
        "gain_reference": "realized gain, RHCP, referred to the accepted power (1 W per port alone; sum of the port powers = 1 W for the weighted cases); anchored to the single-hexagon HFSS run only (OFFSET_DB)",
        "per_port_alone": ports,
        "weighted_examples": examples,
        "scan_tuned_inner_le_30": summ["inner_le_30"],
        "scan_tuned_band_30_45": summ["band_30_45"],
        "scan_mrc3_untuned_for_reference": {k: base[k] for k in ("inner_le_30", "band_30_45")},
        "criterion": tune["criterion"],
        "pass_relaxed_inner": summ["pass_relaxed_inner"],
        "s11": "NOT_AVAILABLE_FROM_CURRENT_DATA (the single 150 mm hexagon measured -15 dB at 10.5 GHz in HFSS; this board is not simulated)",
        "port_coupling": "NOT_AVAILABLE_FROM_CURRENT_DATA (not modelled; the HFSS S-matrix is the first thing to read)",
        "check_items": [
            "each port alone: peak within ~5 deg of the design direction, HPBW 20-30 deg, gain within 3 dB of the table",
            "max |S_ij| between neighbouring ports (> -15 dB would make the uncoupled weights unreliable)",
            "weighted examples: peak direction within ~8 deg of the target, peak gain within 3 dB of the prediction",
            "AR at the peak <= 3 dB",
        ],
    }
    (OUT / "predictions.json").write_text(json.dumps(pred, indent=1, default=float), encoding="utf-8")
    (OUT / "README.md").write_text(readme(R, len(cells), pred), encoding="utf-8")
    print("wrote", OUT, "cells", len(cells), "R", R)


def readme(R, n, p):
    pp = p["per_port_alone"]
    ex = p["weighted_examples"]
    inn, band = p["scan_tuned_inner_le_30"], p["scan_tuned_band_30_45"]
    tab = "\n".join(
        f"| P{q['port']} | ({q['feed_xy_mm'][0]:.1f}, {q['feed_xy_mm'][1]:.1f}) | ({q['beam_design_theta_phi_deg'][0]:g}, {q['beam_design_theta_phi_deg'][1]:g}) | "
        f"{q['peak_rhcp_gain_dbi']:.1f} | ({q['peak_theta_phi_deg'][0]:.0f}, {q['peak_theta_phi_deg'][1]:.0f}) | {q['hpbw_deg']:.0f} | {q['sll_rel_db_outside_1p5_hpbw']:.1f} | {q['axial_ratio_at_peak_db']:.1f} |"
        for q in pp
    )
    etab = "\n".join(
        f"| {k} | {v['active_ports']} | {v['peak_rhcp_gain_dbi']:.1f} | ({v['peak_theta_phi_deg'][0]:.0f}, {v['peak_theta_phi_deg'][1]:.0f}) | "
        f"{v['pointing_error_deg']:.1f} | {v['hpbw_deg']:.0f} | {v['axial_ratio_at_peak_db']:.1f} |"
        for k, v in ex.items()
    )
    return f"""# phase29_compact7: seven hexagonal tiles side by side, one monopole each (Phase 29, compact option)

Seven R {R:g} mm hexagons (vertex at +/-y, flat edges shared, no gap) around a centre tile, one centre monopole + lumped 50 ohm port per tile (7 ports), rotating-slot
tensor-CP cells (circle D 6 mm, slot L x W and angle per cell) on one shared 6.4 mm lattice, RT/duroid 5880 1.575 mm, 10.54 GHz, {n} cells. Port 1 = centre tile (broadside beam);
ports 2..7 = outer tiles at azimuth 0, 60, ..., 300 deg (centre of the tile seen from the board centre), each beam at theta 30 deg pointing outward. Continuous scanning = complex weights on
2-3 ports at a time (`weights_table.csv`), theta 0-45 deg in 3 deg steps, phi in 15 deg steps. **MODEL_ONLY / NOT_HFSS_VALIDATED**: all numbers below are predictions; inter-tile
coupling is NOT modelled; gain is anchored only to the single-hexagon HFSS run (23.03 dBi at 150 mm), so tile size and tilt are extrapolations; S11 and port coupling are
NOT_AVAILABLE_FROM_CURRENT_DATA.

## Files
| File | Content |
|---|---|
| `predictions.json` | everything below, machine readable (per port alone, three weighted examples, scan statistics, check items) |
| `pattern_cuts.csv` | RHCP/LHCP gain (dBi) vs signed theta (-90..90 deg) in the plane of each port's beam (`portK_alone`) and of the three weighted examples (`weighted_*`) |
| `weights_table.csv` | per direction: for ports 1..7 the incident power (W, sum = 1 W) and phase (deg), method, active ports, predicted peak gain / direction / pointing error |
| `build_board.py` | AEDT (IronPython) script: substrate, ground, 7 pins + lumped ports, {n} slotted circles, radiation box, setup 10.54 GHz, far-field sphere |
| `set_weights.py` | AEDT script: sets the port sources (EditSources) from one table row (edit THETA_DEG / PHI_DEG at the top); prints the values if the call fails |
| `layout.csv`, `layout.json`, `layout.png` | every cell (centre, slot, tile/port), tile hexagons, feeds |

## Run order
1. AEDT: Tools > Run Script > `build_board.py` (a few minutes). Check the messages; it saves `phase29_compact7.aedt`. `ADD_SWEEP = True` at the top adds a 10.0-11.1 GHz sweep (slower; not needed for the weights).
2. Validate, then Analyze (Setup1, 10.54 GHz).
3. **Read the S-matrix first**: S11..S77 and the coupling between neighbouring ports (centre-outer, outer-outer). If a neighbouring |S_ij| is above about -15 dB the uncoupled weights are unreliable (note it in the comparison).
4. Port alone: Results > Edit Sources: port k = 1 W, all others 0 W, then the far field (Radiation > Far Fields, Realized Gain, RHCP/LHCP, Axial Ratio). Compare with `pattern_cuts.csv` (`portK_alone`).
5. Weighted cases: run `set_weights.py` after setting THETA_DEG / PHI_DEG (or type the row of `weights_table.csv` by hand), re-run the far field; compare with the weighted examples below.

## How to apply the port weights in HFSS (Edit Sources)
Results (project tree) > Edit Sources... > type **Power** (not Voltage), then for each Port1..Port7: **Magnitude = `P<k>_power_W` from the table in W**, **Phase = `P<k>_phase_deg` in deg**; ports with 0 W stay terminated in 50 ohm. The table's weights are field weights
w_k with |w_k|^2 = incident power (sum = 1 W) and phase = angle(w_k), exp(+j w t) (the project convention, also HFSS's); the phase reference of all seven ports is the identical port plane, so a common offset cancels.
If the beams come out mirrored in phi, the phase sign convention differs: negate all phases. The model assumes uncoupled ports: with strong coupling the realised excitation differs (read the HFSS "port incident / accepted" values). `set_weights.py` automates this
(EditSources call NOT verified against an AEDT version; wrapped, it prints the values to type if it fails).

## Predictions (MODEL_ONLY, 10.54 GHz)
Port alone (1 W accepted; RHCP realised gain; SLL = highest RHCP outside 1.5 x HPBW about the peak, relative):

| Port | Feed (mm) | Design beam (theta, phi) | Peak gain (dBi) | Peak (theta, phi) | HPBW (deg) | SLL (dB) | AR peak (dB) |
|---|---|---|---|---|---|---|---|
{tab}

Weighted examples (sum of port powers = 1 W):

| Case | Active ports | Peak gain (dBi) | Peak (theta, phi) | Pointing error (deg) | HPBW (deg) | AR peak (dB) |
|---|---|---|---|---|---|---|
{etab}

Scan over the table (tuned weights, <= 3 active ports): theta <= 30 deg: worst pointing error {inn["max_pointing_error_deg"]:.1f} deg (median {inn["median_pointing_error_deg"]:.1f}), gain dip {inn["max_gain_dip_db"]:.1f} dB, minimum peak gain {inn["min_gain_dbi"]:.1f} dBi
(criterion: pointing <= 8 deg, dip <= 3 dB: {"MET" if p["pass_relaxed_inner"] else "NOT MET"}). 30-45 deg band, as is: worst pointing error {band["max_pointing_error_deg"]:.1f} deg, only {band["fraction_within_8deg"] * 100:.0f} % of directions within 8 deg,
gain dip {band["max_gain_dip_db"]:.1f} dB, minimum peak gain {band["min_gain_dbi"]:.1f} dBi. Gains are LOW by design (small tiles, weak slot leakage: about 6-8 dBi). The tilted tiles peak near theta 26 deg (not 30) in the model, so the theta 30 deg rows of the table are reached by the combination of ports, not by one tile.

What would falsify the model: a port alone whose peak is more than ~10 deg off its design direction or more than 5 dB from the table; a weighted example whose peak is more than ~15 deg off the target; neighbouring-port coupling above -10 dB.
The weights are a starting point: with coupling in HFSS they should be re-tuned from the measured port patterns (the 7 single-port patterns are all that is needed, the model is linear).
"""


if __name__ == "__main__":
    main()
