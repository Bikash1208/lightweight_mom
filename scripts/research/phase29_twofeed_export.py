"""Phase 29: two-feed HFSS check package `hfss_exports/phase29_twofeed/` (MODEL_ONLY): shared 150 mm hexagon, 2 monopoles 1.35 lambda0 apart on the x axis, two pencil beams one
beamwidth apart (Frey's own N = 2 geometry: feed 1 at +x tilting +x, feed 2 at -x tilting -x), plus a single-feed CONTROL board (same hexagon, centre feed, broadside).
Purpose: measure in HFSS the SLL / gain loss of the shared-aperture superposition (model: SLL -22 dB single feed -> about -11 dB for two feeds).
Predictions of the OLD model (radial reflection cell) are written now; `--new-model` (after the angle-dependent model exists) adds `predictions_new_model.json`.

    python scripts/research/phase29_twofeed_export.py [--new-model]"""

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts/research"))
OUT = ROOT / "hfss_exports/phase29_twofeed"
LAM = 299792458.0 / 10.54e9 * 1e3
R_MM, SPACING, HPBW_SINGLE = 150.0, 1.35 * LAM, 7.3


def cells_of(d):
    return [
        {
            "x_mm": float(x * 1e3),
            "y_mm": float(y * 1e3),
            "slot_L_mm": float(a),
            "slot_W_mm": float(b),
            "slot_angle_deg": float(np.degrees(c)),
        }
        for x, y, a, b, c in zip(d["x_m"], d["y_m"], d["slot_L_mm"], d["slot_W_mm"], d["axis_rad"], strict=True)
    ]


def old_model_designs(table):
    from src.hologram import shared_aperture_frey as fr
    from src.hologram import shared_aperture_multiport as sa
    from src.hologram import tensor_cp_design as t
    from src.hologram.tensor_leaky_cp import march_rays

    mf = lambda j, f, b: march_rays(table, sa.hexagon_m(R_MM), f, t.K0, b[0], b[1], 0.9, taper_c=0.5, n_psi=360)  # noqa: E731
    tilt = HPBW_SINGLE / 2
    two = [(SPACING / 2, 0.0, tilt, 0.0), (-SPACING / 2, 0.0, tilt, 180.0)]
    d2 = fr.design_frey(table, R_MM, two, [0.0, 0.0], 1000.0, 1.0, n_iter=2, n_psi=360, march_fn=mf)
    d1 = fr.design_frey(table, R_MM, [(0.0, 0.0, 0.0, 0.0)], [0.0], 1000.0, 1.0, n_iter=2, n_psi=360, march_fn=mf)
    return d2, d1


def predictions(d, F, th, ph, beams, weights_cases):
    import phase29_compact7_export as ce
    from phase29_pointB_sweep import _hpbw

    from src.hologram import seven_hex_cluster as c
    from src.hologram import shared_aperture_multiport as sa

    out = {"per_port_alone": [], "weighted": {}}
    cuts = []
    for j, f in enumerate(F, start=1):
        m = ce.metrics(*f, th, ph, beams[j - 1], _hpbw, c._unit)
        m.update(
            {
                "port": j,
                "feed_xy_mm": [d["feeds_m"][j - 1][0] * 1e3, d["feeds_m"][j - 1][1] * 1e3],
                "beam_design_theta_phi_deg": list(beams[j - 1]),
            }
        )
        out["per_port_alone"].append(m)
        cuts += [(f"port{j}_alone", beams[j - 1][1], *r) for r in ce.cut(*f, th, ph, beams[j - 1][1])]
    for name, w in weights_cases.items():
        w = np.asarray(w, dtype=complex)
        w = w / np.linalg.norm(w)
        et, ep = sa.combine(F, w)
        m = ce.metrics(et, ep, th, ph, (0.0, 0.0), _hpbw, c._unit)
        m.update(
            {
                "port_power_W": [float(abs(x) ** 2) for x in w],
                "port_phase_deg": [float(np.degrees(np.angle(x))) for x in w],
            }
        )
        out["weighted"][name] = m
        cuts += [(name, 0.0, *r) for r in ce.cut(et, ep, th, ph, 0.0)]
    if len(F) == 2:
        pk = [tuple(p["peak_theta_phi_deg"]) for p in out["per_port_alone"]]
        v = c._unit(*pk[0]) + c._unit(*pk[1])
        v = v / np.linalg.norm(v)
        mid = (float(np.degrees(np.arccos(np.clip(v[2], -1, 1)))), float(np.degrees(np.arctan2(v[1], v[0])) % 360.0))
        G = [sa.gain_db(*f) for f in F]
        i, k = int(round(mid[0])), int(round(mid[1] / (ph[1] - ph[0]))) % ph.size
        out["crossover_db_each"] = [float(G[0][i, k] - G[0].max()), float(G[1][i, k] - G[1].max())]
    return out, cuts


def write_set(d, F, th, ph, beams, cases, name, model_tag, folder):
    import phase29_compact7_export as ce

    pred, cuts = predictions(d, F, th, ph, beams, cases)
    with open(folder / f"pattern_cuts_{model_tag}.csv", "w", newline="", encoding="utf-8") as f:
        w_ = csv.writer(f)
        w_.writerow(["case", "cut_phi_deg", "theta_signed_deg", "rhcp_gain_dbi", "lhcp_gain_dbi"])
        for cs, pc, tsg, gr, gl in cuts:
            w_.writerow([cs, f"{pc:g}", tsg, f"{gr:.2f}", f"{gl:.2f}"])
    with open(folder / f"patterns_3deg_{model_tag}.csv", "w", newline="", encoding="utf-8") as f:
        w_ = csv.writer(f)
        w_.writerow(["port", "theta_deg", "phi_deg", "rhcp_gain_dbi", "lhcp_gain_dbi"])
        from src.hologram import shared_aperture_multiport as sa

        for j, fj in enumerate(F, start=1):
            gr, gl = sa.gain_db(*fj), ce.lhcp_gain_db(*fj)
            for i in range(0, 91, 3):
                for k in range(0, 360, 6):
                    w_.writerow([j, i, k, f"{gr[i, k]:.2f}", f"{gl[i, k]:.2f}"])
    return pred


def main():
    import phase29_B3_debug as dbg

    from src.export import shared7_export as ex
    from src.export.compact7_export import _COMMON
    from src.hologram import shared_aperture_multiport as sa

    ap = argparse.ArgumentParser()
    ap.add_argument("--new-model", action="store_true")
    a = ap.parse_args()
    dbg._init()
    th, ph = np.arange(0.0, 90.0 + 1e-9, 1.0), np.arange(0.0, 360.0, 1.0)
    OUT.mkdir(parents=True, exist_ok=True)
    if not a.new_model:
        d2, d1 = old_model_designs(dbg._S["table"])
        F2, F1 = sa.port_fields(d2, dbg._S["green"], th, ph), sa.port_fields(d1, dbg._S["green"], th, ph)
        tilt = HPBW_SINGLE / 2
        two_cases = {"port1_plus_port2_in_phase": [1, 1], "port1_minus_port2_antiphase": [1, -1]}
        p2 = write_set(d2, F2, th, ph, [(tilt, 0.0), (tilt, 180.0)], two_cases, "two", "old_model", OUT)
        p1 = write_set(d1, F1, th, ph, [(0.0, 0.0)], {}, "control", "old_model_control", OUT)
        (OUT / "build_board.py").write_text(
            ex.aedt_build_script(
                cells_of(d2), [(f[0] * 1e3, f[1] * 1e3) for f in d2["feeds_m"]], name="phase29_twofeed"
            ).replace("7 monopoles, 7 lumped ports", "2 monopoles, 2 lumped ports"),
            encoding="utf-8",
        )
        (OUT / "build_board_control_1feed.py").write_text(
            ex.aedt_build_script(cells_of(d1), [(0.0, 0.0)], name="phase29_twofeed_control").replace(
                "7 monopoles, 7 lumped ports", "1 monopole, 1 lumped port"
            ),
            encoding="utf-8",
        )
        ws = _COMMON.replace("@TITLE@", "port sources for the two-feed checks") + (
            "# Applies port sources (incident power W, phase deg) to the active HFSS design. MODEL_ONLY. EditSources is NOT verified against an AEDT version: wrapped; prints the values to type if it fails.\n"
            "CASE = 0   # 0: port 1 alone, 1: port 2 alone, 2: in phase (sum), 3: antiphase (difference)\n"
            "CASES = [[(1.0, 0.0), (0.0, 0.0)], [(0.0, 0.0), (1.0, 0.0)], [(0.5, 0.0), (0.5, 0.0)], [(0.5, 0.0), (0.5, 180.0)]]\n"
            "oProject = oDesktop.GetActiveProject()\noDesign = oProject.GetActiveDesign()\n"
            "def note(level, text):\n    try:\n        oDesktop.AddMessage(oProject.GetName(), oDesign.GetName(), level, text)\n    except Exception:\n        pass\n"
            '    try:\n        sys.stdout.write(str(text) + "\\n")\n    except Exception:\n        pass\n'
            'row = CASES[CASE]\nports = ["Port1", "Port2"]\n'
            'note(0, "Case %d: " % CASE + "; ".join("%s %.4g W %.1f deg" % (ports[k], row[k][0], row[k][1]) for k in range(2)))\n'
            'try:\n    sol = oDesign.GetModule("Solutions")\n'
            '    sol.EditSources("PowerWatts", ["NAME:Names"] + ports, ["NAME:Magnitudes"] + ["%.6gW" % row[k][0] for k in range(2)],\n'
            '                    ["NAME:Phases"] + ["%.3fdeg" % row[k][1] for k in range(2)], ["NAME:Terminated"], ["NAME:Resistances"], ["NAME:Reactances"])\n'
            '    note(0, "Sources set. Check Results > Edit Sources before analysing.")\n'
            'except Exception as exc:\n    note(1, "EditSources failed (%s): set the sources by hand in Results > Edit Sources with the values above." % exc)\n'
        )
        (OUT / "set_weights.py").write_text(ws, encoding="utf-8")
        pred = {
            "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "SOLVER_SUPERCELL", "INTER_PORT_COUPLING_NOT_MODELLED"],
            "model": "OLD: radial reflection cell (kappa independent of the incidence angle)",
            "frequency_ghz": 10.54,
            "geometry": {
                "hexagon_circumradius_mm": R_MM,
                "feed_spacing_mm": SPACING,
                "feed_spacing_lambda0": 1.35,
                "beam_tilt_deg": HPBW_SINGLE / 2,
                "beam_separation_deg": HPBW_SINGLE,
                "n_cells": len(cells_of(d2)),
                "control_n_cells": len(cells_of(d1)),
            },
            "two_feed": p2,
            "single_feed_control": p1,
            "reading": "model: single feed SLL about -22 dB, two feeds about -11 dB (Frey et al. report -18.5 dB for their two-feed case): HFSS decides whether the superposition loses ~11 dB of SLL.",
            "s11": "NOT_AVAILABLE_FROM_CURRENT_DATA",
            "port_coupling": "NOT_AVAILABLE_FROM_CURRENT_DATA (not modelled)",
            "new_model": "predictions_new_model.json is added when the angle-dependent model exists",
        }
        (OUT / "predictions.json").write_text(json.dumps(pred, indent=1, default=float), encoding="utf-8")
        (OUT / "README.md").write_text(readme(pred), encoding="utf-8")
        print("wrote", OUT)
    else:
        from src.hologram import shared_aperture_tensor as st

        table = dbg._S["table"]
        lut = st.AngleLUT(ROOT / "results/phase29/angle_lut_v1.npz", table)
        tilt = HPBW_SINGLE / 2
        two = [(SPACING / 2, 0.0, tilt, 0.0), (-SPACING / 2, 0.0, tilt, 180.0)]
        d2 = st.design_tensor(table, lut, R_MM, two, [0.0, 0.0], n_psi=360)
        d1 = st.design_tensor(table, lut, R_MM, [(0.0, 0.0, 0.0, 0.0)], [0.0], n_psi=360)
        F2, F1 = sa.port_fields(d2, dbg._S["green"], th, ph), sa.port_fields(d1, dbg._S["green"], th, ph)
        p2 = write_set(
            d2,
            F2,
            th,
            ph,
            [(tilt, 0.0), (tilt, 180.0)],
            {"port1_plus_port2_in_phase": [1, 1], "port1_minus_port2_antiphase": [1, -1]},
            "two",
            "new_model",
            OUT,
        )
        p1 = write_set(d1, F1, th, ph, [(0.0, 0.0)], {}, "control", "new_model_control", OUT)
        old = json.loads((OUT / "predictions.json").read_text())
        new = {
            "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "SOLVER_SUPERCELL", "INTER_PORT_COUPLING_NOT_MODELLED"],
            "model": "NEW: angle-dependent full-tensor cell (results/phase29/angle_lut_v1.npz, src/hologram/shared_aperture_tensor.py)",
            "two_feed": p2,
            "single_feed_control": p1,
            "old_model_for_comparison": {
                "two_feed_per_port": [
                    {k: q[k] for k in ("peak_rhcp_gain_dbi", "hpbw_deg", "sll_rel_db_outside_1p5_hpbw")}
                    for q in old["two_feed"]["per_port_alone"]
                ],
                "control": {
                    k: old["single_feed_control"]["per_port_alone"][0][k]
                    for k in ("peak_rhcp_gain_dbi", "hpbw_deg", "sll_rel_db_outside_1p5_hpbw")
                },
            },
        }
        (OUT / "predictions_new_model.json").write_text(json.dumps(new, indent=1, default=float), encoding="utf-8")
        print("wrote new-model predictions")


def readme(p):
    t, c = p["two_feed"], p["single_feed_control"]
    a, b = t["per_port_alone"]
    return f"""# phase29_twofeed: does the shared-aperture superposition lose SLL? (HFSS check, Phase 29, MODEL_ONLY)

Shared 150 mm-circumradius hexagon (vertex at +/-y), RT/duroid 5880 1.575 mm, 10.54 GHz, {p["geometry"]["n_cells"]} slotted-circle cells (D 6 mm, 6.4 mm lattice), TWO monopoles on the x axis {p["geometry"]["feed_spacing_mm"]:.1f} mm = 1.35 lambda0 apart
(port 1 at +x, beam tilted {p["geometry"]["beam_tilt_deg"]:.2f} deg towards +x; port 2 at -x, beam tilted towards -x): two pencil beams one beamwidth ({p["geometry"]["beam_separation_deg"]:.1f} deg) apart, the geometry of Frey et al.
(IEEE OJAP 2023, 10.1109/OJAP.2023.3261903) Section III. **Control board** `build_board_control_1feed.py`: same hexagon, ONE centre feed, broadside: the single-feed reference of the same cell family.
All numbers are MODEL_ONLY predictions of the OLD model (radial reflection cell); `predictions_new_model.json` (angle-dependent tensor model) is added later. S11 and port coupling: NOT_AVAILABLE_FROM_CURRENT_DATA.

## Files
`build_board.py` (2 ports), `build_board_control_1feed.py` (1 port), `set_weights.py` (CASE 0..3: port 1 alone, port 2 alone, sum, difference), `predictions.json`, `pattern_cuts_old_model*.csv`, `patterns_3deg_old_model*.csv` (RHCP/LHCP gain, hemisphere 3 x 6 deg), `layout`: in the scripts (CELLS).

## Run order
1. AEDT: Tools > Run Script > `build_board.py` (and later the control script in a new project). Messages are printed; projects are saved as `phase29_twofeed.aedt` / `phase29_twofeed_control.aedt`.
2. Validate, Analyze Setup1 (10.54 GHz single frequency; ADD_SWEEP = True at the top adds 10.0-11.1 GHz). One solve gives both excitations: do NOT re-solve per port.
3. Read S11, S22, S21 first (the coupling at 1.35 lambda0).
4. Far field (infinite sphere, theta 0-180, phi -180..180, 1 deg): port 1 alone, port 2 alone (Results > Edit Sources: 1 W on one port, 0 W on the other; or `set_weights.py` CASE), then sum and difference. Export Realized Gain RHCP/LHCP and Axial Ratio.
5. Control board: the same far-field export for its single port.
Compare the SLL of the two-feed beams with the control's: the loss is the point of this check.

## Predictions (old model; RHCP realised gain per 1 W accepted; SLL = highest RHCP lobe outside 1.5 HPBW of the peak, relative)
| Case | Peak gain (dBi) | Peak (theta, phi) | HPBW (deg) | SLL (dB) | AR (dB) |
|---|---|---|---|---|---|
| port 1 alone | {a["peak_rhcp_gain_dbi"]:.1f} | ({a["peak_theta_phi_deg"][0]:.0f}, {a["peak_theta_phi_deg"][1]:.0f}) | {a["hpbw_deg"]:.1f} | {a["sll_rel_db_outside_1p5_hpbw"]:.1f} | {a["axial_ratio_at_peak_db"]:.1f} |
| port 2 alone | {b["peak_rhcp_gain_dbi"]:.1f} | ({b["peak_theta_phi_deg"][0]:.0f}, {b["peak_theta_phi_deg"][1]:.0f}) | {b["hpbw_deg"]:.1f} | {b["sll_rel_db_outside_1p5_hpbw"]:.1f} | {b["axial_ratio_at_peak_db"]:.1f} |
| single-feed control | {c["per_port_alone"][0]["peak_rhcp_gain_dbi"]:.1f} | ({c["per_port_alone"][0]["peak_theta_phi_deg"][0]:.0f}, {c["per_port_alone"][0]["peak_theta_phi_deg"][1]:.0f}) | {c["per_port_alone"][0]["hpbw_deg"]:.1f} | {c["per_port_alone"][0]["sll_rel_db_outside_1p5_hpbw"]:.1f} | {c["per_port_alone"][0]["axial_ratio_at_peak_db"]:.1f} |

Crossover of the two beams (each at the midpoint of the pair, re its own peak): {t["crossover_db_each"][0]:.1f} / {t["crossover_db_each"][1]:.1f} dB. Sum (in phase): {t["weighted"]["port1_plus_port2_in_phase"]["peak_rhcp_gain_dbi"]:.1f} dBi;
difference: {t["weighted"]["port1_minus_port2_antiphase"]["peak_rhcp_gain_dbi"]:.1f} dBi (see `predictions.json`). What would falsify the model: a two-feed SLL near the control's (the loss is absent in HFSS), or a gain drop far from the {c["per_port_alone"][0]["peak_rhcp_gain_dbi"] - a["peak_rhcp_gain_dbi"]:.1f} dB predicted between control and two-feed ports.
"""


if __name__ == "__main__":
    main()
