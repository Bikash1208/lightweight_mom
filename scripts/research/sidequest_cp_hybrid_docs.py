"""Predictions JSON + HFSS_SETUP.md + sync for the CP+LP hybrid package.  MODEL_ONLY."""
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sidequest_cp_hybrid as H

OUT = ROOT / "hfss_exports/sidequest_board_CP_LP_hybrid"
SYNC = Path(r"C:\hmsa_sync\data\sidequest\board_CP_LP_hybrid")
RES = ROOT / "results/sidequest_hybrid"
dz = np.load(RES / "cp_hybrid_design_n180_lb2.0.npz")
cache = H.tensor_design(*H.load_cells()[:2])
out = {}
for rs in (-1.0, +1.0):
    rr, _ = H.run(rs, gb_mask_override=dz["gb"], flip=dz["flip"], cache=cache, verbose=False)
    out[f"rel{int(rs):+d}"] = H.summarize_run(rr)
    out[f"rel{int(rs):+d}"]["detail"] = {k: {a: b for a, b in v.items() if a != "_maps"} for k, v in rr.items() if k in ("GB_OFF", "GB_ON")}
sweep = json.loads((RES / "cp_hybrid_sweep.json").read_text())
P = out["rel-1"]
d_off, d_on = P["detail"]["GB_OFF"], P["detail"]["GB_ON"]
pred = {
    "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "CP part: Phase 26 tensor engine anchored to the HFSS 23.03 dBi (anchor -0.48 dB)", "GB part: SWHB-consistent Born element, A/B board calibrated (launch 0.62 incl.)"],
    "design": {"beam1": "fixed RHCP at (25, 0), tensor rotating-slot cells (Phase 26 table, eta_design 0.9, taper 0.5), 1250 cells of which 50 are rotated by 90 deg (joint optimisation)",
               "beam2": "LP-V (phi-hat) at (25, 180), 180 GB crossed-slit cells (2 diodes each = 360 diodes), ONE group GB, S_w 1e6 (OFF, default) / 6.1 ohm/sq (ON)",
               "joint_optimisation": "GB sites (102 of 180 relocated from the A/B board) + 50 tensor slots rotated by 90 deg so that the tensor background at beam 2 cancels the GB OFF field and adds to the GB ON field; beam-1 amplitude loss budget 2 dB (realised: see below)"},
    "primary_relative_sign": "-1 (physical argument: a slit across E RAISES the sheet susceptance (K2 library: open cross slit 8.69 mS vs solid disc 5.18 mS) -> tensor J = -c M E in the Phase 26 convention)",
    "GB_OFF_default": {"beam1_rhcp_gain_dbi": d_off["beam1"]["rhcp_gain_dbi"], "beam1_peak_dir": d_off["beam1"]["peak_dir"], "AR_peak_db": d_off["beam1"]["ar_peak_db"], "AR_worst_in_minus3dB_beam_db": d_off["beam1"]["ar_worst_in_3dB_beam_db"], "XPD_db": d_off["beam1"]["xpd_db"],
                       "HPBW_deg": d_off["beam1"]["hpbw_deg"], "level_at_25_180_dbi_total": d_off["beam2"]["total_gain_dbi_at_design_dir"], "SLL_db_with_2p5_offset": d_off["sll_total_db_with_offset"], "SLL_db_model": d_off["sll_total_db_model"]},
    "GB_ON": {"beam1_rhcp_gain_dbi": d_on["beam1"]["rhcp_gain_dbi"], "AR_peak_db": d_on["beam1"]["ar_peak_db"], "AR_worst_in_minus3dB_beam_db": d_on["beam1"]["ar_worst_in_3dB_beam_db"],
              "beam2_total_gain_dbi_at_25_180": d_on["beam2"]["total_gain_dbi_at_design_dir"], "beam2_rhcp_dbi": d_on["beam2"]["rhcp_dbi"], "beam2_lhcp_dbi": d_on["beam2"]["lhcp_dbi"], "beam2_pol_tilt_from_thetahat_deg": d_on["beam2"]["pol_tilt_from_theta_hat_deg"],
              "beam2_XPD_phihat_vs_thetahat_db": d_on["beam2"]["xpd_phi_hat_vs_theta_hat_db"], "beam2_HPBW_deg": d_on["beam2"]["hpbw_deg"], "SLL_db_with_2p5_offset": d_on["sll_total_db_with_offset"], "SLL_db_model": d_on["sll_total_db_model"]},
    "beam2_contrast_db": P["contrast_db"], "two_beam_balance_db_beam1_minus_beam2_ON": P["balance_b1_minus_b2_on_db"],
    "comparison_current_LP_board_HFSS": {"beam1_dbi": 13.2, "contrast_db": 6.5, "GB_ON_beam2_dbi": 8.8, "SLL_db": -5},
    "tensor_beta_over_k0": {"min_mean_max": [1.051, 1.062, 1.073], "GB_static_model_beta": "1.056-1.063 (A/B calibration)", "wave_speed_jump": "none (<= 0.01)"},
    "bands": {"beam1_gain_db": "+-1.0 (Phase 26 anchored; model AR is pessimistic: HFSS 1.2-1.3 vs model 2.5 on that board)", "beam2_gain_db": "+-3", "contrast_db": "+-3 and SIGN-CRITICAL (see falsifiers)", "pointing_deg": "+-2", "SLL_db": "+-3 (the +2.5 dB offset is applied)"},
    "falsifiers": ["beam 2 (25, 180) is HIGHER with S_w = 1e6 than with 6.1 ohm/sq (relative sign of tensor and GB cells wrong: the design inverts; the model gives beam 2 at +7.9 dBi OFF and -1.3 dBi ON)",
                   "beam-1 RHCP below 20 dBi or AR at the peak above 3 dB in either state", "beam-2 contrast below 5 dB", "peak direction of beam 1 more than 3 deg from (25, 0)"],
    "relative_sign_alternative_model_result": {k: out["rel+1"][k] for k in ("b2_off_dbi", "b2_on_dbi", "contrast_db", "b1_off_dbi", "ar_off")},
    "diode_reduction_model_only_rel_minus1": {k: {a: v[a] for a in ("ncell", "diodes", "b2_on_dbi", "b2_off_dbi", "contrast_db", "balance_b1_minus_b2_on_db", "b1_off_dbi", "ar_off", "sll_off_with_offset", "flips", "gb_sites_changed")} for k, v in sweep.items()},
}
(OUT / "predictions.json").write_text(json.dumps(pred, indent=1, default=float), encoding="utf-8", newline="\n")
(OUT / "design_summary.json").write_text(json.dumps({"summary": {k: v for k, v in P.items() if k != "detail"}, "model": "sidequest_cp_hybrid.py (rel -1)"}, indent=1, default=float), encoding="utf-8", newline="\n")
g = json.loads((OUT / "geometry_checks.json").read_text())
setup = f"""# HFSS_SETUP - CP + LP hybrid board: fixed RHCP beam 1 (25, 0) from tensor cells; LP-V beam 2 (25, 180) switched by ONE group GB; 150 mm hexagon, RT5880 1.575, 10.5 GHz. MODEL_ONLY. One page.
**1. Assign material `pec` to `METAL_PEC` (patches, ground and probe: ONE named part, {g['checks']['metal_single_part_lumps']['lumps']} DC-isolated lumps; the STEP imports as vacuum).** `SUBSTRATE`: RT5880 eps_r 2.2, tan d 0.0009. `AIR`: vacuum, first in the tree (overlaps the solids on purpose).
**Files** (mm): `board_CP_LP_hybrid.step` ({g['step_mb']} MB: `METAL_PEC`, `SUBSTRATE`, `AIR`, `PORT`, 360 sheets `GB_<ix>_<iy>_DX|DY`), `states.csv`, `sheets.csv`, `cells.csv` (per cell: kind, slot L x W and angle), `predictions.json`, `design_summary.json`, `geometry_checks.json` (all PASS on the re-imported STEP, EXACT_DIMS PASS).
**Geometry**: same hexagon, stack, probe and port as `sidequest_board_k2_lp` (substrate hexagon (0, +-154), (+-133.368, +-77), 1.575 thick, ground 35 um with an opening r 2.05, probe r 0.635 z -1.575..+6.0 above the substrate, `PORT` annulus 0.635..2.05 at z -1.575), 1430 cells on the 6.4 lattice (no cell for rho < 6.5): **1250 tensor cells** = disc r 3.0 (35 um) with ONE centred slot L x W at angle (Phase 26 table slots 4.4 x 1.2, 4.6 x 1.6, 4.8 x 1.6, 4.8 x 2.0; slot corners rounded to 1 um, L and W within +-0.0015 mm; 50 of them rotated by 90 deg by the joint optimisation) and **180 GB cells** = disc r 3.0 with a cross slit 0.4 x 5.6 (rim 0.2) and two sheets `GB_*_DX` (0, +0.8) / `GB_*_DY` (+0.8, 0), 0.4 x 1.0 at z 0.0175. **No sheets and no diodes anywhere else (the former GA region and the static cells are tensor cells).**
**2. Boundaries**: Radiation on all six faces of `AIR` (z -1.61..+20, 14 mm above the probe top, > lambda/4; x/y pad 10 mm; one face per side). PML only if licence and RAM allow; say which was used.
**3. Port**: `PORT` -> Lumped Port 50 ohm, integration line (0.635, 0, -1.575) -> (2.05, 0, -1.575), renormalise 50, no de-embedding.
**4. The two states (same mesh, one project, variable `S_w`)**: ONE Impedance boundary on ALL `GB_*` sheets, reactance 0. **GB OFF (default): S_w = 1e6 ohm/sq** (4e5 ohm across the 0.4 mm gap = open). **GB ON: S_w = 6.1 ohm/sq** (R_on 2.44 ohm x W/L = 1.0/0.4 = 2.5). (The 28.16 fF OFF capacitance is not included, as for the validated v2/addbeam packages; say if you add it.)
**5. Solution**: ONE frequency **10.5 GHz**, Driven Modal, first-order basis, max passes 6, max delta S 0.02, no sweep. Initial mesh: substrate max length 1.0 mm, patch faces 0.4 mm, `GB_*` sheets 0.15 mm, probe and port 0.3 mm, lambda refinement on. Mesh reuse OFF -> ON allowed (say so; S11 must agree within 0.5 dB with a fresh ON).
**6. Export (per state)**: far-field infinite sphere theta 0..180, phi 0..360 step **2 deg**, **complex rETheta and rEPhi in V** (full sphere), `rE_GB_OFF.csv`, `rE_GB_ON.csv`; complex S11 at 10.5 GHz; antenna parameters (realized gain RHCP and LHCP and total, directivity, radiation efficiency, accepted power); passes, delta S, tets, peak RAM, solve time per state. **Also report the axial ratio at the beam-1 peak and its 3 dB-AR beamwidth, and the polarisation tilt at (25, 180).**
**Model predictions (MODEL_ONLY; bands in `predictions.json`)**: GB OFF (default): beam 1 RHCP **{d_off['beam1']['rhcp_gain_dbi']:.1f} dBi** at (25, 0), AR {d_off['beam1']['ar_peak_db']:.1f} dB (worst in the -3 dB beam {d_off['beam1']['ar_worst_in_3dB_beam_db']:.1f}; HFSS AR of the Phase 26 board was 1.2-1.3 against the model 2.5), XPD {d_off['beam1']['xpd_db']:.0f} dB, HPBW {d_off['beam1']['hpbw_deg']:.1f} deg, level at (25, 180) {d_off['beam2']['total_gain_dbi_at_design_dir']:.1f} dBi, SLL {d_off['sll_total_db_with_offset']:.1f} dB (model {d_off['sll_total_db_model']:.1f} + 2.5 offset). GB ON: beam 1 {d_on['beam1']['rhcp_gain_dbi']:.1f} dBi AR {d_on['beam1']['ar_peak_db']:.1f} dB; **beam 2 LP-V {d_on['beam2']['total_gain_dbi_at_design_dir']:.1f} dBi** at (25, 180) (RHCP {d_on['beam2']['rhcp_dbi']:.1f} / LHCP {d_on['beam2']['lhcp_dbi']:.1f} dBi components, tilt {d_on['beam2']['pol_tilt_from_theta_hat_deg']:.0f} deg from theta-hat i.e. close to phi-hat, phi-hat/theta-hat {d_on['beam2']['xpd_phi_hat_vs_theta_hat_db']:.0f} dB), HPBW {d_on['beam2']['hpbw_deg']:.1f} deg, SLL {d_on['sll_total_db_with_offset']:.1f} dB. **Beam-2 contrast {P['contrast_db']:.1f} dB, balance (beam 1 - beam 2) {P['balance_b1_minus_b2_on_db']:.1f} dB** (current LP board HFSS: 13.2 dBi, contrast 6.5 dB).
**KEY RISK / first check**: the relative sign of the tensor cells and the GB cells is taken from a physical argument (a slit across E raises the sheet susceptance). If beam 2 at (25, 180) is HIGHER with S_w = 1e6 than with 6.1 ohm/sq the sign is opposite: tell me, the design is then regenerated (the model gives beam 2 = +7.9 dBi OFF and -1.3 dBi ON in that case).
"""
(OUT / "HFSS_SETUP.md").write_text(setup, encoding="utf-8", newline="\n")
SYNC.mkdir(parents=True, exist_ok=True)
for f in OUT.iterdir():
    shutil.copy2(f, SYNC / f.name)
print(sorted(x.name for x in SYNC.iterdir()))
print(json.dumps({k: pred[k] for k in ("GB_OFF_default", "GB_ON", "beam2_contrast_db", "two_beam_balance_db_beam1_minus_beam2_ON")}, default=float)[:1500])
