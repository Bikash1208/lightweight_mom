"""Finalise the add-a-beam (b) package docs (HFSS_SETUP.md, predictions.json, states.csv) and sync it.  MODEL_ONLY."""
import csv
import json
import re
import shutil
from pathlib import Path

d = Path("hfss_exports/sidequest_board_RO3006_r2p30_addbeam_b")
des = json.loads(Path("results/sidequest_hybrid/ro3006_board_design_r2.3_addbeam_b.json").read_text())
S = des["addbeam_states"]
OFF = {1: -1.78, 2: -2.7}


def g(s, k):
    return s["beam_gain_dbi"][k]


def raw(s, k):
    return s["beam_gain_raw_model_dbi"][k]


p = json.loads(d.joinpath("predictions.json").read_text())  # may contain NaN
states = []
for k, s in enumerate(S):
    hp = s["hpbw"]
    states.append({
        "state": str(k + 1), "name": "default (unpowered, diodes OFF): beam 1 only" if k == 0 else "circuit ON (diodes ON): beam 1 + beam 2",
        "D_sheets": "no boundary (or Impedance Rs = 1e6 ohm/sq)" if k == 0 else "one Impedance boundary, Rs = 6.1 ohm/sq (W/L = 2.5), reactance 0",
        "wanted_beams_theta_phi_deg": [[25.0, 0.0]] if k == 0 else [[25.0, 0.0], [25.0, 180.0]],
        "beam1_25_0_gain_dbi_calibrated": round(g(s, "b1_(25,0)"), 2), "beam2_25_180_gain_dbi_calibrated": round(g(s, "b2_(25,180)"), 2),
        "beam1_gain_dbi_raw_model": round(raw(s, "b1_(25,0)"), 2), "beam2_gain_dbi_raw_model": round(raw(s, "b2_(25,180)"), 2),
        "peak_dir_theta_phi_deg": s["peak_dir"], "sll_first_null_db_rel_stronger_beam": round(s["SLL_firstnull_db"], 2), "sll_first_null_db_rel_weaker_beam": round(s["SLL_firstnull_vs_weaker_db"], 2), "sll_dir_theta_phi": s["SLL_firstnull_dir"],
        "beam_balance_b1_minus_b2_db": None if k == 0 else round(s["balance_db"], 2), "isolation_beam1_minus_level_at_beam2_db": round(s["isolation_db"], 2) if k == 0 else None,
        "expected_hpbw_deg": {("beam1" if int(b) == 0 else "beam2"): {"elevation": round(v["elevation"], 1), "azimuth_angular": round(v["azimuth_angular"], 1)} for b, v in hp.items()},
        "ring_median_rel_peak_db_model": s["ring_median_rel_peak_db"], "nonbeam_floor_db": round(s["floor_db"], 2)})
p["states"] = states
p["falsifiers"] = "beam 1 more than 3 deg from (25, 0) in either state; beam-2 on/off contrast (state 2 minus state 1 level at (25, 180)) below 8 dB; beam balance in state 2 worse than 6 dB; beam 1 more than 4 dB below the prediction after re-calibration; strongest sidelobe above the gap-corrected SLL by more than 3 dB"
p["labels"] = ["MODEL_ONLY", "add-a-beam (b): beam 1 always on, beam 2 switched; w2 1.25, eps_B 0; each state designed with its own real beta/k0 (1.33 / 1.26)",
               "gains calibrated with the board_K2_default offsets (state 1 -1.78 dB, state 2 -2.7 dB): TO BE RE-CALIBRATED with the v2 HFSS offsets before running (see recalibration)"]
p["design"]["id"] = "RO3006 r 2.30, add-a-beam (b), one group D_*"
p["design"]["beam2_on_off_contrast_db"] = round(des["contrast_db"], 2)
p["recalibration"] = {"rule": "offset_k = HFSS gain of the v2 board in state k (design direction) minus the raw model gain of v2 in state k; add the new offset_k to the raw model gains of this board (beam1_gain_dbi_raw_model, beam2_gain_dbi_raw_model) of the same state",
                      "v2_raw_model_gain_at_design_direction_dbi": {"state1": 16.78, "state2": 14.8}, "current_offsets_db": OFF,
                      "example": "calibrated beam 1 in state 2 = raw (%.2f) + offset_2" % raw(S[1], "b1_(25,0)")}
p["sll_definition"] = "highest peak outside the first null of every wanted beam (2 x HPBW if none), relative to the stronger beam; the HFSS gap seen on board_K2_default with this definition is +4.1 dB (model -4.8 / -5.3, HFSS -0.7 / -1.2)"
p["xpd_note"] = "XPD not computed for this variant"
d.joinpath("predictions.json").write_text(json.dumps(p, indent=1, default=str), encoding="utf-8", newline="\n")
with open(d / "states.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["state", "name", "D_sheets", "diodes", "wanted_beams", "polarisation"])
    w.writerow(["1", "default (unpowered)", "NO boundary (or Impedance Rs = 1e6 ohm/sq)", "OFF", "beam 1 (25, 0)", "LP phi-hat (LP-V)"])
    w.writerow(["2", "circuit ON", "one Impedance boundary on all D_*: Rs = 6.1 ohm/sq (R_on 2.44 ohm x W/L = 1.0 / 0.4), Reactance 0", "ON", "beam 1 (25, 0) + beam 2 (25, 180)", "LP phi-hat (LP-V)"])
t = d.joinpath("HFSS_SETUP.md").read_text(encoding="utf-8")
t = re.sub(r"^# HFSS_SETUP.*$", "# HFSS_SETUP - RO3006 board ADD-A-BEAM (b): default (diodes OFF) = beam 1 (25, 0) only; circuit ON = beam 1 (25, 0) + beam 2 (25, 180), ONE diode group `D_*`, LP-V, exactly 10.5 GHz, 150 mm hexagon circumradius. MODEL_ONLY. One page.", t, count=1, flags=re.M)
i = t.index("**Model predictions")
j = t.index("\n", i)
s1, s2 = S
pred = ("**Model predictions (MODEL_ONLY; launch 0.71; calibrated with the board_K2_default offsets, to be re-calibrated with the v2 offsets before running):** "
        f"state 1: beam 1 {g(s1, 'b1_(25,0)'):.1f} dBi (raw {raw(s1, 'b1_(25,0)'):.1f}), level at (25, 180) {g(s1, 'b2_(25,180)'):.1f} dBi, SLL (first-null) {s1['SLL_firstnull_db']:.1f} dB; "
        f"state 2: beam 1 {g(s2, 'b1_(25,0)'):.1f} dBi (raw {raw(s2, 'b1_(25,0)'):.1f}), beam 2 {g(s2, 'b2_(25,180)'):.1f} dBi (raw {raw(s2, 'b2_(25,180)'):.1f}), balance {s2['balance_db']:.1f} dB, SLL {s2['SLL_firstnull_db']:.1f} dB (vs weaker beam {s2['SLL_firstnull_vs_weaker_db']:.1f}); "
        f"**beam-2 on/off contrast {des['contrast_db']:.1f} dB** (level at (25, 180): state 2 minus state 1); 500 diodes (250 switched cells). Falsifiers and the recalibration rule: `predictions.json`.")
t = t[:i] + pred + t[j:]
note = f"""
## Run notes (add-a-beam (b))
**Design: each state's hologram uses its own real beta/k0 (state 1 1.33, state 2 1.26); beam 1 is carried by the static cells in both states, beam 2 only by the switched cells; w2 = 1.25 (beam-2 target weight), eps_B = 0.**
**State runs by Impedance sweep (S_w):** OFF = Impedance boundary Rs = 1e6 ohm/sq on all `D_*` sheets (4e5 ohm across the 0.4 mm gap = open), ON = Rs = 6.1 ohm/sq (reactance 0); one boundary on all `D_*`, one variable `S_w` (1e6, 6.1). 1e6 ohm/sq is equivalent to no boundary.
**SLL definition: the highest peak outside the main lobe's first null (first local minimum along radial cuts from the peak), or 2 x HPBW if no null, of EVERY wanted beam (state 2 has two), relative to the stronger beam** (also given relative to the weaker beam). On board_K2_default the model-to-HFSS gap with this definition is +4.1 dB (model -4.8 / -5.3, HFSS -0.7 / -1.2): expect about 4 dB worse than the model SLL.
**Expected 3 dB beamwidths (model, deg; elevation / azimuth-angular):** state 1 beam 1 {s1['hpbw']['0']['elevation']:.1f} / {s1['hpbw']['0']['azimuth_angular']:.1f}; state 2 beam 1 {s2['hpbw']['0']['elevation']:.1f} / {s2['hpbw']['0']['azimuth_angular']:.1f}, beam 2 {s2['hpbw']['1']['elevation']:.1f} / {s2['hpbw']['1']['azimuth_angular']:.1f} (the model matched the HFSS elevation width of board_K2_default and under-estimated the azimuth width by 1.5-2 deg).
**Re-calibration (before running):** offset_k = HFSS gain of the v2 board in state k at the design direction minus its raw model gain (state 1 16.78 dBi, state 2 14.80 dBi); add it to the raw model gains of this board (`predictions.json`).
**Rim / ring check:** theta 12 deg azimuthal median against `hfss_ring_levels.json` (n = 0, +-6 content).
"""
d.joinpath("HFSS_SETUP.md").write_text(t.rstrip("\n") + "\n" + note, encoding="utf-8", newline="\n")
dst = Path(r"C:\hmsa_sync\data\sidequest\board_RO3006_r2p30_addbeam_b")
dst.mkdir(parents=True, exist_ok=True)
for f in d.iterdir():
    shutil.copy2(f, dst / f.name)
print(sorted(x.name for x in dst.iterdir()))
