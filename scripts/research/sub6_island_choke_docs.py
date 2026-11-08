"""SUB6 island cell, vias / chokes package: predictions.json, states.csv, results_template.csv, HFSS_SETUP.md (reference = the NO-VIA MoM numbers of the island package).  MODEL_ONLY.
    python sub6_island_choke_docs.py
"""
import csv
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sub6_island_choke_package as pk
import sub6_slotdisc_hfss_docs as sd
import sub6_slotdisc_hfss_package as sp

OUT, RES, ROOT = pk.OUT, pk.RES, pk.ROOT
OLD = json.loads((ROOT / "hfss_exports" / "sub6_island_cell" / "predictions.json").read_text())


def rot(p):
    return (-p[1], p[0])  # +90 deg about z: (x, y) -> (-y, x)


def ref_runs():
    """no-via MoM numbers keyed (state, direction, phase) in the NEW frame (across = old y = $ph_x, along = old x = $ph_y)."""
    out = {}
    for r in OLD["runs"]:
        if r["variant"] != "novia":
            continue
        d = "across" if r["direction"] == "across_y" else "along"
        q = r["ph_y_deg"] if d == "across" else r["ph_x_deg"]
        out[(r["state"], d, q)] = (r["f_MoM_GHz"], r["f_MoM_other_branch_GHz"], r["tier"])
    return out


def main():
    fit, chokes = pk.load_chokes()
    ref = ref_runs()
    eq = {"OFF": sp.parallel_equiv(sp.z_off(pk.base.F0), pk.base.F0), "ON": sp.parallel_equiv(sp.z_on(pk.base.F0), pk.base.F0)}
    eq_f = {f"{f:.1f}": {"OFF": sp.parallel_equiv(sp.z_off(f * 1e9), f * 1e9), "ON": sp.parallel_equiv(sp.z_on(f * 1e9), f * 1e9)} for f in (5.2, 5.5, 5.8, 6.1, 6.4)}
    bare = OLD["sanity_bare_substrate_GHz"]
    vi, vd = rot(pk.VIA_I), rot(pk.VIA_D)
    ch_rot = [dict(via=rot(c["via"]), apex=rot(c["apex"]), path=[rot(p) for p in c["path"]], n_runs=c["n_runs"], ls_mm=c["ls"], link_mm=c["link"], lead_mm=c["lead"], meander_gap_mm=c["spacing"], axis_angle_deg=math.degrees(c["axis_angle"]) + 90.0) for c in chokes]
    plan = [("choke", "ISL_choke3.step", "OFF", "across", 95, 1), ("choke", "ISL_choke3.step", "ON", "across", 80, 1),
            ("vias", "ISL_vias2.step", "OFF", "across", 95, 2), ("vias", "ISL_vias2.step", "ON", "across", 80, 2),
            ("choke", "ISL_choke3.step", "OFF", "along", 85, 3), ("choke", "ISL_choke3.step", "ON", "along", 85, 3),
            ("vias", "ISL_vias2.step", "OFF", "along", 85, 4), ("vias", "ISL_vias2.step", "ON", "along", 85, 4)]
    for variant, fn in (("choke", "ISL_choke3.step"), ("vias", "ISL_vias2.step")):
        for st, d, qs in (("OFF", "across", (90, 100)), ("ON", "across", (85, 90)), ("OFF", "along", (80, 90)), ("ON", "along", (80, 90))):
            for q in qs:
                plan.append((variant, fn, st, d, q, 5))
    runs = []
    for variant, fn, st, d, q, pri in plan:
        f0, f1, tier = ref[(st, d, q)]
        runs.append(dict(variant="v3" if variant == "choke" else "v2", state=st, file=fn, direction=d, ph_x_deg=q if d == "across" else 0, ph_y_deg=q if d == "along" else 0, tier=1 if pri <= 4 else 2, priority=pri,
                         f_MoM_noVia_GHz=f0, f_MoM_noVia_other_branch_GHz=f1))
    for n, r in enumerate(runs, 1):
        r["run_id"] = f"C{n:02d}"
    pred = {"labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "PREDICTED_BEFORE_HFSS_RUN"], "frame": "rotated +90 deg about z: slot along y, ACROSS the slot = x = $ph_x, ALONG = $ph_y", "p_mm": pk.P,
            "reference": "the f_MoM values are the NO-VIA MoM numbers of the island package (no vias, no choke, no FR4 layer); the difference HFSS(v2) - MoM and HFSS(v3) - HFSS(v2) at equal phase is the via / choke effect",
            "series_models": {"OFF": "(0.15 pF || 5 kohm) + 0.6 nH", "ON": "3.0 ohm + 0.6 nH"}, "parallel_equivalents_5p8GHz": eq, "parallel_equivalents_vs_f": eq_f,
            "cell": OLD["cell"], "variants": {"v2": "vias (island + disc) with 2.0 x 0.3 mm bias stubs, no choke", "v3": "vias + one 120-deg radial-stub choke per via"},
            "vias_rotated_frame_mm": {"island": vi, "disc": vd, "drill_mm": 0.3, "pad_mm": 0.6, "antipad_mm": 1.0}, "vias_island_frame_mm": {"island": pk.VIA_I, "disc": pk.VIA_D},
            "choke_design": dict(fit["design"], ms=fit["ms"]), "choke_layout_rotated_frame": ch_rot, "choke_gaps_mm": fit["gaps"], "fr4": dict(t_mm=pk.FR4_T, eps_r=pk.FR4_EPS, tan_d=pk.FR4_TAND, assumed=True),
            "runs": runs, "sanity_runs": [dict(run_id="S1", what="bare substrate (ground + RT5880, no metal on top, no FR4), $ph_x = 90, $ph_y = 0", f_model_GHz=bare["90"]["open_air"], priority=1,
                                              note="model 6.291 GHz open air / 6.291 GHz with the 150 mm lid; a different value means a setup error")],
            "modes_to_record": "modes 1-4 for the runs of priority 1 and 2 (f, Q, Ez over the patch decays upward = bound TM-like mode); ON across has two MoM branches (island resonance below, TM mode above)"}
    (OUT / "predictions.json").write_text(json.dumps(pred, indent=1, default=float), encoding="utf-8", newline="\n")
    sh = [[-pk.base.Y_SH0, 0.0, pk.base.Z_SH], [-pk.base.Y_SH1, 0.0, pk.base.Z_SH]]
    with open(OUT / "states.csv", "w", encoding="utf-8", newline="\n") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["variant", "state", "step_file", "DIODE_SHEET", "R_ohm", "C_pF", "L_nH", "series_model", "note"])
        for variant, fnm in (("v2", "ISL_vias2.step"), ("v3", "ISL_choke3.step")):
            w.writerow([variant, "OFF", fnm, "Lumped RLC parallel", f"{eq['OFF']['R_ohm']:.4g}", f"{eq['OFF']['C_pF']:.4g}", "unticked", "(0.15 pF || 5 kohm) + 0.6 nH", f"5.8 GHz parallel equivalent; current line {sh[0]} -> {sh[1]}"])
            w.writerow([variant, "ON", fnm, "Lumped RLC parallel", f"{eq['ON']['R_ohm']:.4g}", "unticked", f"{eq['ON']['L_nH']:.4g}", "3.0 ohm + 0.6 nH", "5.8 GHz parallel equivalent; same geometry as OFF"])
    with open(OUT / "results_template.csv", "w", encoding="utf-8", newline="\n") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["run_id", "priority", "variant", "state", "step_file", "direction", "ph_x_deg", "ph_y_deg", "tier", "f_MoM_noVia_GHz", "f_MoM_noVia_other_branch_GHz", "mode1_f_GHz", "mode1_Q", "mode2_f_GHz", "mode2_Q", "mode3_f_GHz", "mode3_Q", "mode4_f_GHz", "mode4_Q",
                    "bound_TM_mode_no", "Ez_decays_from_patch_YN", "note"])
        w.writerow(["S1", 1, "bare", "none", "(ground + RT5880 only)", "across", 90, 0, 1, f"{bare['90']['open_air']:.3f}"] + [""] * 12)
        for r in runs:
            w.writerow([r["run_id"], r["priority"], r["variant"], r["state"], r["file"], r["direction"], r["ph_x_deg"], r["ph_y_deg"], r["tier"], f"{r['f_MoM_noVia_GHz']:.3f}",
                        "" if r["f_MoM_noVia_other_branch_GHz"] is None else f"{r['f_MoM_noVia_other_branch_GHz']:.3f}"] + [""] * 12)
    setup(pred, eq, bare, ch_rot, fit)
    print("written; runs", len(runs))


def setup(pred, eq, bare, ch_rot, fit):
    d = fit["design"]
    light90 = sd.C0 * math.radians(90) / (pk.P * 1e-3) / (2 * math.pi) / 1e9
    ms = fit["ms"]
    sh0, sh1 = (-pk.base.Y_SH0, 0.0, pk.base.Z_SH), (-pk.base.Y_SH1, 0.0, pk.base.Z_SH)
    vi, vd = pred["vias_rotated_frame_mm"]["island"], pred["vias_rotated_frame_mm"]["disc"]
    cs = "\n".join(f"* {nm}: via ({c['via'][0]:.3f}, {c['via'][1]:.3f}); bias trace 0.3 mm wide, {c['lead_mm']} mm long from the via centre; choke line 0.15 mm wide, {c['n_runs']} meander run(s) of {c['ls_mm']:.3f} mm, meander gap {c['meander_gap_mm']} mm, final straight link {c['link_mm']:.3f} mm; stub apex at ({c['apex'][0]:.3f}, {c['apex'][1]:.3f}), 120 deg sector r {d['r2_mm']:.3f} mm, axis {c['axis_angle_deg']:.0f} deg."
                   for nm, c in zip(("island via", "disc via"), ch_rot))
    txt = f"""# HFSS eigenmode: SUB6 island cell, bias vias without (v2) and with (v3) radial-stub chokes (RT5880 3.175 mm, p 11.6, 5.8 GHz). EXPORT STANDARD (2026-11-05). Predictions MODEL_ONLY (`predictions.json`), NOT_HFSS_VALIDATED.

Cell `ISL_W2.9_q0.8_Li0.5_c0.6` (island 4.348 x 1.300 mm in a 8.696 x 2.900 mm slot, ONE `DIODE_SHEET` across the 1.00 mm gap), same top metal as the no-via package. **Frame: rotated +90 deg about z, slot along y, ACROSS the slot = x = `$ph_x`, ALONG = `$ph_y`.**
**1. Assign material `pec` to `METAL_PEC` (3 solids: disc patch with its via and choke / stub, floating island with its via and choke / stub, ground with two anti-pads). STEP imports as vacuum.**
**2. Air box `AIR`: one vacuum box, ONE face per side wall, z = {pk.Z_AIR} ... +{pk.LID:g} (30 mm below the bottom copper, 150 mm above the patch).** It overlaps all solids on purpose; keep `AIR` first. **Boundaries: lattice (periodic) pairs on the four SIDE faces of `AIR` only (x = -5.8 / +5.8 with `$ph_x`, y = -5.8 / +5.8 with `$ph_y`); PEC (Perfect E) on the TOP face (z = +{pk.LID:g}); RADIATION boundary on the BOTTOM face (z = {pk.Z_AIR}).** (In the earlier no-via package the unassigned bottom face was PEC.)
**3.** `SUBSTRATE`: RT5880 eps_r 2.2, tan d 0.0009 (z -3.175 ... 0). `SUBSTRATE_FR4` (z {pk.Z_FR4B} ... {pk.Z_GB}): **FR-4 {pk.FR4_T} mm, eps_r {pk.FR4_EPS}, tan d {pk.FR4_TAND} (typical values ASSUMED: replace with your laminate).** Copper: patch 35 um (z 0 ... 0.035), ground 35 um (z -3.210 ... -3.175), bottom copper 35 um (z {pk.Z_BOT} ... {pk.Z_FR4B}).
**4.** Project Variables `$ph_x`, `$ph_y` (deg). **ACROSS: `$ph_x = phi`, `$ph_y = 0`; ALONG: `$ph_x = 0`, `$ph_y = phi`.**
**5. Files / states** (`states.csv`): `ISL_vias2.step` (v2) and `ISL_choke3.step` (v3); the geometry differs only below the ground (v2: 2.0 x 0.3 mm bias stubs; v3: chokes).
| state | `DIODE_SHEET` (Lumped RLC, PARALLEL) |
|---|---|
| OFF | **R = {eq['OFF']['R_ohm'] / 1e3:.3f} kohm, C = {eq['OFF']['C_pF']:.4f} pF, L unticked** |
| ON | **R = {eq['ON']['R_ohm']:.1f} ohm, L = {eq['ON']['L_nH']:.4f} nH, C unticked** |

(OFF series model (0.15 pF || 5 kohm) + 0.6 nH = {eq['OFF']['Z_series'][0]:.2f} {eq['OFF']['Z_series'][1]:+.2f}j ohm; ON 3.0 ohm + 0.6 nH = {eq['ON']['Z_series'][0]:.2f} {eq['ON']['Z_series'][1]:+.2f}j ohm; exact at 5.8 GHz only, other frequencies in `predictions.json`.)
**Sheet `DIODE_SHEET`: 1.00 (across, x: -1.450 ... -0.450) x 0.35 (along, y: -0.175 ... +0.175) mm at z = +0.0175.** Current line across the gap: ({sh0[0]:g}, {sh0[1]:g}, {sh0[2]:g}) -> ({sh1[0]:g}, {sh1[1]:g}, {sh1[2]:g}).
**6. Vias** (drill 0.3, pad 0.6, anti-pad 1.0 in the ground; positions = minima of the RF surface current of the MoM, both directions, OFF and ON): island via ({vi[0]:.3f}, {vi[1]:.3f}), disc via ({vd[0]:.3f}, {vd[1]:.3f}) (island-frame values -1.78, -0.53 and -3.63, -3.63). Top pads lie inside the metal (omitted); the barrels (r 0.15) run from the metal down to the bottom copper through the anti-pads (0.35 mm radial clearance to the ground).
**7. Chokes (v3 only)**, analytic design (`sub6_island_choke.py`): choke line 0.15 mm, Z0 {ms['z0']:.1f} ohm, eps_eff {ms['eps_eff']:.3f}, lambda_g/4 {ms['quarter_wave_mm']:.2f} mm, line length {d['l_mm']:.3f} mm meandered inside the cell; ONE 120 deg radial stub r {d['r2_mm']:.3f} mm at the line end (series resonance = RF short at 5.8 GHz); modelled |Z_in| >= {d['zin_min_ohm']:.0f} ohm over 5.725-5.875 GHz (single stub, not HFSS), via barrel {d['via_L_nH']:.2f} nH. Layout (rotated frame, mm):
{cs}
Gaps (mm): between the two chokes {fit['gaps']['between']:.3f}; to the periodic images {fit['gaps']['images']:.3f}; smallest self gap inside a choke {fit['gaps']['self']:.3f}; nothing crosses the cell boundary.
**8.** Eigenmode: min frequency 4.5 GHz, **8 modes** (more modes than the earlier packages: FR4, stubs and chokes add resonances; the radiation boundary makes Q finite), delta f 0.1 %, max 15 passes; initial mesh 1 mm on `SUBSTRATE`, 0.5 mm on the patch faces, 0.1 mm on `DIODE_SHEET`, the gaps, the via barrels and the anti-pads, 0.1 mm on the choke lines, `SUBSTRATE_FR4` 0.5 mm. **Clean Up Solutions before each run.**
**9. What to record: modes 1-4 for the priority 1 and 2 runs** (f, Q, Ez over the patch; the bound TM-like mode has Ez large over the patch and decaying upward; light line at 90 deg {light90:.3f} GHz). Run order = `priority` in `results_template.csv`: **priority 1: C01 v3 OFF across 95 deg and C02 v3 ON across 80 deg; priority 2: C03 / C04 the same for v2; then ALONG (85 deg) for v3 and v2; priority 5 = +-5 deg.**
**Sanity S1:** bare substrate (ground + RT5880 only), `$ph_x = 90`: model {bare['90']['open_air']:.3f} GHz open air / {bare['90']['lid_150mm']:.3f} GHz with a 150 mm lid; with the radiation boundary the value is unchanged to about 0.01 GHz.
**Reference numbers.** `f_MoM_noVia_GHz` in the run list = the NO-VIA MoM prediction at that phase (ON across: upper branch; the lower branch is `f_MoM_noVia_other_branch_GHz`, the island resonance). Compare with the no-via HFSS runs of the island package (I01-I03) if available: HFSS(v2) - HFSS(no via) = via + FR4 + stub effect, HFSS(v3) - HFSS(v2) = choke effect. The MoM geometry is the mesh-64 mask (clearance 0.544, gap 0.906, island 4.350 x 1.450) while HFSS uses nominal 0.600 / 1.000 / 4.348 x 1.300: expect 0.1-0.3 GHz offsets. Not run in AEDT by the author.
"""
    (OUT / "HFSS_SETUP.md").write_text(txt, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
