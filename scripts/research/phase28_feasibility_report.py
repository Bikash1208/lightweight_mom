"""Phase 28 item 2: FEASIBILITY_REPORT.md generated from `results/phase28/feasibility_results.json` ONLY (every number in the report is a value of that
file; `tests/test_phase28_feasibility_report.py` checks it). Prose carries no numbers of its own.

    python scripts/research/phase28_feasibility_report.py [results.json] [out.md]
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results" / "phase28"
CPV, LPV, LPXV = (
    ("cp_cont", "cp_1bit", "cp_2bit"),
    ("lp_cont", "lp_1bit", "lp_2bit"),
    ("lpx_cont", "lpx_1bit", "lpx_2bit"),
)


def f(v, nd=2):
    return "n/a" if v is None else f"{v:.{nd}f}"


def bound(b):
    """An AR / XPD record: the number, or the bound text."""
    return b["text"] if b["value"] is None else f(b["value"], 2)


class Doc:
    def __init__(self, d):
        self.d = d
        self.cases = {(c["pitch_mm"], c["aperture_r_mm"], c["theta0_deg"], c["phi0_deg"]): c for c in d["cases"]}
        self.k = d["constants"]

    def case(self, p, r, t, ph):
        return self.cases[(p, r, t, ph)]

    def pitches(self):
        return sorted({c["pitch_mm"] for c in self.d["cases"]})


def table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return out


def build(d):
    D = Doc(d)
    k = D.k
    L = []
    a = L.append
    a("# Phase 28 item 2 - feasibility of a coded monopole-fed hologram (stop point)")
    a("")
    a(
        "**Labels:** "
        + ", ".join(f"`{x}`" for x in d["labels"])
        + ". **Generated from** `results/phase28/feasibility_results.json` only."
    )
    a("")
    a(
        f"Forward model: the Phase 26 tensor-CP model unchanged (supercell table, ray march, aperture far field on the grounded slab). Frequency {f(k['frequency_ghz'])} GHz, the only frequency of the supercell table; "
        f"{f(k['frequency_ghz_not_available'], 1)} GHz is `{d['frequency_10p0_ghz']}`. Cells are idealised (a state is a tensor or scalar value, no geometry); the table, computed for 6.4 mm cells, is assumed at every pitch, "
        "so the pitch sweep changes only how the hologram is sampled. Gain is referred to the accepted power with the Phase 26 bare launch fraction "
        f"{f(k['launch_fraction_bare'], 3)}; quantisation changes the aperture directivity only (the conversion efficiency is held at the continuous design's value, an assumption)."
    )
    a("")
    # ------------------------------------------------------------------ anchor
    an = d["anchor"]
    a("## Anchor: no quantisation reproduces the Phase 26 single-hexagon prediction")
    a("")
    a(f"Source: {an['source']}. Same grid, same code path; the finer refinement of the peak is reported separately.")
    a("")
    rows = [
        (
            r["quantity"],
            f(r["this_pipeline"], 6),
            f(r["phase26"], 6),
            f"{r['abs_diff']:.1e}",
            f"{r['tolerance']:.0e}",
            "yes" if r["within_tolerance"] else "NO",
        )
        for r in an["rows"]
    ]
    L += table(["quantity", "this pipeline", "Phase 26", "abs diff", "tolerance", "within"], rows)
    a("")
    a(
        f"Refined peak directivity {f(an['refined_d_ap_dbi'], 3)} dBi ({f(an['refinement_gain_over_coarse_db'], 3)} dB above the coarse-grid peak), pointing error {f(an['refined_pointing_error_deg'], 2)} degrees."
    )
    a("")
    # ------------------------------------------------------------------ contrast
    c = d["contrast"]
    a("## Switching contrast, first")
    a("")
    a(
        "**Definition.** The contrast is what the cell must be able to change so that the coded hologram can leak as strongly as the continuous one. Four measures, all at the design frequency:"
    )
    a("")
    a(
        "- `alpha_cap`: the largest local leakage constant [Np/m] the coded surface can realise. The continuous Phase 26 design asks for a profile alpha(x) up to the supercell table's maximum; a coded surface that cannot reach it converts less of the launched surface wave."
    )
    a(
        "- Tensor (CP, rotating axis) cell: the fixed anisotropy `delta` [mS] of `Y = Ybar I + j delta R(a) diag(+,-) R(a)^T`, with alpha equal to a constant times delta squared, from the Phase 26 closed form (`leakage_alpha_first_order_tensor`). `dX_aniso` [ohm] is the difference of the TM surface reactances `X = Im(1/(Y + Y_slab,TM(beta)))` seen along the two principal axes; the switch rotates the axis and does not change it. The largest state-to-state change of the sheet tensor (axes perpendicular) is twice delta."
    )
    a(
        "- Scalar cell (also the amplitude control of the `lpx` bound), states `Ybar +- j dB/2` (`dB` the full swing between the extreme levels): the first harmonic of the quantised surface is c1, one quarter of dB times g (g the first-harmonic gain of the level set), and alpha is a constant times c1 squared (`leakage_alpha_first_order`, TM and TE parts, the larger `dB` reported). `dX_switch` [ohm] is the TM surface-reactance difference of the two extreme states; `d_beta` the difference of the TM0 phase constants of the two uniform states."
    )
    a(
        "- Modulation depth: `delta/|Y_down|` (tensor) or half of dB over `|Y_down|` (scalar), `Y_down = Ybar + Y_slab,TM(beta)`, the Phase 26 `required_map` definition."
    )
    a("")
    a(
        f"**Method.** The supercell table's leakage constants are capped at `alpha_cap` and the single-hexagon CP design is re-marched; `alpha_cap` is the value at which the gain falls {f(k['gain_loss_tolerance_db'], 1)} dB below the uncapped design "
        f"(table maximum {f(c['table_alpha_max_np_per_m'], 3)} Np/m, minimum {f(c['table_alpha_min_np_per_m'], 3)} Np/m). It is converted to physical units with the Phase 26 closed forms at the mean `beta/k0` of the design's cells; "
        "the scalar case is cross-checked with the exact truncated-Floquet complex beta."
    )
    a("")
    rows = []
    for s in c["steerings"]:
        t1, s1, s2 = (
            s["contrast"]["bits_1"]["tensor"],
            s["contrast"]["bits_1"]["scalar"],
            s["contrast"]["bits_2"]["scalar"],
        )
        rows.append(
            (
                f"{f(s['theta0_deg'], 0)} / {f(s['phi0_deg'], 0)}",
                f(s["alpha_cap_for_1db_np_per_m"], 3),
                f(s["check_at_cap_star"]["loss_vs_uncapped_db"], 3),
                f(s["contrast"]["bits_1"]["beta_over_k0_mean"], 4),
                f(t1["delta_ms"], 3),
                f(t1["dx_aniso_ohm"], 1),
                f(t1["d_beta_over_k0_principal_axes"], 4),
                f(t1["modulation_depth_delta_over_ydown"], 3),
                f(s1["d_b_full_swing_ms"], 3),
                f(s1["dx_switch_ohm"], 1),
                f(s1["d_beta_over_k0_uniform_states"], 4),
                f(s2["d_b_full_swing_ms"], 3),
                f(s2["dx_switch_ohm"], 1),
            )
        )
    L += table(
        [
            "beam theta / phi [deg]",
            "alpha_cap [Np/m]",
            "loss re-check [dB]",
            "beta/k0 mean",
            "tensor delta [mS]",
            "tensor dX_aniso [ohm]",
            "tensor d_beta/k0",
            "tensor depth",
            "LP binary dB [mS]",
            "LP dX_switch [ohm]",
            "LP d_beta/k0",
            "LP four-level dB [mS]",
            "LP four-level dX_switch [ohm]",
        ],
        rows,
    )
    a("")
    a(
        "Exact-Floquet cross-check of the scalar closed form (binary sheet at the closed-form full swing found for each beam; TM, collinear):"
    )
    a("")
    rows = []
    for s in c["steerings"]:
        fc = s["floquet_check"]
        rows.append(
            (
                f"{f(s['theta0_deg'], 0)} / {f(s['phi0_deg'], 0)}",
                f(fc["alpha_first_order_target_np_per_m"], 3),
                f(fc["alpha_floquet_exact_np_per_m"], 3),
                f(fc["ratio_exact_over_first_order"], 3),
                fc.get("note", "-"),
            )
        )
    L += table(
        [
            "beam theta / phi [deg]",
            "alpha closed form [Np/m]",
            "alpha Floquet exact [Np/m]",
            "ratio exact/closed",
            "note",
        ],
        rows,
    )
    a("")
    dt = d["diode_terminal"]
    a(
        "For scale only (the cell-level contrast depends on how the diode is embedded, which is item 3): the terminal reactance swing of the diode models at the design frequency."
    )
    a("")
    rows = [
        (n, f(v["r_on_ohm"], 2), f(v["z_on_ohm"][1], 1), f(v["z_off_ohm"][1], 1), f(v["delta_x_off_minus_on_ohm"], 1))
        for n, v in dt.items()
        if n != "note"
    ]
    L += table(
        ["diode (datasheet, PROVISIONAL)", "R_on [ohm]", "X_on [ohm]", "X_off [ohm]", "X_off - X_on [ohm]"], rows
    )
    a("")
    a(dt["note"] + ".")
    a("")
    a(
        "**Table range caveat.** The supercell table covers a limited range of modulation wavenumbers; the share of the design's march points outside it is given in the scan table below and the contrast rows for steering beyond the table are `TABLE_EXTRAPOLATED`."
    )
    a("")
    # ------------------------------------------------------------------ design point
    a("## Design point: 150 mm hexagon, beam theta 30 / phi 225 degrees")
    a("")
    a(
        "Quantisation loss is relative to the continuous design of the same family and pitch. `cp` = rotating-slot tensor, RHCP (the Phase 26 family); `lp` = scalar surface, current along the radial direction; `lpx` = linear polarisation with IDEAL polarisation control (current along x in every cell), an upper bound no cell has been shown to realise. AR and XPD beyond the bound are shown as bounds, never as numbers."
    )
    a("")
    rows = []
    for p in D.pitches():
        cs = D.case(p, 150.0, 30.0, 225.0)
        for v in CPV + LPV + LPXV:
            r = cs["variants"][v]
            xp = bound(r["xpd_at_peak"]) if "xpd_at_peak" in r else "-"
            ar = bound(r["ar_at_peak"]) if "ar_at_peak" in r else "-"
            rows.append(
                (
                    f(p, 1),
                    v,
                    f(r["gain_dbi"], 2),
                    f(r["quantisation_loss_db"], 2),
                    f(r["sll_db"], 1),
                    f(r["pointing_error_deg"], 2),
                    xp,
                    ar,
                )
            )
    L += table(
        [
            "pitch [mm]",
            "variant",
            "gain [dBi]",
            "quantisation loss [dB]",
            "SLL [dB]",
            "pointing error [deg]",
            "XPD [dB]",
            "AR [dB]",
        ],
        rows,
    )
    a("")
    a(
        f"Peak search: {f(k['peak_refine_step_deg'], 2)} degree direction-cosine refinement around the coarse-grid peak; the sidelobe level uses the strongest coarse-grid sidelobe refined the same way. Gain is the aperture directivity plus ten times the base-ten logarithm of (launch fraction x eta_conv); eta_conv of this design is given in the scan table."
    )
    a("")
    # ------------------------------------------------------------------ LP finding
    cs = D.case(6.4, 150.0, 30.0, 225.0)
    lp, lpx = cs["variants"]["lp_cont"], cs["variants"]["lpx_cont"]
    a("## Finding: a scalar surface on a monopole cannot form a linear-polarised pencil beam")
    a("")
    a(
        f"With the induced current along the surface wave's own (radial) direction the continuous scalar design has a pointing error of {f(lp['pointing_error_deg'], 2)} degrees and a strongest sidelobe at {f(lp['sll_db'], 1)} dB: "
        f"the aperture field is radially polarised, so the contributions cancel in the beam direction (a ring or split beam, not a pencil). Gain {f(lp['gain_dbi'], 2)} dBi against {f(lpx['gain_dbi'], 2)} dBi for the same hologram with ideal polarisation control (`lpx`) "
        f"and {f(cs['variants']['cp_cont']['gain_dbi'], 2)} dBi for the CP tensor design. The `lp` rows above are therefore NOT a candidate; the bit-depth conclusion for LP rests on `lpx`, which needs a cell that sets the current direction (not designed)."
    )
    a("")
    # ------------------------------------------------------------------ scan
    a("## Steering: pitch " + f(6.4, 1) + " mm, 150 mm hexagon, CP tensor")
    a("")
    a(
        "Gain by beam polar angle (azimuth 225 degrees; broadside has no azimuth). Scan loss is relative to broadside of the same variant; a negative value means the gain is higher than at broadside (the conversion efficiency rises with the steering angle in this model). `table outside` is the share of the design's march points whose modulation wavenumber lies outside the supercell table: results with a large share are `TABLE_EXTRAPOLATED`."
    )
    a("")
    rows = []
    for t in sorted({c_["theta0_deg"] for c_ in d["cases"] if c_["pitch_mm"] == 6.4 and c_["aperture_r_mm"] == 150.0}):
        ph = 0.0 if t == 0.0 else 225.0
        cs = D.case(6.4, 150.0, t, ph)
        v = cs["variants"]
        spread = [
            D.case(6.4, 150.0, t, q)["variants"]["cp_2bit"]["gain_dbi"]
            for q in (0.0, 45.0, 225.0)
            if (6.4, 150.0, t, q) in D.cases
        ]
        rows.append(
            (
                f(t, 0),
                f(cs["eta_conv"], 3),
                f(cs["table_extrapolated_fraction"], 2),
                f(v["cp_cont"]["gain_dbi"], 2),
                f(v["cp_2bit"]["gain_dbi"], 2),
                f(v["cp_1bit"]["gain_dbi"], 2),
                f(v["cp_2bit"]["scan_loss_db"], 2),
                f(v["cp_2bit"]["sll_db"], 1),
                f(v["cp_2bit"]["pointing_error_deg"], 2),
                f(min(spread), 2) + " to " + f(max(spread), 2),
            )
        )
    L += table(
        [
            "theta0 [deg]",
            "eta_conv",
            "table outside",
            "cp cont [dBi]",
            "cp 2-bit [dBi]",
            "cp 1-bit [dBi]",
            "2-bit scan loss [dB]",
            "2-bit SLL [dB]",
            "2-bit pointing [deg]",
            "2-bit gain over phi0 = 0, 45, 225 [dBi]",
        ],
        rows,
    )
    a("")
    # ------------------------------------------------------------------ pitch and aperture
    a("## Pitch and aperture")
    a("")
    rows = []
    for p in D.pitches():
        cs = D.case(p, 150.0, 30.0, 225.0)
        v = cs["variants"]
        rows.append(
            (
                f(p, 1),
                cs["n_cells"],
                f(v["cp_cont"]["gain_dbi"], 2),
                f(v["cp_2bit"]["quantisation_loss_db"], 2),
                f(v["cp_1bit"]["quantisation_loss_db"], 2),
                f(v["lpx_1bit"]["quantisation_loss_db"], 2),
                f(v["lpx_2bit"]["quantisation_loss_db"], 2),
            )
        )
    L += table(
        [
            "pitch [mm]",
            "cells",
            "cp continuous [dBi]",
            "cp 2-bit loss [dB]",
            "cp 1-bit loss [dB]",
            "lpx 1-bit loss [dB]",
            "lpx 2-bit loss [dB]",
        ],
        rows,
    )
    a("")
    a("The larger reference aperture (hexagon of the larger circumradius), pitch 6.4 mm:")
    a("")
    rows = []
    for t in (0.0, 30.0, 50.0):
        ph = 0.0 if t == 0.0 else 225.0
        cs = D.case(6.4, 225.0, t, ph)
        v = cs["variants"]
        rows.append(
            (
                f(cs["aperture_r_mm"], 0),
                f(t, 0),
                cs["n_cells"],
                f(cs["eta_conv"], 3),
                f(cs["table_extrapolated_fraction"], 2),
                f(v["cp_cont"]["gain_dbi"], 2),
                f(v["cp_2bit"]["gain_dbi"], 2),
                f(v["cp_1bit"]["gain_dbi"], 2),
                f(v["cp_2bit"]["sll_db"], 1),
            )
        )
    L += table(
        [
            "circumradius [mm]",
            "theta0 [deg]",
            "cells",
            "eta_conv",
            "table outside",
            "cp cont [dBi]",
            "cp 2-bit [dBi]",
            "cp 1-bit [dBi]",
            "cp 2-bit SLL [dB]",
        ],
        rows,
    )
    a("")
    # ------------------------------------------------------------------ hardware
    a("## Control hardware (150 mm hexagon, beam theta 30 / phi 225)")
    a("")
    a(
        f"Counting rules (one method for every row): diodes per cell = bits (one diode per bit, an ASSUMPTION; the topology is item 3); bit = 1 is the forward-biased diode at {f(k['diode_forward_current_ma'], 0)} mA, the share of forward bits is counted from the actual coded states; "
        f"DC power = (forward diodes) x I_F x V_rail, V_rail = {f(k['supply_rail_v'], 1)} V (the supply power, series resistor included; the diode alone dissipates {f(k['diode_forward_voltage_v'], 2)} V x I_F). The same formula is evaluated at {f(k['diode_forward_current_low_ma'], 0)} mA. "
        "Independent bias lines = diodes (arbitrary patterns need an independent state per diode); sharing a line among several diodes constrains the pattern (the line-count columns below give the count only; four cells per line is the hologram sampled at twice the pitch, see the pitch table). "
        f"Shift registers: eight- or sixteen-output devices in chains of eight clocked at {f(k['shift_register_clock_mhz'], 0)} MHz; FPGA pins = chains plus the shared clock, latch and enable lines; output current per pin is the diode forward current."
    )
    a("")
    rows = []
    for p in D.pitches():
        cs = D.case(p, 150.0, 30.0, 225.0)
        for v in ("cp_1bit", "cp_2bit"):
            b = cs["bias"][v]
            s8, s16 = b["shift_register"]["outputs_per_chip_8"], b["shift_register"]["outputs_per_chip_16"]
            rows.append(
                (
                    f(p, 1),
                    v,
                    b["n_cells"],
                    b["n_diodes"],
                    b["independent_bias_lines"],
                    b["bias_line_sharing"]["cells_per_line_2"],
                    b["bias_line_sharing"]["cells_per_line_4"],
                    s8["chips"],
                    s16["chips"],
                    s8["fpga_io_pins"],
                    f(s8["update_time_us_at_10mhz"], 1),
                    f(b["fraction_bits_forward"], 3),
                    f(b["dc_current_expected_a"], 1),
                    f(b["dc_power_expected_w"], 1),
                    f(b["dc_at_lower_bias"]["1mA"]["dc_power_expected_w"], 2),
                    f(b["dc_power_all_on_w"], 1),
                )
            )
    L += table(
        [
            "pitch [mm]",
            "variant",
            "cells",
            "diodes",
            "bias lines",
            "lines, 2 per line",
            "lines, 4 per line",
            "eight-output chips",
            "sixteen-output chips",
            "FPGA pins (eight-output)",
            "update [us]",
            "forward share",
            "DC current [A]",
            f"DC power expected, {f(k['diode_forward_current_ma'], 0)} mA [W]",
            f"DC power expected, {f(k['diode_forward_current_low_ma'], 0)} mA [W]",
            f"DC power all on, {f(k['diode_forward_current_ma'], 0)} mA [W]",
        ],
        rows,
    )
    a("")
    b = D.case(6.4, 150.0, 30.0, 225.0)["bias"]["cp_2bit"]
    a(
        f"Driver current per output: {f(b['driver_current_per_output_ma'], 0)} mA (a sixteen-output device with every output on draws {f(b['shift_register']['outputs_per_chip_16']['chip_current_all_on_ma'], 0)} mA, an eight-output device {f(b['shift_register']['outputs_per_chip_8']['chip_current_all_on_ma'], 0)} mA): the output stage must source or sink this per pin; whether common logic-level shift registers are specified for it is NOT checked here (datasheets to be checked in item 3), a current-driver class may be needed. The LP rows are identical in count to the CP rows of the same bit depth (`lp` and `lpx` bias tables are in the JSON)."
    )
    a("")
    a(
        "The lower-bias column is only a power figure; the diode series resistance rises at low bias (item 3 must trade it against insertion loss)."
    )
    a("")
    # ------------------------------------------------------------------ recommendation
    a("## Recommendation and trade-offs")
    a("")
    cs = D.case(6.4, 150.0, 30.0, 225.0)
    v = cs["variants"]
    a(
        "- **Polarisation scope: CP, the rotating-slot tensor family.** It is the family the existing forward model covers, with a single HFSS check of the continuous hexagon (one board, one run); the scalar LP surface does not form a pencil beam (finding above), and the `lpx` bound needs a cell that controls the current direction, which no design exists for."
    )
    a(
        f"- **Bit depth: 2-bit (four slot axes, equally spaced).** At pitch 6.4 mm the 2-bit loss is {f(v['cp_2bit']['quantisation_loss_db'], 2)} dB against {f(v['cp_1bit']['quantisation_loss_db'], 2)} dB for 1-bit; 1-bit costs one diode per cell fewer but loses that much gain and its sidelobe level degrades to {f(v['cp_1bit']['sll_db'], 1)} dB (2-bit: {f(v['cp_2bit']['sll_db'], 1)} dB)."
    )
    b2 = cs["bias"]["cp_2bit"]
    a(
        f"- **Pitch: 6.4 mm** (the pitch of the supercell table and of the HFSS-checked hexagon). In this idealised model the pitch hardly matters (2-bit loss from {f(D.case(4.0, 150.0, 30.0, 225.0)['variants']['cp_2bit']['quantisation_loss_db'], 2)} dB at 4.0 mm to {f(D.case(8.0, 150.0, 30.0, 225.0)['variants']['cp_2bit']['quantisation_loss_db'], 2)} dB at 8.0 mm), so the choice is set by diode count and by what item 3 can build; 8.0 mm is the alternative worth testing with a real cell ({D.case(8.0, 150.0, 30.0, 225.0)['bias']['cp_2bit']['n_diodes']} diodes against {b2['n_diodes']})."
    )
    a(
        f"- **Diode count at the recommendation: {b2['n_diodes']} diodes ({b2['n_cells']} cells x {b2['diodes_per_cell']}), {b2['independent_bias_lines']} independent bias lines**, {b2['shift_register']['outputs_per_chip_16']['chips']} sixteen-output shift registers, {b2['shift_register']['outputs_per_chip_16']['fpga_io_pins']} FPGA pins at the stated chain length; expected DC power {f(b2['dc_power_expected_w'], 1)} W at {f(k['diode_forward_current_ma'], 0)} mA per forward diode ({f(b2['dc_at_lower_bias']['1mA']['dc_power_expected_w'], 2)} W at {f(k['diode_forward_current_low_ma'], 0)} mA). The bias power, not the diode count, is the sharpest cost of this architecture."
    )
    s30 = c["steerings"][1]
    a(
        f"- **Required contrast for the recommendation (beam {f(s30['theta0_deg'], 0)} degrees):** tensor anisotropy {f(s30['contrast']['bits_2']['tensor']['delta_ms'], 3)} mS, i.e. a principal-axis reactance difference of {f(s30['contrast']['bits_2']['tensor']['dx_aniso_ohm'], 1)} ohm, present in every state, with the axis rotating between four orientations; the equivalent scalar swing would be {f(s30['contrast']['bits_2']['scalar']['dx_switch_ohm'], 1)} ohm. These are sheet-level numbers; whether a diode-loaded subwavelength cell reaches them is the question of item 3."
    )
    a(
        f"- **Scan range:** the model is supported by the supercell table up to about the {f(30.0, 0)} degree beam (the share of table-extrapolated march points is in the scan table); beyond it the share of table-extrapolated march points rises (scan table) and the scan losses are `TABLE_EXTRAPOLATED`."
    )
    a("")
    a("## What is not available or not established")
    a("")
    a(
        f"- {f(k['frequency_ghz_not_available'], 1)} GHz: `{d['frequency_10p0_ghz']}`; a table at that frequency is a server job, not run."
    )
    a(
        "- Everything is `MODEL_ONLY`: idealised cells, the table assumed at all pitches, no coupling between neighbouring states, no diode loss, the unconverted surface wave and the feed's direct radiation assumed outside the beam (Phase 26 assumptions)."
    )
    a(
        "- The literature table (item 1) does not support any of the contrast, bit-depth or bias statements here; they are model results."
    )
    a("")
    a("## Stop point")
    a("")
    a(
        "No cell has been designed. Approve or change the pitch, bit depth, polarisation scope and diode count above before item 3 starts."
    )
    return "\n".join(L) + "\n"


def main(src=None, out=None):
    src = Path(src or RES / "feasibility_results.json")
    out = Path(out or RES / "FEASIBILITY_REPORT.md")
    text = build(json.loads(src.read_text(encoding="utf-8")))
    out.write_text(text, encoding="utf-8", newline="\n")
    print("wrote", out)


if __name__ == "__main__":
    main(*sys.argv[1:3])
