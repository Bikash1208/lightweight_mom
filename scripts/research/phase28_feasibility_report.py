"""Phase 28 item 2: FEASIBILITY_REPORT.md generated from `results/phase28/feasibility_results.json` ONLY. Every number printed is a value of that file (the
design point, the recommended pitch and bit depth, the scan azimuth ... are in its `settings` block); `tests/test_phase28_feasibility_report.py` checks this.

    python scripts/research/phase28_feasibility_report.py [results.json] [out.md]
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results" / "phase28"
CPV = ("cp_cont", "cp_1bit", "cp_2bit", "cp_3bit")
LPV = ("lp_cont", "lp_1bit", "lp_2bit")
LPXV = ("lpx_cont", "lpx_1bit", "lpx_2bit")


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
        self.s = d["settings"]
        dp = self.s["design_point"]
        self.dp = dp
        self.dpcase = self.case(dp["pitch_mm"], dp["aperture_r_mm"], dp["theta0_deg"], dp["phi0_deg"])

    def case(self, p, r, t, ph):
        return self.cases[(p, r, t, ph)]

    def at_dp(self, pitch):
        dp = self.dp
        return self.case(pitch, dp["aperture_r_mm"], dp["theta0_deg"], dp["phi0_deg"])

    def pitches(self):
        return sorted({c["pitch_mm"] for c in self.d["cases"]})

    def outside(self, theta, phi):
        dp = self.dp
        return self.case(dp["pitch_mm"], dp["aperture_r_mm"], theta, phi)["table_extrapolated_fraction"]

    def status(self, frac):
        return "TABLE_EXTRAPOLATED" if frac > self.k["table_extrapolated_flag_fraction"] else "inside the table"


def table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return out


def build(d):
    D = Doc(d)
    k, s, dp, rec = D.k, D.s, D.dp, D.s["recommended"]
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
        f"so the pitch sweep changes only how the hologram is sampled. Gain is referred to the accepted power with the Phase 26 bare launch fraction {f(k['launch_fraction_bare'], 3)} "
        f"(the calibrated launch fraction {f(k['launch_fraction_calibrated'], 2)} of the four-hexagon work would add {f(k['launch_calibration_offset_db'], 2)} dB to every gain here, and nothing to any loss)."
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
    # ------------------------------------------------------------------ known HFSS gap
    h = d["hfss_reference"]
    a("## What the model got wrong in the one HFSS run")
    a("")
    a(
        f"Source: {h['source']}. The prediction is the continuous Phase 26 design; the HFSS run is one board at one frequency."
    )
    a("")
    L += table(
        ["quantity", "HFSS", "model"],
        [
            ("peak RHCP gain [dBi]", f(h["gain_hfss_dbi"], 2), f(h["gain_predicted_bare_dbi"], 2)),
            ("RHCP sidelobe level [dB]", f(h["sll_hfss_db"], 1), f(h["sll_predicted_db"], 1)),
            ("XPD at peak [dB]", f(h["xpd_hfss_db"], 1), f(h["xpd_predicted_db"], 1)),
            ("AR at peak [dB]", f(h["ar_hfss_db"], 2), f(h["ar_predicted_db"], 2)),
            ("LHCP share of the radiated power", f(h["lhcp_share_hfss"], 3), f(h["lhcp_share_model"], 3)),
        ],
    )
    a("")
    a(
        "The model sidelobe levels below are therefore optimistic by about the difference above, and the HFSS LHCP share includes the pin's direct radiation, which the model does not radiate. "
        "The bit-depth comparison uses the model's relative numbers; absolute sidelobe levels and cross-polar levels are not predictions."
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
        "- Modulation depth: `delta/|Y_down|` (tensor) or half of dB over `|Y_down|` (scalar), `Y_down = Ybar + Y_slab,TM(beta)`, the Phase 26 `required_map` definition. That map asks for far larger values because it targets a much higher conversion of the launched wave; the continuous Phase 26 design itself converts only the eta_conv of the scan table below."
    )
    a("")
    a(
        f"**Method.** The supercell table's leakage constants are capped at `alpha_cap` and the single-hexagon CP design is re-marched; `alpha_cap` is the value at which the gain falls {f(k['gain_loss_tolerance_db'], 1)} dB below the uncapped design "
        f"(table maximum {f(c['table_alpha_max_np_per_m'], 3)} Np/m, minimum {f(c['table_alpha_min_np_per_m'], 3)} Np/m). It is converted to physical units with the Phase 26 closed forms at the mean `beta/k0` of the design's cells; "
        "the scalar case is cross-checked with the exact truncated-Floquet complex beta. This row is the contrast-only requirement: it ignores the quantisation loss, which is added in the combined table below."
    )
    a("")
    rows = []
    for st in c["steerings"]:
        t1, s1, s2 = (
            st["contrast"]["bits_1"]["tensor"],
            st["contrast"]["bits_1"]["scalar"],
            st["contrast"]["bits_2"]["scalar"],
        )
        fr = D.outside(st["theta0_deg"], st["phi0_deg"])
        rows.append(
            (
                f"{f(st['theta0_deg'], 0)} / {f(st['phi0_deg'], 0)}",
                f(fr, 2) + " " + D.status(fr) + (" (closed form indicative only)" if st["closed_form_note"] else ""),
                f(st["alpha_cap_for_1db_np_per_m"], 3),
                f(st["check_at_cap_star"]["loss_vs_uncapped_db"], 3),
                f(st["contrast"]["bits_1"]["beta_over_k0_mean"], 4),
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
            "march points outside the table",
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
    for st in c["steerings"]:
        fc = st["floquet_check"]
        rows.append(
            (
                f"{f(st['theta0_deg'], 0)} / {f(st['phi0_deg'], 0)}",
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
    a(
        "A ratio away from one means the closed form under- or overstates the leakage for that swing (finite modulation depth); the contrast numbers of that row carry that uncertainty."
    )
    a("")
    a(
        "**Combined requirement: 2-bit coded surface.** The quantisation loss adds to the contrast loss. Each capped design is also quantised to four slot axes; the table gives the smallest `alpha_cap` whose total loss against the uncapped continuous design stays within the tolerance, converted to the tensor and scalar contrast (four-level swing). The 2-bit loss of the uncapped design is the first column. The cap is interpolated linearly in loss between the grid points of the cap sweep, so it can be slightly conservative (a larger cap than strictly needed)."
    )
    a("")
    rows = []
    for st in c["steerings"]:
        fr = D.outside(st["theta0_deg"], st["phi0_deg"])
        fc = st["floquet_check"]
        caveat = D.status(fr)
        if st["closed_form_note"]:
            caveat += "; closed form indicative only"
        if fc["ratio_exact_over_first_order"] is not None and abs(fc["ratio_exact_over_first_order"] - 1.0) > 0.1:
            caveat += "; Floquet exact/closed ratio " + f(fc["ratio_exact_over_first_order"], 3)
        for ce in st["combined_2bit"]:
            if ce["alpha_cap_np_per_m"] is None:
                rows.append(
                    (
                        f"{f(st['theta0_deg'], 0)} / {f(st['phi0_deg'], 0)}",
                        f(st["uncapped"]["q2_fixed_source_loss_db"], 2),
                        f(ce["tolerance_db"], 1),
                        "not attainable: the loss at the table maximum is " + f(ce["loss_at_table_max_db"], 2) + " dB",
                        "-",
                        "-",
                        "-",
                        caveat,
                    )
                )
            else:
                ct = ce["contrast"]
                rows.append(
                    (
                        f"{f(st['theta0_deg'], 0)} / {f(st['phi0_deg'], 0)}",
                        f(st["uncapped"]["q2_fixed_source_loss_db"], 2),
                        f(ce["tolerance_db"], 1),
                        f(ce["alpha_cap_np_per_m"], 3),
                        f(ct["tensor"]["delta_ms"], 3),
                        f(ct["tensor"]["dx_aniso_ohm"], 1),
                        f(ct["scalar"]["dx_switch_ohm"], 1),
                        caveat,
                    )
                )
    L += table(
        [
            "beam theta / phi [deg]",
            "2-bit loss, uncapped [dB]",
            "total loss tolerance [dB]",
            "alpha_cap [Np/m]",
            "tensor delta [mS]",
            "tensor dX_aniso [ohm]",
            "scalar four-level dX_switch [ohm]",
            "caveat",
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
        "**Table range caveat.** The supercell table covers a limited range of modulation wavenumbers; rows with `TABLE_EXTRAPOLATED` above (and in the scan table below) have a large share of march points outside it, and their contrast and scan numbers are extrapolations."
    )
    a("")
    # ------------------------------------------------------------------ design point
    dpc = D.dpcase
    a(
        f"## Design point: {f(dp['aperture_r_mm'], 0)} mm hexagon, beam theta {f(dp['theta0_deg'], 0)} / phi {f(dp['phi0_deg'], 0)} degrees"
    )
    a("")
    a(
        f"This beam has {f(dpc['table_extrapolated_fraction'], 2)} of its march points outside the supercell table: `{D.status(dpc['table_extrapolated_fraction'])}`; the design-point numbers share that caveat."
    )
    a("")
    a(
        "Quantisation loss is relative to the continuous design of the same family and pitch. Two definitions are given because they differ for a phase-only quantisation:"
    )
    a("")
    a(
        "- *fixed source*: the same cell currents, so the peak field scales as the pattern's directivity times its radiated power. A phase-only quantisation (`cp`) leaves every current magnitude unchanged, and the quantisation error is radiated into the non-radiating part of the spectrum; that power does not contribute to the beam (where it goes is not modelled), and this is the loss used as the headline for `cp`."
    )
    a(
        "- *renormalised*: the directivity of the quantised pattern alone, i.e. its own radiated power renormalised to the launched power. It hides the lost power. It is the energy-consistent one for an amplitude quantisation (`lp`, `lpx`), where the currents change magnitude and would change the leakage, and is the headline for those."
    )
    a("")
    a(
        "`cp` = rotating-slot tensor, RHCP (the Phase 26 family); `lp` = scalar surface, current along the radial direction; `lpx` = linear polarisation with IDEAL polarisation control (current along x in every cell), an upper bound no cell has been shown to realise. AR and XPD beyond the bound are shown as bounds, never as numbers. A negative renormalised loss is a quantised pattern whose directivity is higher than the continuous one on this sampling; a negative fixed-source loss for `lp` and `lpx` is not a gain: the plus and minus one levels have larger magnitude than the continuous modulation, so the radiated power rises (radiated power change column)."
    )
    a("")
    rows = []
    for p in D.pitches():
        cs = D.at_dp(p)
        for v in CPV + LPV + LPXV:
            r = cs["variants"][v]
            xp = bound(r["xpd_at_peak"]) if "xpd_at_peak" in r else "-"
            ar = bound(r["ar_at_peak"]) if "ar_at_peak" in r else "-"
            rows.append(
                (
                    f(p, 1),
                    v,
                    f(r["gain_headline_dbi"], 2),
                    f(r["quantisation_loss_db"], 2),
                    f(r["quantisation_loss_fixed_source_db"], 2),
                    f(r["quantisation_loss_renormalised_db"], 2),
                    f(r["radiated_power_change_db"], 2),
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
            "gain, headline [dBi]",
            "loss, headline [dB]",
            "loss, fixed source [dB]",
            "loss, renormalised [dB]",
            "radiated power change [dB]",
            "SLL [dB]",
            "pointing error [deg]",
            "XPD [dB]",
            "AR [dB]",
        ],
        rows,
    )
    a("")
    a(
        f"Peak search: {f(k['peak_refine_step_deg'], 2)} degree direction-cosine refinement around the coarse-grid peak; the sidelobe level uses the strongest coarse-grid sidelobe refined the same way. The headline gain is the continuous design's gain (aperture directivity plus ten times the base-ten logarithm of launch fraction x eta_conv) minus the headline loss; eta_conv of this design is given in the scan table."
    )
    a("")
    # ------------------------------------------------------------------ LP finding
    lp, lpx = dpc["variants"]["lp_cont"], dpc["variants"]["lpx_cont"]
    a("## Finding: a scalar surface on a monopole cannot form a linear-polarised pencil beam")
    a("")
    a(
        f"With the induced current along the surface wave's own (radial) direction the continuous scalar design has a pointing error of {f(lp['pointing_error_deg'], 2)} degrees and a strongest sidelobe at {f(lp['sll_db'], 1)} dB: "
        f"the aperture field is radially polarised, so the contributions cancel in the beam direction (a ring or split beam, not a pencil). Gain {f(lp['gain_dbi'], 2)} dBi against {f(lpx['gain_dbi'], 2)} dBi for the same hologram with ideal polarisation control (`lpx`) "
        f"and {f(dpc['variants']['cp_cont']['gain_dbi'], 2)} dBi for the CP tensor design. The `lp` rows above are therefore NOT a candidate; the bit-depth conclusion for LP rests on `lpx`, which needs a cell that sets the current direction (not designed)."
    )
    a("")
    # ------------------------------------------------------------------ scan
    az = s["scan_azimuth_deg"]
    a(f"## Steering: pitch {f(dp['pitch_mm'], 1)} mm, {f(dp['aperture_r_mm'], 0)} mm hexagon, CP tensor")
    a("")
    a(
        f"Headline gain by beam polar angle (azimuth {f(az, 0)} degrees; broadside has no azimuth). Scan loss is relative to broadside of the same variant; a negative value means the gain is higher than at broadside (the conversion efficiency rises with the steering angle in this model). "
        f"`outside` is the share of the design's march points whose modulation wavenumber lies outside the supercell table; the table supports the beam up to the polar angle {f(s['supported_scan_theta_deg'], 0)} degrees (no point outside), and the rows with a larger share are `TABLE_EXTRAPOLATED`."
    )
    a("")
    rows = []
    for t in sorted(
        {
            c_["theta0_deg"]
            for c_ in d["cases"]
            if c_["pitch_mm"] == dp["pitch_mm"] and c_["aperture_r_mm"] == dp["aperture_r_mm"]
        }
    ):
        ph = 0.0 if t == 0.0 else az
        cs = D.case(dp["pitch_mm"], dp["aperture_r_mm"], t, ph)
        v = cs["variants"]
        spread = [
            D.case(dp["pitch_mm"], dp["aperture_r_mm"], t, q)["variants"]["cp_2bit"]["gain_headline_dbi"]
            for q in sorted(
                {
                    c_["phi0_deg"]
                    for c_ in d["cases"]
                    if c_["theta0_deg"] == t and c_["aperture_r_mm"] == dp["aperture_r_mm"]
                }
            )
        ]
        rows.append(
            (
                f(t, 0),
                f(cs["eta_conv"], 3),
                f(cs["table_extrapolated_fraction"], 2) + " " + D.status(cs["table_extrapolated_fraction"]),
                f(v["cp_cont"]["gain_dbi"], 2),
                f(v["cp_3bit"]["gain_headline_dbi"], 2),
                f(v["cp_2bit"]["gain_headline_dbi"], 2),
                f(v["cp_1bit"]["gain_headline_dbi"], 2),
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
            "outside the table",
            "cp continuous [dBi]",
            "cp 3-bit [dBi]",
            "cp 2-bit [dBi]",
            "cp 1-bit [dBi]",
            "2-bit scan loss [dB]",
            "2-bit SLL [dB]",
            "2-bit pointing [deg]",
            "2-bit gain over the three azimuths [dBi]",
        ],
        rows,
    )
    a("")
    # ------------------------------------------------------------------ pitch and aperture
    a("## Pitch and aperture")
    a("")
    a(
        "Headline quantisation loss at the design point for the CP family, with the renormalised 2-bit value for comparison:"
    )
    a("")
    rows = []
    for p in D.pitches():
        cs = D.at_dp(p)
        v = cs["variants"]
        rows.append(
            (
                f(p, 1),
                cs["n_cells"],
                f(v["cp_cont"]["gain_dbi"], 2),
                f(v["cp_3bit"]["quantisation_loss_db"], 2),
                f(v["cp_2bit"]["quantisation_loss_db"], 2),
                f(v["cp_2bit"]["quantisation_loss_renormalised_db"], 2),
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
            "cp 3-bit loss [dB]",
            "cp 2-bit loss [dB]",
            "cp 2-bit loss, renormalised [dB]",
            "cp 1-bit loss [dB]",
            "lpx 1-bit loss, renormalised [dB]",
            "lpx 2-bit loss, renormalised [dB]",
        ],
        rows,
    )
    a("")
    a(
        f"The larger reference aperture (hexagon of circumradius {f(s['reference_aperture_r_mm'], 0)} mm), pitch {f(dp['pitch_mm'], 1)} mm:"
    )
    a("")
    rows = []
    for t in sorted(
        {
            c_["theta0_deg"]
            for c_ in d["cases"]
            if c_["aperture_r_mm"] == s["reference_aperture_r_mm"] and c_["pitch_mm"] == dp["pitch_mm"]
        }
    ):
        ph = 0.0 if t == 0.0 else dp["phi0_deg"]
        cs = D.case(dp["pitch_mm"], s["reference_aperture_r_mm"], t, ph)
        v = cs["variants"]
        rows.append(
            (
                f(cs["aperture_r_mm"], 0),
                f(t, 0),
                cs["n_cells"],
                f(cs["eta_conv"], 3),
                f(cs["table_extrapolated_fraction"], 2) + " " + D.status(cs["table_extrapolated_fraction"]),
                f(v["cp_cont"]["gain_dbi"], 2),
                f(v["cp_2bit"]["gain_headline_dbi"], 2),
                f(v["cp_1bit"]["gain_headline_dbi"], 2),
                f(v["cp_2bit"]["sll_db"], 1),
            )
        )
    L += table(
        [
            "circumradius [mm]",
            "theta0 [deg]",
            "cells",
            "eta_conv",
            "outside the table",
            "cp continuous [dBi]",
            "cp 2-bit [dBi]",
            "cp 1-bit [dBi]",
            "cp 2-bit SLL [dB]",
        ],
        rows,
    )
    a("")
    # ------------------------------------------------------------------ hardware
    a(
        f"## Control hardware ({f(dp['aperture_r_mm'], 0)} mm hexagon, beam theta {f(dp['theta0_deg'], 0)} / phi {f(dp['phi0_deg'], 0)})"
    )
    a("")
    a(
        f"Counting rules (one method for every row): diodes per cell = bits (one diode per bit, an ASSUMPTION; the topology is item 3); bit = 1 is the forward-biased diode at {f(k['diode_forward_current_ma'], 0)} mA, the share of forward bits is counted from the actual coded states; "
        f"DC power = (forward diodes) x I_F x V_rail, V_rail = {f(k['supply_rail_v'], 1)} V (the supply power, series resistor included; the diode alone dissipates {f(k['diode_forward_voltage_v'], 2)} V x I_F). The same formula is evaluated at {f(k['diode_forward_current_low_ma'], 0)} mA. "
        "Independent bias lines = diodes (arbitrary patterns need an independent state per diode). Sharing a line among several diodes constrains the pattern: the columns give the line count only, they are NOT simulated, and no row/column or per-ring scheme is evaluated (a per-ring scheme would also not carry the azimuth-dependent slot axes of this family). "
        f"Shift registers: eight- or sixteen-output devices in chains of eight clocked at {f(k['shift_register_clock_mhz'], 0)} MHz; FPGA pins = chains plus the shared clock, latch and enable lines (generic counts, not checked against a specific FPGA board); output current per pin is the diode forward current."
    )
    a("")
    rows = []
    for p in D.pitches():
        cs = D.at_dp(p)
        for v in ("cp_1bit", "cp_2bit", "cp_3bit"):
            b = cs["bias"][v]
            s8, s16 = b["shift_register"]["outputs_per_chip_8"], b["shift_register"]["outputs_per_chip_16"]
            rows.append(
                (
                    f(p, 1),
                    v,
                    b["n_cells"],
                    b["n_diodes"],
                    b["independent_bias_lines"],
                    b["bias_line_sharing"]["diodes_per_line_2"],
                    b["bias_line_sharing"]["diodes_per_line_4"],
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
            "lines at 2 diodes per line",
            "lines at 4 diodes per line",
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
    bias_key = "cp_" + str(rec["bits"]) + "bit"
    b = D.at_dp(rec["pitch_mm"])["bias"][bias_key]
    a(
        f"Driver current per output: {f(b['driver_current_per_output_ma'], 0)} mA (a sixteen-output device with every output on draws {f(b['shift_register']['outputs_per_chip_16']['chip_current_all_on_ma'], 0)} mA, an eight-output device {f(b['shift_register']['outputs_per_chip_8']['chip_current_all_on_ma'], 0)} mA): the output stage must source or sink this per pin; whether common logic-level shift registers are specified for it is NOT checked here (datasheets to be checked in item 3), a current-driver class may be needed. The LP rows are identical in count to the CP rows of the same bit depth (`lp` bias tables are in the JSON)."
    )
    a("")
    a(
        "For `lpx` the renormalised loss is a lower bound on the penalty: it hides the power moved out of the beam, and an equal-current-power estimate is not computed here."
    )
    a("")
    a(
        "The lower-bias column is only a power figure; the diode series resistance rises at low bias (item 3 must trade it against insertion loss)."
    )
    a("")
    # ------------------------------------------------------------------ recommendation
    a("## Recommendation and trade-offs")
    a("")
    v = D.at_dp(rec["pitch_mm"])["variants"]
    b3 = D.at_dp(rec["pitch_mm"])["bias"]["cp_3bit"]
    alt_diodes = D.at_dp(rec["alternative_pitch_mm"])["bias"][bias_key]["n_diodes"]
    losses = [D.at_dp(p)["variants"]["cp_2bit"]["quantisation_loss_db"] for p in D.pitches()]
    a(
        "- **Polarisation scope: CP, the rotating-slot tensor family.** It is the family the existing forward model covers, with a single HFSS check of the continuous hexagon (one board, one run); the scalar LP surface does not form a pencil beam (finding above), and the `lpx` bound needs a cell that controls the current direction, which no design exists for."
    )
    a(
        f"- **Bit depth: {rec['bits']}-bit as the baseline, 3-bit as the fallback.** At pitch {f(rec['pitch_mm'], 1)} mm the {rec['bits']}-bit loss is {f(v['cp_2bit']['quantisation_loss_db'], 2)} dB (fixed source; the renormalised figure {f(v['cp_2bit']['quantisation_loss_renormalised_db'], 2)} dB understates it) against {f(v['cp_1bit']['quantisation_loss_db'], 2)} dB for 1-bit and {f(v['cp_3bit']['quantisation_loss_db'], 2)} dB for 3-bit. "
        f"The 2-bit loss uses almost all of the {f(k['gain_loss_tolerance_db'], 1)} dB budget, so the coded surface is within {f(k['gain_loss_tolerance_db'], 1)} dB of the continuous one only if the cell delivers essentially the full table contrast (combined table above); 3-bit costs one more diode per cell ({b3['n_diodes']} diodes against {b['n_diodes']})."
    )
    a(
        f"- **Pitch: {f(rec['pitch_mm'], 1)} mm** (the pitch of the supercell table and of the HFSS-checked hexagon). The 2-bit headline loss is flat in pitch (from {f(min(losses), 2)} to {f(max(losses), 2)} dB over the four pitches), so in this model pitch is not a loss lever and the choice is set by the diode count and by what item 3 can build; {f(rec['alternative_pitch_mm'], 1)} mm is the alternative worth testing with a real cell ({alt_diodes} diodes against {b['n_diodes']})."
    )
    a(
        f"- **Diode count at the recommendation: {b['n_diodes']} diodes ({b['n_cells']} cells x {b['diodes_per_cell']}), {b['independent_bias_lines']} independent bias lines**, {b['shift_register']['outputs_per_chip_16']['chips']} sixteen-output shift registers, {b['shift_register']['outputs_per_chip_16']['fpga_io_pins']} FPGA pins at the stated chain length; expected DC power {f(b['dc_power_expected_w'], 1)} W at {f(k['diode_forward_current_ma'], 0)} mA per forward diode ({f(b['dc_at_lower_bias']['1mA']['dc_power_expected_w'], 2)} W at {f(k['diode_forward_current_low_ma'], 0)} mA). The bias power, not the diode count, is the sharpest cost of this architecture."
    )
    st30 = next(x for x in c["steerings"] if x["theta0_deg"] == dp["theta0_deg"])
    t2 = st30["contrast"]["bits_2"]
    ce_a, ce_b = st30["combined_2bit"][1], st30["combined_2bit"][0]
    a(
        f"- **Required contrast (beam {f(st30['theta0_deg'], 0)} degrees, contrast loss alone within {f(k['gain_loss_tolerance_db'], 1)} dB):** tensor anisotropy {f(t2['tensor']['delta_ms'], 3)} mS, i.e. a principal-axis reactance difference of {f(t2['tensor']['dx_aniso_ohm'], 1)} ohm, present in every state, with the axis rotating between {2 ** rec['bits']} orientations; the equivalent scalar swing would be {f(t2['scalar']['dx_switch_ohm'], 1)} ohm. These are sheet-level numbers; whether a diode-loaded subwavelength cell reaches them is the question of item 3. With the 2-bit loss included the requirement is larger: "
        f"a total loss of {f(ce_a['tolerance_db'], 1)} dB needs a tensor anisotropy of {f(ce_a['contrast']['tensor']['delta_ms'], 3)} mS ({f(ce_a['contrast']['tensor']['dx_aniso_ohm'], 1)} ohm), and a total of {f(ce_b['tolerance_db'], 1)} dB needs {f(ce_b['contrast']['tensor']['delta_ms'], 3)} mS ({f(ce_b['contrast']['tensor']['dx_aniso_ohm'], 1)} ohm), close to the table maximum of the leakage constant."
    )
    a(
        f"- **Scan range:** the supercell table supports the beam up to the polar angle {f(s['supported_scan_theta_deg'], 0)} degrees; the design point has {f(dpc['table_extrapolated_fraction'], 2)} of its march points outside the table and the larger angles more (scan table), so the gain at those angles is `TABLE_EXTRAPOLATED`."
    )
    a(
        f"- **Reliability of the sidelobe and cross-polar numbers:** the one HFSS run has a sidelobe level of {f(h['sll_hfss_db'], 1)} dB against the model's {f(h['sll_predicted_db'], 1)} dB; treat the model's absolute sidelobe levels as optimistic."
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
        "- Bias-line sharing schemes (row/column, per-ring) and the fit to a specific FPGA board are not evaluated; only generic counts are given."
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
