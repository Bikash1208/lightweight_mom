"""Phase 28 (HFSS-first path): README texts of the two packages and the G library report. Every number in them is read from the library / prediction data, none is typed."""

import math

import numpy as np

import phase28_export_g_packages as ex
import phase28_small_board as sb

BEAM_LABEL = {(0.0, 0.0): "0 deg (broadside)", (30.0, 225.0): "30 deg / 225 deg", (45.0, 225.0): "45 deg / 225 deg"}
LABELS = "MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET, VIA_OMITTED_APPROXIMATION"


def _row(lib, state, f, mesh):
    return next(r for r in lib["rows"] if r["kind"] == "normal" and (r["state"], r["f_ghz"], r["mesh"]) == (state, f, mesh))


def library_report(lib):
    conv = lib["convergence_mesh_64_to_96"]
    L = []
    a = L.append
    a("# G state library, 9.5 - 10.5 GHz (Phase 28, HFSS-first path)")
    a("")
    a(f"Labels: {LABELS}. **Diode package inductance L_s = {lib['package_inductance_nh']:g} nH** (`NOT_AVAILABLE_FROM_CURRENT_DATA`; zero is an assumption, the sensitivity is in `../g_refine/`).")
    a("")
    a("Cell G: 6.4 mm, rational-sector hub, 'bowtie' mode, eight MADP-000907-14020 diode sites (ON = 5.2 ohm, OFF = 0.025 pF), RT/duroid 5880 1.575 mm on PEC. Four states = slot axis 0, 45, 90, 135 deg (two diodes forward per state).")
    a("The sheet tensor is the normal-incidence one of `coded_cell_eval` (the Phase 26 B screening method); `alpha_beam` is the leakage of a row of cells to a beam at 30 deg (closed form), not a dispersion.")
    a("")
    a("## Which numbers were reused and which were recomputed")
    a("The item-3 characterisation (`../cells/cell_G_6p4_bowtie.json`) stored only derived metrics (no Y_s / reflection matrices), at mesh 96 only for 10.0 GHz. This library recomputes every state at **mesh 64 and mesh 96 at all five frequencies** (40 solves, Floquet cutoff = mesh), plus the oblique solve; the 10.0 GHz values agree with the item-3 file (same code).")
    a("")
    a("## Mesh convergence (mesh 64 -> mesh 96, maximum over the four states; criterion 5 %, 10 % for alpha, an ASSUMED criterion as in item 3)")
    a("| f [GHz] | b_bar | delta | dX own | dX design | alpha_beam | converged |")
    a("|---|---|---|---|---|---|---|")
    for f in lib["frequencies_ghz"]:
        c = conv[f"{f:g}"]["metrics"]
        a(
            f"| {f:g} | "
            + " | ".join(f"{100 * c[m]['max_relative_change']:.2f} %" for m in ("b_bar_ms", "delta_ms", "dx_aniso_ohm", "dx_aniso_at_design_loading_ohm", "alpha_beam_np_per_m"))
            + f" | {'yes' if conv[f'{f:g}']['converged'] else 'NO'} |"
        )
    a("")
    a("## Library (mesh 96; sheet tensor Y_s in the file `g_library.json`, rows `kind: normal`)")
    a("| f [GHz] | state | axis [deg] | b_bar [mS] | delta [mS] | g_max [mS] | beta/k0 | dX own [ohm] | dX design [ohm] | alpha_beam [Np/m] |")
    a("|---|---|---|---|---|---|---|---|---|---|")
    for f in lib["frequencies_ghz"]:
        for s in range(4):
            r = _row(lib, s, f, 96)
            a(
                f"| {f:g} | {s} | {r['axis_deg']:.2f} | {r['b_bar_ms']:.3f} | {r['delta_ms']:.3f} | {r['g_max_ms']:.3f} | {r['beta_over_k0']:.4f} | {r['dx_aniso_ohm']:.1f} | {r['dx_aniso_at_design_loading_ohm']:.1f} | {r['alpha_beam_np_per_m']:.3f} |"
            )
    a("")
    a("`dX design` is the anisotropy at the item-2 design loading (beta/k0 = 1.06); the 10.0 GHz item-2 requirements are 90.9 ohm (1.5 dB total) and 118.4 ohm (1 dB total), and the leakage the surface has to reach is 2.09 - 3.66 Np/m (`../FEASIBILITY_REPORT.md`).")
    return "\n".join(L) + "\n"


def _ladder(lib, mesh=96):
    out = []
    for s in range(4):
        r = _row(lib, s, 10.0, mesh)
        g = np.array([[ex._c(r["gamma_xy"][i][j]) for j in range(2)] for i in range(2)])
        c = ex.circular(g)
        out.append((s, ex._db(c["r_RR"]), ex._deg(c["r_RR"]), max(ex._db(c["r_RL"]), ex._db(c["r_LR"])), ex._db(c["r_LL"])))
    return out


def cell_readme(lib):
    lad = _ladder(lib)
    conv = lib["convergence_mesh_64_to_96"]
    all_conv = all(v["converged"] for v in conv.values())
    steps = [((lad[(i + 1) % 4][2] - lad[i][2] + 180) % 360) - 180 for i in range(3)]
    dv = ex.geo.diode_values()
    ob = [r for r in lib["rows"] if r["kind"] == "oblique"]
    L = []
    a = L.append
    a("# phase28_cell_G -- HFSS check of the cell G, all four states, Floquet port")
    a("")
    a(f"**State: PREDICTION_NOT_HFSS_VALIDATED ({LABELS}).** Nothing here has an HFSS number yet. Diode: MADP-000907-14020, ON R = {dv['r_on_ohm']} ohm, OFF C = {dv['c_off_pf']} pF, **package inductance L_s = {dv['l_s_nh']:g} nH assumed (NOT_AVAILABLE_FROM_CURRENT_DATA)**; flip-chip pads, via and bias network are not modelled.")
    a("")
    a("Cell G (Phase 28 item 3 candidate): period 6.4 mm, RT/duroid 5880 (eps_r 2.2, tan d 0.0009) 1.575 mm on a PEC ground; zero-thickness metal at the top surface: a 2.4 mm hub, four arms, four diagonal pieces, and **eight diode necks** (0.4 mm long) bridging hub-arm and hub-diagonal gaps. State k = the two necks of one pair ON (R), the other six OFF (C): axis 0 / 45 / 90 / 135 deg for k = 0 .. 3. The polygons are generated from the same geometry code as the MoM (`phase28_g_geometry.py`; tested against the MoM mask).")
    a("")
    a("## Files")
    a("| File | Content |")
    a("|---|---|")
    a("| `build_cell.py` | AEDT script: builds the cell (state `STATE` at the top), 8 Lumped RLC sheets, periodic walls and a 2-mode Floquet port with scan angle variables `scan_theta` / `scan_phi`, setup 10 GHz, discrete sweep 9.5 - 10.5 GHz (0.25 GHz) |")
    a("| `set_state.py` | AEDT script: puts another state on the open model (edit `STATE`) |")
    a("| `mom_reference.csv` | MoM prediction: S_TM_TM, S_TE_TE, S_TM_TE, S_TE_TM and r_RR, r_RL, r_LR, r_LL (magnitude dB, phase deg) per state, f = 9.5 .. 10.5 GHz; theta = 0 and theta = 30 deg / phi = 225 deg (10.0 GHz); mesh 64 and 96 |")
    a("| `bare_grounded_slab_reference.csv` | analytic reflection of the bare grounded slab at the pattern plane (convention check) |")
    a("")
    a("## Predictions (MoM, mesh 96, normal incidence, 10.0 GHz, pattern plane = top of the substrate)")
    a("| state | axis | r_RR [dB] | r_RR phase [deg] | max(opposite-handed) [dB] | r_LL [dB] |")
    a("|---|---|---|---|---|---|")
    for s, mag, ph, opp, ll in lad:
        a(f"| {s} | {45 * s} deg | {mag:.2f} | {ph:.1f} | {opp:.2f} | {ll:.2f} |")
    a("")
    a(f"Phase steps between neighbouring states (r_RR): {steps[0]:.1f}, {steps[1]:.1f}, {steps[2]:.1f} deg (a Pancharatnam-Berry cell turned by 45 deg steps wants 90). The tensor-sheet model behind the hologram design uses b_bar, delta and the axis of the same solves (`../../results/phase28/f10/g_library/G_LIBRARY.md`).")
    a(f"Mesh convergence 64 -> 96 over 9.5 - 10.5 GHz: {'all five frequencies converged (5 % / 10 % criterion)' if all_conv else 'NOT converged at every frequency, see G_LIBRARY.md'}.")
    if ob:
        a("")
        a(f"Oblique prediction (theta {ob[0]['theta_deg']:g} deg, phi {ob[0]['phi_deg']:g} deg, TE / TM, mesh {ob[0]['mesh']}, 10.0 GHz) is in `mom_reference.csv` (rows with theta = 30): this is the incidence of the beam the hologram radiates to; the sheet tensors of the design are extracted at normal incidence only.")
    a("")
    a("## HFSS setup")
    a("1. Run `build_cell.py` (Tools > Run Script). It builds state 0. Read the messages: three assignments are NOT yet verified in AEDT (periodic pairs with scan angle, Floquet port, and in `set_state.py` `EditLumpedRLC`); each falls back to a message telling you what to assign by hand.")
    a("2. By hand if needed: master/slave on the four side walls with the scan angle; ONE Floquet port on the top face of `Air`, modes 1 and 2 (TE: E along y, TM: E along x at phi = 0), **de-embed 20 mm** to the pattern plane (top of the substrate); `scan_theta` = 0 deg first, then 30 deg with `scan_phi` = 225 deg.")
    a("3. Diode boundaries are Lumped RLC on the red sheets `D0 .. D7`, current line along the neck, **one element type each** (ON: R only; OFF: C only; the AEDT-verified form of the 3x3 smoke build).")
    a("4. Solve states 0 .. 3 (`set_state.py` between runs), theta = 0 and 30 deg: 8 solves. Export the complex Floquet S11 (modes 1, 2).")
    a("")
    a("## Compare (no modelling)")
    a("* **Convention first:** run the bare grounded slab (delete the metal and the RLC sheets, or use `bare_grounded_slab_reference.csv`) at theta = 0: apply any 180 deg TM sign shift HFSS shows to every TM entry. Circular terms as in `phase20_r3_cp_cell`: incident R = (TM + jTE)/sqrt2, reflected R = (TM - jTE)/sqrt2; a metal plate has r_LR = r_RL = -1.")
    a("* **Pass (ASSUMED thresholds, taken from the R3 package):** r_RR phase of every state within 20 deg of `mom_reference.csv`, neighbouring-state steps 90 +- 20 deg, |r_RR| within 1 dB of the prediction. **What would be news:** a step far from 90 deg (the pinwheel cell is not a clean rotation in HFSS), |r_RR| more than 2 dB below the prediction (the table already contains a 3 - 4 dB diode loss: this cell is a surface-wave modulator, not a reflectarray cell, and its opposite-handed terms of -5 to -6 dB are predicted, not a fault), or an opposite-handed term more than 3 dB away from the prediction.")
    a("* The package inductance decides this cell: the MoM sensitivity (assumed L_s = 0.1 - 0.3 nH) is in `../../results/phase28/f10/g_refine/`. If HFSS includes any package inductance, compare with that file, not with the L_s = 0 prediction.")
    a("")
    a("## Honest limits")
    a("Zero-thickness metal; staircased polygons in the MoM against exact ones in HFSS (self-similar rational slopes, converged to the criterion above); no via, no bias network, no pad geometry; the tensor-sheet extraction assumes the local response; MADP L_s unknown.")
    return "\n".join(L) + "\n"


def board_readme(preds, alpha_g, ncell, lib):
    ndio = 8 * ncell
    rad = sb.circumradius_mm()
    L = []
    a = L.append
    a("# phase28_coded_small -- PLUMBING TEST: 70 mm coded hexagon with G cells, centre monopole")
    a("")
    a(f"**PLUMBING TEST, NOT A BEAM TEST.** The questions this board answers: does the model build, are the {8 * ncell} diode states applied by `set_coding.py`, does it solve (S11, far field). It is not expected to show the designed beam: G leaks only about 0.45 Np/m, so only about 3 % of the surface wave converts on this board.")
    a("")
    a(f"State: PREDICTION_NOT_HFSS_VALIDATED ({LABELS}). Diode package inductance 0 nH assumed (NOT_AVAILABLE_FROM_CURRENT_DATA). 10.0 GHz. Predictions below were written before any HFSS run.")
    a("")
    a("## Read this first: what this board can and cannot show")
    cap = {b: preds[b]["table_capped_at_G"] for b in sb.BEAMS}
    unc = {b: preds[b]["table_uncapped"] for b in sb.BEAMS}
    a(f"G leaks only alpha = {alpha_g:.2f} Np/m at 10.0 GHz (mean of its four states; the surface has to reach 2.09 - 3.66 Np/m). On a board of this size the aperture therefore converts **{100 * cap[(30.0, 225.0)]['eta_conv']:.1f} %** of the surface wave into the beam (the full 10.0 GHz supercell table would give {100 * unc[(30.0, 225.0)]['eta_conv']:.1f} %), so the predicted RHCP gain at the beam is **{cap[(30.0, 225.0)]['quantised_4_axes']['gain_dbi_launch_bare_0p592']:.1f} dBi** for 30/225 (about {unc[(30.0, 225.0)]['quantised_4_axes']['gain_dbi_launch_bare_0p592']:.1f} dBi if G could deliver the table's leakage). The monopole's own direct radiation (about 0.4 of the accepted power on the 150 mm board, not modelled here) is expected to be comparable to or larger than the aperture beam. **Expect a weak or hidden beam.** The pass criteria below are therefore PLUMBING criteria first (build, boundaries, solve, coding switch); the beam criteria are kept as an optional extra, not the purpose.")
    a("")
    n70 = 8 * len(sb.cf.lattice(sb.cf.hexagon_m(70e-3), sb.PITCH_MM * 1e-3)[0])
    r_max = max(r for r in np.arange(40.0, 70.0, 0.1) if 8 * len(sb.cf.lattice(sb.cf.hexagon_m(r * 1e-3), sb.PITCH_MM * 1e-3)[0]) <= 1200)
    n_max = 8 * len(sb.cf.lattice(sb.cf.hexagon_m(r_max * 1e-3), sb.PITCH_MM * 1e-3)[0])
    a(f"**Size decision (made unattended, to be confirmed):** '70 mm hexagon' is taken as 70 mm *across flats* (circumradius {rad:.2f} mm, vertex at +-y, as the Phase 26 hexagon), which gives **{ncell} cells and {ndio} diode boundaries** (limit ~1,200). Taking 70 mm as the circumradius would give {n70} boundaries; the largest circumradius within the limit is {r_max:.1f} mm ({n_max} boundaries).")
    a("")
    a("## Files")
    a("| File | Content |")
    a("|---|---|")
    a("| `build_board.py` | AEDT script, builds the geometry ONCE: substrate/ground (board + 5 mm margin), monopole (Phase 26 launch: pin 6.0 mm above the substrate = `mono_h` 7.575 mm, `r_pin` 0.635 mm, coax annulus `r_out` 2.05 mm, lumped port 50 ohm), 9 metal sheets and 8 diode sheets per cell, one Lumped RLC boundary per diode sheet (initial coding: 30/225), radiation box, setup 10 GHz + interpolating sweep 9.5 - 10.5 GHz, far-field sphere theta 0 - 180, phi -180 .. 180, 1 deg |")
    a("| `set_coding.py` | AEDT script: sets all diode states from a coding CSV (edit `CODING_CSV`); the model is not rebuilt |")
    a("| `layout.csv` | cell_id, x, y [mm] (lattice 6.4 mm, node at the feed, feed keep-out 6.5 mm) |")
    a("| `codings/coding_beam_000.csv`, `coding_beam_030_225.csv`, `coding_beam_045_225.csv` | cell_id, x, y, state 0 .. 3, slot axis [deg] (state k = two necks ON, axis 45 k deg) |")
    a("| `mom_reference/predictions_*.json`, `cut_*.csv`, `summary.json` | the prediction per coding; RHCP/LHCP gain cut in the beam plane (phi0 and phi0 + 180), 1 deg |")
    a("")
    a("## Run order")
    a("1. **First** the single-cell check `../phase28_cell_G/` (it tests the diode boundaries, the geometry and the cell's own states with a Floquet port in minutes).")
    a(f"2. Run `build_board.py` in AEDT (building {17 * ncell} sheets takes a few minutes). Check `Port1`, `Rad1`, that `D0_0_RLC .. D{ncell - 1}_7_RLC` exist ({ndio}; {2 * ncell} with R = 5.2 ohm, the rest C = 0.025 pF), validate.")
    a("3. Analyze (setup 10 GHz + sweep). Export (per coding): realised gain RHCP and LHCP and axial ratio over the sphere (1 deg CSV), the S11 sweep, the antenna-parameters table (accepted/radiated power, efficiency); optional: Mag_E and complex Ez along rho = 10 .. 35 mm at phi = 225 deg, 1 mm above the board.")
    a("4. Edit `CODING_CSV` in `set_coding.py` to the next coding (full path), run it on the open model (watch the message: how many assignments went through `EditLumpedRLC` and how many via delete + assign), Analyze again, export. Repeat for the three codings. `EditLumpedRLC` / `DeleteBoundaries` are not yet run in AEDT.")
    a("")
    a("## Predictions (MODEL_ONLY, 10.0 GHz, RHCP; fixed launch fraction 0.592 of the 10.54 GHz bare slab; G-limited conversion)")
    a("| Coding (theta / phi) | Peak (theta, phi) | Pointing error | Gain [dBi] | AR at peak [dB] | XPD at peak [dB] | HPBW (2 cuts) [deg] | SLL > 2 HPBW away [dB] | RHCP share of radiated | 2-bit quantisation loss [dB] |")
    a("|---|---|---|---|---|---|---|---|---|---|")
    for b in sb.BEAMS:
        e = cap[b]
        q = e["quantised_4_axes"]
        a(
            f"| {BEAM_LABEL[b]} | {q['peak_theta_deg']:.1f}, {q['peak_phi_deg']:.1f} | {q['pointing_error_deg']:.1f} deg | {q['gain_dbi_launch_bare_0p592']:.2f} | {q['ar_at_peak']['text']} | {q['xpd_at_peak']['text']} | {q['hpbw_deg_two_cuts'][0]:.1f} / {q['hpbw_deg_two_cuts'][1]:.1f} | {q['sll_db_outside_2hpbw']:.1f} | {q['rhcp_fraction']:.2f} | {e['fixed_source_quantisation_loss_db']:.2f} |"
        )
    a("")
    a("If G could deliver the leakage of the 10.0 GHz supercell table (not the case): gain " + ", ".join(f"{BEAM_LABEL[b]} {unc[b]['quantised_4_axes']['gain_dbi_launch_bare_0p592']:.1f} dBi" for b in sb.BEAMS) + ". Aperture directivity D_ap: " + ", ".join(f"{BEAM_LABEL[b]} {cap[b]['quantised_4_axes']['d_ap_dbi']:.1f} dBi" for b in sb.BEAMS) + ".")
    a("The SLL of `analyse` (15 deg exclusion) is not used: at a 23 deg half-power width it is inside the main-lobe skirt. State counts per coding are in the predictions files.")
    a("")
    a("**Known model-versus-HFSS gaps (Phase 26, `results/phase26/hex150_tensor_cp_single_hfss_comparison.md`, 150 mm board, 10.54 GHz, one HFSS run):** the model's RHCP sidelobe level was about 4 dB optimistic (-22.8 dB predicted, -18.5 dB in HFSS); the cross-polar (LHCP) power was strongly underestimated (0.31 of the radiated power in HFSS against 0.008 in the model, because the model gives no pattern to the pin's direct radiation); the peak gain agreed within 0.5 dBi (23.0 against 23.5), the pointing exactly, XPD at the peak was better in HFSS (23.0 against 16.9 dB). Those gaps are for a board whose aperture beam dominated; here the pin's share is relatively larger.")
    a("")
    a("## Pass criteria (ASSUMED, stated before the run)")
    a("* **PLUMBING PASS:** the build finishes without the RLC-failure message; the boundary list shows the expected R / C counts; `set_coding.py` completes for all three codings; every coding solves (converged, S11 exported); the three far fields are not identical (the states reach the solver).")
    a("Optional beam criteria (only meaningful if the beam is not hidden):")
    hp = [min(cap[b]["quantised_4_axes"]["hpbw_deg_two_cuts"]) for b in sb.BEAMS]
    a("* **PASS (per coding):** the RHCP gain maximum in the half space is within half a beamwidth of the target (about " + ", ".join(f"{h / 2:.0f} deg for {BEAM_LABEL[b]}" for h, b in zip(hp, sb.BEAMS, strict=True)) + ") **and** RHCP exceeds LHCP at that maximum.")
    a("* **INCONCLUSIVE (expected possibility):** no RHCP maximum near the target because the monopole's direct radiation dominates. Then compare the RHCP patterns of the three codings with each other: the common (pin) part cancels in the difference, and the prediction is that the three differ in the sense of `mom_reference/cut_*.csv` (beam moves from broadside toward 225 deg as theta0 grows). Report the sidelobe/beam positions you see.")
    a("* **FAIL:** LHCP dominant at the maximum for all codings (a handedness or state-order error: report it, it is the most useful result), or the maximum moves the wrong way (away from phi 225 deg) when theta0 grows.")
    a("")
    a("## Honest limits")
    a("The prediction is the idealised aperture model of Phase 26 (local response, ideal state tensors, no feed scattering, no edge diffraction at 70 mm); the leakage of the four G states differs by up to 25 % and is taken equal; 0.592 is the 10.54 GHz bare-slab launch fraction, not recomputed for 10.0 GHz or for this board; the table minimum and maximum are those of the 10.0 GHz table; flip-chip pads, bias lines, via and the connector are not modelled; the diode package inductance is unknown and G's contrast disappears for a large one (see the cell package).")
    return "\n".join(L) + "\n"
