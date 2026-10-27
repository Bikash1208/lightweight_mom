# phase28_loop_R4 -- HFSS check of the loop cell R4 with the MEASURED MA4AGP907 diode, ON and OFF, Floquet/periodic

**State: PREDICTION_NOT_HFSS_VALIDATED (MODEL_ONLY, MEASURED_VENDOR_SPARAM diode).** Nothing here has an HFSS number yet. Candidate cell, not a design (`results/phase28/f10/ring_cell/`, `loop_R4/`).

**HFSS result (user, 2026-10-27, `hfss_results/phase28_loop_R4/NOTES.md`): the MoM prediction below is HFSS_CHECKED at ROT 0, normal incidence, 9 - 11 GHz** (all sheet-admittance components within 4.3 %, Delta B 4.22 mS against 4.05 mS predicted; the project's ON/OFF labels were swapped). The first HFSS build used manual boundaries because `build_cell.py` failed at AssignMaster; the script is fixed (one master/slave pair per face).

Cell: period 9.0 mm; RT/duroid 5880 (eps_r 2.2, tan d 0.0009) 3.175 mm on PEC; zero-thickness metal: a closed square loop, outer side 6.6 mm, strip width 1.2 mm, ONE diode (0.6 mm long, 1.2 mm wide, current along x) in the middle of the top side.
**Diode (choice stated):** MA4AGP907 = MADP-000907-14020, vendor-measured S-parameters (`results/phase28/diodes/ma4agp907_series_z.json`; series Z = Z0 (S11+S22)/S21). In HFSS: **ON = R only, 2.4 ohm (10 mA; 3.3 ohm at 5 mA, 2.0 ohm at 15 mA), OFF = C only, 28.2 fF (0 V)**, one single-type lumped boundary each (the AEDT-verified form). The measured ON reactance (-3.7 ohm at 10 GHz) is not modelled; the MoM effect of dropping it is quoted below. No series inductance is visible in the measurement (`LS_NH = 0`); `LS_NH = 0.2` is an UNTESTED sensitivity option with no prediction here. Flip-chip pads, via and bias network are not modelled.

## Files
| File | Content |
|---|---|
| `build_cell.py` | AEDT script; variables at the top: `STATE` (ON / OFF), `ROT_DEG` (0, 20), `LS_NH` (0), `R_ON`, `C_OFF`. Builds the cell, the diode sheet(s), periodic walls with scan-angle variables, a 2-mode Floquet port, setup 10 GHz, discrete sweep 9.5 - 10.5 GHz |
| `set_state.py` | AEDT script: puts ON or OFF on the open model |
| `mom_reference.csv` | MoM prediction (tag `HFSS_R10mA_Conly` = exactly what this package builds; `HFSS_R5mA_Conly`; `FULLZ_10mA` = full measured Z(f); `TENSOR_ROTATED_20deg_...` = the unrotated tensor turned by 20 deg): S_xx, S_yy, S_xy and the response of linear polarisation at 0 / 20 / 40 deg (`co_psiN`, `cross_psiN`), dB and deg |
| `predictions_summary.json`, `measured_cell_report.json` | b_bar, delta B, anisotropy, leakage alpha, band, tolerance |
| `bare_grounded_slab_reference.csv` | analytic reflection of the bare grounded slab at the pattern plane |

## Predictions first (MoM, 10.0 GHz)
* mean susceptance b_bar (xx) = 4.27 mS; switching contrast Delta B along x = 4.05 mS (R only / C only, 10 mA).
* anisotropy |b_xx - b_yy|: ON 0.01 mS, OFF 4.01 mS; Delta B along y = 0.03 mS (the switching acts along the diode axis only).
* response vs in-plane polarisation angle psi (Delta B along psi = e^T (Y_on - Y_off) e): psi 0 deg: 4.05 mS, psi 20 deg: 3.58 mS, psi 40 deg: 2.39 mS (about cos^2 psi).
* full measured Z(f) instead of R / C only (10 mA): Delta B = 3.95 mS, b_bar = 4.22 mS (R only / C only changes Delta B by +2.6 %).
* leakage alpha (exact sheet model, 6-cell period): 2.62 Np/m with the full Z(f); 2.69 Np/m with R only / C only (10 mA).

**Rotation.** A MoM solve of the cell turned by 20 deg is NOT available: the rotated-load path of the MoM gave ON = OFF (the rotated load is inactive), so it is not used. The axis-aligned 90 deg rotation reproduces the rotated tensor to 1e-14 (`loop_R4/rotation90_check.json`). The 20 deg rows of `mom_reference.csv` are the unrotated tensor turned by 20 deg (an approximation; a loop turned by 40 deg would overlap its neighbours at 9 mm). `ROT_DEG = 20` in HFSS is therefore a genuine test of the covariance.

**What HFSS checks.** A Floquet cell cannot launch a surface wave. The directionality is checked with a normal-incidence plane wave whose linear polarisation is turned by psi: from the complex S_xx, S_yy, S_xy of each state, e^T S e (co) and e_perp^T S e (cross) at psi = 0, 20, 40 deg must follow `mom_reference.csv`. ON should be isotropic (co(psi) independent of psi); OFF should depend on psi.

## Run order
1. **Convention first:** the bare grounded slab (delete the metal and the RLC sheets, or use `bare_grounded_slab_reference.csv`) at theta = 0: apply any 180 deg TM sign shift HFSS shows to every TM entry.
2. `STATE = "ON"`, `ROT_DEG = 0`: run `build_cell.py`; read the messages (the periodic walls are four master/slave pairs, one pair per face: substrate X, air X, substrate Y, air Y; they and the Floquet port have not run successfully in AEDT yet, and every failure only prints a message naming the pair or boundary to assign by hand). Floquet port: modes 1, 2 (TE: E along y, TM: E along x at phi = 0), de-embed 25 mm to the pattern plane (top of the substrate), `scan_theta` = 0.
3. Analyze; export complex Floquet S11 (modes 1, 2). Then `set_state.py` with `STATE = "OFF"`; Analyze; export.
4. Optional: `ROT_DEG = 20` (rebuild); optional: `R_ON = "3.3ohm"` (5 mA).

## Pass criteria (ASSUMED thresholds, stated before the run)
* Delta B (from the HFSS S, via Y_s = Y_in - Y_slab) within 15 % of the prediction (the MoM mesh convergence is 0.2 %; earlier HFSS-MoM cell gaps were of order 10 deg in phase).
* OFF-state anisotropy |b_xx - b_yy| within 1.5 mS of the prediction and ON isotropic within 1 mS.
* co(psi) for OFF changes with psi as predicted (Delta B along psi falls about as cos^2 psi).
* **News:** Delta B below 3 mS or above 6 mS (the design needs about 4.3 mS), or a visible series inductance (the measurement shows none).

## Honest limits
Zero-thickness metal; the lumped diode without pads; the MoM staircases the loop on a mesh (self-similar for the 0.6 mm multiples); the sheet tensor assumes the local response of an isolated cell in a periodic lattice; the diode data are one wafer/run of the vendor's series-mounted fixture (extraction noise on the OFF real part, set to zero).
