# BASELINE NO_VIAS — HFSS eigenmode (one page)
MODEL_ONLY / NOT_HFSS_VALIDATED. One geometry STEP; OFF/ON differ only in sheet impedance.

**1. Import STEP in mm; assign material pec to METAL_PEC (all three metal solids). STEP imports as vacuum.** SUBSTRATE: RT5880 eps_r=2.2, tanD=0.0009;  AIR: vacuum. Copper thickness 0.035 mm. no vias, FR4, pads or stubs.

2. AIR intentionally overlaps all objects. Create project variables $ph_x, $ph_y in degrees. Lattice pairs on its FOUR complete side faces ONLY: x=-5.8/+5.8, phase $ph_x; y=-5.8/+5.8, phase $ph_y. TOP z=150: Perfect E. BOTTOM z=-33.210: Radiation, **no PEC below the stubs**. Single six-face AIR box retains the calibration cavity fix.

3. Frame is the calibration's +90-degree rotation: slot along y. ACROSS: ($ph_x,$ph_y)=(phi,0); ALONG: (0,phi). p=11.6 mm. ONE Impedance boundary on DIODE_SHEET (1.00 x 0.35 mm): OFF Rs=5.424, Xs=-79.320; ON Rs=1.05, Xs=+7.65 ohm/sq. Current line: (-0.450,0,0.0175) to (-1.450,0,0.0175). Do not multiply sheet values by W/L again. Record actual boundary values; these are the 5.8-GHz values, replacing the old parallel-RLC model.

4. Eigenmode: min f=4.5 GHz, 8 modes, delta f=0.1%, max 15 passes. Initial mesh: substrate 1 mm; top metal 0.5 mm; diode/gaps/vias/anti-pads 0.1 mm; FR4 0.5 mm. Clean Up Solutions for every state/phase. Start phi=80,85,90,95,100 deg; extend until the SAME top-surface TM branch brackets 5.8 GHz. Track Ez above the patch and decay upward; record all nearby modes/Q. Reject underside/air/via modes. ON mode 1 need not be the wanted branch.

5. Read f (Hz), Q and phase for >=3 nearby points per direction/state; refine phase spacing locally if necessary. Interpolate phi_58 and Q_58. beta=phi_58*pi/(180*0.0116) rad/m; v_g=2*pi*0.0116*df/d(phi radians); alpha_x=2*pi*5.8e9/(2*Q_58*v_g) Np/m. Return results_template.csv and mode-Ez identification. This alpha_x is Q-derived attenuation, NOT the MoM holographic leakage alpha. beta/k0=beta/(2*pi*5.8e9/299792458).

6. Sanity: bare RT5880/ground, 90 deg: historical model 6.291 GHz (150-mm lid); identify the top bound mode. No-via/via comparisons require the same state/phase/branch. Predictions and their limitations are in predictions.json; via+FR4 attenuation has no 3D MoM prediction. Method source: sub6_final_cal.py and results/sub6_final_board/hfss_data.json; archived setup at commit 21fe444, hfss_exports/sub6_island_choke/HFSS_SETUP.md.

<!-- MODEL_READOUT_REFERENCE -->
**Expected no-via sheet-model readouts, mesh 96, MODEL_ONLY (along / across):**

| Quantity | OFF | ON |
|---|---:|---:|
| beta/k0 | 1.05580 / 1.14403 E | 1.05239 / 1.00759 E |
| alpha_x loss estimate, Np/m | 0.02221 / 0.22622 E | 0.26183 / 0.47777 E |

E = extrapolated; not acceptance limits. Mesh trends and all roots: model_reference.csv. ON is mesh-unreliable; determine its branch and slope from HFSS.
Check Ez is small at the 150-mm lid; otherwise flag AIR_HEIGHT_NOT_CONVERGED and verify a higher lid before interpreting the mode.
