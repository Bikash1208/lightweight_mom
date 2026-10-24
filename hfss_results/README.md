# `hfss_results/` — External HFSS Simulation Results

Storage for simulation outputs returned from external Ansys HFSS runs. These files provide the independent full-wave ground truth used to validate Method-of-Moments predictions, establish empirical calibration deltas, and test comparison workflows in the HMSA Designer application.

Raw Touchstone files (`*.s4p`, `*.s36p`) are gitignored by repository convention; small representative files, far-field pattern CSVs, and parameter reports are tracked here.

---

## 1. HFSS Ingestion Lifecycle

```
hfss_exports/<case>/          HFSS (External Run)          hfss_results/<case>/
  AEDT script / STEP    ───>    Solve in Ansys     ───>      Gain CSVs, S-params,
  layout / README                                            antenna parameter reports
                                                                    │
                                                                    ▼
                                                            scripts/research/
                                                            src/hologram/compare
                                                            (MoM vs HFSS validation)
```

1. **Export:** A design phase or solver tool writes an export package to `hfss_exports/<case>/` containing layout geometry (STEP/DXF), AEDT Python setup scripts, and a README specifying boundary conditions, excitation ports, and frequency sweeps.
2. **Simulation:** The user or compute resource executes the model in Ansys HFSS.
3. **Ingest:** Output data (Floquet S-parameters, Touchstone matrices, realized gain tables, axial ratio sweeps, and antenna parameters) are placed into `hfss_results/<case>/`.
4. **Comparison:** Research scripts (e.g. `scripts/research/phase20_r4_hfss_check.py`, `scripts/research/phase26_fourhex_design.py`) or the web app (`src/app/routers/tensor_cp.py`, `/compare`) parse these results, align reference planes, and report numerical differences against MoM forward models.

---

## 2. Subfolder Index

| Subfolder | Phase / Source | Contents & Cited File Structure | Used By / Status |
|---|---|---|---|
| `hex150_gap0/` | Phase 23 / 26 (`docs/briefs/PHASE26_LEAKY_TENSOR_CP_HEXAGON.md`) | • `HEX_HFSSDesign1.s4p`: 4-port Touchstone matrix for 150 mm hexagon board.<br>• `port3_antenna_params_8_10p5GHz.txt`: Antenna parameter report at 8.0 GHz and 10.5 GHz (Peak Directivity, Peak Realized Gain, Radiated/Accepted Power, Front-to-Back Ratio, rE field).<br>• `port3_realized_gain_RHCP_LHCP_10p5GHz.csv`: Far-field pattern cut at 10.5 GHz.<br>Header: `"Freq [GHz]","Phi [deg]","Theta [deg]","dB(RealizedGainLHCP) []","dB(RealizedGainRHCP) []"` | Multi-port launcher calibration and leakage analysis (`scripts/research/phase26_hex150_leakage.py`, `scripts/research/phase26_fourhex_scope.py`). |
| `hex150_tensor_cp_single/` | Phase 26 / 27 (`docs/briefs/PHASE27_APP_TENSOR_CP.md`) | • `S11_8_12GHz.csv`: Monopole feed return loss sweep 8–12 GHz.<br>Header: `"Freq [GHz]","dB(S(1,1)) []"`<br>• `axial_ratio_phi224_8_12GHz.csv`: Axial ratio sweep along phi=224° beam cut.<br>Header: `"Freq [GHz]","Phi [deg]","Theta [deg]","dB(AxialRatioValue) []"`<br>• `realized_gain_RHCP_LHCP_10p5GHz.csv` & `_1deg.csv`: Realized gain cuts.<br>Header: `"Freq [GHz]","Phi [deg]","Theta [deg]","dB(RealizedGainLHCP) []","dB(RealizedGainRHCP) []"`<br>• `efield_line_magE_phi225_1mm.csv`: Near-field magnitude line cut 1 mm above board along phi=225°.<br>Header: `"Freq [GHz]","Phase [deg]","Distance [cm]","Mag_E [V_per_meter]"`<br>• `antenna_params_10p5GHz.txt`: Parameter summary at 10.5 GHz (Max U, Peak Directivity, Peak Realized Gain, rE Field at Theta 30°, Phi 225°). | Verified single-hexagon tensor-CP board; app comparison flow (`src/hologram/tensor_cp_compare.py`, `tests/test_app_hfss_compare.py`). |
| `phase20_r3_cp_cell/` | Phase 20 R3 / R4 (`docs/briefs/PHASE20_RIS_ROADMAP.md`) | • `NOTES.md`: Advisor summary of first HFSS run on CP-2 cell with MADP diodes (states 1–4, theta 0° and 30°, Floquet ports).<br>• `NOTES_tolerance.md`: Summary of diode capacitance tolerance sweep (Coff = 0.0225, 0.025, 0.030 pF).<br>• `hfss_bare_slab_theta0_30.csv`: Bare grounded substrate reference.<br>Header: `"scan_theta [deg]","state []","Freq [GHz]","re(S(FloquetPort1:1,FloquetPort1:1)) []",...`<br>• `hfss_cp2_states1-4_theta0_30.csv`: Initial 4-state run 8–12 GHz.<br>• `hfss_cp2_states1-4_theta0_30_run2_fine.csv`: Fine frequency step run.<br>• `hfss_cp2_tolerance_Coff_0225_025_030.csv`: Diode tolerance sweep.<br>Header: `"Coff [pF]","scan_theta [deg]","state []","Freq [GHz]",...`<br>• `advisor_compare.py` & `advisor_compare_output.txt`: Standalone Floquet S-parameter comparison script and run log. | Ground-truth calibration for CP-2 RIS cell at 10.5 GHz; anchors R4/R5 Delta tables (`scripts/research/phase20_r4_hfss_check.py`, `scripts/research/phase20_r5_rebuild_delta.py`, `src/ris/calibrated_library.py`). |
