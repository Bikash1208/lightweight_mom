# C4X SUPERCELL 1 — post-run gate audit

STATUS: STOPPED at the ON surface-wave leakage gate after all five supercell rows completed. The unattended queue initially labelled its two slow/hybrid ON tracks as gate passes. Post-run slope classification corrected that label; the queue's original status/done flag is retained as an audit record. OFF 15/30 had already completed before this correction. No further solves were run.

GATES:
- Step 0 calibration PASS: all four 2025.2 PEC-lid frequencies agree with S2's 2025.1 values within 0.007%, versus the 0.5% limit.
- Step 1 PML40 spurious loss PASS: ON across/along 0.02322/0.00929 Np/m; OFF across/along 0.01929/0.01518 Np/m, all below 0.1.
- Step 2 ON SW leakage FAIL at sampled points. The clean SW track 4→3→3 has α 1.1727, 0.8822, 0.6046 Np/m at φ_s 15/30/45, all below 1.2. Its 16-point |E| field-cosine matches are 0.998/0.999. Other clean SW tracks have still lower α. The tracks 6→6→6 and 7→7→8 give larger apparent α but have mean slopes 0.01683 and 0.01778 GHz per *unit-cell* phase degree: HYBRID under the established <0.02 threshold, so they cannot certify SW leakage.
- Step 2 OFF common-frequency gate NOT CERTIFIED: the ON SW gate fails, so no valid ON operating frequency is established. As a diagnostic only, at the hybrid ON mode-6 frequency 5.863076 GHz, the field-matched OFF track 6→6 brackets it at 5.822254/5.869128 GHz with Q 57,453/13,580 and field cosine 0.960 at φ_s 30°. Its baseline-subtracted α is negative (about −0.095 Np/m by interpolation), below extraction resolution. It is not a measured positive OFF α or a certified 12:1 contrast.

DATA: All 100 raw modes, including the 20 modes returned at each of five phases, are in `C4X_SUPERCELL1_rows.csv`. The CSV gives f, Q, final tets, both kinematic Bloch β/k0 folds, local v_g, phase-matched Q_base, α, radiated-harmonic angles, matched track and 16-point |E| cosine. The unit-cell-phase slope is eight times the supercell-phase slope in magnitude. Q_base is held on the same unrotated PML branch across each matched track.

| ON track | Class | Modes φ_s 15→30→45 | f GHz | Q | α Np/m |
|---|---|---|---|---|---|
| 4 | SW | 4→3→3 | 5.76385 / 5.68789 / 5.61019 | 86.8 / 111.1 / 154.4 | 1.1727 / 0.8822 / 0.6046 |
| 9 | SW | 9→9→11 | 6.05895 / 6.12344 / 6.18871 | 414 / 447 / 734 | 0.3190 / 0.2971 / 0.1816 |
| 11 | SW | 11→12→12 | 6.31464 / 6.25084 / 6.19549 | 202 / 235 / 213 | 0.6923 / 0.6294 / 0.7409 |
| 6 | HYBRID | 6→6→6 | 5.88283 / 5.86308 / 5.81972 | 66.0 / 99.1 / 92.4 | 6.132 / 2.501 / 1.947, **not SW leakage** |
| 7 | HYBRID | 7→7→8 | 5.93730 / 5.94675 / 6.00398 | 193 / 76.4 / 66.8 | 4.099 / 3.156 / 2.135, **not SW leakage** |

OFF resonance inventory: low-Q supercell modes appear at φ_s 15/30 near 5.519/5.503 GHz (Q 5.38/5.40) and 6.074/6.079 GHz (Q 1.64/1.63). A second φ_s 30 mode at 6.0669 GHz has Q 591.8. The low-Q modes have broad linewidths, so their spectral influence at 5.863 GHz cannot be dismissed from these eigenvalues alone. The supplied 5.61 GHz unit-cell resonance does not map to a uniquely identified low-Q supercell mode at exactly that frequency. The field-matched OFF continuation near 5.863 GHz is the high-Q track stated above.

Solver resources and wall time: 48 logical cores requested / 48 granted in each profile. ON φ_s 15/30/45: 1244/1143/1253 s, final tets 84,927/67,729/86,518. OFF φ_s 15/30: 1097/892 s, final tets 64,258/56,319. Each returned 20 modes through at least 7.00 GHz, exceeding the 6.5 GHz span requirement. One local AEDT 2025.2 solve ran at a time, with a 3600 s exact-own-PID timeout. No retry was needed.

FLAGS: The queue's preliminary gate used a minimum v_g and missed the slope classification; its `COMPLETE`/`OFF PASS` status is superseded by this post-run audit. The OFF solves were therefore run although the corrected ON gate says they should have been skipped. α for RES/HYBRID tracks is only an arithmetic diagnostic because division by their small or changing slope can inflate it. The OFF negative subtraction reflects the supercell Q exceeding the unrotated baseline Q; it is not physical negative leakage. Matching uses sampled |E| magnitudes, not full complex-field overlap. A supercell eigenmode does not establish which mode a pin feed excites. No C4X realisation cell or board solve was touched.

FILES: `C4X_SUPERCELL1_rows.csv`, `C4X_SUPERCELL1_analysis.json`, queue status/done flags, raw mode/build/convergence/profile JSON and logs. Projects remain in private `hfss_work/` and are not pushed.

HANDOFF: STOPPED_ON_SW_GATE | C4X ideal | best clean SW α 1.1727 Np/m at 5.76385 GHz | OFF contrast not certified | HFSS 2025.2 | no further solve requested.
