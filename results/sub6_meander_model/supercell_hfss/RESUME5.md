# OFF twin critical-path resume (2026-11-19)

## STATUS

Geometry audit, eight-cell dry build, OFF-twin phi_s=15/30 supercell solves, four unrotated PML40 baseline rows, and two lidded calibration rows completed. The OFF gate is not passed; no driven-strip test was started.

## GATES

- Geometry: PASS. AEDT 2025.2 read back two PEC halves per cell, each 43.539602 mm2, total 87.079204 mm2 against expected 43.5+43.5 / 87.08.
- Lidded calibration: PASS. Own 2025.2 PEC 85-degree branch-1 f differs from the 2025.1 reference by +0.234% across and +0.144% along, both within 0.5%.
- PML spurious alpha: PASS. Branch 1/2 across = 0.0136/0.0458 Np/m; along = 0.0200/0.0288 Np/m, all <= 0.1.
- OFF leakage: **FAIL TO PASS** at the ON x-low 15-degree f = 5.664117 GHz. ON alpha = 1.538821 Np/m, so alpha_OFF must be <= **0.128235 Np/m**. No eligible field-matched OFF surface-wave track brackets that frequency with alpha below the limit. This is not a numerical alpha_OFF at exactly f_op, because no clean corresponding branch exists there.
- Driven strip: NOT RUN because the OFF gate did not pass.

## DATA

The old end-run formula stopped at x=5.230679 mm, while the rim at the slot's inner y edge is x=5.366563 mm. It left a 0.135884 mm copper bridge. The corrected cut extends to x=+/-5.700 mm, beyond the 5.6 mm disc. After subtract, AEDT listed one PEC multi-lump object of volume 3.047772 mm3. `separate_bodies` returned two PEC objects, 43.539602 mm2 each. The full eight-cell dry build read back eight such pairs, 16 separate copper halves, lattice X/Y, PML top, PEC bottom, RT5880 substrate, and 20 requested modes from 5.4 GHz. Design validation passed; no solve was run in the dry build.

OFF-twin phi_s=15: HFSS 2025.2 returned 20 modes, highest 6.696874 GHz; 20 adaptive passes, max 91,679 solved tets, 22m27s solution time, normal completion. Near the 5.664 GHz ON x-low point, OFF modes 2/3 are at 5.662908/5.663298 GHz with Q 8.2/3.0. Their low Q makes frequency-only matching unsafe; the phi_s=30 row and saved fields will determine which, if either, follows the ON surface-wave branch. All 20 raw f/Q values are in the per-row modes JSON.

OFF-twin phi_s=30: 20 modes through 6.718255 GHz; max 131,335 solved tets, 31m41s solution time, normal completion. The near-5.664 modes at phi_s=15 move **up** to 5.677509/5.678460 GHz at 30 degrees (Q 7.6/3.0), unlike the decreasing-frequency ON x-low track. The closest sampled field pattern to ON x-low at both phases is OFF mode 1 (field cosine 0.925/0.899), but its 5.592315 -> 5.602267 GHz track is nearly flat, Q 38.8/40.3. OFF candidates 4/5 have decreasing frequency but weaker field match and higher frequency. None is silently assigned as an OFF surface-wave continuation.

Own 2025.2 unrotated PML40 unit-cell baseline rows (same corrected open geometry) at 85/95 degrees returned branch-1 Q 1905.6/1686.8 across and 1881.9/1787.2 along. They supply the phase-matched Q baseline required for the radiative-alpha subtraction; reference 2025.1 Q values are not substituted. The upper unit-cell branches are also saved in the row JSONs.

**All OFF modes:** `off_all_modes_resume5.csv` contains every raw mode returned at phi_s=15/30 (40 rows, 32 within 5.4-6.5 GHz). Each row gives f, Q, max tets, both possible unfolded beta/k0 values and -1 harmonic angles, a field-matched partner, 16-point |E| field cosine, frequency slope, and |v_g|. The two fold values are kinematic for box/low-Q modes and are not asserted to be propagating branches. v_g is blank for box modes or modes outside the analysis window. Nonbox modes were paired one-to-one by field cosine plus frequency continuity; the pairing is a candidate, not a full complex-field overlap.

The ON x-low modes 3/3 at 15/30 degrees were compared with the OFF fields at the same phases. The strongest sampled-field match is OFF modes 1->1 (cosine 0.925/0.899), but this track moves 5.592315 -> 5.602267 GHz, has Q 38.8/40.3 and |v_g| 22.17 Mm/s: a flat low-Q mode, not a clean surface wave. OFF modes 2->2 and 3->3 bracket 5.664 GHz and are close in f at 15 degrees, but their Q is only 8.2/7.6 and 3.0/3.0, and their ON-field cosine is weaker. The plausible decreasing-frequency OFF tracks 4->5 and 5->4 occur above f_op (5.804->5.781 and 5.832->5.774 GHz), with ON-field cosine 0.611/0.685 and 0.785/0.684. **There is no unambiguous OFF surface-wave continuation of ON x-low at f_op.**

`off_gate_candidates_resume5.csv` gives every candidate's raw f/Q, 15-to-30-degree |v_g|, both phase-matched own Q baselines, ON-field cosine and diagnostic alpha. To minimize alpha under branch ambiguity, it uses the smaller branch-1/2 Q_base: at phi_s=15, 1633.5 versus branch-1 1844.8; at 30, 1543.0 versus 1874.2. If the descending tracks were treated as surface waves, the **smallest** diagnostic alpha is 1.170 Np/m (track 4->5 at 30 degrees), over nine times the 0.128 gate. The low-Q/flat tracks yield much larger diagnostic values, but their alpha formula is not interpreted as surface-wave leakage. No candidate establishes the required OFF suppression at the common ON frequency.

## FLAGS

- The prior 0 mm2 result came from reading a stale/disconnected object handle instead of enumerating the copper bodies after the edge-opening subtract. The end-run geometry also needed the rim-reaching correction above.
- The OFF sheet is an open-gap twin bound; it has no diode impedance. Its Q baseline must be obtained from own 2025.2 unrotated PML40 rows before applying the alpha gate.
- The Q baseline was obtained as described above. Lidded calibration references from commit 242f3fb were used only for the 0.5% geometry/stack check, never as Q_base.
- The 15/30-degree OFF modes include numerous Q<=60 modes. Their beta and harmonic angle columns are Bloch-fold arithmetic, not evidence of a radiating surface wave. With only two supercell phases, mode crossings and v_g remain uncertain; this uncertainty cannot turn any matched candidate into a demonstrated gate pass.
- No B2p2, anchor bracketing, beam-angle, or f_op selection is included.

## FILES

`D3p6_OPEN_area_audit_resume5.json` and its log contain the single-cell object/volume audit. `SUPER_D3p6_OPEN_ph15_PML40_resume5_dry_build.json` and status/log/validation contain the eight-cell dry-build readback. `SUPER_D3p6_OPEN_ph{15,30}_PML40_resume5_*` contain the raw supercell modes, profiles, field samples and watchdog records. `D3p6_OPEN_PML_*_resume5_*` and `D3p6_OPEN_PEC_*_resume5_*` contain own unit-cell baselines/calibration. `off_all_modes_resume5.csv`, `off_gate_candidates_resume5.csv`, and `off_resume5_summary.json` provide the complete mode inventory and gate calculations, reproduced by `scripts/research/sub6_open_twin_resume5_analyze.py`. Private AEDT projects are not pushed.

## HANDOFF

- STATUS: OFF twin complete at phi_s=15/30; gate not passed; driven strip not run.
- CELL JSON: two 43.539602-mm2 halves per cell, 87.079204 mm2 total, read back in all eight cells.
- beta/k0 and alpha: every raw mode/kinematic fold and candidate leakage in the CSVs; no clean OFF SW at 5.664117 GHz.
- M_t, loss ratio: not established by this modal test.
- Evidence: HFSS 2025.2 solved rows and fields; reference frequencies only for calibration.
- Next request: decide whether to change the cell or investigate the low-Q OFF modes before any driven-strip run.
