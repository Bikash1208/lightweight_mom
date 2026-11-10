# STATUS (session handover 2026-11-10). MODEL_ONLY / NOT_HFSS_VALIDATED everywhere. Branch compute/sub6-jobs1-prep (see `git log -1` for the commit; pushed to github + hub).
Python: C:/lightweight_mom/.venv313/Scripts/python.exe. DO NOT switch branches while a job runs (it reverts uncommitted/branch-local fixes, e.g. ISL_GAP_CELLS=8 in diag_full_tensor_lib.py). Use `SUB6_MESH=96`, `ROT180=1` for the HFSS-consistent frame.

## DONE
- Island gap check: design table used the 8-cell gap (0.967 mm); my earlier tensor test used 5 -> superseded; corrected tensors in results/sub6_final_board/diag_full_tensor/ (tensors.json, board_full_tensor_gap8.json; *_SUPERSEDED = old).
- No-FR4 HFSS vs model: rot180 fits (power corr OFF 0.80 / ON 0.94 vs identity 0.71 / -0.02); no sign difference found in model/layout lines (see chat RESULT); cause unknown -> keep the rot180 design convention (A designed at phi 180, B at phi 0). results/jobs1/noFR4/.
- Field maps + via candidates (fieldmaps/), small-board STEP variants V0-V5 (hfss_exports/sub6_jobs1_variants/, STEPs git-ignored; copies on Desktop sub6_board_hfss/jobs1), literature (LITERATURE.md, abstract-level only).
- BAR50 cell MoM (mesh 96, gap 8): 02V C0.11/Rp4k/0.6 nH: alpha_OFF 2.01 / alpha_ON 0.226 = 9.5 dB; 02L (Ls 0.4 nH): 2.00 / 0.0874 = 13.6 dB (diode_mom_bar50.json, diode_mom_bar50_02L.json).
- Iteration-1 study v2 (k = 1 primary, k = 2 option, N 80-150, 220/280/360 mm, inner/outer/sectors/field-weighted, rot180 frame, robustness sw/static alpha x0.85/x1.15 + ON x1.2): iter2_02V.json, iter2_02L.json (+ .log). Summary:
  * 02V: targets (|G_A-G_B|<=3, supp>=10, |dA|<=1) met only marginally (supp 10.0-11.3 nominal; 8.2-9.3 if the ON amplitude is 20 % worse -> no margin). Best k=1 that survives the alpha robustness: 220 mm, N150, sector_perp, G 7.7/8.2 dBi; 280 mm k=1 N150 outer G 10.4/7.7 not robust.
  * 02L: all targets with margin: best k=1 survivors 280 mm N150 sector_perp cap 0.5 G 8.7/8.3 dBi supp 12.8 (>= 11.1 worst); k=2 280 mm N120 G 10.2/10.1 supp 11.8-13.8; 360 mm k2 N150 G 12.6/10.4.
- Minimum board for a clear beam vs the feed (jobs1_min_board.py, min_board.json): all-static G_A = 1.7 (80 mm) / 5.8 (100) / 8.5 (120) / 10.6 (140) / 12.4 (160) / 14.0 (180) / 15.4 (200) / 16.3 (220) / 19.2 (280 mm) dBi; feed+rim reference 7.6 dBi (HFSS 9-cell) -> >= +6 dB needs >= 180 mm (164 cells), +8 dB 200 mm (208 cells). A+B boards need >= 280-360 mm (G 8-12 dBi). Iteration-1 220 mm designs (G 7-8 dBi) are NOT above the feed level.

## PARTIAL / NOT DONE
- 9-cell HFSS package for the best design: NOT exported (waits for V0/V2/V5/V1 results); needs a circular-partition cell_table generator (see hfss_exports/sub6_jobs1_variants/PACKAGE_9CELL_PLAN.md). 9-cell boards cannot validate beams (user HFSS result): use a full iteration-1 board.
- Gathered elements (k >= 2): strength / ON suppression per diode are ASSUMED equal to the single cell; needs an HFSS cell check. Study gains use a simple efficiency model (launch 0.656 x converted power).
- Diode survey (E) dropped by user decision; premium diodes later.
- Cause of the 180 deg beam swap (HFSS vs model) open: check AEDT sphere CS / handedness (user doing).

## RESUME (from C:/lwmom_final, branch compute/sub6-jobs1-prep)
  ROT180=1 SUB6_MESH=96 python scripts/research/jobs1_iter2_study.py results/jobs1/iter2_02L.json 2.0 0.0874      # study, ~2 min (02V: ... 2.01 0.226)
  DIODE_SET=02L python scripts/research/jobs1_diode_mom.py 4 results/jobs1/diode_mom_bar50_02L.json             # cell MoM, ~11 min, 4 workers x 8 threads
  python scripts/research/jobs1_min_board.py results/jobs1/min_board.json                                          # board-size check
  RF_OUT_DIR=... RF_FR4_T=0.2 RF_R_ANTI=0.3 RF_VIA_MODE=base|moved|shield|none RF_OUT_NAME=X.step python scripts/research/sub6_final_rf_std.py test9   # STEP variants (+ sub6_final_rf_std_check.py <step> 9 <json>)
