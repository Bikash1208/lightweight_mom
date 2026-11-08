# sub6-board-model round 2 (same branch compute/sub6-board-model, MODEL_ONLY, PENDING_HFSS_CELL_CHECK)
1. Polarisation-aware synthesis: use each cell's full tensor including A_perp in the far field, and choose psi of the static AND switched cells jointly so beam A is RHCP in BOTH states (correct co-pol, as joint_opt2 in sidequest_cp_hybrid / SIDEQUEST_CP_STEP3).
2. ON-state SLL: amplitude taper (static alpha profile) and keep switched cells out of the first 2 rings around the feed and the outer ring; report SLL_on.
3. Sweep N_switched 180/200/220/240/260, sides 200 and 220 mm (220 fits the 457 mm panel). Also report eta_off/eta_on.
4. Tolerance at the best point: diode C 0.12/0.20 pF (switched-cell OFF leakage from results/sub6_slotdisc_ls_tol) -> B suppression, AR.
Same gates. Commit, push github, RESULT format, STOP.
