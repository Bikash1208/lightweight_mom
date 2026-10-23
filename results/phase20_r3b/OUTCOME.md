# Phase 20 R3b outcome (closed 2026-10-23): baseline CP-2 retained

**Decision (user, 2026-10-23):** the R3b optimisation is NOT resumed. The cell stays the baseline CP-2 (octagon circumradius 4.601 mm, stubs line 1.748 mm (edge-clamped) x width 1.375 mm, period 15 mm, MADP-000907-14020), which the user's HFSS run
checked at 10.5 GHz (`docs/validation/README.md` section 8). The R4 state library (`results/phase20_r4/cp2_library_v1.npz`) therefore stays valid for the RIS cell; the operating frequency is 10.5 GHz.

## What was done
* `Cp2Cell.line_diag` / `stub_width_diag` (additive; default None = unchanged geometry), a checkpointed Nelder-Mead optimiser (`scripts/research/phase20_r3b_opt.py`; every evaluation cached in `opt_cache.jsonl`, replayable), the run-2 HFSS delta
  table (`hfss_delta_table_run2.npz`), a shared-pool multi-start runner, re-score / verification / package scripts (written and smoke-tested, **not run on a real finalist**).
* Objective: nominal conditions at 10.5 GHz (theta 0; theta 30 at phi 0 / 45 / 90), tolerance corners (C_T +20 %, C_T -10 % ASSUMED, L_s +0.1 nH ASSUMED) at (0, 0) and (30, 45), band edges 10.25 / 10.75 GHz and a lower-Q term, with the
  baseline mesh correction and the HFSS delta (phi = 0 only) added.
* 249 evaluations: baseline start 95, diag_long 86, diag_longer_wide 68 (of maxfev 120); best J 26.5 / 49.6 / 74.4; the starting cell scored 27.4.

## Why it stopped (user's reading, consistent with the cache)
* 21.6 of the best J = 26.5 came from the tolerance corners, mostly the theta-30 / phi-45 corners, which are the MoM's least trustworthy condition (staircased 45 deg geometry, no HFSS calibration at phi = 45) and did not respond to the five parameters.
* Splitting the diagonal stubs from the axis stubs made the objective worse (the four states are the same cell rotated by 45 deg only when axis and diagonal stubs are equal).
* No finalist was re-scored and no new HFSS package was generated; `jobs/r3b-continue.md` is withdrawn.

## Open items carried forward
* Steps 98/92/102/68 deg and theta-30 opposite-handed -11.9 dB in HFSS stay as the cell's known properties; R5 must design codings with the **actual** state phases (not an ideal 0/90/180/270 ladder) and treat the opposite-handed term as a specular-lobe contamination.
* Diode-tolerance HFSS cases are being run by the user on the existing project (C_off = 0.0225 / 0.025 / 0.030 pF, states 1-4, theta 0 / 30). MoM comparison to follow when the files arrive (the R4 corner table has C_T +-20 % at 9.5 / 10 / 10.5 GHz).
* The mesh-correction residual of the library (12 deg near 10.3-10.8 GHz) and phi != 0 behaviour (no HFSS data) remain library limitations.
