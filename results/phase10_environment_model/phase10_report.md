# Phase 10 report: environment-aware electromagnetic response model

**Final status: `SURROGATE_INVALIDATED`**

Geometry inversion: **LOCKED**  
Final MTS STEP generation: **LOCKED** (Sec. 43-44 -- unaffected by this phase's findings)

## Q1. How strongly does each cardinal neighbor influence the center cell?

- C3/weak: N: 4.381e-10, E: 1.111e-09, S: 6.108e-10, W: 3.762e-10
- C3/moderate: N: 1.217e-09, E: 1.267e-09, S: 1.217e-09, W: 1.267e-09
- C3/strong_low: N: 3.746e-10, E: 3.961e-10, S: 3.746e-10, W: 3.961e-10
- C3/strong_high: N: 3.944e-09, E: 4.118e-09, S: 3.944e-09, W: 4.118e-09

## Q2. Is the influence directional?

Yes -- individual cardinal (and diagonal) contributions are reported per direction (N/E/S/W/NE/NW/SE/SW) throughout `neighborhood_study/five_cell_cases.csv` and `diagonal_cases.csv`; no rotational symmetry is assumed or imposed, and the measured values are not equal across directions for anisotropic centers.

## Q3. Is the influence approximately additive?

- C1/weak/N+S: nonlinearity error / |S0| = 0.000e+00
- C1/weak/E+W: nonlinearity error / |S0| = 0.000e+00
- C1/weak/N+E+S+W: nonlinearity error / |S0| = 0.000e+00
- C1/moderate/N+S: nonlinearity error / |S0| = 1.797e-03
- C1/moderate/E+W: nonlinearity error / |S0| = 1.952e-03
- C1/moderate/N+E+S+W: nonlinearity error / |S0| = 1.689e-03
- C1/strong_low/N+S: nonlinearity error / |S0| = 0.000e+00
- C1/strong_low/E+W: nonlinearity error / |S0| = 0.000e+00
- C1/strong_low/N+E+S+W: nonlinearity error / |S0| = 0.000e+00
- C1/strong_high/N+S: nonlinearity error / |S0| = 1.679e-01
- C1/strong_high/E+W: nonlinearity error / |S0| = 2.084e+00
- C1/strong_high/N+E+S+W: nonlinearity error / |S0| = 1.512e+00
- C2/weak/N+S: nonlinearity error / |S0| = 4.963e-04
- C2/weak/E+W: nonlinearity error / |S0| = 9.706e-05
- C2/weak/N+E+S+W: nonlinearity error / |S0| = 4.250e-04
- C2/moderate/N+S: nonlinearity error / |S0| = 6.966e-04
- C2/moderate/E+W: nonlinearity error / |S0| = 1.567e-04
- C2/moderate/N+E+S+W: nonlinearity error / |S0| = 4.974e-04
- C2/strong_low/N+S: nonlinearity error / |S0| = 2.017e-03
- C2/strong_low/E+W: nonlinearity error / |S0| = 1.844e-03
- C2/strong_low/N+E+S+W: nonlinearity error / |S0| = 1.981e-03
- C2/strong_high/N+S: nonlinearity error / |S0| = 1.572e-01
- C2/strong_high/E+W: nonlinearity error / |S0| = 2.199e+00
- C2/strong_high/N+E+S+W: nonlinearity error / |S0| = 1.753e+00
- C3/weak/N+S: nonlinearity error / |S0| = 4.371e-03
- C3/weak/E+W: nonlinearity error / |S0| = 8.551e-03
- C3/weak/N+E+S+W: nonlinearity error / |S0| = 6.564e-03
- C3/moderate/N+S: nonlinearity error / |S0| = 1.829e-02
- C3/moderate/E+W: nonlinearity error / |S0| = 2.097e-02
- C3/moderate/N+E+S+W: nonlinearity error / |S0| = 1.701e-02
- C3/strong_low/N+S: nonlinearity error / |S0| = 2.017e-03
- C3/strong_low/E+W: nonlinearity error / |S0| = 1.844e-03
- C3/strong_low/N+E+S+W: nonlinearity error / |S0| = 1.981e-03
- C3/strong_high/N+S: nonlinearity error / |S0| = 1.572e-01
- C3/strong_high/E+W: nonlinearity error / |S0| = 2.199e+00
- C3/strong_high/N+E+S+W: nonlinearity error / |S0| = 1.753e+00
- C4/weak/N+S: nonlinearity error / |S0| = 2.512e-03
- C4/weak/E+W: nonlinearity error / |S0| = 1.098e-03
- C4/weak/N+E+S+W: nonlinearity error / |S0| = 2.793e-03
- C4/moderate/N+S: nonlinearity error / |S0| = 8.814e-02
- C4/moderate/E+W: nonlinearity error / |S0| = 3.854e+00
- C4/moderate/N+E+S+W: nonlinearity error / |S0| = 4.204e+00
- C4/strong_low/N+S: nonlinearity error / |S0| = 2.184e-02
- C4/strong_low/E+W: nonlinearity error / |S0| = 2.726e-02
- C4/strong_low/N+E+S+W: nonlinearity error / |S0| = 2.649e-02
- C4/strong_high/N+S: nonlinearity error / |S0| = 1.443e-01
- C4/strong_high/E+W: nonlinearity error / |S0| = 3.167e+00
- C4/strong_high/N+E+S+W: nonlinearity error / |S0| = 3.272e+00
- C5/weak/N+S: nonlinearity error / |S0| = 0.000e+00
- C5/weak/E+W: nonlinearity error / |S0| = 0.000e+00
- C5/weak/N+E+S+W: nonlinearity error / |S0| = 0.000e+00
- C5/moderate/N+S: nonlinearity error / |S0| = 5.283e+00
- C5/moderate/E+W: nonlinearity error / |S0| = 1.900e-01
- C5/moderate/N+E+S+W: nonlinearity error / |S0| = 2.523e+00
- C5/strong_low/N+S: nonlinearity error / |S0| = 4.918e+00
- C5/strong_low/E+W: nonlinearity error / |S0| = 1.460e-01
- C5/strong_low/N+E+S+W: nonlinearity error / |S0| = 1.702e+00
- C5/strong_high/N+S: nonlinearity error / |S0| = 0.000e+00
- C5/strong_high/E+W: nonlinearity error / |S0| = 0.000e+00
- C5/strong_high/N+E+S+W: nonlinearity error / |S0| = 0.000e+00

## Q4. Are diagonal neighbors significant?

- C1/weak: mean|diagonal|/mean|cardinal| = inf
- C1/moderate: mean|diagonal|/mean|cardinal| = 0.423
- C1/strong_low: mean|diagonal|/mean|cardinal| = inf
- C1/strong_high: mean|diagonal|/mean|cardinal| = 0.416
- C2/weak: mean|diagonal|/mean|cardinal| = 0.584
- C2/moderate: mean|diagonal|/mean|cardinal| = 0.397
- C2/strong_low: mean|diagonal|/mean|cardinal| = 0.442
- C2/strong_high: mean|diagonal|/mean|cardinal| = 0.423
- C3/weak: mean|diagonal|/mean|cardinal| = 0.601
- C3/moderate: mean|diagonal|/mean|cardinal| = 0.431
- C3/strong_low: mean|diagonal|/mean|cardinal| = 0.442
- C3/strong_high: mean|diagonal|/mean|cardinal| = 0.423
- C4/weak: mean|diagonal|/mean|cardinal| = 0.712
- C4/moderate: mean|diagonal|/mean|cardinal| = 0.493
- C4/strong_low: mean|diagonal|/mean|cardinal| = 0.511
- C4/strong_high: mean|diagonal|/mean|cardinal| = 0.481
- C5/weak: mean|diagonal|/mean|cardinal| = inf
- C5/moderate: mean|diagonal|/mean|cardinal| = 0.261
- C5/strong_low: mean|diagonal|/mean|cardinal| = 0.247
- C5/strong_high: mean|diagonal|/mean|cardinal| = inf

Diagonal contributions are consistently a substantial fraction (roughly 25-70%) of the cardinal contributions across every center and contrast level tested -- NOT negligible.

## Q5. Is a 5-cell model sufficient?

## Q6. Is a 9-cell model required?

## Q7. Is a larger neighborhood required?

- C1 ((1.4, 1.4) mm): **9_CELL_REQUIRED** -- worst-case contrast (strong_high): the diagonal-ignoring error (3.070e+00 of the baseline magnitude) exceeds the solver's self-consistency floor, so the cardinal-only (5-cell) model is not sufficient; distance-2 neighbor effects remained at or below the solver's noise floor, so a 3x3 (9-cell) neighborhood is reported as sufficient beyond the immediate ring.
- C2 ((2.0, 2.8) mm): **LARGER_NEIGHBORHOOD_REQUIRED** -- worst-case contrast (strong_high): the diagonal-ignoring error is substantial (3.381 of the baseline magnitude, >= the 0.5 substantiality criterion) AND distance-2 (two-period) neighbors still produce a measurable effect above the solver's noise floor -- a 3x3 (9-cell) neighborhood does not appear to capture the full coupling.
- C3 ((2.8, 2.8) mm): **LARGER_NEIGHBORHOOD_REQUIRED** -- worst-case contrast (strong_high): the diagonal-ignoring error is substantial (3.381 of the baseline magnitude, >= the 0.5 substantiality criterion) AND distance-2 (two-period) neighbors still produce a measurable effect above the solver's noise floor -- a 3x3 (9-cell) neighborhood does not appear to capture the full coupling.
- C4 ((3.6, 2.8) mm): **LARGER_NEIGHBORHOOD_REQUIRED** -- worst-case contrast (strong_high): the diagonal-ignoring error is substantial (4.927 of the baseline magnitude, >= the 0.5 substantiality criterion) AND distance-2 (two-period) neighbors still produce a measurable effect above the solver's noise floor -- a 3x3 (9-cell) neighborhood does not appear to capture the full coupling.
- C5 ((4.2, 4.2) mm): **LARGER_NEIGHBORHOOD_REQUIRED** -- worst-case contrast (moderate): the diagonal-ignoring error is substantial (0.903 of the baseline magnitude, >= the 0.5 substantiality criterion) AND distance-2 (two-period) neighbors still produce a measurable effect above the solver's noise floor -- a 3x3 (9-cell) neighborhood does not appear to capture the full coupling.

## Q8. How does coupling decay with distance?

- C1: distance-2 relative deviation -- N2: 0.000e+00, S2: 0.000e+00, W2: 0.000e+00, E2: 0.000e+00 (significant: False)
- C2: distance-2 relative deviation -- N2: 1.663e-02, S2: 1.663e-02, W2: 3.549e-02, E2: 3.549e-02 (significant: True)
- C3: distance-2 relative deviation -- N2: 1.663e-02, S2: 1.663e-02, W2: 3.549e-02, E2: 3.549e-02 (significant: True)
- C4: distance-2 relative deviation -- N2: 3.262e-02, S2: 3.262e-02, W2: 8.408e-02, E2: 8.408e-02 (significant: True)
- C5: distance-2 relative deviation -- N2: 3.515e-01, S2: 2.416e-01, W2: 2.145e-02, E2: 1.831e-02 (significant: True)

## Q9. What geometry range is represented?

L, W in [1.0, 4.4] mm (the existing validated 18x18, 0.2 mm sweep grid), for the center cell and all 8 neighbor positions independently (a 3x3/9-cell environment representation).

## Q10. What frequency range is represented?

A single design frequency, 12.50 GHz -- explicitly NOT a swept range (a disclosed scope limitation, not a broader claim).

## Q11. What is the validated complex-response accuracy?

- **baseline1** (blind test, n=207): mean complex error 9.338e-08, mean relative error 0.574, mean phase error 30.90 deg
- **baseline2** (blind test, n=207): mean complex error 9.248e-08, mean relative error 1.008, mean phase error 40.56 deg
- **surrogate** (blind test, n=207): mean complex error 9.317e-08, mean relative error 1.350, mean phase error 38.95 deg

## Q12. What is the worst-case error?

- **baseline1**: p95 relative error 1.199, max relative error 2.691, max complex error 6.262e-07
- **baseline2**: p95 relative error 3.309, max relative error 24.072, max complex error 6.255e-07
- **surrogate**: p95 relative error 5.043, max relative error 44.831, max complex error 6.320e-07

## Q13. Does the surrogate generalize to unseen environments?

On the blind test set (which includes center C5, held out entirely from training, Sec. 21), the surrogate's mean relative error (1.350) is NOT better than Baseline 1's isolated-cell error (0.574). The surrogate does NOT outperform the (already-invalidated) isolated-cell lookup on this blind test -- it is not shown to generalize.

## Q14. Does it reproduce holographic geometry gradients?

Smooth L gradient (7x7 grid, 25 interior cells): mean |error| = 1.366e-09 A.m (25/25 predictions in-domain).

## Q15. Does it reproduce abrupt geometry transitions?

Region A | Region B step discontinuity, by distance from the boundary (periods):

- distance 0: mean |error| = 4.843e-08 A.m
- distance 1: mean |error| = 2.627e-08 A.m
- distance 2: mean |error| = 3.778e-08 A.m

Error is markedly larger immediately at the boundary than farther from it: confirmed.

## Q16. Does inverse synthesis work?

5/7 test targets REACHED (MoM-verified within the substantiality criterion); 2 reported TARGET_UNREACHABLE.

## Q17. Does heterogeneous MoM confirm the inverse result?

Yes -- every inversion result above is independently checked by re-solving the proposed geometry in its fixed environment through the real heterogeneous-supercell MoM solver (`grid_solver.solve_grid`), never accepted on the surrogate's own say-so (Sec. 39).

## Q18. Does beam reconstruction remain valid?

MoM-reference vs. surrogate-based aperture far fields: mean complex error 3.095e-13, mean phase error 1.77 deg. Detected beam pointing: MoM (17.999999999999986,), surrogate (17.999999999999986,).

## Q19. What physical control variable is validated?

NONE. Gate 10A remains CLOSED: no HFSS or measurement reference exists anywhere in this repository (the same, unchanged finding as Phase 9A). Neither S11 nor S21 (nor the local-dipole-moment quantity used throughout Phase 10) has been checked against an independent reference.

## Q20. Is geometry inversion unlocked?

**NO -- LOCKED**. Gate status: 10A=CLOSED, 10B=OPEN, 10C=OPEN, 10D=OPEN, 10E=OPEN. All five must be OPEN AND the forward/inverse/beam outcomes favorable (Sec. 42); Gate 10A alone permanently blocks the unlock until real HFSS or measurement data exists.

## Q21. Is final STEP generation unlocked?

NO -- LOCKED (Sec. 43-44 requires geometry inversion to be unlocked first, and at least one complete end-to-end design to be independently checked against HFSS where available; neither condition is met).

## Q22. What assumptions remain unresolved?

- No HFSS or measurement data exists anywhere in this project (Gate 10A, inherited from Phase 9A).
- The center-cell response quantity used throughout (`px`, the local induced-current dipole moment) is a standard, well-defined EM observable, but its relationship to S11/S21 as an actual holographic control variable remains unestablished.
- The dataset and surrogate cover a single design frequency only (10.0-15.0 GHz sweep range is NOT represented).
- Even the 9-cell (3x3) neighborhood does not fully capture the coupling for every center cell tested (`LARGER_NEIGHBORHOOD_REQUIRED` for 4 of 5 centers); the surrogate is built on the 9-cell representation regardless, as the largest tractable option, and this is a known, disclosed source of residual error.
- Sec. 33's finite-array edge/boundary effects are explicitly OUT OF SCOPE for this infinite-periodic-supercell framework (`NOT_FEASIBLE_WITHIN_INFINITE_PERIODIC_FRAMEWORK`).
- A small numerical asymmetry (a few percent) was observed between physically equivalent (translation-related) cells in the largest (7x7) solved grids; it does not affect the qualitative conclusions but is not fully explained.

## Rule compliance

- Rule 1: no HFSS/measurement data fabricated; every unavailable comparison is reported as such, not invented.
- Rule 2: heterogeneous supercells retain full geometric identity (center/neighbor/period/substrate/frequency/polarization all recorded per sample, see `environment_dataset/*.csv`); never collapsed to an isolated-cell equivalent.
- Rule 3: reference simulations -> dataset -> train/val/blind split -> surrogate fit (on TRAIN only) -> hyperparameter selection (on VALIDATION only) -> blind test -> error analysis -> this report, in that order (`dataset.py` built before `surrogate.py` was fit; `fit_surrogate` never touches the blind split).
- Rule 4: symmetry-safe grouping (`dataset.symmetry_canonical_key`) plus an explicit geometry-based holdout (center C5 reserved entirely for blind testing).
- Rule 5: Re(S)/Im(S) fit and reported directly (cartesian); phase is reported separately as a derived, wrapped quantity, never the sole target.
- Rule 6: no amplitude normalization anywhere in the error metrics (`error_metrics.py` reports absolute complex/magnitude error alongside, never instead of, relative error).
- Rule 7: no geometry inversion beyond the first-stage, MoM-verified test harness in `inverse.py`; production geometry inversion remains locked by the gate system.

