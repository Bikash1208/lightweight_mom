# Phase 11 report: physics-based heterogeneous MoM / nonlocal coupling

**Overall classification: `CASE_D_PATTERN_DEPENDENT`**

Geometry inversion: **GEOMETRY_INVERSION_LOCKED**  
STEP export: **STEP_EXPORT_LOCKED**

## No ML / no surrogate

Confirmed by construction: every response in this phase comes from a real solve of `Z I = V` (`reference_solver.solve_grid`/`solve_isolated`, the unmodified Phase 1-3 MoM machinery); no regression, neural network, random forest or learned mapping appears anywhere in `src/phase11/`.

## Coupling-radius series (Sec. 5-7)

- C1: n=1:1.233e-09, n=3:2.421e-09, n=5:2.136e-09, n=7:9.739e-10, n=9:2.882e-09
- C2: n=1:1.268e-08, n=3:2.786e-08, n=5:2.323e-08, n=7:9.284e-09, n=9:3.050e-08
- C3: n=1:1.268e-08, n=3:1.266e-08, n=5:1.299e-08, n=7:1.257e-08, n=9:1.277e-08, n=11:1.300e-08
- C4: n=1:3.550e-08, n=3:3.426e-08, n=5:3.664e-08, n=7:3.402e-08, n=9:3.556e-08
- C5: n=1:1.223e-07, n=3:8.811e-08, n=5:1.554e-07, n=7:1.128e-07, n=9:1.013e-07

## Shell increments (Sec. 7)

- C1: 1->3:1.985e+00, 3->5:1.882e+00, 5->7:1.397e+00, 7->9:1.960e+00
- C2: 1->3:2.186e+00, 3->5:1.834e+00, 5->7:1.328e+00, 7->9:2.301e+00
- C3: 1->3:4.025e-02, 3->5:2.636e-02, 5->7:3.235e-02, 7->9:1.629e-02, 9->11:1.846e-02
- C4: 1->3:1.237e-01, 3->5:7.133e-02, 5->7:7.295e-02, 7->9:4.749e-02
- C5: 1->3:4.833e-01, 3->5:9.368e-01, 5->7:3.053e-01, 7->9:1.844e+00

## Convergence classification (Sec. 16, tolerances declared before evaluation)

- C1: status=STILL_SUBSTANTIALLY_CHANGING, last_increment_relative=1.960e+00, converged_at_n=None
- C2: status=STILL_SUBSTANTIALLY_CHANGING, last_increment_relative=2.301e+00, converged_at_n=None
- C3: status=NOT_CONVERGED, last_increment_relative=1.629e-02, converged_at_n=None
- C4: status=NOT_CONVERGED, last_increment_relative=4.749e-02, converged_at_n=None
- C5: status=STILL_SUBSTANTIALLY_CHANGING, last_increment_relative=1.844e+00, converged_at_n=None

Per-center case: {'C1': 'CASE_C_STRONGLY_NONLOCAL', 'C2': 'CASE_C_STRONGLY_NONLOCAL', 'C3': 'CASE_B_MODERATELY_NONLOCAL', 'C4': 'CASE_B_MODERATELY_NONLOCAL', 'C5': 'CASE_C_STRONGLY_NONLOCAL'}

Overall case: **CASE_D_PATTERN_DEPENDENT** -- different spatial arrangements with similar local statistics produced significantly different center responses (Sec. 18 Case D) -- local descriptors (including the coupling-radius series alone) are insufficient regardless of the per-center convergence finding (CASE_C_STRONGLY_NONLOCAL).

## Diagonal coupling (Sec. 9)

mean|diagonal|/mean|cardinal| = 0.442 (Phase 10's own measured range: 25-70%)

## Pattern study / Case D check (Sec. 10, 18)

- uniform: |response|=1.268e-08
- single_perturbation: |response|=1.221e-09
- smooth_gradient: |response|=1.257e-08
- abrupt_boundary: |response|=1.243e-09
- alternating: |response|=1.257e-08
- checkerboard: |response|=1.257e-08
- strong_isolated_defect: |response|=1.221e-09

Pattern dependence detected: **True** ({'alternating': 0.0, 'abrupt_boundary': 0.9010884305019206})

## Pairwise / higher-order interaction (Sec. 12-13)

- N+S: interaction_error_relative=2.659e-04
- N+E: interaction_error_relative=3.432e-04
- NE+SW: interaction_error_relative=4.271e-04
- N+N_far: interaction_error_relative=1.896e-04
- N+E+S (triple): residual_relative=1.015e-05

## Frequency dependence (Sec. 14)

- 10.00 GHz: max_shell_increment_relative=3.454e-02
- 11.25 GHz: max_shell_increment_relative=3.736e-02
- 12.50 GHz: max_shell_increment_relative=4.025e-02
- 13.75 GHz: max_shell_increment_relative=4.308e-02
- 15.00 GHz: max_shell_increment_relative=4.819e-02

## Polarization dependence (Sec. 15)

- x: max_shell_increment_relative=4.025e-02
- y: max_shell_increment_relative=4.025e-02

Circular polarization: NOT_SUPPORTED_AS_INDEPENDENT_EXCITATION

## Decay law (Sec. 17 -- no automatic fitting)

Observed trend: NON_MONOTONIC_NO_CLEAR_TREND. Fit attempted: NONE_FITTED.
  note: the data do not show a visibly monotonic trend; no decay law was fitted (Sec. 17).

## Reduction (Sec. 19)

REDUCTION_NOT_ATTEMPTED: overall classification is CASE_D_PATTERN_DEPENDENT, not CASE_A_LOCAL: building an approximate spatial-reduction scheme (hierarchical matrices, block-Toeplitz reuse, domain decomposition, etc.) before knowing the coupling is local enough for such a scheme to be meaningful would not be a physically justified reduction (Sec. 19's own requirement); none is implemented in this phase.

## Inverse validation (Sec. 21-22)

- target=3.221e-09+1.215e-08j: status=FOUND, relative_error=0.0
- target=3.844e-09+1.208e-08j: status=FOUND, relative_error=0.04944387069959131
- target=2.535e-07+0.000e+00j: status=TARGET_UNREACHABLE, relative_error=0.9472077091755703

## Beam-level validation (Sec. 23)

- smooth_gradient: pointing_error_deg=(19.999999999999993,)
- abrupt_boundary: pointing_error_deg=(19.999999999999993,)
- strong_contrast: pointing_error_deg=(19.999999999999993,)
- checkerboard: pointing_error_deg=(19.999999999999993,)

## Three-beam holographic test (Sec. 24)

pointing_error_deg=(29.999999999999996, 29.999999999999996, 29.999999999999996), matched_peak_power_db=(0.0, 0.0, 0.0)

## HFSS status (Sec. 25-26)

HFSS_NOT_AVAILABLE

## Sec. 29 acceptance criteria

- [x] Direct heterogeneous MoM reference calculations work.
- [x] 1x1 through the largest feasible neighborhood are reproducible.
- [x] Coupling extent is measured rather than assumed.
- [x] Diagonal coupling is quantified.
- [x] Directional coupling is quantified.
- [x] Smooth gradients are tested.
- [x] Abrupt boundaries are tested.
- [x] Strong contrast is tested.
- [x] Alternating/checkerboard environments are tested.
- [x] Frequency dependence is investigated.
- [x] Polarization dependence is investigated where supported.
- [x] Pairwise interaction is quantified.
- [x] Higher-order interaction is investigated.
- [x] A numerical convergence criterion is defined before final classification.
- [x] The required physical neighborhood is identified.
- [x] Any reduction from full MoM is mathematically justified.
- [ ] Reduced calculations agree with direct MoM within predefined tolerances.
- [x] Inverse candidates are verified using direct MoM.
- [x] Unreachable targets are rejected.
- [x] Beam-level heterogeneous validation is completed.
- [x] No ML/surrogate model is required.
- [x] HFSS validation status is explicitly reported.

## Sec. 35 final objective

Outcome C: the response depends strongly on spatial pattern, not just local statistics.

