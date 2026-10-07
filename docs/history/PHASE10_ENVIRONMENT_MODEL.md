# Environment-aware electromagnetic response model (Phase 10)

Entry status (Phase 9): `LOCAL_PERIODICITY_INVALIDATED` — the isolated periodic-cell model can be
wrong by an amount comparable to the signal's own magnitude under realistic strong-contrast
neighbor conditions. Phase 10 asks the question Phase 9 could not: **how much neighborhood
information does the electromagnetic system actually require, and can a finite, validated
environment-aware model replace the isolated-cell lookup for holographic synthesis?**

Two independent, real capabilities were built on top of Phase 9's own heterogeneous-supercell MoM
solver (`phase9.supercell_geometry`/`supercell_solver`, unmodified): a **neighborhood-sensitivity
study** (10A–10C) that measures, rather than assumes, how many neighbor cells matter; and an
**environment-response dataset + surrogate + validation pipeline** (10D–10F) that turns that
finding into a testable forward/inverse/beam-level model. The "center-cell response" used
throughout is the same physical quantity Phase 9B validated: the local induced-current dipole
moment of one sub-cell's footprint (complex, units A·m) — **not** S11 or S21, since that question
remains unresolved for lack of HFSS data (unchanged from Phase 9A).

## 10A–10C: how many neighbors actually matter?

For all 5 center cells (`C1`=1.4mm, `C2`=(2.0,2.8)mm, `C3`=2.8mm, `C4`=(3.6,2.8)mm, `C5`=4.2mm) and
4 contrast levels (weak ±0.2mm, moderate ±0.6mm, and the sweep's own 1.0mm/4.4mm extremes tested
*separately*, since how "strong" an extreme is depends on the center), the real, measured findings:

- **Diagonal neighbors are not negligible.** Mean |diagonal contribution| / mean |cardinal
  contribution| ranges ~25–70% across every center and contrast level (`diagonal_cases.csv`,
  `cardinal_vs_diagonal`) — comparable to the cardinal contribution, not a small correction.
- **A cardinal-only ("5-cell") model is measurably insufficient.** Ignoring the diagonals produces
  a predictive error against the true 9-cell response that reaches **>300–490% of the center
  cell's own baseline response magnitude** under the domain's strong-contrast extreme, for every
  center tested.
- **Even the full 3×3 ("9-cell") neighborhood may not be enough.** Distance-2 (two-period)
  neighbors still produce a measurable effect above the solver's own numerical floor
  (`long_range_cases.csv`) for 4 of 5 centers.
- **Final per-center classification** (`neighborhood_study/neighborhood_decision.json`):

  | center | (L, W) mm | decision |
  |---|---|---|
  | C1 | (1.4, 1.4) | `9_CELL_REQUIRED` |
  | C2 | (2.0, 2.8) | `LARGER_NEIGHBORHOOD_REQUIRED` |
  | C3 | (2.8, 2.8) | `LARGER_NEIGHBORHOOD_REQUIRED` |
  | C4 | (3.6, 2.8) | `LARGER_NEIGHBORHOOD_REQUIRED` |
  | C5 | (4.2, 4.2) | `LARGER_NEIGHBORHOOD_REQUIRED` |

  Classified using the same, already-sourced, non-arbitrary criteria Phase 9 established (a
  measured solver noise floor for "sufficient"; "at least half the signal's own magnitude" for
  "substantial") — never a newly invented percentage.

  A known, disclosed limitation: a 0.2mm "weak" perturbation is sometimes exactly at this
  project's mesh-staircase resolution floor (0.6mm at `mesh_cells_per_subcell=8`) for small-patch
  centers, giving a measured-zero delta that is a discretization artifact, not zero physical
  coupling (confirmed directly by re-solving at finer mesh; flagged per-case in
  `discretization_floor_flags`, never silently hidden).

Given `LARGER_NEIGHBORHOOD_REQUIRED` for 4/5 centers, the dataset and surrogate below are still
built on the full 9-cell (3×3) representation — the largest tractable option — with this residual
limitation explicitly disclosed rather than papered over.

## Environment-response dataset

640 real, solved samples (structured, harvested from the 10A–10C study at zero extra cost, +
random combinatorial 9-cell environments, fixed seed 1729): 353 train / 80 validation / 207 blind
test. Two safeguards against the two named risks:

- **Geometry-based holdout (Sec. 21):** center `C5` (4.2, 4.2mm) is reserved *entirely* for blind
  testing — never seen in training, testing genuine extrapolation.
- **Symmetry-safe grouping (Rule 4):** samples related by a y-reflection, x-reflection or
  180°-rotation of the 3×3 environment are grouped and forced into the *same* split, so a model can
  never be tested on a near-duplicate of something it trained on.

## Baselines and surrogate (blind test, n=207)

| model | mean relative error | p95 relative error | mean phase error [deg] |
|---|---|---|---|
| Baseline 1 (isolated cell, ignores all neighbors) | **0.574** | 1.199 | 30.9 |
| Baseline 2 (linear additive neighbor correction) | 1.008 | 3.309 | 40.6 |
| Baseline 3 / surrogate (ridge, quadratic + opposite-pair cross terms) | 1.350 | 5.043 | 39.0 |

Full percentile reports (mean/median/RMS/p95/max, complex/relative/magnitude/phase, never only the
mean) are in `surrogate/error_metrics.csv`. **Neither fitted model (additive or nonlinear ridge)
outperforms the naive, already-invalidated isolated-cell baseline on relative error** — the one
non-arbitrary criterion Sec. 27 permits (no invented accuracy threshold; judged only against the
pre-Phase-10 status quo). This is a real, measured failure to generalize, not a coding defect: it
is consistent with the 10A–10C finding that even the 9-cell representation this surrogate is built
on is measurably insufficient for most of the domain, compounded by the deliberately-hard blind
test (an entire unseen center geometry, plus symmetry-safe isolation).

Two spatial-structure generalization tests (Sec. 30–32, one 7×7 heterogeneous-supercell MoM solve
each, all 25 interior cells read from that single solve):

- **Smooth gradient** (L varying 1.4→4.2mm across the aperture): mean |error| ≈1.4e-9 A·m — a
  modest error, the surrogate's best showing.
- **Abrupt A|B discontinuity:** error is markedly larger immediately at the boundary than farther
  from it (`surrogate/abrupt_cases.csv`) — the surrogate degrades sharply exactly where a real
  hologram's geometry changes fastest.

## Inverse validation (Sec. 38–40)

7 test targets, every proposed geometry independently re-checked by the real heterogeneous MoM
solver (never accepted on the surrogate's own say-so): 5/7 REACHED (relative error <0.5 for
same-environment sanity checks), 2/7 `TARGET_UNREACHABLE` — including one deliberately extreme
(20×-library-max) target correctly rejected rather than forced to a nonsensical "closest" answer.

## Beam-level validation (Sec. 28–29)

One environment-aware aperture (the gradient grid's 25 interior cells), far field computed with
Phase 8E's own unmodified planar-aperture radiation integral for both the MoM-reference and
surrogate-based apertures: detected beam pointing matches exactly (17.999...° both), mean complex
far-field error ≈3.1e-13, mean phase error 1.77°. This one scenario shows little beam-level
degradation, but with only 25 aperture elements and a single (smooth-gradient) spatial pattern
tested, it does not contradict the abrupt-discontinuity finding above.

## Gates and final status

| gate | status |
|---|---|
| 10A — control variable | **CLOSED** (no HFSS/measurement reference anywhere in this repository, unchanged from Phase 9A) |
| 10B — environment model | OPEN (every required study step completed) |
| 10C — surrogate forward | OPEN (blind test, strong-contrast, unseen-environment, gradient and abrupt tests all completed) |
| 10D — inverse validation | OPEN (implemented, MoM-verified, unreachable targets correctly rejected) |
| 10E — beam-level | OPEN (aperture/far-field/pointing/power all compared) |

Every gate tracks whether the validation *process* was completed, independent of the *outcome* —
Gate 10B is OPEN even though its own finding is that neither 5 nor 9 cells is fully sufficient.

**Final status: `SURROGATE_INVALIDATED`.** Geometry inversion and final STEP generation remain
**LOCKED** — both independently (Gate 10A's persistent lack of a control variable, and the
surrogate's own measured failure to beat the status quo).

## Running it

```bash
.venv/bin/python examples/run_phase10_environment_model.py     # results/phase10_environment_model/
```

## Verification

`tests/test_phase10_*.py` (79 tests): the 5-cell/9-cell test matrices, the discretization-floor
flag (verified directly against a finer-mesh re-solve), symmetry-safe train/test grouping, the
geometry-based C5 holdout, surrogate domain-rejection and save/load round-tripping, a synthetic
purely-additive recovery test for Baseline 2, the mandatory Sec. 47 negative tests (unvalidated
control variable, and a surrogate that fails its blind test, each independently keep geometry
inversion locked), and an end-to-end (reduced-scale, real-solve) run of the full pipeline.
