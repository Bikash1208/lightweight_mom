# Physics-based heterogeneous MoM / nonlocal coupling (Phase 11)

Phase 10 found that a data-fitted surrogate does not solve the heterogeneous-metasurface problem
(`SURROGATE_INVALIDATED`). Phase 11 does not attempt another surrogate. It asks the physics
question directly, with the existing, unmodified Phase 1-3 MoM solver as the sole computational
engine: **how far through the metasurface does electromagnetic interaction from a center cell
remain significant, and is a finite neighborhood ever enough?**

No machine learning, no data-fitted regression, no learned geometry-to-response mapping appears
anywhere in `src/phase11/` (mechanically checked by `tests/test_phase11_run.py`, which greps the
actual source for forbidden terms).

## Method

Every "environment" is a real, solved heterogeneous supercell (`SupercellUnitCell`, Phase 9/10's
own duck-typed composite geometry, reused unmodified), assembled and solved through the exact same
`solver.assemble_system`/`solver.solve` (`Z I = V`) pipeline every phase of this project has used.
A new Phase-11-owned wrapper (`reference_solver.py`) adds a configurable excitation polarization
(needed for Sec. 15's polarization-dependence study) without touching Phase 10's own frozen
wrapper.

"How far does interaction reach" is measured by **growing one fixed pattern**: a checkerboard
between a center cell's own geometry and the domain's opposite extreme, at sizes 1x1 through 11x11,
where every smaller size is the exact centered truncation of the same larger pattern (verified: this
preserves the checkerboard's parity exactly regardless of the truncation offset, so 1x1/.../11x11
are all genuinely nested views of one structure, not independently regenerated ones).

## Real, measured findings

- **Coupling is not local.** The center response has not settled to within the solver's own
  measured self-consistency floor (`phase9.local_periodicity.NOISE_FLOOR_RELATIVE`, reused, not a
  new number) by 9x9/11x11 for any of the 5 tested center geometries.
- **Diagonal coupling is reproduced.** Literal reuse of Phase 10's own `neighbor_sensitivity.py`/
  `diagonal_study.py` gives mean|diagonal|/mean|cardinal| in the same 25-70%-ish range Phase 10
  measured, confirming that finding independently rather than merely repeating the claim.
- **Spatial arrangement matters beyond neighborhood size (Sec. 18 Case D).** Two patterns with
  IDENTICAL local statistics (the same two geometries, the same 50/50 split) but different
  arrangements -- a distributed checkerboard vs. a single contiguous A|B block -- produce
  substantially different center responses. This means the coupling-radius series alone (which only
  answers "how big a neighborhood") is not sufficient on its own: the actual spatial pattern must be
  retained, not just its size or local density.
- **Near-neighbor pairwise interactions are close to additive** (interaction error a few times
  1e-4 relative to the baseline, for cardinal/diagonal/near-far pairs at strong contrast) and the
  measured 3-body residual beyond pairwise is smaller still (~1e-5 relative) -- additivity holds
  well at the level of a *few* nearby perturbations, even though the *aggregate* effect of many
  cells arranged differently (the Case D finding above) does not reduce to a simple neighborhood-
  size rule. These are not contradictory: many small, nearly-additive pairwise terms can still sum
  to a large, arrangement-dependent total.
- **Frequency**: the coupling extent (measured via the maximum shell-increment relative error)
  increases mildly and monotonically from 10 to 15 GHz -- coupling is somewhat stronger at the high
  end of the band, not constant across it.
- **Polarization**: x- and y-polarized excitation give identical coupling-extent numbers for an
  isotropic (L=W) center -- the expected physical symmetry, not "no polarization dependence" as a
  general claim (an anisotropic center would need its own check).
- **No decay law was assumed or fitted** without first checking the data (Sec. 17): the raw
  distance-vs-interaction series is reported directly, and a candidate exponential is only fit when
  the series is visibly monotonic (with fitted parameters, R^2 and residuals reported alongside).

## Reduction, inverse, beam-level

Because the overall classification is not `CASE_A_LOCAL`, no approximate spatial-reduction scheme
(hierarchical matrices, domain decomposition, etc.) was implemented (Sec. 19: building one before
knowing the coupling is local enough to justify it would not be a physically defensible reduction).
The one reduction this project already relies on -- solving a composite grid once and reading out
many sub-cells from that single solve (`reference_solver.solve_grid`) -- is exact, not approximate,
and is used throughout.

First-stage geometry inversion (Sec. 21-22) uses a deterministic, non-ML method: every candidate
geometry is ranked purely by its ISOLATED-CELL response (a real physics computation, zero fitting),
and the winner is verified against the real heterogeneous MoM solver in its actual fixed
environment. Reachable and deliberately-unreachable targets are both correctly classified.

Beam-level validation (Sec. 23-24, including the three-beam holographic test) is generated entirely
from real, direct heterogeneous MoM solves of small (7x7) apertures -- never from a lookup table or
surrogate.

## Gates and final status

Geometry inversion requires 14 explicit conditions (Sec. 30); two of the most consequential are, and
remain, false: no independent HFSS/measurement reference exists anywhere in this project (unchanged
across Phases 9-11), and the heterogeneous forward model does not converge for the classification
found. **Geometry inversion and final STEP export remain LOCKED.**

## Running it

```bash
.venv/bin/python scripts/research/run_phase11_heterogeneous_mom.py     # results/phase11_heterogeneous_mom/
```

## Verification

`tests/test_phase11_*.py`: pattern generators, the truncation-preserves-parity property, the
convergence criterion (declared from Phase 9's own already-sourced constants, before any Phase 11
result existed), CASE A/B/C/D classification, pairwise/higher-order interaction arithmetic, the
14-/6-condition gates (including the two mandatory Sec. 28 negative tests), and an end-to-end
(reduced-scale, real-solve) run of the full pipeline -- plus a direct source-grep confirming no
forbidden ML term appears anywhere in `src/phase11/`.
