# Resolve the local-periodicity step limit; ingest the real HFSS data (Phase 12D)

Implementation: `src/phase12d/` (`reflection_metric.py`, `step_curve.py`, `floquet_power.py`,
`aperture_test.py`, `hfss_ingest.py`, `gate.py`, `specification.py`, `plots.py`, `report.py`,
`run.py`); core extension `src/phase9/hfss_reference.py` (`read_touchstone`, a general N-port
Touchstone reader). Design point unchanged from Phase 12C: grounded RT/Duroid 5880, h = 1.575 mm,
P = 8.0 mm, square PEC patch, 12.5 GHz. Conventions: PHASE12C_GROUNDED_LOCALITY.md,
`results/phase12c_grounded_locality/task_b.json`, `hfss_exports/phase12_grounded_8mm/README.md`.

## 1. Why

Phase 12C's own "local periodicity holds up to 2.7 deg/cell" came from only four gradient cases,
one at 2.66 deg/cell (passed) and three at 81-107 deg/cell (all failed) -- the 0.8 mm-step L grid
used there cannot produce anything in between near resonance, and the failing cases were not
reported in the one metric that actually matters for a design: reflection phase/amplitude error,
not the raw induced-current dipole-moment deviation. Real HFSS Floquet-port data also arrived for
one case (G1, L = 2.0 mm) between Phase 12C and this phase, and Phase 12C's own Task E HFSS
ingestion never ran. This phase (A) re-measures the step limit with the design-relevant metric and
an honest accounting of what mesh resolution is actually achievable, (B) builds the real
aperture-level test -- a full-supercell MoM solve decomposed into every Floquet order, not the
dipole proxy, which is what actually decides Phase 13 -- and (C) ingests the real HFSS data.

## 2. Task A: the step-size curve, with the design-relevant metric

**The dipole-to-reflection-coefficient conversion** (`reflection_metric.py`) reuses
`sparameters.scattering_from_fields` unmodified: it already accepts a supplied mean-current array,
so plugging in the exact Floquet coefficient of ANY current (isolated or embedded) yields the S11
that current would present if periodic. Validated exactly (machine precision) against
`compute_s_parameters` when fed the true Floquet coefficient; applying it to Phase 9/12C's own
dipole-moment convention (a 48x48 Riemann-sum approximation of the same integral) introduces a
small, separately-measured ~0.1-0.2 deg / <0.1% residual, negligible against this phase's 5-10 deg
tolerances.

**Mesh cost** (measured directly, the reason the brief's 0.05-0.1 mm aspiration is not reached):
a 3x3 embedded supercell's cost scales roughly as the 6th power of per-subcell mesh count once
floquet order tracks it (Phase 12C's own finding, confirmed again):

| cells/subcell | runtime (L=6.0mm case) |
|---|---|
| 24 | 76.8 s |
| 28 | 196.4 s |
| 32 | 409.5 s |
| 40 | did not finish in 14+ min; killed |

**24 cells/subcell** (aligned step = 2 x period / 24 = 0.667 mm) is used throughout Task A -- finer
than Phase 12C's own 20-cell/0.8 mm family, but a 0.05-0.1 mm-aligned family near this design's
~6 mm resonance would need a mesh far beyond 40 cells/subcell, genuinely infeasible for an embedded
solve at this project's validated 2x-total-mesh floquet convention within a practical budget. Cost
also depends strongly on patch size at fixed mesh (fewer metal-adjacent edges for a small patch):
off-resonance cases (L around 2 mm) solved in seconds; near/at-resonance cases (L up to 7.3 mm,
close to filling the whole 8 mm period) took up to ~2-3 minutes of real compute each.

Reference curve (self-discovered resonance, not assumed): resonance falls at **L = 6.0 mm** on this
24-cell/subcell mesh (|S11| = 0.9957, phase -164.5 deg) -- consistent with Phase 12B/12C's own
~6.0-6.4 mm finding at their coarser meshes. Centres: off-resonance L=2.0 mm, near-resonance
L=6.667 mm, at-resonance L=6.0 mm.

20 step cases (1-D row and 2-D gradient, k=1-2 family-steps in both directions, off/near/at
resonance -- see `PHASE12D_LOCAL_PERIODICITY.md`'s companion `task_a.json` for the full table):

| region | achieved steps [deg] | reflection phase error [deg] |
|---|---|---|
| off-resonance | 0.84 - 13.45 | 0.000 - 0.039 |
| near-resonance | 41.70 - 160.32 | 1.1 - 122.4 |
| at-resonance | 77.88 - 124.14 | 19.5 - 173.0 |

**Gate: `LOCAL_PERIODICITY_STEP_LIMIT_MEASURED`** -- holds within the 10 deg/1 dB design tolerance
(and the stricter 5 deg/0.5 dB one) up to **5.8 deg/cell**, the largest tested step that stayed
within it (off-resonance only). The data are explicitly **not monotonic** in step size near
resonance (e.g. a 41.7 deg near-resonance step gives 87-97 deg phase error, but an 82.4 deg step at
the same centre gives only 19.9 deg) -- reported as measured, nothing fit to it. The true
near-resonance crossover step is therefore **bounded, not pinned to a single number**: every step
actually reachable near/at resonance at feasible mesh resolution already exceeds the design
tolerance, so the practical conclusion (a lookup table cannot be trusted through this design's
resonance region) does not depend on knowing the exact crossover value.

## 3. Task B: the aperture-level test (the metric that actually decides Phase 13)

`floquet_power.py` computes, directly from the solved surface current
(`basis_functions.RooftopBasis.current_coefficients` at any Floquet harmonic, fed to
`green_function.LayeredGreenFunction.sheet_fields`), the reflected power in every Floquet order of
the aperture's own (super-period) lattice -- no dipole proxy. Validated on both cases the brief
requires: a full-cell PEC sheet (uniform current -> every order but (0,0) is **exactly** zero to
machine precision, |S11|=1) and a uniform supercell (order (0,0)'s fraction reproduces the true
isolated |S11|^2 to 1e-9, every order not a multiple of n is exactly zero).

Periodic phase-gradient supercells (N=6, 8, 12 cells per assumed 2*pi cycle, cells picked from a
real isolated-cell lookup table, reusing `phase12c.locality_study.build_periodic_phase_gradient_
supercell` unmodified), full supercell solved (not just the centre dipole), reduced floquet order
(`total_cells//2+8`, Phase 12C's own Task-D convention for large qualitative effects) -- **all
three sizes were cheap** (7-125 s; N=12 was NOT infeasible, unlike some of Phase 12C's own Task A
large supercells):

| N | achieved span [deg] | n_unknowns | runtime [s] | fraction in intended (1,0) order | fraction specular (0,0) | largest spurious order |
|---|---|---|---|---|---|---|
| 6 | 186.1 | 2,760 | 7.0 | 0.529 | 0.080 | (-1,0): 0.382 |
| 8 | 186.1 | 4,640 | 18.5 | 0.411 | 0.141 | (-1,0): 0.343 |
| 12 | 186.1 | 10,656 | 124.2 | 0.475 | 0.062 | (-1,0): 0.303 |

**Gate: `LOOKUP_TABLE_INADEQUATE`** -- no tested N reaches even 53% of reflected power in the
intended order, let alone a 90% adequacy bar. The largest single spurious order at every N is
**(-1,0), the mirror-image (backward) diffraction order**, carrying 30-38% of the power -- a real,
physically informative failure mode: the lookup-table-constructed phase ramp reflects a large
fraction of the wave the wrong way. The achieved span (186.1 deg, fixed by the underlying L family
used, itself capped by this design's own real S11 curve) is the root cause, not the specific N: a
blazed grating needs a genuine 360 deg ramp per period, and this design point's own phase coverage
(already measured short of 360 deg by Phases 12A-12C) cannot supply one. More cells subdividing the
same 186 deg span more finely does not fix this -- consistent with N=8's fraction being *lower*
than N=6's, not higher.

## 4. Task C: real HFSS data ingestion

`read_touchstone` (new, general N-port Touchstone reader, `phase9/hfss_reference.py`) parses the
real HFSS export exactly: 36 ports (18 Floquet modes each on the top and bottom air-box ports), 401
frequency points (10-15 GHz), MA format, verified against the file's own token count. Reproduces
the pre-existing `hfss_results/G1_comparison.csv` exactly (independent verification, not
regeneration).

Provenance (recorded, never invented): HFSS 2025 R1; Floquet port at z=+5 mm, no de-embed (matches
the MoM reference plane exactly); 18 modes/port; PEC ground; PEC patch 35 micron (CAD only; MoM
uses zero-thickness PEC); max delta-S 0.001; 10-15 GHz, 401 points; `exp(+j*omega*t)` both sides,
no conjugation. Setup lesson recorded: an earlier 2-mode-port/copper-ground run
(`unit_HFSSDesign1.s4p`, kept on disk) was off by up to -4.7 deg; the 18-mode PEC-ground setup
above is the one that agrees this closely.

**G1 (L=2.0mm) comparison**: mean phase difference **-0.294 deg**, max |dphase| **0.525 deg**
(10->15 GHz), max magnitude difference **7e-5 dB**, max cross-pol **2.14e-5**, top-to-bottom
transmission **exactly 0.0** (grounded, as required) -- **meets** the tolerance fixed before
looking at the data (1.0 deg / 0.05 dB; a generous order of magnitude above the measured agreement,
and one the cited earlier 2-mode/copper-ground run would have failed). Geometric identity check
(`cad_export.import_and_summarize`, re-imports both STEP files independently): Phase 12C's own
earlier single-cell exports (G2/L=4.8mm, G3/L=5.6mm, the two values it happens to share with this
run) are **bit-identical in metal volume and bounding box** to the ones actually sent to and
returned from HFSS.

**Gate: `HFSS_SINGLE_CELL_VALIDATED_PARTIAL`** -- G1 validated; G2, G3, G5, G6 have STEP geometry
but no HFSS Touchstone result yet; **G4 (the resonance case, L=6.0mm) is explicitly still
pending**, per the brief's own note. No supercell HFSS result exists (`HFSS_NOT_AVAILABLE_FOR_
SUPERCELLS`) -- Phase 12C's supercell STEP exports have no MoM-vs-HFSS comparison possible until
one arrives.

## 5. Final statuses

```
Task A (step-size curve):  LOCAL_PERIODICITY_STEP_LIMIT_MEASURED (up to 5.8 deg/cell, design and
                            strict tolerance both; near/at-resonance steps all exceed it)
Task B (aperture test):    LOOKUP_TABLE_INADEQUATE (N=6,8,12 all tested and feasible; none reaches
                            90% -- or even 53% -- power in the intended order)
Task C (HFSS ingestion):   HFSS_SINGLE_CELL_VALIDATED_PARTIAL (G1 only; G4 resonance case pending;
                            supercells HFSS_NOT_AVAILABLE_FOR_SUPERCELLS)
Geometry inversion:        GEOMETRY_INVERSION_LOCKED (unaffected)
STEP export:                STEP_EXPORT_LOCKED (unaffected)
```

## 6. Recommendation for Phase 13's forward model (replaces Phase 12C's "<=2.7 deg/cell" figure)

The step-size curve (Task A) shows the design-relevant reflection-phase error stays under tolerance
for every tested step out to double digits of degrees away from resonance (up to 5.8 deg/cell
measured), but every tested near- or at-resonance step -- the smallest achievable at this project's
validated floquet convention, itself tens of degrees because finer is computationally out of reach
-- already fails badly (19-173 deg phase error). The true near-resonance crossover step is
therefore bounded, not pinned to a single number, but the practical conclusion does not depend on
knowing it exactly: an aperture whose phase profile must pass through this design's resonance
region cannot rely on a lookup table there. **The aperture-level test (Task B), which actually
decides Phase 13, confirms this independently**, from a full-supercell MoM solve: at every tested
size, the lookup-table-constructed anomalous reflector puts less than 53% of its reflected power in
the intended order -- `LOOKUP_TABLE_INADEQUATE`, not merely "adequate above some N". The root cause
is visible in the numbers: this design point's own achievable phase span (186 deg) is far short of
the 360 deg a blazed grating needs, the same phase-coverage shortfall Phases 12A-12C already
measured from a different angle. **A full-aperture integral-equation solver (a la Budhu & Grbic
2021) is therefore the more broadly defensible default for this grounded 8 mm design point** --
though a design confined entirely to the flat, off-resonance part of the curve (steps below ~6
deg/cell) remains a narrow case where a lookup table could still work. This is a recommendation,
not a decision -- Phase 13 makes the actual choice.

## 7. Rules followed

Measured numbers only; the near/at-resonance step cases were computed once and reused (never
recomputed with different, more favourable inputs); no HFSS data fabricated or extrapolated --
Task C's tolerance and provenance were fixed before looking at the comparison numbers, and G4's
absence is reported as pending, not filled in. Earlier phases' code and results are untouched
(`phase12c.locality_study`'s grid builders and `classify_centers` are reused unmodified;
`phase9/hfss_reference.py` is extended additively, matching its own Phase 12C precedent). No
ML/fitted surrogate appears anywhere in this phase. Outputs are written only under
`results/phase12d_local_periodicity/`.

## 8. Running it

```bash
.venv/bin/python scripts/research/run_phase12d.py   # results/phase12d_local_periodicity/
.venv/bin/python -m pytest -q tests/test_phase12d_*.py tests/test_phase9_hfss_reference.py
```
