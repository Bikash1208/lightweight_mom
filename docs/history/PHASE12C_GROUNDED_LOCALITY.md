# Does local periodicity hold for the grounded 8 mm reflectarray cell? (Phase 12C)

Implementation: `src/phase12c/` (`design_point.py`, `grounded_supercell.py`, `locality_study.py`,
`coupling_study.py`, `spot_check.py`, `gate.py`, `exports.py`, `specification.py`, `plots.py`,
`report.py`, `run.py`); core extensions `src/cad_export.py` (`export_supercell_step` and helpers)
and `src/phase9/hfss_reference.py` (one-port ingestion). Design point recorded as
`config_phase12c_grounded_8mm.yaml`, loaded (never hard-coded) by `src/phase12c/design_point.py`.
Conventions: PHASE9_VALIDATION.md, PHASE11_HETEROGENEOUS_MOM.md, PHASE12A_GROUND_PLANE.md,
PHASE12B_MESH_PERIOD_STUDY.md.

## 1. Why

Phases 9-11 found local periodicity invalid and coupling non-local for the **ungrounded** stack,
but only on the 8-cell/side (0.6 mm) mesh that Phase 12B later showed cannot even distinguish two
of Phase 10's own named centers (C2 and C3 mesh-collide there). Phase 12B's own period study
favours a **grounded**, **P = 8.0 mm** design point (279.2 deg coverage at 12.5 GHz, resonance near
L = 6.0 mm, grating onset above band even at 45 deg, full surface-wave margin). This phase asks the
question Phase 13 needs answered before it can choose a forward model: does local periodicity hold
for *this* grounded design, at a mesh fine enough to trust, well enough to build a lookup table --
or does the earlier non-locality finding survive into the grounded case too, requiring a
full-aperture integral-equation solver?

Working design point (recorded, not hard-coded): RT/Duroid 5880, h = 1.575 mm, PEC ground
(`backing: pec`), P = 8.0 mm, rectangular PEC patch, 12.5 GHz design frequency, linear x/y
polarization only.

## 2. Task A: grounded supercells

`GroundedSupercellUnitCell` wraps Phase 9's frozen `SupercellUnitCell` (adds only
`backing = "pec"`, delegates everything else) -- the same duck-typed-wrapper pattern Phase 12A
established (`getattr(unit_cell, "backing", "air")` in the solver); `SupercellUnitCell` itself is
untouched.

**Floquet order matters more than it looks like it should.** A rank-sufficiency-only order
(`2*floquet_order+1 >= max(cells_x, cells_y)`, CONVENTIONS.md's mesh rule) is a well-posedness
condition, not an accuracy guarantee. A deliberate scan at fixed mesh (floquet_embedded = 34 -> 120,
n=3, homogeneous case) showed monotonic, exponential convergence of the correctness-proof relative
error: 1.98e-2, 7.49e-3, 5.35e-3, 3.87e-3, 1.50e-3, then **1.51e-14 at floquet = 120 = 2 x total
mesh cells**. Only `floquet_order = 2 * total_mesh_cells` (Phase 9's own original convention)
reaches machine precision; this is the convention used throughout Tasks A-C.

Homogeneous-supercell correctness proof (n=3, 4 cells/subcell, L=W=4.0 mm, this floquet
convention): relative error **1.411e-15** -- matches Phase 9's own ~4e-15 machine-precision
benchmark exactly, confirming the wrapper adds no numerical error of its own.

Mesh-convergence series (same homogeneous case, aligned mesh doubled each step):

| cells/subcell | relative difference from previous |
|---|---|
| 8 | 0.3318 |
| 16 | 0.1225 |
| 32 | 0.0282 |

Monotonically decreasing, consistent with first-order mesh convergence continuing past the
20-cell/subcell mesh used for Tasks A-C's production runs.

Largest feasible supercell size (same homogeneous case, aligned 20-cell/subcell mesh,
`floquet_order = 2 * total_mesh_cells`):

| n | n_unknowns (embedded) | runtime [s] |
|---|---|---|
| 1 | 180 | 0.13 |
| 3 | 1,620 | 6.3 |
| 5 | 4,500 | 104.6 |

Every case above reproduces the single isolated cell to 1e-14-1e-15 relative error, as it must for
a homogeneous supercell. **n=7 was attempted and killed after 42+ minutes without completing** --
solve cost scales roughly as the 6th power of total per-side cell count once Floquet order tracks
mesh size, and n=7 at this convention (280 cells/side, floquet=560) is infeasible within a
reasonable interactive budget. n=5 (104.6 s) is reported as the largest *routinely* feasible size,
per the brief's own instruction to stop at the largest feasible size rather than coarsen silently.

## 3. Task B: local-periodicity test, grounded 8 mm

Reference reflection-phase curve (S11_xx, aligned 20-cell/subcell mesh, the phase's own curated
`L_FAMILY_MM` = 1.6-7.2 mm in 0.8 mm steps, chosen specifically so its smallest fully-aligned mesh
stays small):

| L [mm] | \|S11\| | phase [deg] |
|---|---|---|
| 1.6 | 0.999887 | -21.28 |
| 2.4 | 0.999861 | -22.83 |
| 3.2 | 0.999787 | -26.59 |
| 4.0 | 0.999582 | -35.13 |
| 4.8 | 0.998944 | -55.43 |
| 5.6 | 0.996855 | -110.46 |
| 6.4 | 0.996414 | 142.54 |
| 7.2 | 0.998770 | 78.31 |

Centres, chosen from this measured curve, never assumed: **off-resonance L=2.4 mm**,
**near-resonance L=5.6 mm** (steepest phase slope between neighbours), **at-resonance L=6.4 mm**
(largest |S11| dip / phase wrap in the family).

Smooth-gradient cases (Phase 9's own dipole-moment embedded-vs-isolated metric; achieved step is
the actual neighbour-to-neighbour phase difference realized by the curated L family, not a
requested one):

| case | achieved step [deg] | phase error [deg] | magnitude error [dB] | within tolerance (<=10 deg / <=0.5 dB) |
|---|---|---|---|---|
| off_resonance_1d_row | 2.7 | 0.11 | 0.049 | **yes** |
| near_resonance_1d_row | 81.0 | 6.70 | 1.72 | no |
| near_resonance_2d_gradient | 107.0 | 0.20 | 9.83 | no |
| at_resonance_1d_row | 85.6 | 9.19 | 12.41 | no |

**Gate: GROUNDED_LOCAL_PERIODICITY_VALID_FOR_SMOOTH_GRADIENTS** -- local periodicity holds within
the stated tolerance only up to the smallest tested step, **2.7 deg/cell**; every larger tested
step (all >= 81 deg, since the curated L family's own resonance region has no finer intermediate
values) already fails on phase, magnitude, or both. Phase 9's own criterion (relative-complex-error
substantiality) separately classifies every one of these cases as INVALIDATED -- the two criteria
agree that gradients through or near resonance break local periodicity; they differ only in how
generously the far-off-resonance, small-step case is treated. Both are reported, per the brief's
instruction never to let one tolerance obscure the other.

Strong-contrast case (Phase 9's own extreme, for direct ungrounded/grounded comparison): center
L=5.6mm, opposite extreme L=1.6mm (farther from center than the alternative extreme, matching
Phase 9/11's own convention), achieved step 89.2 deg, relative complex error **0.400** -- large,
but markedly smaller than Phase 9's own ungrounded 0.997 figure at its coarser mesh (Task D below
shows this gap narrows further once the ungrounded case is itself re-measured on a fine mesh).

Periodic phase-gradient supercell (6 cells, spanning the curated L family in phase order): achieved
span **260.4 deg** -- not a full 2π, an honest finding, not forced. Per-cell embedded-vs-isolated
relative complex error ranges **0.52-1.33** across the 6 cells, worst at the smallest/first cell
(L=1.6mm, where the neighbour jump into L=4.8mm is largest in relative terms) and improving
somewhat but never dropping below 0.5 as L approaches resonance -- i.e. even a *smoothly, exactly
periodic* multi-step phase ramp does not recover local periodicity once individual neighbour steps
are large, confirming the smooth-gradient finding above from an independent construction.

## 4. Task C: coupling extent, grounded

Growing-pattern coupling series (Phase 11's own `coupling_radius`/`convergence`/`decay_law`
functions, reused unmodified) for off-resonance (L=2.4mm) and near-resonance (L=5.6mm) centers,
each with a gradient pattern (not only checkerboard) and a checkerboard pattern, sizes n=1/3/5,
20 cells/subcell:

| series | convergence | last shell increment (relative) |
|---|---|---|
| off_resonance_gradient | NOT_CONVERGED | 9.93e-02 |
| off_resonance_checkerboard | NOT_CONVERGED | 2.01e-03 |
| near_resonance_gradient | NOT_CONVERGED | 1.44e-01 (non-monotonic: 3.25e-02 -> 1.44e-01) |
| near_resonance_checkerboard | NOT_CONVERGED | 8.00e-02 |

**None of the four tested series is converged by n=5** -- matching Phases 9-11's own ungrounded
non-locality finding qualitatively. The near-resonance gradient series is non-monotonic
(increment grows from n=1->3 to n=3->5), so per Phase 11's own rule, no decay law is fitted for it
(`NONE_FITTED`); the other three decrease monotonically and a numerical exponential-candidate curve
is reported for them as a compact summary only, not as an asserted physical law. Grounded slabs
support a TM0 surface wave (unlike the ungrounded stack), and the off-resonance checkerboard series
does decay noticeably faster (2.0e-3 by n=5) than the off-resonance gradient series at the same
center (9.9e-2) -- consistent with, but not conclusive proof of, coupling character depending on
pattern type as well as proximity to resonance, exactly as Phase 11's own checkerboard-vs-block
finding (Case D) already established for the ungrounded case.

**Gate: GROUNDED_COUPLING_NOT_LOCAL_FOR_TESTED_SIZES** -- coupling has not localized by n=5 for any
tested grounded center/pattern combination.

## 5. Task D: spot-check on a finer aligned mesh (read-only)

Re-solved, without touching any earlier phase's code or output files:

| case | mesh (cells/subcell or /side) | verdict | detail |
|---|---|---|---|
| Phase 10 C2 vs C3 | 24 | **EARLIER_CONCLUSIONS_CHANGED** | mesh-**distinct** at this resolution (relative isolated-response difference 0.547); Phase 12B's coarse-mesh mesh-identical finding does not hold here |
| Phase 9 strong contrast | 48 | **EARLIER_CONCLUSIONS_CHANGED** | relative complex error drops from 0.997 (>= 0.5 threshold, INVALIDATED) at 8 cells/subcell to **0.218** (< 0.5) at this finer, aligned mesh -- the INVALIDATED classification itself does not survive; independently verified: embedded self-convergence (floquet 80->120) changed by only 0.41%, isolated self-convergence (floquet 32->96) changed by only 0.29%, so this drop is a genuine finding, not a floquet artifact |
| Phase 11 Case D (checkerboard vs block) | 48, n=3 (reduced from n=7) | EARLIER_CONCLUSIONS_HOLD | relative difference **21.7** -- still overwhelmingly substantial; arrangement, not just local statistics, still matters |
| Phase 11 coupling-radius series | 48, n=1->3 | EARLIER_CONCLUSIONS_HOLD | shell increment 7.83e-02 -- still far from converged at this (necessarily reduced) size range |

**This is the single most consequential result in this phase.** Phase 9's original
`LOCAL_PERIODICITY_INVALIDATED` finding rested on a strong-contrast relative complex error of 0.997
measured at an 8-cell/subcell mesh that Phase 12B separately showed cannot even distinguish some
Phase 10 geometries. At a properly aligned, much finer 48-cell/subcell mesh, that same strong
contrast case drops to 0.218 -- below Phase 9's own 0.5 substantiality threshold. Both the C2/C3
mesh-collision reversal and this drop point the same way: **the coarse unaligned mesh used in
Phases 9-11 systematically overstated inter-cell coupling**, echoing (more dramatically) Phase
12A/12B's own earlier "unaligned mesh overstated an effect" pattern. This does **not** overturn
Phase 11's own Case D or coupling-radius findings, which remain substantial/not-converged even at
this finer mesh (reduced size range) -- so non-locality itself is not disproven, but the specific
magnitude of Phase 9's original strong-contrast evidence for it is now known to have been
mesh-inflated.

## 6. Task E: HFSS-ready exports

5 single-cell STEP exports (L = 1.6, 3.2, 4.8, 5.6, 6.4 mm, spanning the reference curve including
resonance, each with the PEC ground solid included), 1 uniform 3x3, 1 smooth-gradient 3x3, 1
strong-contrast 3x3 supercell STEP export (9 metal solids + 1 ground solid each, `Substrate`,
`Metal_k`, `Ground` named parts, via the new `export_supercell_step`). A synthetic, clearly-labelled
one-port HFSS CSV (`SYNTHETIC` in every `cell_id`, `s21=None` for every record) was ingested by
`phase9.hfss_reference.load_hfss_reference` to confirm one-port (`--hfss-reference`, no S21 columns)
files are now accepted -- **no real HFSS results exist anywhere in this project.**

**Gate: HFSS_EXPORT_READY.**

## 7. Final statuses

```
Task A (grounded supercells):       GROUNDED_SUPERCELLS_VALIDATED
Task B (local periodicity):         GROUNDED_LOCAL_PERIODICITY_VALID_FOR_SMOOTH_GRADIENTS
                                     (holds up to 2.7 deg/cell measured; every larger tested
                                     step, all near/at resonance, fails)
Task C (coupling extent):           GROUNDED_COUPLING_NOT_LOCAL_FOR_TESTED_SIZES (n<=5)
Task D (fine-mesh spot-check):      EARLIER_CONCLUSIONS_CHANGED (Phase 9 strong-contrast
                                     INVALIDATED finding and Phase 10 C2/C3 mesh-identity both
                                     reversed at finer aligned mesh; Phase 11 Case D and
                                     coupling-radius findings hold)
Task E (HFSS-ready exports):        HFSS_EXPORT_READY
HFSS status:                        HFSS_NOT_AVAILABLE (unchanged)
Geometry inversion:                 GEOMETRY_INVERSION_LOCKED (this phase does not attempt to unlock it)
STEP export:                        STEP_EXPORT_LOCKED
```

## 8. Recommendation for Phase 13's forward model

Local periodicity, at this grounded 8 mm design point, holds within the stated tolerance only for
small neighbour phase steps (measured up to ~2.7 deg/cell) and degrades sharply approaching
resonance (Sec. 3); a reflectarray-style lookup-table forward model is therefore defensible **only
if** the aperture's own phase gradient per cell is kept at or below the measured max-step-within-
tolerance everywhere -- a real, checkable design constraint, not an assumption -- and the periodic
phase-gradient supercell (Sec. 3) shows this constraint is hard to satisfy near resonance
specifically, since an aperture whose phase profile must pass through that region cannot rely on a
pure lookup table there. A full-aperture integral-equation solver (a la Budhu & Grbic 2021) remains
the more broadly defensible choice unless the design is constrained to stay in the small-step-safe
region measured here. This is offered as a **recommendation, not a decision** -- Phase 13 makes the
actual choice.

## 9. Rules followed

Measured numbers only; the synthetic HFSS ingestion test data is clearly labelled as such
(`SYNTHETIC` cell IDs); no fabricated HFSS results exist anywhere. No earlier phase's code or
results were modified -- Phase 9's `SupercellUnitCell` is wrapped, not edited; Task D writes only
under `results/phase12c_*`, never touching `results/phase9_*`/`phase10_*`/`phase11_*`. No period or
geometry was chosen by phase coverage alone (Phase 12B already made that choice on its own,
independent grounds). No ML/fitted surrogate appears anywhere in the locality or coupling study
(Phase 11's own rule, reused); the one numerical exponential-candidate fit in Task C is applied only
to a monotonic series and is explicitly labelled a summary, not an asserted decay law, per Phase
11's own convention.

## 10. Running it

```bash
.venv/bin/python examples/run_phase12c.py   # results/phase12c_grounded_locality/, results/phase12c_hfss_export/
.venv/bin/python -m pytest -q tests/test_phase12c_*.py tests/test_hfss_reference_one_port.py \
    tests/test_cad_export_supercell.py
```
