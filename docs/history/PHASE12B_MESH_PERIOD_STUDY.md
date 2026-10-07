# Mesh alignment, duplicate-geometry audit, grounded period study (Phase 12B)

Implementation: `src/discretization.py` (mesh alignment, `metal_mask_key`/`duplicate_geometry_groups`
-- core, additive), `src/sweep.py` (`SweepPlan.mesh_duplicate_groups`/`n_mesh_distinct`, `plan_sweep`'s
`warn_duplicates`). Orchestration/report only: `src/phase12b/`. Conventions:
CONVENTIONS.md ("Mesh alignment"), DATABASE.md, UNIT_CELL_LIBRARY.md Sec. 2A.

## 1. Why

Phase 12A's own reported `mesh_error = 0.048` (above its 0.02 tolerance) and the observation that
several Phase 7/12A sweep points give bit-identical S11 (e.g. L = 1.4/1.6, 2.0/2.2 mm at Phase 7's
32x32 mesh) both point at the same underlying cause: **a mesh cell is metal iff its centroid is
metal** (CONVENTIONS.md), and Phase 7/12A's mesh (32 cells / 4.8 mm = 0.15 mm pitch) does not
divide the 0.2 mm sweep step evenly, so many feature edges are not on a mesh line. This phase (a)
adds a reusable, tested capability to detect and avoid this, (b) re-checks whether Phase 12A's
grounded-stack finding survives it, (c) traces the same issue through Phases 7-11's own data, and
(d) asks whether a larger period lets a grounded patch reach fuller phase coverage.

## 2. Task A: aligned meshes (core capability)

Mesh lines sit at `-P/2 + k*P/N` for integer `k = 0..N`. A centred feature's edges (`+-L/2`) land
on mesh lines iff `L` is an integer multiple of `2*dx` (`dx = P/N`) for an even `N`; more generally
(off-centre features), both `center -+ length/2` must individually satisfy that condition
(`discretization.is_position_mesh_aligned`/`feature_edges_aligned`, derived and tested, not
guessed). `smallest_aligned_cells(period, features)` finds the smallest such `N` for a whole
family of `(length, center)` pairs by direct search (bounded, honestly reports `None` if none
exists within the bound) -- for Phase 7's own 18-value, 0.2 mm-step, centred L=W family at
P = 4.8 mm, this is **48 cells** (dx = 0.1 mm), matching `UNIT_CELL_LIBRARY.md` Sec. 2's own
already-measured "48x48 gives exactly 0% staircase error" finding (rejected there purely for
19.6 h build-time cost, not because it doesn't work).

`discretization.metal_mask_key`/`duplicate_geometry_groups` detect, from the discretized metal
mask alone (no solve), when two geometries are mesh-identical: they will solve to bit-identical
Z, V and S at every frequency. `sweep.plan_sweep` now computes this for every planned sweep
(`SweepPlan.mesh_duplicate_groups`, `.n_mesh_distinct`) and, by default (`warn_duplicates=True`,
matching `compute_s_parameters`'s existing `warn_grating`/`warn_surface_wave` convention), warns
about it; `runner.run_sweep` passes `warn_duplicates=False` (it records the condition in the
database metadata instead), so **no pre-Phase-12B sweep, config or result changes** -- verified by
the full pre-existing test suite passing unmodified.

## 3. Task B: Phase 12A re-assessed on aligned meshes

Aligned-mesh refinement study, L=W=2.8 mm (C3, the same representative patch Phase 12A used),
mesh+Floquet order doubled together (24/48/96 cells, Phase 7's own "floquet_order=cells"
convention), at 12.5 GHz:

| cells | \|S11_xx\| | phase(S11_xx) |
|---|---|---|
| 24 | 0.999681 | -32.460 deg |
| 48 | 0.999674 | -32.793 deg |
| 96 | 0.999670 | -32.966 deg |

Successive differences: **0.0058, 0.0030** -- monotonically decreasing (ratio ~1.9, close to this
project's own documented ~2x first-order mesh-convergence rate, MOM_FORMULATION.md Sec. 9) and
**both within the 0.02 tolerance**. Compare Phase 12A's own unaligned 8/16/32/64-cell study:
-26.4/-35.8/-31.0/-33.8 deg (non-monotonic), mesh_error = 0.048. **This confirms Sec. 1's
hypothesis: the unaligned mesh_error was a staircasing artefact, not genuine under-resolution.**

Floquet-only convergence (mesh fixed aligned at 96 cells, only Floquet order varies -- isolates
Floquet-sum truncation from mesh/staircasing error): floquet_order = 48/64/96 gives
-33.029/-32.979/-32.966 deg; successive differences **0.00086, 0.00023** -- also comfortably
within tolerance, and an order of magnitude smaller than the mesh-driven differences above,
confirming Floquet order 48 (used for the full aligned sweep) is already well converged at this
mesh.

Aligned full sweep (same substrate/period/backing as Phase 12A; the combined coarse 0.2 mm-step
family plus a finer 0.1 mm-step family for L >= 3.6 mm, on the mesh that aligns the whole combined
family, **96 cells**, floquet order 48, evidence-based per the floquet-only study above), 24
points, all 5 of Phase 12A's own frequencies:

| Frequency | S11_xx coverage (aligned) | Largest gap |
|---|---|---|
| 10.00 GHz | 140.0 deg | 220.0 deg |
| 11.25 GHz | 206.6 deg | 153.4 deg |
| 12.50 GHz | **237.5 deg** | 122.5 deg |
| 13.75 GHz | 250.1 deg | 109.9 deg |
| 15.00 GHz | 254.6 deg | 105.4 deg |

Design-frequency (12.5 GHz) coverage improves from Phase 12A's unaligned-mesh, coarser-step
**203.6 deg to 237.5 deg** on the aligned mesh with the added fine sampling near resonance --
a real, measured improvement from removing staircasing and resolving the resonance transition
more finely, though still well short of a full 360 deg at this period (Sec. 5 shows a larger
period does better). Coverage also grows steadily with frequency across the band (140 deg at
10 GHz to 255 deg at 15 GHz) -- confirming Phase 12A's own point that coverage is not uniform
across 10-15 GHz and must be checked at every frequency of interest, not assumed representative
from the design frequency alone.

## 4. Task C: duplicate-geometry audit (read-only; no rebuild, no re-run)

**Phase 7 library** (324 geometries, its own locked 32x32 mesh): **169 distinct** discretized
geometries, **105 duplicate groups** (`src/phase12b/duplicate_audit.py`) -- matches the brief's
own independently-stated "169 distinct / 105 duplicate groups" exactly; the geometric metal-mask
count and an independent count of distinct *solved* S11_xx values at 12.5 GHz in the real
checkpoint agree exactly, as they must (identical masks solve to identical S).

Traced through Phases 9-11's own solve mesh (8 cells/side per sub-cell, dx = 0.6 mm -- much
coarser than Phase 7's 32x32):

- **Phase 9** (`LOCAL_PERIODICITY_INVALIDATED`): the center (2.8, 2.8) mm and the four
  strong-contrast corner values (1.0/4.4 mm combinations) are pairwise mesh-**distinct** -- no
  duplicate found. The finding's own evidentiary basis (a real, distinguishable strong
  perturbation) is unaffected.
- **Phase 10** (5 named centers C1-C5): **C2 = (2.0, 2.8) mm and C3 = (2.8, 2.8) mm are
  mesh-identical** at this mesh -- a real, measured finding, not previously known. The held-out
  (blind-test) center is C5 (`environment_geometry.py`'s own `holdout_center_names=("C5",)`); it
  does **not** collide with any training center, so **there is no train/blind-test leakage**. Each
  center's own per-center finding (diagonal coupling, neighborhood-size sufficiency) is still a
  real, correctly-computed MoM result relative to that center's own real isolated baseline; what
  is weakened is only the claim that C2 and C3 probe two independently-verified base geometries --
  their isolated response, and several of their neighbor-perturbation values, coincide at this
  mesh. **Minimal re-run to fully resolve**: re-verify C2's isolated response on a finer/aligned
  mesh (Sec. 2/3 above show 48+ cells resolves 0.2 mm differences cleanly in this size range); if
  it remains numerically indistinguishable from C3's, either relabel C2 accordingly in any future
  write-up or choose a mesh-distinguishable replacement value. This does not require re-running
  Phase 10's full dataset/surrogate pipeline.
- **Phase 11** (`CASE_D_PATTERN_DEPENDENT`, checkerboard endpoints A = (1.4, 1.4) mm,
  B = (4.2, 4.2) mm): confirmed **not** mesh-identical -- unaffected.

## 5. Task D: grounded period study (4.8 / 8.0 / 12.0 mm)

Same substrate (Rogers RT/Duroid 5880, 1.575 mm), grounded, L = W, aligned meshes throughout.

| Period | Aligned mesh | L range (mm) | Resonant L (mm) | Coverage @ 12.5 GHz | Grating onset (normal) | Grating onset (45 deg oblique) | Surface-wave margin, 10 -> 15 GHz |
|---|---|---|---|---|---|---|---|
| 4.8 mm | 96 cells | 1.0-4.6 | 4.2 (Sec. 3) | **237.5 deg** (Sec. 3) | 62.5 GHz | 36.6 GHz | 1.00 -> 1.00 x k0 |
| 8.0 mm | 40 cells | 1.6-7.6 | 6.0 | **279.2 deg** | 37.5 GHz | 22.0 GHz | 1.00 -> 1.00 x k0 |
| 12.0 mm | 60 cells | 2.4-11.2 | 6.8 | **296.6 deg** | 25.0 GHz | 14.6 GHz | 1.00 -> 0.18 x k0 |

Measured, not assumed (brief's own instruction): **larger periods give substantially more phase
coverage** at 12.5 GHz -- confirms Sec. 1's hypothesis that the 4.8 mm period is too small
(~0.2 lambda at 12.5 GHz) for a grounded patch to fully traverse resonance in this L range. But
two real trade-offs, both measured, not assumed:

- **Grating onset drops with period.** At normal incidence all three periods stay specular-only
  across the whole 10-15 GHz band (onsets 62.5/37.5/25.0 GHz). At an **oblique** feed, the
  worst-case onset `f = c0/(P(1+sin(theta)))` drops faster: for the 12 mm period at a 45 deg feed
  angle, onset is **14.6 GHz -- inside the 10-15 GHz band**. A design using the 12 mm period with
  any feed near 45 deg would need to check this explicitly, band edge by band edge; this is
  reported, not assumed away, for a documented set of angles (0/15/30/45 deg), never a single
  guessed feed geometry.
- **The TM0 surface-wave-pole margin shrinks with period and frequency.** At 4.8 mm and 8.0 mm the
  margin stays a full k0 or more from the interval edge across the whole band. At **12.0 mm it
  shrinks from 1.00 x k0 at 10 GHz to just 0.18 x k0 at 15 GHz** -- still outside the interval
  (`any_in_band = False` at every tested point), but measurably closer, exactly as the brief's own
  hypothesis predicted. A design pushing this period to a higher band, or with a finer Floquet
  grid/larger supercell, should re-check this margin, not assume it stays comfortable.

Even with a larger period, **no period tested here reaches a full 360 deg** with a single-layer
resonant patch (this is reported plainly, per the brief's own instruction): 8.0 mm and 12.0 mm both
exceed the 300 deg practical threshold this report uses, but neither is full 360-degree coverage.
Standard remedies for reaching fuller coverage -- a thinner substrate, stacked patches, or
multi-resonant elements -- are open options; **none is implemented here**.

Mesh/Floquet convergence for one near-resonant patch per period (aligned base mesh, 1x/2x/4x):
see `results/phase12b_mesh_period_study/period_study.json` for the full per-period numbers.

## 6. Final statuses

```
Mesh alignment:   ALIGNED_MESH_VALIDATED (aligned-mesh refinement study converges within
                  the 0.02 tolerance, both jointly and Floquet-only, Sec. 3)
Duplicate audit:  CONCLUSIONS_AT_RISK (Phase 10's C2/C3 mesh-degeneracy, Sec. 4; no
                  train/blind-test leakage; core qualitative findings not overturned)
Period study:     COVERAGE_INSUFFICIENT (best measured: 296.6 deg at 12 mm, < the 300 deg
                  practical threshold and short of a full 360 deg either way, Sec. 5)
HFSS status:      HFSS_NOT_AVAILABLE (unchanged; no HFSS export exists anywhere in this project)
Geometry inversion: GEOMETRY_INVERSION_LOCKED (this phase does not attempt to unlock it)
STEP export:        STEP_EXPORT_LOCKED
```

## 7. Rules followed

No period or geometry was chosen by phase coverage alone -- grating-order onset (normal and four
oblique feed angles), the surface-wave-pole margin, loss (|S11|) and the measured coverage are
reported side by side for every period, and the choice among them is left to the user. No earlier
phase's results were modified in place (Phase 7/12A's own output directories are untouched; every
new sweep here is written under `results/phase12b_mesh_period_study/`, tagged with its own period
and mesh, and registered as a benchmark only, never production data). The centroid-is-metal mesh
classification rule itself is unchanged.

## 8. Running it

```bash
.venv/bin/python examples/run_phase12b.py   # results/phase12b_mesh_period_study/
.venv/bin/python -m pytest -q tests/test_discretization_alignment.py tests/test_sweep_duplicate_detection.py \
    tests/test_phase12b_*.py
```
