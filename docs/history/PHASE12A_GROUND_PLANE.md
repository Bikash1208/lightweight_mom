# PEC ground plane / grounded-substrate stack (Phase 12A)

Implementation: `src/green_function.py` (`backing="pec"`), `src/geometry.py` (`UnitCell.backing`), `src/config.py`
(`substrate.backing`), `src/excitation.py` (`ReferencePlanes.bottom_air_thickness=None`), `src/sparameters.py`
(one-port `SParameters`, `surface_wave_proximity`), `src/database.py` (`Record.port_topology`), `src/sweep.py`/
`src/runner.py` (one-port sweeps), `src/cad_export.py` (`build_ground_solid`), `src/phase8d_rev/library_entry.py`+
`query.py` (backing tag), `src/phase8c/topology_audit.py` (real per-database ground-plane read). Orchestration/report
only: `src/phase12a/`. Derivation: GREEN_FUNCTION.md Sec. 12A. Conventions: CONVENTIONS.md, S_PARAMETERS.md Sec. 5A,
DATABASE.md Sec. 6A.

## 1. Why

The ungrounded rectangular-patch library reaches only 24.8° (S11xx) / 43.9° (S21xx) phase coverage at 12.5 GHz
(`LIBRARY_REPORT.md`), and every grounded holographic/reflectarray formulation this project investigated (Fong 2010,
Minatti, Pozar 1997, Budhu & Grbic 2021 — Phases 8B/8E) was `ARCHITECTURE_INCOMPATIBLE` for lack of a ground plane. With
a ground, S21 ≡ 0 by construction (no port 2 exists), so the S11-vs-S21 control-variable question that has been open
since Phase 8A is settled by the architecture itself, for this stack.

## 2. What changed, and what didn't

The existing solver (`Z I = V`, Green's function, basis functions, mesh, periodic treatment, excitation, field
calculation, S-parameters, polarization conventions) is **unmodified** for the default (`backing="air"`, or no
`backing` key at all) stack — verified bit-identical (`test_green_function_grounded.py::test_default_backing_*`,
full pre-existing test suite unchanged, see Sec. 6). Every touched module was **extended**, never rewritten:

- `LayeredGreenFunction(substrate, frequency, backing="air")`: a new optional third field. `backing="pec"` swaps the
  bottom-face load reflection coefficient from the natural air value to a short circuit (`Γ_L = −1`) everywhere it is
  used — the equivalent impedance seen by the sheet, the substrate standing wave, and (for free, no separate code path)
  the field below the ground, which is exactly zero because the `(1 + Γ_L)` factor that already multiplies it vanishes
  at `Γ_L = −1`. The top-of-stack reflection Γ (needed for the unpatterned-stack excitation) required its own new,
  independently-derived generalization (its existing formula assumes a symmetric, air-both-sides slab). Full algebra:
  GREEN_FUNCTION.md Sec. 12A.
- `UnitCell.backing` / config `substrate.backing` (`"air"`\|`"pec"`): plumbs the choice from YAML down to the kernel.
  `reference_planes.bottom_air_thickness` becomes **forbidden** (not optional) when grounded — strict config, matching
  the project's existing "unknown keys are errors" discipline.
- `ReferencePlanes.bottom_air_thickness: float | None`: `None` marks "no bottom port"; `z_bottom` raises if read.
- `SParameters.grounded`: `True` for a one-port result. `s12`/`s21`/`s22` are `None`, and every accessor
  (`.component()`, `.matrix()`, `.to_record()`) raises for anything but S11 — **undefined, not silently zero**
  (S_PARAMETERS.md Sec. 5A). `compute_s_parameters` solves 2 right-hand sides instead of 4 for a grounded cell.
- `Record.port_topology` (`"two_port"` default, `"one_port"` for grounded): stores a `(2,2)` S11-only matrix. JSON is
  purely additive (old files load unchanged, defaulting to `"two_port"`); `FORMAT_VERSION` is deliberately **not**
  bumped (it is checked with strict equality, so bumping it would break loading every pre-Phase-12A file — the opposite
  of what a version bump is for, DATABASE.md Sec. 6A). **CSV export is not supported** for one-port records (the fixed
  16-linear+16-circular column layout assumes two ports); `Database.to_csv()` raises a clear `DatabaseError`, and
  `runner.run_sweep` writes `database.json` only (no `.csv`) for a grounded sweep — an explicit, documented restriction,
  not a silent gap.
- `sweep.sweep_context()` records `"backing"`, so a grounded and an ungrounded sweep of the same geometry get different
  `cell_id`s and are never confused as "the same cell".
- `cad_export.build_ground_solid`: a third, separate solid (a full-period slab flush below the substrate, the same CAD
  convenience as the existing metal-thickness slab — no effect on any solved quantity) for `backing="pec"`;
  `build_unit_cell_solids` itself is untouched (still returns exactly `(substrate, metal)`), so every existing caller
  is unaffected. `import_and_summarize` extended to recognize 3 solids.
- `phase8d_rev.LibraryEntry.backing` + `query.UnitCellQuery.backing`: a grounded and an ungrounded library are never
  silently mixed by `search_unit_cells` when a query specifies `backing`.
- `phase8c.topology_audit._ground_plane_test`: now a real per-database read of `context["backing"]` (`PRESENT` for a
  grounded database, `ABSENT` — unchanged conclusion — for the real Phase 7 library, which has none) instead of a
  hardcoded conclusion that predates the capability to express a ground plane at all.
- `discretization.refinement_study`: compares `a.s11 if a.grounded else a.matrix()` — generalized, not duplicated.
- `sparameters.surface_wave_proximity`/`_warn_if_surface_wave_near`: the TM0 surface-wave-pole-proximity check and its
  `UserWarning` (`compute_s_parameters`'s new `warn_surface_wave=True`, same convention as `warn_grating`).
  `runner.run_sweep` records it once per grounded cell (`CellEntry.validation["surface_wave"]`, at the sweep's own
  highest frequency) rather than as a new per-record diagnostic key, so `Record.diagnostics`'s strict, fixed 7-key
  schema — and every file already written under it — is unaffected either way.

**Not touched at all** (audited, no change needed): `mom_matrix.py`, `basis_functions.py`, `solver.py` (they consume
`Z`/`V` generically and never inspect backing); `validation_cases.py`, `validation_export.py`, `hfss_comparison.py`
(process only the existing, ungrounded, already-exported cases); `library_schema.py` (the one built production library
remains ungrounded — `backing` is simply never set there); `phase8b/conventions.py`'s `NO_GROUND_PLANE` constant
(still true of every existing config/library; no longer a hard solver limit, but this phase does not change what it
describes).

## 3. Verification

- **Independent reference** (`tests/reference_layered.py`, no shared derivation with `src/green_function.py`):
  `solve_current_sheet_grounded`/`solve_incident_fields_grounded` (PEC-terminated 3-D plane-wave field matching —
  region 2 removed, `E_t(−h)=0` replaces the region-2 matching equations) agree with the kernel to **~1e-14** across
  normal/oblique incidence, lossy/lossless substrates, several periods and thicknesses
  (`tests/test_green_function_grounded.py`).
- **Analytic bare-slab check**: `Γ = (jη_s tan(k_s h) − η0)/(jη_s tan(k_s h) + η0)` vs the kernel, lossless and lossy,
  10–15 GHz: maximum disagreement **4.6e-16** (solver floating-point noise floor). `|Γ|=1` lossless at every tested
  frequency; `|Γ|<1` lossy (0.9999→0.9998, decreasing with frequency); `h→0` gives `Γ → −1` (measured
  `−1.000000+5.2e-7j` at `h=1e-9`m).
- **Surface-wave-pole proximity**: for the 4.8 mm (Phase 7 library) and 8 mm (original Phase 1-6) periods, 10–15 GHz,
  normal incidence: **no Floquet harmonic falls in, or within 2% of k0 of, the TM0 surface-wave interval**
  `(k0, Re(k_substrate))` in any tested case; the nearest harmonic is a full `k0` or more from either edge.
- **Full regression**: the entire pre-existing test suite passes unmodified (see README.md's test count); the default
  (`backing="air"`) code path is verified bit-identical for `spectral_dyad` and `stack_fields`.

## 4. Physics-sanity sweep (measured, not assumed)

A grounded rectangular-PEC-patch sweep on the **same** Rogers RT/Duroid 5880, 1.575 mm, 4.8 mm period, 32×32/floquet-32
mesh Phase 7 used (`library_schema.py`, unmodified), L = W over the same 18 values (1.0–4.4 mm, 0.2 mm step):

| Frequency | S11xx coverage (grounded) |
|---|---|
| 12.5 GHz | **203.6°** of 360° (largest gap 156.4°, 18 points) |

Side by side with the real, unmodified Phase 7 checkpoint at 12.5 GHz: **24.84°** (S11xx) / **43.85°** (S21xx) coverage,
324 cells — confirms the two numbers this brief itself quotes, freshly recomputed by this phase's own code
(`sweep_study.ungrounded_phase7_coverage`), not copied from the doc. The grounded patch's coverage is dramatically
larger because, unlike the ungrounded scattering patch (no resonance anywhere in the 10–15 GHz band for these sizes,
`LIBRARY_REPORT.md`), a grounded patch **can** approach an actual resonance within this size/frequency range: reflection
phase moves smoothly from −21° (L=1.0mm) through −138° (L=4.2mm) and wraps to +135° (L=4.4mm) — a real, physically
expected reflectarray-element effect, not a numerical artifact (the classic Pozar 1997 variable-patch-size phase curve,
cited, not digitized).

Coverage is still short of 360° with only 18 discrete points spanning, and not densely sampling, the resonance region —
a thicker/thinner substrate, finer sampling near resonance, or a different topology are the next open questions, not
this phase's to answer (its own brief's explicit instruction).

**A genuine, honestly-reported mesh-sensitivity finding**: `discretization.refinement_study` on one representative
grounded patch (C3, 2.8×2.8 mm, this project's own already-used center geometry) at Phase 7's own 32×32/floquet-32 base
mesh measured `floquet_error = 4.3e-4` (well converged) but `mesh_error = 4.8e-2` — larger than the 0.02 tolerance that
mesh choice was validated against for the *ungrounded* sweep. A supplementary raw scan (mesh and Floquet order doubled
together: 8, 16, 32, 64 cells/side) shows non-monotonic phase values (−26.4°, −35.8°, −31.0°, −33.8°) rather than the
smooth, ~2x-per-doubling convergence documented for the ungrounded case (MOM_FORMULATION.md Sec. 9) — grounded, near-
resonant patches are measurably more sensitive to staircase discretization than the scattering patches this exact mesh/
sweep-step combination (`library_schema.py`) was chosen around. This is reported as data, not smoothed over.

## 5. Final status

```
Grounded-stack status:  GROUNDED_STACK_PARTIAL   (mesh/Floquet refinement study above; every other
                                                    check -- bare-slab closed form, surface-wave-pole
                                                    avoidance -- passed within tolerance)
HFSS status:            HFSS_NOT_AVAILABLE        (unchanged; no HFSS export exists anywhere in this project)
Geometry inversion:     GEOMETRY_INVERSION_LOCKED (this phase does not attempt to unlock it)
STEP export:            STEP_EXPORT_LOCKED
```

`GROUNDED_STACK_PARTIAL`, not `GROUNDED_STACK_VALIDATED`, is itself the honest result: the kernel physics is verified to
the solver's own noise floor by an independent solver, but the *discretization* this phase reused for a fair
Phase-7-comparable sweep does not resolve a near-resonant grounded patch as tightly as it resolves the ungrounded one.
A finer production mesh for grounded designs is the natural next step, not a reason to distrust the underlying algebra.

## 6. Running it

```bash
.venv/bin/python examples/run_phase12a_ground_plane.py   # results/phase12a_ground_plane/
.venv/bin/python -m pytest -q tests/test_green_function_grounded.py tests/test_geometry_grounded.py \
    tests/test_excitation_grounded.py tests/test_config_grounded.py tests/test_sparameters_grounded.py \
    tests/test_database_grounded.py tests/test_runner_grounded.py tests/test_cad_export_grounded.py \
    tests/test_phase8d_rev_grounded.py tests/test_phase12a_*.py
```

Outputs (`results/phase12a_ground_plane/`): `phase12a_report.md`, `phase12a_gate_status.json`, `sweep.{csv,json}`,
`ungrounded_comparison.json`, `analytic_checks.json`, `surface_wave_pole.json`, `convergence.json`, and
`plots/{phase_vs_length,magnitude_vs_length,grounded_vs_ungrounded,convergence}.png`.

## 7. Limitations

- Not a production library: only the L=W diagonal of Phase 7's 18×18 sweep, at one design frequency for the headline
  number (5 frequencies for the frequency-dependence context) — a full 18×18 grounded sweep was not built.
- Mesh convergence for grounded/near-resonant patches needs a finer production mesh than Phase 7's 32×32 choice
  (Sec. 4); no such finer mesh was built or validated here.
- Geometry inversion, STEP export and holographic synthesis remain exactly as locked as before this phase — nothing
  here changes any of Phases 8-11's own findings or gate outcomes.
- HFSS/measured validation is still `NOT_AVAILABLE` — unchanged, now a four-phases-persistent blocker.
