# Sweeps, discretization validation and the unit-cell database (Phase 5)

Phase 5 wraps the **unchanged** solver (`sparameters.compute_s_parameters`, Phases 2–4) in a sweep driver and stores
the results. No electromagnetic formulation is added: one database record is exactly one call of the Phase 4
solver for one unit cell at one frequency. Everything below is bookkeeping, validation, storage and the
circular-polarization change of basis. Conventions of the S-parameters themselves: [S_PARAMETERS.md](S_PARAMETERS.md),
[CONVENTIONS.md](CONVENTIONS.md). Time factor `exp(+jωt)` throughout.

```
config.yaml + sweep: → plan_sweep (expand + validate every cell) → run_sweep (solve, checkpoint) → Database
                                                                                 ├─ database.json (complete)
                                                                                 └─ database.csv  (flat table)
```

Modules: `sweep.py` (definition, grids, cells), `discretization.py` (mesh/Floquet validation), `runner.py` (driver,
checkpoint), `database.py` (data model, JSON/CSV), `polarization.py` (circular basis, handedness), `phase.py`
(wrap/unwrap). Example: `examples/sweep_config.yaml`, `examples/run_sweep.py`.

## 1. Frequency sweep

`sweep.frequency: {start, stop, step}` (units allowed: `8 GHz`).

    f_i = f_start + i·step,   i = 0 … N−1,   N = (f_stop − f_start)/step + 1

Both end points are included and the last point is set exactly to `f_stop`. `(f_stop − f_start)` must be an integer
multiple of `step` (relative tolerance 1e-9 of the step); otherwise the sweep is **rejected**, so a requested end
point is never silently dropped or overshot. `f_i` is computed from the index, never accumulated, so the grid does not
drift. `start == stop` is a single frequency. Without a `frequency` key the sweep uses the single
`excitation.frequency`. Every frequency is an independent solve: nothing is shared between frequencies (verified:
a sub-grid, a different step, or single-frequency runs in reverse order give bit-identical values).

## 2. Geometry-parameter sweep

Sweepable parameters, all fields of `Geometry`: `length`, `width`, `center_x`, `center_y` (lengths, unit suffixes
allowed) and `rotation_deg` (degrees, counter-clockwise). Internal/database names carry the unit:
`length_m, width_m, center_x_m, center_y_m, rotation_deg`.

* Each axis is a `values:` list or a `start/stop/step` range (same end-point rule as the frequency grid); values must
  be finite and distinct.
* Axes combine as a **Cartesian product, first listed axis slowest**; cell indices follow that order.
* Not sweepable (fixed per database): the lattice periods, geometry `type`, substrate, metal, mesh, reference planes.
  A different lattice or stack is a different database. Non-swept parameters take the base `geometry:` values.
* **Validation happens for every cell before the first solve.** Each cell must satisfy the Phase 1 geometry rules
  (feature inside the cell — touching the edge is allowed; `length, width > 0`; `dipole`/`slot` need `length ≥ width`;
  metal must remain) and the mesh must yield at least one rooftop unknown (a feature too small for the mesh would
  otherwise be solved as a *bare slab*, silently wrong). `on_invalid: error` (default) stops with the list of offending
  cells and parameter values; `on_invalid: skip` records them as `status: invalid` with the reason (JSON only) and solves
  the rest.
* `rotation_deg` is converted to radians only when the `Geometry` is built (`units.from_deg`), from the stored degree
  value, so a cell rebuilt from a database is bit-identical to the one solved. The unswept base rotation passes through
  one `to_deg`/`from_deg` round trip (≤ 1 ulp).

`cell_id` = first 16 hex digits of SHA-256 over the canonical JSON of *(context, parameters)*; the context holds the
lattice, substrate (including name), metal, geometry type, mesh, reference planes and the code version
(`src.__version__`, bumped whenever a change alters the numbers). The frequency is not part of it. Two cells with any
different input have different ids; the id is independent of dict order, platform and `PYTHONHASHSEED`.

## 3. Automatic mesh / Floquet validation

The solver keeps the discretization the user configured; validation **reports** how adequate it is and never changes it
(no adaptive refinement). Thresholds are policy set in `sweep.validation` (defaults in `discretization.py`, recorded in
the database metadata); they are heuristics of this project, not physical constants.

**A priori checks** (every cell and frequency, no extra solve):

| quantity | definition | pass criterion (default) |
|---|---|---|
| `n_unknowns` | rooftop unknowns | > 0 (else the cell is invalid) |
| `metal_area_error` | `|A_mesh − A_true| / A_true`, `A_mesh = (#metal cells)·dx·dy` from the centroid rule | ≤ `max_area_error` (0.05); 0 when the pattern edges lie on mesh lines |
| `cells_per_wavelength` | `λ_sub / max(dx, dy)`, `λ_sub = λ0 / Re√ε_c` | ≥ `min_cells_per_wavelength` (10) |
| `specular_only` | `k0 < 2π/max(Px, Py)`, i.e. `f < c0/max(Px, Py)` | true: no grating order propagates |

`specular_only = false` means the (0,0) S-matrix is still passive but omits power radiated into higher orders (checked:
above onset a lossless cell has `I − SᴴS ⪰ 0` with a strictly positive eigenvalue). Each record stores
`checks_passed` = all of the above. `validation.strict: true` turns a failed a priori check into an error before any solve.

**Refinement study** (opt-in, extra solves; `validation.refinement.enabled`). With integer factor `q ≥ 2`, base mesh
`(nx, ny)`, Floquet order `M`, at each checked frequency (`frequencies: max` = the highest sweep frequency, where the
mesh is coarsest relative to the wavelength; `all` = every frequency):

    A = solve(nx, ny, qM)        B = solve(q·nx, q·ny, qM)
    floquet_error = max|S_A − S_base|        mesh_error = max|S_B − S_A|

(max over the 16 complex entries; absolute; `|S| ≤ 1`). These are **successive differences, not error bars**:
`passed` means "the result moved by less than `tolerance` when the discretization was refined". For a first-order
convergent quantity the error of the base mesh is larger than the successive difference (measured successive-difference
ratio ≈ 1.9 in Phase 3, i.e. the base error is roughly twice `mesh_error`); the Floquet error converges about as `M⁻²`
(measured). The study is stored per cell; the tolerance is not part of its cache key (changing it re-labels `passed`
without re-solving).

Analytic anchor: a cell that is all metal has `S11 = −I`, `S21 = 0` for every mesh and order, so both errors vanish
to rounding (tested). Example (`examples/sweep_config.yaml`, 12 GHz): `dFloquet ≤ 5·10⁻⁴`, `dMesh` 0.0014 … 0.0146,
growing with the metal fill.

## 4. Resume / restart

With `out_dir`, every finished result is appended to `checkpoint.jsonl` (one JSON object per line, flushed and
`fsync`ed). An entry's key is the hash of *everything that determines its value*:

    record key     = hash(cell_id, frequency, min_cells_per_wavelength, max_area_error)
    refinement key = hash(cell_id, refinement frequencies, refinement factor)

On (re)start each planned result whose key is present is reused verbatim; only the missing ones are computed. So:

* an interrupted run continues where it stopped, after `Ctrl-C` **and after a hard kill** (tested with `os._exit`);
* a sweep can be **extended** (more frequencies or parameter values): only the new results are computed;
* changing mesh, materials, geometry, reference planes or the code version changes the key: stale results are never
  reused (they stay in the file, unused). Changing only the refinement tolerance re-labels without re-solving;
* a torn last line or any unparsable/structurally invalid line is dropped (counted in `SweepResult.dropped`) and
  recomputed; the file is rewritten atomically from its valid lines on every start;
* a file at the checkpoint path that is not a version-1 checkpoint is **never overwritten** (`CheckpointError`),
  even with `resume=False`; `resume=False` otherwise discards the stored results;
* results reused from the checkpoint and freshly computed results pass through the same JSON round trip, so they are
  bit-identical; the exported database does not depend on interruption history.

Not guaranteed: durability against power loss beyond what `fsync` of the file gives (the directory entry is not
synced), and two processes writing the same `out_dir` concurrently (unsupported).

## 5. Determinism

The database carries no timestamps, host names or absolute paths. Repeated runs — in one process, in separate
processes, with different `PYTHONHASHSEED`, interrupted-and-resumed, or in memory versus on disk — produce
**byte-identical** `database.json`, `database.csv` and `checkpoint.jsonl` (tested on this platform). The only
platform-dependent part is the floating-point result of the linear algebra library: on another machine or BLAS
thread count expect differences at the 1e-13 level, not bit identity. `environment` (Python and numpy versions) is
provenance and is excluded from database equality.

## 6. Database formats

`Database` = metadata + cells + records. A **record** is the 4×4 complex matrix `[[S11, S12],[S21, S22]]` over
(port 1 x, port 1 y, port 2 x, port 2 y) in the *linear* basis at the stored reference planes, plus diagnostics. Entry
names `S<out><in>_<received><incident>`: `S21_yx` = y received at port 2 for x incident at port 1.

**Diagnostics per record:** `n_unknowns`; `incoming_residual` (radiation-condition check of Phase 4, ≈1e-16);
`reciprocity_error = max|S − Sᵀ| / max|S|`; `passivity_min_eig` (smallest eigenvalue of `I − SᴴS`, ≥ 0 for a passive
cell); `cells_per_wavelength`; `specular_only`; `checks_passed`.

**JSON** (`database.json`, complete): `format`, `version` (1), `metadata` (conventions text, context = lattice /
substrate / metal / mesh / planes, frequency grid, parameter axes, validation policy), `environment`, and per cell
`index, cell_id, status, reason, parameters, validation, records`. Each linear entry is `{re, im, mag, phase_deg}`;
`re`/`im` are the value (Python `repr` floats ⇒ **bit-exact round trip**), `mag`/`phase_deg` are derived and
**checked on import**, as are the 16 circular entries (`S_circular`, `mag`/`phase_deg`) and the `cell_id` of every cell
(recomputed from the parameters and context: an edited file is rejected). Non-finite numbers cannot be written or read.

**CSV** (`database.csv`, one row per record, skipped cells absent): columns `cell_id, cell_index, frequency_index,
frequency_hz`, the five parameters, the seven diagnostics, then `<name>_mag`, `<name>_phase_deg` for the 16 linear and
the 16 circular entries (linear and circular blocks each in family order S11, S12, S21, S22). Values are written with
`repr`, so the columns round-trip exactly as numbers; the complex value rebuilt from `mag·exp(j·phase)` differs from
the JSON value by a few ulp (tested < 1e-14). Metadata is not in the CSV; `Database.from_csv(path, metadata=...)`
accepts it. Import validates the header, field counts, numbers, non-negative magnitudes and the circular columns.

Units in files: hertz, metres, degrees (phase, `rotation_deg`), magnitudes linear. Phase is the principal value in
(−180, 180] (Sec. 8).

**Access:** `db.record(cell, freq_index)`, `Record.component("S21_LR")`, `db.series(cell, name)`,
`db.phase_series_deg(cell, name, unwrap=)`, `db.s_array("linear"|"circular")` (cells × frequencies × 4 × 4, NaN where
missing), `db.parameter_array(name)`, `db.unit_cell(cell)` (rebuilds the `UnitCell` from the metadata; recomputing
the stored S from it reproduces it to < 1e-12, tested), `db.failed_validation()`.

### 6A. One-port (grounded) records (Phase 12A)

A grounded stack (`backing="pec"`, no port 2, S_PARAMETERS.md Sec. 5A) has no 4×4 matrix to store. `Record.port_topology`
(`"two_port"` default, every pre-Phase-12A record; `"one_port"` for a grounded solve) selects which shape `Record.s` must
have: `(4, 4)` as always, or `(2, 2)` — S11 only. The same 7 diagnostics are still computed (from the 2×2 block instead of
the 4×4 one: `reciprocity_error = max|S11 − S11ᵀ|/max|S11|`, `passivity_min_eig` from `I₂ − S11ᴴS11`), so the on-disk
diagnostic schema is unchanged either way. `Record.component()` raises for any family other than `S11` of a one-port
record (undefined, not zero, matching `SParameters`).

**JSON:** each record additionally carries `"port_topology"`; a one-port record's `"S"`/`"S_circular"` only have an
`"S11"` key. This is purely additive — a file written before Phase 12A has no `"port_topology"` key at all, which is read
back as `"two_port"` (the only topology that existed then). **`FORMAT_VERSION` is deliberately not bumped**: it is compared
with strict equality on load, so bumping it would make every pre-Phase-12A file fail to load (the opposite of what a
version bump is for) for a change that is fully backward-*and*-forward compatible without one.

**CSV does not support one-port records**: the fixed 16-linear + 16-circular column layout assumes a two-port S-matrix.
`Database.to_csv()` raises `DatabaseError` if any record is one-port; use `to_json()` for a grounded database (`sweep.py`'s
`run_sweep` does this automatically — a grounded sweep writes `database.json` only, no `database.csv`).

`sweep.sweep_context()` records `"backing"` (`"air"`/`"pec"`) alongside lattice/substrate/metal/mesh/reference-planes, so a
grounded and an ungrounded sweep of the same geometry get different `cell_id`s and are never confused as "the same cell"
by the checkpoint or database machinery; absence (a context written before Phase 12A) means `"air"`.

For a grounded cell, `CellEntry.validation["surface_wave"]` additionally records the Phase 12A TM0 surface-wave-pole-proximity
result (`sparameters.surface_wave_proximity`, at the sweep's own highest frequency) -- `None` for an ungrounded cell. This is
per-*cell* (`validation` is already free-form/optional, GREEN_FUNCTION.md Sec. 12A), not a new per-*record* diagnostic key,
so `Record.diagnostics`'s strict 7-key schema is unchanged either way.

## 7. Circular polarization

**Basis (fixed by this project, not taken from any software package):**

    e_R = (e_x + j e_y)/√2        e_L = (e_x − j e_y)/√2,      T = [e_R  e_L]  (unitary)

These are global x–y vectors, identical at both ports and for incident and scattered waves: **R/L is a label of the
transverse vector, not of the propagation direction.** Amplitudes: `(a_R, a_L) = Tᴴ (E_x, E_y)`; S-matrices:
`S_circ = Tᴴ S_lin T`, entries `[received][incident]` (`S21_LR` = L component at port 2 for R incident at port 1);
the 4×4 two-port matrix is transformed with `blockdiag(T, T)`.

**Handedness under `exp(+jωt)`.** The physical field at a fixed point is `E(t) = Re{e·e^{jωt}}`. For `e_R`:
`E(t) = (cos ωt, −sin ωt)/√2` — along +x at t = 0 and along −y a quarter period later. The average rotation rate about
+z is `ω·Im(E_x E_y*)`, which is `−ω/2` for `e_R`. Therefore

* `e_R` rotates **clockwise seen from +z** (from x toward −y); `e_L` counter-clockwise.

Handedness of a *wave* uses the right-hand rule (IEEE): with the right thumb along the propagation direction, the
fingers curl in the direction E rotates in time ⇒ right-handed (RHCP). (The optics convention names the same waves the
other way round and is **not** used.) With `exp(+jωt)`, a wave `exp(−jkz)` travels toward +z and `exp(+jkz)` toward −z.

| wave | travels | `e_R` is | `e_L` is |
|---|---|---|---|
| incident at port 1 (top); transmitted into port 2 (bottom) | −z | **RHCP** | LHCP |
| reflected at port 1; incident at port 2; transmitted into port 1 | +z | LHCP | **RHCP** |

Consequences (`polarization.component_handedness`, verified against a brute-force evaluation of the physical field):

| entry | incident | outgoing | reading |
|---|---|---|---|
| `S21_RR`, `S21_LL` | RHCP / LHCP | same | transmission **preserves** handedness |
| `S21_LR`, `S21_RL` | RHCP / LHCP | opposite | handedness conversion |
| `S11_RR`, `S11_LL` | RHCP / LHCP | **opposite** | same label, reversed physical handedness |
| `S11_LR`, `S11_RL` | RHCP / LHCP | same | reflection preserving handedness |

A full PEC sheet has `S11 = −I` in the linear basis, hence `S11_RR = S11_LL = −1`, `S11_LR = S11_RL = 0` (the label is
kept, the wave reverses direction, so an RHCP wave returns LHCP) — tested on solver output.

**Co- and cross-polar** are defined by label equality, in each basis separately: linear co = `xx, yy`, cross = `yx, xy`;
circular co = `RR, LL`, cross = `LR, RL`. Total power `Σ|S|²` is basis independent (unitarity); a half-wave plate
`diag(1, −1)` is purely co-polar linearly and purely cross-polar circularly.

**Circular-basis identities** (derived here; each verified on solver output):

* Reciprocity `S12 = S21ᵀ` (linear) becomes, because `conj(e_R) = e_L`,
  `S12_RR = S21_LL`, `S12_LL = S21_RR`, `S12_LR = S21_LR`, `S12_RL = S21_RL`.
* Geometric (Pancharatnam–Berry) phase: rotating the cell counter-clockwise (seen from +z) by ψ, `S' = R(ψ) S R(ψ)ᵀ`,
  gives `RR' = RR`, `LL' = LL`, `LR' = e^{+2jψ} LR`, `RL' = e^{−2jψ} RL`, because `R(ψ) e_R = e^{−jψ} e_R`. Checked with
  random matrices for six angles ψ, and on solver output for ψ = 90° (a 90° rotation of a rectangle on a square mesh is exact).
* Lossless, single propagating order: the circular 4×4 matrix is unitary, and every circular input (port 1) sends
  total power 1 into the four circular outputs.

## 8. Phase wrapping and unwrapping

* Stored phase is `arg(S)` in **degrees on the half-open interval (−180, 180]**; the negative real axis is +180 for
  both signs of a zero imaginary part (numpy would return −180 for `−1 − 0j`). A delay is a negative phase.
* Unwrapping along an ordered axis (frequency, or a swept parameter): `u_0 = p_0`, `u_i = u_{i−1} + wrap(p_i − p_{i−1})`;
  anchored to the first sample, defined only up to a global multiple of 360°.
* It equals the continuous phase **only if the true phase changes by less than 180° between neighbouring samples**,
  and nothing in the data can reveal a violation (an undersampled sweep unwraps cleanly to a wrong answer — demonstrated
  in the tests). `Database.phase_series_deg(..., unwrap=True)` therefore warns when the largest wrapped step exceeds
  90° and refuses a gap in the frequency indices; `phase.max_wrapped_step_deg` gives the number.
* The reference planes dominate the slope: `S21` carries `e^{−jk0 (z_top − z_bottom)}`, about −242°/GHz for planes
  100 mm either side of the cell (tested: the unwrapped change over 8 → 10 GHz is −488.4° vs. −484.1° from the delay alone).
  For unwrapping, keep the planes close to the cell (`examples/sweep_config.yaml` uses 0) or take a fine step; re-referencing
  after the fact is `SParameters.at_planes` in a live solver session (not stored records).

## 9. Config reference (`sweep:` section)

```yaml
sweep:
  frequency: {start: 8 GHz, stop: 12 GHz, step: 0.5 GHz}      # optional
  parameters:                                                # optional; length | width | center_x | center_y | rotation_deg
    length:       {start: 2 mm, stop: 6 mm, step: 1 mm}      #   range, or
    rotation_deg: [0, 45, 90]                                #   list (also {values: [...]})
  on_invalid: error                                          # error | skip
  validation:
    min_cells_per_wavelength: 10
    max_area_error: 0.05
    strict: false
    refinement: {enabled: false, factor: 2, tolerance: 0.02, frequencies: max}   # max | all
```
Unknown keys, unknown parameter names (`period_x`, `type`, `rotation`, …), both or neither of `values` / `start-stop-step`,
a step that does not divide the span, and wrong types are `ConfigError`s that name the offending key.

## 10. Verification (`tests/test_phase.py`, `test_polarization.py`, `test_sweep_definition.py`, `test_discretization.py`, `test_database.py`, `test_sweep_run.py`)

| property | how it is checked |
|---|---|
| frequency sweep consistency | every record = independent single-frequency solve (< 1e-13); values independent of the other frequencies (bit-identical for sub-grids, other steps, reversed order); lossless cell: `SᴴS = I`, `S = Sᵀ` at all 9 frequencies × 4 cells; lossy cell strictly passive; change per step ∝ step |
| geometry parameter validation | every rule in Sec. 2 with the exact message fragment; nothing solved before the error (solver spy); `skip` mode; strict a priori mode; config-level errors |
| database round trip | JSON bit-exact and idempotent byte-for-byte; CSV columns exact, complex to 1e-14; rebuilt unit cells reproduce stored S to 1e-12; tampering (parameters, magnitude, phase, circular column, duplicates, versions, NaN) rejected |
| phase wrap / unwrap | wrap ∈ (−180, 180], idempotent, period 360, ±180 cut; ramps recovered for steps < 180 and shown to alias above; along any axis; real sweep: unwrapped phase independent of sampling and equal to the plane-delay prediction |
| co/cross conversion | closed forms derived by hand; unitary, round trip, norm/singular-value/trace/det preserved; power conservation; wave-plate cases; rotation law (arbitrary ψ, and 90° on solver output); circular reciprocity and unitarity on rotated cells |
| handedness | brute-force right-hand-rule evaluation of `Re{e e^{jωt}}` (and direction of travel from the phase fronts) vs. `handedness()`; PEC reflection reverses handedness on solver output |
| determinism | repeated runs, in-memory vs disk, two `PYTHONHASHSEED`s: byte-identical outputs; no timestamps or paths |
| resume | interrupt at every stage, hard `os._exit`, torn/corrupt lines, foreign file, extension, key sensitivity, refinement caching; result identical to an uninterrupted run |
| mesh/Floquet validation | area error and unknown counts against independent counts; `λ_sub`/mesh formula; grating onset `c0/max(P)` vs. the solver's own warning; refinement = the stated pair of solves; full PEC sheet converged to 1e-10; errors decrease with order (`M⁻²`) and mesh |

Mutation checks (scratch copy, 43 deliberate faults across the new modules — e.g. R/L swapped, handedness direction,
wrap interval, unwrap without wrapping, range end point, axis order, degrees as radians, periods swapped, cache key
ignoring frequency/factor, corrupt lines trusted, foreign file overwritten, CSV phase in radians, derived-field and
`cell_id` checks disabled, wrong S block, refinement pairing, grating onset): all caught. Two were initially missed —
period order (every test cell was square) and the refinement factor in the cache key — and now have tests.

## 11. Limitations and open items

* Normal incidence, (0,0) port only; higher grating orders are flagged, not modelled. Circular *ports* at oblique
  incidence are not defined.
* One frequency-independent mesh per database; no per-cell mesh adaptation; no reuse of `Z` between the four right-hand
  sides of different cells (each record re-assembles `Z`); no parallelism.
* Refinement differences measure convergence, not accuracy against an external solver, which has still not been done
  (magnitude **and** phase, after conjugating if the reference uses `exp(−jωt)`).
* The circular identities and handedness table are derived here from the stated definitions; they were verified
  numerically but not against a published source. Please cross-check the handedness naming against the convention
  your synthesis stage and your measurement/simulation reference use.
* Only PEC metal, one substrate layer, rectangular features; slots as metal-complement currents.
* Not implemented (by request): holographic synthesis, cell selection/matching, interpolation of the database.
