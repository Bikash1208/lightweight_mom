# Phase 14 brief — built-in eigenmode (dispersion) library + a method choice in the app

Read `CLAUDE.md`, `docs/theory/CONVENTIONS.md`, `docs/theory/GREEN_FUNCTION.md`, `docs/theory/MOM_FORMULATION.md`,
`docs/briefs/PHASE13D_POLISH_AND_PHYSICS.md` (§1.1) and `src/dispersion/` (analytic TM0, HFSS table import, MoM
scaffold) first. Existing physics untouched; this adds a capability. Commit per part; fast suite green each time.

## Why
In the live walk-through (4 mm cells, 10.54 GHz, RT/duroid 5880 1.575 mm) the reflection-phase method reaches only
≈ 33° of phase, so 4006/5476 cells clip to the smallest (0.6 mm) or largest (3.6 mm) patch — a near-binary hologram
(`docs/screenshots/patch_layout_example.png`). For a surface-wave (monopole-fed) hologram the control variable is the
cell's **surface-wave wavenumber β / surface reactance X_s**, which small patches do change strongly. The user wants
the library to compute this itself (an eigenmode sweep, like HFSS's), **and** wants the app to offer both methods.

## Part A — MoM eigenmode solver (`src/dispersion/mom_eigenmode.py`, replacing the scaffold)
- Problem: grounded periodic cell (period P, square patch L, substrate εr/h), no incident field. Find real β > k0
  (propagation along x first; ky0 = 0) such that the MoM impedance matrix Z(kx0 = β, ky0 = 0) is singular.
  Reuse the existing Floquet assembly (`floquet_dyads(lattice, kx0, ky0, order)` / `mom_matrix`) with kx0 = β; the
  (0,0) harmonic is then evanescent in air — the kz branch convention already handles it. Don't route through the
  plane-wave excitation (it rejects |kt| ≥ k0).
- Root finding: scan β on (k0, π/P) (first Brillouin zone; stop and flag if the mode reaches the zone edge — a
  stop-band) using the normalized smallest singular value σ_min/σ_max of Z (or a deflated determinant), then refine
  (Brent/golden on log σ). Watch the **bare-slab TM0 pole**: G itself diverges at kt = β_TM0 for the (0,0) term;
  the loaded mode sits above it — handle the pole explicitly (exclude a small bracket around β_TM0, or deflate by
  (β − β_TM0)), document the method, and verify each root by checking the null-vector current is a genuine patch
  current (not a pole artefact). Solve lossless (tanδ = 0) for the table; optionally report the complex root.
- Outputs per patch size: β, n_eff = β/k0, X_s = η0·√((β/k0)² − 1) (TM surface), root quality (σ_min), mesh used,
  stop-band flag. Also the isotropy check: β along x vs along the diagonal for one size (square patches should be
  nearly isotropic; report the spread).
- Mesh: aligned meshes (Phase 12B helper); convergence check on one size (Richardson as in Part 0).
- **Validation (tests):**
  1. L → 0 (or no metal) reproduces the analytic TM0: **224.96 rad/m at 10.54 GHz** for RT/duroid 5880, 1.575 mm —
     the same value measured from the user's HFSS monopole field (`inputs/reference_fields/efl_phase_EzZ_example.fld`).
  2. β increases monotonically with L up to the stop-band; physically sensible range.
  3. If the user supplies HFSS eigenmode spot values (a 2.0 mm and a 3.6 mm patch at 4 mm pitch, 10.54 GHz), compare
     and record; until then, the table is labelled `DISPERSION_FROM_MOM` + "validated against the bare substrate only".
- Library builder: `dispersion sweep` for a board / pitch / frequency / list of sizes (default: 12–16 sizes from
  0.1·P to 0.95·P, denser near the stop-band), cached by config hash, CSV/JSON in the same schema the HFSS-table
  import produces (`src/dispersion/hfss_table.py`), so both sources are interchangeable.

## Part B — app: let the user choose the method
- **Step 2 → "Design method"** (segmented control, plain words):
  - **Surface-wave (eigenmode)** — *Recommended for monopole-fed boards.* Sub-choice **"Unit-cell data"**:
    *Compute automatically* (Part A, default) / *Upload HFSS eigenmode table (CSV)*.
  - **Reflection phase** — *Simpler; limited range with small cells.*
- As soon as board + cell size + frequency are set, show a small preview: β (or X_s) vs patch size for the eigenmode
  method, or phase vs patch size for the reflection method, with the **usable range** stated in one line
  ("Surface-wave index 1.02 – 1.31 across 0.4 – 3.8 mm" / "Phase range 33° — too small for a clean hologram").
  Computing the table shows progress and a time estimate; results are cached.
- For the eigenmode method, set the hologram's mean reactance X̄ and modulation depth M **automatically** from the
  table (centre of the usable range, largest symmetric depth), editable under Advanced — the current fixed
  "100 Ω / 0.3" defaults go away. Rescale the reference wave to β̄(X̄) (already implemented; keep it on).
- Step 3 states which method and data source were used, and the Phase 13D outputs (patch-layout PNG in copper/grey,
  pattern, stats) work for both. Compare the two methods on the walk-through case in the phase report: unreachable
  cells, patch-size histogram, pattern/SLL.
- **More tools → Unit-cell library**: add a "Dispersion (eigenmode) sweep" mode next to the reflection-phase sweep,
  with the β/X_s plot and CSV download.
- HFSS table upload: document the CSV format in the UI help (columns: patch size (mm), β (rad/m) or X_s (Ω) or n_eff,
  frequency) and add short HFSS eigenmode setup notes to the HFSS notes export (master/slave with phase shift,
  eigenmode solver, adequate air height/PML) — the user may still run spot checks.

## Part C — pattern model for the impedance hologram
Check that the aperture/pattern model used for the impedance method is physically consistent (leaky-wave radiation
from the modulated surface, not the reflection-phase aperture); document the approximation; include the
element/obliquity factor from 13D §1.2 so the SLL is correct.

## Done when
Part A validated (test 1 exact to a stated tolerance, tests 2 and convergence pass); the app offers both methods with
automatic table building; the walk-through case (1 beam θ=20°, 4 mm cells, 10.54 GHz) run with the eigenmode method
shows a smooth spread of patch sizes and far fewer unreachable cells than 4006/5476 — report the actual numbers either
way; screenshots updated; docs (`PHASE14_EIGENMODE.md`, USER_GUIDE) written.

## Part D — performance (user concern: "it will take a lot of time")
The sweep must not make the wizard feel slow. Measure first, then apply what's needed:
- **Cache forever** by (substrate, thickness, pitch, frequency, sizes, mesh) — a board is computed once, reused by
  every later design; show "Using saved unit-cell data" when it's a hit.
- **Start early**: begin computing in the background as soon as Step 2 has board + cell size + frequency (don't wait
  for "Design"); Step 3 only waits for whatever is left, with a real progress bar and time estimate.
- **Parallel** across patch sizes (process pool, one size per worker; report cores used).
- **Coarse → fine**: first a quick pass (fewer sizes, coarser mesh, few seconds) so the preview and a first design
  appear fast; refine in the background and let the user re-run with the refined table ("Refined data ready — update").
- **Warm root search**: use each size's β as the starting bracket for the next size (β(L) is smooth), instead of a
  full scan per size.
- **Frequency interpolation**: optionally compute a small frequency grid once per board and interpolate β(L, f)
  for nearby frequencies, flagged as interpolated.
- **Presets**: ship pre-computed tables for the default board (RT/duroid 5880, 1.575 mm) at common pitches
  (3, 4, 5, 6 mm) and a frequency grid covering 8–15 GHz, generated by a script in `examples/` and stored small
  (CSV/JSON).
- Report measured times in `PHASE14_EIGENMODE.md`: per root, per size, full table (cold / cached / preset), on the
  user's Mac. Target: first design on screen in < 30 s for the default board; cold custom board < 2–3 min.
