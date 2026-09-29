# Phase 13B brief — fixes from the first live UI review (2026-09-29)

A reviewer drove every page of the running app (`make app`, http://127.0.0.1:8000) in a real browser with real
solves. The backend physics reproduced the recorded numbers (unit cell L=6.0 grounded → −168.74° @ 12.5 GHz;
library sweep matches Phase 12B exactly; Verify N=6 → 0.529 intended-order, |pred−meas| = 0.0017). The problems are
in how results are wired and shown. Fix in this order; each item gets a regression test that would have caught it.
Physics/numerics unchanged. Commit per group.

## A. Correctness bugs (wrong numbers or silently wrong results)
1. **Richardson mode displays the wrong value.** Unit cell, grounded L=W=6.0, 12.5 GHz, mode=richardson: backend
   `convergence.phase.value_inf = −174.73°` (correct), but the results table and `record.S11_xx_phase_deg` show
   **−169.69°** (the 48-cell finest value) while labelled `MESH_CONVERGED`. The displayed/recorded S-parameters in
   richardson mode must be the extrapolated values (magnitude and phase), with the error bar shown next to them.
   Also: the reported error bar (0.0076°) comes from a 3-parameter fit to 4 points and is not a credible uncertainty —
   use a defensible estimate (e.g. leave-one-out spread of φ∞, or |φ∞ − finest|/some factor) and document it.
2. **HMSA synthesis runs on the raw field grid, not on unit cells.** Synthetic `.fld` upload (151×151, 2 mm pitch),
   period 8 mm → result says **22801 cells** (one per 2 mm field sample). It must resample the reference field
   (complex interpolation, per the Phase 13 brief) onto the cell-centre grid of the chosen period (300 mm / 8 mm →
   37×37) and synthesize per cell.
3. **Frequency mismatch is silent.** Reference field uploaded at 10.54 GHz, library only has 12.5 GHz, synthesis at
   10.54 → **100 % unreachable**, no explanation. Synthesis at 12.5 GHz with a 10.54 GHz field is accepted without
   warning. Required: (a) error `LIBRARY_FREQUENCY_UNAVAILABLE` (with the library's frequencies) when the requested
   frequency isn't in the library; (b) carry the upload's frequency into step 2 and refuse or loudly label a
   mismatch.
4. **"Reachable" ignores the phase error.** With the 12.5 GHz library, required +3.6°…+28° cells all map to L=2 mm
   (−21.95°, errors 25–50°) and are flagged `reachable: true`, 0 % unreachable. Report per-cell achieved-vs-required
   error, summary stats (mean/RMS/max), and a threshold-based unreachable/`PHASE_COVERAGE_SHORTFALL(x°)` label — same
   for the reflectarray designer, which currently shows errors up to 56° with no coverage label.
5. **HMSA aperture field ignores the synthesized cells when lookups fail** (all `1+0j` in the 100 %-unreachable run).
   Never hand a placeholder aperture to Pattern prediction; refuse with an error instead.
6. **Library "Resonant?" is always false** for the grounded sweep, including L=W=6.0 mm at 12.5 GHz (the measured
   resonance, |S11| minimum, phase through ±180°). Fix the classifier for 1-port (grounded) libraries and show
   coverage / largest gap / resonant L for the sweep.
7. **TM0 sanity check never runs.** `upload_field` returns `"tm0_comparison": null` although eps_r, thickness and
   frequency were supplied. Expected ≈ 224.9 rad/m at 10.54 GHz for RT/Duroid 5880, 1.575 mm.
8. **Export summary says "10 cells"** for a 10×10 design (layout CSV and STEP correctly contain 100 patches). Show
   nx×ny and total.
9. **peak_memory_kb is off by ~1000×** (1.4e9 "kB" for a 2.5k-unknown solve): macOS `ru_maxrss` is bytes, Linux is
   kB. Normalize per platform.
10. **Verify page**: "(predicted , |…|" — the predicted value is blank; the specular (0,0) order is missing from the
    order table. "N cells (aperture side length)" is wrong wording — it's cells per 2π gradient period.

## B. Workflow gaps
1. **Plots (largest gap).** Server-rendered matplotlib PNG per job (offline, no JS lib): unit-cell S11 mag/phase vs f;
   library phase-vs-L curves per frequency; reflectarray required/achieved/error maps and layout; far-field cuts
   (E/H plane) + 2-D pattern; HMSA reference phase, hologram/control map, layout; HFSS overlay vs MoM; Verify order
   powers. Replace the 100-row tables with plots + a downloadable CSV.
2. **HMSA → Export** not wired. Carry per-cell geometry from HMSA results in the same shape the exporter uses.
3. **HMSA reference-wave rescaling to β̄ (Phase 13 brief) is missing** — no control on the page. Add it (default on),
   show the accumulated phase error of using the bare-slab field.
4. **Formulation B uses a normal-incidence S11-phase library for a surface-wave hologram.** That is a reflectarray
   control variable, not a surface-wave one. Label B `EXPERIMENTAL — reflectarray-style control variable` on the HMSA
   page and make A (impedance/dispersion) the recommended path; A must also run phase-only with uniform amplitude +
   `REFERENCE_AMPLITUDE_ASSUMED_UNIFORM` (the brief allows it; the page currently says A needs a complex field).
5. **HFSS compare** hard-fails unless the file contains every reference frequency (`no exact frequency match for
   10.25 GHz`). Compare on the overlapping frequencies (report which), accept any sweep/step. Update the stale text
   ("Only G1 has real HFSS data") — G4 (zero-thickness sheet) and G6 are on disk and ingested. Allow comparing
   against an arbitrary unit cell from the Unit cell page, not only G1–G6.
6. **Verify** only runs the canned 12D gradient. Let it take a design job and a row/column cut (or a small
   region) from it.
7. **Inputs that should come from the chosen library/design are re-typed** (reflectarray period, export
   period/substrate). Pre-fill from the library/design and refuse mismatches.
8. **Defaults**: Unit cell defaults to *ungrounded* and *fast* mesh; the design point is grounded and the brief made
   Richardson the app default. Library sweep default mesh (24 cells @ 8 mm) is not aligned with its default L values
   and is ~9° off at resonance — default to an aligned mesh (Phase 12B helper) and warn when unaligned. Library sweep
   only offers full L×W grids — add an "L = W (square)" option.
9. **Diagnostic bundle download** button (infrastructure exists) on every result.
10. **Array export README** gives unit-cell Floquet-port advice; a finite array in HFSS needs a different setup
    (radiation boundary/PML, horn or ideal feed, far-field setup). Split into unit-cell vs array guidance.

## C. Presentation
- Show diagnostics (mesh, unknowns, time, margins) in a collapsible panel instead of only in raw JSON.
- The grounded-lossy "energy_conservation_residual" (0.0086) is absorption, not an error — label it "absorbed power".
- Browser pane width ~800 px: layout works; nav wraps to two rows — acceptable.

## Done when
All A items fixed with tests, B1–B5 done, the rest done or explicitly deferred in `docs/phases/PHASE13_WEBAPP.md`;
fast suite green; a scripted walk-through (API-level) of: grounded unit cell (richardson) → square library sweep →
reflectarray design → pattern → export; `.fld` upload → HMSA (A with dispersion CSV, B labelled) → pattern → export.
