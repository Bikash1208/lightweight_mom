# Changelog

Reconstructed per phase (this project did not keep a changelog before Phase 13; entries below
summarize each phase's own write-up in `docs/history/`/`docs/history/`). No version tag existed
before Phase 13; `0.13.0` (`pyproject.toml`) is the first. See `docs/STATUS.md` for current
validation numbers and open next steps; this file is the per-phase "what changed" record.

## Unreleased - Phase 20 R1 - loads and oblique incidence in the MoM (branch `phase20-ris`, from `phase27-app-tensor-cp` 4b87212; MODEL_ONLY)
- 1: oblique incidence: audit (19 tests; Green function, stack fields, excitation and currents were already right for kt < k0), full 2x2 TE/TM power-wave S-parameters `compute_s_parameters_oblique`, per-frequency grating-order warning, config accepts 0 <= theta < 90.
- 2: lumped RLC loads (`src/loads.py`, `UnitCell.loads`, config `loads`); `open_circuit` / `transverse_conducts=False` (a conducting transverse strip is not the empty gap).
- 3: `CompositeGeometry` (union of metal rectangles minus apertures; one rectangle bit-identical to `rectangle`), mesh-alignment warning, STEP export with `Load_<k>` sheets.
- 4: demo `results/phase20_r1_demo/` and HFSS check package `hfss_exports/phase20_r1_loaded_cell/` (Richardson-extrapolated MoM reference with per-row uncertainty).
- Additive only: theta = 0 without loads is bit-identical to the pre-R1 solver (hex-float reference test) and the Phase 23 snapshot is unchanged.

## Unreleased - Phase 27 - tensor-CP design flow in the app (branch `phase27-app-tensor-cp`)
- 1: the Phase 26 design flow as `src/hologram/tensor_cp_design.py` (single hexagon, boards A / B on the shared lattice, prediction builders); the
  `scripts/research/phase26_*` design scripts are wrappers. Regenerated layouts, cells, results JSON and `predictions.json` are **exactly** the committed ones.
- 2: `src/hologram/cell_family_registry.py`, a cell-family / state-library registry by domain of validity (substrate, period, frequencies with data,
  incidence, data kind, validation label, known biases, fabrication rules); out-of-domain requests give `NOT_AVAILABLE_FROM_CURRENT_DATA`; the schema takes
  discrete N-state libraries (2x2 reflection per state vs f, theta, phi) unchanged -- the slot for RIS (Phase 20 R4).
- 3: "Tensor CP (monopole-fed)" mode: `src/hologram/tensor_cp_service.py`, `/api/tensor_cp/{design,defaults,families}` (background job), page `tensor_cp.html`.
  Boards: single hexagon, four hexagons with one beam per port, four hexagons combined; RHCP/LHCP; pin height; every number source-labelled with the registry's biases.
  The 1-bit HMSA flow is untouched.
- 4: `src/export/tensor_cp_export.py`, `/api/tensor_cp/export` (AEDT script with `DRIVE_PORT`, layout CSV/JSON/DXF/PNG, optional STEP, `predictions.json`, README);
  the single-hexagon export equals the committed `hfss_exports/hex150_tensor_cp_single/` files (DXF compared by entities: ezdxf text is not byte-stable). The committed
  `layout.csv` gained the `cross_slot` column the generator writes (cell data unchanged).
- 5: `src/hologram/tensor_cp_compare.py`, `/api/tensor_cp/compare` (upload the HFSS exports; named states and power budget); with the committed HFSS files it gives 23.03 dBi,
  XPD 23.0 dB, S11 -15.0 dB at 10.54 GHz. Fixed in item 6: the `design_job_id` option of `/compare` (found by the end-to-end test).
- 6: end-to-end test `tests/test_app_phase27_tensor_cp.py`; STATUS / START_HERE updated. Additive hooks only: `march_rays(density_fn=)`, `design_board(beams=, hand=)`.

## Unreleased - Phase 26 - leakage-controlled tensor-CP hexagon
- A: `src/dispersion/leakage_model.py` (first-order and Floquet-balance leakage of the n = -1 harmonic, scalar and tensor sheets),
  `src/hologram/leaky_forward.py`, `src/dispersion/pin_launch.py` (monopole launch fraction, loaded sheet); HFSS data of
  `four_port_hex150_gap0` indexed in `docs/validation/README.md` §6.
- B/B2: `src/dispersion/tensor_sheet.py` (tensor sheet from the 2x2 reflection, Bloch kt > k0 extraction, tensor transverse resonance),
  `slotted_cells.py` (slotted circle/square, cross-slotted circle, `CellRow`), `tensor_floquet.py` (Floquet balance of a modulated tensor
  sheet); solver modules untouched, every pre-existing result bit-identical (Phase 23 snapshot).
- C/D: `src/hologram/tensor_leaky_cp.py` (supercell-table design march, rotating-slot currents, far field) and the exact single-hexagon
  export `hfss_exports/hex150_tensor_cp_single/` (+ V2 graded-pad variant, on hold); HFSS run: 23.03 dBi RHCP, AR 1.3 dB, XPD 23 dB
  (`docs/validation/README.md` §7).
- Decision recorded: board B = combined-beam board (in-phase 1:4); Butler scanning on the 4-hex layout shelved (quantisation lobes, needs ~lambda/2 feed spacing); board B export held until board A HFSS.
- Board B tile-amplitude equalisation study (`phase26_tile_equalise.py`; additive `march_rays(density_fn=)`), 4-hex HFSS comparison script (`phase26_compare_fourhex.py`). Targets not reached; see `results/phase26/README.md` O.
- 4-hex tensor-CP boards: `phase26_fourhex_design.py` (boards A/B on one shared 6.4 mm lattice, edge-power bookkeeping, pickup of the unconverted wave, Butler combinations), `phase26_export_fourhex.py` -> `hfss_exports/hex4_tensor_cp_static/` (board A, DRIVE_PORT option), `results/phase26/fourhex/`; board B scoped only. Model only, not HFSS-validated.
- Studies (model only): substrate study (RT5880 3.175, RO3006, RO3010, RT6010), launcher redesign, aperture-size study, 4-hex scope;
  `scripts/research/phase26_*.py`, results in `results/phase26/`.

## Unreleased - Phase 23 - hexagonal ports and the sidelobe-cancelling optimiser
- A: `PortRegion` may be a convex polygon (`vertices_m`; rectangles unchanged and bit-identical); `layout_preset_hex4` (side or
  flat-to-flat, gap, pointy/flat, shared 6.4 mm lattice, cells outside the hexagons empty); polygon overlap/gap/isolation; DXF
  PORTS layer, outlines in the layout PNG, README table, optional trim-board STEP/AEDT; `hex4` in the UI and `/synthesize_ports`.
- B: `src/pattern/cellwise_pattern.py`: per-cell co-/cross-pol contributions (verified against `_vector_far_field` to 1e-10),
  u-v direction sets denser near the beams, combined multi-port excitation.
- C: `src/hologram/sll_optimizer.py` (greedy q-norm cell flips, optional seeded annealing, gain/pointing/level/XPD constraints,
  per-frequency scenario instances, random port-phase draws, absolute pointing), `sll_board.py` (job result <-> optimiser, held-out
  grid evaluation, x/y-cut SLL as in Phase 18, cells re-solved per frequency), `/api/sidelobes/optimize`, Step 3 "Reduce sidelobes"
  with a Robust option.
- D: `scripts/research/phase23_compare_v2.py`, `results/phase23/`, `hfss_exports/four_port_hex150_*`, `four_port_square243_gap0`.
- New key `two_state.achieved_phase_deg` in coded `/synthesize_ports` results (additive).

## Unreleased - Phase 22 - automatic coded-cell states
- Exact multi-port STEP: `Feed1..FeedN` marker sheets (0.6 x 0.6 mm at z = 0, exact vertices) restored; optional `Monopole1..N` pins
  (checkbox on the Export card; verification then checks 3-D vertices and the pins' z/radius). Feed positions are listed in `feeds.csv`,
  `layout.json`, layer FEEDS of `layout.dxf`, the layout PNG (P1..PN) and the board README; the standard multi-port export gets them too.
- `src/hologram/auto_states.py`: empty + mesh-corrected squares (0.2/0.1 mm MoM + Richardson, interpolated, snapped to 0.05 mm),
  1-/2-/3-bit, disk cache; `auto` two-state source in `/synthesize_ports` and `/synthesize_coded`, `/api/hmsa/auto_states`;
  the near-resonance refusal is skipped for the auto pair unless beta_mode is "state".
- Follow-up: squares up to period - min gap (default 0.2 mm, Advanced setting) via a 0.1 mm-mesh extension past the 6.0 mm coarse-scan limit
  (2-bit level 1 on the reference board: 6.2 mm, -126.0 deg, 5.9 deg error instead of 6.0 mm / 49 deg); the gap is named in UNREACHABLE_LEVEL messages.
- `beta_mode 'paper'` removed (Mishra eq. 1): the API accepts only 'model' and 'state' ('paper' fails with "beta_mode 'paper' was removed; use 'model'"); `paper_surface_wave_beta` and `beta_paper_rad_per_m` are gone; the UI offers Model (default) and Each state's own.
- A level within 10% of the level step counts as reached (UNREACHABLE_LEVEL above that); the Advanced section (minimum gap) stays visible for 2-/3-bit.
- UI: "Automatic" default, results table, "use as uploaded values", fix-L; unit-cell tool patch-plane column; multi-port reachable
  without Step 1; one directivity per port; beta help text.

## Unreleased - Phase 21 - multi-port layouts, coded multi-port, exact board export
- C1: `_cell_centers` floors with a 1e-6 tolerance (300.8 mm / 6.4 mm is 47 cells, not 46); the earlier `four_port_1bit_300mm_*`
  boards were affected. One test fixture that relied on the float noise now builds its 61x61 grid explicitly.
- A: `layout_preset_grid` (N x M, serpentine port order), custom-rectangles table with overlap check, snap-side-to-whole-cells.
- B/C2: coded 1-bit and N-bit in `/synthesize_ports` (`n_bits`, `feed_keepout`), per-port report, shared state table.
- C3/D: per-port exact cell positions, exact verified STEP, AEDT script generator, per-export README, `/api/export/board`.

## Unreleased - Phase 19 - 3D-printed nylon transmitarray (model-only; Parts A-D)
- New `src/dielectric_cell/`: minimal RCWA for air | z-invariant extruded layers (octagonal pillar / pocket) | air, oblique TE/TM,
  normal-vector field factorisation, multi-region supercell strips. MoM core untouched.
- New `src/transmitarray/`: lens layout, cos^q horn, symmetric dual-beam aperture-phase synthesis, direct per-cell octagon-radius
  solve with the full 2x2 Jones response, Huygens far-field model, one-solid STEP/STL export with re-import verification.
- 28 GHz nylon cell library (`results/nylon_cells_28ghz/`), plate-thickness / through-hole study, 150 mm +-30 deg lens design
  (`results/transmitarray_dual30/`), local-periodicity supercell check, package `hfss_exports/nylon_lens_28ghz_dual30/`.
- Part D-bis: Monte Carlo error bar on the lens prediction from the measured neighbour-coupling error (exact 1x4 / alternating / 2x2 supercells, real 9-cell window
  check, `src/transmitarray/uncertainty.py`, `scripts/research/phase19_*`); pitch-lever estimate; the inter-beam 0.00 dB relabelled as a layout identity.
  STEP coordinates made exact (snap to design values within 1e-7 mm, shortest round-trip reals, explicit 1e-6 mm uncertainty, verified: no merged vertices, valid solid,
  volume change 1.7e-16); the old "no float noise" check had a blind spot below 1e-13 (34,764 inexact literals were in the first STEP).
- Docs: `docs/theory/DIELECTRIC_CELL_RCWA.md`. Everything is labelled model-only (HFSS: `NOT_AVAILABLE_FROM_CURRENT_DATA`).

## Unreleased - Phase 18 - 300 mm four-port boards (minimum-gap study)
- `hfss_exports/four_port_1bit_300mm_gap0/` and `_gap28/`: 4-port 1-bit coded boards with 300.8 mm
  sub-apertures for a minimum-gap/isolation study in HFSS (model-only predictions in each README).
- `scripts/adhoc/coded4.py`: optional `port_side_mm` argument; fine-grid (0.5 x 1 deg) per-port
  metrics and x/y-cut figure (`port_cuts.png`). No solver changes.
- `examples/make_paper_board_26ghz.py`: board-layout figure title uses the actual frequency
  (was hard-coded "26 GHz"; unchanged output at 26 GHz).

## 0.17.0 - Phase 17 - handoff-ready repo (no physics/solver changes)
- Onboarding rewrite: `README.md`, `docs/START_HERE.md` (new), `docs/STATUS.md` (new, single
  source of truth for current validation state), `CLAUDE.md` (rewritten to ≤80 lines),
  `docs/ARCHITECTURE.md` (added the three-layer + research-phase-package module map and the app's
  3-step data-flow diagram), `docs/USER_GUIDE.md` (brought up to date: stepper navigation, port
  defaults, level matching, pointing correction, Coded = model-only), `docs/validation/README.md`
  (new, indexes every validation run with numbers and provenance).
- Docs reorganized: completed phase write-ups and briefs moved to `docs/history/` (git history
  preserved via `git mv`); `docs/briefs/` keeps only active work + `_parked/`; added
  `docs/README.md`, `docs/briefs/README.md`, `docs/references/README.md`.
- `examples/` split: user-facing scripts stay; one-off research drivers moved to
  `scripts/research/`. `results/submission_work/`'s reusable scripts moved to `scripts/adhoc/`;
  scratch logs/probes removed (listed in the Phase 17 report, user-approved).
- `hfss_exports/README.md` added, indexing every package (status: awaiting HFSS / validated);
  the four untracked packages from the Phase 16 live-review round
  (`three_port_25deg`, `four_port_1bit_gap0/28`, `fss_square_loop_10mm`) committed (README/CSV/
  JSON/PNG only — `*.step` files stay gitignored alongside the already-tracked phase12/two_beam
  STEP files).
- App footer now shows its own version/commit (`HMSA Designer · <tag or short sha>`) so a
  screenshot or bug report is traceable; `make app` prints a restart-after-backend-changes hint.

## 0.16.0 - Phase 16 - follow-ups and a live-review round
- Default-on level matching (`hologram.pointing_correction.balance_beam_levels`) for any port or
  design with 2+ beams: iterates internal weights so achieved levels approach the requested ones,
  up to 3 passes, best kept; reported before/after in the beam cards.
- Optional pointing correction (off by default) wired into the single-aperture job for all three
  formulations, with the same before/after report as the multi-port path; later fixed to measure
  on a fine 0.1°-resolution local grid (the shared plotting grid was too coarse to resolve a
  small-aperture squint, which stalled correction for some beams).
- Co-pol side-lobe level (SLL) exclusion-region fix: first a grid-seam bug (a beam at phi=0/180
  read its own far shoulder as its sidelobe), then a direction-space rewrite with a nearest-beam
  (Voronoi) split so one requested beam's own main lobe is never counted as another beam's
  sidelobe.
- Multi-port layout fixes: active-aperture-side ports with the guard gap outside; the 2x2+1
  preset's large port centred correctly on the 2x2 block (was offset by gap/2) and defaulting to
  the four small ports' combined beams; the isolation-gap check tolerant of its own default-value
  rounding and reworded as information, not a fault; feed markers written in the STEP export's own
  frame and checked against each port.
- The per-port/per-beam results table's "Achieved" direction and "Level" now use the same fine
  local measurement as pointing correction and level matching, so the table and the refinement
  summary always agree.
- Restored deleted-but-tracked `hfss_exports/` files; noted in `hfss_exports/phase12_grounded_8mm/
  README.md` that the raw `.s36p` exports must be re-exported from HFSS to re-enable their tests.

## 0.15.0 - Phase 15 - polarization, vortex, multi-feed, multi-port, coded cells
- Root-cause fix for a vortex HFSS found: a monopole feed's radial surface current carries an
  orbital-angular-momentum term in circular polarization that the project's scalar pattern model
  couldn't see. Added a per-beam polarization-compensation term (natural/LHCP/RHCP/HP/VP/slant,
  or a deliberate OAM vortex request) and a genuine VECTOR far-field model via this project's own
  layered Green's function (`src/pattern/vector_aperture_pattern.py`), replacing the scalar model
  for HMSA pattern prediction.
- Multi-feed ("multipole") reference-field construction (superposition/segmentation).
- Multi-port sub-apertures (Mishra et al., IEEE TVT 2026 layout): per-port synthesis, layout
  presets, a first-order edge-distance isolation estimate.
- Coded (1-bit/N-bit) design method: quantizes the phase target to a small number of reflective
  states instead of a continuous value.
- Corresponding Step 2/3 app UI for all of the above, verified live in a real browser.

## 0.14.0 - Phase 14 - surface-wave eigenmode dispersion
- A direct MoM eigenmode search over the periodic array's own Z(kx,ky) was tried and found not to
  resolve a loading-dependent mode at any affordable mesh size (kept as a documented research
  item, `src/dispersion/mom_scaffold.py`).
- Working replacement: `src/dispersion/sheet_model.py` extracts a patch-sheet admittance from the
  already-validated normal-incidence S11 and solves the closed-form transverse-resonance
  condition for beta. Added a "Design method" choice (surface-wave/eigenmode vs. reflection-phase)
  to Step 2 — on the walk-through case the eigenmode method reaches all 5476 cells vs. the
  reflection-phase method's 4006 unreachable.

## 0.13.0 - Phase 13 - production repo + local web app
- Part 0 carry-overs (`docs/history/PHASE13_WEBAPP.md`): Richardson mesh-convergence extrapolation
  at resonance; resonance-region error re-stated for Phase 12B/12C/12D's own meshes; G4/G6 real
  HFSS ingestion with a mesh/thickness/residual error budget; the missing lookup-table-predicted-
  vs-measured per-order-power comparison for Phase 12D's Task B.
- Part 1: git repository initialized; `pyproject.toml` packaging (`pip install -e .[dev,app,cad]`);
  docs moved into `docs/{theory,phases,briefs}/` plus `docs/history/PROJECT_HISTORY.md`; `configs/` for
  config-file presets; ruff lint/format; pytest markers (`slow`/`cad`/`checkpoint`/`hfss_data`) and
  a fast (<5 min) suite; GitHub Actions CI (written, not run); `requirements-lock.txt`; `Makefile`.
- Part 2: structured logging (`LWMOM_LOG_LEVEL`, per-run ID), solver diagnostics capture, a typed
  error-code taxonomy, a fixed validation-label vocabulary, and a run manifest (inputs, resolved
  config, code version, timings, labels, warnings, output hashes) attached to every job.
- Part 3: local FastAPI web app (`make app`) covering unit cell / library / reflectarray / pattern
  prediction / verify / export / HFSS compare / HMSA designer pages; a real HFSS field-calculator
  `.fld` reference-field importer; an HFSS-eigenmode-table dispersion-library import
  (`DISPERSION_FROM_HFSS`) with a documented, tested MoM-eigenmode scaffold
  (`DISPERSION_FROM_MOM`, not implemented) left for later; both hologram formulations from Phase
  8B wired into an HMSA designer; job queue with progress/cancel/config-hash caching.

## 0.0.0 (unreleased, pre-Phase-13) - Phases 1-12D
- **Phases 1-5**: core periodic MoM solver -- geometry/config/units, dyadic layered Green's
  function, rooftop MoM (`Z I = V`), S-parameter extraction at configurable reference planes,
  circular-polarization basis and the unit-cell database schema.
- **Phase 6/6A**: STEP export of representative unit cells; independent full-wave-solver
  validation-export packages (no EM code changed, no automatic HFSS agreement claimed).
- **Phase 7**: production 324-geometry x 201-frequency unit-cell library with a query/phase-
  inversion API (resumable, checkpointed).
- **Phase 8A-8E**: holographic-synthesis feasibility gate; object-wave/multi-beam formulation;
  unit-cell scattering architecture validation; universal unit-cell database architecture; the
  physical control-variable/holographic-scattering bridge (reference-wave -> object-wave ->
  modulation -> control-variable pipeline).
- **Phase 9**: control-variable and local-periodicity validation against real heterogeneous
  supercell solves.
- **Phase 10**: environment-aware neighborhood/surrogate/inverse/beam validation.
- **Phase 11**: physics-based heterogeneous-MoM nonlocal-coupling study (no ML/fitted surrogates).
- **Phase 12A**: optional PEC ground plane (`backing="pec"`), additive to the Green's function and
  config schema.
- **Phase 12B**: mesh alignment / duplicate-geometry audit / grounded period study.
- **Phase 12C**: grounded-supercell local-periodicity/coupling/fine-mesh spot-check; HFSS-ready
  STEP exports.
- **Phase 12D** (+ addendum): design-relevant local-periodicity step-limit re-measurement;
  aperture-level Floquet-order-power test (decided Phase 13's forward-model recommendation); real
  HFSS Touchstone (`.sNp`) ingestion (G1); mesh non-convergence at resonance discovered and
  handed to Phase 13 Part 0 as a carry-over.
