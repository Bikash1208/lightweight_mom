# Phase 21: flexible multi-port layouts and the coded (1-bit/2-bit) method for multi-port

## Why
Today the Ports UI only offers two layouts:
- "2×2 small + 1 large", which is Mishra's 5-port board and is always 5 ports.
- "1×N strip".

The Coded method is also hidden whenever Multi-port is selected (`index.html` ~L788: "Coded isn't wired for multi-port designs in this version").

So the user cannot build in the app the board that was hand-built for HFSS in `scripts/adhoc/square4/` (output in `hfss_exports/four_port_1bit_square_300mm_gap{0,28}/`). That board is:
- 4 ports in a 2×2;
- 300.8 mm sub-apertures (47 × 47 cells of 6.4 mm), gap 0 or 28 mm, one monopole per sub-aperture centre;
- coded 1-bit with a square-only pair (empty cell plus a 5.85 mm square), beta_mode "model";
- 4 RHCP beams at θ = 30°, φ = 45/135/225/315°.

The backend can already do most of this: `custom` layouts exist in `_layout_regions`, and `formulation='coded'` exists in `_run_synthesize_ports_job`. The gaps are in the UI, plus 2-bit support and four correctness fixes.

## Reference result to reproduce (acceptance test)
Run `scripts/adhoc/square4/coded4sq.py` with the `fixcells.py` patch, or build the same request through the API: layout 2×2, side 300.8 mm, gap 0 and 28 mm, coded 1-bit, uploaded pair.

| Input | Value |
|---|---|
| Uploaded state 0 | Γ = 1∠138.14° |
| Uploaded state 1 | Γ = 0.9962∠−42.0° |
| Geometry | empty / rectangle 5.85 × 5.85 mm |
| beta_mode | model |
| Polarisation | RHCP |
| Beams | θ 30°, φ 45/135/225/315° |

Expected results (from `design_report.json` in those folders):
- every port has 47 × 47 = 2209 cells;
- beams at θ = 30.0° with pointing error ≤ 1°;
- directivity 27.1–27.8 dBi;
- LHCP −16 to −20 dB;
- 4416 patches.

Write this as a test (mark it slow if needed), plus a fast test on a small board.

## Work items

### A. Layouts (UI and backend)
1. **New preset `grid`: "N×M grid of equal sub-apertures".**
   - Inputs: rows, columns, side (mm), gap (mm).
   - The layout is centred on the origin, with one feed at each sub-aperture centre.
   - Port order: row-major starting top-right, or document whatever order is chosen. For a 2×2 the order must match the reference: P1 at (+,+), P2 at (−,+), P3 at (−,−), P4 at (+,−). That counter-clockwise order is the one used in the reference board.
   - Add `layout_preset_grid(n_rows, n_cols, side_m, gap_m)` to `src/hologram/ports.py`, with tests.
2. **Expose `custom` in the UI.**
   - Add an editable table of rectangles: x_min, x_max, y_min, y_max and an optional feed x/y in mm, with add/remove rows.
   - Validate that there is no overlap and that every side is greater than 0.
3. **Keep the existing presets unchanged.** The 2×2 + 1 large preset stays at 5 ports by definition and the 1×N strip is unchanged. The help text should say which presets are fixed-count.
4. **Add a "Snap side to whole cells" checkbox (default on).**
   - It rounds the side to the nearest k·period and shows "side 300.0 → 300.8 mm (47 cells)" before the run.
   - Off-grid sides make the combined board resample cells by nearest-neighbour, so on is the right default.

### B. Coded method for multi-port
1. **Remove the UI block that hides Coded in multi-port mode.**
   - Wire 1-bit to the existing backend path (`formulation='coded'`).
   - Offer the same two-state sources as single-port: square-only library pair, square + annular, or uploaded values. The uploaded option needs a geometry for state 1, at least `rectangle L × W`.
   - Add a beta_mode selector (state / paper / model), default **model**.
2. **2-bit (and N-bit) for multi-port.**
   - Add an `n_bits` field to `HmsaSynthesizePortsRequest`.
   - For `n_bits > 1`, call `synthesize_coded_n_bit` per port (same signature pattern as the 1-bit path; it needs `library_id`). Mirror the single-port N-bit validation messages.
   - In the UI, the 1/2/3-bit control becomes active in multi-port mode.
3. **Per-port report.**
   - Report the two-state (or N-state) phase table once.
   - Report per port: beam θ/φ, pointing error, directivity, cross-pol, and patch count.
   - Report the isolation flags as today.

### C. Correctness fixes found while building the reference board
1. **Cell-count floor bug** in `src/fields/reference_field_grid.py::_cell_centers`.
   - It uses `n = int(span // period_m)`. For 300.8 mm / 6.4 mm, floating point gives 46.999… → 46, so a whole-cell side silently loses a row and column.
   - Fix it with a relative tolerance, for example `n = int(np.floor(span / period_m + 1e-6))`.
   - Add a regression test: 300.8 mm at 6.4 mm gives 47 cells, and 300.0 mm at 8.0 mm still gives 37 (the existing docstring example).
   - Check that every existing test still passes. Note in STATUS.md that the earlier `four_port_1bit_300mm_*` boards were affected (46 cells per port were stretched to 47 in the combined board).
2. **Feed keep-out.**
   - Add an option (default on) that forces the cell(s) under each feed to state 0 (empty).
   - If the feed is on a cell corner, that means the 4 touching cells; if it is at a cell centre, the 1 cell under it.
   - Report how many cells were changed per port.
3. **Combined board for gapped layouts.**
   - With a gap that is not a whole number of cells (for example 28 mm at 6.4 mm), `combined_board_lw_grid_m` snaps each port onto one shared lattice, which shifts cells by up to half a period.
   - For export, place each port's cells at their own exact positions: `x = x_min_port + (j + 0.5)·period`.
   - Use the per-port `design.cells` directly. See `scripts/adhoc/square4/build_sq_export.py`, which uses Decimal arithmetic.

### D. Export
1. **Exact STEP.** Build every face from exact vertex coordinates (no transforms).
   - Bodies: Substrate, Ground (`Face.makePlane` at the origin, so there are no 1e-15 plane-origin residues) and Patches.
   - Verify after writing: parse every `CARTESIAN_POINT` and check that it lies on the 0.001 mm grid and that none print like `…9999…` or `…0000x`. Fail the export otherwise.
   - Reuse or port the logic in `scripts/adhoc/square4/build_sq_export.py`.
2. **AEDT IronPython script.** Port `scripts/adhoc/square4/aedt_template.py` into `src/export/` as a generator.
   - It builds substrate, ground with coax holes, the patches united into "Patches", PEC monopoles, lumped ports, the air box with a radiation boundary, the setup, the sweep and a far-field sphere.
   - Design variables: mono_h, r_pin, r_out, pads.
   - Add an Export-page button: "HFSS script (.py)".
3. **README per export.** Include the port/feed table, the beam table, the cell-pair phases, and the "excite one port / all ports (Edit Sources)" instructions.

## Constraints
- Don't push or add remotes. Commit per part, and don't move tags.
- Don't rename `src/` packages. Never commit PDFs, .s36p or .s4p files.
- Keep the single-port behaviour and the existing presets byte-identical, except for the C1 fix. Note any changed numbers in the report.
- Run the full test suite before each commit.

## Report back
- What changed per part, and the test counts.
- The acceptance-test numbers versus the expected table above.
- Screenshots of the new Ports UI: grid preset, custom table, Coded with 1-bit/2-bit in multi-port.
