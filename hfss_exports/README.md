# HFSS export packages

Each subfolder is a self-contained HFSS import package (STEP + layout CSV/JSON/DXF + an
`HFSS_README.md` with Floquet-port/PEC/thickness setup notes) for one design, plus whatever
real HFSS results exist for it. Raw HFSS Touchstone exports (`*.s36p`, `*.s4p`) are gitignored —
`docs/DEVELOPING.md` has the full account of what's excluded and why.

| Folder | What | Status |
|---|---|---|
| `phase12_grounded_8mm/` | Single-cell (G1–G6) and supercell (S1–S3) grounded 8 mm unit cells | **Validated**: G1 off-resonance (mean phase diff -0.29°); G4 at-resonance (13.51° diff, decomposed into mesh/thickness/residual — `docs/validation/README.md`). G2/G3/G5 and all supercells await an HFSS run. |
| `two_beam_20deg_phi0_180/` | Two-beam (20°, phi=0/180°) single-aperture board, 150 mm and 300 mm | Awaiting HFSS. |
| `three_port_25deg/` | 1×3 port board, 25° RHCP beams per port (`scripts/adhoc/three_port.py`) | **Awaiting HFSS** (`docs/STATUS.md` "Open next steps"). `*.step` gitignored; everything else committed. |
| `four_port_1bit_gap0/`, `four_port_1bit_gap28/` | 4-port 1-bit coded boards, two guard-gap variants (`scripts/adhoc/coded4.py`) | **Awaiting HFSS** (`docs/STATUS.md` "Open next steps"). `*.step` gitignored; everything else committed. |
| `four_port_1bit_square_gap0/`, `four_port_1bit_square_gap28/` | 4-port 1-bit coded boards, plain square cells only (state 0 empty, state 1 a 5.85 mm square), 140.8 mm sub-apertures (22×22 cells), beams (30°, 45/135/225/315°) RHCP, gap 0 vs 28 mm | **Awaiting HFSS**. Hand-built by `scripts/adhoc/square4/`; the app reproduces them exactly (identical patch sets, 960 patches). Regenerate: Multi-port page, grid 2×2, cell side 140.8 mm, gap 0/28 mm, coded 1-bit, uploaded square pair, `beta_mode=model`, then Export board (same request as `_grid_request(140.8, gap)` in `tests/test_app_phase21_multiport.py`). `*.step` gitignored; everything else committed. |
| `four_port_1bit_square_300mm_gap0/`, `four_port_1bit_square_300mm_gap28/` | Same square-cell 4-port board with 300.8 mm (47×47-cell) sub-apertures | **Awaiting HFSS**. Hand-built by `scripts/adhoc/square4/`; the app reproduces them exactly (4416 patches; `test_acceptance_four_port_square_300mm_board`). Regenerate as above with cell side 300.8 mm. `*.step` gitignored (15.7 MB); everything else committed. |
| `four_port_hex150_gap0/`, `four_port_hex150_opt_single/`, `four_port_hex150_opt_all/`, `four_port_hex150_opt_single_robust/`, `four_port_hex150_opt_all_robust/`, `four_port_square243_gap0/` | Phase 23 comparison boards, 10.54 GHz 4-port 1-bit coded: four 150 mm-side hexagons (gap 0), the same with cell states optimised for single-port excitation / for all ports together (port phases in its README), and the equal-area 243.2 mm square (`scripts/research/phase23_compare.py`; predictions in `results/phase23/comparison_v2_held_out.json`; `*_robust`: worst case over 10.33/10.54/10.75 GHz and, all-ports, port-phase errors) | **Awaiting HFSS** (model-only). `*.step` gitignored; everything else committed. |
| `four_port_1bit_300mm_gap0/`, `four_port_1bit_300mm_gap28/` | Minimum-gap study: same 4-port 1-bit design with 300.8 mm (10.6 λ0) sub-apertures, beams (30°, 45/135/225/315°) RHCP, gap 0 vs 28 mm (`scripts/adhoc/coded4.py ... 300.8`) | **Awaiting HFSS**. Model: 27.6 dBi per port, pointing ≤ 0.5°, co-pol SLL −11 to −12 dB along x/y. `*.step` gitignored; everything else committed. |
| `cell_pair_locality_check/` | Phase 25: neighbour-coupling check of the 1-bit cell pairs — 1×2, 2×2 and 5×5/4×4 square \| empty supercells for the current 6.4 mm / 5.85 mm pair and two candidates, plus MoM reference CSVs (`scripts/research/phase25_locality_package.py`, `phase25_locality_blocks.py`) | **Awaiting HFSS** (model-only; decides whether the Phase 21–23 board numbers can be trusted). The small `*.step` are committed with `git add -f`. |
| `fss_square_loop_10mm/` | Square-loop FSS unit cell, periodic MoM + finite-11×11-panel plane-S (`scripts/adhoc/fss_loop.py`, `scripts/adhoc/fss/plane_s.py`) | **Self-consistency reference, not HFSS-validated** — no HFSS data exists for this geometry (`docs/validation/README.md` item 5). No `*.step` (CSV/PNG only). |
| `nylon_lens_28ghz_dual30/` | 28 GHz all-dielectric (nylon) horn-fed transmitarray lens, ±30° dual beam, 150 mm octagon, 25 mm thick (`examples/export_nylon_lens_28ghz.py`) | **Awaiting HFSS** (model-only; its README carries the Monte Carlo error bar on the prediction, the pitch estimate and the exact-coordinate record). `lens.step`/`lens.stl` gitignored, regenerable. |
| `nylon_cells_28ghz/` | Six optional single nylon cells spanning the phase range with RCWA reference CSV (`examples/export_nylon_cell_checkset_28ghz.py`) | **Awaiting HFSS** (optional spot check). |
| `hex150_tensor_cp_single/` | Phase 26: Single 150 mm hexagon, central monopole feed, rotating-slot tensor hologram for pure RHCP beam at (30°, 225°), 10.54 GHz on RT/duroid 5880 1.575 mm (`scripts/research/phase26_fourhex_design.py`, `src/export/tensor_cp_export.py`) | **HFSS Verified** (measured in `hfss_results/hex150_tensor_cp_single/`; `docs/validation/README.md` section 7). `board.step` gitignored; layout files committed. |
| `hex150_tensor_cp_single_v2_pad/` | Phase 26: Single 150 mm hexagon with V2 graded launch pad (cross-slotted cells) and calibrated launcher pin (`scripts/research/phase26_launcher_package.py`) | **Model-only, V2 run on hold** (`docs/STATUS.md:47`). `board.step` gitignored; layout committed. |
| `hex4_tensor_cp_static/` | Phase 26: 4-hexagon tensor-CP Board A (static beams), 652.8 × 537.6 mm board, shared 6.4 mm lattice with 5684 slotted-circle cells (`scripts/research/phase26_fourhex_design.py`, `src/export/tensor_cp_export.py`) | **Awaiting HFSS** (Board A run parked per `docs/STATUS.md:79`). `board.step` gitignored; layout committed. |
| `hex4_tensor_cp_combined/` | Phase 26: 4-hexagon tensor-CP Board B (combined beam, four hexagons aimed at (30°, 225°), DRIVE_PORT=5). | **Branch finding:** Export files (`README.md`, `hfss_build_board.py`, `layout.csv`) exist in commit `764224a` on branch `phase26-leaky-tensor-cp` only, and are not yet merged into `phase20-ris`. |
| `launcher_calibration/` | Phase 26: Bare grounded RT/duroid 5880 substrate board (no hologram), 6.0 mm coaxial probe pin, matching `four_port_hex150_gap0` outline (`scripts/research/phase26_launcher_package.py`) | **Calibration Package**; measures launcher efficiency and reflection in isolation without hologram leakage. |
| `phase20_r1_loaded_cell/` | Phase 20 R1 item 4: One-cell Floquet package for loaded, oblique MoM verification (15 mm cell, lumped RLC loads, theta 0° and 30°) (`scripts/research/phase20_r1_hfss_package.py`) | **Model demonstration cell**, not a final design (`PREDICTION_NOT_HFSS_VALIDATED`). `cell.step` gitignored. |
| `phase20_r2_diode_cell/` | Phase 20 R2: One-cell Floquet package for PIN diode model verification (SMP1345-079LF ON/OFF states, theta 0° and 30°) (`scripts/research/phase20_r2_hfss_package.py`) | **Model demonstration cell** (`PREDICTION_NOT_HFSS_VALIDATED`, diode values `PROVISIONAL_DATASHEET`). `cell.step` gitignored. |
| `phase20_r3_cp_cell/` | Phase 20 R3 item 5: Candidate 2-bit Pancharatnam-Berry octagon unit cell (CP-2) with MACOM MADP-000907-14020 diodes, 4 discrete states, bare slab reference (`scripts/research/phase20_r3_cp_package.py`) | **HFSS Verified** (user run 2026-10-23 in `hfss_results/phase20_r3_cp_cell/`; confirmed 10.5 GHz operating point). `cell.step` gitignored. |
| `slotted_cell_eigenmode/` | Phase 26 B2.5: Slotted-circle unit cell eigenmode checks on RT/duroid 5880 (`scripts/research/phase26_eigenmode_package.py`) | **Optional package, NOT run**; internal periodic MoM solver relied upon for near-resonant Bloch dispersion. |

## What a new package commits (growth policy, Phase 24)

| Committed | Not committed (gitignored; regenerate) |
|---|---|
| `README.md`, `design_report.json`, PNGs, the small `layout.csv` / `feeds.csv`, the AEDT script | `*.step`, `*.stl`, and `layout.dxf` / `layout.json` (1-4 MB each on a 300 mm board) |

- `hfss_exports/**/*.step`, `**/layout.dxf` and `**/layout.json` are in `.gitignore`. Files **already tracked**
  stay tracked (`phase12_grounded_8mm/`, `two_beam_20deg_phi0_180/`, the Phase 18-23 packages whose layouts
  were committed); gitignore never untracks. A layout under 1 MB may still be committed on purpose with
  `git add -f`.
- Every package README (and its row in the index above) names the command that regenerates the ignored files:
  the app's Export page on the matching `/synthesize_ports` request, or the generating script.
- Nothing tracked was removed; no history was rewritten.
