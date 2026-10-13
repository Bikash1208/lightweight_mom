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
| `fss_square_loop_10mm/` | Square-loop FSS unit cell, periodic MoM + finite-11×11-panel plane-S (`scripts/adhoc/fss_loop.py`, `scripts/adhoc/fss/plane_s.py`) | **Self-consistency reference, not HFSS-validated** — no HFSS data exists for this geometry (`docs/validation/README.md` item 5). No `*.step` (CSV/PNG only). |
| `nylon_lens_28ghz_dual30/` | 28 GHz all-dielectric (nylon) horn-fed transmitarray lens, ±30° dual beam, 150 mm octagon, 25 mm thick (`examples/export_nylon_lens_28ghz.py`) | **Awaiting HFSS** (model-only; its README carries the Monte Carlo error bar on the prediction, the pitch estimate and the exact-coordinate record). `lens.step`/`lens.stl` gitignored, regenerable. |
| `nylon_cells_28ghz/` | Six optional single nylon cells spanning the phase range with RCWA reference CSV (`examples/export_nylon_cell_checkset_28ghz.py`) | **Awaiting HFSS** (optional spot check). |

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

