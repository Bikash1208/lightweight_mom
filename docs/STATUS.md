# Status — single source of truth for *now*

**Last updated:** 2026-10-15, by Claude (Phase 26 on branch `phase26-leaky-tensor-cp`, from `phase25-butler-sector`: single-hexagon tensor-CP board checked against HFSS; 4-hex board scoped; Phase 24 curation merged into `main`, tagged `v1.1-handoff`; the Phase 25 locality caveat below stands). **Rule: every phase updates this
file** (last release tag/commit, capability matrix, open next steps) as part of its own final commit.

**Last milestone tags:** `v1.1-handoff` (this state), `phase23-complete` (hexagon ports + sidelobe optimiser), `phase22-complete`, `phase21-complete`,
`phase19-partD`, `phase16-complete`; `pre-cleanup-2026-10-07` / `pre-cleanup-2026-10-13` mark the two cleanups and
`submission-2026-10-06` is the frozen submission snapshot. Tags are never moved or deleted.

> **Caveat for every Phase 21–23 coded board (added in Phase 25, 2026-10-13): local periodicity is poorly supported
> for the 0.55 mm-gap 1-bit pair.** The boards use empty | 5.85 mm square on a 6.4 mm cell (gap 0.55 mm). In the repo's
> only full-array MoM (periodic supercell, normal-incidence plane wave) the isolated-cell model of that pair is off by
> 42–69 percentage points in the worst Floquet-order power on 5×5 blocks, with the error growing as the gap shrinks
> (3–10 pp at 1.07 mm, 12–30 pp at 0.80 mm, 33–57 pp at 0.64 mm); small 1×2 / 2×2 blocks at 2× finer mesh give the
> same picture, so it is not a mesh artefact. It is not an HFSS result and not a surface-wave-fed board, so the size of
> the effect on the real boards is **unknown**, but every "model-only" pattern number for these boards (directivity,
> pointing, SLL, level matching) rests on an assumption that check does not support. Details, tables and the named
> state `LOCAL_PERIODICITY_POOR_AT_DESIGN_GAP_PLANE_WAVE_CHECK`: [`results/phase25/README.md`](../results/phase25/README.md).
> The error tracks the empty-vs-square phase separation, not the gap: at 6.4 mm no square with a ≥ 1.0 mm gap reaches
> 90° separation (73–90°; its locality error is 4–12 pp), and an 8.0 mm cell with a 6.85 mm square reaches 178.7° but
> is still 11–55 pp off. No pair with ≥ 150° separation and small locality error was found.

## Capability matrix

| Capability | Validation level | Numbers / where |
|---|---|---|
| Dielectric (nylon) unit-cell RCWA, oblique TE/TM (Phase 19 Part A, `src/dielectric_cell/`) | Verified against an independent analytic slab formula (agreement ~1e-14) and internal checks; **not HFSS-validated** | Uniform slab vs Fabry–Pérot/Macleod 1e-14; power balance 1e-14 (φ=0) / <1e-8; Fourier order M=8 converged to <0.2° / 0.01 dB (`docs/theory/DIELECTRIC_CELL_RCWA.md`) |
| 28 GHz nylon cell library + dual-beam transmitarray lens (Phase 19 Parts B–D and D-bis, `src/transmitarray/`) | **Model-only, awaiting the user's HFSS run; local-periodicity error propagated** | **Isolated-cell prediction** (deterministic): 25.6 dBi/beam, continuous-peak beam angle 29.87° (nominal ±30°), SLL −9.7 dB, 0 of 756 cells >1° phase residual. **With the measured neighbour-coupling error** (Monte Carlo, 2000 draws, median (p10…p90); percentiles stable to ≤0.04 dB gain / 0.14 dB SLL): primary phase-constrained model 21.5 (21.1…22.0) dBi, beam 30.26 (30.00…30.46)°, SLL -8.2 (-9.2…-7.1) dB, |inter-beam level difference| 0.4 (0.1…0.9) dB; random-phase model 25.6 (25.3…26.0) dBi, SLL -8.8 (-9.8…-7.8) dB. The primary model is only partly validated on four real 9-cell windows solved exactly (isolated-cell order-spectrum error 0.9–1.4 of the signal; the model halves it in places and under-predicts power by ~0.4–1.4 dB), so the true gain change lies between ≈0 and the primary figure. **The nominal 0.00 dB inter-beam level difference is an identity of the mirror-symmetric layout, not a result.** Pitch λ/3 (150/42 mm, non-decimal) estimate: 1.9× cells, pockets 14:1, median neighbour step 43°→27°, coupling-error rms ×0.57. STEP coordinates are exact (snapped, shortest round-trip; `hfss_exports/nylon_lens_28ghz_dual30/README.md`). Details: `results/transmitarray_dual30/`. |
| Periodic unit-cell solver (ungrounded, Phases 1–9) | HFSS-validated, off-resonance | Case G1: mean phase diff **-0.29°**, max \|S\| diff **~1.6e-5** (`docs/history/PHASE13_WEBAPP.md`) |
| Periodic unit-cell solver, at/near resonance | Mesh-convergence-characterized, residual unexplained | G4 (12.5 GHz): HFSS-vs-0.2mm-mesh diff 13.51° = mesh 5.93° + thickness 3.66° + **residual 3.92° unexplained** (`docs/history/PHASE13_WEBAPP.md`) |
| Grounded supercell / local periodicity (Phase 12C/12D) | Cross-checked against real full-array MoM | Lookup-table-vs-measured per-order power differs by 0.17–1.8 percentage points (`docs/history/PHASE12D_LOCAL_PERIODICITY.md`) |
| Vortex↔pencil physics (Phase 15 Part A/B) | Independently cross-checked (this project's own vector far-field model, not HFSS) | RHCP/LHCP pencil XPD **21.4 dB**; HP/VP pencil XPD 62–341 dB (scalar-cell ceiling); deliberate OAM request: winding **−1 matches charge −1**, null **−31 dB** (`docs/history/PHASE15_POLARIZATION_VORTEX.md`) |
| Hexagonal (polygon) multi-port sub-apertures, `hex4` preset (Phase 23 A, `src/hologram/ports.py`) | Model-only; geometry checked against analytic area/membership, rectangle layouts reproduce **bit-identically** (`results/phase23/rect_snapshot_*.json`); hexagon STEP exact (verified); the trimmed-board AEDT script is generated and unit-tested but not run in AEDT | 4 × 150 mm-side hexagons, gap 0: 58,457 mm² each, 1414–1435 cells, directivity 25.8–26.0 dBi (−1.8 dB vs the 300.8 mm square), co-pol SLL −14.1…−13.7 dB vs −11.5…−11.1 dB for the square and −11.2…−10.9 dB for the equal-area 243.2 mm square (held-out 0.5°×1° grid, pointing-corrected starts) — `results/phase23/comparison_v2_held_out.json` |
| Sidelobe-cancelling 1-bit optimiser: single-port, all-ports (cells + port phases), robust (Phase 23 B/C, `src/pattern/cellwise_pattern.py`, `src/hologram/sll_optimizer.py`, `sll_board.py`) | Model-only; every number on a held-out fine grid the optimiser never saw, cells re-solved (Phase 22) at 10.33/10.54/10.75 GHz; **not HFSS-validated**. The all-ports figure needs ONE phase-locked source split four ways (independent per-port signals: use the single-port metric); the beam squints with frequency (0.5–1.5° over ±2 %), the 0.3° absolute pointing limit holds at the design frequency in single-port mode only (all-ports starts at 0.4–0.7° and is not made worse) | hex4, single-port: SLL −17.2…−16.8 dB (from −14.1…−13.7), −0.1 dB directivity, 47–56 cells/port; at 10.33/10.75 GHz −13.4…−12.5 / −16.6…−16.1 (robust version: −14.2…−13.7 / −15.8…−15.6 at ±2 %, −16.2…−15.6 at 10.54). all-ports: −23.3…−21.7 dB at 10.54, −18.3…−15.8 at 10.33, −18.9…−16.3 at 10.75 (robust −22.5…−22.0 / −20.7…−19.1 / −19.3…−18.0); 20 random ±10° port-phase errors: worst-beam median −21.7 dB (−17.9 robust at 10.75). Chosen port phases 0°/0°/0°/0° |
| Four-port 1-bit boards, 300.8 mm sub-apertures, gap 0 vs 28 mm: minimum-gap (isolation) study (Phase 18, `scripts/adhoc/coded4.py`, `scripts/adhoc/square4/`) | Model-only. The Phase 18 "x/y-cut co-pol SLL" metric (`sll_cut_along_x/y_db`) is the one Phase 23 reports | 27.6 dBi per port, pointing ≤ 0.5°, SLL −11…−12 dB; `hfss_exports/four_port_1bit_{300mm,square_300mm}_gap*` awaiting HFSS |
| Automatic coded-cell states (Phase 22, `src/hologram/auto_states.py`) | Model-only: mesh-corrected MoM (Richardson) interpolated between mesh-aligned solves; the 1-bit reference pair is empty +138.1° + square 5.85 mm −41.4° (separation 179.6°), the 4-port 300 mm run with it gives 27.7–27.8 dBi, beams within 0.3° | `tests/test_hologram_auto_states.py`, `tests/test_app_phase22_*.py` |
| Multi-port layouts: N×M grid, custom rectangles, coded 1-/N-bit, exact STEP + AEDT script (Phase 21) | Model-only; the 4-port 300.8 mm square-cell board is reproduced (47×47 cells/port, beams θ=30° within 0.3°, 27.7–27.8 dBi, 4416 patches, layout identical to the hand-built reference) | `tests/test_app_phase21_multiport.py`, `tests/test_export_exact_board.py` |
| Multi-port sub-apertures (Phase 15 Part C gap-fill) | Model-only; reproduces Mishra et al.'s *target directions*, not yet checked against that paper's measured isolation/gain | Mishra 26 GHz 5-port reproduction: pointing error 0–11°, port-5 beam levels matched to ≤0.1 dB after equalisation (`docs/submission/paper26/comparison.md`); **HFSS packages awaiting simulation**, see below |
| 1-bit / N-bit coded (Phase 15 Part E) | Model-only; cell-pair phase separation is mesh-sensitive | Paper-mode pair at 26 GHz: 173–179° separation across 25–27 GHz; **HFSS check of this exact cell pair not yet done** |
| Level matching / pointing correction (Phase 16 follow-ups) | Self-consistent (one measurement method used everywhere, Phase 16 live-review round) | 4-beam single port: level spread 6.4 dB → 0.5 dB in 2 passes; pointing miss converges to <0.3° on a fine 0.1° local grid |
| Polarization efficiency (scalar, isotropic unit cells) | Documented physical ceiling, not a bug | ≈ 50% co-pol power fraction for a linearly-polarized compensated beam (`docs/history/PHASE15_POLARIZATION_VORTEX.md` "Scalar, isotropic unit cells cap achievable XPD") |
| Surface-wave leakage + launch model for a monopole-fed hologram (Phase 26 A, `src/dispersion/leakage_model.py`, `pin_launch.py`, `src/hologram/leaky_forward.py`) | Model; closed-form first-order α validated against a truncated TM+TE Floquet balance (≤ 0.5 %); the pin launch split against free-space integrals (0.01–0.1 %). Against HFSS: brackets the 1-bit board (16.9–19.9 vs 19.85 dBi, `LEAKAGE_MODEL_BRACKETS_HFSS_NOT_DISCRIMINATING`), and with the tensor-CP board the gain bookkeeping launch 0.59–0.62 × conversion 0.54 agrees to 0.5 dB | `results/phase26/README.md` A–A4, H |
| Tensor (anisotropic) cells: rotating-slot circles, Bloch (kt > k0) sheet tensor, supercell β/α table (Phase 26 B/B2, `src/dispersion/tensor_sheet.py`, `slotted_cells.py`, `tensor_floquet.py`) | Solver-internal (`SOLVER_INTERNAL` / `SOLVER_SUPERCELL`), no solver change; the normal-incidence tensor is off by up to 15 % of β across a slot; a Floquet balance of the cells' local tensors does NOT reproduce the supercell for the strongly anisotropic primary cell (0.013 k0 > 0.01 k0 gate), so the design uses the supercell table (nearest grid, collinear modulation only) | `results/phase26/README.md` B–B2.4f; mesh-fragile cells (bridge < 0.35 mm) excluded |
| **Single-hexagon pure-RHCP tensor-CP board, `hex150_tensor_cp_single` (Phase 26 C/D)** | **HFSS-validated prediction (one run, 10.5 GHz, one board)**: 23.03 dBi RHCP at (30°, 225°), AR 1.2–1.3 dB, XPD 23 dB, SLL −18.2 dB (user figure), S11 −15.0 dB; model predicted 23.5 dBi / AR 2.5 / XPD 16.9 / SLL −22.8. Unconverted surface wave + pin direct radiation = 0.70 of the power (grazing 0.26, back 0.27, LHCP 0.31) | `docs/validation/README.md` §7, `hfss_results/hex150_tensor_cp_single/`, `hfss_exports/hex150_tensor_cp_single/` |
| Launcher on the loaded surface, graded pad V2 (Phase 26 item 5), aperture-size study | Model-only, **V2 not run (on hold)**; feed calibrated to the HFSS S11 sweep (R × 1.30, effective pin radius 1.5 mm, `HFSS_CALIBRATED`, one board); size study HFSS-anchored: 300 mm side ≈ 30.2 dBi (model) | `results/phase26/README.md` I, I2, J, K; `hfss_exports/hex150_tensor_cp_single_v2_pad/` |

## Known limitations

- **Cell-count floor bug fixed (Phase 21 C1).** `fields.reference_field_grid._cell_centers` floored
  300.8 mm / 6.4 mm to 46 cells (46.999… in floating point). The earlier `four_port_1bit_300mm_*`
  boards were affected: 46 cells per port were stretched onto 47 in the combined board. A
  whole-cell side now keeps all its cells; `test_pattern_beam_polarization_metrics`'s two-beam
  fixture had relied on the same float noise and now builds its 61x61 grid explicitly.
- **Phase 23 limits**: the hexagon trim-board AEDT script is generated and unit-tested but has never been run in AEDT; an
  all-ports sidelobe figure holds only for ONE phase-locked source split across the ports (independent per-port signals: use the
  single-port figure); every optimiser number is a model, held-out grid, not HFSS.
- **RIS (Phase 20)** is paused, not started: nothing in `src/` implements PIN-diode loads or oblique incidence in the MoM core.
- Normal incidence only; single substrate layer; PEC zero-thickness metal model (a real HFSS
  patch has finite thickness — see the G4/G6 residual above).
- Transmitarray support is **model-only** (Phase 19, dielectric horn-fed lens; no HFSS comparison yet); the MoM solver itself remains reflectarray/metasurface only. No finite-array (N×N, N not a power of a small factor) full-wave MoM
  beyond what `Verify` already runs for a cost-bounded N.
- Multi-port coded designs (Phase 21): 1-bit and N-bit are wired in the Ports UI and `/synthesize_ports`
  (`n_bits`, shared cell pair, feed keep-out for 1-bit). Model-only; the 4-port square-cell reference
  board's predicted numbers are not yet compared with an HFSS solve.
- Geometry inversion and full-array STEP export remain `GEOMETRY_INVERSION_LOCKED` /
  `STEP_EXPORT_LOCKED` scientific-finding records (Phase 12D) even though the app produces designs
  and exports anyway, per the project's Phase 13 decision — every such output still carries these
  labels.
- No run manifest / diagnostic-bundle download wired into the app yet (`src/manifest.py` exists and
  is tested, nothing attaches it to a job — `docs/ARCHITECTURE.md`).

## Open next steps

Order of work, with what each needs from the user (HFSS runs are the user's):

0. **Phase 26 (branch `phase26-leaky-tensor-cp`, brief `docs/briefs/PHASE26_LEAKY_TENSOR_CP_HEXAGON.md`)**: board A (4-hex static beams, shared 6.4 mm lattice, 5685 cells) is exported to `hfss_exports/hex4_tensor_cp_static/` with predictions and the neighbour-LHCP-lobe check item — run it in HFSS (one port driven, others 50 ohm) and send S-matrix, per-port RHCP/LHCP gain, AR; board B (Butler-ready, all four at (30, 225), global phase reference) is scoped only (`results/phase26/README.md` N; the +/-135 deg Butler inputs do not form a single beam in the model), export waits for the user. The single-hex E-line (|E| only) was non-discriminating; V2 pad stays on hold. Merge/tag of the branch and the move of the brief to `docs/history/briefs/` follow the board A/B HFSS results.
1. **Decision, Phase 19**: is ±30.00° a hard specification for the nylon lens? The error bar puts the beam at ≈30.3° median in the
   primary model (0.4° typical pointing error). Parts E/F of `docs/briefs/PHASE19_DIELECTRIC_TRANSMITARRAY_28GHZ.md` wait on this.
2. **HFSS run, one unit-cell check first**: the 10.54 GHz 1-bit pair (empty +138.1° and the 5.85 mm square −41.4°) as a single cell
   (`scripts/adhoc/one_bit_pair_cells.py` writes the check package; not yet run) — every coded board below depends on this pair.
3. **HFSS runs of the 4-port 1-bit boards**, in this order (each README lists model numbers to compare against):
   1. `hfss_exports/four_port_1bit_square_gap0/` and `_gap28/` (140.8 mm sub-apertures; cheapest);
   2. `hfss_exports/four_port_1bit_square_300mm_gap0/` and `_gap28/` — compare S21/S41 between the two (minimum-gap study);
   3. Phase 23: `four_port_hex150_gap0/` (baseline), then `four_port_hex150_opt_single/`, `_opt_all/` (with the port phases of its README
      from ONE phase-locked source) and the `_robust` pair; `four_port_square243_gap0/` is the equal-area square control. Open follow-up:
      all-ports absolute pointing is 0.2–0.7° (not ≤ 0.3°) because the ports' sidelobes interfere in the beams; the optimiser moves cells
      only (a pointing penalty term or a Phase 16 pass on the all-ports pattern would fix it);
   4. earlier boards: `four_port_1bit_gap0/`, `_gap28/`, `four_port_1bit_300mm_gap0/`, `_gap28/`, `three_port_25deg/`.
4. **HFSS run, Phase 19**: `hfss_exports/nylon_lens_28ghz_dual30/lens.step` (beam angle/gain/level), optionally the six cells in
   `hfss_exports/nylon_cells_28ghz/`; then decide Part E (transmitarray in the app).
5. **HFSS neighbour-coupling check (Phase 25)**: `hfss_exports/cell_pair_locality_check/` — 1×2, 2×2 and 5×5/4×4
   square | empty supercells for the current 6.4 mm pair and two candidates (6.4 mm / 5.33 mm square, locality-best;
   8.0 mm / 6.85 mm square, separation-best), with the MoM reference CSVs. Run this before trusting any Phase 21–23
   board number or choosing a new pair.
6. **Housekeeping**: the Phase 24 branch `phase24-cleanup` can be deleted once the user is happy with `v1.1-handoff`
   (`git branch -d phase24-cleanup`; the tag preserves it). The 2550-test count assumes the gitignored data is present: a fresh
   clone runs 2441 passed / 129 skipped (the Phase 7 library checkpoint, real HFSS exports and playwright are absent) — see
   `docs/DEVELOPING.md`.

## Hygiene (Phase 24 §7 findings — report only, not bulk-fixed)

Measured 2026-10-13 on `phase24-cleanup` (ruff 0.16.9, Python 3.14, macOS):

- **Fast suite**: 2550 passed, 20 skipped, 117 deselected in **6 min 18 s** (a fresh clone with the gitignored data absent: 2441 passed, 129 skipped in 6 min 42 s; the earlier < 5 min target is not met). The time is spread over
  many tests rather than a few slow ones; the slowest, in order: `test_phase8c_data_source::test_load_full_database_has_324_cells...`
  (20 s setup), `test_app_phase23_hex4::test_sidelobe_endpoint_end_to_end` (15 s), `test_pattern_beam_polarization_metrics`
  two-beam fixture (8 s), `test_hologram_multi_port_beam_pointing::test_1x3_rhcp_board...` (6.6 s), then ~15 tests at 4-6 s
  (the 300 mm board acceptance tests, exact-board STEP fixtures, Phase 8E/9 run tests). Left unchanged on purpose.
- **Lint debt**: `ruff check src tests examples scripts` reports **207** findings (the Phase 17 tag, same ruff version, reports 122; the
  "~95" earlier figure came from an older ruff). Growth is in the Phase 18-23 code: `scripts/adhoc/square4/*` (~58, incl. 8
  `invalid-syntax` in `aedt_template.py`, which is an IronPython template, not importable Python), `src/export/board_readme.py` and
  `aedt_script.py` (24, mostly `UP031` %-formatting inside generated text), and tests of Phases 21-22 (~25). Top rules: `B905` bare
  `zip` 41, `UP031` 28, `F841` 22, `F811` 19, `I001` 16. Each needs per-site judgement; no bulk `--fix`.
  `ruff format --check` would reformat 24 files (mostly `scripts/adhoc/`); not applied.
- **Dead code**: the import-graph scan (`ast`, relative imports resolved, over `src/`, `tests/`, `scripts/`, `examples/`) finds
  284 modules and none that is never imported (`src.pattern` is a package only imported through its submodules). The Phase 19
  packages `src/dielectric_cell` and `src/transmitarray` are imported only by `examples/`, `scripts/` and tests, not by the app (by design).
- **Debug prints / TODOs**: none in `src/` outside the CLI report printer `src/phase8c/run.py`; no unowned `TODO`/`FIXME`/`XXX`.
- **Platform-sensitive test fixed**: `tests/test_phase.py::test_branch_cut_negative_real_axis_is_plus_180_for_both_signed_zeros`
  built the negative zero with the literal `-1 - 0j`, whose imaginary sign depends on the Python version (it failed on Linux
  Python 3.12). It now uses `complex(-1.0, -0.0)`; the asserted behaviour of `principal_phase_deg` is unchanged.

## Parked ideas (`docs/briefs/_parked/`)

- `PHASE16_UNIT_CELL_SPARAM_TOOL.md` — a unit-cell/FSS tool with plane-defined S-parameters (the
  measurement-plane technique `scripts/adhoc/fss/plane_s.py` prototyped for the
  square-loop FSS reference, `hfss_exports/fss_square_loop_10mm/`). Not implemented.
- `PHASE20_RIS_ROADMAP.md` and `PHASE20_R1_SOLVER_LOADS_OBLIQUE.md` — the 2-bit PIN-diode reconfigurable-metasurface (RIS)
  roadmap and its first item (solver loads + oblique incidence). **Paused, not started**; nothing in `src/` implements it.
- Also parked, no brief written yet: transmitarray support in the app (Phase 19 Part E), oblique incidence in the MoM core, tensor
  (anisotropic) unit cells (would relax the scalar-cell XPD ceiling above), 2-bit multi-port.
