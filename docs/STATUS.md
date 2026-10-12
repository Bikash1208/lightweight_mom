# Status — single source of truth for *now*

**Last updated:** 2026-10-12, by Claude (Phase 23 A-E, branch `phase23-hex-sll` from `phase22-auto-cell-states`, not yet merged/tagged; Phase 19 Parts A-D on `phase19-transmitarray` as before). **Rule: every phase updates this
file** (last release tag/commit, capability matrix, open next steps) as part of its own final commit.

**Last release tag:** `phase16-complete` (commit `1c5ebfb`). The submission snapshot tag
`submission-2026-10-06` is a separate, frozen checkpoint — never move it.

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
| Automatic coded-cell states (Phase 22, `src/hologram/auto_states.py`) | Model-only: mesh-corrected MoM (Richardson) interpolated between mesh-aligned solves; the 1-bit reference pair is empty +138.1° + square 5.85 mm −41.4° (separation 179.6°), the 4-port 300 mm run with it gives 27.7–27.8 dBi, beams within 0.3° | `tests/test_hologram_auto_states.py`, `tests/test_app_phase22_*.py` |
| Multi-port layouts: N×M grid, custom rectangles, coded 1-/N-bit, exact STEP + AEDT script (Phase 21) | Model-only; the 4-port 300.8 mm square-cell board is reproduced (47×47 cells/port, beams θ=30° within 0.3°, 27.7–27.8 dBi, 4416 patches, layout identical to the hand-built reference) | `tests/test_app_phase21_multiport.py`, `tests/test_export_exact_board.py` |
| Multi-port sub-apertures (Phase 15 Part C gap-fill) | Model-only; reproduces Mishra et al.'s *target directions*, not yet checked against that paper's measured isolation/gain | Mishra 26 GHz 5-port reproduction: pointing error 0–11°, port-5 beam levels matched to ≤0.1 dB after equalisation (`docs/submission/paper26/comparison.md`); **HFSS packages awaiting simulation**, see below |
| 1-bit / N-bit coded (Phase 15 Part E) | Model-only; cell-pair phase separation is mesh-sensitive | Paper-mode pair at 26 GHz: 173–179° separation across 25–27 GHz; **HFSS check of this exact cell pair not yet done** |
| Level matching / pointing correction (Phase 16 follow-ups) | Self-consistent (one measurement method used everywhere, Phase 16 live-review round) | 4-beam single port: level spread 6.4 dB → 0.5 dB in 2 passes; pointing miss converges to <0.3° on a fine 0.1° local grid |
| Polarization efficiency (scalar, isotropic unit cells) | Documented physical ceiling, not a bug | ≈ 50% co-pol power fraction for a linearly-polarized compensated beam (`docs/history/PHASE15_POLARIZATION_VORTEX.md` "Scalar, isotropic unit cells cap achievable XPD") |

## Known limitations

- **Cell-count floor bug fixed (Phase 21 C1).** `fields.reference_field_grid._cell_centers` floored
  300.8 mm / 6.4 mm to 46 cells (46.999… in floating point). The earlier `four_port_1bit_300mm_*`
  boards were affected: 46 cells per port were stretched onto 47 in the combined board. A
  whole-cell side now keeps all its cells; `test_pattern_beam_polarization_metrics`'s two-beam
  fixture had relied on the same float noise and now builds its 61x61 grid explicitly.
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

-1. **Phase 23 (hexagons + sidelobe optimiser)**: run `hfss_exports/four_port_hex150_*` in HFSS (one port at a time first; the
   `*_opt_all*` boards also with the port phases in their README from one phase-locked source) and compare the SLL/directivity
   tables in each README with the held-out model numbers; the equal-area `four_port_square243_gap0` is the square control. Open
   follow-ups: all-ports absolute pointing is 0.2–0.7° (not ≤ 0.3°) because the ports' sidelobes interfere in the beams;
   the optimiser moves cells only (a pointing penalty term or a Phase 16 pass on the all-ports pattern would fix it).

0. **Phase 19 (nylon lens)**: decide whether ±30.00° is a hard specification (the error bar puts the beam at ≈30.3° median in the primary model and 0.4° typical pointing error; Parts E/F wait on this); run `hfss_exports/nylon_lens_28ghz_dual30/lens.step` in HFSS and report beam angle/gain/level difference, and optionally the six single cells in `hfss_exports/nylon_cells_28ghz/`; then decide Part E (app integration) and close out (tag `phase19-complete`, move the brief to `docs/history/briefs/`).

1. **HFSS checks awaiting simulation** (packages committed, `*.step` files present locally but
   gitignored — see `hfss_exports/README.md`):
   - `hfss_exports/three_port_25deg/` — 1×3 port board, 25° RHCP beams.
   - `hfss_exports/four_port_1bit_gap0/` and `four_port_1bit_gap28/` — 4-port 1-bit coded boards,
     two guard-gap variants.
2. **Unit-cell check of the 10.54 GHz 1-bit pair** used by the walk-through board (the two-state
   cell pair's own phase separation/|S11| match at that specific frequency, not yet run through
   HFSS single-cell comparison the way G1–G6 were for the continuous library).
3. **2-bit multi-port wiring** (parked, see below) — if a future phase needs > 1-bit coding on a
   multi-port layout.

## Hygiene (Phase 17 §4 findings — report only, not acted on)

- **Dead code**: an import-graph scan (`ast`-based, resolves relative imports, checked against
  `src/`+`tests/`+`scripts/`+`examples/`) found no module under `src/` that is both never imported
  and never covered by a test. None found, none removed.
- **Lint debt**: `ruff check src tests examples scripts` has ~95 pre-existing findings predating
  Phase 17 (mostly `B905` bare `zip()`, `F841` unused locals, `B023`/`B007` loop-variable closures,
  a handful of `SIM105`/`SIM115`/`E731`/`E701`/`B017`/`B018`/`I001`/`F401`) spread across `src/`
  (40), `tests/` (42) and `scripts/` (12). Each needs per-site judgement (e.g. adding `zip(...,
  strict=True)` can change behaviour if lengths ever mismatch) rather than a bulk `--fix`, so it's
  left as a punch list for the user to triage rather than silently edited. Only genuinely safe,
  zero-behaviour-change items were fixed as part of this phase: 34 intentional `E402` cases
  (`pytest.importorskip`/`sys.path.insert` guard patterns — annotated `# noqa: E402`, or two stray
  imports moved to the top of their file) and a repo-wide `ruff format` pass (whitespace only).
- No stray debug `print()` or unowned `TODO`/`FIXME`/`XXX` found anywhere in `src/`, `tests/`,
  `scripts/`, or `examples/`.

## Parked ideas (`docs/briefs/_parked/`)

- `PHASE16_UNIT_CELL_SPARAM_TOOL.md` — a unit-cell/FSS tool with plane-defined S-parameters (the
  measurement-plane technique `scripts/adhoc/fss/plane_s.py` prototyped for the
  square-loop FSS reference, `hfss_exports/fss_square_loop_10mm/`). Not implemented.
- Also parked, no brief written yet: transmitarray support in the app (Phase 19 Part E), oblique incidence in the MoM core, tensor
  (anisotropic) unit cells (would relax the scalar-cell XPD ceiling above), 2-bit multi-port.
