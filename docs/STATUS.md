# Status — single source of truth for *now*

**Last updated:** 2026-10-07, by Claude (Phase 17 handoff/cleanup). **Rule: every phase updates this
file** (last release tag/commit, capability matrix, open next steps) as part of its own final commit.

**Last release tag:** `phase16-complete` (commit `1c5ebfb`). The submission snapshot tag
`submission-2026-10-06` is a separate, frozen checkpoint — never move it.

## Capability matrix

| Capability | Validation level | Numbers / where |
|---|---|---|
| Periodic unit-cell solver (ungrounded, Phases 1–9) | HFSS-validated, off-resonance | Case G1: mean phase diff **-0.29°**, max \|S\| diff **~1.6e-5** (`docs/history/PHASE13_WEBAPP.md`) |
| Periodic unit-cell solver, at/near resonance | Mesh-convergence-characterized, residual unexplained | G4 (12.5 GHz): HFSS-vs-0.2mm-mesh diff 13.51° = mesh 5.93° + thickness 3.66° + **residual 3.92° unexplained** (`docs/history/PHASE13_WEBAPP.md`) |
| Grounded supercell / local periodicity (Phase 12C/12D) | Cross-checked against real full-array MoM | Lookup-table-vs-measured per-order power differs by 0.17–1.8 percentage points (`docs/history/PHASE12D_LOCAL_PERIODICITY.md`) |
| Vortex↔pencil physics (Phase 15 Part A/B) | Independently cross-checked (this project's own vector far-field model, not HFSS) | RHCP/LHCP pencil XPD **21.4 dB**; HP/VP pencil XPD 62–341 dB (scalar-cell ceiling); deliberate OAM request: winding **−1 matches charge −1**, null **−31 dB** (`docs/history/PHASE15_POLARIZATION_VORTEX.md`) |
| Multi-port sub-apertures (Phase 15 Part C gap-fill) | Model-only; reproduces Mishra et al.'s *target directions*, not yet checked against that paper's measured isolation/gain | Mishra 26 GHz 5-port reproduction: pointing error 0–11°, port-5 beam levels matched to ≤0.1 dB after equalisation (`docs/submission/paper26/comparison.md`); **HFSS packages awaiting simulation**, see below |
| 1-bit / N-bit coded (Phase 15 Part E) | Model-only; cell-pair phase separation is mesh-sensitive | Paper-mode pair at 26 GHz: 173–179° separation across 25–27 GHz; **HFSS check of this exact cell pair not yet done** |
| Level matching / pointing correction (Phase 16 follow-ups) | Self-consistent (one measurement method used everywhere, Phase 16 live-review round) | 4-beam single port: level spread 6.4 dB → 0.5 dB in 2 passes; pointing miss converges to <0.3° on a fine 0.1° local grid |
| Polarization efficiency (scalar, isotropic unit cells) | Documented physical ceiling, not a bug | ≈ 50% co-pol power fraction for a linearly-polarized compensated beam (`docs/history/PHASE15_POLARIZATION_VORTEX.md` "Scalar, isotropic unit cells cap achievable XPD") |

## Known limitations

- Normal incidence only; single substrate layer; PEC zero-thickness metal model (a real HFSS
  patch has finite thickness — see the G4/G6 residual above).
- No transmitarray support, no finite-array (N×N, N not a power of a small factor) full-wave MoM
  beyond what `Verify` already runs for a cost-bounded N.
- 2-bit (or higher) multi-port coded designs are **not wired** — Phase 15 Part E's coded method
  does not combine with Phase 15 Part C's multi-port layout in the app (single-aperture coded
  only; see `src/app/static/index.html`'s own note hiding Coded when "Multiple ports" is chosen).
- Geometry inversion and full-array STEP export remain `GEOMETRY_INVERSION_LOCKED` /
  `STEP_EXPORT_LOCKED` scientific-finding records (Phase 12D) even though the app produces designs
  and exports anyway, per the project's Phase 13 decision — every such output still carries these
  labels.
- No run manifest / diagnostic-bundle download wired into the app yet (`src/manifest.py` exists and
  is tested, nothing attaches it to a job — `docs/ARCHITECTURE.md`).

## Open next steps

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

## Parked ideas (`docs/briefs/_parked/`)

- `PHASE16_UNIT_CELL_SPARAM_TOOL.md` — a unit-cell/FSS tool with plane-defined S-parameters (the
  measurement-plane technique `results/submission_work/fss/plane_s.py` prototyped for the
  square-loop FSS reference, `hfss_exports/fss_square_loop_10mm/`). Not implemented.
- Also parked, no brief written yet: transmitarray support, oblique incidence, tensor
  (anisotropic) unit cells (would relax the scalar-cell XPD ceiling above), 2-bit multi-port.
