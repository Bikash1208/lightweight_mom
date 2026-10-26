# Phase 20 / R3 — 2-bit PIN-diode RIS cells, linear and circular polarisation

**Status:** active brief, written 2026-10-19. Branch `phase20-ris` (after R2, `286e889`). Read
`docs/briefs/PHASE20_RIS_ROADMAP.md`, `docs/history/briefs/PHASE20_R2_DIODE_MODELS.md`, `results/phase20_r2/DEMO.md` and
STATUS open steps 0b/0c first. All repo rules apply (one commit per item, `make test-fast` + snapshot gates green,
physics tests, honest labels, additive only, no ML in the solver core).

## Decision from the user (2026-10-19)
**Polarisation is a user choice: the RIS must support both linear (LP) and circular (CP).** So this phase designs
**two cell families** on the same board technology, and the app (R5) will let the user pick one. Record this in the
roadmap's "Design changes" log (polarisation: single linear → LP and CP as user-selectable families).

## Fixed for both families (starting values, may change per the roadmap rules)
RT/duroid 5880 1.575 mm, PEC ground, 10 GHz centre, X-band characterisation 8–12 GHz, period 10 mm (allowed 6–15 mm),
SMP1345-079LF models from R2 (`PROVISIONAL_DATASHEET`), diodes with `transverse_conducts=False`, metal pads around each
diode chain as in R2 item 3, centre/bias vias omitted and labelled `VIA_OMITTED_APPROXIMATION`, bias layer behind the
ground out of the RF model. Every response is the full 2×2 reflection (TE/TM, R1) at θ = 0°, 30°, 60°.

## Targets (starting values; report misses, do not hide them)
| | LP family | CP family |
|---|---|---|
| States | 4 (2-bit) for one linear pol | 4 (2-bit) for one handedness |
| Phases at 10 GHz | 0/90/180/270° ± 15° (co-pol) | 0/90/180/270° ± 15° (handedness-preserving term) |
| Magnitude | each state ≥ −1.5 dB, spread ≤ 1 dB | handedness-preserving ≥ −1.5 dB, spread ≤ 1 dB |
| Polarisation purity | cross-pol ≤ −15 dB | opposite-handedness term ≤ −15 dB |
| Angle | targets met at θ = 30°; report θ = 60° | same |
| Bandwidth | report the band where steps stay 90° ± 20° | same |
| **Diode tolerance (R2 finding)** | ±20 % C_T and ±20 % L_s one at a time: step error ≤ 20° | same |

The tolerance row is a hard objective, not a report item: R2 found the OFF state near its 13.8 GHz series
resonance, moving 20–100° for ±20 % C_T or L_s. A cell that meets the phases only at nominal values fails R3.

**CP convention:** define handedness relative to each wave's own propagation direction and state it in
`docs/theory/CONVENTIONS.md`. A plain metal plate flips handedness on reflection; a CP RIS normally needs the
handedness-preserving response (anisotropic, |r_xx| ≈ |r_yy|, phase difference ≈ 180°), with the 2-bit phase from the
geometric (Pancharatnam–Berry) phase: switching the effective element orientation by 0/45/90/135° gives 0/90/180/270°.
Build the circular-basis 2×2 from the R1 TE/TM matrix with a test (a rotated half-wave-plate-like cell gives phase =
2 × rotation).

## Work items (one commit each; stop points ■)

**1. Literature and topology selection.** For each family, find published 2-bit (or 1-bit extendable) PIN-diode RIS /
reconfigurable reflectarray cells at X-band (cite by DOI, no PDFs committed): LP (e.g. asymmetric two-diode patches,
patch + delay-line, stacked designs) and CP (e.g. rotated/switchable split rings, switchable cross dipoles, diodes
selecting element orientation). For each candidate: diode count per cell, reported phases/loss/bandwidth, whether it
fits a planar single-layer MoM model (no vertical currents beyond labelled vias). Pick 2 topologies per family.
**■ Stop and report the four candidates before optimising.**

**2. Cell geometries in the solver.** Composite geometries (R1 item 3) with diode loads for the chosen topologies,
parametrised; STEP export with `Load_<k>` sheets and pads as in R2. Tests: state symmetry checks (the asymmetric LP cell
must give 4 distinct states; the CP cell's states must be rotations/mirrors as designed), energy conservation.

**3. Optimisation.** Search geometry parameters against the targets, including the tolerance objective (evaluate each
candidate at nominal and at the four ±20 % corners). Compute budget: search on mesh 40; correct finalists with the
observed-order extrapolation from meshes 40/80/120 (R2: observed order 1.2–1.65, not 2; mesh 120 ≈ 15.6 GB per angle;
mesh 160 does not fit in 36 GiB). Keep all mesh choices edge-aligned. Report per family: best cell, all target rows,
and which targets were missed and why.

**4. Locality test (moved forward from R6).** For each final cell: supercell MoM with (a) a periodic 2-bit gradient
(0-1-2-3 repeating, an anomalous-reflection case) and (b) random 2-bit codings (≥ 5 draws, 4×4 or larger), at θ = 0°
and 30°. Compare the Floquet-order amplitudes/phases with the isolated-cell prediction. Gate: anomalous-order level
within 1 dB and 15° of prediction; report the random-coding error distribution. If the gate fails, report it — do not
tune the gate.

**5. HFSS check packages.** One per family (`hfss_exports/phase20_r3_lp_cell/`, `.../phase20_r3_cp_cell/`): STEP with
load sheets, all 4 states, θ = 0° and 30°, 8–12 GHz, predictions first (`mom_reference.csv`, with uncertainty), and
the R2 RLC-boundary notes. Keep the run count small (Floquet port gives both modes per solve).

**6. ■ Close-out.** Register both families in the Phase 27 cell-family registry as `MODEL_ONLY` (discrete-state, 2×2
per state). Update STATUS, CHANGELOG, roadmap log; move this brief to `docs/history/briefs/`. Stop and report.

## Do not
- Build the R4 library (f/θ/φ sweep at production density) or any app/RIS design mode (R5).
- Model the bias network; fit diode values; relax the tolerance objective to make a cell pass.
- Commit reference PDFs or raw Touchstone exports.
