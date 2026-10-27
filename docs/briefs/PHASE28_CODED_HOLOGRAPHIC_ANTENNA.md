# Phase 28 — Monopole-fed coded (PIN-diode reconfigurable) holographic metasurface antenna

**Status:** active brief, written 2026-10-24. Branch `phase28-coded-holographic` from `phase20-ris` after the curation
merge. All repo rules apply (`CLAUDE.md`): one commit per item, fast suite green, physics tests, honest labels, no
silent solver changes (additive, bit-identical for every existing case), no ML in the solver core, heavy compute only
on the server via `jobs/`.

## Decision from the user (2026-10-24)
The project focuses on the **monopole-fed coded holographic metasurface antenna**: an integrated feed (the monopole
launching a surface wave across the board) and PIN-diode switched cells whose states rewrite the hologram, so the beam
is steered electronically (the proposal's "digital coded metasurface antenna" with an integrated feed, FPGA-controlled).
**Dropped:** reflectarray / plane-wave RIS designs with an external feed, and transmitarrays.
Record this in `docs/ROADMAP.md` and the roadmap log as a design change.

## What carries over (reuse, don't rebuild)
- Phase 26: monopole launch model (launch fraction 0.62 calibrated at a 6.0 mm pin), leaky-wave forward model,
  tensor-CP synthesis (`src/hologram/tensor_cp_design.py`), the single-hexagon board validated in HFSS
  (23.03 dBi RHCP, `docs/STATUS.md`), the 150 mm hexagon and 6.4 mm lattice.
- Phase 14 / Phase 26 B2.5: eigenmode / dispersion tooling for surface-wave cells.
- Phase 20 R1/R2: lumped loads in the MoM, oblique incidence, diode equivalent circuits (`src/ris/diodes.py`,
  MADP-000907-14020 and SMP1345), the HFSS lumped-RLC diode representation (verified in AEDT on 2026-10-24).
- Infrastructure: cell-family registry, server job pipeline, library + held-out validation method, HFSS calibration hook.

## Parked, not deleted
R5 (plane-wave RIS design mode) and the CP-2 cell: move `PHASE20_R5_RIS_DESIGN_MODE.md` to `docs/briefs/_parked/`
with a status note (R5-1/2/3/5 and the R5-6 exporter done; R5-4 stopped). The RIS designer page and CP-2 registry entry
stay as they are (MODEL_ONLY). Jobs `r5-supercell-gradients`, `r5-supercell-random` and `r4-band-refine`: mark
**WITHDRAWN 2026-10-24 (direction change)** in their job files; keep any finished results with their manifests.

## Work items (one commit each; stop point ■)

**0. Re-plan the repo.** ROADMAP (planned track: deliverables, done, next), STATUS open steps, briefs README, the
parking and job withdrawals above. Docs only.

**1. Literature (cite by DOI, no PDFs committed).** Reconfigurable / coded holographic and modulated-surface
leaky-wave antennas with PIN diodes or varactors, surface-wave or monopole/probe feed, X/Ku band. Per paper: cell
pitch vs λ, diodes per cell, bit depth, switching contrast (surface reactance or β/α change), aperture size, gain,
scan range, polarisation, bias approach, measured vs simulated. Short table plus 5–10 lines of conclusions.

**2. Feasibility model (no new cell yet).** Using the Phase 26 forward model with *idealised* switchable cells
(each state = a surface-impedance / tensor value; no geometry):
- Bit depth: 1-bit and 2-bit quantisation of the holographic modulation (scalar for LP; tensor/rotation for CP).
- Pitch: 4, 5, 6.4, 8 mm.
- Aperture: the 150 mm hexagon with the monopole at its centre (and one larger size for reference).
- Steering: θ = 0–50°, several φ; frequency **10.0 GHz (PI decision 2026-10-26; the 10.54 GHz Phase 26 results are history)**.
- Report per case: gain, SLL, pointing error, CP purity (AR / XPD) where applicable, quantisation loss vs continuous,
  scan loss, and the **required switching contrast** (the minimum Δ between states that keeps the gain within 1 dB of
  the continuous hologram).
- Practicality per case: diode count, independent bias lines (and if lines can be shared, e.g. row/column or
  per-ring), and the FPGA/shift-register count for a Basys 3 class board.
Heavy sweeps go to the server as a job. Label everything MODEL_ONLY (idealised cells).
**■ Stop and report:** a recommended pitch, bit depth, polarisation scope, and diode count, with the trade-offs. No
cell design before the user approves.

**3 (after approval). Cell candidates.** Diode-loaded subwavelength cells reaching item 2's required contrast at the
chosen pitch, characterised by dispersion (eigenmode) with diode loads in each state, including the ±20 % diode
tolerance objective. If the eigenmode solver needs lumped loads, that is an additive solver extension (its own commit,
bit-identical for existing cases, with physics tests). HFSS eigenmode check package per candidate, predictions first.

**4 (after 3). State library** (dispersion per state vs f) on the server, held-out validation as in R4.

**5 (after 4). Antenna:** coded hologram synthesis per beam (codebook), forward-model predictions, app mode, HFSS
board package (smaller board first), bias-layout sketch.

## Do not
- Design or optimise a cell before the item 2 stop point is approved.
- Delete anything from the RIS track (park and label only).
- Fit surrogates and call them physics; idealised-cell results must say so.

## PARKED 2026-10-27 (Phase 28, branch `phase28-coded-holographic-v2`; the user switched projects, resume later)

**Architecture decided (design level, MODEL_ONLY except where stated):** side/sector-fed surface-wave hologram with a Yagi-type launcher (driven pin + reflector post + directors), SCALAR 1-bit, LINEAR polarisation, ONE MA4AGP907 (= MADP-000907-14020)
PIN diode per cell, individual control line per cell, 9 mm pitch, RT/duroid 5880 3.175 mm, 10.0 GHz, 150 mm hexagon (724 cells, 724 diodes). Diode model `MEASURED_VENDOR_SPARAM` (`src/ris/diodes_measured.py`, extraction
`scripts/research/phase28_diode_extract.py`; R 3.3 / 2.4 / 2.0 ohm at 5 / 10 / 15 mA, C 28.2 fF at 0 V, no series L visible). Lead cell: loop cell R4 (6.6 mm loop, 1.2 mm strip, diode in the middle of one side), **HFSS_CHECKED at ROT 0, normal incidence,
9 - 11 GHz** (one run: all sheet-admittance components within 4.3 % of the MoM, Delta B 4.22 mS against 4.05 mS; `hfss_results/phase28_loop_R4/NOTES.md`; calibration hook `src/ris/loop_cell_calibration.py`). Cells stay UNROTATED (the diode axis along x fixes a clean
x polarisation; model: single beam, SLL -14 to -21 dB, 21 - 23 dBi with a placeholder launch fraction; 1-bit sampling limits the scan toward the feed to about 30 deg). Earlier cell tracks (G, J, two-diode searches) are superseded; their
files stay as history.

**User note (2026-10-27 HFSS run):** the state files were swapped in the project: `OFF.csv` carries the R = 2.4 ohm state (lossy, isotropic) and `ON.csv` the C = 28.2 fF state (open, anisotropic). Which of the files or the boundaries was
swapped is still TO BE CONFIRMED by the user. (The first HFSS build used hand-set boundaries; `build_cell.py` is fixed since 9ea6897 and tested against an AEDT recording stub, but not yet re-run in AEDT.)

**State of the open items, and the next steps in order:**
1. **Rotated-load solver bug** (`transverse_conducts=False` on a rotated load gives ON = OFF). Cause found: the documented R3 limit (the 1e12 ohm transverse penalty locks the Cartesian rooftops on a rotated region). Done: a `UserWarning`
   in `check_load_geometry` (results bit-identical, fast suite 2900 passed), a minimal reproduction and a rotation-covariance physics test (conducting rotated loads: 2 - 7 % at 30 / 45 deg). NOT done: a real fix of blocked rotated loads
   (would need a basis aligned with the load or a projected constraint); not needed for the unrotated lead cell. Rotated-cell MoM results are good to 10 - 17 % only (staircase of the rotated geometry).
2. **Launcher design + HFSS package.** Model (coupled grounded posts through the TM0 wave, two feed calibrations) and package `hfss_exports/phase28_launcher/` exist with predictions (width 77 - 79 deg, S11 -8.7 / -15.9 dB, F/B 5 - 8 dB, predictions only).
   NOT done: any HFSS run; the model neglects near-field/space-wave mutual impedances and the loaded-surface edge. Next: user's HFSS run (bare slab first, then OFF-state loop surface), then re-tune.
3. **Phase 28 item 5:** app design mode (side-fed scalar 1-bit hologram synthesis, codebook per beam, forward-model predictions) and the full-board HFSS export (150 mm hexagon, 724 cells; smaller board first).
4. **Full-board HFSS check** (user): beam angle, gain, SLL, cross-pol against the forward model.
5. **Oblique / rotated-cell HFSS checks:** polarisation turned by 20 / 40 deg in the Floquet cell (the directionality check of `hfss_exports/phase28_loop_R4/`), the rotated cell, the measured diode reactance, package inductance.

**Scope:** design and solver/tool work verified in HFSS. NO fabrication, bias-network or FPGA work. Nothing here is merged to main, tagged or deleted; branch state is the last commit on `phase28-coded-holographic-v2`.
