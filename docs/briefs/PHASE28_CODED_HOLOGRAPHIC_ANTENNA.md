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
- Steering: θ = 0–50°, several φ; frequency 10.0 and 10.5 GHz (the PI decision is open; keep both).
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
