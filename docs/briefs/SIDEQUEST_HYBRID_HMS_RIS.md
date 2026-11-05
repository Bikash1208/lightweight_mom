# Side quest — hybrid HMS + sparse RIS antenna (quadrant switching)

Branch `sidequest-hybrid-hms-ris` (from `phase30-new-cell`, 474218e). Cloud container: no HFSS, no Windows, limited
CPU. Coarse MoM runs only (< ~10 min each); anything heavier goes to `jobs/` as a job file for the Windows server.
Scope: antenna + unit-cell design and solver work only (no bias-network layout, fabrication, FPGA).
Label for everything below: **MODEL_ONLY** until HFSS-checked.

## Concept
One single-port holographic metasurface antenna (HMS), surface-wave fed by one probe, aperture split into **2 or 4
quadrants**. Default state (all switchable cells OFF): one quadrant radiates a beam at (θ1, φ1), polarisation P1.
Switching PIN-diode (RIS) cells ON activates another quadrant: beam (θ2, φ2), polarisation P2. Polarisation per state
is configurable from {LP-H, LP-V, RHCP, LHCP}. Challenge: switchable cells placed **sparsely, only where needed**, to
keep the diode count small.

Reference (treat as preprint unless a journal DOI is verified): Zhao, Cai, Wang, Wu, Luo, Cheng, Cui, "A Sparse
Reconfigurable Intelligent Surface for Single-Beam Scanning" (Southeast Univ.): 16x16 reflect RIS, 9.5 GHz, 14 mm
cells, 66 % reconfigurable 1-bit cells + 3 static cell types, nested-GA placement, ±50° dual-pol scan, 340 switches.
Taken over: the sparse reconfigurable + static cell + GA placement idea. Not transferable: the plane-wave reflect
geometry (our excitation is a surface wave).

## Architectures (evaluate at model level, then pick one)
- **A. Surface-wave gating.** Passive static holograms in every quadrant, each for its own beam + polarisation. Diodes
  only in narrow gate strips on quadrant borders / around the feed. Gate cell: *pass* (supports the TM surface wave,
  matched to the hologram's mean impedance) vs *block* (stopband / strong mismatch). Open questions: isolation vs
  strip width (rows) and insertion loss; fate of blocked power (reflection to feed -> S11, or redirection into the
  open quadrant); feed position (centre vs offset) and quadrant geometry; can default be "all diodes OFF".
- **B. Sparse switchable hologram.** Quadrant holograms partly made of diode cells toggling between "radiating
  (hologram value)" and "non-radiating (average impedance, bound wave passes)". GA-optimised sparse placement.
  More diodes, but allows mixed beams/polarisations.
- **C. Hybrid A + B** if neither alone meets targets.
Quadrant count (2 vs 4) chosen from item-2 results.

## Baseline to reuse
- Phase 26 tensor cell, RT5880 1.575 mm (εr 2.2, tanδ 0.0009), 6.4 mm lattice, 10.54 GHz, monopole probe (launch
  fraction 0.62 calibrated); single-feed hexagon HFSS-validated at 23.03 dBi RHCP (`docs/STATUS.md`, Phase 26 C/D).
- Diode: MACOM MADP-000907-14020 / MA4AGP907 measured model (Phase 28; ON ~2-3.3 Ω series, OFF ~28 fF),
  `src/ris/diodes.py`; HFSS lumped-RLC representation verified in Phase 28.
- **Known solver bug (Phase 28, parked): rotated lumped loads give ON = OFF at non-90° rotations.** Policy: either fix
  first (additive, bit-identical, physics tests, own commit) or restrict diode cells to 0°/90° and state so.
- Phase 30 RO3010 cell: leakage unresolved in HFSS -> **not used here**.

## Work items (one commit each; ■ = stop and report)
0. Brief + plan (this file). Docs only.
1. Literature, journals only, DOIs verified, no conference/arXiv: reconfigurable holographic / modulated-surface
   antennas with PIN diodes; surface-wave gating (switchable EBG / stopband cells); sparse reconfigurable
   metasurfaces; sector-switched leaky-wave / holographic antennas. Per paper: frequency, cell, diodes/cell,
   isolation/contrast, aperture, gain, beam states, polarisation, measured vs simulated. Short table + ≤ 8 lines.
2. Architecture study (forward model, idealised cells). 150 mm hexagon (+ circle comparison), 10.54 GHz, single probe;
   A and B with 2 and 4 quadrants. Per state: gain, direction, SLL, AR/XPD, aperture efficiency, off-quadrant power,
   power reflected to feed, **required gate isolation / contrast for off-quadrant leakage ≤ -15 dB**, diode count.
   Targets: 2 states (25°, 0°)/(25°, 180°); 4 states θ = 25°, φ = 0/90/180/270°; one mixed polarisation set
   (RHCP/LHCP/LP-H/LP-V). ■ Report: architecture, quadrant count, gate width / sparse fraction, diode count,
   required contrast. No cell design before approval.
3. (after approval) Switchable cell design on the chosen stack: gate cell (A) and/or switchable hologram cell (B),
   one diode per cell, 0°/90° orientation unless the bug is fixed. MoM dispersion ON and OFF with the measured diode
   model, ±20 % tolerance, insertion loss and isolation of a 2-4-row gate supercell.
4. (after 3) Placement optimisation: GA / deterministic placement of diode cells or gate strips, co-optimised with
   static holograms per state, using the forward model. Per-state patterns.
5. (after 4) HFSS package: STEP files only (no AEDT .py), one per state or one board with named `D_ON`/`D_OFF`
   lumped sheets, port table, `HFSS_SETUP.md` (≤ 1 page), predictions JSON. Re-open every STEP and check body counts.

## Rules
One commit per item; fast suite green; solver changes additive and bit-identical (`tests/_platform_tol.py`); no PDFs
or vendor `.s2p` committed; never move/delete tags; delete nothing without approval; utf-8 + LF; terse, coarse grids
first, no unrequested sweeps, detail in repo files not chat.

## Status
- [x] 0 brief
- [~] 1 literature: 8 user-supplied DOIs recorded, BIBLIOGRAPHIC_ONLY (`results/sidequest_hybrid/LITERATURE.md`; doi.org / crossref not reachable)
- [x] 2 architecture study ■ (`results/sidequest_hybrid/ARCH_STUDY.md`; awaiting approval)
- [x] 1b coherent gate-ring reflection (`results/sidequest_hybrid/GATE_REFLECTION.md`)
- [~] 3 gate cell (`results/sidequest_hybrid/GATE_CELL.md`): 16 dB (not 20 dB) isolation, pass state mismatched; stop and report
- [x] round 3 (`results/sidequest_hybrid/GATE_RETUNE_MATCH.md`): retune, transition-row match (RL 22 dB), ring with real pass state; **correction: the ON gate blocks only ~3 dB once the damped quasi-TEM root is included**
- [x] 5a HFSS gate test package `hfss_exports/sidequest_gate_test/` (STEP, 35 sites, D_ON/D_OFF groups, predictions both states)
- [ ] 4 placement (blocked until a gate that blocks exists)   - [ ] 5b full-antenna HFSS package

Approved by the user: architecture A, 2 quadrants, centre feed, 3 gate rows, isolation target 20 dB per strip; 4 quadrants later with the same gate cell.

- Round 4: dog-bone test board ARCHIVED (negative result, not run). Mushroom-type via-diode gate: `results/sidequest_hybrid/MUSHROOM_CELL.md`, package `hfss_exports/sidequest_mushroom_cell/` (eigenmode, ON/OFF, PEC lid 30 mm); thin-wire via MoM extension scoped (planning only).
- Round 5: HFSS mushroom baseline returned (OFF ok, ON failed, resonance 2.2x low); model recalibrated (`results/sidequest_hybrid/MUSHROOM_CALIBRATION.md`); four ON variants + two OFF twins in `hfss_exports/sidequest_mushroom_variants/` (user runs phases 90 and 180 only). Via-MoM extension parked.
- Round 6 (new direction, mushroom gate parked): switched pin launcher (method 1), two-probe fallback (method 2), polarisation note (method 3): `results/sidequest_hybrid/LAUNCHER_STUDY.md`. HFSS launcher test only after approval.
- Round 7: grouped coded 1-bit hologram, K = 2 and 4 (model level): `results/sidequest_hybrid/CODED_HOLOGRAM.md`.
- Round 8 (decision: grouped coded hologram K = 2 is the side-quest design; launcher and gates parked): (1) Phase 28 rotated-load bug fixed (cause: `transverse_conducts=False` default on rotated diode loads freezes the footprint to an open circuit; `loads.default_transverse_conducts`, `tests/test_loads_rotated_diode.py`); (2) K = 2 RHCP (25, 0) / LHCP (25, 180) with rotated 1-diode cells vs 2 diodes per cell: `results/sidequest_hybrid/CODED_POLARISATION.md`; (3) HFSS unit-cell package for the 1-bit pair prepared, NOT SENT: `hfss_exports/sidequest_coded_cell/` (ON/OFF, eigenmode at 90 and 120 deg, PEC lid, EXACT_DIMS PASS).
- Round 9: crossed two-diode rim cell solved and chosen (`CODED_POLARISATION.md` sec. 6); route confirmed: crossed rim cell, two diodes per cell, K = 2 RHCP (25, 0) / LHCP (25, 180), 1 control line + inverter. Fast suite in a venv from `requirements-lock.txt` (fastapi, pytest 9.1.1): 2901 passed, 2 failed, both environmental: `tests/test_app_main.py::test_version_endpoint_reports_the_installed_package_version` needs `pip install -e .` (package metadata for `lwmom`); `tests/test_tensor_cp_service.py::test_app_export_of_the_single_hexagon_matches_the_committed_phase26_folder` fails on a `layout.png` byte difference (rendering-library version issue, not investigated). `test_hologram_ports*` failures were a pytest-version issue (pass with the locked pytest 9.1.1; fail with 9.2.0).
- Round 11 (coded HMS efficiency items I1-I8): `results/sidequest_hybrid/CODED_IMPROVE.md`, paper notes `docs/paper/CODED_HMS_PAPER_NOTES.md`. Result (MODEL_ONLY): 17.2 dBi baseline (rebased) -> 19.4 dBi (18.3 without the flagged launcher), 759 diodes, 1 line + inverter per axis; edge power stays ~80 % (cell alpha 0.7 Np/m); taper and edge rings do not help.
- Round 12 (cell search before K = 4): `results/sidequest_hybrid/CODED_CELL_SEARCH.md`. Correction: cell tensor at kt = beta raises the RT5880 alpha to 1.3-1.9 Np/m (board 21.5 dBi); RO3010 superpixel is out of the first-order regime and gives 1.15x (not 2x) the RT5880 aperture efficiency -> RT5880 retained; leakage alpha is the limiting parameter. Top-hat launcher: HFSS list only.
- Round 13: RT5880 cell HFSS package `hfss_exports/sidequest_rt5880_cell/` (4 states, eigenmode 90 deg x/y, alpha back-out script); K = 4 model `results/sidequest_hybrid/CODED_K4.md`; paper notes updated.
- Round 14 (model only, HFSS CSV pending): K = 4 minimum control lines (3 per axis, -0.56 dB) and a real 4-probe quadrant model (12.9/12.5 dBi before SP4T loss); `results/sidequest_hybrid/CODED_K4.md`, paper notes updated.
- Round 15 (model only): fewer-diodes study D1-D6, K = 2 and K = 4: `results/sidequest_hybrid/DIODE_REDUCTION.md` (Pareto figure, tables); no design adopted; best <= 500 diodes: K = 2 484 diodes 20.09 dBi; K = 4 none (515 diodes 18.0 dBi).
- Round 16: HFSS no-switching diagnosed (patch imported as vacuum, sheets touching); package re-exported (single STEP, named bodies, geometry checks); results on hold: `results/sidequest_hybrid/CELL_VALIDATION.md`.
- Round 17: HFSS (corrected setup) shows the crossed rim-disc cell switches only with both diodes ON; MoM analysis, 4-diode independent-axes candidate + HFSS package `hfss_exports/sidequest_rt5880_cell_4diode/`, scalar fallback; EXPORT STANDARD `hfss_exports/EXPORT_STANDARD.md`; boundary audit `EIGENMODE_BOUNDARY_AUDIT.md`; all coded results on hold.
