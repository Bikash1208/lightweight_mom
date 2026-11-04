# Grouped coded holographic metasurface antenna — paper notes (living document)

## Headline claim (draft)
A flat, single-port, surface-wave-fed holographic metasurface antenna that switches between K predefined
beam/polarisation states (demo: RHCP (25°,0°) ↔ LHCP (25°,180°); next K = 4 with four polarisations), placing
diodes only where the state holograms differ and driving them with a minimal number of control lines
(K = 2: one line + inverter per diode axis).

## Contributions (candidate)
1. State-difference placement: cells whose 1-bit code is identical in all K states are static (no diode).
2. Code-word grouping: diode cells grouped by their K-bit state word → ≤ 2^K − 2 lines (halved with inverters);
   K = 2 needs one line per axis. (Present as a design rule supporting the main claim, not as the headline.)
3. Beam + handedness switching from one feed with crossed two-diode tensor cells (0°/90° diodes only).
4. Efficiency at fixed diode/line count: analog static cells, pre-phased switched cells, radial leakage taper,
   launcher optimisation, diode-loss-aware cell design, edge rings, joint state-phase optimisation (I1–I7).
5. (Method/solver, optional) MoM model of diode-loaded tensor cells incl. the rotated-load fix, HFSS-validated.

## Novelty check (2026-11-04, journals; full text of nearest work still to be read)
Published pieces (cite, do not claim): 1-bit/binary reconfigurable holographic antennas (Hwang, Sci. Rep. 2020,
doi 10.1038/s41598-020-65458-3; Zhao, Pan, Shen, IEEE TAP 2026, doi 10.1109/TAP.2025.3621516); varactor-tuned
multibeam holographic surface (Kampouridou & Feresidis, IEEE AWPL 2022, doi 10.1109/LAWP.2022.3192977);
static + reconfigurable hybrid cells (sparse RIS, Zhao, Cai, Wang, Wu, Luo, Cheng, Cui — reflect-type, GA
placement; Sheng et al., Micromachines 2023, doi 10.3390/mi14081631); polarisation-reconfigurable holographic
antenna (Li, Tang, Xiao, IEEE Access 2019); CP from two orthogonal holograms (Fong et al., IEEE TAP 2010).
Not found: simultaneous beam-direction AND CP-handedness switching on a single-port surface-wave holographic
board with crossed two-diode tensor cells; state-difference static/switched split applied to surface-wave
holography; code-word grouping to one control line. Residual risk: Zhao/Pan/Shen TAP 2026 full text unread;
IEEE Xplore search ("reconfigurable holographic" AND "circular polarization", 2023–) still to do.
(These DOIs are as supplied by the project owner; none was opened or verified from this container.)

## Positioning vs alternatives
- vs reflect-type RIS: integrated feed (no horn at a focal distance), single ~2 mm board, no spill-over or
  blockage, far fewer diodes and control lines, polarisation switching built in. RIS wins on arbitrary
  continuous steering and (today) aperture efficiency → report gain per volume and gain per diode, and close
  the efficiency gap with I1–I7.
- vs feed switching (multi-probe + SPDT/SP4T): every state uses the full aperture, no extra feeds/RF switch, K
  can grow without new feeds. Feed switching wins on diode count (2–4) and currently ~2 dB gain → must be in the
  comparison table.
- vs full 1-bit coded hologram (diode in every cell): same capability with ~27–50 % fewer diodes and ~1 line.

## Validation plan (for a Transactions paper)
Unit cell (4 diode states) MoM vs HFSS eigenmode → K = 2 board MoM vs HFSS → K = 4 design → comparison table and
analyses (loss budget, bandwidth, diode/etch/εr tolerance, limitations) → fabricated prototype measured
(S11, patterns, gain, AR, efficiency, both/all states). Without K = 4 and measurement, target IEEE AWPL instead.

## Results log
(append per item: item | gain A/B | SLL | AR | aperture eff. | edge % | static % | diodes | lines | commit)

All rows MODEL_ONLY (no HFSS). Cumulative chain; I7 selected eps_B = 0, g = 2.7 (final = I5 row: 19.38 dBi, 759 diodes). Detail: `results/sidequest_hybrid/CODED_IMPROVE.md`.

| item | gain A/B [dBi] | SLL [dB] | AR [dB] | aperture eff. | edge [%] | static [%] | diodes | lines | commit |
|---|---|---|---|---|---|---|---|---|---|
| B0 (rebased) | 17.22 / 17.22 | -13.7 | 0.84 | 5.81 % | 82.5 | 73.2 | 766 | 1 | see git log |
| I1 | 17.45 / 17.45 | -13.6 | 0.70 | 6.13 % | 82.7 | 73.7 | 752 | 1 | see git log |
| I2 | 17.47 / 17.47 | -13.5 | 0.67 | 6.15 % | 82.9 | 73.3 | 764 | 1 | see git log |
| I3 (not adopted: taper saturated) | 17.47 / 17.47 | -13.5 | 0.67 | 6.15 % | 82.9 | 73.3 | 764 | 1 | see git log |
| I4 (flagged: S11 not met) | 18.53 / 18.53 | -13.5 | 0.67 | 7.85 % | 82.9 | 73.3 | 764 | 1 | see git log |
| I5 | 19.38 / 19.38 | -14.6 | 0.87 | 9.55 % | 79.1 | 73.5 | 759 | 1 | see git log |
| I6 (not adopted) | 19.13 / 19.13 | -15.1 | 0.80 | 9.02 % | 79.5 | 73.4 | 760 | 1 | see git log |

Findings for the paper: (1) the board is leakage-limited: the crossed rim cell has alpha_beam 0.7 Np/m (mesh-converged estimate), 80 % of the surface-wave power reaches the edge, so taper, edge rings and analog statics cannot move gain; (2) the lever that works is the cell: diode position 0.8 mm raises the swing and cuts the ON conductance (+0.85 dB); (3) launcher: launch fraction 0.80 needs a 10 mm top-hat pin whose S11 is not met in the model; (4) honest baseline is 17.2 dBi (aperture efficiency 5.8 %) not the 18.2 dBi of the mesh-64 cell; (5) state-difference placement holds at 73 % static, 766 diodes, 1 line + inverter per axis. HFSS validation of the rim cell (4 states, `hfss_exports/sidequest_coded_cross_cell/`) is prepared and not sent.
