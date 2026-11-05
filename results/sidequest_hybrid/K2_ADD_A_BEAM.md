# Add-a-beam candidate (beam 1 always on, beam 2 switched) on the r 2.30 RO3006 board, state-dependent beta - MODEL_ONLY (2026-11-06); no package
State 1 (default, diodes OFF) = beam 1 at (25, 0) only; state 2 (diodes ON) = beam 1 + beam 2 at (25, 180). Same cells (r 2.30 library, decoupled depth, pairs at the calibrated contrast), <= 514 diodes in one group (one control line), neutral-static penalty off (the statics carry beam 1 in both states), LP-V, launch 0.71, each state designed with its own real beta/k0 (1.33 / 1.26). Gains calibrated with the HFSS offsets (state 1 -1.78 dB, state 2 -2.7 dB). Script `sidequest_ro3006_v2.py` (mode `add`), data `ro3006_v3_add2.json`, `ro3006_add_abc.json`.
**Port to the strong optimiser**: priced diode allocation (lambda bisected to <= 514 diodes), exact wrong-beam sums (beam 2 in state 1), target weights per (state, beam), amplitude profile iterated. Parity check on the beam swap (v2 settings, same cells): generic routine 15.6 / 12.7 dBi, SLL (first-null) -12.3 / -11.9 dB against the strong v2 routine 15.0 / 12.1 dBi, -12.7 / -12.3 dB (0.4 dB apart). SLL = highest peak outside the first null of every wanted beam (2 x HPBW if no null), relative to the stronger beam; "vs weaker" = relative to the weaker beam of state 2.

## Three designs (beam-2 target weight w2 and beam phase eps_B chosen from a 48-point scan)
| | (a) max beam 1, beam 2 as is | (b) balance <= 3 dB | (c) equal beams | reference: beam swap v2 |
|---|---|---|---|---|
| settings | w2 0.3, eps_B 270 | w2 1.25, eps_B 0 | w2 1.7, eps_B 90 | v2 final |
| diodes (one line) | 340 | 500 | 482 | 514 |
| state 1: beam 1 / level at beam 2 [dBi] | 18.4 / -8.3 | 17.1 / 0.2 | 16.8 / 2.7 | 15.0 / -2.2 (wrong beam) |
| state 2: beam 1 / beam 2 [dBi] | 17.7 / 1.2 | 15.2 / 12.7 | 14.3 / 14.3 | 12.1 (beam 2 only) |
| beam-2 balance in state 2 (b1 - b2) | 16.5 dB (beam 2 is not a beam: 1.2 dBi) | 2.5 dB | 0.0 dB | NA |
| beam-2 on/off contrast (b2 state 2 - state 1) | 9.4 dB (not meaningful) | 12.5 dB | 11.6 dB | NA |
| isolation state 1 (beam 1 - level at beam 2) | 26.7 dB | 16.9 dB | 14.1 dB | 17.2 dB |
| SLL state 1 / state 2 (first-null) [dB] | -13.8 / -13.5 (vs weaker beam: 3.0) | -12.6 / -11.9 (vs weaker -9.4) | -11.9 / -10.0 (vs weaker -10.0) | -12.7 / -12.3 |
| margin to break-up (3 dB loss of the weaker beam), beta/k0 negative / positive: state 1 | 0.13 / 0.08 | 0.14 / 0.09 | 0.14 / 0.09 | 0.13 / 0.10 |
| state 2 | 0.09 / 0.09 | 0.09 / 0.09 | 0.08 / 0.09 | 0.12 / 0.12 |
Beam-2 weight against balance (eps_B 0 / 90 / 180 / 270, state 2 balance in dB): w2 0.3: 18 / 18 / 17 / 17; 0.7: 9 / 9 / 9 / 10; 1.0: 7 / 6 / 5 / 7; 1.15: 3.7 / 4.3 / 3.9 / 5.1; 1.25: 2.5 / 3.2 / 3.4 / 2.9; 1.4: 1.3 / 1.7 / 1.9 / 1.2; 1.5-1.7: 0.9 .. 0.0; 2.0: -1.3 .. -2.2; 3.0: -4.7 .. -5.8. Beam 1 in state 2 falls from 17.9 (w2 0.3) to 14.3 (equal beams): -3.6 dB for the beam-2 power.

## Findings
1. **Balance costs beam 1**: in state 2 beam 1 is 17.7 / 15.2 / 14.3 dBi for (a) / (b) / (c); the equal-beam design has each beam about 3.4 dB below the maximum beam 1 (the power is split) and about 2 dB above the swap design's single beam (12.1), because the statics radiate beam 1 in both states.
2. State 1 is better than the swap in every add-a-beam design (16.8-18.4 against 15.0 dBi), with the same diode cap, because the statics are useful in both states.
3. Beam 2 is made by the switched cells alone: 340-500 diodes (170-250 cells, 13-19 % of the cells); the on/off contrast of beam 2 is 11.6-12.5 dB for (b)/(c). The contrast cannot be larger than the isolation allows: the level at beam 2 in state 1 is 0.2 / 2.7 dBi (a diffuse level of the model).
4. SLL: (a) -13.8 / -13.5 dB; (b) -12.6 / -11.9; (c) -11.9 / -10.0 dB relative to the stronger beam; relative to the weaker beam of state 2 the sidelobes are only -9.4 / -10.0 dB below it in (b) / (c) (the sidelobes of the stronger beam reach up to the level of the weaker beam). Under the HFSS gap seen on board_K2_default with this definition (+4.1 dB) expect about 4 dB worse in HFSS.
5. Margins to the beta break-up are 0.08-0.14 (the positive side and state 2 are at 0.08-0.09, slightly below the 0.10 target of the swap design; the expected real beta lies within +-0.02-0.04).
6. Compromise cells: the static cells carry beam 1 in both states; beam 2 in state 1 is suppressed by the wrong-beam penalty; no neutral-static device is used (it does not apply when one beam is common).
Limits: model only; calibration offsets (-1.78 / -2.7 dB) come from board_K2_default; the model ring (cone) is far below the HFSS ring (10-17 dB) and the cell at r 2.30 is not HFSS-checked; state 2 of add-a-beam shares the power of the aperture between two beams, so the gain per beam also depends on the assumed equal loss mechanisms of the two states.
