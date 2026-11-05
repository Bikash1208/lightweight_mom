# r 2.30 board v2 (state-dependent beta) and the add-a-beam variant - MODEL_ONLY (2026-11-06)
v2 package built locally, **NOT synced to `C:\hmsa_sync`** (waits for the state-1 ffd of v1): `hfss_exports/sidequest_board_RO3006_r2p30_v2/` (same format, 1350 cells, 257 switched = 514 diodes, one line, pad, PinLength 5.0 mm; geometry checks and EXACT_DIMS PASS; STEP 74.6 MB gitignored). Design JSON `ro3006_board_design_r2.3_v2.json`.

## 1. What v1 does at the real beta (`sidequest_eval_v1.py`)
The delivered v1 (designed with one beta/k0 1.404): at 1.404 / 1.404: 15.7 / 13.2 dBi, SLL -13.7 / -11.8 dB (21 deg exclusion; -13.6 / -13.2 with the earlier 1.2-HPBW metric); **at the real 1.33 / 1.26: 14.5 / 8.0 dBi, SLL -11.6 / -4.8 dB** (state 2 collapses: wrong-beam level -1.0 dBi).

## 2. v2 beam swap, each state designed with its own beta/k0 (1.33 / 1.26; A0 = mean of the A0 at the two betas = 19.9 Np/m decoupled, 11.2 coupled)
Same ingredients as v1 (neutral statics kappa 2, inverse taper 0.5, edge ring x1.5, <= 514 diodes, one line); the state-dependent beta enters the coupling of every cell to each beam (`bd.unit` per state) and the far field of each state. Scan (12 settings, `ro3006_v2b_swap.json`; gain at the design direction state 1 / 2, SLL (1.2-HPBW metric), isolation):
| settings | gain | SLL | isolation |
|---|---|---|---|
| decoupled g 1.5 mu 30 | 15.9 / 12.8 | -13.4 / -12.9 | 22.3 / 15.9 |
| decoupled g 2.0 mu 30 | 15.1 / 12.6 | -13.0 / -13.2 | 20.6 / 18.6 |
| **decoupled g 2.0 mu 70 (v2)** | **15.0 / 12.1** | **-13.4 / -13.8** (strict 21 deg exclusion: -12.5 / -11.4) | 17.2 / 17.0 |
| decoupled g 3.0 mu 70 | 14.5 / 11.0 | -13.1 / -13.3 | 20.8 / 20.1 |
| coupled g 1.5 mu 30 | 17.7 / 11.8 | -13.5 / -12.6 | 23.1 / 17.8 |
| coupled g 2.0 mu 30 | 17.1 / 11.7 | -13.3 / -12.0 | 26.5 / 17.3 |
**Decoupling kept**: it improves the state balance (decoupled 15.0 / 12.1 = 2.9 dB apart; coupled 17.1 / 11.7 = 5.4 dB apart) at the same SLL. v2: edge 14 %, SW to the beam 78 %; at the real betas the gain is 15.0 / 12.1 against v1's 14.5 / 8.0 dBi (state 2 +4.1 dB).
**Margin to break-up (each state's real beta varied alone, the other fixed; criterion: 3 dB loss of the gain at the design direction or the peak off the design direction by > 3 deg):**
| state | beta/k0 design | negative side | positive side | gain at -0.10 / -0.05 / +0.05 / +0.10 |
|---|---|---|---|---|
| 1 | 1.33 | **-0.13** | **+0.10** | 13.3 / 14.7 / 14.0 / 11.8 dBi |
| 2 | 1.26 | **-0.12** | **+0.12** | 9.9 / 11.6 / 11.5 / 9.6 dBi |
Target >= 0.10 on both sides and both states is met (the estimated real beta lies within +-0.02-0.04 of the design values: state 1 1.312-1.350, state 2 1.248-1.276). Break-up (peak leaves the design direction) starts at about -0.20 (state 1 peak (27, 342)) and +0.15 (state 1 (24, 340); state 2 (27, 170)).
Note on the SLL numbers: the metric of the earlier design runs excludes 1.2 HPBW (about 45 deg) around the beam; with a fixed 21 deg exclusion v2 gives -12.5 / -11.4 dB (v1: -13.7 / -11.8 at its design beta). Neither is validated: the model underestimates the HFSS ring by 10-17 dB.

## 3. Add-a-beam variant (state 1 = beam 1 only; state 2 = beams 1 + 2 at equal target power), model only
Generic joint design routine `sidequest_ro3006_v2.py` (same cells, decoupled library, <= 514 diodes, one line, state-dependent beta, no neutral-static penalty because the statics carry beam 1 in both states). **This routine is a weaker optimiser than the one behind section 2 (at the old single beta it gives -9.8 / -8.1 dB where the section 2 routine gives -13.6 / -13.2), so swap and add are compared inside the same routine; the strong-routine swap is the last row.** SLL here = fixed 21 deg exclusion around every wanted beam, relative to the stronger beam; compromise cell = static cell whose ideal real m differs by more than 0.5 between the states; "must switch" = cells whose ideal m differs by more than 0.5.
| design (decoupled, 514 diodes) | switched cells | must switch / compromise static cells | state 1: own beam / other beam [dBi], SLL [dB] | state 2: beam 1 / beam 2 [dBi], SLL [dB] (vs weaker beam), balance | isolation S1 |
|---|---|---|---|---|---|
| beam swap, generic routine (g 1.5, mu 30) | 257 (19 %) | 77 % / 783 | 15.8 / 1.7, -10.9 | beam 2 14.5 (own), beam 1 2.0 (wrong), -8.6 | 14.1 |
| **add-a-beam, generic routine (eps_B 45)** | 257 (19 %) | 84 % / 889 | **18.4** / -2.4, **-13.7** | **16.6 / 10.0**, -12.9 (-6.3), **balance 6.6 dB** | 20.7 |
| beam swap, strong routine (v2 final; strict SLL) | 257 (19 %) | NA | 15.0 / -2.2 (isolation 17.2), -12.5 | beam 2 12.1 (own), -11.4 | 17.2 |
Findings: (1) add-a-beam does NOT reduce the number of cells that must switch (84 % against 77 % by the ideal-m criterion, 19 % realised in both because of the 514-diode cap): the beam-2 hologram spans the whole aperture. (2) State 1 is better (+2.6 dB gain, SLL -13.7 against -10.9 in the same routine, isolation 20.7 against 14.1) because the statics radiate beam 1 coherently instead of being neutral. (3) State 2 does not deliver equal power: beam 2 is 5.9-7.7 dB below beam 1 over the eight phase offsets (best 5.9 dB at eps_B 90): with 19 % of the cells switchable the added beam is limited; beam 2 at 10.0 dBi is 4.5 dB below the swap design's beam. (4) More compromise cells (889 against 783). (5) Not compared with the strong optimiser (not implemented for two simultaneous beams): the difference between the routines is about 3 dB in SLL, so the SLL advantage in state 1 is smaller than it looks.
