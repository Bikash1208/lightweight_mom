# Crop testbed: do the switched cells switch? (HFSS 2025.1, direct solver, 32 cores)

**STATUS:** three small boards solved and converged. In HFSS the switched cells' coherent field with the diodes OFF is about 0.8 of the diodes-ON field, not the 0.262 the model uses.

## Setup
Layout = sites with rho <= 64 mm of `layout_FINAL2.json` (84 sites, 21 switched; cell area R 70 mm, 30 mm idealised absorber ring, board R 100 mm). Cases: all sites STATIC (same 21 sites as static cells), switched cells with diodes OFF (-10 V, 21.6 - j297.5 Ohm/sq), switched cells with diodes ON (2.8 + j14.6). Auto mesh, no mesh operations, Setup1 5.8 GHz, 20 passes max / 6 min / dS 0.02, direct solver. validate_simple 1 on all designs.

## GATES
| Case | Passes | Final max dS | Direct solver |
|---|---|---|---|
| all static | 6 | 0.0051 | yes |
| switched OFF | 6 | 0.0090 | yes |
| switched ON | 6 | 0.0091 | yes |

## DATA (HFSS far field, per 1 W incident; subtraction = switched minus static at the same sites, independent meshes, ratios good to about +-1 dB)
- Radiated power: static 0.170 W, OFF 0.196 W, ON 0.202 W. m = 0 Etheta carries 0.162 / 0.178 / 0.180 W.
- RHCP power in the phi 0 beam cone (theta 24, 14 deg): static 0.0059, OFF 0.0153, ON 0.0192 W. phi 180 cone: 0.0053 / 0.0062 / 0.0063 W.
- Switched-cell contribution (state minus static) in the phi 0 cone: ON/OFF amplitude ratio 1.24 (+1.9 dB), phase 5 deg, correlation 0.98. So OFF/ON = 0.81 (-1.9 dB). The model assumes 0.262 (-11.6 dB).
- Model check (MODEL_ONLY, FINAL2 layout, `suppression_dB` = switched beam, two-beam minus one-beam, absolute): r = 0.262 -> 11.5 dB, r = 0.5 -> 5.6, r = 0.74 -> 2.3, r = 0.81 -> 1.6; targets r = 0.20 -> 14.1, 0.15 -> 16.9, 0.10 -> 21.2 dB. HFSS full board, RHCP power in the phi 0 cone, two-beam -> one-beam: 0.0367 -> 0.0331 W = 0.45 dB, so the full board switches even less than r = 0.81 predicts (something state-independent also feeds that beam). CORRECTION: an earlier version of this file compared directivities (1.0 dB vs 1.1 dB); that comparison is withdrawn.
- All-static small board: no beams (RHCP D 4.2 dBi at both beam directions), 95 % of the radiated power is m = 0 Etheta (0.162 W of 0.170 W).

## FLAGS
- phi 0 beam is the switched beam (B); it appears only when switched cells are present. A is the static beam at phi 180.
- m = 0 Etheta(theta) shape of the crop-static board and the full board correlate only 0.31, so the ring is not purely a centre source; not resolved.
- S11 / accepted power not read for the crop cases (GetSetups call fails in this PyAEDT session).
- The crop is a diagnostic board, not a design.

## FILES
static_off.ffd, sw_OFF_ELS10.ffd, sw_ON_ELS.ffd, layout_crop_*.json, readback_*.json, solve_*.json.

## Addendum: accepted power and idealised bounds (HFSS, same settings, 6 passes, dS <= 0.0091)
- S11 / accepted power at 5.8 GHz: static -14.6 dB / 96.5 %, switched OFF -15.7 dB / 97.3 %, switched ON -16.2 dB / 97.6 %, open bound -15.6 dB / 97.2 %, short bound -15.9 dB / 97.4 %. Accepted minus radiated is about 0.77-0.80 W in all five cases (absorber ring + dielectric); case-to-case differences are within about 0.02 W, so the diode-sheet loss is too small to resolve.
- IDEALISED BOUNDS (not a diode sweep): sheets open (Rs 1e5 Ohm/sq) and shorted (Rs 0.001). Switched-cell contribution (state minus all-static) in the phi 0 cone relative to ON: OFF 0.77, open 0.74, short 0.88. Any capacitive (reverse-bias) load lies between open and short, so it cannot go below about 0.74 of ON; the 0.13/0.18 pF variants cannot reach the needed 0.1-0.26. Inductive loads are not bounded by these two cases.
- A single-pole fit E(Z) = E_open + K/(Z + Z_in) to open/short/ON predicts OFF poorly (residual 0.44 of the small OFF-open difference), so it is approximate; it suggests the total amplitude varies only between about 0.8 and 1.2 of ON over passive loads (R >= 2.8 Ohm/sq).
- Idealised-bound solves peak at about 17 GB.

## Addendum 2: where the m = 0 ring comes from (HFSS, 6 passes each, dS <= 0.0032)
Per 1 W incident; "ring" = m = 0 Etheta radiated power.
| Case | S11 | Accepted | Radiated | ring | ring / accepted |
|---|---|---|---|---|---|
| all-static crop, pin top 8.925 mm | -14.6 dB | 96.5 % | 0.170 W | 0.162 W | 0.168 |
| same, pin top 4.5 mm | -1.1 dB | 22.3 % | 0.041 W | 0.039 W | 0.177 |
| same, pin top 2.0 mm | -0.2 dB | 5.2 % | 0.010 W | 0.009 W | 0.182 |
| NO CELLS (pin + slab + ground + bias layer + absorber only) | -16.2 dB | 97.6 % | 0.237 W | 0.237 W | 0.243 |
- The ring is present with no cell at all: 24 % of accepted power (all of it m = 0, as an axisymmetric structure must give, which also validates the m = 0 extraction). Static crop 17 %, full board about 20 %. The ring starts at the launch; the no-cells guide is a different structure, so these numbers are not a decomposition (the m = 0 profile of the no-cells crop correlates only 0.15 with the static crop).
- Pin shortening changes the accepted power by a factor of 20 but not the ring fraction. This does NOT show whether the pin itself radiates (a short probe's own radiation and the wave it launches scale with the same current); it shows only that a shorter pin is not a fix because the match collapses.
- WITHDRAWN: the earlier statement "it is not direct pin radiation".
- The static-crop m = 0 profile is a single smooth hump (peak theta 22 deg, no nulls, smooth phase), so a rim source is unlikely on the crop.
- NOT CHECKED: static-cell leakage in watts against the model's converted fraction (HFSS static crop has 0.008 W outside the ring).

## Addendum 3: noise floor, baselines, matched-dummy test (HFSS, 6-8 passes, direct solver)
- Noise floor: the all-discs-at-the-21-sites reference solved at dS 0.02 and again at dS 0.005 (8 passes) differs by 0.000018 W in the phi 0 cone and 0.000005 W in the phi 180 cone (radiated power 0.191 vs 0.186 W).
- Baseline dependence of r (switched-cell contribution, OFF relative to ON, phi 0 cone, projection): zero = all-static crop 0.77-0.79; zero = unslotted 5.435 mm discs at the 21 sites 0.39-0.45. The two zeros disagree, so r is only stated as 0.4-0.8 from the crop; the full board (reference-free) shows under 1 dB.
- MATCHED-DUMMY TEST (84 sites all carry switched-cell copper; 63 host sites on their A axes with open sheets; 21 B sites with real sheets OFF / ON; zero = same board with the 21 B sites turned to their A axis, all sheets open): r(OFF/ON) = 0.81-0.86 against the matched zero (raw 0.86, ring removed 0.81); raw B-cone power OFF/ON = -0.75 dB. Cone watts of differences: ON-zero B 0.00345 W, OFF-zero B 0.00257 W, ON-OFF B 0.00020 W, A 0.00024 W. An r of 0.3 would need OFF-zero of about 0.0003 W. The neighbour-mismatch hypothesis is REFUTED: with all cells the same copper, the B cell still radiates about 0.8-0.9 of its ON level when OFF. r should be quoted as about 0.8; the ON-OFF effect sits near this board's noise level.
- Full board with the ring removed: phi 0 cone 0.0244 W ON, 0.0238 W OFF (phi 180: 0.0228 / 0.0251 W): the B beam does not respond to the diode state on the full board either.
- Dummy lattice outside-ring radiation: 0.058-0.059 W (OFF/ON), 0.078 W (zero) vs 0.008 W for the all-static crop; below-horizon 9 %, LHCP share about 50 %. Flagged only.
- S11 / accepted: dummy OFF -18.9 dB / 98.7 %, dummy ON -19.9 dB / 99.0 %, zero board -17.5 dB / 98.2 %. Peak memory of the dummy boards about 95 GB (the static-host crop needs about 17 GB).
- The eigenmode-Delta proxy (r = |Delta_OFF|/|Delta_ON| = 0.26) has no support from any HFSS board test.

## Addendum 4: orientation test, step 1 (HFSS, static-host crop, 21 B cells turned 90 deg, direct solver; MODEL_ONLY plan from the advisor)
- Runs: `OPEN_IDEAL` (O90) and `ON_ELS` (ON90) on `layout_crop_sw90.json`; both 6 passes, Normal Completion (target dS 0.02 met, value not read). ON90 S11 -15.0 dB, accepted 0.968. ON90 took 1:43 on 64 free cores (O90 8:38 on shared cores), same solver. Fields: `or90/O90_1.ffd`, `or90/ON90_1.ffd`; scoring `or90/or90.py`, output `or90/or90_out.txt`.
- Oriented part = half the difference of the 0 and 90 deg far fields (open: OPEN_IDEAL 0 deg vs O90; ON: sw_on 0 deg vs ON90). The 0 deg and 90 deg states come from different projects/meshes (noise floor 0.0001-0.0002 W for this pairing).
- RHCP oriented power in a 10 deg cone around design B (theta 25, phi 0): open 0.00102 W, ON 0.00169 W, so quiet/radiating r = open/ON = 0.78 (15 deg cone 0.81). Around the argmax (20, 22) r = 0.95. The centre matters: the number is 0.78-0.95 depending on it. Around design A (phi 180) the oriented powers are 0.00002-0.00009 W, i.e. below the noise floor.
- Raw RHCP in the design-B 10 deg cone: open 0 deg 0.00833 W, open 90 deg 0.00125 W, ON 0 deg 0.01085 W, ON 90 deg 0.00077 W.
- Reading: the calibration gives about 0.8 at the designed B direction, consistent with the earlier switched-minus-static r. Step 2 (bridging bounds) decides whether the copper can be fixed by loading.
