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
