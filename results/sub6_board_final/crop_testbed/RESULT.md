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
- Model check (MODEL_ONLY, FINAL2 layout): with r = 0.262 the switched beam falls 9.4 dB from two-beam to one-beam; with r = 0.81 it falls 1.0 dB. HFSS full board: 1.1 dB (16.5 -> 15.4 dBi).
- All-static small board: no beams (RHCP D 4.2 dBi at both beam directions), 95 % of the radiated power is m = 0 Etheta (0.162 W of 0.170 W).

## FLAGS
- phi 0 beam is the switched beam (B); it appears only when switched cells are present. A is the static beam at phi 180.
- m = 0 Etheta(theta) shape of the crop-static board and the full board correlate only 0.31, so the ring is not purely a centre source; not resolved.
- S11 / accepted power not read for the crop cases (GetSetups call fails in this PyAEDT session).
- The crop is a diagnostic board, not a design.

## FILES
static_off.ffd, sw_OFF_ELS10.ffd, sw_ON_ELS.ffd, layout_crop_*.json, readback_*.json, solve_*.json.
