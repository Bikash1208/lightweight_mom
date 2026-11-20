# v3.1 full-board sidelobe diagnostic — stored HFSS far fields

## STATUS

COMPLETE. Post-processing only at 5.68 GHz; no AEDT process, geometry operation, or new solve. The saved BOTH, A-only-IDEAL (234), B-ALL-OFF-IDEAL, and ALL-PLAIN far fields were used.

## METHOD

For each hand, the peak and its absolute gain came from the converged HFSS metrics. The stored complex FFD gives the relative gain at every 2° grid point. The FFD-derived gains match the direct HFSS gain readbacks at (32°,0°) and (36°,180°) to the printed precision. The existing SLL definition is the highest same-hand gain above the horizon outside the first-null main-lobe mask (2× mean HPBW fallback). The second lobe is the next highest candidate outside that mask and at least 15° from the first, using great-circle separation. `lobes.csv` includes the peak and five separated sidelobe candidates per state and hand.

## PEAKS AND LOBES

Each lobe cell gives **level dBi / dB below that hand's peak at (theta°,phi°)**.

| State | Hand | Peak dBi at (theta,phi) | Highest sidelobe | Second sidelobe (>=15° away) |
|---|---|---|---|---|
| BOTH | RHCP | +9.03 at (32,0) | +0.61 / -8.43 at (10,164) | +0.23 / -8.80 at (40,330) |
| BOTH | LHCP | +9.00 at (36,180) | -0.53 / -9.52 at (18,130) | -0.63 / -9.62 at (36,146) |
| A-only-IDEAL, 234 | RHCP | +8.32 at (32,0) | +0.11 / -8.22 at (44,188) | -0.27 / -8.59 at (24,324) |
| A-only-IDEAL, 234 | LHCP | -1.59 at (44,176) | -2.11 / -0.52 at (44,144) | -3.14 / -1.55 at (48,84) |
| B-ALL-OFF-IDEAL | RHCP | +8.60 at (32,0) | -0.55 / -9.15 at (44,188) | -0.99 / -9.59 at (44,332) |
| B-ALL-OFF-IDEAL | LHCP | -1.91 at (44,144) | -2.82 / -0.91 at (48,358) | -3.19 / -1.28 at (34,176) |
| ALL-PLAIN | RHCP | -3.87 at (48,278) | -3.94 / -0.07 at (48,186) | -3.97 / -0.09 at (48,8) |
| ALL-PLAIN | LHCP | -3.78 at (48,84) | -3.80 / -0.03 at (48,264) | -3.87 / -0.09 at (48,174) |

## GAIN AT FOUR FIXED DIRECTIONS

All values are gain in dBi, directly reconstructed from the stored HFSS FFD and its HFSS peak calibration.

| State | Hand | (44,144) | (44,176) | (32,0) | (36,180) |
|---|---|---:|---:|---:|---:|
| BOTH | RHCP | -2.40 | -1.71 | +9.03 | -2.33 |
| BOTH | LHCP | -2.05 | +5.99 | -19.25 | +9.00 |
| A-only-IDEAL, 234 | RHCP | -4.14 | -0.88 | +8.32 | -1.74 |
| A-only-IDEAL, 234 | LHCP | -2.11 | -1.59 | -26.94 | -2.80 |
| B-ALL-OFF-IDEAL | RHCP | -5.16 | -2.35 | +8.60 | -3.43 |
| B-ALL-OFF-IDEAL | LHCP | -1.91 | -2.78 | -14.50 | -3.48 |
| ALL-PLAIN | RHCP | -6.94 | -5.68 | -18.10 | -11.30 |
| ALL-PLAIN | LHCP | -5.32 | -5.06 | -18.11 | -11.60 |

## BOTH-STATE REGION CHECK

Excluding all points within 20° of **either** (44,144) or (44,176) leaves the highest RHCP SLL unchanged at **+0.61 dBi, -8.43 dB at (10,164)** and the highest LHCP SLL unchanged at **-0.53 dBi, -9.52 dB at (18,130)**. Their directions are **11.16° apart**, so **yes for RHCP and yes for LHCP**, the two highest sidelobes occupy the same ±20° direction region. Separately, **no for RHCP and no for LHCP**, neither highest sidelobe is within 20° of (44,144) or (44,176): nearest distances are 34.27° and 26.82°, respectively.

## FLAGS

For the suppressed and all-plain states, the global residual peaks are not intended beams; “sidelobe” means the next spatially separated residual maximum under the same algorithm. The reported second candidate is at least 15° from the first; a broad lobe can have spatially separated sampled structure. The phi=360° FFD column duplicates phi=0° and was excluded from candidate selection.

## FILES

`lobes.csv` lists the peak and five separated sidelobe candidates per state/hand, plus BOTH's residual-region-excluded sidelobe. `sll_diagnostic.json` contains unrounded gains, lobe directions, target samples and great-circle region checks. The underlying FFD files and converged solver readbacks are in the preceding v3.1 result folders.
