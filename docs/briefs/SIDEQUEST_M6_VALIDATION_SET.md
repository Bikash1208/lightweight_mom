# SIDEQUEST_M6_VALIDATION_SET - HFSS data the new board model must reproduce (2026-11-06)
Status: data record; HFSS values as reported by the user (not re-derived here). Model-side predictions for items 3-5 are **NOT_YET_COMPUTED** (no new cell exists; the retired model's numbers are not to be reused as a pass).

| # | Board / state | HFSS observation | Acceptance for the new model |
|---|---|---|---|
| 1 | board_K2_default, state 1 and 2 | peak gain, lobe positions, HPBW (earlier reports) | +-2 dB peak, +-3 deg lobe, +-3 deg HPBW |
| 2 | add-a-beam (b) ON, far field | broad dome, D 15.1 dBi; beam 1 split: 13.6 dBi at (20, -18) X-pol and 13.2 dBi at (32, +14); radiation efficiency 0.72 | same, +-2 dB, +-3 deg |
| 3 | add-a-beam (b) ON, **near field on the top surface** | \|E\| falls about 30 dB from the probe (about 78 dB) to the edge (about 45 dB), i.e. 15-20 dB beyond cylindrical spreading, mostly within rho = 30-40 mm; faint phase rings at 20-25 mm spacing; -x side weaker | \|E\|(rho) along x and y within +-3 dB over rho = 10-75 mm; ring spacing within +-3 mm; -x/+x asymmetry reproduced in sign |
| 4 | **v2 (beam swap) ON, HFSS** | P_rad = 0.155 W per 1 W incident (i.e. most of the power is not accounted by radiation into the main beams), peak gain -3.4 dBi, D 4.7 dBi, 46 % of the radiated power at theta 60-90 deg | P_rad, D and the theta 60-90 deg fraction within +-2 dB / +-5 % |
| 5 | v2 OFF / ON gain offsets | state 1 -1.78 dB, state 2 -2.7 dB between model and HFSS (old model) | the new model needs no fitted offsets; if it does, M6 fails |

Reading: items 3 and 4 say that in the ON state the surface wave leaks far more strongly than the retired model predicted (alpha larger by roughly 15-20 dB of extra decay within 30-40 mm; most of the radiated power goes to low elevations). The new cell and board model must carry alpha per cell from the Bloch / driven-strip extraction, not from a perturbation formula, and must include the radiation from the transitions (rim, feed region) that produced the 46 % at 60-90 deg.
Model-side values for 3 and 4: NOT_YET_COMPUTED (item 4 of the task).
