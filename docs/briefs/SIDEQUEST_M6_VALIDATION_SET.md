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
Model-side values for 3 and 4: see the section below (computed 2026-11-06, M6 NOT_MET).

## Model-side result (2026-11-06): SWHB Born board model against the old HFSS boards. **M6 = NOT_MET (partial)**
Model: `scripts/research/sidequest_swhb.py` (cell sheet currents J = dY E, exact grounded-slab element factor, incoherent + coherent radiation, SW depletion, loss; launch efficiency not modelled). Old-cell input: MoM mesh-96 library, raw (cal 1.0) and with the old HFSS-derived contrast factor (cal 0.75); both reported. Script `sidequest_swhb_m6.py`.
| item | HFSS | model raw / cal 0.75 | verdict |
|---|---|---|---|
| 3 add-a-beam ON, probe-to-edge \|E\| drop | about 33 dB (78 -> 45), 15-20 dB beyond cylindrical spreading | 35.7 dB / 27.9 dB (13.5 dB spreading + edge power 0.006 / 0.036) | within +2.7 / -5.1 dB: PARTIAL (band), mechanism reproduced: diffuse leakage of the quantised 1-bit pattern is strong |
| 2 add-a-beam ON far field | D 15.1 dBi, split beams (20,-18) 13.6 dBi and (32,+14) 13.2 dBi, eta 0.72 | single beam (25,0), D 17.8-19.8 dBi, G 17.6-19.2 dBi | FAIL (gain +4 to +6 dB, no split) |
| 2b same, cell beta shifted by +0.15 (common susceptance +2 mS) | as above | peak (31,15), D 14.5 dBi | beam break-up reproduced qualitatively: the HFSS board behaves like a cell beta error of about +0.15 |
| 4 v2 ON | P_rad 0.155, peak -3.4 dBi, D 4.7 dBi, 46 % at theta 60-90 | P_rad 0.48 / 0.28, D 21.3 / 22.0 dBi, 18 % / 15 % at 60-90 | FAIL |
| 4b v2 ON, beta shifted by +0.2 | same | P_rad 0.14-0.15 (matches), D 15-16 dBi, 10-15 % at 60-90, edge power 0.55-0.56 | P_rad reproduced, pattern NOT: the model has no rim diffraction (the 0.55 of the power reaching the edge is not radiated) |
Reading: (i) the leakage level and the near-field decay are reproduced within a few dB; (ii) the far field of the old boards is reproduced only if the old cell beta is wrong by about +0.15, which is consistent with the known old-cell errors (OFF branch not mesh-converged, contrast x0.69-0.81), i.e. the validation set mainly tests the cell data, not the board model; (iii) the rim diffraction is missing, so the diffuse 60-90 degree content of v2 ON is not reproduced. **Consequently every SWHB board prediction is MODEL_ONLY with wide bands (gain +-3 dB, SLL +-3 dB, pointing +-3 deg, edge/rim content unknown) and must not be quoted as validated.** Re-test after the new cell passes HFSS (item 3 of the cell task).
