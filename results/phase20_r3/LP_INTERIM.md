# LP interim: 1-bit linear-polarisation codes of the CP-2 cell (Phase 20 R3, decision 2026-10-20)

**Labels: MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET, VIA_OMITTED_APPROXIMATION.** Cell: CP-2 finalist A (`OPTIMISATION.md`), MADP-000907-14020, 10 GHz, MoM mesh 90, data in
`lp_interim.json` (`scripts/research/phase20_r3_cp_package.py`). The 2-bit LP cell targets were **not met** (best: LP-3, 3 diodes, nearest miss; see `OPTIMISATION.md`). All the LP
searches of item 3 used the MADP-000907 diode (`kit` in every `lp*_madp*.json`): the SMP1345 OFF loading was NOT the cause of the stuck all-OFF state; SMP1345 had been used only in the
item-2 starting-cell scans. Interim: LP is offered as **1-bit codes of the CP-2 cell**. 2-bit LP is a future phase needing a multilayer / slot-coupled cell.

The CP-2 states are the same half-wave-plate-like element turned by alpha = 0 / 45 / 90 / 135 deg, so for linear incidence the reflection is a Jones matrix of the form
r [[cos 2a, sin 2a], [sin 2a, -cos 2a]]: states 1 and 3 (alpha = 0, 90) reverse the **co-polar** sign (LP -> co-LP code), states 2 and 4 (alpha = 45, 135) reverse the
**cross-polar** sign (LP -> cross-LP code, i.e. the reflected linear polarisation is turned by 90 deg).

| code (x-polarised / TM incidence) | theta | phase difference | amplitudes (dB) | unwanted term vs wanted (dB) |
|---|---|---|---|---|
| **LP -> co-LP**, states 1 vs 3 | 0 | 163 deg | -1.70 / -0.11 | cross-pol -17.6 / -23.1 |
| LP -> co-LP, states 1 vs 3 (TM) | 30 | 178 deg | -1.91 / -0.12 | cross-pol -16.6 / -22.7 |
| **LP -> cross-LP**, states 2 vs 4 | 0 | 178 deg | -1.04 / -0.95 | co-pol -16.9 / -23.1 |
| LP -> cross-LP, states 2 vs 4 (TM in, TE out) | 30 | 172 deg | -0.99 / -0.98 | co-pol -21.1 / -16.9 |

y-polarised / TE incidence is similar (co code 157 deg at theta 0 and 173 deg at 30; cross code 178 / 176 deg) except that the leftover co-polar term of the cross-LP code for y-incidence
is only -12.7 dB below the wanted term at theta 0 (-12.4 dB at 30 deg).

Reading: the 1-bit LP codes of this cell are close to 180 deg and within ~2 dB amplitude of each other, with cross-pol 17-23 dB below the wanted term in most cases and a worst case of
12.4-12.7 dB (y/TE incidence, cross-LP code). The phase difference of the co-LP code at normal incidence is 163 deg / 157 deg, i.e. 17-23 deg from the ideal 180 deg (model uncertainty ~10 deg).
These are the SAME eight HFSS solves as the CP check (`hfss_exports/phase20_r3_cp_cell/`): no separate LP package.
