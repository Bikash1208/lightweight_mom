# Phase 20 R3 item 4 — locality test of the CP-2 cell

**Labels: MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET, VIA_OMITTED_APPROXIMATION.** CP-2 finalist A (`OPTIMISATION.md`; effective stub line 1.748 mm), MADP-000907-14020 (L_s = 0 assumed),
10 GHz, RT5880 1.575 mm, PEC ground, **R-handed incidence**, theta = 0 and 30 deg (phi = 0). Code: `src/ris/supercell.py` (tiling with per-site diode states, one shared Z0, sparse per-coding
loads, Floquet-order extraction; tests `tests/test_ris_supercell.py`: 1x1 supercell = `compute_s_parameters_oblique` to 1e-9, a lossless non-uniform supercell sums to power 1 over all
orders to 2e-3, a uniform supercell has no anomalous orders), driver `scripts/research/phase20_r3_locality.py`, statistics `phase20_r3_locality_summary.py`, data `results/phase20_r3/locality/`.
Mesh: 48 cells per 15 mm period (pitch 0.3125 mm) for BOTH the isolated cells and the supercells (so the staircase is the same), Floquet cutoff 97 on the supercell mesh; 4x1 supercell 45 s per angle,
4x4: 5 codings 33 / 27 min per angle, 11.4 GB peak. The mesh is coarser than the optimisation mesh (60): isolated r_RR phases differ by a few degrees; the comparison supercell vs isolated is
at equal discretisation.

**Prediction ("isolated cell")** = the project's Phase 13 array-factor convention: every site contributes its isolated-cell r_RR (R in, R out) with no element pattern and no coupling,
b_(m,n) = (1/N) sum_site r_site exp(+j 2 pi (m x + n y)/period_super). Measured order amplitudes are expressed in the fixed Cartesian circular frame (the per-order TE/TM frame rotation exp(+-j phi_o)
is removed). **Gate (brief): the anomalous-order level within 1 dB and 15 deg of the prediction.** Not tuned.

## (a) Periodic gradient 0-1-2-3 (4x1 supercell, the anomalous-reflection case) — handedness-preserving order
| case | intended order | predicted | measured | level error | phase error | gate |
|---|---|---|---|---|---|---|
| theta 0, states 0-1-2-3 | (-1, 0), +-30 deg | -1.0 dB | -1.0 dB | +0.08 dB | -0.7 deg | **pass** |
| theta 0, states 3-2-1-0 | (+1, 0) | -1.0 dB | -0.9 dB | +0.11 dB | +2.4 deg | **pass** |
| theta 30, states 0-1-2-3 | (-1, 0), reflected to normal | -1.0 dB | -0.9 dB | +0.02 dB | 0.0 deg | **pass** |
| theta 30, states 3-2-1-0 | (+1, 0) and (-3, 0) | -1.0 dB each | -11.7 / -9.2 dB | -10.7 / -8.2 dB | -11 / -12 deg | **FAIL** |

The failing case steers into near-grazing orders (|k_x|/k_0 = 0.9997 and 0.999, kz/k_0 ~ 0.02-0.04). The array-factor prediction aliases m = +1 and m = -3 (4 sites), so it assigns -1.0 dB to
BOTH (2 x 0.79 > 1: not a physical prediction); the full-wave solve puts -11.7 dB / -9.2 dB there and ~-3 dB in the mirror order (-2, 0) instead. This case is outside what the prediction
can describe (grazing aliased orders); it is reported as a failure of the gate, not removed. The pass cases are the useful ones: the PB phase ladder behaves locally for a 4-cell gradient to 0.1 dB and 3 deg.

## (b) Random 2-bit codings, 4x4 (5 draws, seed 20261020), strongest predicted non-specular order
| theta | draw | order | predicted | measured | level error | phase error | gate |
|---|---|---|---|---|---|---|---|
| 0 | 1 | (-1, 1) | -10.1 | -8.3 | +1.82 | +1.7 | fail (level) |
| 0 | 2 | (-1,-1) | -8.9 | -8.0 | +0.87 | -4.5 | pass |
| 0 | 3 | (-1,-1) | -10.4 | -9.0 | +1.43 | +5.6 | fail (level) |
| 0 | 4 | (1,-1) | -7.2 | -5.7 | +1.44 | -0.3 | fail (level) |
| 0 | 5 | (1,-1) | -10.7 | -14.4 | -3.71 | +8.1 | fail (level) |
| 30 | 1 | (-1, 1) | -9.7 | -11.0 | -1.27 | +4.5 | fail (level) |
| 30 | 2 | (-1,-1) | -8.8 | -10.2 | -1.39 | +2.7 | fail (level) |
| 30 | 3 | (-1,-2) | -5.7 | -15.4 | -9.75 | -8.0 | fail (level) |
| 30 | 4 | (-2,-1) | -8.1 | -5.6 | +2.52 | -9.7 | fail (level) |
| 30 | 5 | (-2,-1) | -7.9 | -8.1 | -0.16 | +3.4 | pass |

**The gate passes in 2 of 10 random cases.** The **phase** is within 15 deg in all 10 (|error| <= 9.7 deg: the geometric PB phase survives neighbours), the **level** is not: errors of
+-1 to 4 dB on the strongest order and a median |error| of 0.8-5.4 dB over all orders predicted above -25 dB (12 or fewer per draw; maximum single-order errors 9-21 dB, the worst on weak or near-grazing
orders). The isolated-cell array factor is therefore NOT a reliable predictor of the individual order levels of a random 2-bit pattern at this 15 mm (lambda/2) pitch; it is reliable for a
regular 4-cell gradient. (No attribution to a mechanism is made here: element pattern, near-field coupling between stubs of neighbouring cells and the near-grazing orders all contribute and were not separated.)

## (c) Where the opposite-handed term goes
The opposite-handed term (R in, L out) carries no PB phase: its phase does not follow the state, so in an array it should add in the specular order, not in the steered beam. Measured (gradient, 4x1):

| case | specular (0,0) | steered order | total over orders | share in the specular order (measured / predicted) |
|---|---|---|---|---|
| theta 0, 0-1-2-3 | -19.8 dB | -26.8 dB (order -1) | -18.2 dB | 69 % / 60 % |
| theta 0, 3-2-1-0 | -21.9 dB | -24.5 dB (order +1) | -18.7 dB | 48 % / 60 % |
| theta 30, 0-1-2-3 | -20.8 dB | -30.7 dB (order -1) | -19.2 dB | 68 % / 33 % |
| theta 30, 3-2-1-0 | -12.7 dB | -28.1 dB (order +1) | -12.3 dB | 93 % / 33 % |

The isolated cell has -13.6 dB opposite-handed (the target miss). In the steered beam the opposite-handed contamination is **-25 to -31 dB** (11-17 dB better than the isolated-cell figure) and
the larger part, 48-93 %, sits in the specular order — as expected, but not "all": the rest is spread over the other orders (4x4 random draws: only 13-60 % in the specular order). So at array level the -13.6 dB
miss is a **loss of the aperture** (the isolated cell reflects ~4 % of the incident power into the wrong handedness; about half or more of that ends in the specular order) rather than a contamination of the steered beam. Whether that is acceptable for the application is the user's call; the beam purity numbers above are model-only.

## Honest limits
Staircased octagon (mesh 48), coarser than the optimisation mesh; the uncertainty of the isolated-cell phases is ~10 deg (`OPTIMISATION.md`) and is NOT propagated here; L_s = 0 assumed; 4 gradient + 10
random cases are a small sample; the prediction convention is the project's own (no element pattern); the near-grazing orders (theta 30) are numerically delicate (kz/k_0 down to 0.02).
