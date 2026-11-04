# Step 2 — switchable gate cell, architecture A-2q (MODEL_ONLY, solver-internal)

Scripts: `scripts/research/sidequest_gate_cell.py` (cell, mode solver, strip line model), `sidequest_gate_screen.py`, `sidequest_gate_search.py`
(screening), `sidequest_gate_final.py` + `sidequest_gate_report.py` (characterisation) -> `gate_cell_final.json`; tests `tests/test_sidequest_gate_line.py`.
Imported additively from branch `phase28-coded-holographic-v2` (new files only, no solver change): `src/ris/diodes_measured.py`,
`results/phase28/diodes/ma4agp907_series_z.json` (MEASURED_VENDOR_SPARAM, vendor .s2p NOT committed), `src/hologram/coded_cell_eval.py`.

## Setup
RT5880 1.575 mm (eps_r 2.2, tan d 0.0009), PEC ground, 6.4 mm lattice, 10.54 GHz. One MA4AGP907 per cell, ON = 10 mA measured Z(f)
(R ~ 2.4 ohm), OFF = 0 V (C = 28.2 fF), conduction axis **x (0 deg) or y (90 deg) only: the rotated-load bug is avoided**. Reference
(hologram mean) sheet: Phase 26 design loading b = 7.82 mS, beta = 1.065 k0, TM0 characteristic X = 138 ohm.

**Cell ("dog-bone")**: lattice 6.4 mm; two end plates 1.2 x 6.0 mm (x by y) at the cell ends (cell length 5.6 mm, neighbour gap 0.8 mm), a centre strip 0.8 mm
wide, the diode footprint 0.4 x 0.8 mm bridging a 0.4 mm gap in the strip. Orientation y: the same cell turned by 90 deg (not rotated loads: separate
axis-aligned build). Every edge sits on a mesh line (0.2 mm at mesh 32, 0.1 mm at mesh 64; `make_dogbone` refuses other values).
The gate state map: **OFF = pass, ON = block** (a normally-pass cell). Its y_xx / y_yy [mS] at 10.54 GHz, mesh 64:

| state | y_xx | y_yy | TM0 beta/k0 (x propagation) | 3-row strip, ref | strip |
|---|---|---|---|---|---|
| OFF (28 fF) | 0.00 + 1.93j | 5.72j | 1.0243 (alpha 0.011 Np/m) | | S21 = -0.93 dB, S11 = -7.2 dB |
| ON (2.4 ohm) | 99.3 - 71.4j | 5.72j | **no bound TM0 mode** (NO_BOUND_MODE) | | blocked by cut-off, not quantifiable in this model |

## Stopband (ON state, sheet model, mesh 32; mesh 64 spot checks agree to ~0.1 GHz)
Re Y_xx peaks at 153 mS near 10.4 GHz (cell resonance). **No bound TM0 mode from 10.1 to 10.54 GHz** (mesh 64 also at 10.2 GHz, 10.54 GHz), i.e. a ~0.45 GHz
(4 %) cut-off band around the design frequency. Above it the mode is weakly bound (beta ~ 1.0001 k0) and the strip gives S21(3 rows) = -16.3 dB at
10.6 GHz, -14.9 / -13.7 / -12.8 / -12.1 dB at 10.7 / 10.8 / 10.9 / 11.0 GHz (mesh 64 at 10.9: -13.3 dB). The OFF (pass) state is flat: IL 0.6-1.2 dB, RL 8.7-6.3 dB over 10.0-11.0 GHz.

## Result vs the targets
| Quantity | Target | Result |
|---|---|---|
| Isolation of a 3-row strip | 20 dB | **not demonstrated.** 16.3 dB at the nearest bound-mode frequency (10.6 GHz), cut-off (no TM0 mode, unquantified, likely > 16 dB) at 10.1-10.54 GHz |
| Pass-state match | matched to the hologram mean | **not matched**: b_xx 1.9 mS vs 7.8 mS, beta 1.024 vs 1.065 k0 -> 3 rows: IL 0.93 dB, RL 7.2 dB (2 rows: 0.1 dB / 18 dB; 4 rows: 0.3 / 12.4 dB; Fabry-Perot ripple) |
| +-20 % diode tolerance (mesh 64) | | OFF C +-20 %: b_xx 1.75 / 2.12 mS, IL 0.98 / 0.89 dB (insensitive); ON R x0.8 / x1.2: still NO_BOUND_MODE (mesh 32: -17.7 dB / no mode). The ON state is resonance-limited, not diode-limited |
| Direction dependence (x-axis cell, ON, 3 rows) | | blocked (no mode) up to 45 deg off the diode axis, -6.0 dB at 60 deg, -0.8 dB at 75 deg, **-0.15 dB at 90 deg (no block)**. The ring therefore needs x-cells for |psi - {0, 180}| <= 45 deg and y-cells for the rest: two diode orientations, same cell, two bias groups per quadrant |
| Number of diodes | 1 per cell | 1 per cell; ring r_front = 13 mm: ~68 cells (3 rows) |

## Required contrast (what the 20 dB needs, any cell)
Single-mode line model with the reference sheet, 3 rows, pure susceptance b of the gate: S21 = -6.9 / -12 / -17 / -22 / -30 dB for b = -20 / -50 / -100 / -200 / -500 mS
(beta/k0 -> 1.0001; mode reactance X from 138 ohm to ~5 ohm). **20 dB per 3-row strip needs |b_ON| >= ~150 mS (inductive, i.e. E_x almost shorted), loss-free,
with the gate wide enough to be a clean section; a capacitive block never exceeds ~4 dB** (beta saturates at 1.48 k0; there is no TM0 cut-off from capacitive loading
in this single-layer stack: the surface wave is bound for every sheet susceptance, only a resonant, inductive ON state works). The cell reaches |b| ~ 70-90 mS with
g ~ 80-100 mS at 10.5 GHz (resonance), enough for 16-17 dB, not 20.

## Findings / limits
1. **Single-layer, via-less gates cannot be strongly blocking**: the only mechanism is the cell's own resonance, so isolation is narrowband (a 4 % cut-off band, 12-16 dB
   above it), mesh-sensitive (mesh 32 vs 64: Y_ON 81.8-78.0j vs 99.3-71.4j, the resonance moves ~0.1 GHz) and tolerance-sensitive in geometry (not in the diode).
   A more reliable block would short the patch to ground through the diode (via) or add rows: the solver does not model vias (Phase 20 note).
2. The pass state is mismatched (RL 7 dB, 3 rows). Step 1 assumed an ideal pass state, so the open-quadrant S11 and gain of step 1 are optimistic by this reflection;
   the pass reflection can be phase-tuned together with the ring radius but not removed. A static parallel capacitance (interdigital) to lift b_OFF to ~7 mS was not found
   within the 0.4 mm gap rules (b_OFF <= 3 mS in every aligned variant, 160+ screened).
3. **A normally-block complementary cell (needed so that "all diodes OFF" = quadrant 1 radiates) was not found**: with a single normally-pass cell, state 1 needs the
   quadrant-2 gate ON and state 2 the quadrant-1 gate ON; the all-OFF default radiates BOTH beams (two-beam state).
4. Method limits: transmission through the strip is the single-mode TM0 junction/propagation model (`strip_transmission`); power a junction mismatch would scatter into
   radiation is not removed, so |S21| is an UPPER bound on transmission (isolation is conservative where a bound mode exists); where NO_BOUND_MODE the model gives no number.
   HFSS (or a full-wave strip run) is needed to confirm. Screening artefact logged: variants with edges on half-mesh lines (cell length 6.2 mm at mesh 32) gave spurious
   resonances (apparent 20.6 dB, b_OFF 7.8 mS); they disappear when edges are aligned, and `make_dogbone` now refuses unaligned dimensions.
5. 4 quadrants later: the same cell (x and y orientations), 4 bias groups; at 4q the ring is the same 68 cells, S11/gain consequences as in `ARCH_STUDY.md`.
