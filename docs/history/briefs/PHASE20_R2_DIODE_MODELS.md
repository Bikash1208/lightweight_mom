# Phase 20 / R2 — PIN-diode models for the loaded MoM

**Status:** active brief, written 2026-10-17. Starts after R1's close-out (item 5). Same branch `phase20-ris`. Read
`docs/briefs/PHASE20_RIS_ROADMAP.md` and the R1 results (`results/phase20_r1_demo/README.md`) first. All repo rules apply
(one commit per item, `make test-fast` + snapshot gates green, physics tests, honest labels, additive only).

## Why
R1 gave the MoM lumped loads and oblique incidence. A RIS cell needs those loads to be **real PIN diodes**, ON and OFF,
with their frequency dependence, at the size of a real diode, and in a form the user can reproduce in HFSS. R1 also
found three things this phase must handle:
1. A load with `transverse_conducts=True` leaves a floating strip; the physical lumped element is
   `transverse_conducts=False` (or `open_circuit=True` for an open).
2. The R1 demo load was a full-width 8 mm strip, so TE drove circulating currents through it (TE absorption up to 0.68
   at 60°). A real diode is ≈ 0.3–0.6 mm wide; that behaviour should largely disappear — verify it.
3. HFSS lumped RLC boundaries put R, L and C **in parallel**; a PIN diode is a **series** circuit (ON: R_s + L_s;
   OFF: C_T + L_s, optional parallel R_p). The HFSS representation must be worked out, not assumed.

## Work items (one commit each; stop point ■)

**1. `src/ris/diodes.py`.** Diode = named equivalent circuit → Z_L(f). Models: PIN ON (R_s + L_s), PIN OFF
(C_T ‖ R_p, in series with L_s), plus generic series RLC. Values for Skyworks SMP1345-079LF from the manufacturer's
datasheet, each value cited (document, revision, page, bias condition); package parasitics separately. If a value is
not in the datasheet, report `NOT_AVAILABLE_FROM_CURRENT_DATA` and ask the user — do not invent it. All values labelled
`PROVISIONAL_DATASHEET` (no measured de-embedding yet). Tests: passivity (Re Z ≥ 0 over 1–20 GHz), limits (L_s → 0,
C_T → ∞ give the R1 short; C_T → 0 gives open), series resonance frequency of OFF state equals 1/(2π√(L_s C_T)).

**2. Diode → R1 load.** Feed Z_L(f) into `LumpedLoad` per frequency, with `transverse_conducts=False` as the diode default
and the load region sized to the diode/pad footprint (parameters, not constants). Tests: frequency sweep 8–12 GHz runs the
load per frequency; a frequency-independent model reproduces the R1 fixed-R/L/C load bit-identically.

**3. HFSS representation.** Work out and document how the user should model ON and OFF in HFSS given parallel-only RLC
boundaries — e.g. two RLC sheets in series across the gap, or a lumped port terminated in a circuit, or another
standard practice; cite the source. Then show in the MoM that the chosen HFSS representation (as geometry + loads)
gives the same S11 as the ideal series model within mesh tolerance at 10 GHz, θ = 0/30°, TE/TM.

**4. Demo: one diode in a gapped patch.** Reuse the R1 demo cell (10 mm, RT5880 1.575 mm, PEC ground) with a narrow
diode-sized load: ON vs OFF phase and |S11| over 8–12 GHz at θ = 0°, 30°, 60°, TE and TM. Report: ON/OFF phase
difference vs f, loss in each state, and how the TE behaviour changed vs the R1 8 mm strip. Sensitivity: ±20 % on each
diode parameter one at a time → phase spread at 10 GHz. Mesh: use the R1 convergence practice (extrapolation with the
observed order reported; widen the stated uncertainty where the order < 2). **■ Stop and report.** Write an HFSS check
package for the ON and OFF cell (`hfss_exports/phase20_r2_diode_cell/`, predictions first) only after the user agrees.

## Do not
- Design the 2-bit cell (that is R3) or pick the final topology.
- Model the bias network or vias beyond the roadmap's labelled approximations.
- Fit diode values to make results look good; datasheet values only, labelled.
