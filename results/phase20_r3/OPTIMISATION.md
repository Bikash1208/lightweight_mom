# Phase 20 R3 item 3 — optimisation of the 2-bit cells

**Labels: MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET, VIA_OMITTED_APPROXIMATION, ASSUMED_SENSITIVITY_NOT_DEVICE_DATA (MADP L_s).** 10 GHz, RT5880 1.575 mm, PEC ground,
zero-thickness metal, bias via and bias network omitted. Nothing here is an HFSS result. All data: `results/phase20_r3/opt/*.json`; drivers: `scripts/research/phase20_r3_opt*.py`,
`phase20_r3_finalist.py`, `phase20_r3_lp_report.py`, `phase20_r3_band.py` (optimiser kept OUTSIDE the solver core; Nelder–Mead, no surrogate, no fitted diode values).

## Verdict per family (targets from `docs/briefs/PHASE20_R3_TWO_BIT_CELLS.md`; misses are reported, nothing relaxed)

| family | best cell | verdict |
|---|---|---|
| **CP-2** (MADP-000907-14020, period 15 mm) | octagon circumradius 4.60 mm, stub line 2.48 mm, stub width 1.37 mm (other dimensions as `Cp2Cell` defaults: diode span 0.75, pad 0.25, overlap 0.1) | **meets phase, amplitude and tolerance targets at θ = 0 and 30° within the extrapolation uncertainty; misses opposite-handed isolation by ~1.4 dB; degrades at 60°** |
| **LP-3** (3 diodes, MADP) | patch 6.53 × 3.37 mm, period 14.97 mm; codes (000),(010),(011),(110) | nearest LP miss: phase deviation 21.7° (target 15°), amplitude −1.74 dB (target −1.5), spread 1.7 dB (target 1), tolerance step errors 30–49° (target 20°) |
| **LP-2** (2 diodes, MADP) | best of 17 starts | **misses**: phase deviation 48° (period 12 mm) / 58° (period free), tolerance step errors 58–71° |
| **LP-1** (2 diodes, MADP) | best of 16 starts | **misses**: phase deviation 79°, ladder order destroyed under the diode corners |

## CP-2 (lead CP candidate)

Stage 1 (single-state PB condition, stub-0 diode ON, three OFF: x′/y′ reflection-phase difference 180° ± 10°, |r| ≥ −1 dB, widest band) → stage 2 (full 4-state objective incl. θ = 30° and
the tolerance corners, Nelder–Mead on circumradius / line / stub width at period 15 mm) → mesh-90 polish → finalist correction (aligned meshes 60/90/120, pitch 0.25/0.167/0.125 mm,
observed-order extrapolation, R2 practice; corners on meshes 60/90 corrected by the nominal extrapolation ratio per state). The period went to the top of the allowed 6–15 mm range:
stub extension is the limit (an axis-aligned stub must stay clear of the cell boundary), so the period-15 optimum sits against a constraint.

**Finalist A** (stage-2 optimum; `cp2_madp_final.json`), corrected values, r_RR (handedness-preserving) in the propagation-relative basis:

| quantity | target | θ=0 | θ=30° | θ=60° (reported) |
|---|---|---|---|---|
| r_RR phases (states 1–4) | 0/90/180/270° ±15° | −104 / −10 / 98 / 178° — deviation **12.1°** | −101 / −9 / 94 / 176° — **8.8°** | −94 / −10 / 69 / 151° — 12.7° |
| min amplitude | ≥ −1.5 dB | −1.09 dB | **−1.60 dB (miss, 0.1 dB)** | −7.1 dB |
| amplitude spread | ≤ 1 dB | 0.16 dB | 0.43 dB | 5.6 dB |
| opposite-handed max(\|r_LR\|,\|r_RL\|) | ≤ −15 dB | **−13.6 dB (miss)** | **−14.0 dB (miss)** | −3.1 dB |
| step error, C_T −20 % / +20 % / L_s +0.1 nH (ASSUMED) | ≤ 20° | 15.7° / 17.0° / 18.0° | not run (cost) | — |

Bandwidth (nominal, θ = 0, search mesh): the **step error stays ≤ 17° over the whole 8.5–11.5 GHz sweep** (the phase is geometric, a PB rotation); amplitude ≥ −1.5 dB limits the usable band
to 10.0–10.75 GHz on the 0.25 GHz grid (**0.75 GHz, ~7 %**; the true edges are about 9.8 and 10.9 GHz). The limit is diode loss in the ON stub, not the ladder.

**Uncertainty — read before trusting the margins.** The observed Richardson order was unusable (non-monotone convergence on the staircased octagon) for 10 of 12 nominal values; those fall back to
first order. The resulting phase uncertainty is 3.5–8.3° per state for A. The raw fine-mesh values differ from the extrapolated ones by about as much as the margins (A, nominal θ=0 lattice
deviation: mesh 60 → 9.0°, 90 → 7.3°, 120 → 8.6° raw; 12.1° extrapolated). **The polished finalist B** (`cp2_madp_final2.json`, circumradius 4.605 / line 2.485 / width 1.372 mm — a 0.01 mm
change) has the same raw numbers but extrapolates to deviation 15.2° and tolerance step errors 25.6–34.3°: the two points are indistinguishable at the mesh uncertainty, and the choice of A
over B rests on the extrapolation noise. **Do not read the CP-2 margins as better than about ±10°.** HFSS (item 5) is the real test. The opposite-handed miss is robust: −13.5 to −15.2 dB on every
mesh and corner (the three OFF stubs break mirror symmetry).

**SMP1345 comparison** (same geometry, NOT re-optimised for SMP1345, mesh 60): lattice deviation 35°, amplitudes −2.8 to −7.9 dB, opposite-handed −1 to −3.6 dB, tolerance step errors 55–79°. As
found in item 2 the SMP1345 OFF state (0.19 pF + 0.7 nH, ≈ −40 Ω at 10 GHz) is not an open; this cell needs a different design for it, not a retune. MADP is the CP-family diode (roadmap log).

**Caveats for MADP:** it is a flip-chip part, flip-chip assembly is not yet confirmed with fabrication; if impractical, re-optimise CP-2 around a packaged low-capacitance PIN diode
(C_T ≤ ~0.05 pF, low L_s) — logged in the roadmap. L_s is NOT_AVAILABLE_FROM_CURRENT_DATA: nominal L_s = 0 is an assumption; the L_s corner is the one-sided ASSUMED offset +0.1 nH. The tolerance
requirement on C_T is exactly the datasheet maximum (0.030 pF = +20 %).

## LP families

Search effort: LP-1 5 + 8 starts (+3 early 150-evaluation starts), LP-2 6 + 8 (+3 early), LP-3 8 + 7 refinement restarts; Nelder–Mead, 280–500 evaluations each, the second round with the period free
(10–15 mm) and a phase-spreading term (the max-deviation objective is flat when the states coincide). The first LP-1/LP-2 rounds (random starts, tolerance corners in the objective) ended with
the four states at nearly the same phase; they are kept in `lp1_madp.json` / `lp2_madp.json` as a record.

| cell (MADP, θ=0 nominal) | phases (deg) | amplitudes (dB) | xpol (dB) | lattice dev θ=0 / 30° | step errors C_T−20 / +20 / L_s | 
|---|---|---|---|---|---|
| LP-2, period free (`lp2_madp_report.json`) | 116 / 86 / 50 / −104 | −0.01 / −0.25 / −0.88 / −1.39 | −23 | 58° / 73° | 71° / 58° / 66° |
| LP-2, period 12 mm (`lp2_madp_stage1.json`, nominal only) | 97 / 12 / 39 / −114 | −0.02 / −0.77 / −0.45 / −0.78 | −23 | 48.5° / — | not run |
| LP-1 (`lp1_madp_report.json`) | 129 / −158 / 84 / 168 | −0.00 / −1.24 / −1.24 / −0.28 | −53 | 79° / 86° | 166° / 185° / 183° (ladder order lost) |
| **LP-3**, codes 000,010,011,110 (`lp3_madp_report.json`) | 120 / 69 / −131 / −49 | −0.01 / −0.56 / −1.56 / −1.74 | −161 (mirror-symmetric) | 21.7° / 28.4° | 49° / 30° / 37° |

(Search mesh only; no finalist mesh correction was run for the LP cells — they miss by far more than the mesh uncertainty.)

**What the runs show (an observation about these searches, not a proof of impossibility):** in nearly every LP run the all-OFF state sits at 120–137° (a bare patch with open diodes) and the ON states
move only part of the 270° required, so the switched stubs do not move the patch resonance far enough. (Rough scale from the CP-2 x-polarisation slope, ~65° per 0.5 GHz near resonance: 270° needs the resonance
to move by roughly ±1 GHz, ~10 %; this is an estimate, not a derived limit.) For the LP-2 cell the C_T ±20 % corners barely move the states (OFF stubs are effectively open), so its tolerance
errors are the nominal miss, not tolerance sensitivity. The 2-diode cells (4 states, no choice) fail by 48–80°; the 3-diode cell (8 codes, best 4 picked) comes within
~7° of the phase target but misses the amplitude, spread and tolerance targets. A fourth diode or a different topology might close it; that has not been tried.

## Decision needed

LP is the open point. Options for the user: (a) continue LP-3 (more starts, 4 diodes, wider topology — each is another optimisation round), (b) accept LP-3 as the best 2-bit LP cell with its
reported misses, (c) revisit the scope decision (LP from the CP-2 cell via LP→co/cross codes was rejected on 2026-10-19; the CP-2 cell is the only one that meets the phase targets so far).
