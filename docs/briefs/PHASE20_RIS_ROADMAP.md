> **RE-ACTIVATED 2026-10-16** (user: "our ultimate goal is to move to RIS"). Parked 2026-10-13 in Phase 24, never
> started. Order: Phase 27 (tensor-CP flow + cell-family registry in the app, `PHASE27_APP_TENSOR_CP.md`) first, then
> **R1** (`PHASE20_R1_SOLVER_LOADS_OBLIQUE.md`). See "What Phase 26/27 give RIS" below before starting R1.

# Phase 20 roadmap — 2-bit PIN-diode reconfigurable metasurface (RIS), X-band

**Status:** roadmap, written 2026-10-09, updated 2026-10-16. Each item below becomes its own brief; only **R1** is active
(`PHASE20_R1_SOLVER_LOADS_OBLIQUE.md`). Scope is deliberately limited to **designing the product and the solver**:
no FPGA/codebook/MATLAB export, no fabrication files beyond what verification needs.

## What Phase 26/27 give RIS (added 2026-10-16)
1. **Transverse wavenumber in the MoM is already exercised.** Phase 26 B2 extracted Bloch tensors at real kt > k0 with
   the impressed-field J/E_tot approach and *no solver change* (`src/dispersion/tensor_sheet.py`; equals the
   normal-incidence tensor at kt = 0 to 1e-8). R1 item 1 should first confirm what of the (kx0, ky0) path through
   assembly/excitation/currents already works, then add only what is missing (oblique TE/TM port modes, the (0,0)
   field decomposition at oblique incidence, config). Verify, don't assume.
2. **Full 2×2 reflection is already handled** (tensor extraction from the 2×2 Γ). RIS state libraries (R4) should store
   the full 2×2 Jones matrix per state (co- and cross-pol, TE/TM at oblique incidence), not a scalar phase.
3. **App registry.** Phase 27 adds a cell-family / state-library registry whose schema accepts N-state switchable
   libraries; R4 registers its library there and R5 builds the RIS design mode on it.
4. **Locality lesson (Phase 25).** Abrupt neighbour changes (the 1-bit 0.55 mm-gap pair) broke the isolated-cell model;
   a 2-bit RIS coding map is abrupt by nature. R6's supercell/coupling check is therefore **moved forward**: run a
   locality test on the R3 cell (random and periodic 2-bit codings, supercell MoM vs isolated-cell prediction) before
   R5 relies on the library.
5. **HFSS loop.** The user runs HFSS; every phase that predicts a number writes the prediction before the HFSS run and
   provides a comparison script (Phase 26 pattern).

## Target product — starting values, not locked
Values come from the ISRO/IIT Indore technical proposal and are **starting points**. Any phase may change them
when its results call for it (e.g. no asymmetric cell reaches four states ~90° apart at 10 mm), provided it
(1) shows the result that forced the change, (2) records old → new value and the reason in a "Design changes" log
at the bottom of this file, and (3) re-checks anything downstream that depended on the old value. Only the
frequency band (X-band, ~10 GHz) and the technology (PIN-diode, 2-bit, PCB on a grounded substrate) are fixed
unless the user says otherwise.

| Item | Starting value | May change within |
|---|---|---|
| Frequency | 10 GHz (λ0 = 30 mm); characterise 8–12 GHz | fixed (X-band) |
| Cell period | 10 mm (λ0/3) | ~λ0/5 … λ0/2 (6–15 mm); upper bound for grating-lobe-free 60° steering at 10 GHz is p < λ/(1+sin60°) = 16.1 mm |
| Aperture | 300 × 300 mm (30 × 30 cells at 10 mm) | follows period and gain/beamwidth needs; cell count = aperture / period |
| Substrate | Rogers RT/duroid 5880, εr 2.2, tanδ 0.0009, 1.575 mm, PEC ground | other standard 5880 thicknesses (0.787, 3.175 mm) or another Rogers laminate if phase range/bandwidth requires |
| Cell | patch with two PIN diodes across gaps at the two ends, one via patch-centre → ground (DC bias return); bias lines on a layer behind the ground | topology may change (diode count/positions, slots) if 2-bit isn't reachable |
| Diode | Skyworks SMP1345-079LF (SC-79); ON = R_s + L_s, OFF = C_T + L_s from the datasheet, labelled provisional | any X-band PIN diode in a pick-and-place package |
| States | 2-bit: four diode combinations → four phases ≈ 0/90/180/270° | fixed (2-bit) |
| Polarisation | **LP and CP as user-selectable cell families** (decided 2026-10-19) | dual-pol later |
| Steering | beam 0°–60°; illumination: plane wave and offset horn ≈ 25°, F ≈ 240–300 mm (R5) | horn offset/distance free |

## Engineering constraints carried into every brief
1. **Symmetric two-diode cells give 3 states, not 4** (ON/OFF ≡ OFF/ON by mirror symmetry). The cell must be
   deliberately asymmetric so all four combinations are distinct.
2. **Centre via**: at the patch's fundamental-mode voltage null, expected RF-benign. The MoM is planar and cannot model
   vertical currents, so the via is omitted and labelled `VIA_OMITTED_APPROXIMATION`; one HFSS solve (run by the user)
   quantifies it.
3. **Bias layer behind the ground** is shielded and outside the RF model — state it, don't model it.
4. Repo rules apply: no silent changes to frozen solver modules (additive, documented, bit-identical for every
   pre-existing case, like Phase 12A), physics tests, honest labels, no ML in the solver core.

## Phases
| # | Phase | Output |
|---|---|---|
| **R1** | Solver: lumped RLC loads across gaps + oblique incidence (TE/TM, θ ≤ 60°) + a composite-rectangle geometry for gapped patches | extended MoM, tests |
| R2 | Diode models (PIN ON/OFF; varactor C(V) optional) | `src/ris/diodes.py` |
| R3 | 2-bit cell design: asymmetric two-gap patch on 5880, 10 mm; optimise for four states ~90° apart, balanced \|Γ\|, stable to 60° both pols; **and robust to diode tolerance: the objective must include ±20 % on C_T and L_s (R2 finding: the OFF state sits near its series resonance, 13.8 GHz, and moves 20–100° for ±20 %), not only four states 90° apart**; **locality test of the cell under 2-bit codings (moved forward from R6)**; one HFSS single-cell check package (user runs it) | cell geometry + state responses |
| R4 | State library vs f (8–12 GHz), θ, φ, pol (full 2×2 per state); metrics: phase error, magnitude spread, 90°-step bandwidth, angular stability; registered in the Phase 27 registry | library + report |
| R5 | RIS design mode in the app (on the Phase 27 registry): plane-wave and offset-horn illumination, 2-bit coding map from each cell's real oblique response, quantisation-lobe report (+ random-offset suppression), predicted pattern/gain; HFSS export of a coded board | design module + app path |
| R6 | Coupling check on real 2-bit coding patterns at array scale (supercell / full-array tools), building on the R3 locality test | report |

## Design changes log
| Date | Phase | Item | Old → new | Result that forced it |
|---|---|---|---|---|
| 2026-10-16 | roadmap | R6 locality check | after R5 → first part moved into R3 | Phase 25: isolated-cell model fails for abrupt neighbour changes |
| 2026-10-19 | R3 | polarisation | single linear → **LP and CP as user-selectable families** (two cell families on the same technology; R5 lets the user pick) | user decision at the R3 decision point (the HMSA work requires pure CP; the RIS must serve both) |
| 2026-10-20 | R3 | CP-family diode | SMP1345-079LF (R2 comparison) → **MADP-000907-14020** for CP-2 (SMP1345 kept as the comparison case in every report) | R3 item 2: the SMP1345 OFF state is not an open at 10 GHz (0.19 pF + 0.7 nH ~ -40 ohm) and breaks the CP-2 PB ladder (steps 62/197/59/106 deg); MADP (0.025 pF) gives 91/99/92/78 deg unoptimised. **Caveat:** MADP is a flip-chip part; flip-chip assembly is NOT yet confirmed with fabrication. If impractical, re-optimise CP-2 around a packaged low-capacitance PIN diode (target C_T <= ~0.05 pF, low L_s) |
| 2026-10-21 | R3 | LP scope | 2-bit LP cell (LP-1 / LP-2 / LP-3) → **LP as 1-bit codes of the CP-2 cell** (LP→co-LP states 1/3, LP→cross-LP states 2/4) for now; **2-bit LP deferred to a future phase needing a multilayer / slot-coupled cell** | R3 item 3: no 2-diode LP cell comes within 48° of the 15° phase target and the 3-diode LP-3 misses amplitude, spread and tolerance; the MADP-000907 diode (OFF ~ open) was used in ALL item-3 LP searches (SMP1345 only in the item-2 starting scans), so the SMP1345 OFF loading was not the cause |
| 2026-10-21 | R3 | for R5 | R5 as planned → **R5 must verify codings at array level (supercell / coupling-aware), not trust isolated-cell predictions** | R3 item 4: isolated-cell predictions are reliable for regular gradients (0.1 dB, 3 deg on a 0-1-2-3 gradient) but not for random codings at the 15 mm pitch (strongest-order levels off 1-4 dB, single orders up to ~10 dB, phases within 10 deg; gate passed 2 of 10 random cases); `src/ris/supercell.py` is the starting tool |
| 2026-10-21 | R3 | open item | pitch vs stub length | the CP-2 optimum sits ON the cell-edge limit (stub line parameter 2.48 mm clamped to 1.748 mm at the 15 mm maximum period): longer stubs would need a larger pitch (quantisation lobes, locality) or a meandered stub; to be settled if HFSS confirms CP-2 |
