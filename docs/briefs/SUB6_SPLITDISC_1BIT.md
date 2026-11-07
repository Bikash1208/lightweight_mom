# SUB6_SPLITDISC_1BIT - sub-6 GHz redesign: split-disc tensor cell, 1-bit CP board (2 beams OFF / 1 beam ON)

Status: **active brief, written 2026-11-07 (decision session).** MODEL_ONLY until HFSS returns. Read `CLAUDE.md`,
`docs/START_HERE.md`, `docs/STATUS.md` first. All repo rules apply (one commit per item, fast suite green, physics tests,
honest named states, never fabricate data, additive-only solver changes, heavy compute on the server only).

Two sessions: the **decision session** (this brief, decisions, reviews) and the **server session** (all MoM runs, scripts,
commits on `compute/<job>`). The server runs exactly what a `jobs/<job>.md` says and reports. It does not decide the next stage.

Branch: create `sub6-splitdisc` from the current head of `sidequest-hybrid-hms-ris` (it holds the rotated-load fix
`fcc6691`, which rotated diodes need). Record the base commit hash in the first job report.

## 1. Requirement (user, 2026-11-07)
- Frequency: sub-6 GHz, 4-6 GHz; exact band chosen from the application study (section 2).
- Cell: **two half-circles (split disc) joined by one PIN diode across the gap.** Degrees of freedom: cell orientation
  psi and half-circle radius, either one radius (variant S) or two radii r1, r2 (variant T). The results decide which.
- First deliverable: a **1-bit board**. Control circuit OFF gives **two CP beams (A and B)**. ON gives **one CP beam = beam A**
  (same direction as OFF beam A; B suppressed). All beams CP.
- Board substrate, size and feed: open; proposal in section 3, user decides.
- Rule from the user: **no assumptions carried as results**. Every architecture choice below waits for its gate.

## Decisions log
- 2026-11-07: **frequency 5.8 GHz (ISM 5.725-5.875) confirmed** by the user.
- 2026-11-07: **beams A (theta 25, phi 0) and B (theta 25, phi 180)** (user: 25 deg instead of 30 deg).
- 2026-11-07: **bias = ideal**: one lumped RLC sheet between the two halves in HFSS, same as the MoM lumped load; no vias or bias lines
  in the cell or board models. The HFSS sheet has the **same footprint as the MoM load (gap g x 1.0 mm along the gap, at the centre)**
  so MoM and HFSS describe the same device. Real bias hardware is a later item; until then results carry `IDEAL_BIAS`.
- 2026-11-07: **gap and sheet footprint come from the diode data sheet** (user), not from the design: MA4AGP907 g 0.33 mm, sheet 0.33 x 0.35 mm
  (from the MA4GP907 outline, AGP outline UNVERIFIED); SMP1345-079LF g 0.6, sheet 0.6 x 0.3 mm (ESTIMATE, no land pattern in the data sheet).
  Details and sources in `jobs/sub6-splitdisc-cell-scan.md`. The earlier "g x 1.0 mm" footprint is superseded.

- 2026-11-07: **diode budget Rs 20,000 for the whole project (taken as INR), diode must be purchasable in India** (user). See section 2b.

## 2b. Diode choice under the budget (decision session, 2026-11-07; prices from distributor listings via Octopart / Findchips, USD,
## converted at about Rs 95/USD (Aug 2026 close) with an ASSUMED x1.3 for customs duty + GST + shipping; verify on DigiKey/Mouser India before ordering)
| diode | package | lifecycle | price (qty 100-1000) | landed each (est.) | diodes for Rs 20k (15 % spares) |
|---|---|---|---|---|---|
| MA4AGP907 | flip chip | active | $6.92-8.30 (Mouser) | Rs 855-1025 | **17-20: ruled out** (and needs flip-chip assembly) |
| SMP1345-079LF | SC-79 | active, large stock at DigiKey/Mouser (both sell in INR in India) | $1.77-2.26 | Rs 219-279 | **62-80** |
| BAR64-02V | SC-79 | OBSOLETE (aftermarket stock: Arrow, Avnet, Rochester) | $0.054-0.080 | Rs 7-10 | 1760-2600 |
| BAR63-02V | SC-79 | OBSOLETE, last-time-buy 2026-06-30 (Rochester 15k stock) | $0.052-0.090 | Rs 6-11 | 1565-2700 |
| BAP64-02 (Nexperia) | SOD-523 | OBSOLETE, LTB 2026-03-31, little stock | $0.11-0.17 | Rs 14-21 | 820-1235 |
Electrical at 5.8 GHz is similar for all SC-79 parts (C 0.19-0.3 pF, 1.5-2.1 ohm at 10 mA, Ls 0.6-0.7 nH); the cell scan covers SMP1345 and
BAR64-02V explicitly. **The budget, not the physics, sets the diode count:**
- with SMP1345 the board can switch at most about 70 cells. On a 1039-cell board that is 7 % of the aperture, too few for a beam B
  of similar strength to A (the 10.5 GHz hybrid needed 90-180 switched cells for a much weaker beam 2);
- with a BAR63/64-class part (one-lot purchase of aftermarket stock, acceptable for a prototype if the lot is bought at once) the budget
  covers the whole B group (about 500 diodes) with margin.
Open (user): accept an end-of-life part bought as one lot? Diode count cap for the board model: **<= 70 (SMP1345)** or **<= 600 (BAR6x)**.

**Decision 2026-11-07 (user): Infineon BAR50-02V (SC-79).** Data sheet (Infineon BAR50 series, typ / max): C_T 0.15 pF at 0 V, 1-1.8 GHz
(0.2 pF at 0 V / 100 MHz; 0.24 / 0.5 pF at 1 V, 1 MHz); r_f 3 / 4.5 ohm at 10 mA, 16.5 / 25 ohm at 1 mA (100 MHz); L_S 0.6 nH. Lower 0 V capacitance
than BAR64 (good for an unpowered OFF state), higher ON resistance (more ON loss; bias at >= 10 mA). Footprint values in the extracted data sheet
text (0.35 / 1.35 / 0.35 mm) suggest pads 0.35 mm with a 0.65 mm gap: **g = 0.65 mm, sheet 0.65 x 0.35 mm, ESTIMATE until read from the
drawing.** Procurement (Octopart / Findchips, 2026-11-07): **OBSOLETE, last-time-buy 2026-06-30; authorised distributors (Avnet, Arrow,
Future) show 0 stock**; brokers only: Win Source 227k at $1.17-1.76 (landed about Rs 145-215, so about 80-120 diodes for Rs 20k), Component
Stockers 2182 at $0.05-0.26 (Rs 6-32, so about 600-700 diodes fit). Broker parts carry a counterfeit risk and India shipping is unconfirmed:
**the diode count cap waits for the actual quote**. SMP1345-079LF stays the in-production fallback.

## 2. Frequency: 5.8 GHz ISM (5.725-5.875 GHz), confirmed 2026-11-07
| option | band | for | against |
|---|---|---|---|
| **5.8 GHz ISM (proposed)** | 5.725-5.875 (2.6 %) | delicensed in India and ISM worldwide; standard band for microwave power transfer and point-to-point links, where CP and switchable beams are used; narrow band, so beam squint is about 1 deg (29.5-30.5 deg for a 30 deg beam, 1.06 index, model estimate) | crowded with Wi-Fi |
| Wi-Fi 5 GHz | 5.15-5.85 (12 %) | large market | 12 % band: squint 27.6-32.1 deg, about half the 8.7-11.7 deg beamwidth; one hologram cannot serve the whole band well |
| 5G n79 | 4.4-5.0 | 5G | not allocated in India (India uses 3.3-3.67 GHz); 13 % band |
| lower 6 GHz | 5.925-6.425 | delicensed in India from Jan 2026 (Wi-Fi 6E/7) | outside the 4-6 GHz window the user gave |

Effects of the move from 10.5 to 5.8 GHz (computed, `sizing.py` in the decision session):
- Diode OFF reactance at 5.8 GHz: MA4AGP907 chip -974 ohm (was -538), MADP-000907 -1098 ohm, **SMP1345-079LF (packaged) -119 ohm**
  (was -34, series resonance 13.8 GHz): **a packaged PIN may now be usable**, which the via-cell study ruled out at 10.5 GHz.
  The scan includes it as a state.
- The OFF diode capacitance competes with the gap capacitance of the split, so the OFF state is not "open". The scan uses
  the real OFF models, not only an open circuit.

## 3. Board proposal (user decides; numbers scale the HFSS-measured hybrid board, 150 mm side hexagon, 20.5 dBi, aperture efficiency 12.5 %)
| hexagon side | across corners / flats | Dmax at 5.8 | expected, one beam (ON) | each of two beams (OFF, -3 dB split) | cells at p = 10 mm |
|---|---|---|---|---|---|
| 150 mm (same board) | 300 / 260 mm | 24.4 dBi | about 15.3 dBi | about 12.3 dBi | 585 |
| **200 mm (proposed)** | 400 / 346 mm | 26.9 | about 17.8 | about 14.8 | 1039 |
| 230 mm | 460 / 398 mm | 28.1 | about 19.1 | about 16.1 | 1374 |
| 272 mm (same electrical size as 10.5 GHz) | 544 / 471 mm | 29.6 | about 20.5 | about 17.5 | 1922 |
"Expected" assumes the 12.5 % aperture efficiency HFSS gave at 10.5 GHz. That is an extrapolation, NOT a prediction for this cell.
- 200 mm side fits a standard 18 x 24 in (457 x 610 mm) laminate panel; 272 mm does not (flats 471 > 457 mm).
- Substrate: **RT5880 3.175 mm** keeps the electrical thickness of the validated 1.575 mm @ 10.5 GHz stack (h/lambda 0.061 vs 0.055);
  1.575 mm @ 5.8 GHz is thin (0.030 lambda) and binds the surface wave weakly. The scan runs **both** thicknesses; the scan decides.
- Feed: one centre monopole probe (SMA) as before, length scaled; launcher redesigned after the cell is fixed.
- Cell period: about 10 mm (lambda0/5.2) baseline; scan p = 10 mm (9 and 11 mm only if the scan's resonance check asks for it).

## 4. The physics point that decides the architecture (hypothesis H1, to be tested, not assumed)
Variant S and variant T both keep a mirror line perpendicular to the gap. If that symmetry holds in the solver, the 2x2 sheet tensor of the
cell is diagonal in the gap frame **in both diode states**: principal axes along / across the gap, fixed by psi. Switching
changes only the two principal values (Y_par, Y_perp), not the axes. For a CP tensor hologram the cell's radiated contribution is set
by its anisotropy Delta = Y_par - Y_perp, and its phase by 2 psi. Then for every cell **m_ON = rho * m_OFF** with rho = Delta_ON / Delta_OFF real.
The diode can only scale a cell's contribution or flip its sign; it cannot rotate its phase.
Consequence (if H1 holds): OFF = A + B -> ON = A is not a per-cell phase change. It must be an **amplitude change by cell groups**:
- **Architecture I (interleaved, preferred if available):** A cells with rho near 1 (diode has little effect, or no diode at all) phased for beam A;
  B cells with **rho near 0** (diode ON removes their anisotropy) phased for beam B. OFF: both groups radiate (2 beams). ON: B cells go
  quiet and their surface-wave power continues to the A cells (beam A should get stronger). Diodes only on B cells.
- **Architecture II (sign flip):** if some geometry gives rho near -1, cells can encode A + B in OFF and A - B (or A) in ON with a per-cell fit.
- If rho stays between about 0.5 and 1 for every geometry, neither works well; then the cell must change (decision session).
Also check: the **mean** admittance changes with the state (beta shift); B-cell switching shifts the local beta seen by beam A (phase error in ON).
Gates in the job report rho range, mean shift and symmetry. **No board design until the decision session reads them.**

## 5. Bias (DECIDED 2026-11-07: ideal RLC sheet, see the decisions log; the text below is kept as the record of the later real-bias item)
The two half-discs are DC-isolated by the gap, so each half needs its own DC path: (a) a via from one half to ground and a via
through an anti-pad to a bias layer for the other (as in `SIDEQUEST_VIA_BIAS`; **the MoM has no vertical currents, so vias are
HFSS-only**); or (b) thin high-impedance top-layer lines with ballast resistors (in-plane, so the MoM can model them, but they perturb the tensor).
The first scan is the bare cell (no bias hardware) and is labelled `BARE_CELL_UPPER_BOUND`.
Lesson from the 10.5 GHz cell: HFSS did not reproduce the MoM's single-diode switching (`SIDEQUEST_HANDOVER.md` sec. 2, unresolved).
So an **HFSS cell eigenmode check comes before any board** here too.

## 6. Stages (stop points marked ■)
1. **Cell scan (server job `jobs/sub6-splitdisc-cell-scan.md`).** ■ Decision session reads the gates, picks variant S/T, thickness, geometry
   window, architecture I/II, and confirms frequency + board size with the user.
2. HFSS cell package (eigenmode OFF/ON, diode as one lumped RLC sheet between the halves, MoM footprint), user runs HFSS. ■ Compare, calibrate.
3. Board model: group assignment / joint optimisation with the Phase 26 tensor engine (`src/hologram/tensor_leaky_cp`); check first whether
   `src/hologram/shared_aperture_tensor.py` and `src/hologram/two_state_cells.py` already cover architecture I; launcher at 5.8 GHz. ■
4. Board HFSS package, user runs it. ■

## 7. Targets (first draft, the user may change them)
OFF: beams A and B both CP, AR <= 3 dB, balance |G_A - G_B| <= 3 dB. ON: beam A, B suppressed by >= 12 dB relative to A.
SLL <= -13 dB both states. Beam directions (user, 2026-11-07): **A (theta 25, phi 0), B (theta 25, phi 180)**; handedness RHCP for both
unless the user says otherwise.
