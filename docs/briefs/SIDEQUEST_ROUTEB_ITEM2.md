# SIDEQUEST_ROUTEB_ITEM2 - candidate B (back-switched mushroom, 1 diode/cell): cell + board design, MODEL_ONLY (2026-11-06)
**Everything here is MODEL_ONLY.** Cell: MoM Bloch sheet tensor of the patch grid (grounded RO4003C, no via) + a lumped via/diode branch per square cell (Sievenpiper/Luukkonen mushroom model, L_via = mu0 h; the via is NOT in the MoM). Board: SWHB Born model (`sidequest_swhb.py`); **M6 is NOT_MET** (`SIDEQUEST_M6_VALIDATION_SET.md`), so board numbers carry wide bands (gain +-3 dB, SLL +-3 dB, pointing +-3 deg) and rim diffraction is not modelled. Nothing is HFSS-checked. Scripts: `scripts/research/sidequest_cellB_*.py`, `sidequest_swhb*.py`; data `results/sidequest_hybrid/cellB_*.json`, `swhb_scan_*.json`.

## 1. Correction of the earlier concept numbers
The lumped estimate in the item-1 write-up (patch 7 mS, via 0.7 nH, resonance 18.5 GHz) was wrong. MoM patch grid at p 6.4: Y_patch = +30 mS at gap 0.4, +9-10 mS at 0.8; the dense-lattice via inductance is mu0 h (1.9 nH at 1.524 mm), not 0.7 nH. With the ON via the top sheet is capacitive-minus-inductive (+20 mS) for narrow gaps: the TM wave is bound only when the net sheet susceptance is small (4-10 mS), which fixes the gap (0.6-1 mm) and makes the cell a **high-impedance surface operated below its stop band, not far from it** (see 3).

## 2. Recommended cell (D2) and results
RO4003C (er 3.38 +-0.05, tan d 0.0027), h 1.524 mm, p 6.4 mm, square patch 5.76 mm, gap 0.64 mm (= 0.1 p, allows meshes 40/60/80), central via to a back island (diameter free in this model; the +-0.05 mm drill is bounded as +-8 % of L_via), diode MA4AGP907 measured (ON 2.44 ohm, OFF 28.16 fF) on the back between island and ground, **knob: series inductance of the back trace L_ser = 1.0 nH** (about 1.7 mm of 0.2 mm line at an assumed 0.6 nH/mm; ESTIMATE). Knob scan (mesh 60): L_ser 0 / 0.5 / 1.0 / 1.5 / 2.0 nH -> A0 26.8 / 17.9 / 12.2 / 8.9 / 7.0 Np/m, beta_ON 1.059 / 1.086 / 1.112 / 1.135 / 1.154: the contrast knob trades alpha against beta_ON; the window 12-18 Np/m and beta_ON >= 1.10 overlap only at 0.75-1.0 nH.
| quantity | value |
|---|---|
| beta/k0 OFF / ON (x, mesh 80) | 1.395 / 1.113 (45 deg, mesh 60: 1.418 / 1.130) |
| convergence mesh 40 -> 60 -> 80 | beta_ON 1.1093 / 1.1120 / 1.1134 (0.24 %, 0.12 %); beta_OFF 1.3905 / 1.3934 / 1.3948 (0.21 %, 0.10 %); A0 12.56 / 12.15 / 11.95 Np/m (3.2 %, 1.7 %): M2 met |
| isotropy x vs 45 deg (mesh 60) | beta_ON +1.6 %, beta_OFF +1.8 % (limit 3 %); A0 12.15 vs 11.33 (-6.8 %); y = x by symmetry (checked at meshes 32/64 for the first cell) |
| contrast | dY = 8.09 mS; A0 (radiative alpha of the fully 1-bit modulated lattice, first order, theta0 25 deg, with the MoM patch Y at kx = k0 sin 25) = 12.0 Np/m, at the bottom of the 12-18 target |
| absorption | alpha_loss ON 0.71, OFF 0.62 Np/m (diode R + dielectric; copper not included) -> about 5 % of A0 |
| cell cross-pol | Y_uv/Y_uu = 0 (C4v mirror symmetry); NOT a verified number (diode / pad asymmetry not modelled) |
| resonance distance | ON: normal-incidence cell resonance (Im Y_tot = 0) at **11.27 GHz (+0.77 GHz)**; OFF: the TM root disappears between 11.0 and 11.5 GHz (stop band edge, **+0.5-1.0 GHz**), series resonance 7.98 GHz (-2.5 GHz): **the 1.5-2 GHz requirement is NOT met** (dbeta_OFF/df about +0.2 per GHz) |
| mean-alpha target | ln10/(2R): 15.4 Np/m (R 75), 11.5 (100), 9.2 (125), 7.7 (150); the board designs below reach the edge targets with partial modulation, A0 is the full-modulation ceiling |
| back (anti-pad) leakage | ESTIMATES: Bethe small hole (b 0.9 mm) radiating into the half space: alpha_back about 1e-4 Np/m; monopole (island + diode, 0.2-1 mm) 0.02-0.5 ohm against R_on 2.44 ohm -> <= 0.1 Np/m; dielectric loss of a 25 fF pad C: 5 uS equivalent. Negligible against A0 12; the back-cavity (bias plane) coupling and bias-line radiation are NOT estimated (needs HFSS) |

### Tolerance (mesh 60, x; single factor; d = change from nominal)
| factor | d beta_ON | d beta_OFF | d A0 |
|---|---|---|---|
| etch +0.03 / edge (gap +0.06) | -0.0125 | -0.0145 | +15.5 % |
| etch -0.03 / edge (gap -0.06) | +0.0189 | +0.0178 | -20.5 % |
| er 3.33 / 3.43 | -0.0073 / +0.0079 | -0.0153 / +0.0154 | +3.2 / -3.6 % |
| diode R +-20 % | 0 | 0 | -+0.3 % |
| **diode C +-20 %** | 0 | **-0.0230 / +0.0261** | -7.8 / +5.0 % |
| via L +-8 % (drill / plating bound) | -0.0117 / +0.0110 | -0.0034 / +0.0036 | +18.6 / -14.3 % |
| L_ser +-15 % (trace etch) | -0.0075 / +0.0072 | -0.0022 / +0.0023 | +11.4 / -9.6 % |
| **pad capacitance +-5 fF (20 % of an assumed 25 fF; NOT in the diode model)** | 0 | **+0.0230 / -0.0206** | +4.7 / -6.8 % |
| geometry +-0.05 mm (gap) / h +-0.05 mm (patch Y held) | -0.0105 .. +0.0156 / -0.0145 .. +0.0159 | -0.0120 .. +0.0149 / -0.0191 .. +0.0179 | +13 / -17 % ; +13.5 / -13.3 % |
| via-L model (kappa 0.7 instead of 1) | **-0.0300** | -0.0079 | **+56 %** |
Combined over the spec factors (etch, er, diode R, diode C, via L, L_ser): worst case (sum) d beta_ON 0.046, d beta_OFF 0.065, d A0 62 %; RSS 0.025 / 0.035 / 31 %. **Worst-case beta_ON = 1.073 (>= 1.05 met).** **R5 / criterion 4 NOT met** (dominant: diode C and pad C in the OFF state, etch; the OFF state sits on the steep side of the dispersion). The via-inductance model (kappa) is the largest unknown; the HFSS cell check (item 3) decides it.
Important: the pad / island capacitance to ground is **not** in the measured diode model (28.16 fF); a 0.3 mm pad over a 0.1 mm bond ply is about 25 fF, i.e. as large as the diode. With +20 fF the OFF root leaves the window (g 0.8 cell). The geometry must be re-centred once the pad is drawn.

### Static cells (diode-free)
Via + series trace, no diode: beta/k0 = 1.044 (hard via, L_tot 1.9 nH), 1.084 (3 nH), 1.154 (6 nH), 1.201 (12 nH), 1.226 (24 nH), 1.25 (no via) for the g 0.8 cell. **The static cell with the same mean beta as the pair (1.23) needs L_ser about 27 nH (about 45 mm of 0.15 mm trace: a meander on the back layer) or an open via (beta 1.25, +0.02).** Practical static family: shorted via (ON-like), open (OFF-like), and meander values in between. Polarity: type A (ON lowers the sheet susceptance) is the cell above. Type B (ON raises) is a CONCEPT not simulated: permanent inductive via plus a diode in series with a MIM capacitor to ground (ON adds a capacitive branch); it is used in the board designs only as a modelled option (susceptance mirror of A).

### Diodes, laminates, bias
* **MA4AGFCP910** (MACOM datasheet V4, copy at datasheet.chipsmall.com/MA4AGFCP910/334350.pdf; not the macom.com original): C_T 0.018 pF typ / 0.021 max (-10 V, 1 MHz and -5 V, 10 GHz), R_S 5.2 ohm typ / 6.0 ohm max (+10 mA, 10 GHz), V_F 1.33 V typ / 1.45 V max (+10 mA), reverse voltage -75 V, carrier lifetime 4 ns. In the same cell geometry (L_ser 1.6 nH): beta_OFF 1.311, A0 9.0 (typ) / 10.0 (max) Np/m, alpha_loss ON 0.78 / 0.86 Np/m: lower contrast (smaller C_off) and 1.6x the ON loss; usable only with L_ser -> 0 and a gap re-tune. Baseline stays the measured MA4AGP907.
* **Laminates** (Dk, tan d from rogerscorp.com product pages: RO4003C 3.38 +-0.05 / 0.0027, RO3003 3.00 +-0.04-0.047 / 0.0010, RO3006 6.15 +-0.15 / 0.0020 at 10 GHz; thickness lists from third-party mirrors of the Rogers data sheets, to be checked: RO4003C 60 mil = 1.524 mm, RO3003 60 mil, RO3006 50 mil = 1.28 mm): **RO4003C 1.524 chosen.** RO3006 1.28: no TM root for gaps <= 2.0 mm, at 2.4 mm roots exist but A0 5.6 Np/m, and Dk +-0.15 (2.4 %) is the worst tolerance. RO3003 1.52: roots, lowest loss, but A0 <= 7 Np/m at beta_ON >= 1.1 (gap 0.4, 1 nH) -> alpha target missed. Multilayer bonding: RO4003C with RO4450F/T prepreg is the standard RO4000 combination (used in the read journal paper Pan et al., IEEE IoT-J 2022, RO4350B + RO4450F); RO3000 laminates have a separate RO3000 bondply (Rogers processing guide, not read).
* **Stack (3 metal layers):** L1 patches (top), L2 ground (bottom of the 1.524 mm core), a thin bond ply (about 0.1 mm), L3 back layer: islands, diodes, chokes, bias tree. Total about 1.75 mm + copper. A 2-metal-layer version (islands and bias traces on the ground layer, slots cut in the ground) is not recommended (ground slots leak, continuity).
* **Bias (feasibility only):** 32 parallel strings of 16 series diodes = 512 diodes (<= 514), one supply line. V_F 1.33 V typ (MA4AGP907/FCP910 datasheet, +10 mA) gives 21.3 V per string typ, 23.2 V max: a 24 V supply has only 0.8 V headroom at max V_F -> per-string constant-current / series resistor ((24-21.3)/10 mA = 270 ohm typ, but the V_F spread makes 80-270 ohm): choose 26-28 V compliance. Current 32 x 10 mA = 0.32 A (<= 0.5 A met), diode dissipation about 6.8 W. Default = unpowered = OFF (model state 0 V). **Topology consequence:** a series chain needs the diode cathode NOT at ground: per cell a cathode pad RF-grounded by a bypass capacitor (overlap with L2 across the 0.1 mm ply: >= 2 pF -> pad about 2.6 x 2.6 mm, X = -7.6 ohm in series with the ON branch, absorbed into L_ser) and a DC link to the next cell's anode island through an RF choke. Choke estimate: lambda/4 high-Z line (0.15 mm, Z0 about 120 ohm, 4.2 mm) into the bypass: Z_in = Z0^2/|X_C| about 1.9 kohm, i.e. an inductive shunt of about -0.5 mS on every cell's branch (ON branch -4.3 mS, OFF +1.9 mS): **not negligible (7-25 %), must be put into the cell model and its tolerance; isolation per cell is limited to about 1.9 kohm** unless a higher-Z / longer / two-stage choke is drawn. R3 met: every diode terminal has a DC path; no top bias lines.

## 3. Board results (SWHB, MODEL_ONLY, 1-bit coding with static levels; cell D2 unless stated; launch efficiency not modelled -> gains are per launched power)
Design method: far-field-matrix coordinate descent (beam power, sidelobe power outside 1.3 HPBW + 3 deg of every beam, wrong-beam power), SW depletion iterated, diode budget by a price bisection, switched cell types A (ON lowers) and B (ON raises, concept), 14 static levels.
**Recommended: swap, circular board D = 300 mm (R 150), 514 diodes of 1644 cells (31 % switched), p 6.4:**
| | state 1 (OFF) beam (25, 0) | state 2 (ON) beam (25, 180) |
|---|---|---|
| gain / D [dBi] (band +-3) | 23.9 / 24.3 | 24.3 / 24.6 |
| pointing | (25, 0) | (25, 180) (grid 1 deg) |
| SLL, first-null definition [dB] | **-12.2** | **-12.5** |
| HPBW [deg] | 8.7 | 8.6 |
| XPD [dB] (no diode asymmetry) | 44 | 43 |
| level at the other beam / isolation | 4.5 dBi / 19.4 dB | -0.1 dBi / 24.4 dB |
| P_rad / edge / absorbed (fractions of the launched power) | 0.93 / 0.04 / 0.04 | 0.92 / 0.02 / 0.06 |
| S11 | NA (probe not modelled); edge reflection <= -14 dB | NA |
Other seeds: 512 diodes SLL -11.4 / -12.5; R 125 (D 250): 514 diodes G 22.7 / 23.1, SLL -10.5 / -9.7, edge 6 / 2 %. **Add-a-beam (b), D = 300 mm, 514 diodes:** state 1 G 24.2 dBi, SLL -11.3, edge 7 %; state 2 beams 19.5 / 18.6 dBi (balance 0.9 dB), SLL -13.2 (outside both lobes), level at beam 2 in state 1 -4.2 dBi -> beam-2 contrast 22.8 dB.

### Criteria (swap, D 300, 514 diodes; add-a-beam in brackets)
1 clean single beam, pointing <= 2 deg: **pass in the model** (band +-3 deg; the old-board HFSS split was a cell-beta error of +0.15, which this cell has not been checked against).
2 SLL <= -13: **FAIL** (-12.2 / -12.5; [state 1 -11.3, state 2 -13.2]).
3 gain >= 17 dBi: pass (23.9 / 24.3); radiation efficiency >= 85 %: pass (0.93 / 0.92 of the launched power, before the unknown launch efficiency); XPD >= 15: pass (43-44, model); isolation >= 15 dB: pass (19.4 / 24.4); S11 <= -10: UNKNOWN.
4 robustness: **FAIL** (cell: combined d beta 0.025 RSS / 0.065 worst, d alpha 31 % / 62 %; worst beta_ON 1.073 passes).
If no design meets them at <= 514 diodes: SLL and robustness fail. **Doubling the diodes does not fix the SLL in this model** (swap D 300 with 1016 diodes: -12.4 / -11.2; D 400 with 1026: -11.8 / -12.4; first cell g 0.8, info only), so the limit is the 1-bit quantisation of the hologram and my optimiser, not the diode count; I cannot name a diode count that reaches -13. Robustness is a cell property, independent of the diode count.

### Open questions
a. **Pitch.** p 6.4 (D2): converged root, beta window above. p 4.5 (patch Y at mesh 60 only, gap 0.3 mm = etch +-0.06 is +-20 % of the gap): roots exist, beta_ON 1.06-1.09, beta_OFF 1.29-1.33, A0 13-15 Np/m, but boards at <= 514 diodes: 22-35 % of the cells switched, edge power 38-95 %, SLL -8.3 .. -9.9 dB, three of six runs collapsed to 0 diodes. **p 4.5 is worse at 514 diodes.** Grating lobes: lambda/p = 4.5 at p 6.4 -> none; the quantised 1-bit pattern (3.8 cells per SW wavelength) puts its harmonic content into the diffuse floor (included in the SLL above), not into grating lobes.
b. **Full vs partial** (p 6.4, only possible for R <= about 80 mm: 396 cells, 514 diodes available): R 75 swap, full 396 diodes G 9.8 / 12.1 dBi, SLL -6.9 / -7.2; partial 211 diodes G 18.3 / 18.6, SLL -8.7 / -8.6; add-a-beam full G 14.3 / 17.3, partial 16.6 / 16.9. **Partial wins**: the analog static levels matter more than a diode in every cell (first cell g 0.8).
c. **Diameter** (g 0.8 cell, 514 diodes, swap, partial; G S1/S2 dBi, SLL S1/S2 dB, edge S1/S2, diodes/cells): R 75: 18.3 / 18.6, -8.7 / -8.6, 0.23 / 0.20, 211/396; R 100: 21.1 / 21.4, -9.6 / -9.2, 0.11 / 0.08, 431/712; R 125: 22.9 / 23.1, -10.0 / -9.7, 0.09 / 0.04, 511/1132; R 150: 24.3 / 24.4, -11.7 / -12.3, 0.08 / 0.04, 507/1644; R 200: 26.3 / 26.1, -10.0 / -11.5, 0.08 / 0.05, 500/2960. Switched fraction: 100 % (R 75), 61 % (R 100), 45 % (R 125), 31 % (R 150), 17 % (R 200) (the old 2-diode cells: 19 % of 1350). Edge power <= 15 % needs R >= 100. Gain and SLL improve up to R 150; beyond that the 514 diodes thin out. Add-a-beam needs R >= 150 for state-1 edge <= 7 % (R 100-125: 22-32 %).
d. **Worst-case beta_ON** over the tolerance box: 1.073.

## 4. Risks
1. The OFF state is on the steep side of the dispersion, 0.5-1 GHz below a stop-band edge, and the pad capacitance is not in the diode model: +-5 fF moves beta_OFF by +-0.023.
2. The via is a lumped mushroom inductance (L = mu0 h); if the true value is 0.7x, beta_ON drops by 0.03 and A0 rises by 56 %: the HFSS cell check decides.
3. M6 NOT_MET: the board model reproduces the leakage level but not the old HFSS far fields (they need a cell-beta error of +0.15) and has no rim diffraction.

## 5. Next decision
Approve an HFSS cell package (item 3: eigenmode OFF/ON at two phases along x, y, 45 deg with the real via, anti-pad, pad and bond ply, plus a radiation-top run for alpha) for D2, or first re-centre the cell for robustness (lower beta_OFF, accept smaller alpha) and put the pad capacitance, bypass and choke into the cell model.
