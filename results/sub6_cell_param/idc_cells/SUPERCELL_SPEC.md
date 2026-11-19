# Supercell leakage spec (results only; no tool code)

Goal: measure the radiative leakage of rotated interdigital (IDC) cells for the 5.8 GHz RHCP/LHCP holographic board. Unit-cell eigenmode tuning is finished (see README.md); near-resonance cells have non-unique beta pairs, so leakage has to be measured on a rotated supercell.

## Stack (all unit-cell rows)
- RT5880: eps_r 2.2, tan_d 0.0009, h = 3.175 mm; copper 0.035 mm, modelled as PEC; ground plane under the substrate; a 0.254 mm RT5880 bias layer below the ground (electrically irrelevant).
- Unit-cell rows: period 11.6 mm; PEC boxes ("C40B"): air_top = +40.035 mm (lid above the substrate top at z = 0, copper z = 0..0.035), air_bottom = -43.499 mm; eigenmode search from 3.0 GHz; master/slave lattice pairs.
- Direction convention: `across` = phase on x, `along` = phase on y. beta/k0 = (phase/80.8) * (5.8/f_GHz) at 11.6 mm pitch. Stored beta_over_k0 in the pushed rows uses 80.8 (ignore it for the 14.5 mm rows).

## Cell geometry (top view, copper = disc minus slot; all slot cuts 0.4 mm wide unless W given; axis of the meander = x)
Common IDC (N = 2 teeth, leg pitch T = 1.2 mm, x0 = -N*T):
- legs at x = x0 + i*T (i = 0..4), each spanning y = -(D/2 + w/2) .. +(D/2 + w/2);
- horizontal bars i = 0..3 at y = +D/2 (i even) or -D/2 (i odd), spanning x = x0 + i*T - w/2 .. x0 + (i+1)*T + w/2;
- end runs: left at y = -D/2 from the rim to x0 + w/2; right at y = +D/2 from x_e - w/2 to the rim (x_e = x0 + 2*N*T = 2.4).
- bridge b: each end run stops at |x| = rimx - b, rimx = sqrt(r^2 - (D/2 + w/2)^2); no bridge = run to the rim (two floating halves).
| Name | r | w | D | b | copper area (build log) |
|---|---|---|---|---|---|
| MN_R5p6_W0p4_N2_T1p2_D3p0_E_B2p2 (B2p2) | 5.6 | 0.4 | 3.0 | 2.2 | 90.0 mm2 |
| MN_R5p6_W0p4_N2_T1p2_D3p6_E_B2p2 (D3p6 bridged) | 5.6 | 0.4 | 3.6 | 2.2 | 88.9 mm2 |
| MN_R5p6_W0p4_N2_T1p2_D3p6_E (D3p6 open twin, OFF bound) | 5.6 | 0.4 | 3.6 | none | 87.08 mm2 (two halves 43.5 + 43.5, 0.4 mm apart) |
Anchor (not interdigital, reference): ST_R5p45_L8p72_W2p91 = disc r 5.45 with one centred rectangular slot 8.72 (y) x 2.91 (x) mm, slot along y; copper ~67.9 mm2 (computed).
Clearance: all discs are round (rotation safe), gap to the neighbouring cell >= 0.4 mm (r 5.6) or 0.7 mm (r 5.45); minimum feature 0.3 mm.

## Calibration rows (pushed JSON in this folder; reproduce lidded, same stack, within 0.5 % in f before any supercell run)
- B2p2: `MN_R5p6_W0p4_N2_T1p2_D3p0_E_B2p2_OFF_across_85.json` (SW 5.685 GHz) and `..._along_88/95` (y).
- D3p6 bridged: `MN_R5p6_W0p4_N2_T1p2_D3p6_E_B2p2_OFF_across_80/85/95.json`, `..._along_88/95/102.json`.
- D3p6 open twin: `MN_R5p6_W0p4_N2_T1p2_D3p6_E_OFF_across_85/95.json`, `..._along_85/95.json`.
- Anchor: `ST_R5p45_L8p72_W2p91_OFF_across_85/95.json`, `..._along_85/95.json`.
Each row lists modes (f, Q). Q > 1e5 modes are box/lid artefacts; ignore them.

## Supercell
- 8 cells along x, rotations psi_n = 0, 22.5, ..., 157.5 deg (in-plane, about each disc centre), pitch 11.6 mm: supercell 92.8 x 11.6 mm. y periodic at 0 deg; x Bloch phase phi_s is a parameter.
- Rotation of the polarisation tensor by psi gives a 2*psi = 45 deg step per cell, so K = 2*pi/92.8 mm, K/k0 = 0.557 (Phase 26: 0.5555). The -1 harmonic (beta - K)/k0 ~ 0.52 radiates at ~31.5 deg (beam design 30 deg).
- Mode identification: unit-cell branch f(phase) with beta(f) from the rows; in the supercell a branch appears where beta(f)*92.8 mm = phi_s + 360*m. Use the unit-cell tables to predict mode frequencies near 5.8 GHz; match by Q and field pattern.
- Top boundary: the lidded unit-cell Q includes only dielectric loss (PEC copper). The supercell needs an open top (radiation boundary or PML, air >= lambda/4 ~ 13 mm, bound-wave field decay length 14-29 mm, so use ~40 mm unless a PML spurious-loss check (unrotated lidless unit cell, alpha_spur <= 0.1 Np/m) passes at a lower height), or a driven strip (TM0 wave port, 1-2 rotation periods, absorber far end, alpha from decay).
- Extraction: alpha = pi*f/(Q_rad*v_g), v_g = d(2*pi*f)/d(beta) = 360*p*(df/dphase) with p = 0.0116 m (unit cell) or from the supercell phase scan; 1/Q_rad = 1/Q - 1/Q_lidded (Q_lidded from the unit-cell row of the matched branch). Dielectric baseline 0.01-0.09 Np/m vs a 1.2 gate.

## Expectations (MODEL_ONLY maps from the other agent; calibration, not truth)
- Anchor 1.16-1.48 Np/m (y pairing ambiguous); B2p2 ~0.42 Np/m; D3p6 bridged 1.09 Np/m (x mode 1.081) to 2.38 Np/m (x mode ~1.063); OFF open gap: 0.015 or ~4.6 Np/m depending on y pairing (1.139 vs 1.038).
- Gates: ON alpha >= 1.2 Np/m and alpha_ON/alpha_OFF >= 12 (~10 dB, field-ratio r <= 0.3).

## Not in this round
- The diode-sheet cell: replacing a 2.2 mm copper bridge by a ~0.4 mm sheet changes the inductance that sets the resonance, so the sheet position and series L need design first (and a land-pattern check against 0.3 mm). The open-gap twin is the OFF bound for now. Crop-testbed diode OFF (-10 V) sheet: 21.6 - j297.5 ohm; tool OFF sheet 22.23 - j335.4 ohm.
- Which mode a pin launch excites: the supercell gives each mode's leakage only. If one D3p6 x mode fails, flag it for a driven-strip run.

## Corrections after the first run (other agent, AEDT 2025.2; local edit, not yet pushed)
- Eigenmode accepts a PML top (finite Q), rejects a radiation top. Calibration passed: anchor 5.56231 vs 5.55317 GHz (0.165 %), B2p2 5.69323 and D3p6 bridged 5.68256 GHz matched. The open twin returned 0 mm2 copper (readback/Boolean issue in their build; 87.08 mm2 in 2025.1): it blocks only the OFF bound.
- Air height >= lambda/4 was too small; use the PML spurious-loss check and Q_base = mean of unrotated lidless-PML across/along Q at the matched per-cell phase.
- Supercell Bloch phases 15/30/45 deg (avoid 0 degeneracy): per-cell phases (720 +/- phi_s)/8 = 91.9/88.1, 93.8/86.3, 95.6/84.4 deg.

## Update 2026-11-18
5.8 GHz is not strict (user). Gates evaluated at one common f_op shared by ON/OFF and both beams; hologram period re-derived at f_op.
Beam angle theta free in {25,30,35,40} deg, both beams equal, chosen jointly with f_op at the end; alpha(theta) from the first-order chain.
