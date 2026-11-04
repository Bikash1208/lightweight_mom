# New direction: switched surface-wave launcher (method 1), two-probe fallback (method 2), polarisation switching (method 3) — MODEL_ONLY

Scripts: `sidequest_launcher_model.py` (network), `_scan.py`, `_gain.py`, `sidequest_twoprobe_model.py`; tests `tests/test_sidequest_launcher.py`; data `launcher_scan2q.json`, `launcher_board_gain.json`, `twoprobe_model.json`.
References (bibliographic-only, none opened here): Zhao/Pan/Shen, IEEE TAP 2026, doi 10.1109/TAP.2025.3621516 (user-verified; Yagi launcher + coded hologram); Wang et al., IEEE AWPL 2024 (multifeed beam-switchable holographic); Li/Tang/Xiao, IEEE Access 2019 (polarisation-reconfigurable holographic, polarisation-insensitive hologram). No DOIs for the last two (not supplied, not verified).

## 1. Element type — what our tools can and cannot model
- **MoM:** has no vias or free-standing pins (the thin-wire via extension is parked). A planar-only element (printed strip + diode) is what the MoM sheet model covers. Our own ON-strip row data (dog-bone, 1 row): |S11| = 0.73, |S21| = 0.64 for an infinite row,
  i.e. an isolated element scattering amplitude |t| ~ 0.5 (t ~ r beta p / 2). Ring of 12 such real scatterers: front-to-back **6.5 dB at |t| = 0.5, 8 dB at 0.9** (best radius): **a planar-only ring does not reach 10 dB**.
- **A post inside the substrate** (shorted, top at the surface): the repo's pin-launch model gives R_sw = **0.26 ohm** into the surface wave (vs 28 ohm for the 6 mm probe): with the 2.44 ohm diode it is an |t| < 0.1 scatterer. Same reason the mushroom gate failed.
- **Raised parasitic pins** (as the Yagi launcher of the reference): pin top 4-9 mm above the surface, diode at the foot (between the pin and the ground, on the back). Coupling from `PinLaunch`: R_sw / R_rad / X_est =
  top 4: 8.9 / 14.6 / -27 ohm; 5: 16 / 26 / +8; **6: 28 / 44.6 / +48**; 7: 49 / 77 / +97; 8: 89 / 137 / +163 ohm (resonance ~ top 4.7 mm). t = -R_sw/(R_rad + jX + Z_diode): max |t| = R_sw/R_rad = 0.63 (direct space radiation of the pin is the loss).
  Diode: ON 2.44 - 3.86j ohm (measured), OFF -536j ohm (28 fF): OFF pin |t| ~ 0.05 (transparent), ON pin strong. Network model: normalised impedance matrix Z~_mn = H0^(2)(beta rho_mn) (surface-wave mutual), self (R_rad + jX + Z_diode)/R_sw;
  **power balance (surface wave + space radiation + diode loss = accepted) closes to 1.0000**; bare probe reproduces the PinLaunch launch fraction 0.627; the probe is the HFSS-calibrated feed (R x 1.30, S11 -15 dB bare).

## 2. Launcher result (2 quadrants, centre probe; ON = pins in the back half, OFF = transparent; each pin is ON in exactly one of the two states, so diodes = N)
Scan N 4-12, top 3-9 mm, ring radius 4-35 mm (896 cases): **only 7 meet F/B >= 10 dB and S11 <= -10 dB**. Recommended **N = 12, top 6.0 mm (same as the probe), ring radius 10.5 mm (pin spacing 5.5 mm), 12 diodes**:
| quantity | value |
|---|---|
| surface-wave front/back (half-plane vs half-plane) | **10.6 dB** (93 % of the launched SW power into the chosen half: 0.585 of the accepted power; bare probe: 0.31 per half) |
| S11 | **-10.3 dB** (bare calibrated probe -15 dB) |
| launched SW / accepted | 0.636 (bare 0.627); space radiation 0.35, diode loss 0.014 |
| tolerance | diode R +-20 %, C +-20 %: F/B 10.1-10.2 dB, S11 -10.8 to -11.0; pin top -0.5 / +0.5 mm: 10.5 / 9.4 dB; ring radius -1 / +1 mm: 11.0 / 8.8 dB (**thin margin**) |
| N = 8, top 8 mm, rho 6 mm | F/B 10.8 dB, S11 -11.1, but only 0.43 of the power launched into SW (tall pins radiate: 0.57 into space) |
| 4 quadrants (ON = pins outside the chosen 90-deg window), N <= 16 | best F/B 10.0 dB at S11 -8.0 (N 16, top 8 mm, rho 12 mm): **no design meets both** |
On the 2-quadrant board (idealised holograms, alpha capped at the available 3.4 Np/m, SW weighted by the launcher's angular pattern): **gain 20.3 dBi** (D 24.6 dBi, conversion 0.58), **other-quadrant beam only 9.9 dB down**
(quadrant isolation ~ launcher F/B). Compare the gate approach (ideal 20.75 dBi, but the real gate leaks 3 dB).

## 3. Method 2 (fallback): one probe per quadrant, ideal SPDT (150 mm hexagon, probes at +-60.2 mm)
Gain **20.5 dBi** per state (D 25.2 dBi, SLL -17.5 dB, conversion 0.55 at the available alpha); **isolation between quadrants 22.9 dB** (the other half's hologram sees the wrong phase front);
**probe-probe coupling S21 = -30.4 dB** at 120 mm (SW mutual impedance, hologram attenuation included), bare-probe S11 -15.0 dB (neighbour terminated in 50 ohm, no effect at this distance).
Costs: 2 probes and an RF switch (SPDT insertion loss not included), no parasitic diodes.

## 4. Method 3 (polarisation switching with a polarisation-insensitive hologram + switched feed) — not worth it on this stack
1. RT5880 1.575 mm has **no TE surface wave** (TE1 cutoff c/(4 h sqrt(eps_r - 1)) = 43.4 GHz); every feed (probe, loop, slot) excites only the TM0 wave, whose local E is radial.
2. The radiated polarisation then comes from the cells' anisotropy relative to the radial direction, not from the feed; a switched feed can move the radial pattern but not select LP-H/V or RHCP/LHCP.
3. Polarisation switching would need a second TE-capable (thicker, ~4x) substrate or per-state tensor cells: outside the fixed stack. Recommendation: keep polarisation per quadrant (static tensor holograms).

## Recommendation
Method 1 with 12 raised pins (12 diodes) meets 10 dB / -10 dB only marginally (F/B 10.6 dB, tolerance +-1 mm ring radius = +-1.1 dB), gives 20.3 dBi but only ~10 dB quadrant isolation; method 2 gives 22.9 dB isolation and 20.5 dBi with two probes and a switch.
If the quadrant isolation target is 15 dB, method 1 does not reach it in this model. One small HFSS launcher test (probe + 12 pins on a bare 80 mm board, STEP + one page) is the next step **after approval**.
Model limits: surface-wave-only mutual coupling (no pin-pin space-wave coupling), thin-wire reactance estimate from PinLaunch, hologram scattering back to the ring not modelled, phase-centre shift of the ring ignored, free-standing pins are a fabrication question (out of scope), beta of the bare slab (1.018 k0) around the launcher.
