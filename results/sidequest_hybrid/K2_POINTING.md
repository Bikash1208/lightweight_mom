# Beam pointing over tolerance, probe without series C, bound-mode flag - MODEL_ONLY (2026-11-06); no package built
Bare launch 0.71 for every gain; gains calibrated to the HFSS state-1 point (offset -1.78 dB); <= 514 diodes, one line. r 2.25 is the tentative nominal pending the HFSS cell check.

## 0. Caution first: the MoM cell numbers are pixel-limited at these radii (`sidequest_mesh_check.py`, `k2_mesh_check.json`)
RO3006 1.27, eps_r 6.15, pair (OFF / ON), mesh 64 (0.1 mm pixel) against mesh 96 (0.067 mm pixel); Ybar [mS], A0 [Np/m], beta/k0:
| r [mm] | 2.19 | 2.22 | 2.25 | 2.28 | 2.31 |
|---|---|---|---|---|---|
| mesh 64 | 6.82, 4.6, 1.199 | 7.18, 6.0, 1.220 | 7.36, 6.7, 1.232 | 8.96, 19.3, 1.346 | 9.22, 19.3, 1.366 |
| mesh 96 | 6.71, 3.9, 1.192 | 6.96, 4.6, 1.207 | 7.94, 11.7, 1.271 | 8.24, 13.1, 1.293 | 9.77, 22.5, 1.407 |
The jump of beta/k0 sits between 2.25 and 2.28 at mesh 64 but between 2.22-2.25 and 2.28-2.31 at mesh 96: the disc edge is staircased, so radius steps of 0.03 mm are below the resolution. The trend (A0 about x5 and beta/k0 +0.2 over 0.12 mm, i.e. d beta/dr about 1.7 per mm) is real, the exact position of the steps is not, and A0 at one radius is uncertain by about x2 (beta/k0 by +-0.06). A converged curve needs mesh >= 128 (hours per solve): the HFSS parametric sweep r 2.15 / 2.20 / 2.25 / 2.30 (OFF/ON eigenmode at 90 deg) would give d f/d r and replace the MoM slope. All tolerance numbers below inherit this limit.

## 1. Pointing and design-direction gain over the cases (`sidequest_pointing.py`, `k2_pointing_*.json`, `k2_pointing_table.json`)
Method: the hologram (cell type per position) is designed ONCE for the nominal (beta_nom, nominal library); each case is the same board with that case's levels, A0, beta/k0 and background loss. Peak = maximum on the 1 deg x 2 deg grid; gain at the design directions (25, 0) / (25, 180).
**Nominal r 2.25 / eps_r 6.15 (beta/k0 1.232): 16.0 / 16.7 dBi.** Columns: beta/k0 | peak (theta, phi) state 1 / state 2 | pointing error [deg] | gain at the design directions [dBi]
| case | beta/k0 | peak s1 / s2 | error | gain at design dir s1 / s2 |
|---|---|---|---|---|
| r 2.20, 6.00 | 1.175 | (25,0) / (25,180) | 0 / 0 | 13.0 / 13.9 |
| r 2.22, 6.00 / 6.15 / 6.30 | 1.19 / 1.22 / 1.26 | (25,0) / (25,180) | 0 | 14.2 / 15.0; 15.7 / 16.5; 16.5 / 17.1 |
| r 2.25, 6.00 / 6.30 | 1.20 / 1.27 | (25,0) / (25,180) | 0 | 14.7 / 15.5; 16.4 / 17.0 |
| r 2.28, 6.00 | 1.293 | (25,0) / (25,180) | 0 | 16.5 / 13.8 |
| r 2.28, 6.15 | 1.346 | (25,0) / (25,180) | 0 | 14.0 / 11.3 |
| **r 2.28, 6.30** | **1.406** | **(29,14) / (24,164)** | **7.5 / 6.7** | **9.7 / 7.2** |
| r 2.30, 6.15 | 1.359 | (25,0) / (25,180) | 0 | 13.6 / 12.2 |
| **r 2.30, 6.30** | **1.420** | **(29,16) / (25,162)** | **8.3 / 7.6** | **8.4 / 6.9** |
| r 2.35, 6.30 | 1.585 | (21,302) / (22,238) | 22 / 22 | 1.2 / -1.2 |
**The linear estimate sin(theta) = beta/k0 - 0.807 (1.36 -> 33 deg, 1.42 -> 38 deg) is not what the centre-fed model gives.** For a radial surface wave the beta mismatch is a phase error -(beta - beta_nom) rho, an EVEN function of position: the two halves of the aperture tilt the beam in opposite directions, so the peak stays at (25, 0) / (25, 180) (resolution 1 deg) and the penalty is gain at the design direction (-2 dB at beta/k0 1.35, -6 .. -9 dB at 1.40-1.42, where the beam breaks up and the peak moves to (29, +-15) = 7.5-8.3 deg off). The linear rule applies to a plane-wave (side-fed) launch.
**Pointing error > 3 deg occurs only at beta/k0 >= 1.40 (r 2.28 / eps_r 6.30, r 2.30 / eps_r 6.30)**; inside the box r 2.22-2.28 x eps_r +-0.15 around r 2.25: maximum error 7.5 deg, worst gain at the design direction 7.2 dBi (-8.8 dB against 16.0 / 16.7). The loss is on the high-beta side only.
**Fix, trade-off: lower the nominal to r 2.22** (design beta/k0 1.220, 15.7 / 16.6 dBi nominal = -0.3 / -0.1 dB against r 2.25). Box r 2.20-2.25 x eps_r +-0.15 (beta/k0 1.175-1.27): pointing error 0 deg (all peaks at the design directions), worst gain at the design direction 13.3 dBi (r 2.20, eps_r 6.00) = **-2.4 dB** against the nominal (the low-A0 side). Nominal r 2.20 (beta 1.201) is similar (15.2 / 16.1 nominal, same pattern, worst -2.4 dB on the low side, -4.6 .. -6.8 dB if r 2.28 / eps_r 6.15-6.30 is reached). The r 2.22 nominal is therefore the robust choice in the model; it moves the corner of the box away from the A0 shelf (beta/k0 1.35+). Gain is almost flat in the nominal radius (15.7 against 16.0 dBi), so the lower nominal costs 0.3 dB. Subject to the HFSS recentring (the cell offset shifts every radius above) and the pixel limit of section 0.

## 2. Probe without the series capacitance (thin-pin FD, `sidequest_probe_disc.py`, `k2_probe_disc.json`, bare-launch geometry)
* The FD is **not valid for pins shorter than about 3 mm above the surface** (it returns negative R: the thin-pin current / flux balance is not resolved; the spectral model gives R 4 ohm at t 1 mm). A short pin with a printed top disc therefore cannot be evaluated in this model for t < 3 mm; the disc runs shown there are discarded.
* Valid range (t >= 3 mm, pin through a PEC top disc of radius R_d): the disc moves X **up**, by +1 ohm (R_d 1.5) / +3.3 (2.5) / +8.6 (3.5 mm) at t 6 and leaves R unchanged (35 -> 37 ohm): the wrong direction for the inductive X of the 6-7 mm pin. At t = 5 mm: R 24.3 -> 25.5 ohm, X +4.7 -> +14.8 ohm (R_d 0 -> 3.5). At t = 3-4 mm the disc raises a capacitive X towards zero (t 4: -18 -> -6 ohm for R_d 3.5) but R stays 16 ohm (S11 -5.9 dB).
* Result: **|X| < 14 ohm AND S11 <= -15 dB (R >= 35 ohm) is not reachable by pin length plus a printed disc in this model**: X crosses zero at t about 4.8-5.0 mm where R is only 23-25 ohm (S11 about -9 dB); R reaches 35 ohm only at t 6 mm, where X is +29 .. +38 ohm. The two needs (R from pin length, X from the series element) are coupled by the pin physics.
* Starting point for the HFSS sweep (existing parametric SMA model): **PinLength (above the surface) 5.0 mm with no disc** (model R 24 ohm, X +5 ohm, S11 -9 dB), tuning direction: **each +1 mm of pin length gives about +11 ohm in R and +25 ohm in X**; to raise R without X, enlarge the pin diameter (a 1.0 mm: X -6 ohm at t 6, R 34 ohm) and increase the loading near the pin (the pad: rho_in 2.6 mm gives R 43 ohm, X +9 ohm at t 5, S11 -18 dB; UNVALIDATED). The baseline check (model S11 -12.7 dB against HFSS -15 dB at t 6 mm) says the X offset is a few ohm, so the sweep should cover t 4.5-6.5 mm and the pin diameter 1.27-2.0 mm.
* A capacitively coupled pin tip (pin ending g below the top layer, C about eps0 eps_r pi a^2 / g: 0.35 pF at g about 0.2 mm in RO3006, ignoring fringing) would give the series C without an extra part, but it is an idea only (not modelled here).

## 3. Flagged risk: bound-mode and model validity for beta/k0 > 1.3 (NOT resolved)
* The cell Y is extracted at NORMAL incidence (k_t = 0) but the surface wave travels at beta = 1.3-1.4 k0: 1.83-1.97 rad = 105-113 deg per 6.4 mm cell. The HFSS eigenmode check is at 90 deg; the RT5880 pipeline used the Bloch tensor at 90 deg. Spatial dispersion of the cell is NOT included in the A0 / beta used for the boards; it can move beta/k0 and A0 as much as the tolerances. The final library should use the Bloch tensor at beta (as in `C.bloch_loaded`).
* First-order leakage needs alpha << beta: at A0 20 Np/m, alpha / beta = 20 / 300 = 7 %; the SW is down to e^-3 in amplitude across 150 mm, so the aperture efficiency of the outer cells is low (engine: SW at the edge 2-9 %).
* Brillouin zone: beta p / pi = 0.63 at 1.42 (below 1.0), modulation period 29 mm = 4.5 cells (above Nyquist).
* The engine phase constant is one number per case (mean loading of the pair); the real board has a position-dependent mean loading (statics). A fixed-design / varied-beta check is in section 1.
