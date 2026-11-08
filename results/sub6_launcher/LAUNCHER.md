# SUB6 centre-probe launcher, 5.8 GHz, RT5880 3.175 mm (stage: sub6-launcher)

**MODEL_ONLY / NOT_HFSS_VALIDATED.** Code `scripts/research/sub6_launcher.py` (uses `src/dispersion/pin_launch.py` unmodified; reference: `phase26_launcher_sweep.py`); data `sweep_coarse.csv` (1080 pts), `sweep_refined.csv` (274), `extend_coarse.csv` (900), `extend_refined.csv`, `launcher_summary.json`.

## Method and its limits
* Launch fraction L = P_sw/(P_sw + P_space) from `PinLaunch.split()` (thin filament from the ground to z_top above the surface, optional top-hat of radius a; bare slab, beta/k0 = 1.0228). Independent check `power_balance()` (lossy numerical integral) agrees with the split to 0.1-0.8 %.
* **Impedance estimate:** Z(f) = R_rad(f) + j X(f); R_rad = 2 (P_sw + P_space) is the radiated-power resistance (no dielectric, conductor or coax loss); X is the thin-wire transmission-line estimate -Z_w cot(Phi) of `pin_launch.py` (documented there as good to tens of ohms); no SMA flange, coax aperture or PTFE bead. Gamma = (Z - 50)/(Z + 50) at 5.725, 5.7625, 5.8, 5.8375, 5.875 GHz; the match criterion is the worst of the five. Pin length from the ground = 3.175 mm + z_top.
* Sweep: z_top -2.9..+12 mm (0.25 step, refined to 0.05), hat radius 0..6 mm (cap 6.0 mm: the first-ring cell metal starts at 11.6 - 5.44 = 6.2 mm from the feed; the lattice already keeps rho < 10 mm clear, and round-2 keeps switched cells out of rings 1-2), pin radius 0.25 / 0.4 / 0.635 mm (SMA centre pins), then an extension with 0.8-1.6 mm (non-standard thick rods, flagged) and tops up to 16 mm.
* L is for a bare slab around the feed; the real board has the cell layer there (beta 1.05-1.2), so L is a model number.

## Results
| pin | z_top (mm) | pin from ground (mm) | hat (mm) | radius (mm) | L (5.8) | L min band | R / X at 5.8 (ohm) | S11 5.8 / worst band (dB) | accepted share L(1-G^2) |
|---|---|---|---|---|---|---|---|---|---|
| used so far (Phase 26 pin scaled by h) | -0.72 | 2.46 | 0 | 0.615 | 0.6263 | 0.622 | 0.35 / -136.5 | -0.01 / -0.01 | 0.002 |
| best with a standard SMA pin (0.635 mm) | 8.95 | 12.1 | 0 | 0.635 | 0.656 | 0.650 | 27.3 / +14.2 | -9.4 / -9.3 (**misses -10 dB**) | 0.580 |
| best meeting -10 dB (non-standard 3.2 mm rod) | 10.00 | 13.2 | 0 | 1.6 | 0.6615 | 0.656 | 37.6 / +23.9 | -10.6 / -10.3 (meets, on the estimate) | 0.604 |
| max L at any match (unmatched) | 11.85 | 15.0 | 5.0 | 0.25 | 0.714 | 0.706 | 4861 / -1991 | -0.15 | 0.025 |

* **No standard-pin geometry reaches -10 dB** over 5.725-5.875 GHz in the model (closest -9.3 dB). A -10 dB match appears only with a 3.2 mm diameter rod (16 geometries at the edge of the extension grid; thicker is still better, so this is a boundary, not an optimum). The X of +24 ohm at that point is inside the estimate's own uncertainty: the result is **marginal, not a demonstrated match**.
* **Hat does not help:** best S11 worst-band by hat radius (r 0.635): 0 mm -9.2 dB, 1 mm -8.2, 2 mm -6.4, 3 mm -5.0, 4 mm -4.1, 5 mm -3.5, 6 mm -3.1. A larger hat raises L (to 0.71 at 3-5 mm) but R_rad collapses or the pin must be very tall; none is matched. Best matched designs therefore have no hat, which keeps the first cell ring clear.
* L varies only 0.63-0.71 over the whole sweep: the launch fraction is a weak function of the geometry; the geometry mostly sets the match.
* **Board gain change vs 0.6263:** the board model takes L as a scalar, so G changes by 10 log10(L_new / 0.6263) = **+0.20 dB** (standard pin, 0.656) or **+0.24 dB** (thick rod, 0.6615), referred to accepted power. Verified on the round-3 point 200 mm / N 135 / a: OFF G_A 9.82 -> 10.02 / 10.05 dBi, ON G_A 14.7 -> 14.9 / 14.93 dBi; gate unchanged (MISSES_SLL).
* **Important side finding:** the pin used so far (inside the substrate, 2.46 mm from the ground) has R_rad = 0.35 ohm and X = -137 ohm: against 50 ohm that is |S11| = -0.01 dB, i.e. in this model only about 0.2 % of the source power would be accepted. All board gains so far (rounds 1-3) are referred to accepted power and do not include this mismatch. Including it, the old pin would cost about 27 dB; the matched pin recovers an accepted share of 0.58-0.60 (this uses the same crude estimate, so treat the order of magnitude, not the number).

## Flags
MODEL_ONLY; X estimate good to tens of ohms; no coax/SMA detail; pin heights of 9-10 mm above the board are long probes (SMA centre pins protrude about 3-5 mm); the thick-rod match needs a 3.2 mm conductor (or an equivalent wide strap/cone), not a standard SMA pin; an HFSS calibration (as `hfss_exports/launcher_calibration` at 10.5 GHz) is needed before building. Not done: a matching section, a cone/strap probe, the effect of the cell layer on L.
