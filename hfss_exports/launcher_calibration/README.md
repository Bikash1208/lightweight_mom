# launcher_calibration — bare grounded board, no hologram (Phase 26)

Purpose: measure what the surface-wave launcher alone does, so the hologram's leakage and the launcher's
launch fraction are no longer entangled in one HFSS number. Same substrate (RT/duroid 5880, eps_r 2.2,
tan d 0.0009, 1.575 mm), same board outline as `four_port_hex150_gap0` (652.8 x 537.6 mm, x -326.4..326.4,
y -268.8..268.8), ONE feed at port 3's position (-64.952, -112.5) mm, **no patches**. 10.54 GHz setup,
+-5 {'kappa': np.float64(42.63957340998781), 'b': 1.0184589177682728, 'w': 2.400984883806226, 'ad': 0.0067262750972242194, 'table': '| Variant | launch fraction | R_rad [ohm] | X est [ohm] | est. |S11| [dB] | |E_z| 20 / 50 / 100 mm [V/m per sqrt(W accepted)] |\n|---|---|---|---|---|---|\n| A_current_pin | 0.594 | 1.7 | -115 | -0.1 | 515.6 / 326.1 / 230.5 |\n| B_matched_quarter_wave_pin | 0.617 | 26.0 | 8 | -9.6 | 525.6 / 332.3 / 234.9 |\n| C_tall_pin_top_hat | 0.821 | 43.8 | -161 | -1.3 | 606.4 / 383.4 / 271.0 |'}weep, far-field sphere theta 0..180 deg (1 deg), phi -180..180 deg (2 deg), as in the main export.

## Which pin is which

| Script | Pin | Why |
|---|---|---|
| `build_A_current_pin.py` | `mono_h` = 2.794 mm from the ground (110 mil), radius 12 mil | your current pin, read as PinLength measured from the ground plane (the `mono_h` convention of the exported scripts). **If your PinLength is measured above the substrate surface, set `mono_h = 4.369mm` in HFSS** (design variable) -- please tell me which it is. |
| `build_B_matched_quarter_wave_pin.py` | `mono_h` = 6.575 mm (5.0 mm above the surface), radius 12 mil | the model's best *launched share of accepted power* (launch fraction x match) |
| `build_C_tall_pin_top_hat.py` | `mono_h` = 11.575 mm with a 4 mm-radius, 0.1 mm top-hat | the model's best launch fraction (0.82) -- estimated unmatched, so this is a physics probe, not a design |

Run in AEDT: Tools > Run Script. Each script builds substrate, ground (with the coax opening), pin, one 50 ohm
lumped port, air box/radiation boundary, setup, sweep, far-field sphere, and a non-model polyline
`RadialLine225` 1 mm above the surface from rho = 10 to 140 mm at phi = 225 deg from the feed.

## What I need back (per script you run; A first, then B)

1. **Antenna parameters** table (port 1 excited, 1 W): radiated / accepted power, efficiencies, peak
   directivity, as in `hfss_results/hex150_gap0/port3_antenna_params_*.txt`.
2. **Far-field CSV**: realized gain RHCP and LHCP (or Etheta/Ephi, whichever is easier) over the full sphere,
   at 10.5 or 10.54 GHz, same format as `port3_realized_gain_RHCP_LHCP_10p5GHz.csv`.
3. **E-field along `RadialLine225`** at the same frequency: Fields calculator > Input > Quantity > E; Complex
   > Mag (for Mag_E) and separately the **complex Ez** (Real and Imag, or Mag and Phase); Input > Geometry >
   Line > `RadialLine225`; Value (eval) at 141+ points; Results > Create Fields Report > Rectangular Plot with
   X = Distance, then File > Export to a CSV. Export both Mag_E and the complex Ez (Re, Im) vs distance.
4. S11 vs frequency (the touchstone `.s1p` or a CSV of dB(S11) and the real/imag input impedance, 8-12 GHz) --
   this is what checks the model's thin-wire match estimate.

## How I will use them

* The surface wave above a grounded slab: E_z falls off as exp(-kappa z) with kappa = 42.6 Np/m
  (bound mode beta/k0 = 1.0185), so |E_z| at 1 mm gives the field at the surface. The cylindrical wave then
  gives the power through a circle of radius rho: P_sw(rho) = 2 pi rho w |H0|^2 with w = 2.401 ohm m and
  |H0| = |E_z0| omega eps0 / beta. Fitting |E_z|(rho) to rho^-1/2 exp(-alpha_d rho) separates the launched power
  from the (small, predicted alpha_d = 0.007 Np/m) dielectric loss, and `launch fraction = P_sw / accepted power`.
* The far-field integral gives the directly radiated power, a second, independent measurement of the same split.
* The far-field *pattern* of A (donut-like, azimuth-uniform, peak towards the horizon) vs the patterned
  board's HFSS pattern shows how much of the hex150 board's grazing power is the launcher itself.

## Model predictions to compare (MODEL_ONLY; thin-wire match estimate, good to tens of ohms)

| Variant | launch fraction | R_rad [ohm] | X est [ohm] | est. |S11| [dB] | |E_z| 20 / 50 / 100 mm [V/m per sqrt(W accepted)] |
|---|---|---|---|---|---|
| A_current_pin | 0.594 | 1.7 | -115 | -0.1 | 515.6 / 326.1 / 230.5 |
| B_matched_quarter_wave_pin | 0.617 | 26.0 | 8 | -9.6 | 525.6 / 332.3 / 234.9 |
| C_tall_pin_top_hat | 0.821 | 43.8 | -161 | -1.3 | 606.4 / 383.4 / 271.0 |

|E_z| is a peak amplitude per sqrt(W) of ACCEPTED power (HFSS antenna-parameter table: Accepted Power; scale by sqrt(accepted [W])), 1 mm above the
surface, rho = 10, 20, 50, 100, 140 mm, including the predicted dielectric loss and cylindrical spreading.
A near-field radial line passes close to the pin at rho = 10 mm: the 10 mm value includes non-TM0 near field.
Model launch fraction is the TM0 share of the power radiated by the pin (`src/dispersion/pin_launch.py`).
