# phase28_launcher -- Yagi-type side launcher on RT/duroid 5880 3.175 mm, 10.0 GHz (driven pin + reflector post + 2 director posts)

**State: PREDICTION_NOT_HFSS_VALIDATED (MODEL_ONLY, ESTIMATE_COUPLED_MONOPOLES).** The launcher after [78] (a monopole turned into a Yagi-Uda with a reflector and directors, all grounded posts) for the first board. Target: half-power width of the launched surface wave about 78 deg, S11 <= -10 dB. This is a launcher check, not a hologram: the surface is the uniform OFF-state loop cell, so the launched wave sees the real surface.

## Geometry (from the model; top heights are above the substrate surface, post radius = pin radius 0.635 mm)
| element | x [mm] | top above surface [mm] | top from ground [mm] |
|---|---|---|---|
| driven pin (coax lumped port at the base) | 0 | 3.73 | 6.90 |
| reflector post | -9.03 | 8.36 | 11.54 |
| director 1 | 13.14 | 5.52 | 8.70 |
| director 2 | 27.10 | 5.11 | 8.29 |

Surface: hexagon of circumradius 90 mm (vertex at +-y), loop cells (period 9.0 mm, outer 6.6 mm, strip 1.2 mm, OFF diode = 28.2 fF only) outside a keep-out of radius 38 mm: 200 cells. The wave goes out along +x (directors) and is damped toward -x (reflector). `LOADED_SURFACE = False` builds the bare slab (run this first).

## Files
| File | Content |
|---|---|
| `build_launcher.py` | AEDT script (variables at the top), setup 10 GHz + interpolating sweep 9.5 - 10.5 GHz, a non-model polyline `AzimuthCircle60` (R = 60 mm, 1 mm above the surface) for the azimuth pattern |
| `mom_reference_pattern.csv` | predicted TM0 power pattern versus azimuth, normalised, for the two feed calibrations |
| `mom_reference_s11.csv` | predicted S11, half-power width, front-to-back, launch fraction at 9.5 - 10.5 GHz |

## Predictions first (model, 10.0 GHz; A = thin wire r 0.635 mm, B = Phase 26 HFSS calibration; they bracket the feed, neither is validated at 3.175 mm)
| quantity | A | B |
|---|---|---|
| input impedance [ohm] | 29.3 + j22.2 | 39.3 + j9.6 |
| S11 [dB] | -8.7 | -15.9 |
| half-power width of the TM0 wave [deg] | 77 | 79 |
| front-to-back [dB] | 8.3 | 5.4 |
| share of the surface-wave power in |psi| < 90 deg | 0.70 | 0.63 |
| launch fraction (TM0 share of the accepted power) | 0.83 | 0.82 |

**Model limits (stated, not hidden).** The pins are coupled through the TM0 cylindrical wave only (Z_mn = sqrt(R_sw,m R_sw,n) H0^(2)(beta d), the pins are 9 - 14 mm = 0.3 - 0.5 lambda apart); the near-field and space-wave mutual impedances are neglected. The posts' own radiation is in their self terms. The loaded surface and the keep-out edge are not in the model (bare slab, beta/k0 = 1.07). The feed calibrations differ by 7 dB in S11, so S11 is an estimate.

## Run order
1. Run `build_launcher.py` with `LOADED_SURFACE = False` (bare slab, 4 cylinders, 1 lumped port); read the messages. Analyze. Export S11 (port 1, 9.5 - 10.5 GHz).
2. Export the field on `AzimuthCircle60`: Fields Calculator, E (complex magnitude) and the complex Ez, versus the polyline distance (the polyline starts at azimuth -180 deg and runs counter-clockwise in 5 deg steps; 73 points).
3. Set `LOADED_SURFACE = True`, rebuild, Analyze, export the same. The loaded surface adds the wave's refraction at the keep-out edge: compare the half-power width with the bare run.

## Pass criteria (ASSUMED, stated before the run)
* S11 <= -10 dB at 10.0 GHz (model range -9 to -16 dB: a miss of up to 4 dB is within the calibration uncertainty; > -6 dB means re-tune the driven-pin height `TD`).
* Half-power width of |Ez| (field, so -3 dB of power at the 0.707 level) on the 60 mm circle within +-15 deg of 78 deg; peak azimuth within +-10 deg of 0 (toward +x).
* Front-to-back of the power pattern >= 5 dB (model 5 - 8 dB). **News:** a width above 110 deg or a second lobe within 6 dB of the main lobe.

## Honest limits
The pattern is the TM0 part only; the space-wave radiation of the pins (the launch fraction is 0.8 in the model) is not in the azimuth pattern; the keep-out edge and the loaded surface scatter; zero-thickness metal; post and pin radius equal.
