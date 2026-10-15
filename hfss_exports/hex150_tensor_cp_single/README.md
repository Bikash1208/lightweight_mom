# hex150_tensor_cp_single

Single 150 mm-circumradius hexagon, ONE monopole at the centre, rotating-slot tensor hologram for a **pure RHCP** pencil
beam at theta 30 deg, phi 225 deg, 10.54 GHz. Board RT/duroid 5880 (eps_r 2.2, tan d 0.0009), 1.575 mm. Phase 26 item D.
Everything below the "Predictions" heading is a MODEL prediction (SOLVER_SUPERCELL, written before any HFSS run); nothing here
is HFSS-validated. See `results/phase26/README.md` (sections B2.4b-f and C) for how it was derived and what it assumes.

## Files
| File | Content |
|---|---|
| `hfss_build_board.py` | AEDT (IronPython) script: builds substrate, ground, pin, lumped port, 1430 slotted circles (exact: circle minus a rectangle with exactly rotated corners, mm), radiation box, setup 10.54 GHz (+ 10.01-11.07 GHz interpolating sweep), far-field sphere theta 0-180 / phi -180..180, 1 deg |
| `layout.json`, `layout.csv`, `layout.dxf`, `layout.png` | Every cell: centre, circle D 6.0 mm, slot length x width, slot angle (long axis, from +x, deg). DXF layers PATCH, SLOT, OUTLINE |
| `board.step` | Thin (35 um) copper solids, NOT committed (35 MB); regenerate with `python scripts/research/phase26_export_single_hex.py` |
| `predictions.json` | The numbers below, machine readable |

Geometry: hexagon vertices (0, +/-150), (+/-129.904, +/-75) mm (vertex at +/-y, the orientation of port 3 of `four_port_hex150_gap0`), feed at (0, 0); lattice
6.4 mm centred on the feed cell (which, and everything within r < 6.5 mm, is empty); board 270 x 310 mm centred on the feed. Slot sizes in
use: 4.8x1.6 (869 cells), 4.6x1.6 (339), 4.8x2.0 (195), 4.4x1.2 (27); smallest feature: slot-end bridge 0.43 mm, gap between circles 0.4 mm.

## HFSS setup
- One lumped port (50 ohm) on the coax annulus r_pin .. r_out = 0.305 .. 2.05 mm at z = 0; radiation boundary on the padded air box (15 mm sides/bottom, 30 mm top).
- **Feed:** `mono_h` = 110 mil = 2.794 mm and `r_pin` = 12 mil = 0.305 mm, the HFSS variables of the `hex150_gap0` run. The real pin height above the
  ground is still to be confirmed (the old exporter default was 7.1 mm; the launcher-calibration package carries both readings). Edit `mono_h` after the build if it differs.
  Section A2 of the results README showed the launch fraction barely depends on the pin height (0.59 bare slab) — the pin's direct space wave (~41 %) is a property of the monopole.
- Export: antenna parameters table at 10.54 GHz and 10.01-11.07 GHz, realized-gain RHCP/LHCP and axial ratio over the sphere (1 deg CSV), the S11 sweep. For the
  surface-wave decay: Mag_E and complex Ez along rho = 10..140 mm at phi = 225 deg, 1 mm above the board (as in the calibration package).
- Run AEDT: Tools > Run Script > `hfss_build_board.py`; building 1430 circle/slot subtractions takes a few minutes.

## Predictions (MODEL_ONLY, 10.54 GHz; compare with the HFSS run)
| Quantity | Prediction | Brief target |
|---|---|---|
| Beam peak | theta 30 deg, phi 224-225 deg | 30 / 225 |
| Peak realized gain RHCP | **23.5 dBi** with the bare-slab launch fraction 0.592; 24.0-25.2 dBi if the loaded-slab launch is 0.66-0.88 | >= 24 dBi |
| Surface-wave power converted in the hexagon | 0.54 (rays 0.44-0.64) | >= 0.85 (unreachable on this board: B2.4f) |
| Axial ratio at the peak / worst over the -3 dB beam | 2.5 dB / 3.1 dB | <= 3 dB |
| XPD at the peak | 16.9 dB | >= 25 dB |
| RHCP / LHCP fraction of the aperture pattern | 0.975 / 0.025; RHCP within 8 deg: 0.93 of aperture power | LHCP <= 10 % |
| RHCP SLL | -22.8 dB | <= -15 dB |
| Power budget (of accepted) | main beam 0.30, aperture-radiated 0.32, unconverted at the edge 0.27, pin direct space wave 0.41 (assumed outside the beam; not patterned) | grazing + back <= 15 % |
| Frequencies 10.33 / 10.75 GHz | NOT_AVAILABLE_FROM_CURRENT_DATA | |

What would falsify the model: peak gain below ~19 dBi (the 1-bit board's HFSS value), peak not within 2 deg of (30, 225), RHCP fraction in the beam < 0.8,
or a decay of the surface wave along the 225 deg line much slower than exp(-alpha rho) with alpha ~ 2-4 Np/m (fit from the Ez cut).
The aperture part of the model is idealised (no feed scattering, perfect phases); the previous HFSS comparison says the leakage/launch parts under-predicted, so
a value above 23.5 dBi would not be a surprise and a value below 20 dBi would point at the cell/locality model.
