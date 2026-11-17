# Item 0: ideal aperture estimate (MODEL_ONLY)

At 5.8 GHz, a circular, uniformly illuminated aperture was divided into two half-discs or four quadrants. `D150` means diameter 150 mm; `R150` means radius 150 mm. Sector phase ramps point at theta 30 degrees and phi 0/180 (two sectors) or phi 0/90/180/270 (four sectors). Hands alternate RHCP/LHCP. A 5 mm aperture grid and 2 degree upper-hemisphere angular grid were used. Each sector is assigned ideal circular polarization and a cos(theta) field element factor. Directivity uses total radiated power of both hands. AR and XPD include the opposite-hand field from other sectors. SLL is the highest **local maximum** outside approximate first-null exclusion cones centred on all intended same-hand beams. The cone radius uses 1.22 lambda / equal-area circular diameter, so this SLL definition is approximate for sectors and overlapping beams. The earlier global-maximum-outside-cone method selected a main-lobe shoulder; its SLL values are superseded here.

| Aperture | Partition | Beam directivity (dBi) | Worst SLL (dB) | Worst AR (dB) | Minimum XPD (dB) | Meets all targets? |
|---|---:|---:|---:|---:|---:|---|
| D150 | 2 sectors | 13.84 | -15.04 | 2.77 | 16.02 | Yes under ideal assumptions |
| D150 | 4 sectors | 8.98–9.91 | undefined | 5.93 | 9.66 | No: AR, XPD; lobes overlap |
| R150 | 2 sectors | 19.12 | -13.62 | 1.49 | 21.36 | No: SLL |
| R150 | 4 sectors | 13.10–13.56 | -9.03 | 4.94 | 11.14 | No: AR, XPD, SLL |

Targets are AR <= 3 dB, XPD >= 15 dB and SLL <= -15 dB. For the R150 two-sector case, the worst true local maximum is near boresight: RHCP at theta 2°, phi 180° and its LHCP counterpart at theta 2°, phi 0°, both -13.62 dB relative to their intended beam. Values per beam and actual peak angles are in [aperture.csv](aperture.csv), with full parameters in [aperture.json](aperture.json). No feed, real cell response, coupling, taper, dielectric/conductor loss, absorber ring or HFSS validation is included. This estimate does not justify building a board or claiming a switchable antenna.

## HANDOFF

- STATUS: Item 0 corrected, MODEL_ONLY; SLL uses actual local maxima.
- Cell JSON: pending 5.8 GHz baseline screen; proposed baseline p 11.6 mm, r 5.6 mm, 3 teeth, 1.2 mm tooth width, 0.4 mm slot, 3.0 mm run along x.
- Beta/k0 along/across: pending MoM screen.
- Alpha: pending MoM screen; aperture estimate has no leakage parameter.
- M_t: pending MoM screen.
- Loss ratio: pending MoM screen.
- HFSS or MODEL_ONLY: MODEL_ONLY; no HFSS check.
- Next request: check whether aperture taper/sector boundaries can meet SLL, AR and XPD after the cell response is measured.
