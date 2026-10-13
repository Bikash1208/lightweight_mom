# four_port_hex150_opt_all

4-port multi-port HMSA board, 10.54 GHz. Substrate eps_r 2.2, tan d 0.0009, 1.575 mm thick, PEC ground. Cells 6.4 mm. The board is 652.8 x 537.6 mm (x -326.4..326.4, y -268.8..268.8), one monopole per sub-aperture.
2867 patches in total. Every coordinate is exact (rounded to 0.001 mm; the largest rounding applied was 0.0004736 mm), and each port's cells sit at that port's own positions, not on a shared lattice.

## Ports and feeds

| Port | Feed (x, y) mm | Sub-aperture x mm | Sub-aperture y mm | Cells | Patches |
|---|---|---|---|---|---|
| 1 | (64.952, 112.5) | -70.4 .. 198.4 | -38.4 .. 268.8 | 1414 | 697 |
| 2 | (-194.856, 112.5) | -326.4 .. -64 | -38.4 .. 268.8 | 1435 | 714 |
| 3 | (-64.952, -112.5) | -198.4 .. 70.4 | -268.8 .. 38.4 | 1414 | 715 |
| 4 | (194.856, -112.5) | 64 .. 326.4 | -268.8 .. 38.4 | 1435 | 741 |

## Sub-aperture outlines

Polygon sub-apertures (cells whose centre is inside an outline carry the design; every other cell is empty, state 0). Outlines are on layer PORTS of `layout.dxf`, dashed in `layout.png`. The substrate and ground are the bounding rectangle of the cluster.

| Port | Vertices (x, y) mm |
|---|---|
| 1 | (64.952, 262.5); (-64.952, 187.5); (-64.952, 37.5); (64.952, -37.5); (194.856, 37.5); (194.856, 187.5) |
| 2 | (-194.856, 262.5); (-324.76, 187.5); (-324.76, 37.5); (-194.856, -37.5); (-64.952, 37.5); (-64.952, 187.5) |
| 3 | (-64.952, 37.5); (-194.856, -37.5); (-194.856, -187.5); (-64.952, -262.5); (64.952, -187.5); (64.952, -37.5) |
| 4 | (194.856, 37.5); (64.952, -37.5); (64.952, -187.5); (194.856, -262.5); (324.76, -187.5); (324.76, -37.5) |

## Beams (model prediction, not an HFSS result)

| Port | Beam | Target (θ, φ) | Predicted (θ, φ) | Pointing error ° | Directivity dBi | Cross-pol dB |
|---|---|---|---|---|---|---|
| 1 | 1 | (30.0°, 45.0°) | 29.7°, 44.6° | 0.36 | 25.1 | - |
| 2 | 1 | (30.0°, 135.0°) | 29.9°, 134.8° | 0.14 | 25.3 | - |
| 3 | 1 | (30.0°, 225.0°) | 29.8°, 224.6° | 0.28 | 25.1 | - |
| 4 | 1 | (30.0°, 315.0°) | 29.9°, 315.0° | 0.10 | 25.4 | - |

## Cell states

| State | Cell | abs(Γ) dB | Phase at the patch plane |
|---|---|---|---|
| 0 | empty (bare substrate) | -0.001 | 138.1° |
| 1 | square patch 5.85 mm | -0.034 | -41.4° |

Separation 179.6°. Check the pair with one HFSS unit cell first (master/slave boundaries, 6.4 mm period, Floquet port, PEC ground, 10.54 GHz) and nudge the patch size if the separation is far from 180°.

## Files

- `board.step`
- `layout.csv`
- `feeds.csv`
- `layout.json`
- `layout.dxf`
- `layout.png`
- `hfss_build_board.py`

The STEP carries `Feed1`..`Feed4`: 0.6 x 0.6 mm marker sheets at z = 0 centred on each feed (listed above and in `feeds.csv`, `layout.json`, and on layer FEEDS of `layout.dxf`) (monopole pins were not requested).

`hfss_build_board.py` is a native AEDT script (HFSS: Tools > Run Script). It builds the substrate, ground with a coax opening at each feed, all patches united into one body "Patches" (Perfect E), the PEC monopoles with 4 lumped ports (50 ohm), an air box with a radiation boundary, a 10.54 GHz setup, a sweep and a far-field sphere. Design variables: `mono_h` (default 7.1 mm: set it to your monopole), `r_pin`, `r_out`, `pad_xy`, `pad_top`, `pad_bot`. Set `BUILD_MONOPOLES = False` at the top for geometry only. Creating a few thousand rectangles takes minutes. The script targets AEDT / HFSS 2019 R1 through 2025 R1 (sheets are created with `SolveInside` enabled because 2025 R1 rejects a zero-conductivity material such as vacuum with `SolveInside` off; Perfect E is assigned afterwards).

## Exciting the ports in HFSS

- **One port at a time**: HFSS > Excitations > Edit Sources, set port N to 1 W and every other port to 0 W (or simply read the S-parameters: S(N,N) is the match, S(M,N) the isolation to port M).
- **All ports together**: Edit Sources, give every port 1 W at the phase offsets you want (all 0° is the default); the far-field sphere then shows the superposed beams.
- The board is large (about 23 wavelengths across): use one port at a time first, an IE/FE-BI region or a mixed-order mesh, and compare with the model's beam table above.

Model-only: none of these numbers has been compared with an HFSS solve of this board.

## Phase 23 model predictions on the held-out fine grid (theta 0.5 deg x phi 1 deg)

Model-only; none of these numbers has been compared with HFSS.

**port1**: directivity 25.15 dBi, pointing error 0.36 deg, co-pol SLL -13.15 dB, x-cut -14.77 dB, y-cut -13.3 dB, XPD 18.6 dB

**port2**: directivity 25.35 dBi, pointing error 0.14 deg, co-pol SLL -11.86 dB, x-cut -11.87 dB, y-cut -16.04 dB, XPD 19.1 dB

**port3**: directivity 25.11 dBi, pointing error 0.28 deg, co-pol SLL -13.03 dB, x-cut -14.98 dB, y-cut -13.15 dB, XPD 17.1 dB

**port4**: directivity 25.42 dBi, pointing error 0.1 deg, co-pol SLL -12.47 dB, x-cut -12.5 dB, y-cut -15.13 dB, XPD 17.2 dB

**all**: directivity 23.17, 22.8, 23.13, 23.13 dBi, pointing error 0.22, 0.71, 0.41, 0.73 deg, co-pol SLL -23.32, -23.08, -22.76, -21.66 dB, x-cut -24.56, -25.88, -24.64, -24.17 dB, y-cut -23.09, -23.36, -22.99, -20.67 dB, XPD 17.0, 22.1, 16.5, 21.5 dB

## Port phases for the all-ports case (HFSS > Excitations > Edit Sources)

| Port | Magnitude | Phase |
|---|---|---|
| 1 | 1 W | 0.0 deg |
| 2 | 1 W | 0.0 deg |
| 3 | 1 W | 0.0 deg |
| 4 | 1 W | 0.0 deg |

Port 1 is the 0 deg reference. For a one-port-at-a-time run set that port to 1 W and the others to 0 W (the cell states were optimised for the mode named in this file).

Cell states optimised for ALL FOUR ports excited together.

**All-ports SLL needs ONE phase-locked source split four ways** (the cell states and the port phases above were optimised for the four ports radiating coherently with fixed relative phases). If the four ports are driven by independent signals the all-ports sidelobe figures do not apply: use the single-port metric (each port on its own) instead.
