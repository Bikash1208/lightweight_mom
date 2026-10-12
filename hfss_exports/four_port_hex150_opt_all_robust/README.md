# four_port_hex150_opt_all_robust

4-port multi-port HMSA board, 10.54 GHz. Substrate eps_r 2.2, tan d 0.0009, 1.575 mm thick, PEC ground. Cells 6.4 mm. The board is 652.8 x 537.6 mm (x -326.4..326.4, y -268.8..268.8), one monopole per sub-aperture.
2837 patches in total. Every coordinate is exact (rounded to 0.001 mm; the largest rounding applied was 0.0004736 mm), and each port's cells sit at that port's own positions, not on a shared lattice.

## Ports and feeds

| Port | Feed (x, y) mm | Sub-aperture x mm | Sub-aperture y mm | Cells | Patches |
|---|---|---|---|---|---|
| 1 | (64.952, 112.5) | -70.4 .. 198.4 | -38.4 .. 268.8 | 1414 | 700 |
| 2 | (-194.856, 112.5) | -326.4 .. -64 | -38.4 .. 268.8 | 1435 | 683 |
| 3 | (-64.952, -112.5) | -198.4 .. 70.4 | -268.8 .. 38.4 | 1414 | 715 |
| 4 | (194.856, -112.5) | 64 .. 326.4 | -268.8 .. 38.4 | 1435 | 739 |

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
| 1 | 1 | (30.0°, 45.0°) | 29.6°, 44.8° | 0.41 | 25.1 | - |
| 2 | 1 | (30.0°, 135.0°) | 29.8°, 135.0° | 0.20 | 25.2 | - |
| 3 | 1 | (30.0°, 225.0°) | 29.6°, 224.8° | 0.41 | 25.0 | - |
| 4 | 1 | (30.0°, 315.0°) | 29.7°, 315.0° | 0.30 | 25.2 | - |

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

**port1**: directivity 25.12 dBi, pointing error 0.41 deg, co-pol SLL -12.75 dB, x-cut -15.05 dB, y-cut -13.16 dB, XPD 18.6 dB

**port2**: directivity 25.21 dBi, pointing error 0.2 deg, co-pol SLL -12.26 dB, x-cut -12.26 dB, y-cut -12.65 dB, XPD 18.8 dB

**port3**: directivity 25.03 dBi, pointing error 0.41 deg, co-pol SLL -12.59 dB, x-cut -14.68 dB, y-cut -12.97 dB, XPD 17.7 dB

**port4**: directivity 25.24 dBi, pointing error 0.3 deg, co-pol SLL -11.25 dB, x-cut -11.25 dB, y-cut -14.29 dB, XPD 17.2 dB

**all**: directivity 23.43, 23.06, 23.38, 23.31 dBi, pointing error 0.2, 0.72, 0.2, 0.67 deg, co-pol SLL -22.54, -22.28, -22.34, -22.0 dB, x-cut -24.27, -24.57, -23.2, -23.94 dB, y-cut -23.12, -23.53, -23.52, -23.83 dB, XPD 16.9, 20.5, 16.9, 19.8 dB

## Port phases for the all-ports case (HFSS > Excitations > Edit Sources)

| Port | Magnitude | Phase |
|---|---|---|
| 1 | 1 W | 0.0 deg |
| 2 | 1 W | 0.0 deg |
| 3 | 1 W | 0.0 deg |
| 4 | 1 W | 0.0 deg |

Port 1 is the 0 deg reference. For a one-port-at-a-time run set that port to 1 W and the others to 0 W (the cell states were optimised for the mode named in this file).

Cell states optimised for ALL FOUR ports excited together. ROBUST: the worst case over 10.33 / 10.54 / 10.75 GHz (cells re-solved at each) and 20 seeded random port-phase errors of +-10 deg.

**All-ports SLL needs ONE phase-locked source split four ways** (the cell states and the port phases above were optimised for the four ports radiating coherently with fixed relative phases). If the four ports are driven by independent signals the all-ports sidelobe figures do not apply: use the single-port metric (each port on its own) instead.
