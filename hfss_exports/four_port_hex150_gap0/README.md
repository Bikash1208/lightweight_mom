# four_port_hex150_gap0

4-port multi-port HMSA board, 10.54 GHz. Substrate eps_r 2.2, tan d 0.0009, 1.575 mm thick, PEC ground. Cells 6.4 mm. The board is 652.8 x 537.6 mm (x -326.4..326.4, y -268.8..268.8), one monopole per sub-aperture.
2847 patches in total. Every coordinate is exact (rounded to 0.001 mm; the largest rounding applied was 0.0004736 mm), and each port's cells sit at that port's own positions, not on a shared lattice.

## Ports and feeds

| Port | Feed (x, y) mm | Sub-aperture x mm | Sub-aperture y mm | Cells | Patches |
|---|---|---|---|---|---|
| 1 | (64.952, 112.5) | -70.4 .. 198.4 | -38.4 .. 268.8 | 1414 | 696 |
| 2 | (-194.856, 112.5) | -326.4 .. -64 | -38.4 .. 268.8 | 1435 | 719 |
| 3 | (-64.952, -112.5) | -198.4 .. 70.4 | -268.8 .. 38.4 | 1414 | 717 |
| 4 | (194.856, -112.5) | 64 .. 326.4 | -268.8 .. 38.4 | 1435 | 715 |

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
| 1 | 1 | (30.0°, 45.0°) | 29.6°, 45.0° | 0.40 | 25.9 | -17.1 |
| 2 | 1 | (30.0°, 135.0°) | 29.6°, 135.0° | 0.40 | 26.0 | -17.8 |
| 3 | 1 | (30.0°, 225.0°) | 29.6°, 225.0° | 0.40 | 25.9 | -18.1 |
| 4 | 1 | (30.0°, 315.0°) | 29.6°, 315.0° | 0.40 | 26.0 | -17.0 |

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

**port1**: directivity [25.91] dBi, pointing error [0.4] deg, co-pol SLL [-13.86] dB, x-cut [-16.17] dB, y-cut [-13.87] dB, XPD [17.1] dB

**port2**: directivity [25.98] dBi, pointing error [0.4] deg, co-pol SLL [-14.23] dB, x-cut [-16.43] dB, y-cut [-15.47] dB, XPD [17.0] dB

**port3**: directivity [25.9] dBi, pointing error [0.4] deg, co-pol SLL [-13.73] dB, x-cut [-16.01] dB, y-cut [-13.74] dB, XPD [17.1] dB

**port4**: directivity [25.99] dBi, pointing error [0.4] deg, co-pol SLL [-14.3] dB, x-cut [-16.51] dB, y-cut [-15.58] dB, XPD [17.0] dB

**all**: directivity [20.05, 19.97, 20.02, 20.0] dBi, pointing error [0.36, 0.41, 0.45, 0.4] deg, co-pol SLL [-12.99, -12.58, -13.13, -12.48] dB, x-cut [-15.91, -15.33, -15.62, -15.24] dB, y-cut [-13.03, -16.66, -13.17, -16.85] dB, XPD [16.8, 18.1, 16.9, 17.9] dB

Starting design (no sidelobe optimisation).
