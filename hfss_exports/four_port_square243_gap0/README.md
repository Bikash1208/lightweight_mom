# four_port_square243_gap0

4-port multi-port HMSA board, 10.54 GHz. Substrate eps_r 2.2, tan d 0.0009, 1.575 mm thick, PEC ground. Cells 6.4 mm. The board is 486.4 x 486.4 mm (x -243.2..243.2, y -243.2..243.2), one monopole per sub-aperture.
2880 patches in total. Every coordinate is exact (rounded to 0.001 mm; the largest rounding applied was 0 mm), and each port's cells sit at that port's own positions, not on a shared lattice.

## Ports and feeds

| Port | Feed (x, y) mm | Sub-aperture x mm | Sub-aperture y mm | Cells | Patches |
|---|---|---|---|---|---|
| 1 | (121.6, 121.6) | 0 .. 243.2 | 0 .. 243.2 | 1444 | 717 |
| 2 | (-121.6, 121.6) | -243.2 .. 0 | 0 .. 243.2 | 1444 | 704 |
| 3 | (-121.6, -121.6) | -243.2 .. 0 | -243.2 .. 0 | 1444 | 723 |
| 4 | (121.6, -121.6) | 0 .. 243.2 | -243.2 .. 0 | 1444 | 736 |

## Beams (model prediction, not an HFSS result)

| Port | Beam | Target (θ, φ) | Predicted (θ, φ) | Pointing error ° | Directivity dBi | Cross-pol dB |
|---|---|---|---|---|---|---|
| 1 | 1 | (30.0°, 45.0°) | 30.0°, 45.0° | 0.00 | 26.0 | -16.6 |
| 2 | 1 | (30.0°, 135.0°) | 30.0°, 135.2° | 0.10 | 26.0 | -18.3 |
| 3 | 1 | (30.0°, 225.0°) | 30.0°, 225.0° | 0.00 | 26.0 | -17.7 |
| 4 | 1 | (30.0°, 315.0°) | 30.0°, 315.2° | 0.10 | 26.0 | -16.8 |

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
- The board is large (about 17 wavelengths across): use one port at a time first, an IE/FE-BI region or a mixed-order mesh, and compare with the model's beam table above.

Model-only: none of these numbers has been compared with an HFSS solve of this board.

## Phase 23 model predictions on the held-out fine grid (theta 0.5 deg x phi 1 deg)

Model-only; none of these numbers has been compared with HFSS.

**port1**: directivity 25.97 dBi, pointing error 0.0 deg, co-pol SLL -10.93 dB, x-cut -11.29 dB, y-cut -10.95 dB, XPD 16.4 dB

**port2**: directivity 26.03 dBi, pointing error 0.1 deg, co-pol SLL -11.13 dB, x-cut -11.63 dB, y-cut -11.13 dB, XPD 16.5 dB

**port3**: directivity 25.97 dBi, pointing error 0.0 deg, co-pol SLL -10.9 dB, x-cut -11.29 dB, y-cut -10.91 dB, XPD 16.4 dB

**port4**: directivity 26.03 dBi, pointing error 0.1 deg, co-pol SLL -11.2 dB, x-cut -11.69 dB, y-cut -11.2 dB, XPD 16.6 dB

**all**: directivity 20.19, 20.52, 20.17, 20.57 dBi, pointing error 0.6, 0.14, 0.6, 0.14 deg, co-pol SLL -8.65, -9.19, -8.64, -9.26 dB, x-cut -8.99, -11.16, -9.07, -11.23 dB, y-cut -9.59, -9.83, -9.55, -9.89 dB, XPD 15.2, 17.0, 15.0, 17.3 dB

Starting design (Phase 16 pointing-corrected, no sidelobe optimisation).
