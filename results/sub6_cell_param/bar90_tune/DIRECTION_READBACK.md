# BAR90 tune — direction and via read-back (solved model, AEDT 2025.1)

Read back from the saved BASE model (project copy, no solve) on 2026-11-12. Geometry is the BASE config; all
seven configs share the same construction.

## Phase axes
| Label | Project variable | Lattice pair | Faces (global) | Phase runs along |
|---|---|---|---|---|
| across | `$ph_x` | `Lattice_X` | x = -5.8 / +5.8 mm (normals ±X) | global **X** |
| along  | `$ph_y` | `Lattice_Y` | y = -5.8 / +5.8 mm (normals ±Y) | global **Y** |

The boundary property `Phase` reads `$ph_x` on `Lattice_X` and `$ph_y` on `Lattice_Y`. A row sets `($ph_x,$ph_y)` = (phi, 0)
for across and (0, phi) for along; each row's JSON stores both readbacks in SI units.

## Slot long axis
Slot cut is |x| < 1.45 mm by |y| < 4.348 mm (slot_w = 2.9, slot_l = 8.696 mm), so the long axis is global **Y**.
Point probes at z = 17.5 um: (0, 4.0) is air only, (1.3, 0) is air only, (1.6, 0) and (0, 4.6) are DISC metal.
**"Along" = phase parallel to the slot long axis (Y), "across" = phase perpendicular to it (X).** This matches the
calibration frame in `scripts/research/sub6_island_cell_packages.py` ("slot along y. ACROSS: (ph_x,ph_y)=(phi,0);
ALONG: (0,phi)").

## Diode gap and current direction
`DIODE_SHEET` is a rectangle in the XY plane at z = 17.5 um: x from -1.45 to -1.15 mm (0.3 mm), y from -0.3 to +0.3 mm
(0.6 mm). Its -X end touches DISC metal (probe at x = -1.46: DISC); its +X end touches ISLAND metal (probe at x = -1.14:
ISLAND). So current crosses the gap **along the X axis**. By the build's naming convention the island is the anode side and
the disc the cathode side, so forward current flows from the island to the disc, i.e. in **-X**. Note the boundary
`BAR90_02LRH` is an isotropic surface impedance (Rs, Xs) with no polarity of its own; the direction comes from the geometry
only. OFF: Rs = 22.23 ohm, Xs = -335.4 ohm. ON: Rs = 1.60 ohm, Xs = 29.2 ohm.

## Via bottom pads (neither grounded)
Island via at (x, y) = (-0.15, -1.34) mm (the local-frame x = -1.34 mm point after the +90 degree rotation); disc via at
(3.63, -3.63) mm.
- Each via barrel (r = 0.15 mm) runs from the bias-layer bottom to the top metal and is unioned with its patch: the island
  barrel and pad are part of ISLAND (probe at the pad centre: ISLAND), the disc barrel and pad are part of DISC (DISC).
- Each bottom pad (r = 0.3 mm, 35 um thick, z from -3.499 to -3.464 mm) has only air below it (probe 0.5 mm below the pad:
  air only). The AIR box extends to z = -153.5 mm where the PEC_BOTTOM_150 boundary sits, 150 mm away.
- GROUND (z -3.21 to -3.175 mm) has an antipad of r = 0.5 mm around each via: probe at r = 0.4 mm is air, at r = 0.6 mm is
  GROUND. The barrel crosses the ground plane through the antipad without touching it.
- Result: the island pad connects only to the island, the disc pad only to the disc, and neither is connected to GROUND or
  to the other.

Boundaries present: BAR90_02LRH (impedance), Lattice_X, Lattice_Y, PEC_LID (z = 150.035 mm), PEC_BOTTOM_150 (z = -153.499 mm).
