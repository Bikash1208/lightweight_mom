# sub6_edge_diode_eigenmode: HFSS eigenmode package, best 90 deg edge-diode row

**MODEL_ONLY, IDEAL_BIAS, NOT_HFSS_VALIDATED. Predictions were written before any HFSS run.**

> **GATED / PROVISIONAL.** This package was generated before the resonance + mesh 64/96 check (`results/sub6_edge_diode_check/`) finished. It must not be used until that check passes; the row and phases may change.

Row: `EDGE_th90_q0.7_W2.4` (best 90 deg row on mesh 64: alpha_OFF 2.96 Np/m, alpha_ON 0.21 Np/m, ratio 14, class OK_OFF_ONLY).
RT5880 (eps_r 2.2, tan d 0.0009) 3.175 mm on a PEC ground, p = 11.6 mm, disc r = 5.30 mm, slot 7.42 x 2.4 mm along y (rim intact), 5.8 GHz.

## Geometry (cell window shifted by p/2 in x)
- Window x, y in [-5.8, 5.8] mm. The disc centres are on the window edges x = -+5.8, so the window holds half of each of two discs, the 1.00 mm gap and the diode, at the window centre. No vias.
- Each half disc carries a half slot (width 2.4 mm in x, length 7.42 mm in y) on the window edge.
- Diode sheets at z = 3.175 (patch plane): `DiodeA` x -0.50..0.00 and `DiodeB` x 0.00..0.50, both 0.35 mm wide in y. Two 0.02 mm PEC tabs (`TabL`, `TabR`) straddle the disc edges, because the circular edge is up to 0.003 mm off the straight sheet end over +-0.175 mm.
- AIR box from the substrate top to 60 mm with a PEC lid. There is no air below the ground, which is a PEC sheet at z = 0.
- `cell_edge_diode.step`: exact-coordinate STEP (substrate, ground, patches, DiodeA, DiodeB). The 40 3-D vertex coordinates sit on the 0.001 mm grid (checked by `verify_step_exact`, vertices only). The arcs are analytic circles, so their pcurve parameters are not on the grid, as in the other circle exports.

## Diode model (as in the MoM): (C || R) + L
HFSS has no series-L lumped element, so the diode is two Lumped RLC sheets in series (current line along x), which is exact:

| State | RLC_A (`DiodeA`) | RLC_B (`DiodeB`) |
|---|---|---|
| OFF | R 5000 ohm, C 0.15 pF (parallel) | L 0.6 nH |
| ON | R 3 ohm | L 0.6 nH |

Per-sheet width/length is 0.35 / 0.50 mm. HFSS applies the RLC value to the whole boundary, so the series total is the stated R, C, L regardless of sheet aspect ratio.

## Runs (crucial first)
Master/slave in x and y, a PEC lid, and the bound mode = the lowest eigenmode whose field decays away from the patch plane. The Bloch phase is beta(state, direction) * p, with beta from the MoM sheet admittance (Y_uu for kt along the slot axis y = "along", Y_vv for kt along x = "across") at 5.8 GHz, mesh 64.

| Run | Script | phase_x | phase_y | beta/k0 | Predicted f |
|---|---|---|---|---|---|
| 1 | `build_EDGE_OFF_across.py` | 103.997 | 0 | 1.287 | 5.8 GHz |
| 2 | `build_EDGE_ON_across.py` | 81.42 | 0 | 1.009 | 5.8 GHz |
| 3 | `build_EDGE_OFF_along.py` | 0 | 84.82 | 1.050 | 5.8 GHz |
| 4 | `build_EDGE_ON_along.py` | 0 | 84.82 | 1.050 | 5.8 GHz |

Runs 1 and 2 (across, kt along the diode axis) are the crucial ones: that is the only direction where the diode state changes the dispersion. Along the slot axis the diode carries no current, so ON and OFF give the same phase (84.82 deg) and should give the same eigenfrequency.

Each script is a copy of the proven `phase26_eigenmode_package.py` template: run it in AEDT (Tools > Run Script), check the master/slave dialogs, then Analyze. Record f and Q (the diode adds loss, so the eigenvalue is complex) in a copy of `predictions.json`.

## Caveats
- The author did not run these scripts in AEDT.
- The line of diodes is periodic: every cell has a diode to its neighbour.
- Phases and the row choice come from mesh 64. See `results/sub6_edge_diode_check/` for the mesh 96 and resonance check.
- Eigenmode f is a model-vs-HFSS check of the dispersion, not an S-parameter validation.
