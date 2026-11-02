# Phase 30 - refined screen of the approved cell (MODEL_ONLY, NOT_HFSS_VALIDATED)

Cell: Rogers RO3010 0.635 mm (eps_r 10.2, tan d 0.0022; `DATASHEET_RO3010.md`), 3.2 mm square lattice, full-fill super-ellipse
n = 4, major 3.05 mm along x (gap 0.15 mm), minor = ratio x major, centre slot along x (L x W). Scripts: `scripts/research/phase30_refine.py`
(`grid`, `stage`), `phase30_rotation_check.py`. Data: `results/phase30/refine/` (`refine_grid.json` 87 cells at mesh 48, `refine_stage.json`
7 picks x 3 frequencies / 3 meshes / gap +-0.025 / rotation, `rotation_check.json`, `eps11p2_sensitivity.json`, `edge_fraction_tapered.json`).

## Reference cell E1 (ratio 1.0, slot 1.89 x 0.58 mm), mesh 48
| f GHz | beta/k0 | alpha/k0 first-order | alpha/k0 Floquet (ratio) | M_t | X_hi/X_lo | X_bar ohm |
|---|---|---|---|---|---|---|
| 10.30 | 1.214 | 0.0364 | 0.0345 (0.95) | 0.47 | 2.8 | 259 |
| 10.54 | 1.302 | 0.0584 | 0.0310 (0.53) | 0.64 | 4.6 | 314 |
| 10.80 | 1.420 | 0.0583 | 0.0098 (0.17) | 0.88 | 15 | 380 |
The cell is a narrow-band resonant element: 5 % in frequency moves beta/k0 by ~0.1 and M_t by ~0.2; at 10.8 GHz (M_t 0.88) the first-order alpha is
invalid (Floquet/first-order 0.17). Use the Floquet value where the ratio is 0.5-1.5.

## M(size) down to M ~ 0 (ratio 1.0, mesh 48; slot width 0.10 mm, slot length sweep)
| slot L mm | 0 | 1.22 | 1.53 | 1.89 | 2.29 |
|---|---|---|---|---|---|
| M_t | 0 | 0.064 | 0.117 | 0.266 | 0.621 |
| alpha/k0 (1st) | 0 | 0.0006 | 0.0020 | 0.0113 | 0.0592 |
| X_bar ohm | 215 | 226 | 234 | 256 | 301 |
Monotone and smooth in slot length, continuous to M = 0 (no slot). Caveat: slot widths 0.10 and 0.20 mm give identical results at mesh 48
(mask quantisation, cell = 0.067 mm), and the weak-M values shrink at finer meshes (M_t 0.064 -> 0.040 at mesh 64): the low-M end is
quantisation-limited. X_bar rises with M (215 -> 314 ohm), so a taper needs a compensating size change to keep X_bar constant (not done).
Slot lengths above 0.75 x major leave < 0.3 mm bridges and are invalid (a first test at 0.9 x major gave unphysical admittances).

## Mesh convergence: NOT converged (converged = n)
| cell | mesh 48 | 64 | 80 | spread (max-min)/mean |
|---|---|---|---|---|
| E1 beta/k0 | 1.302 | 1.208 | 1.265 | 7.5 % |
| E1 alpha/k0 (1st) | 0.0584 | 0.0419 | 0.0629 | 38 % |
| E1 X_bar ohm | 314 | 255 | 292 | 21 % |
| E2 (L2.29 W0.10) alpha/k0 | 0.0592 | 0.0197 | 0.0295 | 100 % |
Non-monotonic: the rooftop mask samples the outline and the slot at mesh-cell centres, so features (gap 0.15 mm, slot 0.58 mm) change by whole cells
(0.067 / 0.050 / 0.040 mm) from mesh to mesh. Edges align with cell boundaries only for gap 0.10 or 0.20 mm and slot dimensions that are multiples of 0.1 mm
(meshes 64 and 128); mesh 128 needs ~13 GB per process and ~3 h and was not run. Quantitative values above are therefore +-10 % in X_bar,
+-4 % in beta, a factor ~1.5-2 in alpha; HFSS is the arbiter. A conforming-mask or sub-cell treatment of the MoM is NOT part of this task.

## Gap +-0.025 mm (E1, mesh 48; staircase-limited: the shift is below one mesh cell)
| gap | X_bar ohm | M_t | beta/k0 |
|---|---|---|---|
| 0.125 | 325 (+11) | 0.677 (+0.036) | 1.321 (+1.5 %) |
| 0.175 | 277 (-37) | 0.528 (-0.113) | 1.241 (-4.7 %) |
Both signs and sizes are uncertain (asymmetric because of quantisation); the coarse grid slope (gap 0.20 -> 0.15 -> 0.10 mm: X_bar 180 -> 314 -> 351 ohm for the
same slot) suggests ~35-70 ohm per 0.025 mm near 0.15 mm.  Read as: an etch tolerance of +-0.025 mm changes X_bar by tens of ohm (10-20 %) and M_t by 0.05-0.1.

## Substrate epsilon_r sensitivity (E1; `eps11p2_sensitivity.json`)
eps_r 10.2 -> 11.2 (the datasheet 'design Dk'): beta/k0 1.30 -> 1.57, X_bar 314 -> 456 ohm, M_t 0.64 -> 1.22 (leaves the bound-mode regime). The
+-0.30 process tolerance alone is worth ~+-0.08 k0 in beta and ~+-45 ohm. A measured eps_r of the actual lot is a prerequisite for a build.

## Tensor rotation (E1, pick-5, pick-3; mesh 48, kt = 1.1 k0)
Slot-only rotation (outline fixed, as in the board): relative error 2.4 % / 3.4 % (22.5 / 45 deg) for E1, 1.5 / 2.1 % and 1.0 / 1.7 % for the others, 0 at 90 deg.
Whole-cell rotation: 1.14-1.30 (invalid: a 4-fold super-ellipse is not rotation invariant). Designs must keep the outline fixed and turn only the slot.
(The `rot` entries of `refine_stage.json` are the whole-cell version from the first run; `rotation_check.json` is the correct one; `phase30_refine.py`
now runs the slot-only version.)

## Losses (E1, 10.54 GHz, mesh 48)
alpha_rad(1st) 12.9 Np/m (Floquet 6.8); alpha_diel 0.22 Np/m (tan d 0.0022); alpha_cond 0.43 Np/m = patch 0.20 + ground 0.23 (Cu, sigma 5.8e7 S/m, smooth, Rs 25.6 mohm/sq; thickness 18 um >> 0.64 um skin depth).
Ratios alpha_rad/alpha_loss: dielectric 59 (Floquet-based 31), conductor 30 (16). Conductor loss is an ESTIMATE and a lower bound (rooftop currents under-resolve edge
crowding at the 0.15 mm gaps, roughness ignored; a rough-copper factor of 1.2-1.5 would be typical) -> total loss <= ~10 % of the leaked power over a 70 mm path.

## Edge power
Uniform M_max everywhere, R = 70 mm: exp(-2 alpha R) = 16 % (first-order alpha 12.9 Np/m) ... 38 % (Floquet alpha 6.8 Np/m).
Tapered (ray-march design, table of 22 cells, M_t <= 0.66; `edge_fraction_tapered.json`): 150 mm hexagon 28-30 % (first-order alpha) ... 45-46 % (alpha x 0.53);
90 mm test board 39-40 % ... 58 %. The < 10 % goal is NOT met by this cell in the model (needs ~1.6-3x more alpha).

## Other picks (mesh 48; Floquet-checked at 10.54 GHz)
E2 L2.29 W0.10: alpha/k0 0.059 (Floquet 0.036), M_t 0.62; E3 L1.89 W0.10: 0.0113 (0.0117), M_t 0.27; E4 L1.53 W0.10: 0.0020 (0.0025), M_t 0.12;
pick-5 L1.89 W0.45: 0.0416 (0.0320), M_t 0.51; pick-6 ratio 0.97 L2.29 W0.58: 0.0534 (0.0409) but mesh 64/80 give M_t 1.2/0.96 (unconverged, leaves the valid regime).
