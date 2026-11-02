# Phase 30 item 2 - cell screening report (MODEL_ONLY, NOT_HFSS_VALIDATED)

10.54 GHz, grounded single-layer patch cells, our MoM (Bloch sheet tensor at kt = beta0, `src/dispersion/tensor_sheet.py`),
mesh 48, 696 resolved cells (gap >= 1.5 mesh cells; unresolved gaps skipped - they gave unphysical negative admittance).
Data: `results/phase30/screen/cell_screen.json` (+ `.jsonl`), `stageb.json`. Scripts: `scripts/research/phase30_cell_screen.py`,
`phase30_cell_stageb.py`. New geometry: `src/dispersion/elliptical_cells.py` (ellipse / super-ellipse, slot, cross slot; 14 tests).

Metrics: M_t = |X_par - X_perp| / |X_par + X_perp| (surface-wave reactances along / across the cell axis), anisotropy ratio
X_hi/X_lo, alpha from the first-order rotating-axis tensor formula (`leakage_alpha_first_order_tensor`, theta 17 deg), Floquet
harmonic-balance alpha (`tensor_floquet`) as cross-check. Edge power on a hexagon of effective radius 70 mm:
eps_uniform = exp(-2 alpha R) (M_max everywhere); eps_taper_bound = 1/(1 + alpha R) is the best possible for UNIFORM aperture
illumination (a hard lower bound for that goal).

## Ranked (best per geometry x substrate, pitch in brackets)

| # | Cell | M_t | aniso X_hi/X_lo | X_bar ohm / beta/k0 | alpha/k0 first-order (Floquet) | eps_uniform 1st (Floq) | eps_taper_bound | rad/diel |
|---|---|---|---|---|---|---|---|---|
| 1 | RO3010 0.635 mm [3.2 mm], squircle (n=4) + centre slot 0.58 mm along the long axis, patch ratio 1.0, gap 0.15 mm | 0.64 | 4.6 | 314 / 1.30 | 0.058 (0.043, ratio 0.73 OK) | 0.16 (0.27) | 0.53 | ~50 |
| 1b | same, gap 0.10 mm | 0.76 | 7.5 | 351 / 1.37 | 0.059 (Floquet 0.16, unreliable) | 0.16 (<=0.01?) | 0.52 | ~50 |
| 2 | RO3010 [3.2 mm], rounded rectangle (n=8) ratio 0.92, gap 0.10 mm | 0.52 | 3.2 | 188 / 1.12 | 0.024 (Floquet root jumped, unreliable) | 0.48 | 0.73 | 73 |
| 3 | RO3010, circle + slot (the current cell's geometry), gap 0.10 | 0.27 | 1.75 | 152 / 1.08 | 0.0035 (0.009) | 0.90 (0.75) | 0.95 | 19 |
| 4 | RO3006 0.635 mm [4.0], squircle + slot, gap 0.15 | 0.28 | 1.8 | 146 / 1.07 | 0.0058 (0.020) | 0.84 (0.55) | 0.92 | 28 |
| 5 | RT5880 1.575 mm [6.4] rounded rectangle ratio 0.75, gap 0.20 | 0.30 | 1.9 | 118 / 1.05 | 0.0089 (0.0065) | 0.76 (0.82) | 0.88 | 75 |
| 6 | RT5880 [6.4] circle + slot (baseline-type), gap 0.20 | 0.15 | 1.35 | 134 / 1.06 | 0.0025 (0.0012) | 0.93 (0.96) | 0.96 | 15 |

Other geometries (best alpha/k0 over all substrates/sizes): ellipse 0.0007 (brief candidate a), ellipse + slot 0.0014 (brief
candidate b), circle + cross slot 0 (isotropic: AR = 1, no tensor, so useless for CP), squircle 0.0047, rounded rectangle 0.024,
squircle + slot 0.059.

## What this says
1. The Phase 26-type cell is confirmed weak: baseline-type cells give alpha/k0 ~ 0.001-0.003 and anisotropy 1.2-1.4:1 (matches the
   HFSS finding that the surface wave reaches the edge).
2. Brief candidates (a) elliptical patch and (b) elliptical patch + slot are poor at every substrate here: an ellipse fills only
   ~79 % of the cell, so it cannot get near the capacitive cancellation of the slab inductance. The cell shape that works is a
   FULL-fill patch (super-ellipse / rounded rectangle) with a centre slot; the tensor swing comes from operating near the
   X_par -> large resonance, where X_par is large and X_perp stays small.
3. Substrate: only thin high-er (RO3010, er 10.2, 0.635 mm) reaches the strongly bound regime (beta/k0 1.3-1.4, X_bar 300-350 ohm,
   close to the -450 ohm of Faenzi's wideband design). RO3006 and RT5880 do not get there at printable gaps (>= 0.15-0.2 mm).
4. Dielectric loss is small for the top cells (radiated/dielectric alpha ~ 50); conductor loss NOT estimated (PEC model).
5. Directional selectivity: no spatial-dispersion selectivity found. The deviation from the local rotated-tensor model is
   0.07-0.13 of the anisotropy for the RO3010 cells (0.33 for the baseline-type cell). Beta does vary with arrival angle
   by 0.25-0.43 k0 for the top cells (that is the tensor itself, not selectivity) - a large anisotropy means a strongly
   direction-dependent local wave speed, which the existing slow-modulation synthesis does not handle; this needs care.

## Does a top cell meet the goal on a 150 mm hexagon?
- alpha/k0 >= 0.03: yes for #1 (0.043 Floquet-checked at gap 0.15 mm), plausibly #1b.
- < 10 % edge power: NOT met by a Floquet-checked cell. With M_max everywhere the edge keeps 16 % (first order) to 27 % (Floquet).
  An ideal taper that radiates the power UNIFORMLY across the aperture cannot go below ~53 % edge power with this alpha
  (bound 1/(1 + alpha R)). So "< 10 % at the edge" and "uniform aperture illumination" are incompatible at R = 70 mm for any
  screened cell: a board of this size must accept a strongly tapered aperture (lower aperture efficiency) or an absorbing/
  leaky rim. Reaching 10 % needs about alpha/k0 >= 0.07 with uniform M (edge) - roughly 1.6x the best verified cell.
- Anisotropy >= 3:1: yes for #1, #1b, #2. M_t up to 0.64-0.76 (target ~0.5).
- Single layer, printable, pitch <= lambda_g/5: pitch 3.2 mm = lambda0/8.9; lambda_sw/5 at beta/k0 = 1.3 is 4.4 mm, OK.

## Cautions (do not skip)
- The best cells sit at gap 0.10-0.15 mm. X_bar for the same cell goes 180 -> 314 -> 351 ohm for gap 0.20 -> 0.15 -> 0.10 mm:
  etch tolerance of +-0.05 mm changes the surface reactance by tens of percent. Process capability must be checked first
  (fabrication is out of this session's scope).
- At M_t 0.6-0.8 the first-order alpha is outside its tested range (<= 3 % at Delta Y/Y 0.4); the Floquet cross-check agrees (0.73)
  at gap 0.15 mm but its root jumped to another mode for 4 of 8 cells (alpha/k0 > 0.5 or ratio > 40: discarded as unreliable,
  not used). Single-cell tensors are used for all harmonics (no spatial dispersion between harmonics).
- Mesh 48: resolved gaps only; mesh convergence of the near-resonant cells was NOT checked (32 vs 48 differs ~2 % in Y for a mid cell;
  near resonance X is far more sensitive). Item 3 (HFSS eigenmode check) is the real test.
- M(size) near zero: not mapped. Only two slot widths were sampled, so smoothness of M from 0 is NOT_AVAILABLE_FROM_CURRENT_DATA.
- Alternative cause not excluded: part of the near-surface power in the HFSS boards may be the monopole's direct radiation / edge
  diffraction rather than an un-leaked surface wave (no launch-fraction model in the repo, see results/phase26/README.md).
  A better cell fixes only the leakage part.
- Rogers loss tangents are datasheet-nominal values entered from memory (RO3006 0.0020, RO3010 0.0022), not re-verified.

## Recommendation
Go with candidate #1 (RO3010, 3.2 mm pitch, squircle + centre slot) to item 3, but first agree: (i) whether the 0.10-0.15 mm gap is
acceptable to fabricate, (ii) accept ~16-27 % edge power with a tapered aperture (or add a rim), and (iii) a refined
(finer slot/gap, mesh-converged) screen of the #1 family before the HFSS package. Runner-up: #2 (rounded rectangle, no slot; simpler,
alpha/k0 0.024 unverified). The original three substrates/geometries of the brief other than RO3010 + full-fill + slot do not reach the goal.
