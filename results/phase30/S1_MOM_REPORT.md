# S1 direct MoM supercell solve (MODEL_ONLY, 2026-11-04)

Question: HFSS gave S1 Q = 403 (alpha/k0 ~ 0.0078) against the Floquet-tensor model's 45-91 (alpha/k0 0.0373). Is the locality assumption (periodic-neighbour tensors while the slot turns 25.7 deg/cell with 0.15 mm gaps) the cause?

Method (`src/dispersion/leaky_supercell_mom.py`, additive; `scripts/research/phase30_s1_supercell_mom.py`; tests `tests/test_leaky_supercell_mom.py`): the whole 7-cell supercell (22.4 x 3.2 mm, E1 outline fixed, slot angle 180 i/7 deg, grounded RO3010 0.635 mm, eps_r 11.2, F = 2 x mesh per cell, mesh = rooftop cells per 3.2 mm cell) is one periodic MoM problem. The impedance matrix is continued analytically (outgoing branch for radiating harmonics) and the singular point is found in the complex plane: f-mode = real kx0 (0.4 deg/supercell, 0.31 rad/m), complex f (Q = f_r / 2 f_i); k-mode = real f (= f_r of the f-mode), complex kx0 (alpha = |Im kx0|). Z fill verified equal to the production assembler to 1e-12; the uniform E1 cell reproduces the tensor-model frequencies (m48: y 9.870 vs model 9.876 GHz; x 10.610 vs 10.528; HFSS 9.862 / 10.478).

| mesh | n unknowns | f_r GHz (tan d 0) | Q_rad (tan d 0) | Q_tot (0.0022) | alpha_rad/k0 | alpha_tot/k0 |
|---|---|---|---|---|---|---|
| 24 | 5.4k | 10.442 | 44.4 | 41.1 | 0.0536 | 0.0574 |
| 32 | 9.8k | 10.240 | 36.2 | 34.0 | 0.0797 | 0.0836 |
| **48** | 23.1k | **9.981** | **17.4** | **16.7** | **0.0752** | **0.0790** |

(m24/m32 JSON `derived` blocks come from the first script version, signed v_g; use the stage entries.) Mesh 48 is the mesh of the model's 0.0373 and of the 0.5 % HFSS agreement on E1; Q_rad moves a factor 2-3 between m32 and m48 (gap staircase 0.2 -> 0.133 mm), alpha/k0 only 6 %, so alpha is the trustworthy number.

Comparison at m48 (alpha/k0): direct MoM 0.0752 (tan d 0) / 0.0790 (0.0022); Floquet-tensor model 0.0373 (first order 0.0665); HFSS-implied 0.0078. The direct solve gives **2.0x the Floquet value and 9.6x the HFSS-implied value**.

LOCALITY: ruled out as the cause of the HFSS deficit. Dropping locality makes the leakage stronger, not weaker (the Floquet staircase underestimates it by ~2x). The mode radiates broadside (harmonic kt = beta - 2pi/L ~ 0), runs toward -x for this rotation sense (v_g|MoM| ~ 0.29 c at m48 from the 0.4 -> 6.74 deg pair; strongly dispersive, so the Q/v_g conversion is only indicative).

Q split expected from HFSS (m48 MoM): Q_rad 17 (Im f 0.29 GHz), Q_total 16.7, Q_dielectric ~ 400 (m24/m32: 555-563; 1/Q_d = 1/Q_tot - 1/Q_rad is a small difference of large terms at Q ~ 17, so 400-560 +-). Dielectric energy fraction p_e = 1/(Q_d tan d) ~ 0.8 (m24/m32), >0.8 and not meaningful beyond 1 at m48. Hence HFSS Q = 403 equals the pure dielectric Q: **the HFSS S1 mode behaves as if it did not radiate**. Falsifiable predictions for the pending HFSS runs: tan d = 0 with the same radiation top -> Q stays ~ 400-600 or higher if the boundary does not absorb the broadside leakage; PML top -> Q ~ 15-40, f_i ~ 0.15-0.3 GHz; PEC lid -> a cavity (no decay). Suspect HFSS setup (radiation-boundary distance / reflection, the lid) before the model.

Labels: MODEL_ONLY, NOT_HFSS_VALIDATED. Not done: m64 (41k unknowns, ~1 h per stage), ky0 != 0, mesh-converged Q.

## Addendum (2026-11-04): vs the HFSS PML run, mesh and F convergence (MODEL_ONLY)

HFSS S1 PML top (eps_r 11.2, ph_x 6.74 deg): f = 9.394/9.655/10.558/10.954/12.747/12.875 GHz, Q = 830/286/282/477/169/176; near 10 GHz Q ~ 280, alpha/k0 ~ 0.011 total.

Q_rad of the bright x-mode (tan d = 0, kx0 = 6.74 deg/supercell, F = 2n), logs `s1_mom_Q_vs_mesh_F_kx6p74.log`:

| mesh | N | f (GHz) | Q_rad |
|---|---|---|---|
| 24 | 5174 | 10.424 | 49 |
| 32 | 9756 | 10.15 | 32 |
| 40 | 15712 | 10.06 | 30 |
| 48 | 23092 | 10.05 | 19 |
| 64 | 40182 | 10.196 | 38.5 |

- Q_rad 17 does NOT survive m64: Q_rad(m) scatters 19-49 (non-monotone; gap staircase), best estimate 30-40 (alpha/k0 ~ 0.04 by the alpha ~ 1/Q scaling at fixed v_g: not re-solved). Still 7-9x below the HFSS Q ~ 280.
- F: m24 F = 1n/2n/3n/4n -> Q 48.9/48.9/49.1/49.1; m32 F = 1n/2n/3n -> 32.4/32.5/32.6. The Floquet truncation is not the cause (< 1 %).
- Mode match (m24 scan at Im f = 0.08, `s1_mom_m24_modescan_imf0p08.log`): bright x-mode 10.4-10.55 <-> HFSS 10.558 (Q 282); MoM 10.8 <-> HFSS 10.954; MoM 12.74 (Q 402) <-> HFSS 12.747 (Q 169); 12.95 <-> 12.875. The other HFSS modes (9.39, 9.655) lie below the scanned MoM peaks at m24 (converged-mesh modes shift down 0.2-0.5 GHz). The mode-by-mode match is by frequency only, no eigenvector comparison.
- Where the MoM radiation could be overestimated (untested): (1) staircase of the 0.15 mm gap and of the rotated slot (Q scatter +-40 % with mesh; the near-field of a narrow gap sets the n = 0 harmonic coupling); (2) infinite lateral periodicity (no edge, no truncation of the leaky wave by a finite supercell count: Bloch real kx0 with the n = 0 harmonic always in the light cone); (3) ideal PEC ground and zero-thickness PEC pattern (HFSS metal and finite-thickness ground); (4) the n = 0 harmonic is the only radiating order, so any error in the lowest-harmonic substrate/air Green function at small kt directly scales Q_rad (checked equal to the production assembler for real k, 1e-12); (5) HFSS side: PML distance and eigen-solver damping accuracy, 6-mode request (a Q ~ 20-40 mode may not be returned), mesh of the 0.15 mm gap.
- Not done: eigenvector (field) comparison, an independent far-field power check of Q_rad, a finer-than-m64 mesh, a staircase-free (curved-edge) basis.
