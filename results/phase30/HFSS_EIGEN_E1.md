# Phase 30 - HFSS eigenmode of E1 vs model (E1: RO3010 0.635 mm, 3.2 mm, squircle n=4, slot 1.89 x 0.58 mm; phase 53 deg/cell; air box, auto lattice pairs)

HFSS (user-run, mode 1 bound, Mag_E cut checked): eps_r 10.2: x 10.861, y 10.232 GHz (mode 2: x 11.039, y 11.915); diag (53/53) 13.189. eps_r 11.2: x 10.478, y 9.862 GHz.
Model (mesh 48 / 64): eps_r 10.2 x 10.90 (10.912 root) / 11.14, y 10.26 (10.269) / 10.43; mode 2 x 11.064, y 11.965 (mesh 48). eps_r 11.2 (mesh 48, `eps11p2_hfss_followup.json`): x 10.528, y 9.876.
Mesh 48 vs HFSS: x +0.5 %, y +0.4 % (eps_r 10.2); x +0.5 %, y +0.1 % (eps_r 11.2); mode 2 +0.2 % / +0.4 %. Mesh 64 is 2.5-2 % high: the HFSS result favours mesh 48 (luck of the staircase, not convergence).
x/y split: HFSS 0.629 GHz, model 0.64 (m48) / 0.71 (m64). Pass criterion of the package (within 0.3 GHz of mesh-64, split within 30 %): PASSED.
**DIAGONAL (53/53) NOT reproduced**: HFSS 13.189 GHz vs model 11.19 (m48, TM fraction 0.89) / 11.36 (m64): +18 %. Unresolved; the model's oblique-azimuth Bloch tensor is unvalidated.

HFSS-derived (transverse resonance at each root, own f and beta = 53 deg/3.2 mm = 289.1 rad/m): x: f 10.861, beta/k0 1.270, X = kappa/(omega eps0) = 294.9 ohm, Y_sheet 17.10 mS (model m48 16.81, ratio 1.017; m64 ratio 1.10). y: f 10.232, beta/k0 1.348, X 340.5 ohm, Y 19.54 mS (model 19.27, ratio 1.014; m64 1.085).
At the design f the mean of the two roots is 10.547 GHz -> beta/k0 = 1.308 (model 1.302).
Anisotropy in the SAME metric as 4.6:1 (X_hi/X_lo of Im[1/(jY_par,perp + Y_slab)] at common f = 10.54 GHz and common beta = 1.302 k0): model 4.72:1 (X_par 192, X_perp 907 ohm, M_t 0.65); with the model Y rescaled by the HFSS/model ratios above (assumed constant in f): 6.0:1 (203/1226 ohm, M_t 0.72). Near-resonance quantity, very sensitive; the own-root X (295 vs 341 ohm, 1.16:1) is a different metric and not comparable.

eps_r sensitivity: HFSS f(fixed phase) -3.6 % (x) / -3.6 % (y) for eps_r 10.2 -> 11.2; model -3.5 % / -3.8 %: AGREE. The earlier "beta/k0 1.30 -> 1.57" is the shift at FIXED FREQUENCY (10.54 GHz, +21 % in beta); at fixed beta the same dispersion (v_g 0.16-0.18 c) gives -3.5 % in f (dbeta 60 rad/m x v_g / 2 pi = -0.5 GHz). Not a discrepancy; the package README line "drops roughly 0.8-1 GHz" was WRONG (-0.38 GHz measured). Consequence: at eps_r 11.2 the E1 geometry works at ~10.1 GHz (model beta/k0 1.326 at 10.1 GHz, M_t 0.70, X_bar 328 ohm), NOT at 10.54 GHz (M_t 1.22 there).
S1 at eps_r 11.2 (model, E1 geometry, 10.1 GHz, v_g from 9.9/10.3 GHz: 0.160 c): alpha/k0 0.0687 (first order) / 0.0342 (Floquet, ratio 0.50); Q = 1/(2 (alpha/k0)(v_g/c)) = 45 / 91; Im f = f/(2Q) = 0.112 / 0.055 GHz; f_r ~ 10.07 GHz (+-0.15). alpha[Np/m] = 2 pi Im f / v_g = pi f_r /(Q v_g); alpha/k0 = c/(2 Q v_g).
