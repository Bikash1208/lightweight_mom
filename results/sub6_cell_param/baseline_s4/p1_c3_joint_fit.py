"""Phase 1 c: JOINT fit of (beta_tA, beta_tB) to the HFSS beam angles, then the layout beta_d per hand that puts both beams at the design angle. MODEL_ONLY."""
import json, math
import numpy as np
from scipy.optimize import least_squares
from model_core import *
isb0 = assign(30.0, 0); psi0 = layout(isb0, 1.10, 1.10); assert np.allclose(psi0, m.build(0, True)["psi"])
HF = (31.58, 35.94)            # HFSS BOTH sub-grid peaks (same estimator), from p1_c_pointing.py
(pa0, pb0, _, _) = beam_peaks(psi0, isb0, 1.10, 1.10); NOM = (pa0[0], pb0[0])
print("nominal model peaks (estimator bias included):", [round(v, 2) for v in NOM])
def res_fit(b):
    pa, pb, _, _ = beam_peaks(psi0, isb0, b[0], b[1]); return [pa[0]-HF[0], pb[0]-HF[1]]
sol = least_squares(res_fit, [1.15, 1.25], bounds=([1.0, 1.0], [1.32, 1.32]), diff_step=0.01, xtol=1e-4)
bA, bB = sol.x; print(f"JOINT FIT beta_t/k0: A {bA:.3f}  B {bB:.3f}   residual (deg) {np.round(sol.fun, 3)}")
# +-0.4 deg uncertainty corners
for dA, dB in ((0.4, 0.4), (-0.4, -0.4), (0.4, -0.4), (-0.4, 0.4)):
    t = (HF[0]+dA, HF[1]+dB)
    s2 = least_squares(lambda b: [beam_peaks(psi0, isb0, b[0], b[1])[0][0]-t[0], beam_peaks(psi0, isb0, b[0], b[1])[1][0]-t[1]], sol.x, bounds=([1.0, 1.0], [1.32, 1.32]), diff_step=0.01, xtol=1e-4)
    print(f"   target {t}: beta_t A {s2.x[0]:.3f} B {s2.x[1]:.3f}")
pa, pb, R, L = beam_peaks(psi0, isb0, bA, bB)
print(f"VERIFY: layout 1.10/1.10 + fitted beta_t -> RHCP th {pa[0]:.2f} ph {pa[1]:.1f} ; LHCP th {pb[0]:.2f} ph {pb[1]:.1f}  (HFSS 31.58 / 35.94)")
# design beta_d per hand so that, with the fitted physical beta_t, both beams sit at the nominal-model angle
def res_des(d):
    psi = layout(isb0, d[0], d[1]); pa, pb, _, _ = beam_peaks(psi, isb0, bA, bB); return [pa[0]-NOM[0], pb[0]-NOM[1]]
sd = least_squares(res_des, [bA, bB], bounds=([1.0, 1.0], [1.35, 1.35]), diff_step=0.01, xtol=1e-4)
dA, dB = sd.x; psiD = layout(isb0, dA, dB); pa, pb, R, L = beam_peaks(psiD, isb0, bA, bB)
print(f"DESIGN beta_d/k0 that puts both beams at the design angle (model, beta_t fixed at fit): A {dA:.3f}  B {dB:.3f}; residual {np.round(sd.fun,3)}")
print(f"   check: RHCP th {pa[0]:.2f} ph {pa[1]:.1f}  LHCP th {pb[0]:.2f} ph {pb[1]:.1f}")
json.dump(dict(HFSS_peaks_deg=HF, nominal_model=NOM, beta_t_fit=[bA, bB], beta_d_design=[dA, dB], target_design_deg=NOM), open("p1_c3_fit.json", "w"), indent=1)
