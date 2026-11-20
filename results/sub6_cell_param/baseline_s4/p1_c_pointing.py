"""Phase 1 c: effective beta per hand from the HFSS beam pointing, using the layout phase formula + af_jones-style model.
MODEL_ONLY for the model parts. Same peak estimator (3x3 quadratic on a 2-deg grid, log-power) applied to HFSS and model."""
import sys, math, json
import numpy as np
sys.path.insert(0, r"C:\work\lwmom_s4\results\sub6_cell_param\idc_cells\full_board_layout")
import make_layout_v3 as m
from common import *

k0 = m.k0; P = m.P
TH2 = np.arange(0, 91, 2.0); PH2 = np.arange(0, 360, 2.0)

def est_peak(th, ph, pw, near):
    """sub-grid peak (theta,phi) near 'near'=(theta,phi) by quadratic fit of dB power in direction cosines (u,v) over a 5x5 window."""
    i0 = int(np.argmin(abs(th-near[0]))); j0 = int(np.argmin(abs(ph-near[1])))
    # refine to local max in 3x3
    best = None
    for di in range(-2, 3):
        for dj in range(-2, 3):
            i = i0+di; j = (j0+dj) % len(ph)
            if 0 <= i < len(th) and (best is None or pw[i, j] > pw[best]): best = (i, j)
    i0, j0 = best
    U = []; V = []; Y = []
    for di in range(-2, 3):
        for dj in range(-2, 3):
            i = i0+di; j = (j0+dj) % len(ph)
            if not (0 <= i < len(th)): continue
            t = math.radians(th[i]); p = math.radians(ph[j0]+dj*(ph[1]-ph[0]))
            U.append(math.sin(t)*math.cos(p)); V.append(math.sin(t)*math.sin(p)); Y.append(10*math.log10(max(pw[i, j], 1e-300)))
    U = np.array(U); V = np.array(V); Y = np.array(Y)
    A = np.column_stack([np.ones_like(U), U, V, U*U, U*V, V*V])
    c = np.linalg.lstsq(A, Y, rcond=None)[0]
    H = np.array([[2*c[3], c[4]], [c[4], 2*c[5]]]); g = -np.array([c[1], c[2]])
    u, v = np.linalg.solve(H, g)
    th_pk = math.degrees(math.asin(min(1, math.hypot(u, v)))); ph_pk = math.degrees(math.atan2(v, u)) % 360
    return th_pk, ph_pk, u, v

print("== HFSS beam positions (sub-grid, same estimator)")
hf = {}
for s in ("BOTH", "B_ALL_OFF_IDEAL", "A_ONLY_IDEAL", "B_ALL_SHORT"):
    d = load(s)
    for hand, pw, near in (("RHCP", d["rh"], (32, 0)), ("LHCP", d["lh"], (36, 180))):
        r = est_peak(d["theta"], d["phi"], pw, near)
        hf[(s, hand)] = r
        print(f"{s:16s} {hand}: grid-peak -> sub-grid (theta {r[0]:6.2f}, phi {r[1]:6.1f}) u={r[2]:.4f} v={r[3]:.4f}")

# ---- model (af_jones physics): radial launch rho_hat*exp(-j beta_t r)*amp; each cell half-wave plate R(psi)diag(1,-1)R(-psi)
d0 = m.build(0, True)
x, y, r, phi, psi0, isb = d0["x"], d0["y"], d0["r"], d0["phi"], d0["psi"], d0["isb"]
amp = np.exp(-np.interp(r, [0, 75, 150], [1.11, 1.45, 1.83])*r*1e-3)
eR = np.array([1, -1j])/math.sqrt(2); eL = np.array([1, 1j])/math.sqrt(2)
TT, PP = np.meshgrid(np.deg2rad(TH2), np.deg2rad(PH2), indexing="ij")
DX = np.sin(TT)*np.cos(PP); DY = np.sin(TT)*np.sin(PP)

def layout_psi(beta_dA, beta_dB, az=True):
    PhiA = -k0*math.sin(m.TH)*x + beta_dA*r
    PhiB = +k0*math.sin(m.TH)*x + beta_dB*r
    return 0.5*(np.where(isb, -PhiB, PhiA) + (phi if az else 0))

def model_far(psi, beta_t, amp_=None):
    a = amp if amp_ is None else amp_
    c, s = np.cos(2*psi), np.sin(2*psi)
    Ox = c*np.cos(phi)+s*np.sin(phi); Oy = s*np.cos(phi)-c*np.sin(phi)
    w = a*np.exp(-1j*beta_t*r)
    Fx = np.zeros(TT.shape, complex); Fy = np.zeros(TT.shape, complex)
    for k in range(len(x)):
        ph_ = np.exp(1j*k0*(DX*x[k]+DY*y[k]))*w[k]
        Fx += ph_*Ox[k]; Fy += ph_*Oy[k]
    R = np.conj(eR[0])*Fx+np.conj(eR[1])*Fy; L = np.conj(eL[0])*Fx+np.conj(eL[1])*Fy
    return np.abs(R)**2, np.abs(L)**2

def peaks(psi, beta_t):
    R, L = model_far(psi, beta_t)
    return est_peak(TH2, PH2, R, (30, 0)), est_peak(TH2, PH2, L, (30, 180)), R, L

psi_ref = layout_psi(1.10*k0, 1.10*k0)
assert np.allclose(psi_ref, psi0)
print("\n== model, layout beta_d=1.10k0, scan actual beta_t/k0 (peak theta of RHCP / LHCP, sub-grid)")
scan = {}
for bt in np.arange(0.95, 1.40001, 0.025):
    (a, b, R, L) = peaks(psi_ref, bt*k0)[:4] if False else (None,)*4
    pa, pb, R, L = peaks(psi_ref, bt*k0)
    scan[round(bt, 3)] = (pa, pb)
    print(f"beta_t/k0={bt:5.3f}: RHCP peak (th {pa[0]:5.2f}, ph {pa[1]:6.1f})   LHCP peak (th {pb[0]:5.2f}, ph {pb[1]:6.1f})")
json.dump({str(k): [list(v[0]), list(v[1])] for k, v in scan.items()}, open("p1_c_scan.json", "w"))
