"""MODEL_ONLY af_jones-physics core with per-hand beta_d (layout) and beta_t (physical effective), dither width, seed. Reused by Phase 2 generator."""
import sys, math
import numpy as np
sys.path.insert(0, r"C:\work\lwmom_s4\results\sub6_cell_param\idc_cells\full_board_layout")
import make_layout_v3 as m
k0 = m.k0; P = m.P
TH2 = np.arange(0, 91, 2.0); PH2 = np.arange(0, 361, 2.0)
TT, PP = np.meshgrid(np.deg2rad(TH2), np.deg2rad(PH2), indexing="ij"); DX = np.sin(TT)*np.cos(PP); DY = np.sin(TT)*np.sin(PP)
ix_, iy_ = np.meshgrid(np.arange(-12, 13), np.arange(-12, 13)); ix_ = ix_.ravel(); iy_ = iy_.ravel()
X_, Y_ = ix_*P, iy_*P; keep = (X_**2+Y_**2 <= m.R**2) & ((ix_ != 0) | (iy_ != 0)); IX, IY, X, Y = ix_[keep], iy_[keep], X_[keep], Y_[keep]
R_ = np.hypot(X, Y); PHI = np.arctan2(Y, X); AMP = np.exp(-np.interp(R_, [0, 75, 150], [1.11, 1.45, 1.83])*R_*1e-3)
eR = np.array([1, -1j])/math.sqrt(2); eL = np.array([1, 1j])/math.sqrt(2)

def assign(W, seed):
    pb = 0.5*(1-np.tanh(X/W)) if W > 0 else (X < 0).astype(float)
    isb = np.zeros(len(X), bool); rng = np.random.default_rng(seed)
    for c in np.unique(IX):
        ids = np.flatnonzero(IX == c); e = pb[ids].sum(); n = int(math.floor(e)+(rng.random() < e % 1))
        isb[rng.choice(ids, size=n, replace=False)] = True
    return isb

def layout(isb, bdA, bdB, thA=30.0, thB=30.0, az=True):
    """rotation psi (rad) = 0.5*(Phi + az*phi); same formula as make_layout_v3 but beta per hand (units of k0) and theta per hand."""
    PhiA = -k0*math.sin(math.radians(thA))*X + bdA*k0*R_
    PhiB = +k0*math.sin(math.radians(thB))*X + bdB*k0*R_
    return 0.5*(np.where(isb, -PhiB, PhiA) + (PHI if az else 0))

def far(psi, isb, btA, btB, amp=None):
    a = AMP if amp is None else amp
    bt = np.where(isb, btB, btA)*k0
    c, s = np.cos(2*psi), np.sin(2*psi)
    Ox = c*np.cos(PHI)+s*np.sin(PHI); Oy = s*np.cos(PHI)-c*np.sin(PHI); w = a*np.exp(-1j*bt*R_)
    Fx = np.zeros(TT.shape, complex); Fy = np.zeros(TT.shape, complex)
    for k in range(len(X)):
        e = np.exp(1j*k0*(DX*X[k]+DY*Y[k]))*w[k]; Fx += e*Ox[k]; Fy += e*Oy[k]
    return Fx, Fy        # complex Cartesian transverse field

def cp_power(Fx, Fy):
    return np.abs(np.conj(eR[0])*Fx+np.conj(eR[1])*Fy)**2, np.abs(np.conj(eL[0])*Fx+np.conj(eL[1])*Fy)**2

def est_peak(th, ph, pw, near):
    i0 = int(np.argmin(abs(th-near[0]))); j0 = int(np.argmin(abs(ph-near[1]))); best = None
    for di in range(-2, 3):
        for dj in range(-2, 3):
            i = i0+di; j = (j0+dj) % (len(ph)-1 if ph[-1] >= 360 else len(ph))
            if 0 <= i < len(th) and (best is None or pw[i, j] > pw[best]): best = (i, j)
    i0, j0 = best; U = []; V = []; Y = []; n = (len(ph)-1 if ph[-1] >= 360 else len(ph)); dph = ph[1]-ph[0]
    for di in range(-2, 3):
        for dj in range(-2, 3):
            i = i0+di
            if not (0 <= i < len(th)): continue
            j = (j0+dj) % n; t = math.radians(th[i]); p = math.radians(ph[j0]+dj*dph)
            U.append(math.sin(t)*math.cos(p)); V.append(math.sin(t)*math.sin(p)); Y.append(10*math.log10(max(pw[i, j], 1e-300)))
    U = np.array(U); V = np.array(V); Y = np.array(Y)
    A = np.column_stack([np.ones_like(U), U, V, U*U, U*V, V*V]); c = np.linalg.lstsq(A, Y, rcond=None)[0]
    u, v = np.linalg.solve(np.array([[2*c[3], c[4]], [c[4], 2*c[5]]]), -np.array([c[1], c[2]]))
    return math.degrees(math.asin(min(1, math.hypot(u, v)))), math.degrees(math.atan2(v, u)) % 360, u, v

def beam_peaks(psi, isb, btA, btB):
    Fx, Fy = far(psi, isb, btA, btB); R, L = cp_power(Fx, Fy)
    return est_peak(TH2, PH2, R, (30, 0)), est_peak(TH2, PH2, L, (30, 180)), R, L
