"""Phase 1 c (fit + verification) and e-support (MODEL_ONLY lever screen). Run after p1_c_pointing.py definitions (imported).
Model = af_jones physics; per-hand beta_t is an EFFECTIVE parameter (A cells use beta_tA, B cells beta_tB)."""
import io, contextlib, math, json, sys
import numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    import p1_c_pointing as pc          # reuses est_peak etc. (re-runs its scan quietly)
from p1_c_pointing import *
from scipy.optimize import brentq
from sidequest_sll_def import analyze

PH2 = np.arange(0, 361, 2.0); pc.PH2 = PH2
TT, PP = np.meshgrid(np.deg2rad(TH2), np.deg2rad(PH2), indexing="ij"); DX = np.sin(TT)*np.cos(PP); DY = np.sin(TT)*np.sin(PP)
ix_, iy_ = np.meshgrid(np.arange(-12, 13), np.arange(-12, 13)); ix_ = ix_.ravel(); iy_ = iy_.ravel()
X_, Y_ = ix_*P, iy_*P; keep = (X_**2+Y_**2 <= m.R**2) & ((ix_ != 0) | (iy_ != 0)); IX, IY, X, Y = ix_[keep], iy_[keep], X_[keep], Y_[keep]
R_ = np.hypot(X, Y); PHI = np.arctan2(Y, X); AMP = np.exp(-np.interp(R_, [0, 75, 150], [1.11, 1.45, 1.83])*R_*1e-3)

def assign(W, seed):
    pb = 0.5*(1-np.tanh(X/W)) if W > 0 else (X < 0).astype(float)*1.0
    isb = np.zeros(len(X), bool); rng = np.random.default_rng(seed)
    for c in np.unique(IX):
        ids = np.flatnonzero(IX == c); e = pb[ids].sum(); n = int(math.floor(e)+(rng.random() < e % 1))
        isb[rng.choice(ids, size=n, replace=False)] = True
    return isb

def layout(isb, bdA, bdB, thA=30.0, thB=30.0, az=True):
    PhiA = -k0*math.sin(math.radians(thA))*X + bdA*k0*R_
    PhiB = +k0*math.sin(math.radians(thB))*X + bdB*k0*R_
    return 0.5*(np.where(isb, -PhiB, PhiA) + (PHI if az else 0))

def far(psi, isb, btA, btB):
    bt = np.where(isb, btB, btA)*k0
    c, s = np.cos(2*psi), np.sin(2*psi)
    Ox = c*np.cos(PHI)+s*np.sin(PHI); Oy = s*np.cos(PHI)-c*np.sin(PHI)
    w = AMP*np.exp(-1j*bt*R_)
    Fx = np.zeros(TT.shape, complex); Fy = np.zeros(TT.shape, complex)
    for k in range(len(X)):
        e = np.exp(1j*k0*(DX*X[k]+DY*Y[k]))*w[k]; Fx += e*Ox[k]; Fy += e*Oy[k]
    return np.abs((Fx-1j*Fy)/math.sqrt(2))**2, np.abs((Fx+1j*Fy)/math.sqrt(2))**2   # RHCP, LHCP (eR=(1,-j): conj -> (Fx + j Fy)?) checked below

# convention check against af_jones: with beta_t=beta_d the designed beam hand/positions must come out as RHCP@(30,0), LHCP@(30,180)
isb0 = assign(30.0, 0); assert np.array_equal(isb0, m.build(0, True)["isb"])
psi0 = layout(isb0, 1.10, 1.10); assert np.allclose(psi0, m.build(0, True)["psi"])
def far_cp(psi, isb, btA, btB):
    bt = np.where(isb, btB, btA)*k0
    c, s = np.cos(2*psi), np.sin(2*psi)
    Ox = c*np.cos(PHI)+s*np.sin(PHI); Oy = s*np.cos(PHI)-c*np.sin(PHI); w = AMP*np.exp(-1j*bt*R_)
    Fx = np.zeros(TT.shape, complex); Fy = np.zeros(TT.shape, complex)
    for k in range(len(X)):
        e = np.exp(1j*k0*(DX*X[k]+DY*Y[k]))*w[k]; Fx += e*Ox[k]; Fy += e*Oy[k]
    eR = np.array([1, -1j])/math.sqrt(2); eL = np.array([1, 1j])/math.sqrt(2)
    return np.abs(np.conj(eR[0])*Fx+np.conj(eR[1])*Fy)**2, np.abs(np.conj(eL[0])*Fx+np.conj(eL[1])*Fy)**2
far = far_cp
R0, L0 = far(psi0, isb0, 1.10, 1.10)
pa = est_peak(TH2, PH2, R0, (30, 0)); pb = est_peak(TH2, PH2, L0, (30, 180))
print(f"nominal model (beta_t=beta_d=1.10): RHCP peak th {pa[0]:.2f} ph {pa[1]:.1f}; LHCP peak th {pb[0]:.2f} ph {pb[1]:.1f}  (design 30/30)")

def sll(Rp, Lp):
    out = {}
    for name, pw, near in (("R", Rp, (30, 0)), ("L", Lp, (30, 180))):
        pk = est_peak(TH2, PH2, pw, near)
        i = int(np.argmin(abs(TH2-near[0])));
        # integer grid peak for analyze()
        ii, jj = np.unravel_index(int(np.argmax(np.where((TH2[:, None] > 6), pw, 0)*((np.abs(((PH2[None, :]-near[1]+180) % 360)-180) < 30)))), pw.shape)
        r = analyze(TH2, PH2, pw, peak=(TH2[ii], PH2[jj]))
        out[name] = (r["SLL_db"], TH2[ii], PH2[jj])
    return out
s0 = sll(R0, L0)
print(f"nominal model SLL: RHCP {s0['R'][0]:.2f} dB, LHCP {s0['L'][0]:.2f} dB   (HFSS BOTH: -8.43 / -9.52)")

# --- fit effective beta per hand: A cells' beta_tA from RHCP peak theta (31.58), B cells' beta_tB from LHCP peak (35.94)
tgt_A, tgt_B = 31.58, 35.94
fA = lambda b: est_peak(TH2, PH2, far(psi0, isb0, b, 1.10)[0], (30, 0))[0]-tgt_A
fB = lambda b: est_peak(TH2, PH2, far(psi0, isb0, 1.10, b)[1], (30, 180))[0]-tgt_B
bA = brentq(fA, 1.10, 1.30, xtol=1e-3); bB = brentq(fB, 1.10, 1.315, xtol=1e-3)
print(f"FIT: beta_effA/k0 = {bA:.3f}, beta_effB/k0 = {bB:.3f}")
for tol in (0.4, ):   # +-0.4 deg pointing-estimator uncertainty -> beta interval
    try:
        print(f"  A interval for +-{tol} deg:", round(brentq(lambda b: fA(b)+tol, 1.10, 1.32, xtol=1e-3), 3), round(brentq(lambda b: fA(b)-tol, 1.10, 1.32, xtol=1e-3), 3))
        print(f"  B interval for +-{tol} deg:", round(brentq(lambda b: fB(b)+tol, 1.10, 1.315, xtol=1e-3), 3), round(brentq(lambda b: fB(b)-tol, 1.10, 1.32, xtol=1e-3), 3))
    except Exception as e: print("  interval failed", e)
Rf, Lf = far(psi0, isb0, bA, bB)
pa = est_peak(TH2, PH2, Rf, (30, 0)); pb = est_peak(TH2, PH2, Lf, (30, 180))
print(f"VERIFY (layout beta_d=1.10, model with fitted beta_t): RHCP th {pa[0]:.2f} ph {pa[1]:.1f}  LHCP th {pb[0]:.2f} ph {pb[1]:.1f}   (HFSS 31.58 / 35.94)")
sf = sll(Rf, Lf); print(f"   model SLL with the pointing error: RHCP {sf['R'][0]:.2f} LHCP {sf['L'][0]:.2f}")

# --- retarget: layout beta_dA=bA, beta_dB=bB, same physical beta_t
psiR = layout(isb0, bA, bB)
Rr, Lr = far(psiR, isb0, bA, bB)
pa = est_peak(TH2, PH2, Rr, (30, 0)); pb = est_peak(TH2, PH2, Lr, (30, 180)); sr = sll(Rr, Lr)
print(f"RETARGET (layout beta_dA={bA:.3f}, beta_dB={bB:.3f}; model beta_t = same): RHCP th {pa[0]:.2f}, LHCP th {pb[0]:.2f}; SLL R {sr['R'][0]:.2f} L {sr['L'][0]:.2f}")
json.dump(dict(betaA=bA, betaB=bB), open("p1_c_fit.json", "w"))

print("\n== MODEL_ONLY lever screen at nominal pointing (beta_t=beta_d=1.10 both): SLL R/L (dB)")
print("dither width W (mm), seed 0:")
for W in (0, 10, 20, 30, 45, 60, 90):
    isb = assign(W, 0) if W > 0 else assign(0, 0); psi = layout(isb, 1.10, 1.10); Rp, Lp = far(psi, isb, 1.10, 1.10); s = sll(Rp, Lp)
    print(f"  W={W:3d}: nB={isb.sum():3d}  RHCP {s['R'][0]:6.2f}  LHCP {s['L'][0]:6.2f}")
print("seed scan, W=30:")
sl = []
for sd in range(8):
    isb = assign(30.0, sd); psi = layout(isb, 1.10, 1.10); Rp, Lp = far(psi, isb, 1.10, 1.10); s = sll(Rp, Lp); sl.append((s['R'][0], s['L'][0]))
    print(f"  seed {sd}: RHCP {s['R'][0]:6.2f} LHCP {s['L'][0]:6.2f}")
sl = np.array(sl); print(f"  mean R {sl[:,0].mean():.2f} (min {sl[:,0].min():.2f}) L {sl[:,1].mean():.2f} (min {sl[:,1].min():.2f})")
