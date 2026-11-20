"""Phase 1 e: lever ranking support. (1) HFSS-derived: SLL if the structure/ring component S is scaled (coherent superposition of stored fields).
(2) MODEL_ONLY: 16 gate discs removed vs anchors; pointing error cost; with model_core."""
import math, json
import numpy as np
from common import *
from sidequest_sll_def import analyze
assert all(load(s)["rhcp_conv"] == "et_plus_j_ep" for s in STATES)
k = 10**(load("BOTH")["scale_db"]/20)*math.sqrt(load("BOTH")["m"]["accepted_power_fraction"])
def fld(s):
    d = load(s); a = k*10**(-d["scale_db"]/20); Et, Ep = d["Et"]*a, d["Ep"]*a; return (Et+1j*Ep, Et-1j*Ep)
th, ph = load("BOTH")["theta"], load("BOTH")["phi"]
S, Ao, Bt = fld("ALL_PLAIN"), fld("B_ALL_OFF_IDEAL"), fld("BOTH")
cal = 0.0926
print("== HFSS-derived (stored fields, linear superposition): structure component S scaled by g; A+B kept. SLL = first-null definition at the BOTH beam peaks")
res = {}
for g_db in (0, -3, -6, -10, -20, "removed"):
    g = 0.0 if g_db == "removed" else 10**(g_db/20)
    row = []
    for h, name, pk in ((0, "RHCP", (32, 0)), (1, "LHCP", (36, 180))):
        f = (Bt[h]-S[h]) + g*S[h]; pw = np.abs(f)**2
        sep = direction_separation(th, ph, pk[0], pk[1]); reg = np.where((sep <= 10) & (th[:, None] <= 90), pw, 0)
        i, j = np.unravel_index(int(np.argmax(reg)), reg.shape)
        r = analyze(th, ph, pw, peak=(float(th[i]), float(ph[j])))
        row.append((r["SLL_db"], r["SLL_dir"], 10*math.log10(pw[i, j])+cal))
    res[str(g_db)] = row
    print(f"S {str(g_db):>8s} dB: RHCP SLL {row[0][0]:6.2f} at {row[0][1]} (peak {row[0][2]:.2f} dBi) | LHCP SLL {row[1][0]:6.2f} at {row[1][1]} (peak {row[1][2]:.2f} dBi)")
json.dump(res, open("p1_e_ring_scaling.json", "w"), indent=1)

from model_core import *
isb0 = assign(30.0, 0)
fit = json.load(open("p1_c3_fit.json")); bA, bB = fit["beta_t_fit"]; dA, dB = fit["beta_d_design"]
def sll_of(psi, btA, btB, amp=None, label=""):
    Fx, Fy = far(psi, isb0, btA, btB, amp); R, L = cp_power(Fx, Fy); out = []
    for pw, near in ((R, (30, 0)), (L, (30, 180))):
        pk = est_peak(TH2, PH2, pw, near); i = int(round(pk[0]/2)); j = int(round(pk[1]/2))
        out.append(analyze(TH2, PH2, pw, peak=(TH2[i], PH2[j % 180 if j < 181 else 0]))["SLL_db"])
    print(f"{label:60s} SLL R {out[0]:6.2f} L {out[1]:6.2f}"); return out
print("\n== MODEL_ONLY (cell hologram only; no feed ring, no absorber)")
sll_of(layout(isb0, 1.10, 1.10), 1.10, 1.10, label="nominal (beta_t=beta_d=1.10)")
sll_of(layout(isb0, 1.10, 1.10), bA, bB, label="as built: beta_d 1.10, beta_t fitted (pointing error)")
sll_of(layout(isb0, dA, dB), bA, bB, label="retargeted beta_d per hand, beta_t fitted")
gate = (R_ >= 5.8) & (R_ < 5.8+2*P)
print("n gate-region cells (rho 5.8..29 mm):", int(gate.sum()))
a2 = AMP.copy(); a2[gate] = 0
sll_of(layout(isb0, 1.10, 1.10), 1.10, 1.10, amp=a2, label="nominal, 16-ish near-feed cells plain (gate discs, no radiation)")
sll_of(layout(isb0, dA, dB), bA, bB, amp=a2, label="retargeted, near-feed cells plain")
