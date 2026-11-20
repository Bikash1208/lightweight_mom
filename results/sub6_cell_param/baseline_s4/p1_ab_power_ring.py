"""Phase 1 a+b: power budget per state and ring diagnostics, from stored FFD only (no AEDT started)."""
import json, math
import numpy as np
from common import *

out = {}
print("== (a) power budget  (Prad/Pacc from integrated far field; Pacc = 1-|S11|^2)")
print(f"{'state':16s} {'S11dB':>7s} {'Pacc%':>6s} {'Prad/Pacc%':>10s} {'RH%':>6s} {'LH%':>6s} {'th>90%':>7s} {'rest%Pacc':>9s} {'Prad/Pinc%':>10s} {'rest%Pinc':>9s}")
for s in STATES:
    d = load(s); m = d["m"]
    th, ph = d["theta"], d["phi"]; w = sphere_weights(th, ph)
    Prh = float((d["rh"]*w).sum()/(4*math.pi)); Plh = float((d["lh"]*w).sum()/(4*math.pi))
    up = (th <= 90)[:, None]
    Pdown = float(((d["rh"]+d["lh"])*w*(~up)).sum()/(4*math.pi))
    S11 = m["s11_db"]; acc = m["accepted_power_fraction"]
    eff = Prh + Plh
    out[s] = dict(s11_db=S11, accepted_frac=acc, rad_over_acc=eff, rh_over_acc=Prh, lh_over_acc=Plh, below_horizon_over_acc=Pdown,
                  rest_over_acc=1-eff, rad_over_inc=eff*acc, rest_over_inc=1-eff*acc, wsum=float(w.sum()/(4*math.pi)))
    print(f"{s:16s} {S11:7.2f} {acc*100:6.2f} {eff*100:10.2f} {Prh*100:6.2f} {Plh*100:6.2f} {Pdown*100:7.3f} {(1-eff)*100:9.2f} {eff*acc*100:10.2f} {(1-eff*acc)*100:9.2f}")
print("sphere weight sum/(4pi) =", out["BOTH"]["wsum"])

print("\n== (b) ALL_PLAIN ring diagnostics")
d = load("ALL_PLAIN"); th, ph = d["theta"], d["phi"]; P = ph[:-1]
Et, Ep = d["Et"][:, :-1], d["Ep"][:, :-1]
Gtot = (np.abs(Et)**2+np.abs(Ep)**2)*2   # == rh+lh (|Et+jEp|^2+|Et-jEp|^2 = 2(|Et|^2+|Ep|^2))
print("theta  Gtot_dBi(phi-avg power)  phi-std(dB)  |Ep|/|Et|  m0_frac_of_Et")
rows = []
for i in range(0, 91, 2):
    g = Gtot[i]; avg = g.mean()
    m0 = abs(Et[i].mean())**2 / max(np.mean(abs(Et[i])**2), 1e-300)
    rows.append((th[i], db(avg), np.std(db(g)), np.sqrt(np.mean(abs(Ep[i])**2)/max(np.mean(abs(Et[i])**2),1e-300)), m0))
for r in rows: print("%5.0f %8.2f %8.2f %9.3f %9.3f" % r)
# power by theta band, fraction of ALL_PLAIN radiated power
w = sphere_weights(th, ph)[:, :]; Gall = d["rh"]+d["lh"]; tot = (Gall*w).sum()
for a, b in ((0, 20), (20, 40), (40, 56), (56, 90), (90, 180)):
    msk = ((th >= a) & (th < b))[:, None]
    print(f"band {a:3d}-{b:3d} deg: {(Gall*w*msk).sum()/tot*100:6.2f}% of ALL_PLAIN radiated power")

print("\n== ring m=0 amplitude per state (phi-mean of Et), relative to ALL_PLAIN, theta 40-52")
ref = load("ALL_PLAIN")
for s in STATES:
    dd = load(s)
    e0 = dd["Et"][:, :-1].mean(axis=1); r0 = ref["Et"][:, :-1].mean(axis=1)
    sel = (th >= 40) & (th <= 52)
    ratio = np.sqrt(np.mean(abs(e0[sel])**2)/np.mean(abs(r0[sel])**2))
    dif = np.sqrt(np.mean(abs(e0[sel]-r0[sel])**2)/np.mean(abs(r0[sel])**2))
    print(f"{s:16s} |m0(Et)|/|m0 plain| = {20*math.log10(ratio):6.2f} dB ; |m0 - m0plain|/|m0plain| = {20*math.log10(max(dif,1e-9)):6.2f} dB")
json.dump(out, open("p1_power_budget.json", "w"), indent=1)
