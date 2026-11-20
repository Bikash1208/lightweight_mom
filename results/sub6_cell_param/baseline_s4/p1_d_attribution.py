"""Phase 1 d (+ radial profile for a): lobe attribution by linear superposition of stored HFSS complex far fields.
S = ALL_PLAIN (structure), A = B_ALL_OFF_IDEAL - ALL_PLAIN (A cells), B = BOTH - B_ALL_OFF_IDEAL (B cells + interaction).
Linear superposition is an ASSUMPTION (cell-cell interaction terms folded into A and B). Labels: STORED_HFSS_2026R1 + MODEL_SUPERPOSITION."""
import csv, math, json
import numpy as np
from common import *
from sidequest_sll_def import analyze

def raw(s):
    d = load(s); a = 10**(-d["scale_db"]/20)
    return d["Et"]*a, d["Ep"]*a, d["scale_db"]
# scale check: scale_s should track 1/accepted-power only if FFD is incident-normalised
for s in STATES:
    sc = load(s)["scale_db"]; acc = load(s)["m"]["accepted_power_fraction"]
    print(f"{s:16s} gain-scale {sc:8.3f} dB ; 10log10(Pacc) {10*math.log10(acc):7.3f} dB ; scale+10log10(Pacc) = {sc+10*math.log10(acc):8.3f}")
# raw-field common reference: use ALL_PLAIN scale (gain referenced to incident power, constant C/Pacc folded)
ref_scale = load("BOTH")["scale_db"]
def cp(Et, Ep, state_conv="et_plus_j_ep"):
    # RHCP/LHCP projections, convention identical in all 7 states (checked)
    return (Et + 1j*Ep), (Et - 1j*Ep)          # RHCP, LHCP (conv et_plus_j_ep)
assert all(load(s)["rhcp_conv"] == "et_plus_j_ep" for s in STATES)
F = {}
for s in STATES:
    Et, Ep, sc = raw(s)
    # fields in 'gain units referenced to incident 1 W' using BOTH's constant: gain_ref = C*|E|^2, C=10^(scale_BOTH/10)*Pacc_BOTH
    k = 10**((load("BOTH")["scale_db"])/20)*math.sqrt(load("BOTH")["m"]["accepted_power_fraction"])
    F[s] = cp(Et*k, Ep*k)
th, ph = load("BOTH")["theta"], load("BOTH")["phi"]
S = F["ALL_PLAIN"]; Aoff = F["B_ALL_OFF_IDEAL"]; Bth = F["BOTH"]
comp = {}
for h, name in enumerate(("RHCP", "LHCP")):
    comp[name] = dict(S=S[h], A=Aoff[h]-S[h], B=Bth[h]-Aoff[h], Bleak=F["B_ALL_SHORT"][h]-Aoff[h], B24=F["A_ONLY_IDEAL"][h]-Aoff[h])
def gdb(x): return 10*np.log10(np.maximum(np.abs(x)**2, 1e-300)/2)   # |Et±jEp|^2 = 2*(...) ; dB in 'gain referenced to Pinc' (relative ok)
# calibrate to dBi of BOTH peak so numbers match: BOTH RHCP peak 9.03 dBi
cal = {}
for h, name in enumerate(("RHCP", "LHCP")):
    p = load("BOTH")["m"][name]; i = int(np.argmin(abs(th-p["peak_theta_deg"]))); j = int(np.argmin(abs(ph-p["peak_phi_deg"])))
    cal[name] = p["peak_gain_db"] - 10*math.log10(abs(Bth[h][i, j])**2)
print("cal", cal)
def gd(name, key, i, j): return 10*math.log10(abs(comp[name][key][i, j])**2) + cal[name]
rows = []
print("\nLobe attribution (levels in dBi, per hand; |.| of each complex component at the lobe direction)")
print(f"{'hand':5s} {'(th,ph)':>9s} {'BOTH':>6s} | {'|S|':>6s} {'|A|':>6s} {'|B|':>6s} | {'S+A':>6s} {'S+B':>6s} | {'|Bleak|':>7s} {'|B24|':>6s} | dominant")
for r in csv.DictReader(open(GATE/"v31_sll/lobes.csv")):
    if r["state"] != "BOTH" or r["kind"] in ("peak", "sll_excluding_44deg_regions"): continue
    name = r["hand"]; t, p = float(r["theta_deg"]), float(r["phi_deg"])
    i = int(np.argmin(abs(th-t))); j = int(np.argmin(abs(ph-p)))
    c = comp[name]; h = 0 if name == "RHCP" else 1
    vals = dict(S=gd(name,"S",i,j), A=gd(name,"A",i,j), B=gd(name,"B",i,j))
    both = 10*math.log10(abs(Bth[h][i,j])**2)+cal[name]
    sa = 10*math.log10(abs(c["S"][i,j]+c["A"][i,j])**2)+cal[name]; sb = 10*math.log10(abs(c["S"][i,j]+c["B"][i,j])**2)+cal[name]
    dom = max(vals, key=vals.get)
    print(f"{name:5s} ({t:3.0f},{p:3.0f}) {both:6.2f} | {vals['S']:6.2f} {vals['A']:6.2f} {vals['B']:6.2f} | {sa:6.2f} {sb:6.2f} | {gd(name,'Bleak',i,j):7.2f} {gd(name,'B24',i,j):6.2f} | {dom}  (stored {float(r['gain_dbi']):.2f})")
    rows.append(dict(hand=name, kind=r["kind"], theta=t, phi=p, both_dbi=both, S_dbi=vals["S"], A_dbi=vals["A"], B_dbi=vals["B"], SA_dbi=sa, SB_dbi=sb, dominant=dom))
json.dump(rows, open("p1_lobe_attribution.json", "w"), indent=1)

print("\nSLL of component sums (first-null mask, above horizon; same analyze() as the stored scoring), peak mask centred on BOTH peak")
for name, key in (("RHCP", "RHCP"), ("LHCP", "LHCP")):
    h = 0 if name == "RHCP" else 1; p = load("BOTH")["m"][name]
    pk = (p["peak_theta_deg"], p["peak_phi_deg"])
    for lab, fld in (("BOTH (S+A+B)", Bth[h]), ("S+A  (B off ideal)", Aoff[h]), ("S+B  (A removed)", comp[name]["S"]+comp[name]["B"]),
                     ("A only", comp[name]["A"]), ("B only", comp[name]["B"]), ("S only (plain)", comp[name]["S"])):
        pw = np.abs(fld)**2
        try:
            r = analyze(th, ph, pw, peak=pk)
            print(f"{name} {lab:20s} SLL {r['SLL_db']:7.2f} dB at {r['SLL_dir']}  (peak at BOTH peak dir: {10*math.log10(pw[int(np.argmin(abs(th-pk[0]))), int(np.argmin(abs(ph-pk[1])))])+cal[name]:6.2f} dBi)")
        except Exception as e: print(name, lab, "analyze failed", e)

print("\n== radial near-field profile (BOTH, 2 mm above copper): mean r*Ez^2 over 8 azimuths (surface-wave power-flux proxy), normalised to r=8..14 mm")
rp = {}
for r in csv.DictReader(open(GATE/"v31/field_diagnostic/radial_profiles.csv")):
    if r["state"] != "BOTH": continue
    try: rp.setdefault(int(r["r_mm"]), []).append(float(r["e_z_sq"]))
    except: pass
rs = np.array(sorted(rp)); ez = np.array([np.mean(rp[k]) for k in rs]); pf = rs*ez
ref = pf[(rs >= 8) & (rs <= 14)].mean()
for a, b in ((8,14),(20,30),(40,50),(60,70),(80,90),(100,110),(120,130),(130,140),(140,150)):
    m = (rs>=a)&(rs<=b)
    if m.any(): print(f"r {a:3d}-{b:3d} mm: r*Ez^2 = {10*math.log10(pf[m].mean()/ref):6.2f} dB vs r~11 mm")
