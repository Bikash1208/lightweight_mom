"""Phase 2 per-solve analysis. Usage: python p2_analyze.py TAG   (reads C:/work/s4_work/<TAG>/{extract,profile,status.json})
Reports: peak gain+direction per hand, SLL stored definition and fixed-cone SLL max(1.5 HPBW,10deg), theta cut through each beam peak,
m=0 ring level at theta 44-48, accepted/radiated fractions, memory, wall time, near-surface field band levels."""
import sys, json, math, csv
from pathlib import Path
import numpy as np
sys.path.insert(0, r"C:\work\lwmom_s4\scripts\research"); sys.dont_write_bytecode = True
from full_board_gate_v31_extract import ffd, direction_separation
from sidequest_sll_def import analyze

def sphere_w(theta, phi):
    t = np.deg2rad(theta); p = np.deg2rad(phi)
    def tw(a):
        w = np.zeros_like(a); w[1:-1] = (a[2:]-a[:-2])/2; w[0] = (a[1]-a[0])/2; w[-1] = (a[-1]-a[-2])/2; return w
    return (np.sin(t)*tw(t))[:, None]*tw(p)[None, :]

def run(tag, work=Path(r"C:\work\s4_work")):
    d = work/tag; m = json.loads((d/"extract"/"metrics.json").read_text()); st = json.loads((d/"status.json").read_text())
    th, ph, et, ep = ffd(Path(m["ffd_path"]))
    plus, minus = np.abs(et+1j*ep)**2, np.abs(et-1j*ep)**2
    rh, lh = (plus, minus) if m["ffd_rhcp_convention"] == "et_plus_j_ep" else (minus, plus)
    i = int(np.argmin(abs(th-m["RHCP"]["peak_theta_deg"]))); j = int(np.argmin(abs(ph-m["RHCP"]["peak_phi_deg"])))
    sc = m["RHCP"]["peak_gain_db"]-10*math.log10(rh[i, j]); k = 10**(sc/10)
    G = {"RHCP": rh*k, "LHCP": lh*k}
    res = {"tag": tag, "s11_db": m["s11_db"], "accepted_frac": m["accepted_power_fraction"], "wall_s": st.get("wall_s"),
           "hfss_max_memory_gb": (st.get("profile") or {}).get("hfss_max_memory_gb"), "adaptive_passes": len((st.get("profile") or {}).get("delta_s", [])),
           "tets_max": max((st.get("profile") or {}).get("tetrahedra", [0])), "converged": (st.get("profile") or {}).get("converged")}
    w = sphere_w(th, ph); Prh = float((G["RHCP"]*w).sum()/(4*math.pi)); Plh = float((G["LHCP"]*w).sum()/(4*math.pi))
    res.update(rad_over_acc=Prh+Plh, rh_over_acc=Prh, lh_over_acc=Plh, rad_over_inc=(Prh+Plh)*m["accepted_power_fraction"])
    # m=0 ring (phi-mean of complex Et,Ep), gain per hand at theta 44,46,48
    k_amp = math.sqrt(k)
    ring = {}
    for t in (44, 46, 48):
        ii = int(np.argmin(abs(th-t))); e0t = et[ii, :-1].mean(); e0p = ep[ii, :-1].mean()
        a, b = abs(e0t+1j*e0p)**2*k, abs(e0t-1j*e0p)**2*k
        r_, l_ = (a, b) if m["ffd_rhcp_convention"] == "et_plus_j_ep" else (b, a)
        ring[str(t)] = {"RHCP_m0_dbi": 10*math.log10(max(r_, 1e-300)), "LHCP_m0_dbi": 10*math.log10(max(l_, 1e-300))}
    res["ring_m0"] = ring
    res["ring_peak_dbi_any_phi"] = {h: float(10*np.log10(G[h][(th >= 40) & (th <= 52)].max())) for h in G}
    for hand in ("RHCP", "LHCP"):
        p = m[hand]; pk = (p["peak_theta_deg"], p["peak_phi_deg"]); pw = G[hand]
        an = analyze(th, ph, pw, peak=pk)
        hp = an["HPBW_mean_deg"]; cone = max(1.5*hp, 10.0)
        sep = direction_separation(th, ph, pk[0], pk[1]); vis = (th[:, None] <= 90)
        mm = np.where((sep > cone) & vis, pw, 0); kk = np.unravel_index(int(np.argmax(mm)), mm.shape)
        # theta cut through peak at the peak phi
        jj = int(np.argmin(abs(ph-pk[1]))); cut = {int(t): round(float(10*np.log10(max(pw[int(np.argmin(abs(th-t))), jj], 1e-300))), 2) for t in range(0, 62, 4)}
        res[hand] = dict(peak_dbi=p["peak_gain_db"], peak_theta=pk[0], peak_phi=pk[1], sll_stored_db=an["SLL_db"], sll_stored_dir=an["SLL_dir"],
                         hpbw_mean=hp, cone_deg=cone, sll_cone_db=float(10*np.log10(mm[kk]) - p["peak_gain_db"]), sll_cone_dir=[float(th[kk[0]]), float(ph[kk[1]])],
                         at_30_0=m[hand]["at_30_phi0_db"], at_30_180=m[hand]["at_30_phi180_db"], theta_cut_dbi=cut)
    # near-surface profile bands
    pr = d/"profile"/f"{tag.split('_',1)[1] if False else ''}"
    csvs = list((d/"profile").glob("*_radial.csv")); bands = {}
    if csvs:
        rows = list(csv.DictReader(open(csvs[0])))
        for a_, b_ in ((6, 29), (40, 100), (100, 119), (120, 145)):
            e = [float(r["e_total_sq"]) for r in rows if a_ <= int(r["r_mm"]) <= b_ and r["e_total_sq"] not in ("nan",) and math.isfinite(float(r["e_total_sq"]))]
            z = [float(r["e_z_sq"]) for r in rows if a_ <= int(r["r_mm"]) <= b_ and math.isfinite(float(r["e_z_sq"]))]
            bands[f"r{a_}-{b_}"] = {"E2_mean": float(np.mean(e)), "Ez2_mean": float(np.mean(z))}
    res["near_surface"] = bands
    (d/"p2_metrics.json").write_text(json.dumps(res, indent=1))
    return res

if __name__ == "__main__":
    r = run(sys.argv[1])
    for hand in ("RHCP", "LHCP"):
        h = r[hand]; print(f"{hand}: peak {h['peak_dbi']:.2f} dBi at ({h['peak_theta']:.0f},{h['peak_phi']:.0f}); SLL stored {h['sll_stored_db']:.2f} at {h['sll_stored_dir']}; cone({h['cone_deg']:.1f}) SLL {h['sll_cone_db']:.2f} at {h['sll_cone_dir']}; HPBW {h['hpbw_mean']:.1f}")
        print("   theta cut at peak phi:", h["theta_cut_dbi"])
    print("S11 %.2f dB, accepted %.2f%%, radiated/acc %.2f%% (RH %.2f LH %.2f), ring m0 %s, ring max(any phi,40-52) %s" % (r["s11_db"], r["accepted_frac"]*100, r["rad_over_acc"]*100, r["rh_over_acc"]*100, r["lh_over_acc"]*100, {k: (round(v['RHCP_m0_dbi'], 2), round(v['LHCP_m0_dbi'], 2)) for k, v in r["ring_m0"].items()}, {k: round(v, 2) for k, v in r["ring_peak_dbi_any_phi"].items()}))
    print("wall", r["wall_s"], "mem GB", r["hfss_max_memory_gb"], "tets", r["tets_max"], "converged", r["converged"], "passes", r["adaptive_passes"])
    print("near-surface", r["near_surface"])
