"""Conservative C4X raw-mode inventory and automated ON/OFF screen."""

import csv
import json
import math
import re
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[2]
R = ROOT / "results" / "sub6_meander_model" / "supercell_hfss"
C0, P, PS = 299792458.0, .0116, .0928


def read(cell, phase):
    return json.loads((R / f"SUPER_C4X_{cell}_ph{phase}_PML40_{tag(cell,phase)}_modes.json").read_text())["modes"]


def tag(cell,phase):
    if cell=="ON" and phase==15:
        return "c4xrun"
    second=R/f"SUPER_C4X_{cell}_ph{phase}_PML40_c4xauto2_modes.json"
    return "c4xauto2" if second.exists() else "c4xauto1"


def field(cell, phase, mode):
    path = R / f"SUPER_C4X_{cell}_ph{phase}_PML40_{tag(cell,phase)}_mode{mode}_MagE.fld"
    if not path.exists():
        return None
    try:
        samples = np.loadtxt(path, skiprows=2)[:, 3].reshape(8, 2)
        vec = np.sqrt(np.mean(samples**2, axis=1))
        return vec / np.linalg.norm(vec) if np.linalg.norm(vec) else None
    except (ValueError, OSError):
        return None


def cosine(a, b):
    return float(a @ b) if a is not None and b is not None else None


def unit(cell, direction, phase):
    return json.loads((R / f"C4X_{cell}_PML_{direction}_{phase}_air40.035_c4xbase_modes.json").read_text())["modes"]


def baseline(cell, phase, f, fixed_branch=None):
    phase_cell = (720-phase)/8
    candidates = []
    for index in range(4):
        fs, qs = [], []
        for direction in ("across", "along"):
            a, b = unit(cell, direction, 85)[index], unit(cell, direction, 95)[index]
            t = (phase_cell-85)/10
            fs.append(a["f_ghz"]+t*(b["f_ghz"]-a["f_ghz"]))
            qs.append(a["Q"]+t*(b["Q"]-a["Q"]))
        candidates.append((abs(sum(fs)/2-f), index+1, sum(qs)/2))
    _, branch, qbase = (next(x for x in candidates if x[1]==fixed_branch)
                        if fixed_branch else min(candidates))
    return branch, qbase


def tets(cell, phase):
    path = R / f"SUPER_C4X_{cell}_ph{phase}_PML40_{tag(cell,phase)}_convergence.conv"
    matches = re.findall(r"^\s*\d+\|\s*(\d+)\|", path.read_text(), re.M)
    return int(matches[-1]) if matches else None


def analyze(cell, phases):
    raw = {p: read(cell, p) for p in phases}
    groups = {p: [m for m in raw[p] if 5.3 <= m["f_ghz"] <= 6.5 and m["Q"] < 1e5]
              for p in phases}
    links = {}
    for pa, pb in zip(phases[:-1], phases[1:]):
        a, b = groups[pa], groups[pb]
        if not a or not b:
            continue
        cost = np.zeros((len(a), len(b)))
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                c = cosine(field(cell, pa, x["mode"]), field(cell, pb, y["mode"]))
                cost[i, j] = (1-c if c is not None else 2) + .3*((x["f_ghz"]-y["f_ghz"])/.15)**2
        ia, ib = linear_sum_assignment(cost)
        links[(pa,pb)] = {a[i]["mode"]: b[j]["mode"] for i,j in zip(ia,ib)}
    tracks = []
    if len(phases) > 1:
        for start in groups[phases[0]]:
            ids = [start["mode"]]
            for pa,pb in zip(phases[:-1],phases[1:]):
                nxt = links.get((pa,pb),{}).get(ids[-1])
                if nxt is None:
                    break
                ids.append(nxt)
            if len(ids) == len(phases):
                tracks.append(ids)
    owner = {(p,mid): (index,ids) for index,ids in enumerate(tracks,1)
             for p,mid in zip(phases,ids)}
    track_branch={index:baseline(cell,phases[0],raw[phases[0]][ids[0]-1]["f_ghz"])[0]
                  for index,ids in enumerate(tracks,1)}
    track_slope={index:abs((raw[phases[-1]][ids[-1]-1]["f_ghz"]-
                            raw[phases[0]][ids[0]-1]["f_ghz"])/
                           (phases[-1]-phases[0])*8)
                 for index,ids in enumerate(tracks,1)}
    rows = []
    for p in phases:
        for m in raw[p]:
            mode,f,q = m["mode"],m["f_ghz"],m["Q"]
            track = owner.get((p,mode))
            beta_m=(720-p)/8*C0/(360*P*f*1e9)
            beta_p=(720+p)/8*C0/(360*P*f*1e9)
            angle=lambda beta: math.degrees(math.asin(beta-C0/(PS*f*1e9))) if abs(beta-C0/(PS*f*1e9))<=1 else None
            others=[abs(x["f_ghz"]-f) for x in raw[p]
                    if x["mode"]!=mode and 60<x["Q"]<1e5]
            row=dict(cell=cell,phase_deg=p,mode=mode,f_ghz=f,Q=q,max_solved_tets=tets(cell,p),
                     beta_minus_over_k0=beta_m,beta_plus_over_k0=beta_p,
                     harmonic_minus_angle_deg=angle(beta_m),harmonic_plus_angle_deg=angle(beta_p),
                     nearest_other_mode_ghz=min(others) if others else None,
                     track_id=track[0] if track else None,matched_modes="->".join(map(str,track[1])) if track else None,
                     track_slope_ghz_per_cell_deg=track_slope[track[0]] if track else None,
                     track_class=("SW" if track_slope[track[0]]>.02 else
                                  "RES" if track_slope[track[0]]<.005 else "HYBRID") if track else None,
                     field_cosine=None,vg_m_s=None,Q_base=None,Q_base_unit_branch=None,alpha_np_m=None,
                     gate_candidate=False)
            if track:
                idx=phases.index(p)
                ids=track[1]
                neighbor=idx+1 if idx<len(phases)-1 else idx-1
                pn=phases[neighbor]
                m2=raw[pn][ids[neighbor]-1]
                row["field_cosine"]=cosine(field(cell,p,mode),field(cell,pn,ids[neighbor]))
                # Adjacent local phase slope; for middle phase average both neighbours.
                if len(phases)==3 and idx==1:
                    m0=raw[phases[0]][ids[0]-1]; m3=raw[phases[2]][ids[2]-1]
                    slope=(m3["f_ghz"]-m0["f_ghz"])/(phases[2]-phases[0])
                else:
                    slope=(m2["f_ghz"]-f)/(pn-p)
                vg=abs(slope*1e9)*360*PS
                branch,qbase=baseline(cell,p,f,track_branch[track[0]])
                row.update(vg_m_s=vg,Q_base=qbase,Q_base_unit_branch=branch)
                if vg>0:
                    row["alpha_np_m"]=math.pi*f*1e9*(1/q-1/qbase)/vg
                row["gate_candidate"]=(cell=="ON" and row["track_class"]=="SW" and q>60 and q<1e5 and
                    vg>50e6 and row["field_cosine"] is not None and row["field_cosine"]>=.8 and
                    row["nearest_other_mode_ghz"]>=.05 and row["alpha_np_m"]>=1.2)
            rows.append(row)
    return rows,tracks


def run(cells=("ON",)):
    all_rows=[]; tracks={}
    for cell in cells:
        phases=(15,30,45) if cell=="ON" else (15,30)
        rows,track=analyze(cell,phases)
        all_rows+=rows; tracks[cell]=track
    out=R/"C4X_SUPERCELL1_rows.csv"
    with out.open("w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh,fieldnames=list(all_rows[0]))
        w.writeheader();w.writerows(all_rows)
    candidates=[r for r in all_rows if r["gate_candidate"]]
    summary={"on_gate_candidates":candidates,"on_pass":bool(candidates),"tracks":tracks,
             "mode_rows":len(all_rows),"method":"16-point |E| cosine plus frequency Hungarian pairing; kinematic folds shown for all raw modes"}
    (R/"C4X_SUPERCELL1_analysis.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    return summary


if __name__=="__main__":
    import sys
    print(json.dumps(run(("ON","OFF") if "--off" in sys.argv else ("ON",)),indent=2))
