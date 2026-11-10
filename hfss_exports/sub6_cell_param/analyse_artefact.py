"""Compare nominal modes across PEC floor depths and a 377-ohm floor."""
import csv
import itertools
import json
from pathlib import Path

import cell_param_pyaedt as cp

HERE=Path(__file__).resolve().parent
cp.set_profile("bar50")   # artefact test runs on the BAR50 calibration cell
DEPTH={"PEC_30mm":30.0,"PEC_60mm":60.0,"Z377_30mm":30.0}
RUNS=HERE/"artefact_runs_nm"
DIRS={"PEC_30mm":"PEC_30mm","PEC_60mm":"PEC_60mm","Z377_30mm":"Z377_30mm"}


def load(label,state,direction):
    p=RUNS/DIRS[label]/f"ARTEFACT_{label}_{state}_{direction}.json"
    r=json.loads(p.read_text())
    ok,why=cp.readback_ok(r)
    if not ok:
        raise SystemExit(f"ANALYSIS REFUSED {p.name}: {why}")
    # Intended floor depth for this case must be what AEDT actually had (the 2026-11-10 air_below bug).
    if abs(r["variables"]["air_below"]-DEPTH[label])>1e-9 or abs(r["readback_si"]["air_below"]*1e3-DEPTH[label])>1e-6:
        raise SystemExit(f"ANALYSIS REFUSED {p.name}: air_below {r['variables']['air_below']} / read back "
                         f"{r['readback_si']['air_below']*1e3} mm, intended {DEPTH[label]} mm")
    return r["modes"]


def match(reference,other):
    """Minimum total relative-frequency assignment; never treated as field branch proof."""
    best=None
    for perm in itertools.permutations(range(len(other)),len(reference)):
        cost=sum(abs(other[j]["f_ghz"]/reference[i]["f_ghz"]-1) for i,j in enumerate(perm))
        if best is None or cost<best[0]:
            best=(cost,perm)
    return [other[j] for j in best[1]]


def main():
    rows=[]
    for state in ("OFF","ON"):
        for direction in ("across","along"):
            p30=load("PEC_30mm",state,direction)
            p60=match(p30,load("PEC_60mm",state,direction))
            z377=match(p30,load("Z377_30mm",state,direction))
            for a,b,c in zip(p30,p60,z377):
                dp=100*abs(b["f_ghz"]/a["f_ghz"]-1)
                dz=100*abs(c["f_ghz"]/a["f_ghz"]-1)
                verdict="REAL" if max(dp,dz)<=1 else "ARTEFACT" if max(dp,dz)>=2 else "UNRESOLVED"
                rows.append({"state":state,"direction":direction,"pec30_mode":a["mode"],
                    "pec30_f_ghz":a["f_ghz"],"pec30_q":a["q"],
                    "pec60_mode":b["mode"],"pec60_f_ghz":b["f_ghz"],"pec60_q":b["q"],
                    "z377_mode":c["mode"],"z377_f_ghz":c["f_ghz"],"z377_q":c["q"],
                    "pec_depth_shift_pct":dp,"absorbing_shift_pct":dz,"verdict":verdict,
                    "branch_identity":"FREQUENCY_ASSIGNMENT_ONLY_EZ_REVIEW_REQUIRED"})
    out=HERE/"ARTEFACT_TEST.csv"
    with out.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    summary={"status":"COMPLETE_WITH_FIELD_REVIEW_REQUIRED",
       "rule":"REAL if PEC30-vs-PEC60 and PEC30-vs-Z377 matched f shifts <=1%; ARTEFACT if either >=2%; otherwise UNRESOLVED",
       "counts":{v:sum(r["verdict"]==v for r in rows) for v in ("REAL","ARTEFACT","UNRESOLVED")},
       "warning":"Assignments minimize total frequency shift and do not prove branch identity at crossings; Ez review is mandatory.",
       "rows":len(rows)}
    (HERE/"ARTEFACT_TEST.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    main()
