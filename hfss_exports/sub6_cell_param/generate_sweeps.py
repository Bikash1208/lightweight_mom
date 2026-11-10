"""Generate reproducible, fabrication-preflighted OAT screening CSVs.

Usage: python generate_sweeps.py [bar90|bar50]   (default bar90; bar50 writes *_bar50 files)
"""
import csv
import sys
from pathlib import Path
import cell_param_pyaedt as cp

HERE=Path(__file__).resolve().parent
PROFILE_NAME=sys.argv[1] if len(sys.argv)>1 else "bar90"
cp.set_profile(PROFILE_NAME)
# BAR90: island keeps its width, is shifted to the 0.30 mm diode gap, and its via sits at the island
# CENTROID in x and y for EVERY row (island_w, island_l, slot_w, island_offset all move it).
SCREEN_STATES=["OFF","ON"]
PHASES=[90,80,100]   # 80/100 added 2026-11-10 (criterion fix): slope across 80/90/100 separates resonance from surface wave
LATE=None   # OFF_C018 rows DROPPED 2026-11-10 (user): C sensitivity is re-run for finalists only
GROUPS_BAR90=[
    {"island_l":[4.16,4.348,4.7,5.2], "island_w":[1.0,1.15,1.3,1.45,1.6],
     "slot_l":[7.8,8.2,8.696,9.2], "slot_w":[2.6,2.75,2.9,2.94],
     "diode_pos":[-1.2,-.6,0,.6,1.2], "island_offset":[-.08,0,.08]},
    # stub 0 mm is a separate project (nostub); 0.3/1/2 mm on the full project. Anti-pad stays 0.6 mm (0.4 mm dropped:
    # not fabricable with the 0.3 mm drill and the 0.15 mm rule - user decision 2026-11-10).
    {"stub_l":[.3,1,2], "fr4_t":[.4,.8,1.2], "disc_r":[5.435,5.55,5.65,5.725],   # stub_l first: key bias-path spread result
     "via_drill":[.2,.25,.3], "pad_dia":[.6,.7,.8]},
    {"disc_via_x":[3.3,3.45,3.63], "disc_via_y":[-3.63,-3.45,-3.3],
     "stub_w":[.3,.4,.5]},
]
# Projects without stubs: nostub (stub 0 mm) and trace (periodic 0.25 mm bias trace on the FR4 bottom). Island size, anti-pad, FR4.
GROUP_VARIANT={"island_l":[4.16,4.348,4.7,5.2], "island_w":[1.0,1.15,1.3,1.45,1.6],
               "fr4_t":[.4,.8,1.2]}
# (file stem, groups, primary states, phases, late rows). SCREEN_4_nostub (90 only) was already running when the
# criterion changed, so its 80/100 partners are a separate supplement file.
SPECS_BAR90=[("SCREEN_4_nostub",GROUP_VARIANT,["OFF","ON"],[90],None),
             ("SCREEN_4_nostub_extra",GROUP_VARIANT,["OFF","ON"],[80,100],None),
             ("SCREEN_5_trace",GROUP_VARIANT,["OFF","ON"],PHASES,None),
             ("SCREEN_2",GROUPS_BAR90[1],SCREEN_STATES,PHASES,LATE),
             ("SCREEN_1",GROUPS_BAR90[0],SCREEN_STATES,PHASES,LATE),
             ("SCREEN_3",GROUPS_BAR90[2],SCREEN_STATES,PHASES,LATE)]
GROUPS_BAR50=[
    {"island_l":[4.16,4.348,4.7,5.2], "island_w":[1.28,1.3,1.4,1.45],
     "slot_l":[7.8,8.2,8.696,9.2], "slot_w":[2.6,2.75,2.9,2.94],
     "diode_pos":[-1.2,-.6,0,.6,1.2], "island_offset":[-.08,0,.08]},
    {"disc_r":[5.435,5.55,5.65,5.725], "stub_l":[.8,1.4,2],
     "fr4_t":[.4,.8,1.2], "antipad_dia":[.6,.8,1.0],
     "via_drill":[.2,.25,.3], "pad_dia":[.6,.7,.8]},
    {"island_via_x":[.52,.53,.55], "island_via_y":[-1.85,-1.78,-1.5], **GROUPS_BAR90[2]},
]
GROUPS=GROUPS_BAR90 if PROFILE_NAME=="bar90" else GROUPS_BAR50
SKIPPED=[]


def centroid_via(v):
    """Island via at the island centroid (x and y); BAR50 keeps its original via location."""
    if PROFILE_NAME=="bar90":
        v["island_via_x"]=round(-v["slot_w"]/2+v["diode_gap"]+v["island_w"]/2,6)
        v["island_via_y"]=round(v["island_offset"],6)
    return v


def write_cases(path, geometries, states, phases, late=None):
    """Rows geometry by geometry with all phases adjacent (classification completes per geometry, nominal first);
    `late` = (states, phases) rows are appended after every primary row of the file."""
    base=cp.nominal()
    columns=["*","$ph_x","$ph_y","Rs","Xs"]+[k for k in base if k not in {"$ph_x","$ph_y","Rs","Xs"}]
    out=[]
    plans=[(states,phases)]+([late] if late else [])
    for st_list,ph_list in plans:
        for name,v in geometries:
            gate=cp.fab_gate(v,name)
            if gate["status"]!="PASS":
                SKIPPED.append({"file":path.name,"geometry":name,"failures":str(gate["failures"])})
                continue
            for phi in ph_list:
                for state in st_list:
                    rs,xs=cp.state_values(state,v)
                    for direction,(px,py) in [("across",(phi,0)),("along",(0,phi))]:
                        out.append({"*":f"{path.stem}_{name}_{phi}_{state}_{direction}",**v,"$ph_x":px,"$ph_y":py,"Rs":rs,"Xs":xs})
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows(out)
    return len(list(cp.read_rows(path)))


def main():
    base=centroid_via(cp.nominal())
    tag="" if PROFILE_NAME=="bar90" else "_bar50"
    screen_states=SCREEN_STATES if PROFILE_NAME=="bar90" else list(cp.states())
    specs = SPECS_BAR90 if PROFILE_NAME=="bar90" else [(f"SCREEN_{i}",g,screen_states,[90],None) for i,g in enumerate(GROUPS,1)]
    for stem,group,states,phases,late in specs:
        geometries=[("NOMINAL",base)]
        for name,values in group.items():
            for value in values:
                if value==base[name]:
                    continue
                v=base.copy();v[name]=value
                label=f"{name}_{value:g}"
                geometries.append((label,centroid_via(v)))
        n=write_cases(HERE/f"{stem}{tag}.csv",geometries,states,phases,late)
        assert n<=400
        print(f"{stem}{tag}: {n} rows")
    if PROFILE_NAME=="bar90":
        # User-required short-island rows (2026-11-10): island_l 3.5/3.0/2.5/2.0 mm x island_w 1.3/1.0, via at centroid,
        # nostub AND trace, phases 90/80/100, OFF and ON. Shorter island = higher island resonance.
        short=[]
        for W in (1.3,1.0):
            for L in (3.5,3.0,2.5,2.0):
                v=base.copy(); v["island_l"]=L; v["island_w"]=W
                short.append((f"island_l_{L:g}_island_w_{W:g}",centroid_via(v)))
        for stem in ("SCREEN_6_short_nostub","SCREEN_6_short_trace"):
            n=write_cases(HERE/f"{stem}.csv",short,["OFF","ON"],PHASES,None)
            print(f"{stem}: {n} rows")
    # Executable nominal full-phase reference template, not solved finalists.
    write_cases(HERE/f"FULL_nominal{tag}.csv",[("NOMINAL",base)],screen_states,[90,85,95])
    if PROFILE_NAME=="bar90":
        # Diode-model sensitivity at nominal geometry: OFF C 0.13/0.18 pF and Rp 2/5 kohm.
        n=write_cases(HERE/"SENS_nominal.csv",[("NOMINAL",base)],list(cp.SENS_STATES),PHASES)
        print(f"SENS_nominal: {n} rows")
    if SKIPPED:
        with (HERE/f"SKIPPED_FAB{tag}.csv").open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=["file","geometry","failures"]);w.writeheader();w.writerows(SKIPPED)
        print(f"{len(SKIPPED)} geometries failed fabrication preflight and were skipped (SKIPPED_FAB{tag}.csv)")


if __name__=="__main__":
    main()
