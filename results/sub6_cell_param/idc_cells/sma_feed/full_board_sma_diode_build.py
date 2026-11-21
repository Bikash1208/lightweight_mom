"""Build corrected-frame DIODE V2 from converged B-ALL-SHORT; no solve."""

import json
import math
import os
from pathlib import Path
import shutil
import threading
import traceback
from datetime import datetime, timezone
import zipfile

import psutil

import argparse
VERSION = "2025.1"   # set from --aedt-version in main()
from ansys.aedt.core import Hfss

SOURCE = GEOMETRY = WORK = PROJECT = OUT = ARCHIVE = STATUS = None
ZR = ZX = "0ohm"
SRC_DESIGN = "SMA_B_ALL_SHORT"
NEW_DESIGN = "SMA_B_ALL_DIODE"
EXPECTED_METAL = 517   # 516 discs + ground; the SMA exporter drops the probe pin (the old probe board had 518)


def utc(): return datetime.now(timezone.utc).isoformat()


def save(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+".tmp")
    tmp.write_text(json.dumps(data,indent=2,default=str)+"\n",encoding="utf-8")
    os.replace(tmp,path)


def need(ok,label):
    if not ok: raise RuntimeError(label)
    return ok


def row_xy(row):
    return tuple(round(float(v),3) for v in row["xy_mm"])


def main():
    global VERSION, SOURCE, GEOMETRY, WORK, PROJECT, OUT, ARCHIVE, STATUS, ZR, ZX, SRC_DESIGN, NEW_DESIGN
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", type=Path, required=True, help="base SMA B_ALL_SHORT project (.aedt) built by full_board_sma_hfss.py --build-only")
    ap.add_argument("--manifest", type=Path, required=True, help="manifest_B_ALL_SHORT_SMA.json from full_board_sma_export.py")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--tag", required=True, help="state tag used in the project name, e.g. ONC")
    ap.add_argument("--zr", required=True, help="sheet resistance ohm/sq, e.g. 24.0")
    ap.add_argument("--zx", required=True, help="sheet reactance ohm/sq, e.g. -234.0")
    ap.add_argument("--aedt-version", default="2025.1")
    ap.add_argument("--src-design", default="SMA_B_ALL_SHORT")
    ap.add_argument("--new-design", default="SMA_B_ALL_DIODE")
    a = ap.parse_args()
    VERSION = a.aedt_version
    root_env = "ANSYSEM_ROOT" + VERSION.replace(".", "")
    if root_env not in os.environ:
        os.environ[root_env] = str(Path("C:/Program Files/ANSYS Inc") / ("v" + VERSION.replace(".", "")) / "AnsysEM")
    SOURCE = a.source.resolve(); GEOMETRY = a.manifest.resolve(); OUT = a.out.resolve()
    WORK = OUT / "work"; PROJECT = WORK / f"BOARD_R150_V31_B_ALL_DIODE_SMA_{a.tag}.aedt"
    ARCHIVE = OUT / f"BOARD_R150_V31_B_ALL_DIODE_SMA_{a.tag}.aedtz"; STATUS = OUT / "build_status.json"
    ZR = f"{float(a.zr)}ohm"; ZX = f"{float(a.zx)}ohm"
    SRC_DESIGN = a.src_design; NEW_DESIGN = a.new_design
    OUT.mkdir(parents=True,exist_ok=True)
    WORK.mkdir(parents=True,exist_ok=True)
    need(SOURCE.is_file() and GEOMETRY.is_file(),"SOURCE_MISSING")
    need(not PROJECT.exists() and not ARCHIVE.exists(),"REFUSE_EXISTING_OUTPUT")
    temp=WORK/"aedt_temp"
    temp.mkdir(exist_ok=True)
    os.environ.update(TEMP=str(temp),TMP=str(temp))
    status={"phase":"START","controller_pid":os.getpid(),"start_utc":utc(),
            "source_project":str(SOURCE),"output_project":str(PROJECT)}
    save(STATUS,status)
    stop=threading.Event()
    def beat():
        while not stop.wait(60):
            status["heartbeat_utc"]=utc();save(STATUS,status)
    threading.Thread(target=beat,daemon=True).start()
    h=None
    try:
        # Copy only the project file. Its CAD solids are embedded in the AEDT file.
        shutil.copy2(SOURCE,PROJECT)
        h=Hfss(project=str(PROJECT),design=SRC_DESIGN,solution_type="Modal",
               version=VERSION,non_graphical=True,new_desktop=True,close_on_exit=False)
        status["aedt_pid"]=h.desktop_class.aedt_process_id
        status["aedt_version"]=str(h.desktop_class.aedt_version_id)
        need(status["aedt_version"]==VERSION,"AEDT_VERSION")
        h.desktop_class.temp_directory=temp
        need(Path(h.desktop_class.temp_directory).resolve()==temp.resolve(),"TEMP_READBACK")
        m=h.modeler
        m.model_units="mm"
        imported=[n for n in m.object_names if n.startswith("METAL_PEC")]+[n for n in ("SUBSTRATE","BIAS_RT5880","AIR","Pin","Dielectric","Ground_Plane") if n in m.object_names]
        need(m.purge_history(imported),"PURGE_HISTORY")   # self-contained geometry before the booleans, as in the probe-feed V2 pipeline
        need(h.design_name==SRC_DESIGN,"SOURCE_DESIGN_READBACK")
        setup=h.get_setup("Setup1")
        setup_before={k:str(setup.props.get(k)) for k in ("Frequency","MaximumPasses","MinimumPasses",
                     "MaxDeltaS","BasisOrder","DrivenSolverType")}
        boundary_before={b.name:dict(b.props) for b in h.boundaries}
        sphere_before=[v.name for v in h.field_setups]
        geometry=json.loads(GEOMETRY.read_text(encoding="utf-8"))
        targets=[row for row in geometry["rows"] if row["beam"]=="B" and not row["gate"]]
        need(len(targets)==250,"TARGET_COUNT")
        need(sum(bool(row["gate"]) for row in geometry["rows"])==16,"GATE_COUNT")
        metal_names=[n for n in m.object_names if n.startswith("METAL_PEC")]
        need(len(metal_names)==EXPECTED_METAL,"SOURCE_METAL_COUNT")
        by_xy={}
        for name in metal_names:
            bb=[float(x) for x in m[name].bounding_box]
            if abs(bb[5]-bb[2]-0.035)>1e-4 or abs(bb[3]-bb[0]-10.9)>.02 or abs(bb[4]-bb[1]-10.9)>.02:
                continue
            xy=(round((bb[0]+bb[3])/2,3),round((bb[1]+bb[4])/2,3))
            need(xy not in by_xy,f"DUPLICATE_DISC_{xy}")
            by_xy[xy]=name
        need(len(by_xy)==516,f"DISC_BODY_COUNT_{len(by_xy)}")
        need(all(row_xy(row) in by_xy for row in geometry["rows"]),"SITE_BODY_MAP")
        non_targets={}
        for row in geometry["rows"]:
            if row["beam"]=="B" and not row["gate"]:
                continue
            obj=m[by_xy[row_xy(row)]]
            non_targets[obj.name]={"volume_mm3":float(obj.volume),
                                   "bounding_box_mm":[float(v) for v in obj.bounding_box],
                                   "beam":row["beam"],"gate":row["gate"]}
        need(len(non_targets)==266,"NON_TARGET_COUNT")
        status.update(phase="CUTTING",mapped_discs=len(by_xy),target_sheets=len(targets))
        save(STATUS,status)
        cut_audit=[]
        sheets=[]
        for idx,row in enumerate(targets):
            x,y=row_xy(row)
            angle=float(row["rotation_deg"])
            name=by_xy[(x,y)]
            body=m[name]
            volume_before=float(body.volume)
            cutter=need(m.create_box([-.30,-.15,-.035],[.6,.3,.105],
                                     name=f"DIODE_GAP_CUT_{idx:03d}",material="vacuum"),
                        f"CREATE_GAP_{idx}")
            need(cutter.rotate("Z",angle),f"ROTATE_CUT_{idx}")
            need(cutter.move([x,y,0]),f"MOVE_CUT_{idx}")
            need(m.subtract(name,cutter.name,keep_originals=False),f"SUBTRACT_{idx}")
            if name not in m.object_names:     # pyaedt's object cache can lag behind AEDT: refresh once, then re-check
                try: m.refresh_all_ids()
                except Exception: pass
            need(name in m.object_names,f"DISC_NAME_AFTER_CUT_{idx}: {sorted(set(m.object_names)-set(by_xy.values()))[:8]}")
            volume_after=float(m[name].volume)
            delta=volume_before-volume_after
            # Imported STEP facets perturb the analytic 0.00525 mm3 by up to
            # ~1e-5 mm3 across rotations. The probes below check the actual gap.
            need(abs(delta-.5*.3*.035)<5e-5,f"GAP_VOLUME_{idx}: {delta}")
            probes=[]
            ca,sa=math.cos(math.radians(angle)),math.sin(math.radians(angle))
            for u in (-.2,0,.2):
                px,py=x+ca*u,y+sa*u
                bodies=[n for n in m.get_bodynames_from_position([px,py,.0175])
                        if n.startswith("METAL_PEC")]
                probes.append({"local_uv_mm":[u,0],"metal_bodies":bodies})
            need(all(not p["metal_bodies"] for p in probes),f"GAP_PROBE_{idx}: {probes}")
            for v in (-.3,.3):
                px,py=x-sa*v,y+ca*v
                bodies=[n for n in m.get_bodynames_from_position([px,py,.0175])
                        if n.startswith("METAL_PEC")]
                need(name in bodies,f"STRIP_HALF_{idx}_{v}")
            sheet=need(m.create_rectangle("XY",[-.25,-.15,.035],[.5,.3],
                                          name=f"DIODE_{idx:03d}"),f"CREATE_SHEET_{idx}")
            need(sheet.rotate("Z",angle),f"ROTATE_SHEET_{idx}")
            need(sheet.move([x,y,0]),f"MOVE_SHEET_{idx}")
            area=float(sheet.faces[0].area)
            need(abs(area-.15)<1e-5,f"SHEET_AREA_{idx}: {area}")
            sheets.append(sheet.name)
            cut_audit.append({"i":row["i"],"j":row["j"],"beam":row["beam"],
                              "copper_object":name,"sheet":sheet.name,"xy_mm":[x,y],
                              "rotation_deg":angle,"removed_copper_mm3":delta,
                              "sheet_area_mm2":area,"sheet_z_mm":.035,
                              "gap_probes":probes})
            if (idx+1)%25==0:
                status.update(cut_count=idx+1,heartbeat_utc=utc());save(STATUS,status)
        need(len([n for n in m.object_names if n.startswith("METAL_PEC")])==EXPECTED_METAL,
             "COPPER_OBJECT_COUNT_CHANGED")
        need(len(set(sheets))==250,"SHEET_COUNT")
        for name,old in non_targets.items():
            need(name in m.object_names,f"NON_TARGET_LOST_{name}")
            need(abs(float(m[name].volume)-old["volume_mm3"])<1e-6,
                 f"NON_TARGET_VOLUME_{name}")
            need(all(abs(float(a)-b)<1e-4 for a,b in
                     zip(m[name].bounding_box,old["bounding_box_mm"])),
                 f"NON_TARGET_POSITION_{name}")
        h["Zr"]=ZR
        h["Zx"]=ZX
        diode=need(h.assign_impedance_to_sheet(sheets,name="DIODE",resistance="Zr",
                                               reactance="Zx"),"DIODE_BOUNDARY_ASSIGN")
        props=dict(diode.props)
        need(props.get("Resistance")=="Zr" and props.get("Reactance")=="Zx",
             "IMPEDANCE_READBACK")
        need(set(props.get("Objects",[]))==set(sheets),"BOUNDARY_SHEET_MEMBERSHIP")
        need(len(props.get("Objects",[]))==250,"BOUNDARY_SHEET_COUNT")
        def _ohm(v): return float(str(v).lower().replace("ohm","").replace(" ",""))
        status["zr_zx_readback"]=[str(h["Zr"]),str(h["Zx"])]
        need(abs(_ohm(h["Zr"])-_ohm(ZR))<1e-9,f"ZR_READBACK {h['Zr']}")
        need(abs(_ohm(h["Zx"])-_ohm(ZX))<1e-9,f"ZX_READBACK {h['Zx']}")
        need(h.rename_design(NEW_DESIGN,save=False),"DESIGN_RENAME")
        setup_after=h.get_setup("Setup1")
        need({k:str(setup_after.props.get(k)) for k in setup_before}==setup_before,"SETUP_CHANGED")
        boundary_after={b.name:dict(b.props) for b in h.boundaries if b.name!="DIODE"}
        need(boundary_after==boundary_before,"OTHER_BOUNDARY_CHANGED")
        need([v.name for v in h.field_setups]==sphere_before,"SPHERE_CHANGED")
        valid=h.validate_simple(str(OUT/"validation.log"))
        need(valid==1 or valid is True,"DESIGN_VALIDATION")
        need(h.save_project(),"SAVE_PROJECT")
        audit={"source_project_bytes":SOURCE.stat().st_size,"project_path":str(PROJECT),
               "project_bytes":PROJECT.stat().st_size,"archive_path":str(ARCHIVE),
               "aedt_version":status["aedt_version"],"design":h.design_name,
               "source_metal_objects":EXPECTED_METAL,"final_metal_objects":EXPECTED_METAL,
               "diode_sheet_count":len(sheets),"diode_boundary":props,
               "non_target_cells_unchanged":len(non_targets),
               "non_target_gate_discs_unchanged":16,
               "diode_removed_expected_mm3":.00525,
               "diode_removed_tolerance_mm3":5e-5,
               "all_gap_probes_pass":all(not p["metal_bodies"] for r in cut_audit for p in r["gap_probes"]),
               "solve_launched":False,"mesh_launched":False,
               "variables":{"Zr":str(h["Zr"]),"Zx":str(h["Zx"])},
               "other_states_ohm_sq":{"ONC":[24.0,-234.0],"OFF_BAR90":[2.5,12.0]},
               "setup_before_after":setup_before,"other_boundaries":list(boundary_before),
               "infinite_spheres":sphere_before,"validation":valid,"cuts":cut_audit}
        save(OUT/"audit.json",audit)
        status.update(phase="BUILT",end_utc=utc(),
                      project_bytes=audit["project_bytes"])
        save(STATUS,status)
        (OUT/"build.done").write_text("BUILT\n",encoding="utf-8")
    except BaseException:
        status.update(phase="FAILED",traceback=traceback.format_exc(),end_utc=utc())
        save(STATUS,status)
        (OUT/"build.done").write_text("FAILED\n",encoding="utf-8")
        raise
    finally:
        stop.set()
        if h is not None:
            own_pid=int(h.desktop_class.aedt_process_id)
            try:
                own=psutil.Process(own_pid);created=own.create_time()
            except psutil.NoSuchProcess: own=None
            h.release_desktop(close_projects=True,close_desktop=True)
            if own is not None:
                try: own.wait(timeout=10)
                except psutil.TimeoutExpired:
                    if own.create_time()==created: own.terminate()
                except psutil.NoSuchProcess: pass


if __name__=="__main__": main()
