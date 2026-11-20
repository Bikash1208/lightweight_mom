"""Build-only bias-via variant from a copy of the unsolved v3.1 diode board."""

import json
import math
import os
from pathlib import Path
import shutil
import threading
import traceback
import zipfile
from datetime import datetime, timezone

import psutil

os.environ["ANSYSEM_ROOT261"]=r"C:\Program Files\ANSYS Inc\v261\AnsysEM"
from ansys.aedt.core import Hfss

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/"hfss_work/full_board_gate_v31_diode/BOARD_R150_V31_B_ALL_DIODE.aedt"
SOURCE_AUDIT=ROOT/"results/sub6_meander_model/full_board_gate/v31_diode_build/audit.json"
GEOM=ROOT/"results/sub6_meander_model/full_board_gate/v31_short/geometry.json"
WORK=ROOT/"hfss_work/full_board_gate_v31_biasvia"
PROJECT=WORK/"BOARD_R150_V31_B_ALL_DIODE_BIASVIA.aedt"
SCRATCH=WORK/"BIASVIA_ISOLATION_SCRATCH.aedt"
OUT=ROOT/"results/sub6_meander_model/full_board_gate/v31_biasvia_build"
ARCHIVE=OUT/"BOARD_R150_V31_B_ALL_DIODE_BIASVIA.aedtz"
STATUS=OUT/"build_status.json"
AUDIT=OUT/"audit_biasvia.json"
R=5.45
CU=.035
P=11.6
Z_TOP=.035


def utc(): return datetime.now(timezone.utc).isoformat()


def write(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+".tmp")
    tmp.write_text(json.dumps(data,indent=2,default=str)+"\n",encoding="utf-8")
    os.replace(tmp,path)


def need(ok,label):
    if not ok: raise RuntimeError(label)
    return ok


def world(x,y,angle,u,v):
    a=math.radians(angle)
    return (x+math.cos(a)*u-math.sin(a)*v,
            y+math.sin(a)*u+math.cos(a)*v)


def rim_loss_per_cell():
    a=.15
    primitive=lambda x:.5*(x*math.sqrt(R*R-x*x)+R*R*math.asin(x/R))
    one_end_area=2*primitive(a)-.3*4.36
    return 2*one_end_area*CU


def own_session(path,design,temp):
    temp.mkdir(parents=True,exist_ok=True)
    os.environ.update(TEMP=str(temp),TMP=str(temp))
    h=Hfss(project=str(path),design=design,solution_type="Modal",
           version="2026.1",non_graphical=True,new_desktop=True,close_on_exit=False)
    need(str(h.desktop_class.aedt_version_id)=="2026.1","AEDT_VERSION")
    h.desktop_class.temp_directory=temp
    need(Path(h.desktop_class.temp_directory).resolve()==temp.resolve(),"PRIVATE_TEMP")
    return h


def close_own(h):
    if h is None:return
    pid=int(h.desktop_class.aedt_process_id)
    try:own=psutil.Process(pid);created=own.create_time()
    except psutil.NoSuchProcess:own=None
    h.release_desktop(close_projects=True,close_desktop=True)
    if own is not None:
        try:own.wait(timeout=10)
        except psutil.TimeoutExpired:
            if own.create_time()==created:own.terminate()
        except psutil.NoSuchProcess:pass


def rotated_box(m,name,origin,size,x,y,angle,material="vacuum"):
    obj=need(m.create_box(origin,size,name=name,material=material),"BOX_"+name)
    need(obj.rotate("Z",angle),"ROTATE_"+name)
    need(obj.move([x,y,0]),"MOVE_"+name)
    return obj


def rotated_sheet(m,name,plane,origin,size,x,y,angle,expected_area):
    obj=need(m.create_rectangle(plane,origin,size,name=name),"SHEET_"+name)
    need(obj.rotate("Z",angle),"ROTATE_"+name)
    need(obj.move([x,y,0]),"MOVE_"+name)
    area=float(obj.faces[0].area)
    need(abs(area-expected_area)<1e-5,f"SHEET_AREA_{name}_{area}")
    return obj,area


def build(h,rows,source_audit,state):
    m=h.modeler
    m.model_units="mm"
    metals=[n for n in m.object_names if n.startswith("METAL_PEC")]
    need(len(metals)==518,"SOURCE_METAL_COUNT")
    target_names={r["copper_object"] for r in rows}
    need(len(target_names)==250 and all(n in metals for n in target_names),"TARGET_COPPER_MAPPING")
    diode=next((b for b in h.boundaries if b.name=="DIODE"),None)
    need(diode is not None and len(diode.props.get("Objects",[]))==250,"DIODE_SOURCE")
    boundaries_before={b.name:dict(b.props) for b in h.boundaries}
    setup=h.get_setup("Setup1")
    setup_before={k:str(setup.props.get(k)) for k in
                  ("Frequency","MaximumPasses","MinimumPasses","MaxDeltaS","BasisOrder","DrivenSolverType")}
    need(setup_before["Frequency"]=="5.68GHz" and "direct" in setup_before["DrivenSolverType"].lower(),
         "SOURCE_SETUP")
    sphere_before=[v.name for v in h.field_setups]
    need(sphere_before==["Global"],"SOURCE_SPHERE")
    # Identify the original ground, not a cell or feed object, by its z and span.
    ground=[]
    by_xy={}
    for name in metals:
        bb=[float(v) for v in m[name].bounding_box]
        if bb[3]-bb[0]>300 and bb[4]-bb[1]>300 and abs(bb[2]+3.21)<.01:
            ground.append(name)
        if abs(bb[5]-bb[2]-CU)<1e-4 and abs(bb[3]-bb[0]-10.9)<.02 and abs(bb[4]-bb[1]-10.9)<.02:
            xy=(round((bb[0]+bb[3])/2,3),round((bb[1]+bb[4])/2,3))
            need(xy not in by_xy,"DUPLICATE_SITE")
            by_xy[xy]=name
    need(len(ground)==1 and len(by_xy)==516,"GROUND_OR_SITE_MAPPING")
    ground=ground[0]
    geometry=json.loads(GEOM.read_text(encoding="utf-8"))
    other_rows=[r for r in geometry["rows"] if not (r["beam"]=="B" and not r["gate"])]
    need(len(other_rows)==266 and sum(bool(r["gate"]) for r in other_rows)==16,"NON_TARGET_ROWS")
    others={}
    for row in other_rows:
        xy=tuple(round(float(v),3) for v in row["xy_mm"])
        name=by_xy[xy]
        others[name]={"xy_mm":xy,"beam":row["beam"],"gate":row["gate"],
                      "volume_mm3":float(m[name].volume),
                      "bounding_box_mm":[float(v) for v in m[name].bounding_box]}
    metal_volume_before=sum(float(m[n].volume) for n in metals)
    substrate_volume_before=float(m["SUBSTRATE"].volume)
    cap_names=[];choke_names=[];hole_names=[];lower_names=[];cells=[]
    state.update(phase="BUILDING",ground_object=ground,non_target_count=len(others))
    write(STATUS,state)
    for index,row in enumerate(rows):
        x,y=(float(v) for v in row["xy_mm"])
        angle=float(row["rotation_deg"])
        disc=row["copper_object"]
        original_volume=float(m[disc].volume)
        cuts=[]
        for side,ya,yb in (("P",4.30,5.76),("M",-5.76,-4.30)):
            cuts.append(rotated_box(m,f"RIM_CUT_{index:03d}_{side}",
                                    [-.15,ya,-CU],[.3,yb-ya,3*CU],x,y,angle))
        need(m.subtract(disc,[c.name for c in cuts],keep_originals=False),f"RIM_SUBTRACT_{index}")
        need(disc in m.object_names,f"RIM_BODY_NAME_{index}")
        rim_loss=original_volume-float(m[disc].volume)
        need(abs(rim_loss-rim_loss_per_cell())<.0001,f"RIM_VOLUME_{index}_{rim_loss}")
        caps=[]
        for side,cy in (("P",4.96),("M",-4.96)):
            sheet,area=rotated_sheet(m,f"CAP_{index:03d}_{side}","XY",
                                     [-.15,cy-.25,Z_TOP],[.3,.5],x,y,angle,.15)
            cap_names.append(sheet.name);caps.append({"name":sheet.name,"area_mm2":area})
        uppers=[];lowers=[];chokes=[];holes=[]
        for side,cx in (("P",3.3),("M",-3.3)):
            holes.append(rotated_box(m,f"SUB_HOLE_{index:03d}_{side}",
                                     [cx-.15,-.15,-3.175],[.3,.3,3.175],x,y,angle))
            upper=rotated_box(m,f"UPPER_{index:03d}_{side}",
                              [cx-.15,-.15,-1.6],[.3,.3,1.635],x,y,angle,"pec")
            need(m.unite([disc,upper.name]),f"UPPER_UNITE_{index}_{side}")
            need(disc in m.object_names,f"UPPER_BODY_NAME_{index}_{side}")
            uppers.append(upper.name)
            lower=rotated_box(m,f"LOWER_{index:03d}_{side}",
                              [cx-.15,-.15,-3.21],[.3,.3,1.31],x,y,angle,"pec")
            lower_names.append(lower.name);lowers.append(lower.name)
            sheet,area=rotated_sheet(m,f"CHOKE_{index:03d}_{side}","XZ",
                                     [cx-.15,0,-1.9],[.3,.3],x,y,angle,.09)
            choke_names.append(sheet.name);chokes.append({"name":sheet.name,"area_mm2":area})
        hole_names.extend(q.name for q in holes)
        cells.append({"i":row["i"],"j":row["j"],"site_xy_mm":[x,y],
                      "rotation_deg":angle,"copper_object":disc,"diode_sheet":row["sheet"],
                      "rim_loss_mm3":rim_loss,"cap_sheets":caps,"choke_sheets":chokes,
                      "upper_barrel_count":2,"lower_barrel_names":lowers,
                      "copper_volume_after_upper_mm3":float(m[disc].volume)})
        if (index+1)%25==0:
            state.update(cells_built=index+1,heartbeat_utc=utc());write(STATUS,state)
    need(len(cap_names)==500 and len(choke_names)==500 and len(hole_names)==500 and len(lower_names)==500,
         "PART_COUNTS")
    state.update(phase="SUBTRACTING_SUBSTRATE",cells_built=250);write(STATUS,state)
    need(m.subtract("SUBSTRATE",hole_names,keep_originals=False),"SUBSTRATE_HOLES")
    need("SUBSTRATE" in m.object_names,"SUBSTRATE_NAME")
    expected_sub_loss=500*.3*.3*3.175
    actual_sub_loss=substrate_volume_before-float(m["SUBSTRATE"].volume)
    need(abs(actual_sub_loss-expected_sub_loss)<.1,"SUBSTRATE_HOLE_VOLUME")
    state.update(phase="UNITING_GROUND");write(STATUS,state)
    for start in range(0,len(lower_names),25):
        need(m.unite([ground]+lower_names[start:start+25]),f"GROUND_UNITE_{start}")
        need(ground in m.object_names,"GROUND_NAME")
    need(len([n for n in m.object_names if n.startswith("METAL_PEC")])==518,
         "METAL_OBJECT_COUNT")
    for name,item in others.items():
        need(name in m.object_names,f"NON_TARGET_LOST_{name}")
        need(abs(float(m[name].volume)-item["volume_mm3"])<1e-6,
             f"NON_TARGET_VOLUME_{name}")
        need(all(abs(float(a)-b)<1e-4 for a,b in zip(m[name].bounding_box,item["bounding_box_mm"])),
             f"NON_TARGET_POSITION_{name}")
    h["Zrc"]="0.5ohm";h["Xrc"]="0ohm"
    h["Zvr"]="1ohm";h["Zvx"]="1000ohm"
    cap=need(h.assign_impedance_to_sheet(cap_names,name="CAP",resistance="Zrc",reactance="Xrc"),
             "CAP_BOUNDARY")
    choke=need(h.assign_impedance_to_sheet(choke_names,name="CHOKE",resistance="Zvr",reactance="Zvx"),
               "CHOKE_BOUNDARY")
    for b,names,rv,xv in ((cap,cap_names,"Zrc","Xrc"),(choke,choke_names,"Zvr","Zvx")):
        props=dict(b.props)
        need(props.get("Resistance")==rv and props.get("Reactance")==xv,"BOUNDARY_VARIABLES")
        need(len(props.get("Objects",[]))==500 and set(props["Objects"])==set(names),
             "BOUNDARY_MEMBERSHIP")
    need(len(next(b for b in h.boundaries if b.name=="DIODE").props.get("Objects",[]))==250,
         "DIODE_COUNT_AFTER")
    need({b.name:dict(b.props) for b in h.boundaries if b.name not in ("CAP","CHOKE")}==boundaries_before,
         "OLD_BOUNDARIES_CHANGED")
    need({k:str(h.get_setup("Setup1").props.get(k)) for k in setup_before}==setup_before,
         "SETUP_CHANGED")
    need([v.name for v in h.field_setups]==sphere_before,"SPHERE_CHANGED")
    metal_volume_after=sum(float(m[n].volume) for n in metals)
    expected_per_cell=2*(.3*.3*1.6+.3*.3*1.275)-rim_loss_per_cell()
    expected_delta=250*expected_per_cell
    actual_delta=metal_volume_after-metal_volume_before
    need(abs(actual_delta-expected_delta)<.1,f"METAL_VOLUME_{actual_delta}_EXPECTED_{expected_delta}")
    validation=h.validate_simple(str(OUT/"validation.log"))
    need(validation==1 or validation is True,"VALIDATION")
    need(h.save_project(),"SAVE_PROJECT")
    state.update(phase="SAVED_BEFORE_ISOLATION",validation=validation)
    write(STATUS,state)
    return {"cells":cells,"source_metal_count":518,"final_metal_count":518,
            "source_metal_volume_mm3":metal_volume_before,"final_metal_volume_mm3":metal_volume_after,
            "expected_rim_loss_per_cell_mm3":rim_loss_per_cell(),
            "expected_metal_delta_per_cell_mm3":expected_per_cell,
            "expected_metal_delta_total_mm3":expected_delta,
            "actual_metal_delta_total_mm3":actual_delta,
            "volume_tolerance_mm3":.1,"substrate_expected_hole_loss_mm3":expected_sub_loss,
            "substrate_actual_hole_loss_mm3":actual_sub_loss,
            "ground_object":ground,"non_target_cells_unchanged":len(others),
            "non_target_gate_discs_unchanged":16,"diode_sheets":250,
            "cap_sheets":len(cap_names),"choke_sheets":len(choke_names),
            "cap_area_mm2":.15,"choke_area_mm2":.09,
            "boundaries":{"DIODE":dict(diode.props),"CAP":dict(cap.props),"CHOKE":dict(choke.props)},
            "variables":{"Zrc":str(h["Zrc"]),"Xrc":str(h["Xrc"]),
                         "Zvr":str(h["Zvr"]),"Zvx":str(h["Zvx"])},
            "alternates_ohm_sq":{"CAP":[.5,18.333],"CHOKE":[1000,0]},
            "setup":setup_before,"sphere":sphere_before,
            "source_boundary_names":list(boundaries_before),"validation":validation,
            "bias_layer_traces_modelled":False,"ground_clearance_holes_modelled":False,
            "solve_launched":False,"mesh_launched":False}


def isolation(rows,audit,state):
    need(not SCRATCH.exists(),"REFUSE_EXISTING_SCRATCH")
    shutil.copy2(PROJECT,SCRATCH)
    h=None
    failures=[]
    try:
        h=own_session(SCRATCH,"GATE_B_ALL_DIODE",WORK/"scratch_temp")
        state.update(phase="ISOLATION_SCRATCH",scratch_aedt_pid=h.desktop_class.aedt_process_id)
        write(STATUS,state)
        m=h.modeler
        for index,(row,record) in enumerate(zip(rows,audit["cells"])):
            name=row["copper_object"]
            x,y=(float(v) for v in row["xy_mm"])
            angle=float(row["rotation_deg"])
            reason=[]
            try:
                parts=need(m.separate_bodies(name),f"SEPARATE_{index}")
                names=[p.name for p in parts]
                if len(parts)!=2:reason.append(f"pieces={len(parts)}")
                for sign in (-1,1):
                    sx,sy=world(x,y,angle,sign*.4,0)
                    vx,vy=world(x,y,angle,sign*3.3,0)
                    strip=[n for n in m.get_bodynames_from_position([sx,sy,.0175]) if n in names]
                    barrel=[n for n in m.get_bodynames_from_position([vx,vy,-1.0]) if n in names]
                    if len(strip)!=1 or len(barrel)!=1 or strip[0]!=barrel[0]:
                        reason.append(f"side{sign}:strip={strip},barrel={barrel}")
                    record[f"side_{'minus' if sign<0 else 'plus'}_body"] = strip[0] if len(strip)==1 else None
                if record.get("side_minus_body")==record.get("side_plus_body"):
                    reason.append("halves_connected")
                if any(abs(float(p.bounding_box[2])+1.6)>.002 for p in parts):
                    reason.append("upper_barrel_z_missing")
                record["scratch_piece_names"]=names
                record["scratch_piece_count"]=len(parts)
            except Exception as exc:
                reason.append(f"{type(exc).__name__}: {exc}")
            record["isolation_pass"]=not reason
            record["isolation_failure_reason"]="; ".join(reason) if reason else None
            if reason:failures.append({"i":row["i"],"j":row["j"],"reason":record["isolation_failure_reason"]})
            if (index+1)%25==0:
                state.update(isolation_checked=index+1,heartbeat_utc=utc());write(STATUS,state)
        audit["isolation_passed"]=250-len(failures)
        audit["isolation_failures"]=failures
        write(AUDIT,audit)
        need(not failures,f"ISOLATION_FAILED_{len(failures)}")
    finally:
        close_own(h)


def embed_archive(h,audit,state):
    m=h.modeler
    imported=[n for n in m.object_names if n.startswith("METAL_PEC")]+[
        "SUBSTRATE","BIAS_RT5880","AIR","PORT"]
    need(len(imported)==522,"IMPORTED_COUNT")
    state["phase"]="PURGING_CAD_HISTORY";write(STATUS,state)
    need(m.purge_history(imported),"PURGE_HISTORY")
    need(len([n for n in m.object_names if n.startswith("METAL_PEC")])==518,"METAL_AFTER_PURGE")
    need(h.validate_simple(str(OUT/"validation_embedded.log")),"VALIDATION_EMBEDDED")
    need(h.save_project(),"SAVE_EMBEDDED")
    steprefs=PROJECT.read_bytes().lower().count(b".step")
    need(steprefs==0,f"STEP_REFS_{steprefs}")
    need(not ARCHIVE.exists(),"REFUSE_EXISTING_ARCHIVE")
    need(h.archive_project(ARCHIVE,include_external_files=False,include_results_file=False),
         "ARCHIVE_FAILED")
    with zipfile.ZipFile(ARCHIVE) as z:
        need(z.testzip() is None,"ARCHIVE_CRC")
        members=z.namelist()
        model=next(n for n in members if n.lower().endswith(".aedt"))
        archive_refs=z.read(model).lower().count(b".step")
    need(archive_refs==0 and not any(n.lower().endswith(".step") or "aedtresults" in n.lower()
                                     for n in members),"ARCHIVE_EXTERNAL_OR_RESULTS")
    audit.update(imported_history_purged=522,project_step_references=steprefs,
                 archive_step_references=archive_refs,archive_members=members,
                 project_path=str(PROJECT),project_bytes=PROJECT.stat().st_size,
                 archive_path=str(ARCHIVE),archive_bytes=ARCHIVE.stat().st_size,
                 validation_embedded=True)
    write(AUDIT,audit)


def main():
    need(SOURCE.is_file() and SOURCE_AUDIT.is_file() and GEOM.is_file(),"INPUT_MISSING")
    need(not PROJECT.exists() and not SCRATCH.exists() and not ARCHIVE.exists(),"REFUSE_EXISTING_OUTPUT")
    WORK.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
    source_audit=json.loads(SOURCE_AUDIT.read_text(encoding="utf-8"))
    rows=source_audit["cuts"]
    need(len(rows)==250,"SOURCE_TARGET_COUNT")
    state={"phase":"START","controller_pid":os.getpid(),"started_utc":utc(),
           "source_project":str(SOURCE),"output_project":str(PROJECT)}
    write(STATUS,state)
    stop=threading.Event()
    def beat():
        while not stop.wait(60):
            state["heartbeat_utc"]=utc();write(STATUS,state)
    threading.Thread(target=beat,daemon=True).start()
    h=None
    try:
        shutil.copy2(SOURCE,PROJECT)
        h=own_session(PROJECT,"GATE_B_ALL_DIODE",WORK/"build_temp")
        state.update(aedt_pid=h.desktop_class.aedt_process_id,aedt_version="2026.1")
        write(STATUS,state)
        audit=build(h,rows,source_audit,state)
        write(AUDIT,audit)
        close_own(h);h=None
        isolation(rows,audit,state)
        h=own_session(PROJECT,"GATE_B_ALL_DIODE",WORK/"embed_temp")
        state.update(phase="EMBEDDING",embed_aedt_pid=h.desktop_class.aedt_process_id)
        write(STATUS,state)
        embed_archive(h,audit,state)
        state.update(phase="COMPLETE",finished_utc=utc(),archive_bytes=ARCHIVE.stat().st_size)
        write(STATUS,state)
        (OUT/"build.done").write_text("COMPLETE\n",encoding="utf-8")
    except BaseException:
        state.update(phase="FAILED",traceback=traceback.format_exc(),finished_utc=utc())
        write(STATUS,state)
        (OUT/"build.done").write_text("FAILED\n",encoding="utf-8")
        raise
    finally:
        stop.set();close_own(h)


if __name__=="__main__":main()
