"""One-cell corrected-frame dry run. Stops before the 250-cell builds."""

import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import traceback

import psutil

os.environ["ANSYSEM_ROOT261"]=r"C:\Program Files\ANSYS Inc\v261\AnsysEM"
from ansys.aedt.core import Hfss

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/"hfss_work/full_board_gate_v31_short/solve_B_ALL_SHORT_001/BOARD_R150_GATE_V31_B_ALL_SHORT.aedt"
SOURCE_STATUS=ROOT/"hfss_work/full_board_gate_v31_short/solve_B_ALL_SHORT_001/status.json"
ROW_AUDIT=ROOT/"results/sub6_meander_model/full_board_gate/v31_diode_build/audit.json"
WORK=ROOT/"hfss_work/full_board_gate_v31_diode_v2_dryrun"
PROJECT=WORK/"DRYRUN_B_CELL_M4_M12.aedt"
OUT=ROOT/"results/sub6_meander_model/full_board_gate/v31_diode_v2_dryrun"
STATUS=OUT/"dryrun_status.json"


def write(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,default=str)+"\n",encoding="utf-8")


def need(value,label):
    if not value:raise RuntimeError(label)
    return value


def xy(x,y,angle,u,v):
    a=math.radians(angle)
    return [x+math.cos(a)*u-math.sin(a)*v,
            y+math.sin(a)*u+math.cos(a)*v]


def cutter(m,name,origin,size,x,y,angle):
    q=need(m.create_box(origin,size,name=name,material="vacuum"),"CREATE_"+name)
    need(q.rotate("Z",angle),"ROTATE_"+name)
    need(q.move([x,y,0]),"MOVE_"+name)
    return q


def rim_loss():
    r=5.45;a=.15
    f=lambda x:.5*(x*math.sqrt(r*r-x*x)+r*r*math.asin(x/r))
    return 2*(2*f(a)-.3*4.36)*.035


def close(h):
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


def main():
    status=json.loads(SOURCE_STATUS.read_text(encoding="utf-8"))
    need(status.get("phase")=="COMPLETE","SOURCE_NOT_CONVERGED")
    row=json.loads(ROW_AUDIT.read_text(encoding="utf-8"))["cuts"][0]
    need((int(row["i"]),int(row["j"]))==(-4,-12),"WRONG_DRYRUN_SITE")
    need(SOURCE.is_file() and not PROJECT.exists(),"SOURCE_OR_OUTPUT_STATE")
    WORK.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
    original_hash=hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    shutil.copy2(SOURCE,PROJECT)
    temp=WORK/"temp";temp.mkdir(exist_ok=True)
    os.environ.update(TEMP=str(temp),TMP=str(temp))
    state={"phase":"START","controller_pid":os.getpid(),"source_hash_sha256":original_hash}
    write(STATUS,state)
    h=None
    try:
        h=Hfss(project=str(PROJECT),design="GATE_B_ALL_SHORT",solution_type="Modal",
               version="2026.1",non_graphical=True,new_desktop=True,close_on_exit=False)
        state["aedt_pid"]=h.desktop_class.aedt_process_id
        need(str(h.desktop_class.aedt_version_id)=="2026.1","VERSION")
        h.desktop_class.temp_directory=temp
        need(Path(h.desktop_class.temp_directory).resolve()==temp.resolve(),"TEMP")
        m=h.modeler;m.model_units="mm"
        name=row["copper_object"]
        x,y=(float(q) for q in row["xy_mm"])
        angle=float(row["rotation_deg"])
        need(name in m.object_names,"TARGET_BODY_MISSING")
        need(len([n for n in m.object_names if n.startswith("METAL_PEC")])==518,"SOURCE_METAL_COUNT")
        before=float(m[name].volume)
        q=cutter(m,"DRY_DIODE_GAP",[-.30,-.15,-.035],[.6,.3,.105],x,y,angle)
        need(m.subtract(name,q.name,keep_originals=False),"DIODE_SUBTRACT")
        diode_loss=before-float(m[name].volume)
        need(abs(diode_loss-.00525)<1e-5,f"DIODE_VOLUME_{diode_loss}")
        gap_probes=[]
        for u in (-.2,0,.2):
            p=xy(x,y,angle,u,0)
            found=[n for n in m.get_bodynames_from_position([*p,.0175])
                   if n.startswith("METAL_PEC")]
            gap_probes.append({"local_uv_mm":[u,0],"metal_bodies":found})
        need(all(not p["metal_bodies"] for p in gap_probes),"DIODE_GAP_PROBE")
        strip_before_rims=[]
        for v in (-.3,.3):
            p=xy(x,y,angle,0,v)
            found=[n for n in m.get_bodynames_from_position([*p,.0175])
                   if n.startswith("METAL_PEC")]
            strip_before_rims.append({"local_uv_mm":[0,v],"metal_bodies":found})
            need(name in found,"STRIP_HALF_MISSING")
        sheet=need(m.create_rectangle("XY",[-.25,-.15,.035],[.5,.3],name="DRY_DIODE"),
                   "DIODE_SHEET")
        need(sheet.rotate("Z",angle),"ROTATE_SHEET")
        need(sheet.move([x,y,0]),"MOVE_SHEET")
        area=float(sheet.faces[0].area)
        need(abs(area-.15)<1e-5,"DIODE_AREA")
        h["Zr"]="3.333ohm";h["Zx"]="-232.667ohm"
        b=need(h.assign_impedance_to_sheet(sheet.name,name="DIODE",resistance="Zr",
                                           reactance="Zx"),"BOUNDARY")
        need(dict(b.props).get("Objects")==[sheet.name],"DIODE_MEMBERSHIP")
        before_rim=float(m[name].volume)
        cuts=[cutter(m,"DRY_RIM_P",[4.30,-.15,-.035],[1.46,.3,.105],x,y,angle),
              cutter(m,"DRY_RIM_M",[-5.76,-.15,-.035],[1.46,.3,.105],x,y,angle)]
        need(m.subtract(name,[q.name for q in cuts],keep_originals=False),"RIM_SUBTRACT")
        measured_rim=before_rim-float(m[name].volume)
        expected_rim=rim_loss()
        need(abs(measured_rim-expected_rim)<1e-4,f"RIM_VOLUME_{measured_rim}")
        parts=need(m.separate_bodies(name),"SEPARATE_PARTS")
        names=[p.name for p in parts]
        need(len(names)==2,f"ISOLATION_PIECES_{len(names)}")
        strip_after=[]
        for v in (-.3,.3):
            p=xy(x,y,angle,0,v)
            found=[n for n in m.get_bodynames_from_position([*p,.0175]) if n in names]
            strip_after.append({"local_uv_mm":[0,v],"copper_piece":found})
            need(len(found)==1,"STRIP_PIECE_PROBE")
        need(strip_after[0]["copper_piece"]!=strip_after[1]["copper_piece"],
             "HALVES_STILL_CONNECTED")
        # The impedance sheet remains separate copper-wise and spans the two cut edges.
        need(abs(float(sheet.faces[0].area)-.15)<1e-5,"SHEET_AFTER_SEPARATE")
        validation=h.validate_simple(str(OUT/"validation.log"))
        need(validation==1 or validation is True,"VALIDATION")
        need(h.save_project(),"SAVE_DRYRUN")
        need(hashlib.sha256(SOURCE.read_bytes()).hexdigest()==original_hash,"SOURCE_CHANGED")
        audit={"status":"PASS_STOP_AFTER_ONE_CELL","site":{"i":row["i"],"j":row["j"],
                 "xy_mm":[x,y],"rotation_deg":angle,"original_copper_object":name},
               "frame":{"local_x":"along 8.72mm slot","local_y":"across 2.91mm slot"},
               "diode_cut_local_mm":{"x":[-.30,.30],"y":[-.15,.15]},
               "diode_removed_expected_mm3":.00525,"diode_removed_readback_mm3":diode_loss,
               "gap_probes":gap_probes,"strip_probes_before_rims":strip_before_rims,
               "diode_sheet_area_mm2":area,"diode_boundary":dict(b.props),
               "rim_cut_expected_mm3":expected_rim,"rim_cut_readback_mm3":measured_rim,
               "separate_bodies_count":len(names),"separate_bodies_names":names,
               "strip_probes_after_rims":strip_after,"validation":validation,
               "source_project_sha256":original_hash,"dryrun_project":str(PROJECT),
               "solve_launched":False,"mesh_launched":False,"output_A_built":False,
               "output_B_built":False}
        write(OUT/"dryrun_audit.json",audit)
        state.update(phase="PASS_STOP_AFTER_ONE_CELL",aedt_version="2026.1")
        write(STATUS,state)
        (OUT/"dryrun.done").write_text("PASS_STOP_AFTER_ONE_CELL\n",encoding="utf-8")
    except BaseException:
        state.update(phase="FAILED",traceback=traceback.format_exc())
        write(STATUS,state)
        (OUT/"dryrun.done").write_text("FAILED\n",encoding="utf-8")
        raise
    finally:
        close(h)


if __name__=="__main__":main()
