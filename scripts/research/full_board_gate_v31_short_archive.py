"""Archive converged v3.1 B-ALL-SHORT without results or external CAD links; no solve."""

import json
import math
import os
from pathlib import Path
import shutil
import traceback
import zipfile

import psutil

os.environ["ANSYSEM_ROOT261"]=r"C:\Program Files\ANSYS Inc\v261\AnsysEM"
from ansys.aedt.core import Hfss

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/"hfss_work/full_board_gate_v31_short/solve_B_ALL_SHORT_001/BOARD_R150_GATE_V31_B_ALL_SHORT.aedt"
PROFILE=ROOT/"hfss_work/full_board_gate_v31_short/solve_B_ALL_SHORT_001/status.json"
WORK=ROOT/"hfss_work/full_board_gate_v31_short_archive"
PROJECT=WORK/"BOARD_R150_V31_B_ALL_SHORT.aedt"
OUT=ROOT/"results/sub6_meander_model/full_board_gate/v31_short_archive"
ARCHIVE=OUT/"BOARD_R150_V31_B_ALL_SHORT.aedtz"
STATUS=OUT/"archive_status.json"


def save(data):
    OUT.mkdir(parents=True,exist_ok=True)
    STATUS.write_text(json.dumps(data,indent=2,default=str)+"\n",encoding="utf-8")


def need(value,label):
    if not value: raise RuntimeError(label)
    return value


def main():
    need(SOURCE.is_file() and PROFILE.is_file(),"CONVERGED_SOURCE_MISSING")
    solved=json.loads(PROFILE.read_text(encoding="utf-8"))
    profile=solved.get("profile",{})
    need(solved.get("phase")=="COMPLETE" and profile.get("converged") and profile.get("normal_completion"),
         "SOURCE_NOT_CONVERGED")
    need(not PROJECT.exists() and not ARCHIVE.exists(),"REFUSE_EXISTING_OUTPUT")
    WORK.mkdir(parents=True,exist_ok=True)
    OUT.mkdir(parents=True,exist_ok=True)
    temp=WORK/"aedt_temp";temp.mkdir(exist_ok=True)
    os.environ.update(TEMP=str(temp),TMP=str(temp))
    state={"phase":"START","controller_pid":os.getpid(),"source_project":str(SOURCE)}
    save(state)
    h=None
    try:
        shutil.copy2(SOURCE,PROJECT)
        h=Hfss(project=str(PROJECT),design="GATE_B_ALL_SHORT",solution_type="Modal",
               version="2026.1",non_graphical=True,new_desktop=True,close_on_exit=False)
        state["aedt_pid"]=h.desktop_class.aedt_process_id
        need(str(h.desktop_class.aedt_version_id)=="2026.1","VERSION")
        h.desktop_class.temp_directory=temp
        need(Path(h.desktop_class.temp_directory).resolve()==temp.resolve(),"TEMP")
        m=h.modeler
        metal=[n for n in m.object_names if n.startswith("METAL_PEC")]
        need(len(metal)==518 and len(m.object_names)==523,"SOURCE_OBJECT_COUNTS")
        imported=metal+["SUBSTRATE","BIAS_RT5880","AIR","PORT"]
        need(all(n in m.object_names for n in imported),"IMPORTED_OBJECT_MISSING")
        volume_before=sum(float(m[n].volume) for n in metal)
        port=m["PORT"]
        bb=[float(x) for x in port.bounding_box]
        area=sum(float(f.area) for f in port.faces)
        ro=max(abs(bb[0]),abs(bb[3]),abs(bb[1]),abs(bb[4]))
        ri=math.sqrt(ro**2-area/math.pi)
        need(abs(ro-2.05)<.001 and abs(ri-.635)<.002,"PORT_GEOMETRY")
        need(abs(bb[2]+3.175)<1e-4 and abs(bb[5]+3.175)<1e-4,"PORT_Z")
        boundaries={b.name:dict(b.props) for b in h.boundaries}
        need("Port1" in boundaries,"PORT_BOUNDARY_MISSING")
        port_props=boundaries["Port1"]
        need("50" in str(port_props),"PORT_IMPEDANCE")
        setup=h.get_setup("Setup1")
        setup_props={k:str(setup.props.get(k)) for k in
                     ("Frequency","MaximumPasses","MinimumPasses","MaxDeltaS","BasisOrder","DrivenSolverType")}
        need(setup_props["Frequency"]=="5.68GHz","FREQUENCY")
        need("direct" in setup_props["DrivenSolverType"].lower(),"SOLVER_TYPE")
        spheres=[v.name for v in h.field_setups]
        need(spheres==["Global"],"SPHERE")
        state.update(phase="PURGING_CAD_HISTORY",imported_objects=len(imported))
        save(state)
        need(m.purge_history(imported),"PURGE_HISTORY")
        need(len(m.object_names)==523 and len([n for n in m.object_names if n.startswith("METAL_PEC")])==518,
             "OBJECT_COUNTS_CHANGED")
        volume_after=sum(float(m[n].volume) for n in metal)
        need(abs(volume_after-volume_before)<1e-3,"METAL_VOLUME_CHANGED")
        need({b.name:dict(b.props) for b in h.boundaries}==boundaries,"BOUNDARIES_CHANGED")
        need({k:str(h.get_setup("Setup1").props.get(k)) for k in setup_props}==setup_props,
             "SETUP_CHANGED")
        need([v.name for v in h.field_setups]==spheres,"SPHERE_CHANGED")
        validation=h.validate_simple(str(OUT/"validation.log"))
        need(validation==1 or validation is True,"VALIDATION")
        need(h.save_project(),"SAVE")
        step_refs=PROJECT.read_bytes().lower().count(b".step")
        need(step_refs==0,f"EXTERNAL_STEP_REFERENCES_{step_refs}")
        need(h.archive_project(ARCHIVE,include_external_files=False,include_results_file=False),"ARCHIVE")
        with zipfile.ZipFile(ARCHIVE) as z:
            members=[(info.filename,info.file_size) for info in z.infolist()]
            need(z.testzip() is None,"ZIP_CRC")
            project_member=next(n for n,_ in members if n.lower().endswith(".aedt"))
            archive_step_refs=z.read(project_member).lower().count(b".step")
        need(archive_step_refs==0,"ARCHIVE_STEP_REFERENCES")
        need(not any(n.lower().endswith(".step") or "aedtresults" in n.lower()
                     for n,_ in members),"ARCHIVE_EXTERNAL_OR_RESULTS")
        audit={"source_project":str(SOURCE),"converged_evidence":solved,
               "saved_project":str(PROJECT),"project_bytes":PROJECT.stat().st_size,
               "archive":str(ARCHIVE),"archive_bytes":ARCHIVE.stat().st_size,
               "archive_members":members,"aedt_version":"2026.1",
               "metal_objects":len(metal),"imported_cad_history_purged":len(imported),
               "metal_volume_before_after_mm3":[volume_before,volume_after],
               "step_references_project_archive":[step_refs,archive_step_refs],
               "port":{"name":"Port1","type":"Lumped Port","impedance_ohm":50,
                       "sheet_z_mm":bb[2],"outer_radius_mm":ro,"inner_radius_mm":ri,
                       "area_mm2":area,"properties":port_props},
               "boundaries":list(boundaries),"sphere":spheres,"setup":setup_props,
               "validation":validation,"solve_launched":False}
        (OUT/"audit.json").write_text(json.dumps(audit,indent=2,default=str)+"\n",encoding="utf-8")
        state.update(phase="COMPLETE",project_bytes=PROJECT.stat().st_size,
                     archive_bytes=ARCHIVE.stat().st_size)
        save(state)
        (OUT/"archive.done").write_text("COMPLETE\n",encoding="utf-8")
    except BaseException:
        state.update(phase="FAILED",traceback=traceback.format_exc());save(state)
        (OUT/"archive.done").write_text("FAILED\n",encoding="utf-8")
        raise
    finally:
        if h is not None:
            pid=int(h.desktop_class.aedt_process_id)
            try: own=psutil.Process(pid);created=own.create_time()
            except psutil.NoSuchProcess: own=None
            h.release_desktop(close_projects=True,close_desktop=True)
            if own is not None:
                try: own.wait(timeout=10)
                except psutil.TimeoutExpired:
                    if own.create_time()==created: own.terminate()
                except psutil.NoSuchProcess: pass


if __name__=="__main__": main()
