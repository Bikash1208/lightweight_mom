"""Read-only comparison of B-ALL-SHORT source and corrected DIODE V2 settings."""
import json
import os
from pathlib import Path
import traceback

import psutil

os.environ["ANSYSEM_ROOT261"]=r"C:\Program Files\ANSYS Inc\v261\AnsysEM"
from ansys.aedt.core import Hfss

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/"hfss_work/full_board_gate_v31_short/solve_B_ALL_SHORT_001/BOARD_R150_GATE_V31_B_ALL_SHORT.aedt"
TARGET=ROOT/"hfss_work/full_board_gate_v31_diode_v2/BOARD_R150_V31_B_ALL_DIODE_V2.aedt"
OUT=ROOT/"results/sub6_meander_model/full_board_gate/v31_diode_v2"


def own_close(h):
    pid=int(h.desktop_class.aedt_process_id)
    try:p=psutil.Process(pid);created=p.create_time()
    except psutil.NoSuchProcess:p=None
    h.release_desktop(close_projects=True,close_desktop=True)
    if p is not None:
        try:p.wait(timeout=10)
        except psutil.TimeoutExpired:
            if p.create_time()==created:p.terminate()
        except psutil.NoSuchProcess:pass


def read(path,design,label):
    temp=OUT/f"readback_temp_{label}";temp.mkdir(exist_ok=True)
    os.environ.update(TEMP=str(temp),TMP=str(temp))
    h=Hfss(project=str(path),design=design,solution_type="Modal",
           version="2026.1",non_graphical=True,new_desktop=True,close_on_exit=False)
    try:
        h.desktop_class.temp_directory=temp
        assert str(h.desktop_class.aedt_version_id)=="2026.1"
        assert Path(h.desktop_class.temp_directory).resolve()==temp.resolve()
        m=h.modeler
        objects={}
        for name in ("PORT","ABSORBER","SUBSTRATE","BIAS_RT5880","AIR"):
            q=m[name]
            objects[name]={"bbox_mm":[float(v) for v in q.bounding_box],
                           "volume_mm3":float(q.volume),"material":str(q.material_name),
                           "area_mm2":sum(float(f.area) for f in q.faces) if name=="PORT" else None}
        setup=h.get_setup("Setup1")
        setup_props={k:str(setup.props.get(k)) for k in
                     ("Frequency","MaximumPasses","MinimumPasses","MaxDeltaS","BasisOrder","DrivenSolverType")}
        mesh={"operation_names":list(h.mesh.meshoperation_names),
              "initial_settings":str(h.mesh.initial_mesh_settings.props)}
        return {"objects":objects,"boundaries":{b.name:dict(b.props) for b in h.boundaries},
                "setup":setup_props,"mesh":mesh,
                "sphere":{v.name:dict(v.props) for v in h.field_setups}}
    finally:own_close(h)


def main():
    src=read(SOURCE,"GATE_B_ALL_SHORT","source")
    dst=read(TARGET,"GATE_B_ALL_DIODE_V2","v2")
    checks={"port_absorber_substrate_bias_air":src["objects"]==dst["objects"],
            "old_boundaries":src["boundaries"]=={k:v for k,v in dst["boundaries"].items() if k!="DIODE"},
            "setup":src["setup"]==dst["setup"],
            "mesh":src["mesh"]==dst["mesh"],"sphere":src["sphere"]==dst["sphere"]}
    record={"checks":checks,"source":src,"v2":dst,"solve_launched":False}
    (OUT/"unchanged_readback.json").write_text(json.dumps(record,indent=2,default=str)+"\n",encoding="utf-8")
    assert all(checks.values()),checks
    print("PASS",checks,flush=True)


if __name__=="__main__":
    try:main()
    except BaseException:
        (OUT/"unchanged_readback_error.log").write_text(traceback.format_exc(),encoding="utf-8")
        raise
