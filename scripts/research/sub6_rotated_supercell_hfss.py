"""Fresh local AEDT sessions for IDC supercell leakage; never attaches to an existing session.

The initial command is a single-cell open-top support probe. Projects/logs are
created in this worktree's private hfss_work directory. No process-wide search,
attachment, or wildcard deletion is performed.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import sys
import time
import traceback
import threading
from datetime import datetime, timezone

from ansys.aedt.core import Hfss


ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "hfss_work" / "rotated_supercell"
RESULTS = ROOT / "results" / "sub6_meander_model" / "supercell_hfss"
AEDT_ROOT = Path(r"C:\Program Files\ANSYS Inc\v252\AnsysEM")
VERSION = "2025.2"  # user approved after v251 executable was found absent
P = 11.6
CU = .035
RT = 3.175
BIAS = .254

CELLS = {
    "MOAT_ON": dict(kind="moat", on=True, area=None),
    "MOAT_OFF": dict(kind="moat", on=False, area=None),
    "B2p2_D3p0": dict(kind="mn", r=5.6, w=.4, depth=3.0, bridge=2.2, area=90.0),
    "B2p2_D3p6": dict(kind="mn", r=5.6, w=.4, depth=3.6, bridge=2.2, area=88.9),
    "D3p6_OPEN": dict(kind="mn", r=5.6, w=.4, depth=3.6, bridge=None, area=87.08),
    "ANCHOR": dict(kind="st", r=5.45, slot_l=8.72, slot_w=2.91, area=67.9),
}


def now():
    return datetime.now(timezone.utc).isoformat()


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, default=str, allow_nan=False) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def need(value, label):
    if not value:
        raise RuntimeError(label)
    return value


def launch(project):
    temp = WORK / "temp" / project.stem
    temp.mkdir(parents=True, exist_ok=True)
    os.environ.update(TEMP=str(temp), TMP=str(temp), ANSYSEM_ROOT252=str(AEDT_ROOT))
    h = Hfss(project=str(project), design="EigenCell", solution_type="Eigenmode",
             version=VERSION, non_graphical=True, new_desktop=True, close_on_exit=False)
    h.desktop_class.temp_directory = temp
    need(str(h.desktop_class.aedt_version_id) == VERSION, "AEDT_VERSION_READBACK")
    need(Path(h.desktop_class.temp_directory).resolve() == temp.resolve(), "PRIVATE_TEMP_READBACK")
    write(RESULTS / "active_pid.json", {"controller_pid": os.getpid(),
                                        "aedt_pid": h.desktop_class.aedt_process_id,
                                        "project": str(project), "started_utc": now()})
    return h


def build(h, cell, top, phase=85, direction="across", air_mm=20.035):
    spec = CELLS[cell]
    h.modeler.model_units = "mm"
    mat = h.materials.add_material("RT5880_SUPERCELL")
    mat.permittivity = 2.2
    mat.dielectric_loss_tangent = .0009
    h["$ph_x"] = f"{phase if direction == 'across' else 0}deg"
    h["$ph_y"] = f"{phase if direction == 'along' else 0}deg"
    for key, value in dict(period=P, cu_t=CU, rt_t=RT, bias_t=BIAS,
                           air_top=40.035 if top == "PEC" else air_mm,
                           air_bottom=-43.499).items():
        h[key] = f"{value:.12g}mm"
    m = h.modeler

    def box(name, pos, size, material="pec"):
        return need(m.create_box(pos, size, name=name, material=material), f"CREATE_{name}")

    if spec["kind"] == "moat":
        def cylinder(name, radius):
            return need(m.create_cylinder("Z", ["0mm", "0mm", "0mm"],
                                          f"{radius}mm", "cu_t", name=name,
                                          material="pec"), name)
        ring = cylinder("RING", 5.25)
        core_cut = cylinder("RING_CORE_CUT", 4.20)
        need(m.subtract(ring.name, core_cut.name, keep_originals=False), "RING_SUBTRACT")
        island = cylinder("ISLAND", 2.75)
        def sector(name, ri, ro, deg):
            angles = [deg-12+24*i/48 for i in range(49)]
            xy = [(ro*math.cos(math.radians(a)), ro*math.sin(math.radians(a))) for a in angles]
            xy += [(ri*math.cos(math.radians(a)), ri*math.sin(math.radians(a))) for a in reversed(angles)]
            points = [[f"{x:.12f}mm", f"{y:.12f}mm", "0mm"] for x,y in xy]
            poly = need(m.create_polyline(points, name=name, cover_surface=True,
                                          close_surface=True, material="pec"), name)
            need(m.thicken_sheet(poly.name, "cu_t", both_sides=False), f"THICKEN_{name}")
            return poly
        for label, angle in (("E",0),("N",90),("W",180),("S",270)):
            paddle = sector(f"PADDLE_{label}", 3.20, 3.85, angle)
            if spec["on"] and label in ("E","W"):
                bridge = sector(f"JOIN_{label}", 2.74, 3.21, angle)
                need(m.unite([island.name,paddle.name,bridge.name]), f"UNITE_{label}")
        copper_solids = [m[name] for name in m.object_names
                         if m[name].material_name.lower() == "pec"]
        expected_count = 4 if spec["on"] else 6
        need(len(copper_solids) == expected_count,
             f"MOAT_COPPER_BODY_COUNT {len(copper_solids)} expected {expected_count}")
        copper_area = sum(float(obj.volume) for obj in copper_solids) / CU
        expected_area = math.pi*(5.25**2-4.2**2+2.75**2) + 4*math.pi/15*(3.85**2-3.2**2)
        if spec["on"]:
            expected_area += 2*math.pi/15*(3.2**2-2.75**2)
        need(abs(copper_area-expected_area)/expected_area < .01,
             f"MOAT_COPPER_AREA {copper_area} expected {expected_area}")
    else:
        disc = need(m.create_cylinder("Z", ["0mm", "0mm", "0mm"],
                                      f"{spec['r']}mm", "cu_t", name="DISC", material="pec"), "DISC")
    cuts = []
    def cut_box(name, x0, y0, dx, dy):
        cuts.append(box(name, [f"{x0:.9f}mm", f"{y0:.9f}mm", "-cu_t"],
                        [f"{dx:.9f}mm", f"{dy:.9f}mm", "3*cu_t"], "vacuum"))

    if spec["kind"] == "moat":
        pass
    elif spec["kind"] == "st":
        cut_box("SLOT", -spec["slot_w"]/2, -spec["slot_l"]/2, spec["slot_w"], spec["slot_l"])
    else:
        w, depth, bridge = spec["w"], spec["depth"], spec["bridge"]
        x0 = -2.4
        for i in range(5):
            cut_box(f"LEG_{i}", x0+i*1.2-w/2, -depth/2-w/2, w, depth+w)
        for i in range(4):
            yc = depth/2 if i%2 == 0 else -depth/2
            cut_box(f"BAR_{i}", x0+i*1.2-w/2, yc-w/2, 1.2+w, w)
        rimx = math.sqrt(spec["r"]**2-(depth/2+w/2)**2)
        # For the open twin, the old rimx end stopped 0.136 mm short of the
        # rim at the inner y edge. Extend past the disc to open the full cut.
        end_abs = spec["r"] + 0.1 if bridge is None else rimx-bridge
        if end_abs <= 2.2:
            raise RuntimeError("END_RUN_LENGTH_NONPOSITIVE")
        cut_box("END_LEFT", -end_abs, -depth/2-w/2, end_abs-2.2, w)
        cut_box("END_RIGHT", 2.2, depth/2-w/2, end_abs-2.2, w)
    if spec["kind"] != "moat":
        need(m.subtract(disc.name, [c.name for c in cuts], keep_originals=False), "SLOT_SUBTRACT")
    if cell == "D3p6_OPEN":
        audit = {"end_run": {"left_x": [-end_abs, -2.2], "right_x": [2.2, end_abs],
                             "rim_x_outer_y": rimx,
                             "rim_x_inner_y": math.sqrt(spec["r"]**2-(depth/2-w/2)**2),
                             "width_mm": w},
                 "objects_after_subtract": [{"name": name,
                     "material": m[name].material_name, "volume_mm3": float(m[name].volume)}
                    for name in m.object_names]}
        audit["separate_result"] = [x.name for x in need(m.separate_bodies(disc.name),
                                                          "SEPARATE_TWIN_BODIES")]
        audit["objects_after_separation"] = [{"name": name,
                     "material": m[name].material_name, "volume_mm3": float(m[name].volume)}
                    for name in m.object_names]
        write(RESULTS / "D3p6_OPEN_copper_audit.json", audit)
    # An edge-to-edge cut splits DISC into two AEDT solids; the original handle
    # then reports zero volume. Sum the two read-back PEC solids instead.
    copper_solids = [m[name] for name in m.object_names
                     if m[name].material_name.lower() == "pec"]
    copper_area = sum(float(obj.volume) for obj in copper_solids) / CU
    if spec["kind"] != "moat" and not math.isclose(copper_area, spec["area"], rel_tol=.02):
        raise RuntimeError(f"COPPER_AREA_READBACK {cell}: {copper_area:.4f} vs expected {spec['area']}")
    ground = box("GROUND", ["-period/2", "-period/2", "-rt_t-cu_t"], ["period", "period", "cu_t"])
    substrate = box("SUBSTRATE", ["-period/2", "-period/2", "-rt_t"],
                    ["period", "period", "rt_t"], "RT5880_SUPERCELL")
    bias = box("BIAS_RT5880", ["-period/2", "-period/2", "-rt_t-cu_t-bias_t"],
               ["period", "period", "bias_t"], "RT5880_SUPERCELL")
    air = box("AIR", ["-period/2", "-period/2", "air_bottom"],
              ["period", "period", "air_top-air_bottom"], "vacuum")
    topz = 40.035 if top == "PEC" else air_mm
    zm = (topz - 43.499) / 2
    face = lambda pos: need(m.get_faceid_from_position(pos, assignment="AIR"), "AIR_FACE")
    need(h.assign_lattice_pair([face([-P/2,0,zm]),face([P/2,0,zm])],
                               phase_delay="InputPhaseDelay", phase_delay_param1="$ph_x",
                               name="Lattice_X"), "LATTICE_X")
    need(h.assign_lattice_pair([face([0,-P/2,zm]),face([0,P/2,zm])],
                               phase_delay="InputPhaseDelay", phase_delay_param1="$ph_y",
                               name="Lattice_Y"), "LATTICE_Y")
    topface = face([0,0,topz])
    if top == "PEC":
        need(h.assign_perfect_e([topface], name="PEC_LID"), "PEC_LID")
    elif top == "RADIATION":
        need(h.assign_radiation_boundary_to_faces([topface], name="RAD_TOP"), "RAD_TOP")
    elif top in ("PML", "PML2"):
        h.oboundary.CreatePML(["NAME:PMLCreationSettings", "UserDrawnGroup:=", False,
                               "PMLFaces:=", [topface], "CreateJoiningObjs:=", False,
                               "Thickness:=", "10mm", "RadDist:=", "20mm",
                               "UseFreq:=", True, "MinFreq:=", "5.8GHz"])
    else:
        raise ValueError(top)
    need(h.assign_perfect_e([face([0,0,-43.499])], name="PEC_BOTTOM"), "PEC_BOTTOM")
    setup = need(h.create_setup("EigenCell"), "SETUP")
    setup.props.update(MinimumFrequency="3.0GHz", NumModes=12, MaxDeltaFreq=.2,
                       MinimumPasses=4, MinimumConvergedPasses=2, MaximumPasses=20,
                       ConvergeOnRealFreq=True)
    need(setup.update(), "SETUP_UPDATE")
    # validate_full_design calls GetMessages, which this run explicitly forbids.
    validation_file = RESULTS / "validation" / f"{cell}_{top}.log"
    validation_file.parent.mkdir(parents=True, exist_ok=True)
    need(h.validate_simple(validation_file) == 1, "DESIGN_VALIDATION")
    validation = validation_file.read_text(errors="replace") if validation_file.exists() else ""
    boundaries = {b.name: dict(b.props) for b in h.boundaries}
    write(RESULTS / f"{cell}_{top}_boundary_probe.json",
          {"boundaries": boundaries, "objects": m.object_names,
           "boundary_names": [b.name for b in h.boundaries]})
    required = ["Lattice_X", "Lattice_Y", "PEC_BOTTOM"]
    if top == "PEC": required.append("PEC_LID")
    if top == "RADIATION": required.append("RAD_TOP")
    need(all(k in boundaries for k in required), "BOUNDARY_READBACK")
    if top in ("PML", "PML2"):
        need(any("pml" in n.lower() for n in m.object_names + list(boundaries)), "PML_READBACK")
    need(abs(float(substrate.bounding_box[5])-float(substrate.bounding_box[2])-RT)<1e-6,
         "SUBSTRATE_THICKNESS_READBACK")
    need(substrate.material_name.lower() == "rt5880_supercell", "SUBSTRATE_MATERIAL_READBACK")
    need(h.save_project(), "SAVE_PROJECT")
    return {"cell":cell,"top":top,"copper_area_mm2":copper_area,
            "copper_body_count":len(copper_solids),
            "copper_solids":[{"name": obj.name,"volume_mm3":float(obj.volume)} for obj in copper_solids],
            "objects":m.object_names,
            "boundaries":boundaries,"setup":dict(setup.props),"validation":validation,
            "aedt_version":h.desktop_class.aedt_version_id,"aedt_pid":h.desktop_class.aedt_process_id,
            "air_top_mm":topz}


def modes(h):
    fs = h.post.available_report_quantities(quantities_category="Eigen Modes")
    qs = h.post.available_report_quantities(quantities_category="Eigen Q")
    need(fs and len(fs) == len(qs), f"MODE_LIST_MISMATCH {len(fs)} {len(qs)}")
    out = []
    for i,(fq,qq) in enumerate(zip(fs,qs),1):
        fd = need(h.post.get_solution_data(expressions=fq,setup_sweep_name="EigenCell : LastAdaptive",
                                           report_category="Eigenmode"), "MODE_F")
        qd = need(h.post.get_solution_data(expressions=qq,setup_sweep_name="EigenCell : LastAdaptive",
                                           report_category="Eigenmode"), "MODE_Q")
        f = float(fd.get_expression_data()[1][0])
        q = float(qd.get_expression_data()[1][0])
        unit = fd.units_data.get(fq, "")
        fghz = f*({"Hz":1e-9,"MHz":1e-3,"GHz":1}.get(unit,1e-9))
        out.append({"mode":i,"f_ghz":fghz,"Q":q,"raw_f":f,"raw_unit":unit})
    return out


def one(cell, top, cores, tag="", phase=85, direction="across", air_mm=20.035, build_only=False, timeout_s=1500):
    RESULTS.mkdir(parents=True, exist_ok=True)
    label = f"{cell}_{top}_{direction}_{phase:g}_air{air_mm:g}{('_'+tag) if tag else ''}"
    project = WORK / "projects" / f"{label}.aedt"
    project.parent.mkdir(parents=True, exist_ok=True)
    if project.exists():
        raise RuntimeError(f"REFUSE_EXISTING_PROJECT {project}")
    start = time.monotonic()
    h = None
    state = {"status":"STARTING","cell":cell,"top":top,"controller_pid":os.getpid(),"started_utc":now()}
    status_path = RESULTS / f"{label}_status.json"
    write(status_path,state)
    finished = threading.Event()
    def heartbeat():
        while not finished.wait(60):
            state["heartbeat_utc"] = now()
            if state.get("aedt_pid"):
                try:
                    import psutil
                    state["own_child_pids"] = [p.pid for p in psutil.Process(int(state["aedt_pid"])).children(recursive=True)]
                except psutil.NoSuchProcess:
                    state["own_child_pids"] = []
            write(status_path, state)
    threading.Thread(target=heartbeat, daemon=True).start()
    try:
        h = launch(project)
        state.update(status="BUILDING",aedt_pid=h.desktop_class.aedt_process_id)
        write(status_path,state)
        built = build(h,cell,top,phase=phase,direction=direction,air_mm=air_mm)
        write(RESULTS / f"{label}_build.json", built)
        if build_only:
            state.update(status="BUILT", build_end_utc=now())
            write(status_path,state)
            return
        state.update(status="SOLVING",solve_start_utc=now())
        write(status_path,state)
        # The watchdog targets only the AEDT PID launched above.
        def timeout_own_aedt():
            import psutil
            state.update(status="TIMEOUT", error=f"{timeout_s}s solve timeout", timeout_utc=now())
            write(status_path, state)
            pid = int(h.desktop_class.aedt_process_id)
            try:
                p = psutil.Process(pid)
                targets = p.children(recursive=True) + [p]
                state["terminated_own_pids"] = [q.pid for q in targets]
                write(status_path, state)
                for q in targets:
                    if q.is_running():
                        q.terminate()
            except psutil.NoSuchProcess:
                pass
        watchdog = threading.Timer(timeout_s, timeout_own_aedt)
        watchdog.start()
        try:
            need(h.analyze_setup("EigenCell",cores=cores,use_auto_settings=True,blocking=True), "SOLVER_FAILURE")
            need(state["status"] != "TIMEOUT", "SOLVE_TIMEOUT")
        finally:
            watchdog.cancel()
        result = {"status":"COMPLETE","cell":cell,"top":top,"aedt_version":VERSION,
                  "phase_deg":phase,"direction":direction,"air_top_mm":built["air_top_mm"],
                  "aedt_pid":h.desktop_class.aedt_process_id,"modes":modes(h),
                  "seconds":time.monotonic()-start,"solve_end_utc":now(),"project":str(project)}
        try:
            result["profile_path"] = h.export_profile("EigenCell", output_file=str(RESULTS / f"{label}_profile.prof"))
            result["convergence_path"] = h.export_convergence("EigenCell", output_file=str(RESULTS / f"{label}_convergence.conv"))
        except Exception as exc:
            result["profile_export_error"] = f"{type(exc).__name__}: {exc}"
        write(RESULTS / f"{label}_modes.json",result)
        state.update(status="COMPLETE",solve_end_utc=now(),seconds=result["seconds"])
        write(status_path,state)
        print(json.dumps({"cell":cell,"top":top,"aedt_version":VERSION,"modes":result["modes"][:5],
                          "seconds":result["seconds"]},indent=2),flush=True)
    except BaseException as exc:
        state.update(status="FAILED",error=f"{type(exc).__name__}: {exc}",traceback=traceback.format_exc(),failed_utc=now())
        write(status_path,state)
        traceback.print_exc()
        raise
    finally:
        finished.set()
        if h is not None:
            h.release_desktop(close_projects=True,close_desktop=True)


if __name__ == "__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--cell",choices=CELLS,required=True)
    ap.add_argument("--top",choices=("PEC","RADIATION","PML","PML2"),required=True)
    ap.add_argument("--cores",type=int,default=8)
    ap.add_argument("--tag",default="")
    ap.add_argument("--phase",type=float,default=85)
    ap.add_argument("--direction",choices=("along","across"),default="across")
    ap.add_argument("--air-mm",type=float,default=20.035)
    ap.add_argument("--build-only",action="store_true")
    ap.add_argument("--timeout-s",type=int,default=1500)
    args=ap.parse_args()
    one(args.cell,args.top,args.cores,args.tag,args.phase,args.direction,args.air_mm,args.build_only,args.timeout_s)
