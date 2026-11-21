"""Step 2 (mesh only, no field solve) on a COPY of the stored BOTH archive: geometry/assignment checks, GenerateMesh, tet count, AEDT messages. Saves + exports before release."""
import os, sys, json, traceback, time
from pathlib import Path
os.environ["ANSYSEM_ROOT261"] = r"C:\Program Files\ANSYS Inc\v261\AnsysEM"; os.environ["ANSYSLMD_LICENSE_FILE"] = "1055@10.100.105.13"
W = Path(r"C:\work\s4_work\both_step2"); temp = W/"aedt_temp"; temp.mkdir(exist_ok=True); os.environ.update(TEMP=str(temp), TMP=str(temp))
from ansys.aedt.core import Hfss
res = {}; out = W/"step2.json"
def save(): out.write_text(json.dumps(res, indent=1, default=str))
h = None
try:
    h = Hfss(project=str(W/"BOARD_R150_V31_BOTH.aedt"), design="GATE_BOTH", solution_type="Modal", version="2026.1", non_graphical=True, new_desktop=True, close_on_exit=False)
    m = h.modeler; res["pid"] = h.desktop_class.aedt_process_id
    names = list(m.object_names); res["n_objects"] = len(names)
    def bb(n): return [round(x, 4) for x in m[n].bounding_box]
    res["non_metal"] = {n: (m[n].material_name, bb(n)) for n in names if not n.startswith("METAL_PEC")}
    mp = [n for n in names if n.startswith("METAL_PEC")]; res["metal_objects"] = {n: (m[n].material_name, bb(n), round(m[n].volume, 4)) for n in mp[:3]}
    res["metal_materials"] = sorted(set(m[n].material_name for n in mp))
    # the pin is one of the metal solids: find solids touching r<=0.7 around x=y=0
    pins = [n for n in mp if abs(m[n].bounding_box[0]+0.635) < 1e-3 and abs(m[n].bounding_box[3]-0.635) < 1e-3]
    res["pin_candidates"] = {n: (bb(n), m[n].material_name) for n in pins}
    # any metal in z>=0 within 6 mm of the axis (centre cell is the feed, expected empty)
    near = [(n, bb(n)) for n in mp if m[n].bounding_box[2] > -0.01 and abs(m[n].bounding_box[0]) < 6 and abs(m[n].bounding_box[1]) < 6 and m[n].bounding_box[3] < 6]
    res["metal_near_axis_above_z0"] = near
    res["port"] = [(p.name, p.props.get("Impedance"), p.props.get("DeembedDist") if hasattr(p, "props") else None) for p in h.boundaries if p.type in ("Port", "Lumped Port") or "Port" in p.name]
    res["boundaries"] = [(b.name, b.type) for b in h.boundaries]
    res["mesh_ops"] = [mo.name for mo in h.mesh.meshoperations]
    res["setups"] = [(s.name, s.props.get("DrivenSolverType")) for s in h.setups]; save()
    t0 = time.monotonic()
    ok = h.mesh.generate_mesh("Setup1"); res["generate_mesh_ok"] = ok; res["mesh_wall_s"] = round(time.monotonic()-t0, 1); save()
    ms = h.export_mesh_stats("Setup1", mesh_path=str(W/"mesh_stats.ms"))
    res["mesh_stats_file"] = str(ms); res["mesh_stats_text"] = Path(ms).read_text(errors="replace")[:1500] if ms else None; save()
    msgs = {}
    for sev in (0, 1, 2, 3):
        try: msgs[sev] = list(h.odesktop.GetMessages("BOARD_R150_V31_BOTH", "GATE_BOTH", sev))
        except Exception as e: msgs[sev] = "ERR " + str(e)[:100]
    res["aedt_messages_by_severity"] = msgs; save()
except BaseException:
    res["error"] = traceback.format_exc(); save()
finally:
    try:
        if h is not None:
            h.save_project(); res["saved"] = True; res["results_dir_exists"] = (W/"BOARD_R150_V31_BOTH.aedtresults").exists(); save()
            h.release_desktop(close_projects=True, close_desktop=True)
    except Exception: pass
