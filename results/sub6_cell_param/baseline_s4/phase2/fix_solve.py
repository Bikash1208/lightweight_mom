"""Feed-fix one-pass test on a COPY of the BOTH archive (min=max=1, 16 cores). mode c1: length mesh ops <=0.2 mm on PORT sheet, pin, ground-hole wall.
mode c2: (on top of c1 if arg 'c1+c2') redraw annulus as new sheet, reassign lumped port by object, integration line by coordinates. Saves + exports before release."""
import os, sys, json, math, re, time, traceback, shutil
from pathlib import Path
import numpy as np
os.environ["ANSYSEM_ROOT261"] = r"C:\Program Files\ANSYS Inc\v261\AnsysEM"; os.environ["ANSYSLMD_LICENSE_FILE"] = "1055@10.100.105.13"
mode = sys.argv[1]; tag = sys.argv[2]; cores = 16
W = Path(r"C:\work\s4_work") / tag; W.mkdir(exist_ok=True, parents=True)
shutil.copy(r"C:\work\s4_work\both_step3\BOARD_R150_V31_BOTH.aedt", W/"BOARD_R150_V31_BOTH.aedt")
temp = W/"aedt_temp"; temp.mkdir(exist_ok=True); os.environ.update(TEMP=str(temp), TMP=str(temp))
from ansys.aedt.core import Hfss
res = {"mode": mode, "cores": cores}; out = W/"result.json"
def save(): out.write_text(json.dumps(res, indent=1, default=str))
h = None
try:
    h = Hfss(project=str(W/"BOARD_R150_V31_BOTH.aedt"), design="GATE_BOTH", solution_type="Modal", version="2026.1", non_graphical=True, new_desktop=True, close_on_exit=False)
    res["pid"] = h.desktop_class.aedt_process_id; m = h.modeler
    st = h.setups[0]; res["setup_before"] = {k: str(st.props.get(k)) for k in ("MinimumPasses", "MaximumPasses", "DrivenSolverType")}
    assert str(st.props["MinimumPasses"]) == "1" and str(st.props["MaximumPasses"]) == "1"
    if "c1" in mode:
        ops = []
        ops.append(h.mesh.assign_length_mesh(["PORT"], inside_selection=False, maximum_length=0.2, maximum_elements=None, name="ML_PORT"))
        ops.append(h.mesh.assign_length_mesh(["METAL_PEC"], inside_selection=False, maximum_length=0.2, maximum_elements=None, name="ML_PIN"))
        hole = [f for f in m["METAL_PEC_1"].faces if abs(math.hypot(*f.center[:2]) - 2.05) < 1e-3 and abs(f.center[2] + 3.1925) < 1e-3]
        res["hole_face_ids"] = [f.id for f in hole]
        ops.append(h.mesh.assign_length_mesh([f.id for f in hole], inside_selection=False, maximum_length=0.2, maximum_elements=None, name="ML_HOLE"))
        res["mesh_ops"] = [(o.name if o else None) for o in ops]; assert all(ops), res["mesh_ops"]
    if "c2" in mode:
        old = [b for b in h.boundaries if b.name == "Port1"][0]; old.delete()
        m.delete("PORT")
        outer = m.create_circle("XY", [0, 0, -3.175], 2.05, name="PORT_NEW"); inner = m.create_circle("XY", [0, 0, -3.175], 0.635, name="PORT_IN")
        m.subtract(outer, inner, keep_originals=False)
        ok = h.lumped_port(assignment="PORT_NEW", integration_line=[[0.635, 0, -3.175], [2.05, 0, -3.175]], impedance=50, name="Port1", renormalize=True, deembed=False)
        res["port_reassigned"] = bool(ok); res["new_port_area"] = [f.area for f in m["PORT_NEW"].faces]
    res["boundaries"] = [(b.name, b.type) for b in h.boundaries]; res["objects_count"] = len(m.object_names); save()
    h.save_project(); save()
    t0 = time.monotonic(); ok = h.analyze_setup("Setup1", cores=cores, use_auto_settings=False, blocking=True)
    res["wall_s"] = round(time.monotonic()-t0, 1); res["analyze_ok"] = bool(ok); save()
    h.save_project(); res["saved"] = True; res["results_dir_exists"] = (W/"BOARD_R150_V31_BOTH.aedtresults").exists(); save()
    prof = Path(h.export_profile("Setup1", output_file=str(W/"solve.prof"))); txt = prof.read_text(errors="replace")
    res["tets"] = re.findall(r"Tetrahedra:\s*(\d+)", txt); res["max_mem_gb"] = re.findall(r"Max memory/process:\s*([\d.]+)\s*GB", txt); res["solve_types"] = sorted(set(re.findall(r"Type:\s*(\w+), Cores", txt))); save()
    sw = "Setup1 : LastAdaptive"
    def val(expr):
        try:
            d = h.post.get_solution_data(expressions=expr, setup_sweep_name=sw)
            return float(d.get_expression_data(expression=expr, formula="real")[1][0]) if d else None
        except Exception as e: return "ERR " + str(e)[:100]
    for e in ("dB(S(Port1,Port1))", "mag(S(Port1,Port1))", "ang_deg(S(Port1,Port1))", "re(Z(Port1,Port1))", "im(Z(Port1,Port1))"): res[e] = val(e)
    if isinstance(res["mag(S(Port1,Port1))"], float): res["accepted_power"] = 1 - res["mag(S(Port1,Port1))"]**2
    save()
    f = W/"port_face.fld"; pts = [[r*math.cos(math.radians(a)), r*math.sin(math.radians(a)), -3.176] for r in (0.7, 1.0, 1.4, 1.8) for a in (0, 90, 180, 270)]
    assert h.post.export_field_file(quantity="E", solution=sw, output_file=str(f), sample_points=pts, export_in_si_system=False, intrinsics={"Freq": "5.68GHz", "Phase": "0deg"})
    with f.open() as s:
        s.readline(); s.readline(); raw = np.loadtxt(s)
    e = np.sqrt((raw[:, 3::2]**2 + raw[:, 4::2]**2).sum(axis=1)); res["E_port_face_by_r"] = {r: round(float(e[i*4:(i+1)*4].mean()), 2) for i, r in enumerate((0.7, 1.0, 1.4, 1.8))}
    res["aedt_warnings"] = [str(x) for x in h.odesktop.GetMessages("BOARD_R150_V31_BOTH", "GATE_BOTH", 1)][:6]; save()
    h.save_project(); res["results_dir_exists_final"] = (W/"BOARD_R150_V31_BOTH.aedtresults").exists(); save()
except BaseException:
    res["error"] = traceback.format_exc(); save()
finally:
    try:
        if h is not None: h.save_project(); h.release_desktop(close_projects=True, close_desktop=True)
    except Exception: pass
