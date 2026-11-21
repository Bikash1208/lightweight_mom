"""Step 3: known-answer coax-fed patch in a radiation box, built with the SAME construction as the board (pin r0.635, port annulus r0.635-2.05 at z=-H, ground opening r2.05, PEC, RT5880, 50 ohm lumped port, radiation box), same Direct solver / 30 cores. Saves the project and exports S11."""
import os, sys, json, time, traceback, re
from pathlib import Path
os.environ["ANSYSEM_ROOT261"] = r"C:\Program Files\ANSYS Inc\v261\AnsysEM"; os.environ["ANSYSLMD_LICENSE_FILE"] = "1055@10.100.105.13"
cores = int(sys.argv[1]) if len(sys.argv) > 1 else 30; tag = sys.argv[2] if len(sys.argv) > 2 else "patch30"
W = Path(r"C:\work\s4_work") / tag; W.mkdir(exist_ok=True, parents=True)
temp = W/"aedt_temp"; temp.mkdir(exist_ok=True); os.environ.update(TEMP=str(temp), TMP=str(temp))
from ansys.aedt.core import Hfss
res = {"cores": cores}; out = W/"result.json"
def save(): out.write_text(json.dumps(res, indent=1, default=str))
H, T, PIN_R, PORT_R = 3.175, 0.035, 0.635, 2.05
PL, PW, XF, GND, PAD = 15.2, 20.9, 5.0, 60.0, 13.2
h = None
try:
    h = Hfss(project=str(W/"patch.aedt"), design="PATCH", solution_type="Modal", version="2026.1", non_graphical=True, new_desktop=True, close_on_exit=False)
    res["pid"] = h.desktop_class.aedt_process_id
    rt = h.materials.add_material("RT5880_K"); rt.permittivity = 2.2; rt.dielectric_loss_tangent = 0.0009
    m = h.modeler
    sub = m.create_box([-GND/2, -GND/2, -H], [GND, GND, H], name="SUBSTRATE", material="RT5880_K")
    gnd = m.create_box([-GND/2, -GND/2, -H-T], [GND, GND, T], name="GROUND", material="pec")
    hole = m.create_cylinder("Z", [XF, 0, -H-T], PORT_R, T, name="HOLE"); m.subtract(gnd, hole, keep_originals=False)
    pin = m.create_cylinder("Z", [XF, 0, -H], PIN_R, H, name="PIN", material="pec")
    m.subtract(sub, m.create_cylinder("Z", [XF, 0, -H], PIN_R, H, name="PINCUT"), keep_originals=False)
    patch = m.create_box([-PL/2, -PW/2, 0], [PL, PW, T], name="PATCH", material="pec")
    outer = m.create_circle("XY", [XF, 0, -H], PORT_R, name="PORT"); inner = m.create_circle("XY", [XF, 0, -H], PIN_R, name="PIN_SHEET")
    m.subtract(outer, inner, keep_originals=False)
    air = m.create_box([-GND/2-PAD, -GND/2-PAD, -H-T-PAD], [GND+2*PAD, GND+2*PAD, H+T+T+2*PAD], name="AIR", material="vacuum")
    assert h.lumped_port(assignment="PORT", integration_line=[[XF+PIN_R, 0, -H], [XF+PORT_R, 0, -H]], impedance=50, name="Port1", renormalize=True, deembed=False)
    assert h.assign_radiation_boundary_to_objects("AIR")
    st = h.create_setup("Setup1", Frequency="5.68GHz", MaximumPasses=8, MinimumPasses=2, MaxDeltaS=0.02)
    st.props["BasisOrder"] = 1; st.props["DrivenSolverType"] = "Direct"; assert st.update()
    res["setup"] = {k: str(st.props.get(k)) for k in ("Frequency", "MaximumPasses", "MinimumPasses", "MaxDeltaS", "DrivenSolverType", "BasisOrder")}
    h.create_linear_count_sweep("Setup1", "GHz", 5.0, 6.4, 57, name="Sweep1", sweep_type="Interpolating", save_fields=False)
    h.save_project(); res["saved_before"] = True; save()
    t0 = time.monotonic(); ok = h.analyze_setup("Setup1", cores=cores, use_auto_settings=False, blocking=True)
    res["wall_s"] = round(time.monotonic()-t0, 1); res["analyze_ok"] = bool(ok); save()
    h.save_project()
    prof = Path(h.export_profile("Setup1", output_file=str(W/"solve.prof"))); txt = prof.read_text(errors="replace")
    res["delta_s"] = re.findall(r"Max Mag\. Delta S:\s*([\d.]+)", txt); res["tets"] = re.findall(r"Tetrahedra:\s*(\d+)", txt)
    res["solve_types"] = sorted(set(re.findall(r"Type:\s*(\w+), Cores", txt)))
    def get(expr, sw):
        d = h.post.get_solution_data(expressions=expr, setup_sweep_name=sw)
        return list(d.primary_sweep_values), [float(x) for x in d.get_expression_data(expression=expr, formula="real")[1]]
    f, s11 = get("dB(S(Port1,Port1))", "Setup1 : Sweep1"); _, ph = get("ang_deg(S(Port1,Port1))", "Setup1 : Sweep1")
    k = min(range(len(s11)), key=lambda i: s11[i])
    res["sweep"] = {"min_s11_db": s11[k], "f_at_min_ghz": f[k], "s11_at_5p68": s11[min(range(len(f)), key=lambda i: abs(f[i]-5.68))], "max_s11_db": max(s11), "phase_at_min_deg": ph[k]}
    f2, s2 = get("dB(S(Port1,Port1))", "Setup1 : LastAdaptive"); res["lastadaptive_s11_5p68_db"] = s2[0]
    res["setups_solved"] = [(s.name, s.is_solved) for s in h.setups]; save()
except BaseException:
    res["error"] = traceback.format_exc(); save()
finally:
    try:
        if h is not None: h.release_desktop(close_projects=True, close_desktop=True)
    except Exception: pass
