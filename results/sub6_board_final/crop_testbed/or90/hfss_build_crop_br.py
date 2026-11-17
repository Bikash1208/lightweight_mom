"""Build BOARD_FINAL_ELS.aedt (UNSOLVED): four designs OFF_ELS10_ABS / ON_ELS_ABS (idealised edge absorber) and OFF_ELS10 / ON_ELS (references).  AEDT 2025.1, validate_simple only,
automatic meshing only (no mesh operations).  usage: hfss_build_final.py <step> <project.aedt> <readback.json> <cell_radius_mm> <board_radius_mm> <n_metal_objects> <n_diodes>"""
import hashlib, json, math, os, sys, time
from pathlib import Path
os.environ["ANSYSEM_ROOT251"] = r"C:\Program Files\ANSYS Inc\v251\AnsysEM"
from ansys.aedt.core import Hfss

STEP, PROJ, REP = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
R_IN, R_OUT, N_METAL = float(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6])
N_D = int(sys.argv[7]); MODE = sys.argv[8]
ABS_H = 30.0
temp = Path(r"C:/Users/Bikash/Desktop/sub6_board_hfss/aedt_temp_bar90/board_crop_build_br_" + sys.argv[8]); temp.mkdir(parents=True, exist_ok=True)
os.environ.update(TEMP=str(temp), TMP=str(temp))
PROJ.parent.mkdir(parents=True, exist_ok=True)
SHEETS = {"OFF": (21.6, -297.5), "ON": (2.8, 14.6)}
DESIGNS = [("OFF_ELS10_ABS", "OFF", True, "-10 V reverse bias, idealised edge absorber; run first"), ("ON_ELS_ABS", "ON", True, "10 mA forward, idealised edge absorber; run second"),
           ("OFF_ELS10", "OFF", False, "-10 V, no absorber (reference)"), ("ON_ELS", "ON", False, "10 mA, no absorber (reference)")]
rep = {"step": str(STEP), "step_sha256": hashlib.sha256(STEP.read_bytes()).hexdigest(), "step_bytes": STEP.stat().st_size, "project": str(PROJ), "started": time.strftime("%Y-%m-%d %H:%M:%S"), "errors": {},
       "absorber": dict(inner_radius_mm=R_IN, outer_radius_mm=R_OUT, height_mm=ABS_H, z_from_mm=0.0, z_to_mm=ABS_H + 0.035, er="1.5 - j1.0", ur="1.5 - j1.0", idealised=True)}


def need(c, m):
    if not c:
        raise RuntimeError(m)
    return c


h = Hfss(project=str(PROJ), design="OFF_ELS10", solution_type="Modal", version="2025.1", non_graphical=True, new_desktop=True, close_on_exit=False)
rep["aedt_pid"] = h.desktop_class.aedt_process_id
print("OWN_AEDT_PID", rep["aedt_pid"], flush=True)
try:
    need(str(h.desktop_class.aedt_version_id) == "2025.1", "VERSION")
    m = h.modeler; m.model_units = "mm"
    rt = h.materials.add_material("RT5880_V2"); rt.permittivity = 2.2; rt.dielectric_loss_tangent = 0.0009
    ab = h.materials.add_material("ABSORBER_15_10")
    ab.permittivity = 1.5; ab.dielectric_loss_tangent = 1.0 / 1.5; ab.permeability = 1.5; ab.magnetic_loss_tangent = 1.0 / 1.5
    need(m.import_3d_cad(str(STEP)), "IMPORT")
    names = list(m.object_names)
    metal = [n for n in names if n.startswith("METAL_PEC")]
    diodes = sorted([n for n in names if len(n) > 1 and n[0] == "D" and n[1:].isdigit()], key=lambda n: int(n[1:]))
    need(len(metal) > 0 and len(diodes) == N_D, f"COUNTS metal={len(metal)} diodes={len(diodes)}"); print("COUNTS", len(metal), len(diodes), flush=True)
    for n in metal:
        m[n].material_name = "pec"
    m["SUBSTRATE"].material_name = "RT5880_V2"; m["BIAS_RT5880"].material_name = "RT5880_V2"; m["AIR"].material_name = "vacuum"
    h["Rs"] = f"{SHEETS['OFF'][0]}ohm"; h["Xs"] = f"{SHEETS['OFF'][1]}ohm"
    if diodes: need(h.assign_impedance_to_sheet(diodes, name="BAR90_ELS", resistance="Rs", reactance="Xs"), "IMPEDANCE")
    brs = [n for n in names if n.startswith("BR") and "x" in n]
    print("BRIDGE_SHEETS", len(brs), flush=True)
    need(len(brs) > 0 and h.assign_impedance_to_sheet(brs, name="BRIDGE_SHORT", resistance="0.001ohm", reactance="0ohm"), "BRIDGE")
    need(h.lumped_port(assignment="PORT", integration_line=[[0.635, 0, -3.175], [2.05, 0, -3.175]], impedance=50, name="Port1", renormalize=True, deembed=False), "PORT")
    need(h.assign_radiation_boundary_to_objects("AIR"), "RADIATION")
    setup = need(h.create_setup("Setup1", Frequency="5.8GHz", MaximumPasses=20, MaxDeltaS=0.02, MinimumPasses=6), "SETUP")
    setup.props["MinimumPasses"] = 6; setup.props["MaximumPasses"] = 20; setup.props["MaxDeltaS"] = 0.02; setup.props["BasisOrder"] = 1
    need(setup.update(), "SETUP_UPDATE")
    need(h.insert_infinite_sphere(theta_start=0, theta_stop=180, theta_step=2, phi_start=0, phi_stop=360, phi_step=2, name="Sphere1"), "SPHERE")
    need(h.insert_infinite_sphere(theta_start=90, theta_stop=180, theta_step=2, phi_start=0, phi_stop=360, phi_step=2, name="Lower"), "SPHERE_LOWER")
    need(h.save_project(), "SAVE_FIRST")
    def add_absorber():
        o = m.create_cylinder("Z", [0, 0, 0], R_OUT, ABS_H + 0.035, name="ABS_OUT", material="ABSORBER_15_10")
        i = m.create_cylinder("Z", [0, 0, -1], R_IN, ABS_H + 2.0, name="ABS_IN", material="vacuum")
        need(o and i, "ABS_CYL"); need(m.subtract("ABS_OUT", "ABS_IN", keep_originals=False), "ABS_SUBTRACT")
        m["ABS_OUT"].name = "ABSORBER"; m["ABSORBER"].material_name = "ABSORBER_15_10"
    names = []
    if MODE == "sw":
        h.rename_design("OFF_ELS10_ABS") if False else None
        add_absorber(); h.save_project()
        h.set_active_design("OFF_ELS10"); need(h.duplicate_design("ON_ELS"), "DUP_ON"); h.set_active_design("ON_ELS")
        h["Rs"] = f"{SHEETS['ON'][0]}ohm"; h["Xs"] = f"{SHEETS['ON'][1]}ohm"
        names = ["OFF_ELS10", "ON_ELS"]
    else:
        add_absorber(); names = ["OFF_ELS10"]
    need(h.save_project(), "SAVE_DESIGNS")
    rep["designs"] = {}
    for name in names:
        h.set_active_design(name)
        d = {"mesh_ops": [x.name for x in h.mesh.meshoperations], "n_objects": len(m.object_names)}
        if diodes:
            d["Rs_read"] = h.get_evaluated_value("Rs"); d["Xs_read"] = h.get_evaluated_value("Xs")
            b = next(x for x in h.boundaries if x.name == "BAR90_ELS"); d["n_sheets"] = len(list(b.props.get("Objects", [])))
        s = next(x for x in h.setups if x.name == "Setup1")
        d["setup"] = {k: str(s.props.get(k)) for k in ("Frequency", "MaximumPasses", "MinimumPasses", "MaxDeltaS", "BasisOrder", "DrivenSolverType")}
        d["spheres"] = [f.name for f in h.field_setups]
        vf = REP.parent / f"validate_simple_{name}.log"
        d["validate_simple"] = h.validate_simple(str(vf)); rep["designs"][name] = d
        print("READBACK", name, d, flush=True)
    need(h.save_project(), "SAVE_FINAL")
    rep["design_list"] = list(h.design_list)
    mat = h.materials.material_keys["absorber_15_10"] if "absorber_15_10" in h.materials.material_keys else None
    if mat is not None:
        rep["absorber_material_readback"] = dict(permittivity=str(mat.permittivity), dielectric_loss_tangent=str(mat.dielectric_loss_tangent), permeability=str(mat.permeability), magnetic_loss_tangent=str(mat.magnetic_loss_tangent))
    REP.write_text(json.dumps(rep, indent=1, default=str))
    print("BUILT_FINAL", rep["design_list"], flush=True)
finally:
    try:
        REP.write_text(json.dumps(rep, indent=1, default=str))
    except Exception:
        pass
    h.release_desktop(close_projects=True, close_desktop=True)
print("DONE", flush=True)
