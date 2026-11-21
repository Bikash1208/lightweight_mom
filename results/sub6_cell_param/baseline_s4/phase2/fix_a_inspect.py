"""(a)(b) read-only on a COPY of the saved one-pass solution: faces near the feed, mesh OBJ, |E| on the port face. Saves nothing that changes geometry."""
import os, json, math, traceback
from pathlib import Path
import numpy as np
os.environ["ANSYSEM_ROOT261"] = r"C:\Program Files\ANSYS Inc\v261\AnsysEM"
W = Path(r"C:\work\s4_work\fix_a"); temp = W/"aedt_temp"; temp.mkdir(exist_ok=True); os.environ.update(TEMP=str(temp), TMP=str(temp))
from ansys.aedt.core import Hfss
res = {}; out = W/"fix_a.json"
def save(): out.write_text(json.dumps(res, indent=1, default=str))
h = None
try:
    h = Hfss(project=str(W/"BOARD_R150_V31_BOTH.aedt"), design="GATE_BOTH", solution_type="Modal", version="2026.1", non_graphical=True, new_desktop=True, close_on_exit=False)
    res["pid"] = h.desktop_class.aedt_process_id; m = h.modeler
    res["setups_solved"] = [(s.name, s.is_solved) for s in h.setups]
    rows = []
    for n in ["BIAS_RT5880", "SUBSTRATE", "METAL_PEC", "METAL_PEC_1", "PORT"]:
        o = m[n]; bb = o.bounding_box
        if bb[2] > -3.14 or bb[5] < -3.25 or abs(bb[0]) > 1e5: pass
        # object must reach the slab z in [-3.21,-3.175] near the axis
        if bb[5] < -3.215 or bb[2] > -3.17: continue
        for f in o.faces:
            try: c = f.center
            except Exception: continue
            r = math.hypot(c[0], c[1])
            if -3.215 <= c[2] <= -3.17 and r < 2.1:
                rows.append({"object": n, "material": o.material_name, "face_id": f.id, "center": [round(x, 4) for x in c], "area_mm2": round(f.area, 4), "is_planar": bool(f.is_planar)})
    res["faces_slab_r2p1"] = rows
    res["boundaries"] = [(b.name, b.type, getattr(b, "props", {}).get("Objects") if hasattr(b, "props") else None, getattr(b, "props", {}).get("Faces") if hasattr(b, "props") else None) for b in h.boundaries]
    res["mesh_ops"] = [mo.name for mo in h.mesh.meshoperations]; save()
    # PORT object details
    p = m["PORT"]; res["PORT"] = {"pid_obj": getattr(p, "id", None), "object_type": str(getattr(p, "object_type", None)), "bbox": p.bounding_box, "area_mm2": [f.area for f in p.faces], "expected_area": math.pi*(2.05**2-0.635**2)}
    pn = m["METAL_PEC"]; res["pin"] = {"bbox": pn.bounding_box, "faces": [(f.id, [round(x, 3) for x in f.center], round(f.area, 4)) for f in pn.faces]}
    g = m["METAL_PEC_1"]; res["ground_faces_near_axis"] = [(f.id, [round(x, 3) for x in f.center], round(f.area, 3)) for f in g.faces if math.hypot(f.center[0], f.center[1]) < 2.2 or f.area < 100]
    # integration line endpoints
    for b in h.boundaries:
        if b.name == "Port1":
            try: res["Port1_props"] = {k: str(v) for k, v in b.props.items()}
            except Exception as e: res["Port1_props"] = str(e)
    save()
    # mesh OBJ (surface mesh of all objects incl. air)
    try:
        objp = h.post.export_mesh_obj(setup="Setup1 : LastAdaptive", export_air_objects=True, on_surfaces=True); res["mesh_obj"] = str(objp)
    except Exception as e: res["mesh_obj_err"] = str(e)[:300]
    save()
    # |E| on the port face
    def fld(name, pts):
        f = W/f"{name}.fld"
        assert h.post.export_field_file(quantity="E", solution="Setup1 : LastAdaptive", output_file=str(f), sample_points=pts, export_in_si_system=False, intrinsics={"Freq": "5.68GHz", "Phase": "0deg"})
        with f.open() as s:
            s.readline(); s.readline(); raw = np.loadtxt(s)
        return raw[:, :3], np.sqrt((raw[:, 3::2]**2 + raw[:, 4::2]**2).sum(axis=1)), raw
    pts = [[r*math.cos(math.radians(a)), r*math.sin(math.radians(a)), z] for z in (-3.19, -3.176, -3.16) for r in (0.7, 1.0, 1.4, 1.8, 2.0) for a in (0, 90, 180, 270)]
    xyz, e, raw = fld("port_face", pts)
    res["E_port_face"] = {f"z={z}": {f"r={r}": round(float(np.nanmean(e[(np.isclose(xyz[:, 2], z)) & np.isclose(np.hypot(xyz[:, 0], xyz[:, 1]), r)])), 3) for r in (0.7, 1.0, 1.4, 1.8, 2.0)} for z in (-3.19, -3.176, -3.16)}
    # radial vs axial/tangential components at z=-3.176
    k = np.isclose(xyz[:, 2], -3.176) & np.isclose(np.hypot(xyz[:, 0], xyz[:, 1]), 1.4) & np.isclose(xyz[:, 1], 0) & (xyz[:, 0] > 0)
    res["E_vec_r1p4_x_axis_z-3.176_ReIm_xyz"] = [round(float(v), 3) for v in raw[k][0, 3:]] if k.any() else None
    save()
except BaseException:
    res["error"] = traceback.format_exc(); save()
finally:
    try:
        if h is not None: h.release_desktop(close_projects=True, close_desktop=True)
    except Exception: pass


