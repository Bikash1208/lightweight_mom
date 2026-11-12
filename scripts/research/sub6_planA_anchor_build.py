"""Build one Stage 1 anchor variant from the audited imported STEP project.

Always launches its own AEDT 2025.2 process. No solve is performed by this script.
"""
import argparse
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path

from ansys.aedt.core import Hfss

THICK = {"V": 0.8, "T8": 0.8, "T5": 0.508, "T2": 0.254, "T1": 0.127}
RSXS = {"OFF": (3.885, -58.70), "ON": (0.28, 5.11)}


def now():
    return datetime.now(timezone.utc).isoformat()


def box(o):
    return [round(float(x), 6) for x in o.bounding_box]


def require(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--base", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--variant", choices=THICK, required=True)
    p.add_argument("--state", choices=RSXS, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    run = f"ANCHOR_{a.variant}_{a.state}"
    dst = a.out / f"{run}.aedt"
    report = a.out / f"{run}_build.json"
    shutil.copy2(a.base, dst)
    d = {"run_id": run, "started_utc": now(), "own_python_pid": os.getpid(),
         "target_bias_thickness_mm": THICK[a.variant], "state": a.state,
         "status": "BUILDING"}
    report.write_text(json.dumps(d, indent=2))
    h = None
    try:
        h = Hfss(project=str(dst), design="AnchorPreflight", solution_type="Modal",
                 version="2025.2", non_graphical=True, new_desktop=True,
                 close_on_exit=False)
        if os.environ.get("AEDT_PRIVATE_TEMP"):
            h.desktop_class.temp_directory = os.environ["AEDT_PRIVATE_TEMP"]
            require(Path(h.desktop_class.temp_directory).resolve() == Path(os.environ["AEDT_PRIVATE_TEMP"]).resolve(),
                    "PRIVATE_TEMP_READBACK")
        d["aedt_version"] = str(h.desktop_class.aedt_version_id)
        d["own_aedt_pid"] = h.desktop_class.aedt_process_id
        report.write_text(json.dumps(d, indent=2))
        require(d["aedt_version"] == "2025.2", "AEDT_VERSION")
        m = h.modeler
        names = list(m.object_names)
        metal = [n for n in names if n.startswith("METAL_PEC")]
        sheets = [n for n in names if n.startswith("D") and n not in {"SUBSTRATE", "SUBSTRATE_FR4"}]
        require(len(metal) == 1000 and len(sheets) == 240 and len(names) == 1244,
                f"IMPORT_COUNTS metal={len(metal)} diode={len(sheets)} all={len(names)}")
        require(all(n in names for n in ("AIR", "SUBSTRATE", "SUBSTRATE_FR4", "PORT")), "NAMED_PARTS")
        via = [n for n in metal if abs(box(m[n])[2] + 4.045) < 1e-4]
        require(len(via) == 480, f"VIA_COUNT_{len(via)}")
        d["counts_before"] = {"metal_bodies": len(metal), "diode_sheets": len(sheets),
                              "via_bodies": len(via), "substrates": 2, "total": len(names)}
        d["bias_before"] = box(m["SUBSTRATE_FR4"])
        d["top_substrate_before"] = box(m["SUBSTRATE"])
        require(abs(d["bias_before"][2] + 4.010) < 1e-4 and
                abs(d["bias_before"][5] + 3.175) < 1e-4, "BIAS_ORIGINAL_BBOX")
        shift = round(0.8 - THICK[a.variant], 6)
        if shift:
            bottom = []
            for f in m["SUBSTRATE_FR4"].faces:
                c = f.center
                if c and abs(float(c[2]) + 4.01) < 1e-4:
                    bottom.append(f.id)
            require(len(bottom) == 1, f"SUBSTRATE_BOTTOM_FACES_{len(bottom)}")
            require(m.move_face(bottom, offset=-shift), "BIAS_FACE_MOVE")
            require(abs(box(m["SUBSTRATE_FR4"])[2] - (-4.010 + shift)) < 1e-4,
                    "BIAS_MOVE_READBACK")
            for i, name in enumerate(via, 1):
                o = m[name]
                faces = []
                for f in o.faces:
                    c = f.center
                    if not c:
                        continue
                    z, area = float(c[2]), float(f.area)
                    if ((abs(z + 4.045) < 1e-4 and abs(area - 0.2827433388) < 0.01)
                            or (abs(z + 4.010) < 1e-4 and abs(area - 0.2120575041) < 0.01)):
                        faces.append(f.id)
                require(len(faces) == 2, f"PAD_FACE_IDENTIFICATION {name}: {faces}")
                m.oeditor.MoveFaces(
                    ["NAME:Selections", "Selections:=", name, "NewPartsModelFlag:=", "Model"],
                    ["NAME:Parameters", ["NAME:MoveFacesParameters",
                      "MoveAlongNormalFlag:=", False, "OffsetDistance:=", "0mm",
                      "MoveVectorX:=", "0mm", "MoveVectorY:=", "0mm",
                      "MoveVectorZ:=", f"{shift}mm", "FacesToMove:=", faces]])
                bb = box(m[name])
                require(abs(bb[2] - (-4.045 + shift)) < 1e-4 and abs(bb[5] - 0.035) < 1e-4,
                        f"VIA_MOVE_READBACK {name}: {bb}")
                if i % 40 == 0:
                    d["vias_moved"] = i
                    report.write_text(json.dumps(d, indent=2))
        require(box(m["SUBSTRATE"]) == d["top_substrate_before"], "TOP_SUBSTRATE_CHANGED")
        require(len(m.object_names) == 1244, "OBJECT_COUNT_CHANGED")
        d["bias_after"] = box(m["SUBSTRATE_FR4"])
        d["via_zmin_after"] = sorted({box(m[n])[2] for n in via})
        require(d["via_zmin_after"] == [round(-4.045 + shift, 6)], "VIA_READBACK_DISTRIBUTION")
        d["geometry_gate"] = "PASS"
        report.write_text(json.dumps(d, indent=2))
        # Assign the standard materials only after the geometry gate.
        rt = h.materials.add_material("RT5880_STAGE1")
        rt.permittivity = 2.2
        rt.dielectric_loss_tangent = 0.0009
        require(h.assign_material(metal, "pec"), "PEC_ASSIGNMENT")
        require(h.assign_material("SUBSTRATE", "RT5880_STAGE1"), "TOP_SUBSTRATE_MATERIAL")
        bias_mat = "vacuum" if a.variant == "V" else "RT5880_STAGE1"
        require(h.assign_material("SUBSTRATE_FR4", bias_mat), "BIAS_MATERIAL")
        require(h.assign_material("AIR", "vacuum"), "AIR_MATERIAL")
        d["material_readback"] = {n: m[n].material_name for n in ("SUBSTRATE", "SUBSTRATE_FR4", "AIR")}
        require(d["material_readback"]["SUBSTRATE_FR4"].lower() == bias_mat.lower(),
                "BIAS_MATERIAL_READBACK")
        require(h.assign_radiation_boundary_to_objects("AIR", name="RAD_AIR"), "RADIATION")
        rs, xs = RSXS[a.state]
        imp = h.assign_impedance_to_sheet(sheets, name="BAR90_02LRH",
                                           resistance=rs, reactance=xs)
        require(imp, "DIODE_IMPEDANCE")
        d["impedance_readback"] = dict(imp.props)
        require(len(d["impedance_readback"].get("Objects", [])) == 240, "IMPEDANCE_SHEET_COUNT")
        port = h.lumped_port("PORT", integration_line=[[0.635, 0, -3.175], [2.05, 0, -3.175]],
                             impedance=50, name="Port1")
        require(port, "LUMPED_PORT")
        d["port_readback"] = dict(port.props)
        setup = h.create_setup("Setup1")
        setup.props["Frequency"] = "5.8GHz"
        setup.props["MaximumPasses"] = 4
        setup.props["MaxDeltaS"] = 0.02
        require(setup.update(), "SETUP_UPDATE")
        d["setup_readback"] = {k: setup.props.get(k) for k in
                               ("Frequency", "MaximumPasses", "MaxDeltaS")}
        sphere = h.insert_infinite_sphere(name="Global", theta_start=0, theta_stop=180,
                                          theta_step=2, phi_start=-180, phi_stop=180, phi_step=2)
        require(sphere, "INFINITE_SPHERE")
        d["sphere_readback"] = dict(sphere.props)
        h.save_project()
        d["status"] = "BUILT_UNSOLVED"
    except Exception as exc:
        d["status"] = "STOPPED"
        d["error"] = repr(exc)
    finally:
        d["ended_utc"] = now()
        report.write_text(json.dumps(d, indent=2, default=str))
        if h is not None:
            h.release_desktop(close_projects=True, close_desktop=True)


if __name__ == "__main__":
    main()
