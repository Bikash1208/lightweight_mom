"""Test one via/pad translation while preserving the top copper in a private AEDT copy."""
import json
import os
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

from ansys.aedt.core import Hfss


def now():
    return datetime.now(timezone.utc).isoformat()


def main():
    src, out = [Path(x).resolve() for x in sys.argv[1:3]]
    out.mkdir(parents=True, exist_ok=True)
    dst = out / "anchor_via_trial.aedt"
    shutil.copy2(src, dst)
    report = out / "via_trial.json"
    d = {"started_utc": now(), "own_python_pid": os.getpid()}
    report.write_text(json.dumps(d, indent=2))
    h = None
    try:
        h = Hfss(project=str(dst), design="AnchorPreflight", solution_type="Modal",
                 version="2025.2", non_graphical=True, new_desktop=True,
                 close_on_exit=False)
        d["aedt_version"] = str(h.desktop_class.aedt_version_id)
        d["own_aedt_pid"] = h.desktop_class.aedt_process_id
        report.write_text(json.dumps(d, indent=2))
        m = h.modeler
        name = json.loads((src.parent.parent / "edit_trial" / "edit_trial.json").read_text())["via_body"]
        o = m[name]
        d["body"] = name
        d["bbox_before"] = [float(x) for x in o.bounding_box]
        faces = []
        for f in o.faces:
            c = f.center
            faces.append({"id": f.id, "z": float(c[2]), "area": f.area})
        pad_faces = [f["id"] for f in faces if
                     (abs(f["z"] + 4.045) < 1e-5 and abs(f["area"] - 0.2827433388) < 0.01) or
                     (abs(f["z"] + 4.010) < 1e-5 and abs(f["area"] - 0.2120575041) < 0.01)]
        d["faces_before"] = faces
        d["selected_faces"] = pad_faces
        report.write_text(json.dumps(d, indent=2))
        if len(pad_faces) != 2:
            d["gate"] = "STOP_PAD_FACE_IDENTIFICATION"
            return
        m.oeditor.MoveFaces(
            ["NAME:Selections", "Selections:=", name, "NewPartsModelFlag:=", "Model"],
            ["NAME:Parameters", ["NAME:MoveFacesParameters",
              "MoveAlongNormalFlag:=", False, "OffsetDistance:=", "0mm",
              "MoveVectorX:=", "0mm", "MoveVectorY:=", "0mm", "MoveVectorZ:=", "0.292mm",
              "FacesToMove:=", pad_faces]])
        d["bbox_after"] = [float(x) for x in m[name].bounding_box]
        d["faces_after"] = [{"id": f.id, "z": float(f.center[2]), "area": f.area}
                            for f in m[name].faces]
        d["gate"] = "PASS_ONE_VIA" if (
            abs(d["bbox_after"][2] + 3.753) < 0.002 and
            abs(d["bbox_after"][5] - 0.035) < 0.002) else "STOP_VIA_READBACK"
        h.save_project()
    except Exception as exc:
        d["error"] = repr(exc)
        d["gate"] = "STOP_EXCEPTION"
    finally:
        d["ended_utc"] = now()
        report.write_text(json.dumps(d, indent=2))
        if h is not None:
            h.release_desktop(close_projects=True, close_desktop=True)


if __name__ == "__main__":
    main()
