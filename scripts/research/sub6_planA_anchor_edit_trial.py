"""Test a single bias-layer thickness edit in a private copy of the imported STEP."""
import json
import os
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

from ansys.aedt.core import Hfss


def now():
    return datetime.now(timezone.utc).isoformat()


def snap(obj):
    return [round(float(x), 6) for x in obj.bounding_box]


def main():
    src, out = [Path(x).resolve() for x in sys.argv[1:3]]
    out.mkdir(parents=True, exist_ok=True)
    dst = out / "anchor_edit_trial.aedt"
    shutil.copy2(src, dst)
    report = out / "edit_trial.json"
    d = {"started_utc": now(), "source_project": str(src), "own_python_pid": os.getpid()}
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
        metals = [m[n] for n in m.object_names if n.startswith("METAL_PEC")]
        vias = [o for o in metals if abs(float(o.bounding_box[2]) + 4.045) < 0.0001]
        d["via_body_count"] = len(vias)
        v = vias[0]
        d["via_body"] = v.name
        d["via_before"] = snap(v)
        d["substrate_before"] = snap(m["SUBSTRATE_FR4"])
        vf = [{"id": f.id, "center": f.center, "normal": f.normal,
               "area": f.area} for f in v.faces]
        d["via_faces"] = vf
        sf = []
        for f in m["SUBSTRATE_FR4"].faces:
            c = f.center
            if c and abs(float(c[2]) + 4.01) < 0.0001:
                sf.append({"id": f.id, "center": c, "normal": f.normal,
                           "area": f.area})
        d["substrate_bottom_faces"] = sf
        report.write_text(json.dumps(d, indent=2))
        if len(vias) != 480 or len(sf) != 1:
            d["gate"] = "STOP_COUNT_OR_BOTTOM_FACE_MISMATCH"
            return
        # Move only the bias substrate bottom face from -4.010 to -3.718 mm.
        # The copper test below is on one representative via body only.
        bottom = sf[0]
        if abs(float(bottom["normal"][2]) + 1) > 0.01:
            d["gate"] = "STOP_BOTTOM_FACE_NORMAL"
            return
        d["substrate_move_ok"] = bool(m.move_face([bottom["id"]], offset=-0.292))
        d["substrate_after"] = snap(m["SUBSTRATE_FR4"])
        report.write_text(json.dumps(d, indent=2))
        if abs(d["substrate_after"][2] + 3.718) > 0.002:
            d["gate"] = "STOP_SUBSTRATE_READBACK"
            return
        d["gate"] = "SUBSTRATE_EDIT_PROVED_VIA_EDIT_PENDING"
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
