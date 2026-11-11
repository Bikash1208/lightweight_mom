"""Private AEDT import audit for the N240 Stage 1 anchor STEP."""
import json
import os
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from ansys.aedt.core import Hfss


def stamp():
    return datetime.now(timezone.utc).isoformat()


def main():
    step = Path(sys.argv[1]).resolve()
    out = Path(sys.argv[2]).resolve()
    out.mkdir(parents=True, exist_ok=True)
    project = out / "anchor_import_preflight.aedt"
    report = out / "preflight.json"
    data = {"started_utc": stamp(), "step": str(step), "own_python_pid": os.getpid()}
    report.write_text(json.dumps(data, indent=2))
    h = None
    try:
        h = Hfss(project=str(project), design="AnchorPreflight", solution_type="Modal",
                 version="2025.2", non_graphical=True, new_desktop=True,
                 close_on_exit=False)
        data["aedt_version"] = str(h.desktop_class.aedt_version_id)
        data["own_aedt_pid"] = h.desktop_class.aedt_process_id
        report.write_text(json.dumps(data, indent=2))
        m = h.modeler
        before = set(m.object_names)
        data["import_ok"] = bool(m.import_3d_cad(step, create_group=False,
                                                  import_free_surfaces=True,
                                                  separate_disjoints_lumped_object=False))
        names = [n for n in m.object_names if n not in before]
        data["object_count"] = len(names)
        data["names"] = names
        data["name_prefix_counts"] = dict(Counter("D_sheet" if n.startswith("D") else n for n in names))
        data["object_summary"] = {}
        for n in names:
            if not n.startswith("D") or n in {"SUBSTRATE", "SUBSTRATE_FR4"}:
                o = m[n]
                data["object_summary"][n] = {
                    "type": str(o.object_type),
                    "bbox": [float(v) for v in o.bounding_box],
                    "face_count": len(o.faces) if o.object_type == "Solid" else None,
                }
        data["ended_utc"] = stamp()
        h.save_project()
    except Exception as exc:
        data["error"] = repr(exc)
        data["ended_utc"] = stamp()
    finally:
        report.write_text(json.dumps(data, indent=2))
        if h is not None:
            h.release_desktop(close_projects=True, close_desktop=True)


if __name__ == "__main__":
    main()
