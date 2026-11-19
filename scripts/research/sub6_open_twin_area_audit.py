"""Geometry-only readback of the D3p6 open twin in a fresh local AEDT session."""

import json
import math
from pathlib import Path
import traceback

from sub6_rotated_supercell_hfss import CELLS, CU, RESULTS, WORK, launch, need, write


def main():
    spec = CELLS["D3p6_OPEN"]
    project = WORK / "projects" / "OPEN_TWIN_AREA_AUDIT_resume5.aedt"
    need(not project.exists(), f"REFUSE_EXISTING_PROJECT {project}")
    h = None
    audit = {"status": "STARTING", "project": str(project)}
    path = RESULTS / "D3p6_OPEN_area_audit_resume5.json"
    write(path, audit)
    try:
        h = launch(project)
        audit["aedt_pid"] = h.desktop_class.aedt_process_id
        m = h.modeler
        m.model_units = "mm"
        r, w, depth = spec["r"], spec["w"], spec["depth"]
        disc = need(m.create_cylinder("Z", ["0mm", "0mm", "0mm"], f"{r}mm",
                                      f"{CU}mm", name="DISC", material="pec"), "CREATE_DISC")
        cuts = []

        def cut(name, x0, y0, dx, dy):
            cuts.append(need(m.create_box([f"{x0:.9f}mm", f"{y0:.9f}mm", f"{-CU:.9f}mm"],
                                          [f"{dx:.9f}mm", f"{dy:.9f}mm", f"{3*CU:.9f}mm"],
                                          name=name, material="vacuum"), f"CREATE_{name}"))

        for i in range(5):
            cut(f"LEG_{i}", -2.4 + 1.2*i - w/2, -depth/2-w/2, w, depth+w)
        for i in range(4):
            yc = depth/2 if i % 2 == 0 else -depth/2
            cut(f"BAR_{i}", -2.4 + 1.2*i - w/2, yc-w/2, 1.2+w, w)
        original_end = math.sqrt(r*r - (depth/2+w/2)**2)
        inner_rim = math.sqrt(r*r - (depth/2-w/2)**2)
        corrected_end = r + 0.1
        audit["end_run_extents_mm"] = {
            "original_formula_end_x": original_end,
            "actual_rim_at_inner_y": inner_rim,
            "uncut_bridge_at_inner_y": inner_rim - original_end,
            "corrected_end_x": corrected_end,
            "left": [-corrected_end, -2.2], "right": [2.2, corrected_end],
            "bottom_y": [-depth/2-w/2, -depth/2+w/2],
            "top_y": [depth/2-w/2, depth/2+w/2],
        }
        cut("END_LEFT", -corrected_end, -depth/2-w/2, corrected_end-2.2, w)
        cut("END_RIGHT", 2.2, depth/2-w/2, corrected_end-2.2, w)
        need(m.subtract(disc.name, [c.name for c in cuts], keep_originals=False),
             "SLOT_SUBTRACT")

        def objects():
            return [{"name": name, "material": m[name].material_name,
                     "volume_mm3": float(m[name].volume),
                     "bounding_box_mm": list(m[name].bounding_box)}
                    for name in m.object_names]

        audit["objects_after_subtract"] = objects()
        copper = [o for o in audit["objects_after_subtract"]
                  if o["material"].lower() == "pec" and o["volume_mm3"] > 0]
        audit["positive_copper_bodies_after_subtract"] = len(copper)
        if len(copper) == 1:
            result = m.separate_bodies(copper[0]["name"])
            audit["separate_result"] = [x.name for x in result] if result else []
        else:
            audit["separate_result"] = "already_separate"
        audit["objects_after_separation"] = objects()
        copper = [o for o in audit["objects_after_separation"]
                  if o["material"].lower() == "pec" and o["volume_mm3"] > 0]
        audit["copper_bodies"] = copper
        audit["areas_mm2"] = [o["volume_mm3"] / CU for o in copper]
        audit["total_area_mm2"] = sum(audit["areas_mm2"])
        need(len(copper) == 2, f"COPPER_BODIES_{len(copper)}")
        need(abs(audit["total_area_mm2"] - spec["area"]) / spec["area"] <= .02,
             f"COPPER_AREA_{audit['total_area_mm2']:.6f}")
        need(all(abs(a - spec["area"]/2) / (spec["area"]/2) <= .02
                 for a in audit["areas_mm2"]), f"HALF_AREA_{audit['areas_mm2']}")
        need(h.save_project(), "SAVE_AUDIT_PROJECT")
        audit["status"] = "PASS"
    except BaseException as exc:
        audit.update(status="FAILED", error=f"{type(exc).__name__}: {exc}",
                     traceback=traceback.format_exc())
        raise
    finally:
        write(path, audit)
        if h is not None:
            h.release_desktop(close_projects=True, close_desktop=True)


if __name__ == "__main__":
    main()
