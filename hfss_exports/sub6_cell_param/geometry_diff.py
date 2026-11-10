"""Exact geometry diff: PyAEDT-built nominal cell vs the user's BASELINE STEP, by boolean differences.

Runs on a SCRATCH COPY of a built project (never the original): imports the STEP next to the built objects and, per
object pair, computes volume(A-B) and volume(B-A) on clones plus their bounding boxes. Geometry only, no solve.

python geometry_diff.py --project <scratch copy of a built .aedt> --step BASELINE_vias_fr4_stubs.step --output geometry_diff.json
"""
import argparse
import json
import re
from pathlib import Path
import sys

from ansys.aedt.core import Hfss

PAIRS = ["METAL_PEC", "SUBSTRATE", "SUBSTRATE_FR4", "AIR"]


def vol(obj):
    try:
        return float(obj.volume)
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", type=Path, required=True)
    ap.add_argument("--step", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    h = Hfss(project=str(a.project.resolve()), design="CellParam", solution_type="Eigenmode", version="2025.1",
             non_graphical=True, new_desktop=True, close_on_exit=False)
    out = {"project": str(a.project), "step": str(a.step), "objects_before": list(h.modeler.object_names)}
    try:
        m = h.modeler
        before = set(m.object_names)
        ok = m.import_3d_cad(str(a.step.resolve()), create_group=False)
        new = [n for n in m.object_names if n not in before]
        out["import_ok"] = bool(ok); out["imported"] = new
        out["bbox_step_objects"] = {n: [round(x, 6) for x in m[n].bounding_box] for n in new}
        out["bbox_mine"] = {n: [round(x, 6) for x in m[n].bounding_box] for n in before if n in m.object_names}
        out["volume_step"] = {n: vol(m[n]) for n in new}
        out["volume_mine"] = {n: vol(m[n]) for n in before if n in m.object_names}
        # imported name for product X is X or X<suffix>; match by prefix
        res = {}
        for name in PAIRS:
            cand = [n for n in new if re.fullmatch(re.escape(name) + r"(_\d+)?", n)]
            if len(cand) > 1:
                # the STEP stores disconnected metal components as separate solids; mine is one united object
                united = m.unite(cand)
                cand = [united if isinstance(united, str) else getattr(united, "name", cand[0])]
                res.setdefault("_united", {})[name] = cand[0]
            if name not in before or len(cand) != 1:
                res[name] = {"status": "UNMATCHED", "candidates": cand}; continue
            B = cand[0]; r = {"imported_name": B, "volume_mine": vol(m[name]), "volume_step": vol(m[B])}
            for tag, blank, tool in (("mine_minus_step", name, B), ("step_minus_mine", B, name)):
                c1 = m[blank].clone(); c2 = m[tool].clone()
                try:
                    done = m.subtract(c1.name, c2.name, keep_originals=False)
                    if done and c1.name in m.object_names:
                        r[tag] = {"volume": vol(m[c1.name]), "bbox": [round(x, 6) for x in m[c1.name].bounding_box]}
                    else:
                        r[tag] = {"volume": 0.0, "note": "empty result / subtract returned False"}
                except Exception as exc:
                    r[tag] = {"error": str(exc)}
                for n in (c1.name, c2.name):
                    try:
                        if n in m.object_names:
                            m[n].delete()
                    except Exception:
                        pass
            res[name] = r
        out["diff"] = res
        sheet = [n for n in new if n.startswith("DIODE_SHEET")]
        out["diode_sheet_step_bbox"] = {n: out["bbox_step_objects"][n] for n in sheet}
        out["diode_sheet_mine_bbox"] = out["bbox_mine"].get("DIODE_SHEET")
        out["status"] = "DONE"
    except Exception as exc:
        out["status"] = "ERROR"; out["error"] = repr(exc)
    finally:
        a.output.write_text(json.dumps(out, indent=2, default=str) + "\n", encoding="utf-8")
        try:
            h.release_desktop()
        except Exception:
            pass
    print(out.get("status"), out.get("error", ""))


if __name__ == "__main__":
    main()
