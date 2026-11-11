"""Reopen one saved anchor build in a private AEDT session and verify it."""
import json
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

from ansys.aedt.core import Hfss


def utc():
    return datetime.now(timezone.utc).isoformat()


def main():
    project = Path(sys.argv[1]).resolve()
    build = json.loads(Path(sys.argv[2]).read_text())
    report = Path(sys.argv[3]).resolve()
    result = {"project": str(project), "start_utc": utc(), "status": "VALIDATING"}
    report.write_text(json.dumps(result, indent=2))
    h = None
    try:
        assert build["status"] == "BUILT_UNSOLVED" and build["geometry_gate"] == "PASS"
        h = Hfss(project=str(project), design="AnchorPreflight", solution_type="Modal",
                 version="2025.2", non_graphical=True, new_desktop=True, close_on_exit=False)
        result["own_aedt_pid"] = h.desktop_class.aedt_process_id
        result["aedt_version"] = str(h.desktop_class.aedt_version_id)
        assert result["aedt_version"] == "2025.2"
        m = h.modeler
        names = list(m.object_names)
        metal = [n for n in names if n.startswith("METAL_PEC")]
        diode = [n for n in names if n.startswith("D") and n not in ("SUBSTRATE", "SUBSTRATE_FR4")]
        assert (len(names), len(metal), len(diode)) == (1244, 1000, 240)
        for name in ("SUBSTRATE", "SUBSTRATE_FR4", "AIR"):
            assert m[name].material_name.lower() == build["material_readback"][name].lower()
        bbox = [round(float(v), 6) for v in m["SUBSTRATE_FR4"].bounding_box]
        assert bbox == build["bias_after"], (bbox, build["bias_after"])
        b = next(x for x in h.boundaries if x.name == "BAR90_02LRH")
        assert len(b.props["Objects"]) == 240
        for key in ("Resistance", "Reactance"):
            assert str(b.props[key]) == str(build["impedance_readback"][key])
        result["setup_names"] = [s.name for s in h.setups]
        result["setup_settings"] = {s.name: {k: str(s.props.get(k)) for k in
                                    ("Frequency", "MaximumPasses", "MaxDeltaS")} for s in h.setups}
        report.write_text(json.dumps(result, indent=2, default=str))
        assert "Setup1" in result["setup_names"], result["setup_names"]
        target_setup = next(s for s in h.setups if s.name == "Setup1")
        for key in ("Frequency", "MaximumPasses", "MaxDeltaS"):
            assert str(target_setup.props[key]) == str(build["setup_readback"][key])
        messages, ok = h.validate_full_design(output_dir=str(report.parent), ports=1)
        result["validation_messages"] = messages
        assert ok, "AEDT_DESIGN_VALIDATION"
        h.save_project()
        result["status"] = "PASS"
        result["readback"] = {"total_objects": len(names), "metal": len(metal),
                              "diodes": len(diode), "bias_bbox_mm": bbox,
                              "resistance": str(b.props["Resistance"]),
                              "reactance": str(b.props["Reactance"])}
    except BaseException:
        result["status"] = "FAILED"
        result["traceback"] = traceback.format_exc()
        raise
    finally:
        result["end_utc"] = utc()
        report.write_text(json.dumps(result, indent=2, default=str))
        if h is not None:
            h.release_desktop(close_projects=True, close_desktop=True)


if __name__ == "__main__":
    main()
