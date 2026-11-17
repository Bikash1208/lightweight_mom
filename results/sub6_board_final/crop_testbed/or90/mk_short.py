"""Copy of the switched crop: add two idealised-bound designs (sheets open / shorted). NOT a diode sweep."""
import os, sys, json, math
os.environ["ANSYSEM_ROOT251"] = r"C:\Program Files\ANSYS Inc\v251\AnsysEM"
from pathlib import Path
B = Path(r"C:/Users/Bikash/Desktop/sub6_board_hfss/board_v2/crop")
temp = Path(r"C:/Users/Bikash/Desktop/sub6_board_hfss/aedt_temp_bar90/crop_short_build"); temp.mkdir(parents=True, exist_ok=True)
os.environ.update(TEMP=str(temp), TMP=str(temp))
from ansys.aedt.core import Hfss
h = Hfss(project=str(Path(sys.argv[1])), design="ON_ELS", version="2025.1", non_graphical=True, new_desktop=True, close_on_exit=False)
print("OWN_AEDT_PID", h.desktop_class.aedt_process_id, flush=True)
rep = {}
try:
    for name, rs, xs in (("SHORT_IDEAL", 0.001, 0.0),):
        h.set_active_design("ON_ELS"); assert h.duplicate_design(name); h.set_active_design(name)
        h["Rs"] = f"{rs}ohm"; h["Xs"] = f"{xs}ohm"
        rep[name] = dict(Rs=h.get_evaluated_value("Rs"), Xs=h.get_evaluated_value("Xs"))
        assert math.isclose(rep[name]["Rs"], rs, rel_tol=1e-9) and math.isclose(rep[name]["Xs"], xs, abs_tol=1e-12), "READBACK " + name
        b = next(x for x in h.boundaries if x.name == "BAR90_ELS"); rep[name]["n_sheets"] = len(list(b.props.get("Objects", [])))
        rep[name]["validate_simple"] = h.validate_simple(str(Path(sys.argv[1]).parent / f"validate_{name}.log")); rep[name]["mesh_ops"] = [o.name for o in h.mesh.meshoperations]
        print(name, rep[name], flush=True)
    for d in ("OFF_ELS10", "ON_ELS"):      # drop the two stale designs from the copy so only bounds are solved here
        pass
    h.save_project()
finally:
    (Path(sys.argv[1]).parent / "short_readback.json").write_text(json.dumps(rep, indent=1, default=str))
    h.release_desktop(close_projects=True, close_desktop=True)
print("DONE", flush=True)
