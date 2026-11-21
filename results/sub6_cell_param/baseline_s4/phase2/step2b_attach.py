import os, json, traceback
from pathlib import Path
os.environ["ANSYSEM_ROOT261"] = r"C:\Program Files\ANSYS Inc\v261\AnsysEM"
W = Path(r"C:\work\s4_work\both_step2"); res = {}; out = W/"step2b.json"
def save(): out.write_text(json.dumps(res, indent=1, default=str))
from ansys.aedt.core import Hfss
h = None
try:
    h = Hfss(project="BOARD_R150_V31_BOTH", design="GATE_BOTH", solution_type="Modal", version="2026.1", non_graphical=True, new_desktop=False, aedt_process_id=20564, close_on_exit=False)
    res["attached_pid"] = h.desktop_class.aedt_process_id
    ms = h.export_mesh_stats("Setup1", output_file=str(W/"mesh_stats.ms")); res["ms"] = str(ms)
    res["ms_text"] = Path(ms).read_text(errors="replace")[:3000] if ms and Path(ms).exists() else None; save()
    msgs = {}
    for sev in (0, 1, 2, 3):
        try: msgs[sev] = [str(x) for x in h.odesktop.GetMessages("BOARD_R150_V31_BOTH", "GATE_BOTH", sev)]
        except Exception as e: msgs[sev] = "ERR " + str(e)[:100]
    res["msgs"] = msgs; save()
    h.save_project(); res["saved"] = True; res["results_exists"] = (W/"BOARD_R150_V31_BOTH.aedtresults").exists(); save()
except BaseException:
    res["error"] = traceback.format_exc(); save()
finally:
    try:
        if h is not None: h.release_desktop(close_projects=True, close_desktop=True)
    except Exception: pass
