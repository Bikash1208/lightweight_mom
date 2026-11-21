"""Diagnostic re-run of the failed SMA test-2 solve on a COPY of the saved project (no change): capture AEDT messages + profile. Saves/exports before release."""
import os, json, shutil, time, traceback, re
from pathlib import Path
os.environ["ANSYSEM_ROOT261"] = r"C:\Program Files\ANSYS Inc\v261\AnsysEM"; os.environ["ANSYSLMD_LICENSE_FILE"] = "1055@10.100.105.13"
W = Path(r"C:\work\s4_work\sma_t2_diag"); W.mkdir(exist_ok=True, parents=True)
shutil.copy(r"C:\work\s4_work\sma_t2\run\BOARD_R150_V31_BOTH_SMA.aedt", W/"BOARD_R150_V31_BOTH_SMA.aedt")
temp = W/"aedt_temp"; temp.mkdir(exist_ok=True); os.environ.update(TEMP=str(temp), TMP=str(temp))
from ansys.aedt.core import Hfss
res = {}; out = W/"result.json"
def save(): out.write_text(json.dumps(res, indent=1, default=str))
h = None
try:
    h = Hfss(project=str(W/"BOARD_R150_V31_BOTH_SMA.aedt"), design="SMA_BOTH", solution_type="Modal", version="2026.1", non_graphical=True, new_desktop=True, close_on_exit=False)
    res["pid"] = h.desktop_class.aedt_process_id; st = h.setups[0]
    res["setup"] = {k: str(st.props.get(k)) for k in ("Frequency", "MinimumPasses", "MaximumPasses", "DrivenSolverType")}; res["setups"] = [s.name for s in h.setups]; save()
    t0 = time.monotonic(); ok = h.analyze_setup("Setup1", cores=16, use_auto_settings=False, blocking=True); res["analyze_ok"] = bool(ok); res["wall_s"] = round(time.monotonic()-t0, 1); save()
    msgs = {}
    for sev in (0, 1, 2, 3):
        try: msgs[sev] = [str(x) for x in h.odesktop.GetMessages("BOARD_R150_V31_BOTH_SMA", "SMA_BOTH", sev)]
        except Exception as e: msgs[sev] = "ERR " + str(e)[:100]
    res["messages"] = msgs; save()
    try:
        prof = Path(h.export_profile("Setup1", output_file=str(W/"solve.prof"))); res["profile"] = str(prof)
    except Exception as e: res["profile_err"] = str(e)[:200]
    save(); h.save_project(); res["saved"] = True; save()
except BaseException:
    res["error"] = traceback.format_exc(); save()
finally:
    try:
        if h is not None: h.save_project(); h.release_desktop(close_projects=True, close_desktop=True)
    except Exception: pass
