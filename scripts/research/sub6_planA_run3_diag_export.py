"""Read-only far-field export diagnostic for an already solved Run 3 project."""
import json
import os
import sys
import traceback
from pathlib import Path

import psutil

from sub6_planA_run3 import ROOT, TEMP, WORK, write
from sub6_planA_restart_extract import ffd_data, point


def main():
    run = sys.argv[1]
    project = WORK / run / f"{run}.aedt"
    out = WORK / f"{run}_export_diag"
    out.mkdir(exist_ok=True)
    temp = TEMP / run
    os.environ.update(TEMP=str(temp), TMP=str(temp), AEDT_PRIVATE_TEMP=str(temp))
    report = out / "report.json"
    d = {"run_id": run, "status": "STARTING"}
    write(report, d)
    h = None
    try:
        from ansys.aedt.core import Hfss
        h = Hfss(project=str(project), design="AnchorPreflight", solution_type="Modal",
                 version="2025.2", non_graphical=True, new_desktop=True, close_on_exit=False)
        d["own_aedt_pid"] = h.desktop_class.aedt_process_id
        h.desktop_class.temp_directory = temp
        assert Path(h.desktop_class.temp_directory).resolve() == temp.resolve()
        d["all_sources"] = h.get_all_sources()
        d["all_source_modes"] = h.get_all_source_modes()
        d["categories"] = h.post.available_quantities_categories()
        d["s11_db"] = point(h, "dB(S(Port1,Port1))")
        for th, ph in ((25, 0), (25, 180), (24, 180)):
            d[f"th{th}_ph{ph}"] = {
                "rhcp_db": point(h, "dB(GainRHCP)", th, ph),
                "lhcp_db": point(h, "dB(GainLHCP)", th, ph)}
        direct = out / "direct.ffd"
        d["direct_export_ok"] = bool(h.export_element_pattern(
            frequencies=[5.8e9], setup="Setup1 : LastAdaptive", sphere="Global",
            variations={}, element_name="direct", output_dir=str(out)))
        if direct.is_file():
            _, _, et, ep = ffd_data(direct)
            d["direct_ffd_max_abs"] = max(abs(et).max().item(), abs(ep).max().item())
            d["direct_ffd_nonzero"] = int((et != 0).sum() + (ep != 0).sum())
        d["status"] = "COMPLETE"
    except BaseException:
        d["status"] = "FAILED"
        d["traceback"] = traceback.format_exc()
        raise
    finally:
        write(report, d)
        if h is not None:
            h.release_desktop(close_projects=True, close_desktop=True)
        pid = d.get("own_aedt_pid")
        if pid:
            try:
                proc = psutil.Process(pid)
                for child in proc.children(recursive=True):
                    child.terminate()
                proc.terminate()
            except psutil.Error:
                pass


if __name__ == "__main__":
    main()
