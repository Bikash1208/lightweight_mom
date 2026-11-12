"""Extract an already solved Run 3 row without repeating its HFSS solve."""
import json
import os
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

import psutil

from sub6_planA_run3 import OUT, ROOT, TEMP, WORK, extract, write


def utc():
    return datetime.now(timezone.utc).isoformat()


def main():
    run = sys.argv[1]
    project = WORK / run / f"{run}.aedt"
    result = OUT / f"{run}.json"
    state = json.loads((WORK / "controller_status.json").read_text())
    row = dict(state["rows"][run])
    assert row["status"] == "SOLVED_EXTRACTING", "ROW_WAS_NOT_SOLVED"
    assert project.is_file() and (WORK / run / f"{run}.aedtresults").is_dir()
    temp = TEMP / run
    temp.mkdir(parents=True, exist_ok=True)
    os.environ.update(TEMP=str(temp), TMP=str(temp), AEDT_PRIVATE_TEMP=str(temp))
    row["status"] = "RECOVERING_EXTRACTION"
    row["recovery_started_utc"] = utc()
    write(result, row)
    h = None
    try:
        from ansys.aedt.core import Hfss
        h = Hfss(project=str(project), design="AnchorPreflight", solution_type="Modal",
                 version="2025.2", non_graphical=True, new_desktop=True,
                 close_on_exit=False)
        row["recovery_own_aedt_pid"] = h.desktop_class.aedt_process_id
        assert str(h.desktop_class.aedt_version_id) == "2025.2"
        h.desktop_class.temp_directory = temp
        assert Path(h.desktop_class.temp_directory).resolve() == temp.resolve()
        assert [s.name for s in h.setups] == ["Setup1"]
        assert "Setup1 : LastAdaptive" in h.existing_analysis_sweeps
        extract(h, run, OUT / run, row)
        try:
            monitor = json.loads((WORK / "controller_monitor.json").read_text())
            row["independent_monitor_peak_ram_gb"] = monitor["peak_ram_gb_by_run"][run]
        except (OSError, KeyError, ValueError):
            row["independent_monitor_peak_ram_gb"] = None
        row["peak_ram_gb"] = row["aedt_profile_peak_ram_gb"]
        row["status"] = "COMPLETE"
        row["ended_utc"] = utc()
    except BaseException:
        row["status"] = "EXTRACTION_FAILED"
        row["traceback"] = traceback.format_exc()
        with (WORK / f"{run}_recovery.log").open("a", encoding="utf-8") as f:
            f.write(row["traceback"])
        raise
    finally:
        write(result, row)
        if h is not None:
            h.release_desktop(close_projects=True, close_desktop=True)
        pid = row.get("recovery_own_aedt_pid")
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
