"""Export sampled eigenmode field magnitudes from one solved, owned project."""

import argparse
import json
import math
import os
from pathlib import Path
import traceback

from ansys.aedt.core import Hfss

from sub6_rotated_supercell_hfss import AEDT_ROOT, RESULTS, VERSION, WORK, write


def main(cell, phase, tag, selected_modes):
    label = f"SUPER_{cell}_ph{phase}_PML40{('_'+tag) if tag else ''}"
    data = json.loads((RESULTS / f"{label}_modes.json").read_text(encoding="utf-8"))
    project = Path(data["project"])
    temp = WORK / "temp" / f"{label}_fields"
    temp.mkdir(parents=True, exist_ok=True)
    os.environ.update(TEMP=str(temp), TMP=str(temp), ANSYSEM_ROOT252=str(AEDT_ROOT))
    h = None
    log = {"cell": cell, "phase_deg": phase, "project": str(project), "exports": []}
    try:
        h = Hfss(project=str(project), design="EigenCell", solution_type="Eigenmode",
                 version=VERSION, non_graphical=True, new_desktop=True, close_on_exit=False)
        h.desktop_class.temp_directory = temp
        log["aedt_pid"] = h.desktop_class.aedt_process_id
        points = []
        for n in range(8):
            psi = math.radians(n * 22.5)
            cx = (n - 3.5) * 11.6
            for offset in (-3.5, 3.5):
                points.append([cx + offset * math.cos(psi),
                               offset * math.sin(psi), 1.0])
        log["sample_points_mm"] = points
        for mode in data["modes"]:
            if not 5.3 <= mode["f_ghz"] <= 6.5 or mode["Q"] >= 100000:
                continue
            if selected_modes and mode["mode"] not in selected_modes:
                continue
            path = RESULTS / f"{label}_mode{mode['mode']}_MagE.fld"
            record = {"mode": mode["mode"], "f_ghz": mode["f_ghz"], "Q": mode["Q"],
                      "path": str(path)}
            try:
                sources = {f"mode_{i}": ("1" if i == mode["mode"] else "0")
                           for i in range(1, len(data["modes"]) + 1)}
                record["edit_sources_ok"] = bool(h.edit_sources(sources,
                                            eigenmode_stored_energy=True))
                record["export_ok"] = bool(h.post.export_field_file(
                    quantity="Mag_E", solution="EigenCell : LastAdaptive",
                    sample_points=points, export_in_si_system=False,
                    output_file=str(path)))
                record["file_exists"] = path.exists()
            except Exception as exc:
                record["error"] = f"{type(exc).__name__}: {exc}"
            log["exports"].append(record)
            write(RESULTS / f"{label}_field_exports.json", log)
    except BaseException:
        log["fatal_traceback"] = traceback.format_exc()
        write(RESULTS / f"{label}_field_exports.json", log)
        raise
    finally:
        if h is not None:
            own_pid = int(h.desktop_class.aedt_process_id)
            h.release_desktop(close_projects=True, close_desktop=True)
            import psutil
            if psutil.pid_exists(own_pid):
                proc = psutil.Process(own_pid)
                if proc.exe().lower() == str(AEDT_ROOT / "ansysedt.exe").lower():
                    proc.terminate()
                    try:
                        proc.wait(timeout=10)
                    except psutil.TimeoutExpired:
                        pass
                    log["terminated_own_aedt_pid_after_release"] = own_pid
                    write(RESULTS / f"{label}_field_exports.json", log)
            lock = project.with_suffix(".aedt.lock")
            if lock.exists() and not psutil.pid_exists(own_pid):
                lock.unlink()
                log["removed_own_stale_lock"] = str(lock)
                write(RESULTS / f"{label}_field_exports.json", log)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", required=True)
    ap.add_argument("--phase", type=int, required=True)
    ap.add_argument("--tag", default="")
    ap.add_argument("--modes", nargs="*", type=int, default=[])
    args = ap.parse_args()
    main(args.cell, args.phase, args.tag, set(args.modes))
