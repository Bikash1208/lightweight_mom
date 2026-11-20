"""Read the four requested gain samples from saved v3.1 HFSS projects, sequentially."""
import json
import os
import sys
import traceback
from pathlib import Path

os.environ["ANSYSEM_ROOT261"] = r"C:\Program Files\ANSYS Inc\v261\AnsysEM"
from ansys.aedt.core import Hfss

from full_board_gate_v31_extract import solution_point


def main():
    out = Path(sys.argv[1]).resolve()
    out.mkdir(parents=True, exist_ok=True)
    jobs = [
        ("BOTH", Path(sys.argv[2]).resolve()),
        ("A_ONLY_IDEAL", Path(sys.argv[3]).resolve()),
        ("B_ALL_OFF_IDEAL", Path(sys.argv[4]).resolve()),
        ("ALL_PLAIN", Path(sys.argv[5]).resolve()),
    ]
    report = {"status": "RUNNING", "samples": {}, "solves_started": 0}
    path = out / "gain_at_reference_directions.json"
    path.write_text(json.dumps(report, indent=2))
    for state, project in jobs:
        temp = out / f"temp_{state}"
        temp.mkdir(exist_ok=True)
        os.environ.update(TEMP=str(temp), TMP=str(temp))
        h = None
        try:
            h = Hfss(project=str(project), design=f"GATE_{state}", solution_type="Modal",
                     version="2026.1", non_graphical=True, new_desktop=True, close_on_exit=False)
            h.desktop_class.temp_directory = temp
            assert str(h.desktop_class.aedt_version_id) == "2026.1"
            assert "Setup1 : LastAdaptive" in h.existing_analysis_sweeps
            entry = {"own_aedt_pid": h.desktop_class.aedt_process_id, "project": str(project)}
            for theta, phi in ((32, 0), (36, 180)):
                entry[f"th{theta}_ph{phi}"] = {
                    "RHCP_gain_dbi": solution_point(h, "dB(GainRHCP)", theta, phi),
                    "LHCP_gain_dbi": solution_point(h, "dB(GainLHCP)", theta, phi),
                }
            report["samples"][state] = entry
            path.write_text(json.dumps(report, indent=2))
        except BaseException:
            report["status"] = "FAILED"
            report["traceback"] = traceback.format_exc()
            path.write_text(json.dumps(report, indent=2))
            raise
        finally:
            if h is not None:
                h.release_desktop(close_projects=True, close_desktop=True)
    report["status"] = "COMPLETE"
    path.write_text(json.dumps(report, indent=2))
    print("COMPLETE", path, flush=True)


if __name__ == "__main__":
    main()
