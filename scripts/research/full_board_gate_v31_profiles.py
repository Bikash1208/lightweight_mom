"""Export stored complex E fields on eight radial lines from solved v3.1 projects.

Read-only AEDT reopen; no geometry, setup, or solve call. Coordinates are in mm.
"""
import argparse
import csv
import json
import math
import os
import sys
import traceback
from pathlib import Path

import numpy as np

os.environ["ANSYSEM_ROOT261"] = r"C:\Program Files\ANSYS Inc\v261\AnsysEM"
from ansys.aedt.core import Hfss

PHIS = (0, 45, 90, 135, 180, 225, 270, 315)
RADII = tuple(range(151))
Z_MM = 2.035  # Top copper spans z=0..0.035 mm; sample 2 mm above its top.


def extract(state, project, out):
    out.mkdir(parents=True, exist_ok=True)
    temp = out / "aedt_temp"
    temp.mkdir(exist_ok=True)
    os.environ.update(TEMP=str(temp), TMP=str(temp))
    status = {"state": state, "project": str(project), "status": "OPENING", "resolves": 0}
    status_file = out / "status.json"
    status_file.write_text(json.dumps(status, indent=2))
    h = None
    try:
        h = Hfss(project=str(project), design=f"GATE_{state}", solution_type="Modal",
                 version="2026.1", non_graphical=True, new_desktop=True, close_on_exit=False)
        status["own_aedt_pid"] = h.desktop_class.aedt_process_id
        status["version"] = str(h.desktop_class.aedt_version_id)
        assert status["version"] == "2026.1"
        h.desktop_class.temp_directory = temp
        assert Path(str(h.desktop_class.temp_directory)).resolve() == temp.resolve()
        assert "Setup1 : LastAdaptive" in h.existing_analysis_sweeps
        assert str(h.setups[0].props["Frequency"]) == "5.68GHz"
        assert len([n for n in h.modeler.object_names if n.startswith("METAL_PEC")]) == 518
        points = [[r * math.cos(math.radians(phi)), r * math.sin(math.radians(phi)), Z_MM]
                  for phi in PHIS for r in RADII]
        field = out / f"{state}_E_complex.fld"
        status["status"] = "EXPORTING"
        status_file.write_text(json.dumps(status, indent=2))
        assert h.post.export_field_file(quantity="E", solution="Setup1 : LastAdaptive",
                                        output_file=str(field), sample_points=points,
                                        export_in_si_system=False,
                                        intrinsics={"Freq": "5.68GHz", "Phase": "0deg"}), "FIELD_EXPORT_FAILED"
        with field.open(encoding="utf-8") as stream:
            unit = stream.readline().strip()
            header = stream.readline().strip()
            raw = np.loadtxt(stream)
        assert unit == 'Unit: "mm"' and header == 'Complex Vector data "<Ex,Ey,Ez>"'
        assert raw.shape == (len(points), 9), raw.shape
        missing = np.flatnonzero(~np.isfinite(raw[:, 3:]).all(axis=1))
        # The r=0 sample lies inside the PEC feed pin; AEDT correctly returns NaN.
        assert set(missing.tolist()) == {k * len(RADII) for k in range(len(PHIS))}, missing
        assert np.isfinite(raw[:, :3]).all()
        status["nan_samples_inside_feed_pin"] = len(missing)
        assert np.allclose(raw[:, :3], points, rtol=0, atol=1e-9)
        # Three interleaved real/imaginary component pairs after x,y,z.
        e2 = np.sum(raw[:, 3:] ** 2, axis=1)
        ez2 = raw[:, 7] ** 2 + raw[:, 8] ** 2
        rows = []
        for n, (phi, radius) in enumerate((p, r) for p in PHIS for r in RADII):
            rows.append({"state": state, "phi_deg": phi, "r_mm": radius,
                         "x_mm": raw[n, 0], "y_mm": raw[n, 1], "z_mm": raw[n, 2],
                         "e_total_v_per_m": math.sqrt(float(e2[n])),
                         "e_z_v_per_m": math.sqrt(float(ez2[n])),
                         "e_total_sq": float(e2[n]), "e_z_sq": float(ez2[n])})
        with (out / f"{state}_radial.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        status.update(status="COMPLETE", rows=len(rows), field_path=str(field), csv_path=str(out / f"{state}_radial.csv"))
        status_file.write_text(json.dumps(status, indent=2))
        print("COMPLETE", state, len(rows), flush=True)
    except BaseException:
        status.update(status="FAILED", traceback=traceback.format_exc())
        status_file.write_text(json.dumps(status, indent=2))
        print(status["traceback"], file=sys.stderr, flush=True)
        raise
    finally:
        if h is not None:
            h.release_desktop(close_projects=True, close_desktop=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--state", required=True)  # S4 patch: any state name (was choices BOTH/A_ONLY/B_ONLY)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    extract(args.state, args.project.resolve(), args.out.resolve())


if __name__ == "__main__":
    main()
