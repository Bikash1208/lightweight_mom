"""Write the pre-build Option 1 geometry and local copper DC report."""
import json
from pathlib import Path

from opt1_geometry import all_geometries, DISC_R, PERIOD

here = Path(__file__).resolve().parent
out = here.parents[1] / "results" / "sub6_cell_param" / "lab_B_opt1"
geometries = all_geometries(wide_channel=0.5)
report = {
    "status": "PASS" if len(geometries) == 10 and all(g["status"] == "PASS" for g in geometries) else "CONFLICT",
    "geometry_count": len(geometries),
    "rule": "0.10 line / 0.40 channel and 0.20 line / 0.50 channel; both leave 0.15 mm per side",
    "local_dc_rule": "island and local disc copper remain disjoint; only the diode sheet bridges their 0.30 mm gap",
    "periodic_receiver": "opposite-edge top-copper tab connects the outbound line across the periodic seam to the next disc rim",
    "disc_rim_to_cell_boundary_mm": PERIOD/2-DISC_R,
    "geometries": geometries,
}
out.mkdir(parents=True, exist_ok=True)
(out / "geometry_precheck.json").write_text(json.dumps(report, indent=2) + "\n")
print(f"{report['status']}: {len(geometries)} geometries; local DC contacts "
      f"{sum(g['checks']['island_local_disc_copper_dc_contact'] for g in geometries)}; "
      f"clearance 0.15 mm for both line widths; receiver bridge {report['disc_rim_to_cell_boundary_mm']:.3f} mm")
if report["status"] != "PASS":
    raise SystemExit(1)
