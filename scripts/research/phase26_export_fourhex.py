"""Phase 26: export the 4-hexagon tensor-CP board A (static beams phi 45/135/225/315, theta 30, all RHCP) for HFSS.

    python scripts/research/phase26_export_fourhex.py <design_dir> [out_dir]

`design_dir` holds `fourhex_design.json` + `fourhex_A_cells.json` written by `phase26_fourhex_design.py`.
Writes `hfss_build_board.py` (AEDT, native exact geometry, 4 monopoles + 4 lumped ports, DRIVE_PORT option),
`layout.csv/.dxf/.png/.json`, `predictions.json`, `README.md`. No STEP (the AEDT script is the build path).
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.export.tensor_cp_export import README, aedt_script_fourhex, write_layout_fourhex  # noqa: E402,F401
from src.hologram.tensor_cp_design import PHI, board_a_predictions, fourhex_a_predictions  # noqa: E402

OUT = ROOT / "hfss_exports/hex4_tensor_cp_static"


def main(design_dir, out=OUT):
    design_dir, out = Path(design_dir), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    design = json.loads((design_dir / "fourhex_design.json").read_text())
    cells = json.loads((design_dir / "fourhex_A_cells.json").read_text())["cells"]
    beams = "port 1 (theta 30, phi 45), port 2 (30, 135), port 3 (30, 225), port 4 (30, 315): all RHCP"
    (out / "hfss_build_board.py").write_text(aedt_script_fourhex(cells, beams=beams))
    write_layout_fourhex(cells, out, f"hex4_tensor_cp_static: {len(cells)} slotted circles (D 6.0 mm), shared 6.4 mm lattice")
    scope = json.loads((ROOT / "results/phase26/fourhex_scope.json").read_text())
    cond = board_a_predictions(design)
    pred = fourhex_a_predictions(design, scope, len(cells))
    (out / "predictions.json").write_text(json.dumps(pred, indent=1, default=float))
    lines = ["| port | beam (theta, phi) | model peak | gain dBi (s=1 / s=1.7) | XPD dB | AR peak / worst-in-beam dB | neighbour LHCP lobe (rel. main peak, direction) |", "|---|---|---|---|---|---|---|"]
    for j in ("1", "2", "3", "4"):
        a, b = cond["s=1.0"][j], cond["s=1.7"][j]
        lines.append(f"| {j} | (30, {PHI[int(j)]:g}) | {a['model_peak_theta_phi_deg']} | {a['realized_gain_rhcp_dbi_anchored']} / {b['realized_gain_rhcp_dbi_anchored']} | "
                     f"{a['xpd_at_peak_db']} / {b['xpd_at_peak_db']} | {a['ar_at_peak_db']} / {a['ar_worst_in_minus3db_beam_db']} (s=1) | "
                     f"{a['check_item_neighbour_lhcp_lobe']['rel_to_main_rhcp_peak_db']} / {b['check_item_neighbour_lhcp_lobe']['rel_to_main_rhcp_peak_db']} dB at "
                     f"{a['check_item_neighbour_lhcp_lobe']['theta_phi_deg']} ({a['check_item_neighbour_lhcp_lobe']['deg_from_beam']} deg off) |")
    (out / "README.md").write_text(README.format(n_cells=len(cells), table="\n".join(lines)))
    return pred


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], a[1] if len(a) > 1 else OUT)
