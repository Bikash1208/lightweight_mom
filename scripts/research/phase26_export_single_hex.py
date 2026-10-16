"""Phase 26 D: export the single 150 mm hexagon tensor-CP board (rotating slotted circles) for HFSS.

    python scripts/research/phase26_export_single_hex.py [cells.json results.json out_dir]

Writes `layout.json/.csv/.dxf/.png`, `board.step` (CadQuery, if available), `hfss_build_board.py` (AEDT,
native exact geometry, rotated corners, no STEP import needed) and `predictions.json`. The generators live in
`src/export/tensor_cp_export.py` (Phase 27); this script reads the Phase 26 C results file and calls them.
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.export.tensor_cp_export import (  # noqa: E402,F401
    BOARD_X,
    BOARD_Y,
    D_MM,
    H_MM,
    R_HEX,
    aedt_script,
    hexagon_vertices,
    slot_corners,
    write_single_board,
)
from src.hologram.tensor_cp_design import single_hex_predictions  # noqa: E402

CELLS = ROOT / "results/phase26/tensor_cp_hexagon_cells.json"
RESULTS = ROOT / "results/phase26/tensor_cp_hexagon.json"
OUT = ROOT / "hfss_exports/hex150_tensor_cp_single"


def write_all(cells_path=CELLS, results_path=RESULTS, out=OUT, name="hex150_tensor_cp_single", note="", predictions=None, cells_override=None):
    res = json.loads(Path(results_path).read_text())
    if cells_override is not None:
        cells = cells_override
    else:
        c = json.loads(Path(cells_path).read_text())
        cells = [
            {"x_mm": c["x_mm"][i], "y_mm": c["y_mm"][i], "slot_L_mm": c["slot_L_mm"][i], "slot_W_mm": c["slot_W_mm"][i],
             "slot_angle_deg": c["slot_angle_deg"][i], "alpha_np_per_m": c["alpha_np_per_m"][i]} for i in range(len(c["x_mm"]))
        ]
    pred = predictions if predictions is not None else single_hex_predictions(res)
    write_single_board(cells, out, pred, name=name, note=note)
    return (cells, None, res) if predictions is not None else (cells, res["selected"], res)


if __name__ == "__main__":
    a = sys.argv[1:]
    write_all(*(a[:3] if len(a) >= 3 else ()))
