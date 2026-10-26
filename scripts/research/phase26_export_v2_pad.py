"""Phase 26 V2: export the single-hex tensor-CP board with the graded launch pad (see `phase26_v2_pad_design.py`).

    python scripts/research/phase26_export_v2_pad.py
"""

import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
_p = Path(__file__).resolve().parent / "phase26_export_single_hex.py"
_s = importlib.util.spec_from_file_location("p26_export_single", _p)
ex = importlib.util.module_from_spec(_s)
sys.modules["p26_export_single"] = ex
_s.loader.exec_module(ex)

OUT = ROOT / "hfss_exports/hex150_tensor_cp_single_v2_pad"


def main():
    cells = json.loads((ROOT / "results/phase26/v2_pad_design_cells.json").read_text())
    design = json.loads((ROOT / "results/phase26/v2_pad_design.json").read_text())
    note = "# V2: the %d cells with r < 30 mm are cross-slotted circles (cross = 1) with the width graded 2.0 mm (r <= 15 mm) to ~1.1 mm (r = 30 mm); the feed cell is empty." % design["pad"]["n_pad_cells"]
    ex.write_all(results_path=ROOT / "results/phase26/tensor_cp_hexagon.json", out=OUT, name="hex150_tensor_cp_single_v2_pad",
                 note=note, predictions=design, cells_override=cells)
    return design


if __name__ == "__main__":
    main()
