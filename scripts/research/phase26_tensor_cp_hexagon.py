"""Phase 26 C: pure-RHCP rotating-slot tensor hologram on ONE 150 mm hexagon (RT/duroid 5880 1.575 mm).

    python scripts/research/phase26_tensor_cp_hexagon.py <out.json> [--quick]

Design from the SOLVER_SUPERCELL table (`supercell_table_*.json`, nearest grid, mesh-stable slots only: bridges
>= 0.43 mm), mesh pitch 6.4 mm, D 6.0 mm circle, hexagon circumradius 150 mm (vertex at +/-y, the orientation of
port 3 of `four_port_hex150_gap0`), feed at the centre, beam theta 30 deg / phi 225 deg, 10.54 GHz.
A scan over the design conversion `eta_design` and amplitude taper picks the design of the highest predicted
gain whose RHCP SLL is acceptable; every scanned design is reported. Everything is MODEL_ONLY.
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.hologram.tensor_cp_design import *  # noqa: E402,F401,F403  (the design flow lives in the module; this script is a wrapper)
from src.hologram.tensor_cp_design import design_single_hexagon  # noqa: E402


def main(out_path, quick=False):
    res, cells = design_single_hexagon(quick=quick)
    Path(out_path).with_name(Path(out_path).stem + "_cells.json").write_text(json.dumps(cells))
    Path(out_path).write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main(sys.argv[1], "--quick" in sys.argv)
