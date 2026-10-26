"""Phase 26: compare the user's HFSS run of board A (`hfss_exports/hex4_tensor_cp_static`) with `predictions.json`.

    python scripts/research/phase26_compare_fourhex.py --s4p board.s4p
        --gain 3=gain_port3.csv --gain 2=gain_port2.csv [--params 3=params3.txt --params 2=params2.txt]
        [--pred predictions.json] [--out report.json]

The logic lives in `src/hologram/tensor_cp_compare.py` (Phase 27 item 5); this script is a wrapper.
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.hologram import tensor_cp_compare as _m  # noqa: E402

globals().update({k: v for k, v in vars(_m).items() if not k.startswith("__")})
single = _m  # the old four-hex script reached the single-hex helpers as `single.X`
main = _m.main_fourhex
compare = _m.compare_fourhex
markdown = _m.markdown_fourhex

if __name__ == "__main__":
    _m.main_fourhex()
