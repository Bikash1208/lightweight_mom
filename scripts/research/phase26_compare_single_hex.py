"""Phase 26: compare the user's HFSS run of `hfss_exports/hex150_tensor_cp_single` with `predictions.json`.

    python scripts/research/phase26_compare_single_hex.py --gain gain.csv [--s11 file.s1p|.s4p [--port N]]
        [--params antenna_params.txt] [--ar ar_cut.csv] [--ez ez_line.csv] [--pin scripted|connector|<label>]
        [--pred predictions.json] [--out report.json]

The logic lives in `src/hologram/tensor_cp_compare.py` (Phase 27 item 5); this script is a wrapper.
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.hologram import tensor_cp_compare as _m  # noqa: E402

globals().update({k: v for k, v in vars(_m).items() if not k.startswith("__")})  # keep the old script API (tests, other scripts)
main = _m.main_single

if __name__ == "__main__":
    _m.main_single()
