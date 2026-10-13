"""Phase 25 follow-up 2: what 1-bit separation does the empty|square pair reach when neighbouring patches must
keep a minimum gap, and how sensitive is it?

    python scripts/research/phase25_gap_pair_study.py <out.json>

Phase 22 `auto_states` (1-bit: empty + the square closest to empty + 180 deg, mesh-aligned solves,
Richardson-corrected, interpolated) at min gap 0.2 mm (the app default; the current pair), 1.0 and 1.1 mm,
at 10.33 / 10.54 / 10.75 GHz, period 6.4 mm; and, as the only other lever that keeps the structure, the
same at periods 7.0 and 8.0 mm with the 1.0 mm gap. Also the largest-separation square the gap allows
(the pair `auto_states` returns when the 180 deg level is unreachable).
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.app.storage import app_data_dir  # noqa: E402
from src.hologram.auto_states import auto_states, max_square_side_mm  # noqa: E402
from src.materials import Substrate  # noqa: E402
from src.phase import wrapped_difference_deg  # noqa: E402

SUB = Substrate("user", 2.2, 0.0009, 1.575e-3)
FREQS = (10.33e9, 10.54e9, 10.75e9)


def one(period_mm, gap, f, cache):
    a = auto_states(1, SUB, f, period_mm * 1e-3, min_gap_mm=gap, cache_dir=cache)
    s0, s1 = a.states
    sep = abs(wrapped_difference_deg(s0.phase_deg, s1.phase_deg))
    side = s1.geometry["L_mm"]
    return {
        "empty_phase_deg": s0.phase_deg,
        "square_side_mm": side,
        "square_phase_deg": s1.phase_deg,
        "square_magnitude_db": s1.magnitude_db,
        "separation_deg": float(sep),
        "gap_mm": period_mm - side,
        "max_side_allowed_mm": max_square_side_mm(period_mm * 1e-3, gap),
        "d_phase_per_0p01mm_deg": s1.sensitivity_deg_per_0p01mm,
        "d_phase_per_mm_deg": None if s1.sensitivity_deg_per_0p01mm is None else s1.sensitivity_deg_per_0p01mm * 100.0,
        "mesh_uncertainty_deg": s1.mesh_uncertainty_deg,
        "labels": list(s1.labels),
        "message": s1.message,
    }


def main(out):
    cache = app_data_dir() / "auto_states_cache"
    res = {"labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED"], "period_6p4": {}, "other_periods_gap_1p0": {}}
    for gap in (0.2, 1.0, 1.1):
        res["period_6p4"][f"min_gap_{gap}"] = {}
        for f in FREQS:
            r = one(6.4, gap, f, cache)
            res["period_6p4"][f"min_gap_{gap}"][f"{f / 1e9:.2f}GHz"] = r
            print(
                f"P6.4 gap>={gap} {f / 1e9:.2f} GHz: side {r['square_side_mm']:.2f} (gap {r['gap_mm']:.2f}) "
                f"sep {r['separation_deg']:.1f} deg, d/dL {r['d_phase_per_mm_deg']:.0f} deg/mm",
                flush=True,
            )
        Path(out).write_text(json.dumps(res, indent=1))
    for period in (7.0, 8.0):
        res["other_periods_gap_1p0"][f"P{period}"] = {}
        for f in FREQS:
            r = one(period, 1.0, f, cache)
            res["other_periods_gap_1p0"][f"P{period}"][f"{f / 1e9:.2f}GHz"] = r
            print(
                f"P{period} gap>=1.0 {f / 1e9:.2f} GHz: side {r['square_side_mm']:.2f} (gap {r['gap_mm']:.2f}) "
                f"sep {r['separation_deg']:.1f} deg, d/dL {r['d_phase_per_mm_deg']:.0f} deg/mm",
                flush=True,
            )
        Path(out).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
