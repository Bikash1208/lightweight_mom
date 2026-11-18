"""Interpolate unrotated PML Q to the two folded per-cell phases."""

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
R = ROOT / "results" / "sub6_meander_model" / "supercell_hfss"


def modes(cell, direction, phase):
    path = R / f"{cell}_PML_{direction}_{phase}_air40.035_modes.json"
    return [m for m in json.loads(path.read_text(encoding="utf-8"))["modes"]
            if m["Q"] < 100000]


def lerp(a, b, phase):
    return a + (b - a) * ((phase - 85) / 10)


out = {"method": "arithmetic mean of across/along unrotated PML Q, interpolated 85-95 deg",
       "air_top_mm": 40.035, "rows": []}
for cell in ("ANCHOR", "B2p2_D3p0", "B2p2_D3p6"):
    data = {(d, p): modes(cell, d, p) for d in ("across", "along") for p in (85, 95)}
    for phase_s in (15, 30, 45):
        for sign in (-1, 1):
            point = (720 + sign * phase_s) / 8
            for branch in (0, 1):
                q, f = {}, {}
                for d in ("across", "along"):
                    a, b = data[d, 85][branch], data[d, 95][branch]
                    q[d] = lerp(a["Q"], b["Q"], point)
                    f[d] = lerp(a["f_ghz"], b["f_ghz"], point)
                out["rows"].append({"cell": cell, "super_phase_deg": phase_s,
                                    "fold_sign": sign, "per_cell_phase_deg": point,
                                    "branch_index": branch + 1,
                                    "Q_across": q["across"], "Q_along": q["along"],
                                    "Q_base": (q["across"] + q["along"]) / 2,
                                    "f_across_ghz": f["across"],
                                    "f_along_ghz": f["along"],
                                    "phase_extrapolated": not 85 <= point <= 95})

(R / "folded_qbase.json").write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
print(f"wrote {len(out['rows'])} folded baseline rows")
