"""Phase 29: STEP of the compact 7-tile board `hfss_exports/phase29_compact7/compact7_board.step` (same exporter, `phase29_twofeed_step.build`, and settings as `sevenport_board.step`: 10.54 GHz design, 6 mm circle cells,
RT5880 1.575 mm, pins r 0.635 mm to z = 7.575 mm), geometry parsed from the package's own `build_board.py`; re-opened and checked; copied with `ports_and_excitations.csv` to C:\\hmsa_sync\\data\\phase29_compact7\\.
Also the raw-model per-port table and the theta <= 30 deg coverage (switched = best single port, blended = MRC over all 7 ports, uncoupled), no calibration margin. Writes results/phase29/compact7_step_summary.json."""

import csv
import json
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts/research"))
PKG = ROOT / "hfss_exports/phase29_compact7"
SYNC = Path(r"C:\hmsa_sync\data\phase29_compact7")


def main():
    import phase29_compact7_export as ce
    import phase29_compact7_tune as ct
    import phase29_twofeed_step as ts
    from phase29_pointB_sweep import _hpbw

    from src.hologram import seven_hex_cluster as c
    from src.hologram import shared_aperture_multiport as sa
    from src.hologram import tensor_cp_design as t
    from src.hologram.tensor_leaky_cp import SupercellTable

    ts.PKG = PKG
    cells, feeds, wx, wy = ts.parse("build_board.py")
    out = PKG / "compact7_board.step"
    ts.build(cells, feeds, wx, wy, out)
    names, ok = ts.check(out, len(cells), len(feeds))
    pins = {k: v for k, v in names.items() if k and k.startswith("Pin_")}
    # pin positions from the re-opened STEP are checked via the shared parse (same data written); here: counts and names
    SYNC.mkdir(parents=True, exist_ok=True)
    shutil.copy2(out, SYNC / out.name)

    R = float(json.loads((ROOT / "results/phase29/compact7_tune.json").read_text())["chosen_R_mm"])
    table, green = SupercellTable.from_files(t.TABLE_GLOB), c.green_function()
    th, ph = np.arange(0.0, 90.0 + 1e-9, 1.0), np.arange(0.0, 360.0, 1.0)
    _, F = ct.tile_fields(table, green, R, th, ph)
    rows = []
    for j, f in enumerate(F, start=1):
        beam = (0.0, 0.0) if j == 1 else (30.0, 60.0 * (j - 2))
        m = ce.metrics(*f, th, ph, beam, _hpbw, c._unit)
        rows.append(
            {
                "port": j,
                "feed_x_mm": round(feeds[j - 1][0], 3),
                "feed_y_mm": round(feeds[j - 1][1], 3),
                "design_beam_theta_deg": beam[0],
                "design_beam_phi_deg": beam[1],
                "pred_peak_theta_phi_deg": [round(x, 1) for x in m["peak_theta_phi_deg"]],
                "hpbw_deg": round(m["hpbw_deg"], 1),
                "gain_raw_dbi": round(m["peak_rhcp_gain_dbi"], 1),
                "sll_raw_db": round(m["sll_rel_db_outside_1p5_hpbw"], 1),
            }
        )
    with open(PKG / "ports_and_excitations.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(list(rows[0]) + [f"excite_P{k}" for k in range(1, 8)])
        for r in rows:
            w.writerow(list(r.values()) + [("1.0 W / 0 deg" if k == r["port"] else "0 (matched)") for k in range(1, 8)])
    shutil.copy2(PKG / "ports_and_excitations.csv", SYNC / "ports_and_excitations.csv")
    G = np.array([sa.gain_db(*f) for f in F])
    dirs = sa.scan_directions(30.0, 3.0, 15.0)
    sw, bl = [], []
    for d in dirs:
        it, ip = int(round(d[0])), int(round(d[1])) % 360
        gj = G[:, it, ip]
        sw.append(float(gj.max()))
        bl.append(float(10 * np.log10(np.sum(10 ** (gj / 10)))))
    sw, bl = np.array(sw), np.array(bl)
    cov = {
        "n_directions": len(dirs),
        "switched_min_mean_dbi": [round(float(sw.min()), 1), round(float(sw.mean()), 1)],
        "blended_min_mean_dbi": [round(float(bl.min()), 1), round(float(bl.mean()), 1)],
        "worst_switched_dir": list(dirs[int(np.argmin(sw))]),
        "worst_blended_dir": list(dirs[int(np.argmin(bl))]),
    }
    S = {"check_ok": ok, "cells": len(cells), "feeds": feeds, "names": {k: v for k, v in names.items() if k and not k.startswith(("Pin", "Port"))}, "n_pins": len(pins), "ports": rows, "coverage_le_30": cov}
    (ROOT / "results/phase29/compact7_step_summary.json").write_text(json.dumps(S, indent=1, default=float))
    print(json.dumps({k: S[k] for k in ("check_ok", "cells", "n_pins", "names", "coverage_le_30")}))
    for r in rows:
        print(r)
    return 0 if ok and len(pins) == 7 else 1


if __name__ == "__main__":
    sys.exit(main())
