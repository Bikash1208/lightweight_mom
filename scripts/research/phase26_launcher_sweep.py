"""Phase 26 (launcher design): sweep pin height and top-hat radius for the TM0 launch fraction at 10.54 GHz.

    python scripts/research/phase26_launcher_sweep.py <out.json>

For every `(top above the substrate surface, disk radius)`: launch fraction `P_sw / (P_sw + P_space)`, the
radiated-power resistance and thin-wire transmission-line reactance at the feed, the estimated reflection
against 50 ohm, and the launched share of the ACCEPTED power `launch * (1 - |Gamma|^2)`. Negative tops are
pins fully inside the substrate. The match is a thin-wire ESTIMATE (`x_est_ohm`), good to tens of ohms, not a
solved input impedance -- the HFSS calibration package is where it is measured.
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.dispersion.pin_launch import Launcher, PinLaunch  # noqa: E402

F, H, EPS_R, R_PIN = 10.54e9, 1.575e-3, 2.2, 12 * 25.4e-6
TOPS_MM = [-1.2, -0.6, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0]
DISKS_MM = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]


def evaluate(top_mm, disk_mm):
    pl = PinLaunch(F, EPS_R, H, Launcher(top_mm * 1e-3, R_PIN, disk_mm * 1e-3), n_seg=16)
    r = pl.split()
    z = complex(r["r_rad_ohm"], r["x_est_ohm"])
    gamma = abs((z - 50.0) / (z + 50.0))
    return {
        "top_mm": top_mm,
        "disk_mm": disk_mm,
        "launch_fraction": float(r["launch_fraction"]),
        "r_ohm": float(r["r_rad_ohm"]),
        "x_est_ohm": float(r["x_est_ohm"]),
        "gamma_db": float(20 * np.log10(max(gamma, 1e-9))),
        "launched_share_of_accepted": float(r["launch_fraction"] * (1.0 - gamma**2)),
    }


def resonant_top_mm(disk_mm):
    """Pin top (above the surface) at which the thin-wire reactance estimate crosses zero, or None."""
    from scipy.optimize import brentq

    def x(t):
        return PinLaunch(F, EPS_R, H, Launcher(t * 1e-3, R_PIN, disk_mm * 1e-3), n_seg=4)._reactance_estimate()

    grid = np.linspace(-1.0, 12.0, 53)
    vals = [x(t) for t in grid]
    out = []
    for a, b, fa, fb in zip(grid, grid[1:], vals, vals[1:], strict=False):
        if fa * fb < 0 and abs(fa) < 1e3 and abs(fb) < 1e3:  # a zero, not the cot() pole
            out.append(brentq(x, a, b))
    return out


def main(out):
    rows = [evaluate(t, a) for t in TOPS_MM for a in DISKS_MM if not (a > 0 and t < -0.5)]
    for a in (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0):
        for t in resonant_top_mm(a):
            r = evaluate(round(t, 3), a)
            r["resonant_tuned"] = True
            rows.append(r)
    cur_top = 110 * 25.4e-3 - H * 1e3  # PinLength from the ground
    cur = evaluate(round(cur_top, 3), 0.0)
    ok = [r for r in rows if r["gamma_db"] <= -10.0]
    best_match = max(ok, key=lambda r: r["launched_share_of_accepted"]) if ok else None
    best_fraction = max(rows, key=lambda r: r["launch_fraction"])
    best_share = max(rows, key=lambda r: r["launched_share_of_accepted"])
    res = {
        "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "MATCH_ESTIMATE_THIN_WIRE_TL"],
        "frequency_ghz": F / 1e9,
        "current_pin_from_ground_reading": cur,
        "best_launch_fraction_any_match": best_fraction,
        "best_launched_share_of_accepted_power": best_share,
        "best_launched_share_with_estimated_match_below_minus10_db": best_match,
        "rows": rows,
    }
    res["gain_recovery_db_vs_current_pin_launch_fraction"] = {
        "best_launch_fraction_any_match": float(10 * np.log10(best_fraction["launch_fraction"] / cur["launch_fraction"])),
        "best_launched_share_of_accepted_power": float(
            10 * np.log10(best_share["launch_fraction"] / cur["launch_fraction"])
        ),
    }
    Path(out).write_text(json.dumps(res, indent=1))
    print("current", cur)
    print("best fraction", best_fraction)
    print("best matched(<=-10dB)", best_match)
    print("best share", best_share)
    print(res["gain_recovery_db_vs_current_pin_launch_fraction"])
    for a in DISKS_MM:
        sub = [r for r in rows if r["disk_mm"] == a]
        b = max(sub, key=lambda r: r["launch_fraction"])
        print("disk", a, "max launch", round(b["launch_fraction"], 3), "at top", b["top_mm"], "G", round(b["gamma_db"], 1))
    for r in rows:
        if r.get("resonant_tuned"):
            print("resonant", r["disk_mm"], "top", r["top_mm"], "launch %.3f R %.1f X %.1f G %.1f" % (r["launch_fraction"], r["r_ohm"], r["x_est_ohm"], r["gamma_db"]))


if __name__ == "__main__":
    main(sys.argv[1])
